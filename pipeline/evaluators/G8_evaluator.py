"""
G8 evaluator — reporting traceability and release lock.

Implements contracts/G8_reporting.yaml. This is the last gate before anything
leaves the system, so it is deliberately the strictest: every figure, table and
numerical claim must resolve to a registered artifact whose SHA-256 still
matches and whose contract passed (or was explicitly adjudicated by a human).

The hash comparison is what makes a stale figure or a hand-edited results table
detectable. A figure regenerated from different data, or a value edited after
computation, changes the file digest and therefore stops matching the artifact
it claims to come from — the manifest alone cannot reveal that.
"""

from __future__ import annotations

import hashlib
import re
from pathlib import Path
from typing import Any

from ._assay_blindness import BLOCKED_SUBSTRINGS, is_blocked_column
from ._confidence_tier import assess_lbdd, assess_sbdd, combine

# Metric names that may never be reported: these are the SNP-Lib-internal
# figures derived from the excluded in-house assay (ADR-AIIM-004). Only the
# ChEMBL external validation metrics are reportable, under G5b caveat rules.
_FORBIDDEN_METRIC_PATTERNS = (
    r"internal[_\s-]*auc",
    r"in[-\s]?house[_\s-]*auc",
    r"camp[_\s-]*auc",
    r"snp[_\s-]*lib[_\s-]*auc",
    r"assay[_\s-]*auc",
    r"retrospective[_\s-]*auc",
)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest().upper()


def _registry_artifacts(registry: Any) -> dict[str, dict[str, Any]]:
    """artifact_id -> latest artifact record, from a registry or a snapshot list."""
    records: list[dict[str, Any]] = []
    if registry is None:
        return {}
    if isinstance(registry, list):
        records = registry
    elif hasattr(registry, "get_all_records"):
        try:
            records = registry.get_all_records()
        except Exception:
            return {}
    out: dict[str, dict[str, Any]] = {}
    for rec in records:
        if rec.get("record_type") == "artifact":
            aid = rec.get("artifact_id") or rec.get("record_id")
            if aid:
                out[str(aid)] = rec
    return out


def _released_artifact_ids(registry: Any) -> set[str]:
    """Artifacts with a PASS contract outcome or a human adjudication."""
    records: list[dict[str, Any]] = []
    if registry is None:
        return set()
    if isinstance(registry, list):
        records = registry
    elif hasattr(registry, "get_all_records"):
        try:
            records = registry.get_all_records()
        except Exception:
            return set()
    released: set[str] = set()
    for rec in records:
        rtype = rec.get("record_type")
        if rtype == "contract_outcome" and rec.get("outcome") == "PASS":
            released.update(str(a) for a in rec.get("artifact_ids", []))
        elif rtype == "human_intervention":
            released.update(str(a) for a in rec.get("artifact_ids", []))
    return released


def _verify_manifest(
    manifest: list[dict[str, Any]],
    artifacts: dict[str, dict[str, Any]],
    id_key: str,
) -> tuple[list[dict[str, Any]], int]:
    """Return (problems, n_verified) for a figure or table manifest."""
    problems: list[dict[str, Any]] = []
    verified = 0
    for item in manifest or []:
        item_id = item.get(id_key, "<unnamed>")
        artifact_id = item.get("artifact_id")
        path_str = item.get("path") or item.get("source_path")

        if not artifact_id:
            problems.append({id_key: item_id, "problem": "no artifact_id declared"})
            continue

        expected = (item.get("expected_sha256")
                    or artifacts.get(str(artifact_id), {}).get("sha256"))
        if not expected:
            problems.append({id_key: item_id, "artifact_id": artifact_id,
                             "problem": "artifact not registered / no recorded sha256"})
            continue

        if not path_str:
            problems.append({id_key: item_id, "artifact_id": artifact_id,
                             "problem": "no path to verify"})
            continue

        path = Path(str(path_str))
        if not path.is_file():
            problems.append({id_key: item_id, "artifact_id": artifact_id,
                             "problem": f"file missing: {path}"})
            continue

        actual = _sha256(path)
        if actual != str(expected).upper():
            problems.append({
                id_key: item_id,
                "artifact_id": artifact_id,
                "problem": "sha256 mismatch — content differs from registered artifact",
                "expected": str(expected).upper(),
                "actual": actual,
            })
            continue
        verified += 1
    return problems, verified


def evaluate(
    contract: dict[str, Any],
    inputs: dict[str, Any],
    registry: Any,
) -> dict[str, Any]:
    checks = contract.get("checks", [])
    passed, failed, warned, evidence = [], [], [], {}
    critical_failures: list[str] = []

    def fail(check_id: str, severity: str, detail: dict[str, Any]) -> None:
        failed.append(check_id)
        if severity == "critical":
            critical_failures.append(check_id)
        evidence[check_id] = detail

    snapshot = inputs.get("artifact_registry_snapshot", registry)
    artifacts = _registry_artifacts(snapshot)
    released = _released_artifact_ids(snapshot)

    figures = inputs.get("figure_manifest") or inputs.get("figure_manifest_json") or []
    tables = inputs.get("table_manifest") or inputs.get("table_manifest_json") or []
    snippets = (inputs.get("result_snippets")
                or inputs.get("manuscript_result_snippets") or [])
    metadata = inputs.get("manuscript_metadata") or {}

    for check in checks:
        cid = check["id"]
        sev = check.get("severity", "medium")

        if cid == "number_traceability":
            untraceable = []
            for snip in snippets:
                claim = snip.get("claim", "<unnamed claim>")
                aid = snip.get("artifact_id")
                if not aid:
                    untraceable.append({"claim": claim, "problem": "no artifact_id"})
                elif str(aid) not in artifacts:
                    untraceable.append({"claim": claim, "artifact_id": aid,
                                        "problem": "artifact not in registry"})
            if untraceable:
                fail(cid, sev, {"n_untraceable": len(untraceable),
                                "untraceable": untraceable[:20]})
            else:
                passed.append(cid)
                evidence[cid] = {"n_claims": len(snippets), "all_traceable": True}

        elif cid == "contract_coverage":
            cited = {str(s["artifact_id"]) for s in snippets if s.get("artifact_id")}
            cited |= {str(f["artifact_id"]) for f in figures if f.get("artifact_id")}
            cited |= {str(t["artifact_id"]) for t in tables if t.get("artifact_id")}
            uncovered = sorted(cited - released)
            if uncovered:
                # No PASS record and no human adjudication: the report lock must
                # not release these.
                fail(cid, sev, {"n_uncovered": len(uncovered),
                                "uncovered_artifacts": uncovered[:20]})
            else:
                passed.append(cid)
                evidence[cid] = {"n_cited_artifacts": len(cited), "all_covered": True}

        elif cid == "figure_source_verification":
            problems, verified = _verify_manifest(figures, artifacts, "figure_id")
            if problems:
                fail(cid, sev, {"n_problems": len(problems), "problems": problems[:20]})
            else:
                passed.append(cid)
                evidence[cid] = {"n_figures": len(figures), "n_verified": verified}

        elif cid == "table_source_verification":
            problems, verified = _verify_manifest(tables, artifacts, "table_id")
            if problems:
                fail(cid, sev, {"n_problems": len(problems), "problems": problems[:20]})
            else:
                passed.append(cid)
                evidence[cid] = {"n_tables": len(tables), "n_verified": verified}

        elif cid == "no_unregistered_claims":
            # A claim whose recorded value disagrees with the value recomputed
            # from its source artifact indicates a post-computation edit.
            mismatches = []
            for snip in snippets:
                if "value" not in snip or "expected_value" not in snip:
                    continue
                try:
                    reported = float(snip["value"])
                    source = float(snip["expected_value"])
                except (TypeError, ValueError):
                    if str(snip["value"]) != str(snip["expected_value"]):
                        mismatches.append({"claim": snip.get("claim"),
                                           "reported": snip["value"],
                                           "source": snip["expected_value"]})
                    continue
                tol = float(snip.get("tolerance", 1e-9))
                if abs(reported - source) > tol:
                    mismatches.append({"claim": snip.get("claim"),
                                       "reported": reported, "source": source,
                                       "delta": round(reported - source, 10)})
            if mismatches:
                fail(cid, sev, {"n_mismatches": len(mismatches),
                                "mismatches": mismatches[:20]})
            else:
                passed.append(cid)
                evidence[cid] = {"n_checked": len(snippets), "n_mismatches": 0}

        elif cid == "internal_auc_not_cited":
            hits = []
            for snip in snippets:
                text = " ".join(
                    str(snip.get(k, "")) for k in ("claim", "metric", "name", "label")
                ).lower()
                for pattern in _FORBIDDEN_METRIC_PATTERNS:
                    if re.search(pattern, text):
                        hits.append({"claim": snip.get("claim"), "matched": pattern})
                        break
                else:
                    # Also block any claim keyed on an assay-derived column name.
                    metric_name = str(snip.get("metric") or snip.get("name") or "")
                    if metric_name and is_blocked_column(metric_name):
                        hits.append({"claim": snip.get("claim"),
                                     "matched": f"blocked metric name: {metric_name}"})
            if hits:
                fail(cid, sev, {"n_forbidden_metrics": len(hits), "hits": hits[:20],
                                "policy": "ADR-AIIM-004: in-house assay metrics are "
                                          "excluded from all reported results"})
            else:
                passed.append(cid)
                evidence[cid] = {"n_claims_scanned": len(snippets),
                                 "forbidden_patterns": len(_FORBIDDEN_METRIC_PATTERNS)}

        elif cid == "human_intervention_count":
            count = metadata.get("human_intervention_count")
            if count is None and hasattr(snapshot, "count_human_interventions"):
                try:
                    count = snapshot.count_human_interventions()
                except Exception:
                    count = None
            if count is None:
                (failed if sev == "critical" else warned).append(cid)
                if sev == "critical":
                    critical_failures.append(cid)
                evidence[cid] = {"error": "human_intervention_count not reported"}
            else:
                passed.append(cid)
                evidence[cid] = {"human_intervention_count": int(count)}

        elif cid == "confidence_tier_assessment":
            ct = inputs.get("confidence_tier_inputs")
            if not ct or not ct.get("paths_attempted"):
                warned.append(cid)
                evidence[cid] = {"warning": "no confidence_tier_inputs supplied; "
                                            "no technique path was assessed"}
            else:
                # Frozen contract values (ADR-AIIM-046) are the authoritative
                # thresholds, not the Python module's defaults -- mirrors
                # G3's control_representativeness check.get(...) pattern.
                sbdd_kwargs = dict(ct.get("sbdd") or {})
                lbdd_kwargs = dict(ct.get("lbdd") or {})
                paths = ct["paths_attempted"]
                sbdd_report = (
                    assess_sbdd(
                        **sbdd_kwargs,
                        pose_pass_threshold=check.get("sbdd_pose_pass_threshold_angstrom", 2.0),
                        pose_margin_threshold=check.get("sbdd_pose_margin_threshold_angstrom", 1.0),
                        min_similarity_coverage=check.get("sbdd_min_similarity_coverage", 0.50),
                        min_scaffold_coverage=check.get("sbdd_min_scaffold_coverage", 0.10),
                        external_auc_min=check.get("sbdd_external_auc_min", 0.60),
                        external_ef1pct_min=check.get("sbdd_external_ef1pct_min", 3.0),
                        external_bedroc20_min=check.get("sbdd_external_bedroc20_min", 0.40),
                    ) if "SBDD" in paths and ct.get("sbdd") else None
                )
                lbdd_report = (
                    assess_lbdd(
                        **lbdd_kwargs,
                        ad_coverage_high=check.get("lbdd_ad_coverage_high", 0.80),
                        ad_coverage_low=check.get("lbdd_ad_coverage_low", 0.40),
                        cv_auc_high=check.get("lbdd_cv_auc_high", 0.75),
                        cv_auc_low=check.get("lbdd_cv_auc_low", 0.60),
                        external_auc_min=check.get("lbdd_external_auc_min", 0.60),
                        imprecision_min_recommended=check.get("lbdd_imprecision_min_recommended", 30),
                        imprecision_min_floor=check.get("lbdd_imprecision_min_floor", 10),
                    ) if "LBDD" in paths and ct.get("lbdd") else None
                )
                overall = combine(sbdd_report, lbdd_report)
                passed.append(cid)
                evidence[cid] = {
                    "overall_tier": overall["overall_tier"],
                    "sbdd": sbdd_report, "lbdd": lbdd_report,
                }

        elif cid == "confidence_tier_threshold_immutability":
            registered = inputs.get("registered_confidence_tier_thresholds")
            # Thresholds live on the confidence_tier_assessment check block,
            # not on this check's own (threshold-free) block.
            assessment_check = next(
                (c for c in checks if c.get("id") == "confidence_tier_assessment"), {}
            )
            current = {k: v for k, v in assessment_check.items()
                       if k not in ("id", "description", "severity")}
            if registered is None:
                warned.append(cid)
                evidence[cid] = {"warning": "no registered anchor to compare against yet "
                                            "(first run) — record this contract's values "
                                            "before the next run"}
            elif registered != current:
                fail(cid, sev, {"registered": registered, "current": current,
                                "error": "confidence-tier thresholds changed since registration"})
            else:
                passed.append(cid)
                evidence[cid] = {"thresholds": current, "match": True}

        elif cid == "ai_use_declaration_present":
            decl = metadata.get("ai_use_declaration")
            if not decl or not str(decl).strip():
                (failed if sev == "critical" else warned).append(cid)
                if sev == "critical":
                    critical_failures.append(cid)
                evidence[cid] = {"error": "ai_use_declaration missing or empty"}
            else:
                passed.append(cid)
                evidence[cid] = {"declaration_length": len(str(decl))}

    if critical_failures:
        state = "HALT"
    elif failed:
        state = "ESCALATE"
    else:
        state = "PASS"

    return {
        "state": state,
        "passed": passed,
        "failed": failed,
        "warned": warned,
        "critical_failures": critical_failures,
        "evidence": evidence,
        "summary": (
            f"{len(passed)} checks passed, {len(failed)} failed "
            f"({len(critical_failures)} critical), {len(warned)} warned"
        ),
    }

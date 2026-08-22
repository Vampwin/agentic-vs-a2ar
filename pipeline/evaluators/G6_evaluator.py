"""
G6 evaluator — selection assay-blindness and frozen rule verification.

Scans column names of every dataframe read during the selection step.
Any blocked column name (or substring match) is a critical failure.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import pandas as pd

from . import _domain_gate_bias_probe

# Exact column names that must never appear
_BLOCKED_EXACT = frozenset({
    "known_active", "a_sim_actives", "pct_inhibition",
    "cAMP", "assay_label", "is_active",
})

# Substrings that must not appear in any column name (case-insensitive)
_BLOCKED_SUBSTRINGS = ["camp", "assay", "inhibit"]

# SNP-Lib cAMP-positive compound IDs (ADR-AIIM-004).
# If any selection code path references these as explicit targets (not incidentally
# as high-ranked compounds), it is a failure.
_CAMP_ACTIVE_IDS = frozenset({
    "SNP00183", "SNP00182", "SNP00063", "SNP00227", "SNP00235",
})


def evaluate(
    contract: dict[str, Any],
    inputs: dict[str, Any],
    registry: Any,
) -> dict[str, Any]:
    passed, failed, warned, evidence = [], [], [], {}
    critical_failures = []

    ranking_csv = Path(inputs.get("final_ranking_csv", ""))
    selected_csv = Path(inputs.get("selected_candidates_csv", ""))
    column_scan_log = inputs.get("column_scan_log", [])  # List of {file, columns} dicts
    selection_script_hash = inputs.get("selection_script_hash", "")
    registered_script_hash = inputs.get("registered_selection_script_hash", "")

    for check in contract.get("checks", []):
        cid = check["id"]
        sev = check.get("severity", "medium")

        if cid == "assay_blindness_column_scan":
            blocked = _scan_columns(column_scan_log, ranking_csv, selected_csv)
            ok = len(blocked) == 0
            evidence[cid] = {"blocked_columns_found": blocked, "files_scanned": len(column_scan_log)}
            (passed if ok else failed).append(cid)
            if not ok and sev == "critical":
                critical_failures.append(cid)

        elif cid == "eligible_set_completeness":
            # A selection rule may reorder the eligible set but may not quietly
            # shrink it. An undeclared filter applied before ranking removes
            # candidates that the frozen rule never authorised excluding, and
            # because the surviving rows are internally consistent nothing
            # downstream can notice. Every eligible compound must therefore
            # appear in the ranking or carry a recorded exclusion reason.
            expected = inputs.get("expected_eligible_count", check.get("expected_count"))
            excluded = inputs.get("selection_exclusions") or {}
            if not ranking_csv.is_file():
                failed.append(cid)
                if sev == "critical":
                    critical_failures.append(cid)
                evidence[cid] = {"error": f"ranking file not found: {ranking_csv}"}
            elif expected is None:
                warned.append(cid)
                evidence[cid] = {"warning": "no expected eligible count supplied"}
            else:
                try:
                    n_ranked = len(pd.read_csv(ranking_csv))
                except Exception as exc:
                    n_ranked = None
                    evidence[cid] = {"error": f"could not read ranking: {exc}"}
                if n_ranked is None:
                    failed.append(cid)
                    if sev == "critical":
                        critical_failures.append(cid)
                else:
                    accounted = n_ranked + len(excluded)
                    ok = accounted == int(expected)
                    evidence[cid] = {
                        "expected_eligible": int(expected),
                        "n_ranked": n_ranked,
                        "n_recorded_exclusions": len(excluded),
                        "n_unaccounted": int(expected) - accounted,
                    }
                    (passed if ok else failed).append(cid)
                    if not ok and sev == "critical":
                        critical_failures.append(cid)

        elif cid == "scaffold_uniqueness":
            # The panel exists to sample distinct chemotypes. Two members
            # sharing a Bemis-Murcko scaffold spend a slot of a deliberately
            # small panel on a near-duplicate, which silently narrows the
            # chemical space actually carried forward.
            if not selected_csv.is_file():
                failed.append(cid)
                if sev == "critical":
                    critical_failures.append(cid)
                evidence[cid] = {"error": f"selected panel not found: {selected_csv}"}
            else:
                dupes, n_panel, unparsed = _scaffold_duplicates(selected_csv)
                if dupes is None:
                    warned.append(cid)
                    evidence[cid] = {"warning": "RDKit unavailable; scaffolds not computed"}
                else:
                    ok = len(dupes) == 0
                    evidence[cid] = {"n_panel": n_panel,
                                     "n_duplicate_scaffolds": len(dupes),
                                     "duplicate_scaffolds": dupes[:10],
                                     "n_unparseable": unparsed}
                    (passed if ok else failed).append(cid)
                    if not ok and sev == "critical":
                        critical_failures.append(cid)

        elif cid == "frozen_rule_applied":
            if selection_script_hash and registered_script_hash:
                ok = selection_script_hash.upper() == registered_script_hash.upper()
                evidence[cid] = {
                    "registered_hash": registered_script_hash,
                    "actual_hash": selection_script_hash,
                    "match": ok,
                }
                (passed if ok else failed).append(cid)
                if not ok and sev == "critical":
                    critical_failures.append(cid)
            else:
                warned.append(cid)
                evidence[cid] = {"warning": "Script hashes not provided; manual verification required"}

        elif cid == "no_assay_ids_in_selected":
            # Check whether selection code explicitly targeted cAMP-active IDs
            explicit_target = inputs.get("selection_explicitly_targeted_assay_ids", False)
            if explicit_target:
                failed.append(cid)
                critical_failures.append(cid)
                evidence[cid] = {
                    "explicit_target": True,
                    "note": "Selection code referenced cAMP-positive IDs by ID — protocol violation",
                }
            else:
                # Check if any cAMP-active IDs appear in the selected set (acceptable if by rank)
                selected_ids = _read_id_column(selected_csv)
                camp_selected = _CAMP_ACTIVE_IDS & set(selected_ids)
                passed.append(cid)
                evidence[cid] = {
                    "explicit_target": False,
                    "camp_ids_in_selected": sorted(camp_selected),
                    "note": "These IDs may appear because they rank highly by Vina; acceptable",
                }

        elif cid == "diversity_seed_recorded":
            seed = inputs.get("diversity_seed")
            ok = seed is not None
            evidence[cid] = {"diversity_seed": seed}
            (passed if ok else warned).append(cid)

        elif cid == "selection_count_reported":
            count_report = inputs.get("selection_count_report")
            ok = count_report is not None
            evidence[cid] = {"count_report": count_report}
            (passed if ok else warned).append(cid)

        elif cid == "applicability_domain":
            ad_flag = inputs.get("applicability_domain_annotation_present", False)
            coverage = inputs.get("applicability_domain_coverage_fraction")
            evidence[cid] = {"annotation_present": ad_flag}
            if coverage is not None:
                # Real number for the G8 confidence-tier assessment to rate,
                # per the reference metric/threshold declared in this
                # contract (ADR-AIIM-045); this check's own pass/warn status
                # is still about whether the annotation exists, not the
                # coverage value itself -- the tier judgment happens in G8.
                evidence[cid]["coverage_fraction"] = coverage
                evidence[cid]["reference_metric"] = check.get("ad_reference_metric")
                evidence[cid]["reference_threshold"] = check.get("ad_reference_threshold")
            (passed if ad_flag else warned).append(cid)

        elif cid == "applicability_domain_multi_space":
            coverage_by_space = inputs.get("applicability_domain_by_space")
            min_spaces = int(check.get("min_spaces", _domain_gate_bias_probe.DEFAULT_MIN_SPACES))
            report = _domain_gate_bias_probe.assess(coverage_by_space or {}, min_spaces=min_spaces)
            evidence[cid] = report
            if not report.get("available"):
                warned.append(cid)
            elif not report.get("enough_spaces"):
                # Fewer spaces than required is a coverage gap in the evidence
                # itself, not a finding about the library -- warn, don't fail.
                warned.append(cid)
            else:
                # This check exists to surface the divergence, never to hide
                # it: reaching here with enough spaces is a pass regardless of
                # whether divergence was found, because the report itself
                # (in evidence) carries that finding forward to the human
                # decision at G6/G8.
                passed.append(cid)

        elif cid == "drug_likeness_not_used_in_ranking":
            # Check that ranking_csv columns used for selection do not include DL metrics
            ranking_cols = _get_column_names(ranking_csv)
            dl_in_ranking = [
                c for c in ranking_cols
                if any(kw in c.lower() for kw in ["lipinski", "veber", "qed", "pains"])
                and "score" in c.lower()
            ]
            ok = len(dl_in_ranking) == 0
            evidence[cid] = {"dl_ranking_columns": dl_in_ranking}
            (passed if ok else failed).append(cid)
            if not ok and sev == "critical":
                critical_failures.append(cid)

    # Critical breaches halt; non-critical failures still require a human
    # decision. Returning PASS while `failed` is non-empty would let a
    # high-severity finding through silently.
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


def _scaffold_duplicates(
    panel_csv: Path,
) -> tuple[list[dict[str, Any]] | None, int, int]:
    """
    Bemis-Murcko scaffolds shared by more than one panel member.

    Returns (duplicates, n_panel, n_unparseable); duplicates is None when RDKit
    is unavailable, so the caller can warn rather than assert a clean panel it
    could not actually verify.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem.Scaffolds import MurckoScaffold
    except ImportError:
        return None, 0, 0

    try:
        df = pd.read_csv(panel_csv)
    except Exception:
        return [], 0, 0

    smiles_col = next((c for c in df.columns if "smiles" in c.lower()), None)
    id_col = next((c for c in df.columns if "id" in c.lower()), None)
    if smiles_col is None:
        return [], len(df), len(df)

    groups: dict[str, list[str]] = {}
    unparsed = 0
    for _, row in df.iterrows():
        smi = str(row.get(smiles_col, "") or "")
        mol = Chem.MolFromSmiles(smi) if smi else None
        if mol is None:
            unparsed += 1
            continue
        try:
            scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
        except Exception:
            unparsed += 1
            continue
        if not scaffold:
            continue
        label = str(row.get(id_col, row.name)) if id_col else str(row.name)
        groups.setdefault(scaffold, []).append(label)

    duplicates = [
        {"scaffold": sc, "members": members}
        for sc, members in groups.items()
        if len(members) > 1
    ]
    return duplicates, len(df), unparsed


def _scan_columns(
    scan_log: list[dict[str, Any]],
    *additional_csvs: Path,
) -> list[dict[str, str]]:
    """Return list of {file, column} for all blocked column names found."""
    findings = []
    all_entries = list(scan_log)
    for p in additional_csvs:
        if p and p.exists():
            all_entries.append({"file": str(p), "columns": _get_column_names(p)})

    for entry in all_entries:
        file_name = entry.get("file", "unknown")
        for col in entry.get("columns", []):
            col_lower = col.lower()
            if col_lower in {b.lower() for b in _BLOCKED_EXACT}:
                findings.append({"file": file_name, "column": col, "match_type": "exact"})
            elif any(sub in col_lower for sub in _BLOCKED_SUBSTRINGS):
                findings.append({"file": file_name, "column": col, "match_type": "substring"})
    return findings


def _get_column_names(csv_path: Path) -> list[str]:
    if not csv_path or not csv_path.exists():
        return []
    try:
        import pandas as pd
        df = pd.read_csv(csv_path, nrows=0)
        return list(df.columns)
    except Exception:
        return []


def _read_id_column(csv_path: Path) -> list[str]:
    if not csv_path or not csv_path.exists():
        return []
    try:
        import pandas as pd
        df = pd.read_csv(csv_path)
        id_col = next(
            (c for c in df.columns if c.lower() in {"snp_id", "compound_id", "id"}),
            None,
        )
        return df[id_col].astype(str).tolist() if id_col else []
    except Exception:
        return []

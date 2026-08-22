"""
G2 evaluator — receptor and ligand preparation.

Implements contracts/G2_preparation.yaml: receptor integrity, reference-ligand
extraction, docking-box provenance, ligand embedding accounting, and the
requirement that the box be derived from crystallographic geometry alone.

The box checks matter more than they look. A docking box silently recentred on
a compound of interest turns the screen into a self-fulfilling prophecy, and
nothing downstream can detect it — the scores are internally consistent and the
poses look reasonable. So the box must declare which structure it came from and
must not mention any library compound.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

# Any reference to a library compound ID inside a box definition means the box
# was not derived from crystallographic geometry alone.
_SNP_ID_PATTERN = re.compile(r"SNP\s*_?\d{3,}", re.IGNORECASE)

_REQUIRED_BOX_FIELDS = ("center", "size")
_PROVENANCE_FIELDS = ("pdb_used", "source_pdb", "box_center_residue_ref",
                      "centroid_reference", "derived_from")


def _load_json(path_str: Any) -> tuple[dict[str, Any] | None, str | None]:
    if not path_str:
        return None, "no path supplied"
    path = Path(str(path_str))
    if not path.is_file():
        return None, f"file not found: {path}"
    try:
        with open(path, encoding="utf-8") as f:
            return json.load(f), None
    except Exception as exc:
        return None, f"{type(exc).__name__}: {exc}"


def _count_residues(pdbqt_path: Path) -> tuple[int, list[str]]:
    """Return (residue count, ordered chain+resseq keys) from a PDBQT/PDB file."""
    seen: list[str] = []
    unique: set[str] = set()
    with open(pdbqt_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                chain = line[21:22]
                resseq = line[22:27].strip()
                key = f"{chain}:{resseq}"
                if key not in unique:
                    unique.add(key)
                    seen.append(key)
    return len(unique), seen


def _count_structural_waters(receptor_path: Path) -> int:
    """Count distinct HOH/WAT residues (by chain+resseq) in a PDB/PDBQT file."""
    seen: set[tuple[str, str]] = set()
    with open(receptor_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(("ATOM  ", "HETATM")):
                resn = line[17:20].strip()
                if resn in {"HOH", "WAT"}:
                    seen.add((line[21:22], line[22:26].strip()))
    return len(seen)


def _has_chain_break(residue_keys: list[str]) -> list[str]:
    """Detect numbering discontinuities within a chain (missing residues)."""
    breaks: list[str] = []
    prev_chain, prev_num = None, None
    for key in residue_keys:
        chain, _, resseq = key.partition(":")
        try:
            num = int(re.sub(r"[^\d-]", "", resseq) or 0)
        except ValueError:
            continue
        if chain == prev_chain and prev_num is not None and num - prev_num > 1:
            breaks.append(f"{chain}:{prev_num}->{num}")
        prev_chain, prev_num = chain, num
    return breaks


def _check_receptor(
    path_str: Any,
    check: dict[str, Any],
) -> tuple[bool, dict[str, Any]]:
    if not path_str:
        return False, {"error": "no receptor path supplied"}
    path = Path(str(path_str))
    if not path.is_file():
        return False, {"error": f"file not found: {path}"}

    n_res, keys = _count_residues(path)
    detail: dict[str, Any] = {"path": str(path), "residue_count": n_res}

    expected = check.get("expected_residue_count")
    if expected is not None and n_res != int(expected):
        detail["expected_residue_count"] = int(expected)
        detail["error"] = "residue count differs from expected"
        return False, detail

    if n_res == 0:
        detail["error"] = "no ATOM/HETATM records parsed"
        return False, detail

    breaks = _has_chain_break(keys)
    detail["chain_breaks"] = breaks[:20]
    detail["n_chain_breaks"] = len(breaks)
    if breaks and check.get("allow_chain_breaks") is False:
        detail["error"] = "chain breaks present"
        return False, detail
    return True, detail


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

    box, box_error = _load_json(inputs.get("docking_box_json"))

    for check in checks:
        cid = check["id"]
        sev = check.get("severity", "medium")

        if cid == "receptor_integrity_4eiy":
            ok, detail = _check_receptor(inputs.get("prepared_receptor_4eiy_pdbqt"), check)
            if ok:
                passed.append(cid)
                evidence[cid] = detail
            else:
                fail(cid, sev, detail)

        elif cid == "receptor_integrity_3eml":
            ok, detail = _check_receptor(inputs.get("prepared_receptor_3eml_pdbqt"), check)
            if ok:
                passed.append(cid)
                evidence[cid] = detail
            else:
                fail(cid, sev, detail)

        elif cid == "reference_ligand_present":
            ref = inputs.get("reference_ligand_4eiy_sdf")
            path = Path(str(ref)) if ref else None
            if not path or not path.is_file():
                fail(cid, sev, {"error": f"reference ligand not found: {ref}"})
            elif path.stat().st_size == 0:
                fail(cid, sev, {"error": "reference ligand file is empty", "path": str(path)})
            else:
                passed.append(cid)
                evidence[cid] = {"path": str(path), "size_bytes": path.stat().st_size}

        elif cid == "box_derivation_recorded":
            if box is None:
                fail(cid, sev, {"error": box_error or "docking box JSON unavailable"})
                continue
            missing = [f for f in _REQUIRED_BOX_FIELDS if f not in box]
            provenance = [f for f in _PROVENANCE_FIELDS if box.get(f)]
            detail = {"fields_present": sorted(box.keys()),
                      "missing_required": missing,
                      "provenance_fields": provenance}
            if missing:
                detail["error"] = "box definition missing required fields"
                fail(cid, sev, detail)
            elif not provenance:
                # Without a recorded source structure the box cannot be shown to
                # be crystallographically derived rather than chosen.
                detail["error"] = ("no provenance field recording which structure "
                                   "the box was derived from")
                fail(cid, sev, detail)
            else:
                passed.append(cid)
                evidence[cid] = detail

        elif cid == "no_assay_constants_in_box":
            if box is None:
                fail(cid, sev, {"error": box_error or "docking box JSON unavailable"})
                continue
            blob = json.dumps(box)
            hits = sorted(set(_SNP_ID_PATTERN.findall(blob)))
            if hits:
                fail(cid, sev, {"compound_ids_referenced": hits[:20],
                                "error": "box definition references library compounds; "
                                         "it must derive from crystallographic geometry only"})
            else:
                passed.append(cid)
                evidence[cid] = {"compound_ids_referenced": []}

        elif cid == "embedding_count":
            n_embedded = inputs.get("n_embedded_successfully")
            expected = inputs.get("count_valid_for_docking", check.get("expected_count"))
            if n_embedded is None or expected is None:
                warned.append(cid)
                evidence[cid] = {"warning": "embedding or expected count not supplied",
                                 "n_embedded": n_embedded, "expected": expected}
            elif int(n_embedded) == int(expected):
                passed.append(cid)
                evidence[cid] = {"n_embedded": int(n_embedded), "expected": int(expected)}
            else:
                fail(cid, sev, {"n_embedded": int(n_embedded), "expected": int(expected),
                                "delta": int(n_embedded) - int(expected)})

        elif cid == "embedding_failure_accounting":
            failures = inputs.get("embedding_failures")
            if failures is None:
                warned.append(cid)
                evidence[cid] = {"warning": "embedding_failures not supplied"}
                continue
            unexplained = [
                str(k) for k, v in dict(failures).items()
                if not str(v or "").strip()
            ]
            if unexplained:
                fail(cid, sev, {"n_unexplained": len(unexplained),
                                "unexplained_sample": unexplained[:20],
                                "error": "embedding failures lack a recorded reason"})
            else:
                passed.append(cid)
                evidence[cid] = {"n_failures": len(dict(failures)), "all_explained": True}

        elif cid == "cofactor_decision_recorded":
            decisions = inputs.get("cofactor_decisions")
            valid = {"retained", "removed", "modified"}
            if not decisions:
                (failed if sev == "critical" else warned).append(cid)
                if sev == "critical":
                    critical_failures.append(cid)
                evidence[cid] = {"warning" if sev != "critical" else "error":
                                 "no cofactor decisions recorded"}
                continue
            problems = [
                {"group": str(k), "decision": str(v)}
                for k, v in dict(decisions).items()
                if str(v).split(":")[0].strip().lower() not in valid
            ]
            if problems:
                (failed if sev == "critical" else warned).append(cid)
                if sev == "critical":
                    critical_failures.append(cid)
                evidence[cid] = {"invalid_decisions": problems[:20],
                                 "allowed": sorted(valid)}
            else:
                passed.append(cid)
                evidence[cid] = {"n_groups": len(dict(decisions))}

        elif cid == "structural_waters_recorded":
            # Held-out finding H02 (2026-08-08): this check was declared in the
            # contract but had no branch here at all -- it never ran, so a water
            # silently stripped after the retention decision was recorded went
            # undetected. The fix compares the actual retained-water count in the
            # prepared receptor against the recorded decision, the same
            # count-reconciliation pattern used by embedding_count above.
            decision = inputs.get("structural_water_decision")
            receptor_path_str = (inputs.get("prepared_receptor_4eiy_pdbqt")
                                  or inputs.get("prepared_receptor_3eml_pdbqt"))
            if not decision:
                (failed if sev == "critical" else warned).append(cid)
                if sev == "critical":
                    critical_failures.append(cid)
                evidence[cid] = {"warning" if sev != "critical" else "error":
                                  "no structural-water retention decision recorded"}
                continue
            required = ("cutoff", "reference_ligand", "n_waters_retained")
            missing_fields = [f for f in required if f not in decision]
            if missing_fields:
                fail(cid, sev, {"error": "water decision record incomplete",
                                 "missing_fields": missing_fields})
                continue
            receptor_path = Path(str(receptor_path_str)) if receptor_path_str else None
            if not receptor_path or not receptor_path.is_file():
                fail(cid, sev, {"error": f"prepared receptor not found: {receptor_path}"})
                continue
            actual_n = _count_structural_waters(receptor_path)
            recorded_n = int(decision["n_waters_retained"])
            detail = {"recorded_n_waters_retained": recorded_n, "actual_n_waters_in_receptor": actual_n,
                       "cutoff": decision["cutoff"], "reference_ligand": decision["reference_ligand"]}
            if actual_n != recorded_n:
                detail["error"] = ("prepared receptor's actual retained-water count does not "
                                    "match the recorded decision -- waters were added or removed "
                                    "after the decision was recorded")
                fail(cid, sev, detail)
            else:
                passed.append(cid)
                evidence[cid] = detail

        elif cid == "water_shell_selection_criterion_registered":
            # Choosing a water shell because it admits more of the library
            # risks selecting the validation that yields the desired answer
            # (ADR-AIIM-073). Only applies once more than one shell was
            # evaluated; a single frozen shell has nothing to select among.
            candidates = inputs.get("water_shell_candidates") or []
            criterion = inputs.get("water_shell_selection_criterion")
            registered_before = inputs.get("water_shell_criterion_registered_before_selection")
            if len(candidates) <= 1:
                passed.append(cid)
                evidence[cid] = {"n_candidates": len(candidates),
                                 "note": "not applicable -- single water shell used"}
                continue
            allowed_bases = set(check.get("allowed_criterion_bases", []))
            criterion_basis = str((criterion or {}).get("basis", "")).strip()
            detail = {"n_candidates": len(candidates), "candidates": candidates,
                       "criterion": criterion, "registered_before_selection": registered_before}
            if not criterion:
                fail(cid, sev, {**detail, "error": "no selection criterion recorded"})
                continue
            if allowed_bases and criterion_basis not in allowed_bases:
                fail(cid, sev, {**detail,
                                 "error": f"selection criterion basis {criterion_basis!r} is not "
                                          f"an external structure-based recovery criterion -- "
                                          f"library admission rate must not drive the choice"})
                continue
            if registered_before is not True:
                fail(cid, sev, {**detail,
                                 "error": "selection criterion was not confirmed registered "
                                          "before the admission-rate outcome was observed"})
                continue
            passed.append(cid)
            evidence[cid] = detail

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

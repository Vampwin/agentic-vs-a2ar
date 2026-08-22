"""
G0 evaluator — source / provenance checks.

Verifies SHA-256 hashes of all frozen inputs, compound count, ID uniqueness,
SMILES column presence, and absence of assay-derived column headers.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import Chem


from ._assay_blindness import activity_annotations_in, blocked_columns_in


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest().upper()


def _verify_hash(
    cid: str,
    sev: str,
    path: Path,
    check: dict[str, Any],
    inputs: dict[str, Any],
    passed: list[str],
    failed: list[str],
    warned: list[str],
    critical_failures: list[str],
    evidence: dict[str, Any],
) -> None:
    """
    Compare a file's SHA-256 against its registered anchor.

    The anchor may come from the contract (`expected_sha256`) or from the
    caller's registered-hash map. If no anchor exists the check is *unverified*,
    which for a critical check is a failure: reporting PASS would mean an
    unregistered input silently satisfies provenance, and a corrupted file
    would flow downstream unnoticed. Previously a missing anchor raised
    KeyError, which the harness recorded as an evaluator crash rather than a
    provenance finding.
    """
    anchors = inputs.get("registered_input_hashes") or {}
    expected = check.get("expected_sha256") or anchors.get(cid) or anchors.get(str(path))
    actual = _sha256(path)

    if not expected:
        detail = {"actual": actual, "expected": None,
                  "error": "no registered SHA-256 anchor; hash could not be verified"}
        if sev == "critical":
            failed.append(cid)
            critical_failures.append(cid)
        else:
            warned.append(cid)
        evidence[cid] = detail
        return

    expected = str(expected).upper()
    ok = actual == expected
    evidence[cid] = {"expected": expected, "actual": actual, "match": ok}
    (passed if ok else failed).append(cid)
    if not ok and sev == "critical":
        critical_failures.append(cid)


def evaluate(
    contract: dict[str, Any],
    inputs: dict[str, Any],
    registry: Any,
) -> dict[str, Any]:
    checks = contract.get("checks", [])
    passed, failed, warned, evidence = [], [], [], {}
    critical_failures = []

    snp_csv = Path(inputs.get("snp_library_csv", ""))
    snp_sdf = Path(inputs.get("snp_library_sdf", ""))
    snp_xlsx = Path(inputs.get("snp_library_raw_xlsx", ""))
    pdb_4eiy = Path(inputs.get("receptor_pdb_4eiy", ""))
    pdb_3eml = Path(inputs.get("receptor_pdb_3eml", ""))
    chembl_csv = Path(inputs.get("chembl_benchmark_csv", ""))
    chembl_raw_csv = Path(inputs.get("chembl_a2ar_raw_csv", ""))

    # ── File existence shortcuts ──────────────────────────────────────────────
    def file_exists(p: Path, check_id: str, sev: str) -> bool:
        if not p or not p.is_file():
            (failed if sev == "critical" else warned).append(check_id)
            if sev == "critical":
                critical_failures.append(check_id)
            evidence[check_id] = {"error": f"File not found: {p}"}
            return False
        return True

    # ── Defensive CSV access ──────────────────────────────────────────────────
    _csv_cache: dict[Path, Any] = {}

    def read_csv(p: Path, check_id: str, sev: str, **kwargs: Any):
        """
        Parse a CSV, recording a parse failure as a finding rather than raising.

        A corrupted input is exactly what this gate exists to detect, so an
        undecodable or malformed file must surface as a failed provenance check.
        Letting the exception propagate would abort the whole evaluator, and the
        harness scores an evaluator crash as an escape — turning a caught fault
        into a missed one.
        """
        key = (p, tuple(sorted(kwargs.items())))
        if key in _csv_cache:
            cached = _csv_cache[key]
        else:
            try:
                cached = pd.read_csv(p, **kwargs)
            except Exception as exc:
                cached = f"{type(exc).__name__}: {exc}"
            _csv_cache[key] = cached
        if isinstance(cached, str):
            (failed if sev == "critical" else warned).append(check_id)
            if sev == "critical":
                critical_failures.append(check_id)
            evidence[check_id] = {"error": f"Could not parse {p.name}: {cached}"}
            return None
        return cached

    # ── Hash checks ───────────────────────────────────────────────────────────
    for check in checks:
        cid = check["id"]
        sev = check.get("severity", "medium")

        if cid == "hash_snp_library_csv" and file_exists(snp_csv, cid, sev):
            _verify_hash(cid, sev, snp_csv, check, inputs,
                         passed, failed, warned, critical_failures, evidence)

        elif cid == "hash_snp_library_sdf" and file_exists(snp_sdf, cid, sev):
            _verify_hash(cid, sev, snp_sdf, check, inputs,
                         passed, failed, warned, critical_failures, evidence)

        elif cid == "hash_raw_workbook" and file_exists(snp_xlsx, cid, sev):
            _verify_hash(cid, sev, snp_xlsx, check, inputs,
                         passed, failed, warned, critical_failures, evidence)

        elif cid == "compound_count" and file_exists(snp_csv, cid, sev):
            df = read_csv(snp_csv, cid, sev)
            if df is None:
                continue
            count = len(df)
            expected = check["expected_count"]
            ok = count == expected
            evidence[cid] = {"expected": expected, "actual": count}
            (passed if ok else failed).append(cid)
            if not ok and sev == "critical":
                critical_failures.append(cid)

        elif cid == "id_uniqueness" and file_exists(snp_csv, cid, sev):
            df = read_csv(snp_csv, cid, sev)
            if df is None:
                continue
            id_col = _find_id_column(df)
            if id_col is None:
                failed.append(cid)
                critical_failures.append(cid)
                evidence[cid] = {"error": "No SNP_ID or ID column found"}
            else:
                n_dupes = int(df[id_col].duplicated().sum())
                ok = n_dupes == 0
                evidence[cid] = {"id_column": id_col, "duplicate_count": n_dupes}
                (passed if ok else failed).append(cid)
                if not ok and sev == "critical":
                    critical_failures.append(cid)

        elif cid == "smiles_column_present" and file_exists(snp_csv, cid, sev):
            df = read_csv(snp_csv, cid, sev)
            if df is None:
                continue
            smiles_col = _find_smiles_column(df)
            if smiles_col is None:
                failed.append(cid)
                critical_failures.append(cid)
                evidence[cid] = {"error": "No SMILES column found"}
            else:
                null_count = int(df[smiles_col].isna().sum())
                ok = null_count == 0
                evidence[cid] = {"smiles_column": smiles_col, "null_count": null_count}
                (passed if ok else failed).append(cid)
                if not ok and sev == "critical":
                    critical_failures.append(cid)

        elif cid == "receptor_files_present":
            ok = pdb_4eiy.exists() and pdb_3eml.exists()
            evidence[cid] = {"4eiy": pdb_4eiy.exists(), "3eml": pdb_3eml.exists()}
            (passed if ok else failed).append(cid)
            if not ok and sev == "critical":
                critical_failures.append(cid)

        elif cid == "receptor_hash_4eiy" and file_exists(pdb_4eiy, cid, sev):
            expected = check.get("expected_sha256", "").upper()
            actual = _sha256(pdb_4eiy)
            ok = not expected or actual == expected
            evidence[cid] = {"expected": expected or "not_set", "actual": actual, "match": ok}
            (passed if ok else warned).append(cid)

        elif cid == "receptor_hash_3eml" and file_exists(pdb_3eml, cid, sev):
            expected = check.get("expected_sha256", "").upper()
            actual = _sha256(pdb_3eml)
            ok = not expected or actual == expected
            evidence[cid] = {"expected": expected or "not_set", "actual": actual, "match": ok}
            (passed if ok else warned).append(cid)

        elif cid == "chembl_benchmark_present":
            ok = bool(chembl_csv) and chembl_csv.is_file()
            evidence[cid] = {"exists": ok, "path": str(chembl_csv)}
            (passed if ok else failed).append(cid)
            if not ok and sev == "critical":
                critical_failures.append(cid)

        elif cid == "chembl_benchmark_hash" and file_exists(chembl_csv, cid, sev):
            _verify_hash(cid, sev, chembl_csv, check, inputs,
                         passed, failed, warned, critical_failures, evidence)

        elif cid == "chembl_a2ar_raw_present":
            ok = bool(chembl_raw_csv) and chembl_raw_csv.is_file()
            evidence[cid] = {"exists": ok, "path": str(chembl_raw_csv)}
            (passed if ok else failed).append(cid)
            if not ok and sev == "critical":
                critical_failures.append(cid)

        elif cid == "hash_chembl_a2ar_raw_csv" and file_exists(chembl_raw_csv, cid, sev):
            _verify_hash(cid, sev, chembl_raw_csv, check, inputs,
                         passed, failed, warned, critical_failures, evidence)

        elif cid == "no_assay_columns":
            # Scope: the source library only. The public ChEMBL benchmark is
            # deliberately excluded — it is a bioactivity dataset and its labels
            # are supposed to exist. Isolating them is G5a's job (they must not
            # be *read* before the validation gate), not a matter of pretending
            # the columns are absent.
            header = read_csv(snp_csv, cid, sev, nrows=0) if snp_csv.is_file() else None
            if header is None:
                # read_csv already recorded the parse failure; "unreadable" must
                # never be silently treated as "clean".
                continue
            blocked_found = blocked_columns_in(header.columns)
            literature = activity_annotations_in(header.columns)
            ok = len(blocked_found) == 0
            evidence[cid] = {
                "scanned": snp_csv.name,
                "blocked_columns_found": blocked_found,
                "literature_activity_annotations": literature,
                "note": (
                    "Literature activity annotations are permitted in the curated "
                    "source library but must not survive standardization; G1 "
                    "enforces their absence from the ledger."
                ) if literature else "",
                "benchmark_excluded_from_scan": (
                    str(chembl_csv.name) if chembl_csv else None
                ),
            }
            (passed if ok else failed).append(cid)
            if not ok and sev == "critical":
                critical_failures.append(cid)

        elif cid == "name_structure_mapping":
            # Cross-checks the CSV's structure column against the independently
            # generated SDF, matched by compound ID -- not just column presence.
            # A row-level identity/structure misalignment (e.g. a positional
            # shift introduced by an independent sort or merge) leaves the
            # "name" column looking perfectly normal; only comparing the two
            # separately-sourced structures for the same ID can catch it.
            # Held-out finding H01 (2026-08-08): the previous version of this
            # check verified only that a "name" column existed and never
            # compared structures at all, so a full-library row shift passed
            # silently. Falls back to the column-presence check only when no
            # SDF is supplied, since the cross-check needs a second, structure-
            # bearing artifact independent of the CSV to compare against.
            if not snp_csv.exists():
                warned.append(cid)
                continue
            df = read_csv(snp_csv, cid, sev, nrows=None)
            if df is None:
                continue
            has_name = any("name" in c.lower() for c in df.columns)
            id_col = _find_id_column(df)
            smiles_col = _find_smiles_column(df)
            if smiles_col is None and "structure" in df.columns:
                smiles_col = "structure"

            if not snp_sdf.exists() or id_col is None or smiles_col is None:
                evidence[cid] = {
                    "name_column_present": has_name,
                    "cross_check": "not_run",
                    "reason": "SDF, ID column, or structure column unavailable"
                              " -- falling back to column-presence check only",
                }
                (passed if has_name else warned).append(cid)
                continue

            try:
                sdf_structures = _sdf_structures_by_id(snp_sdf)
            except Exception as exc:
                evidence[cid] = {"name_column_present": has_name, "cross_check": "not_run",
                                  "reason": f"could not parse SDF: {type(exc).__name__}: {exc}"}
                (passed if has_name else warned).append(cid)
                continue

            n_checked = 0
            mismatches: list[dict[str, Any]] = []
            for _, row in df.iterrows():
                rid = str(row.get(id_col, "")).strip()
                if rid not in sdf_structures:
                    continue
                n_checked += 1
                csv_canon = _canonical(row.get(smiles_col))
                sdf_canon = sdf_structures[rid]
                if csv_canon is None or sdf_canon is None or csv_canon != sdf_canon:
                    mismatches.append({"id": rid, "csv_smiles": row.get(smiles_col),
                                        "csv_canonical": csv_canon, "sdf_canonical": sdf_canon})

            resolution_rate = (n_checked - len(mismatches)) / n_checked if n_checked else 0.0
            ok = has_name and n_checked > 0 and resolution_rate >= 0.98
            evidence[cid] = {
                "name_column_present": has_name,
                "cross_check": "ran",
                "n_ids_matched_between_csv_and_sdf": n_checked,
                "n_structure_mismatches": len(mismatches),
                "resolution_rate": round(resolution_rate, 4),
                "threshold": 0.98,
                "mismatch_sample": mismatches[:20],
            }
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


def _find_id_column(df: pd.DataFrame) -> str | None:
    for candidate in ["SNP_ID", "snp_id", "Compound", "compound", "ID", "id", "compound_id"]:
        if candidate in df.columns:
            return candidate
    return None


def _find_smiles_column(df: pd.DataFrame) -> str | None:
    for candidate in ["SMILES", "smiles", "Smiles", "canonical_smiles"]:
        if candidate in df.columns:
            return candidate
    return None


def _canonical(smiles: str) -> str | None:
    mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    return Chem.MolToSmiles(mol, canonical=True) if mol is not None else None


def _sdf_structures_by_id(sdf_path: Path) -> dict[str, str | None]:
    """id (mol block title / _Name) -> canonical SMILES, or None if unparseable."""
    out: dict[str, str | None] = {}
    for mol in Chem.SDMolSupplier(str(sdf_path), removeHs=False):
        if mol is None:
            continue
        name = mol.GetProp("_Name").strip() if mol.HasProp("_Name") else None
        if not name:
            continue
        try:
            out[name] = Chem.MolToSmiles(mol, canonical=True)
        except Exception:
            out[name] = None
    return out

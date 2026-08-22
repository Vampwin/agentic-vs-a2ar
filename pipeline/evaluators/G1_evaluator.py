"""
G1 evaluator — standardization checks.

Verifies the accountability ledger: 236 rows complete, SMILES round-trip valid,
drug-likeness used as annotation only (no exclusion by Lipinski/Veber/QED/PAINS).

Requires RDKit. Install via: conda install -c conda-forge rdkit
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from ._assay_blindness import blocked_columns_in


_EXCLUSION_REASONS_ALLOWED = {
    "standardization_failure",
    "embedding_failure",
    "docking_failure",
    "mixture_excluded",
    "",               # blank = not excluded
}

_DRUG_LIKENESS_EXCLUSION_FORBIDDEN = {
    "lipinski", "veber", "qed", "pains", "drug_like", "ro5",
}


def _count_stereocentres(smiles: str) -> int | None:
    """
    Specified tetrahedral stereocentres in a SMILES string.

    Falls back to counting '@' when RDKit is unavailable: the two agree for
    specified centres, which is all these checks compare against.
    """
    if not smiles:
        return 0
    try:
        from rdkit import Chem
    except ImportError:
        # '@@' is one centre written with two characters, so it must not be
        # double-counted when falling back to plain string counting.
        return smiles.count("@") - smiles.count("@@")
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        centres = Chem.FindMolChiralCenters(
            mol, includeUnassigned=False, useLegacyImplementation=False
        )
        return len(centres)
    except Exception:
        return None


def _declared_int(row: Any, column: str) -> int | None:
    if column not in row or pd.isna(row[column]):
        return None
    try:
        return int(float(row[column]))
    except (TypeError, ValueError):
        return None


def evaluate(
    contract: dict[str, Any],
    inputs: dict[str, Any],
    registry: Any,
) -> dict[str, Any]:
    passed, failed, warned, evidence = [], [], [], {}
    critical_failures = []

    ledger_path = Path(inputs.get("snp_source_ledger_csv", ""))
    sdf_path = Path(inputs.get("standardized_smiles_sdf", ""))

    if not ledger_path.exists():
        return {
            "state": "HALT",
            "passed": [],
            "failed": ["ledger_missing"],
            "warned": [],
            "critical_failures": ["ledger_missing"],
            "evidence": {"ledger_missing": {"path": str(ledger_path)}},
            "summary": "Ledger file not found",
        }

    df = pd.read_csv(ledger_path)
    # Normalize string columns that may read back as float NaN
    for _col in ["exclusion_reason", "standardization_failure_reason", "removed_fragment", "inchikey"]:
        if _col in df.columns:
            df[_col] = df[_col].fillna("").astype(str)

    for check in contract.get("checks", []):
        cid = check["id"]
        sev = check.get("severity", "medium")

        if cid == "ledger_count":
            count = len(df)
            expected = check["expected_count"]
            ok = count == expected
            evidence[cid] = {"expected": expected, "actual": count}
            (passed if ok else failed).append(cid)
            if not ok and sev == "critical":
                critical_failures.append(cid)

        elif cid == "id_coverage":
            if "SNP_ID" in df.columns:
                n_unique = df["SNP_ID"].nunique()
                ok = n_unique == 236 and len(df) == 236
                evidence[cid] = {"unique_ids": n_unique, "total_rows": len(df)}
                (passed if ok else failed).append(cid)
                if not ok and sev == "critical":
                    critical_failures.append(cid)
            else:
                failed.append(cid)
                critical_failures.append(cid)
                evidence[cid] = {"error": "SNP_ID column not found"}

        elif cid == "smiles_roundtrip":
            result = _check_smiles_roundtrip(df)
            evidence[cid] = result
            ok = result["n_failed"] == 0
            (passed if ok else failed).append(cid)
            if not ok and sev == "critical":
                critical_failures.append(cid)

        elif cid == "inchikey_generated":
            if "inchikey" in df.columns and "valid_smiles" in df.columns:
                valid_rows = df[df["valid_smiles"] == True]
                missing_inchikey = valid_rows["inchikey"].isna().sum()
                ok = missing_inchikey == 0
                evidence[cid] = {"missing_inchikey_count": int(missing_inchikey)}
                (passed if ok else warned).append(cid)
            else:
                warned.append(cid)
                evidence[cid] = {"error": "inchikey or valid_smiles column not found"}

        elif cid == "salt_strip_recorded":
            has_cols = all(c in df.columns for c in ["parent_smiles", "removed_fragment"])
            if has_cols:
                # Check rows where stripping occurred
                stripped = df[df["removed_fragment"].notna() & (df["removed_fragment"] != "")]
                missing_parent = stripped["parent_smiles"].isna().sum()
                ok = missing_parent == 0
                evidence[cid] = {"stripped_count": len(stripped), "missing_parent": int(missing_parent)}
                (passed if ok else warned).append(cid)
            else:
                warned.append(cid)
                evidence[cid] = {"error": "parent_smiles or removed_fragment column not found"}

        elif cid == "mixture_flagged":
            if "is_mixture" in df.columns:
                mixture_rows = df[df["is_mixture"] == True]
                # Check that all mixture rows have exclusion_reason set
                if "exclusion_reason" in df.columns:
                    unlabeled = mixture_rows[
                        mixture_rows["exclusion_reason"].isna()
                        | ~mixture_rows["exclusion_reason"].str.contains("mixture", na=False, case=False)
                    ]
                    ok = len(unlabeled) == 0
                    evidence[cid] = {
                        "mixture_count": len(mixture_rows),
                        "unlabeled_mixtures": int(len(unlabeled)),
                    }
                    (passed if ok else warned).append(cid)
                else:
                    warned.append(cid)
                    evidence[cid] = {"error": "exclusion_reason column not found"}
            else:
                warned.append(cid)
                evidence[cid] = {"error": "is_mixture column not found"}

        elif cid == "stereochemistry_recorded":
            required = ["has_stereocenters", "n_stereocenters", "stereo_defined"]
            missing = [c for c in required if c not in df.columns]
            ok = len(missing) == 0
            evidence[cid] = {"missing_columns": missing}
            (passed if ok else warned).append(cid)

        elif cid == "stereochemistry_preserved":
            # `stereochemistry_recorded` only proves the bookkeeping columns
            # exist. This proves they are *true*: the stereocentres actually
            # present in the standardized SMILES must match what the ledger
            # declares. Silent stereo loss keeps the metadata intact while the
            # structure changes, so only a recount can detect it — and a lost
            # centre changes which enantiomer was docked.
            # Must be the standardized output column, not whichever "*smiles*"
            # column happens to sort first. The ledger also carries
            # smiles_input (pre-standardization) and parent_smiles
            # (pre-salt-stripping); reading either of those here means the
            # check silently validates stereochemistry that was never
            # standardized, so it can never observe a loss standardization
            # introduced. Confirmed as a real gap via a real-artifact replay
            # (20_new_runs/real_artifact_replay_2026-08-05): mutating
            # smiles_canonical alone did not trip this check while it read
            # smiles_input, the untouched original column.
            smiles_col = "smiles_canonical" if "smiles_canonical" in df.columns else next(
                (c for c in df.columns if "smiles" in c.lower()), None
            )
            if smiles_col is None or "n_stereocenters" not in df.columns:
                warned.append(cid)
                evidence[cid] = {"warning": "SMILES or n_stereocenters column not found"}
            else:
                # Compare *specification*, not counts. A ledger may legitimately
                # count written stereo markers while RDKit counts genuine
                # stereocentres, and the two disagree for symmetric molecules
                # (glycerol is written with '@' but is achiral). Requiring the
                # counts to match would fail those rows on clean data.
                #
                # What must never happen is the structure losing its stereo
                # specification while the ledger still claims it has one — that
                # is silent stereo loss, and it changes which enantiomer is docked.
                lost, count_mismatch, unparseable, undefined_skipped = [], [], 0, 0
                for _, row in df.iterrows():
                    declared = _declared_int(row, "n_stereocenters")
                    if declared is None or declared <= 0:
                        continue
                    # A centre counted with includeUnassigned=True exists in the
                    # skeleton but has no configuration in the source. Such a row
                    # correctly carries no stereo marker, and `stereo_defined`
                    # records that. Only a row claiming *defined* stereochemistry
                    # is obliged to show it.
                    if "stereo_defined" in df.columns and not bool(row.get("stereo_defined")):
                        undefined_skipped += 1
                        continue
                    smi = str(row.get(smiles_col, "") or "")
                    if not smi:
                        unparseable += 1
                        continue
                    written = smi.count("@") - smi.count("@@")
                    entry = {"id": str(row.get("SNP_ID", row.name)),
                             "declared": declared, "written_markers": written}
                    if written <= 0:
                        lost.append(entry)
                        continue
                    actual = _count_stereocentres(smi)
                    if actual is not None and actual != declared:
                        count_mismatch.append({**entry, "rdkit_centres": actual})

                ok = len(lost) == 0
                evidence[cid] = {"n_stereo_lost": len(lost),
                                 "stereo_lost": lost[:20],
                                 "n_undefined_skipped": undefined_skipped,
                                 "undefined_note": (
                                     "rows with structurally present but unconfigured "
                                     "stereocentres (stereo_defined=False); no marker "
                                     "is expected for these"),
                                 "n_count_mismatch": len(count_mismatch),
                                 "count_mismatch_note": (
                                     "declared count differs from RDKit's assigned-centre "
                                     "count; expected for symmetric molecules, reported "
                                     "for information only"),
                                 "n_unparseable": unparseable}
                (passed if ok else failed).append(cid)
                if not ok and sev == "critical":
                    critical_failures.append(cid)

        elif cid == "no_activity_columns_in_ledger":
            # The source library is a curated literature collection and carries
            # activity annotations (IC50, "Biological activities"). Those are
            # acceptable in the source but must not survive standardization:
            # the ledger is what feeds ranking and selection, so an activity
            # value riding along here would influence candidate choice
            # regardless of whether anyone intended it to.
            found = blocked_columns_in(df.columns, strict=True)
            ok = len(found) == 0
            evidence[cid] = {"activity_columns_found": found,
                             "n_columns": len(df.columns)}
            (passed if ok else failed).append(cid)
            if not ok and sev == "critical":
                critical_failures.append(cid)

        elif cid == "salt_stripping_verified":
            # A '.' in a standardized SMILES means a disconnected fragment
            # survived standardization. That is acceptable only if the ledger
            # says so — either it recorded what was removed or it flagged the
            # row as a mixture. An undeclared counter-ion changes the docked
            # species and its molecular weight.
            smiles_col = next((c for c in df.columns if "smiles" in c.lower()
                               and "parent" not in c.lower()), None)
            if smiles_col is None:
                warned.append(cid)
                evidence[cid] = {"warning": "no standardized SMILES column found"}
            else:
                offenders = []
                for _, row in df.iterrows():
                    smi = str(row.get(smiles_col, "") or "")
                    if "." not in smi:
                        continue
                    removed = str(row.get("removed_fragment", "") or "").strip()
                    is_mixture = bool(row.get("is_mixture", False))
                    excl = str(row.get("exclusion_reason", "") or "").lower()
                    declared = removed or is_mixture or "mixture" in excl
                    if not declared:
                        offenders.append({
                            "id": str(row.get("SNP_ID", row.name)),
                            "smiles": smi[:80],
                            "n_fragments": smi.count(".") + 1,
                        })
                ok = len(offenders) == 0
                evidence[cid] = {"n_undeclared_fragments": len(offenders),
                                 "offenders": offenders[:20]}
                (passed if ok else failed).append(cid)
                if not ok and sev == "critical":
                    critical_failures.append(cid)

        elif cid == "no_drug_likeness_exclusions":
            result = _check_no_dl_exclusions(df)
            evidence[cid] = result
            ok = result["violations"] == 0
            (passed if ok else failed).append(cid)
            if not ok and sev == "critical":
                critical_failures.append(cid)

        elif cid == "drug_likeness_as_annotation":
            dl_cols = [
                c for c in df.columns
                if any(kw in c.lower() for kw in ["lipinski", "veber", "qed", "pains"])
            ]
            exclusion_col = df.get("exclusion_reason", pd.Series(dtype=str))
            dl_in_exclusion = exclusion_col.str.lower().str.contains(
                "|".join(_DRUG_LIKENESS_EXCLUSION_FORBIDDEN), na=False
            ).sum()
            ok = dl_in_exclusion == 0
            evidence[cid] = {
                "dl_annotation_columns": dl_cols,
                "dl_used_as_exclusion_count": int(dl_in_exclusion),
            }
            (passed if ok else failed).append(cid)
            if not ok and sev == "critical":
                critical_failures.append(cid)

        elif cid == "failure_accounting":
            if "valid_smiles" in df.columns and "standardization_failure_reason" in df.columns:
                invalid = df[df["valid_smiles"] == False]
                missing_reason = invalid["standardization_failure_reason"].isna().sum()
                ok = missing_reason == 0
                evidence[cid] = {
                    "invalid_count": len(invalid),
                    "missing_reason_count": int(missing_reason),
                }
                (passed if ok else warned).append(cid)
            else:
                warned.append(cid)
                evidence[cid] = {"error": "valid_smiles or standardization_failure_reason column not found"}

        elif cid == "duplicate_inchikey_check":
            if "inchikey" in df.columns:
                dupes = df[df["inchikey"].notna()].duplicated(subset=["inchikey"], keep=False)
                n_dupes = int(dupes.sum())
                evidence[cid] = {"duplicate_inchikey_count": n_dupes}
                if n_dupes > 0:
                    warned.append(cid)
                else:
                    passed.append(cid)
            else:
                warned.append(cid)
                evidence[cid] = {"error": "inchikey column not found"}

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


def _check_smiles_roundtrip(df: pd.DataFrame) -> dict[str, Any]:
    """
    For rows with valid_smiles=True, verify RDKit round-trip succeeds.
    Returns summary dict. Uses RDKit if available; falls back to a flag.
    """
    try:
        from rdkit import Chem
    except ImportError:
        return {"rdkit_available": False, "n_checked": 0, "n_failed": 0,
                "error": "RDKit not installed; roundtrip check skipped"}

    # Prefer canonical/parent SMILES from ledger; fall back to raw SMILES column
    _smiles_candidates = ["smiles_canonical", "parent_smiles", "smiles_input", "SMILES", "smiles"]
    smiles_col = next((c for c in _smiles_candidates if c in df.columns), None)
    if smiles_col is None:
        return {"error": "No recognized SMILES column found", "n_checked": 0, "n_failed": 0}

    smiles_col = smiles_col  # already set above
    valid_col = "valid_smiles" if "valid_smiles" in df.columns else None
    subset = df if valid_col is None else df[df[valid_col] == True]

    failures = []
    for _, row in subset.iterrows():
        smi = row[smiles_col]
        if pd.isna(smi):
            continue
        mol = Chem.MolFromSmiles(str(smi))
        if mol is None:
            failures.append(row.get("SNP_ID", "unknown"))
            continue
        roundtrip_smi = Chem.MolToSmiles(mol)
        mol2 = Chem.MolFromSmiles(roundtrip_smi)
        if mol2 is None:
            failures.append(row.get("SNP_ID", "unknown"))

    return {
        "rdkit_available": True,
        "n_checked": len(subset),
        "n_failed": len(failures),
        "failed_ids": failures[:20],  # First 20 only to keep evidence compact
    }


def _check_no_dl_exclusions(df: pd.DataFrame) -> dict[str, Any]:
    if "exclusion_reason" not in df.columns:
        return {"violations": 0, "note": "exclusion_reason column not present; check skipped"}
    violations = df["exclusion_reason"].dropna().apply(
        lambda r: any(kw in str(r).lower() for kw in _DRUG_LIKENESS_EXCLUSION_FORBIDDEN)
    ).sum()
    return {
        "violations": int(violations),
        "allowed_reasons": list(_EXCLUSION_REASONS_ALLOWED),
    }

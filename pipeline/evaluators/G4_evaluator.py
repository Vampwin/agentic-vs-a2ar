"""
G4 evaluator — docking / scoring verification.

Implements contracts/G4_scoring.yaml: output completeness, score direction,
seed repeatability, missing-output accounting, score-distribution sanity,
molecular-weight bias, and assay-column exclusion.

Severity follows the contract, not convenience. score_direction,
output_completeness, seed_repeatability and no_assay_columns_in_output are
critical and halt the workflow. score_distribution_sanity and mw_bias_check are
declared medium and warn only — a size-correlated ranking can be legitimate, so
the contract deliberately surfaces it for human inspection rather than blocking.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from ._assay_blindness import scan_csv_header, scan_dataframe

_SCORE_CANDIDATES = ("vina_score", "docking_score", "vina_minimize", "vina_box_best", "score")
_ID_CANDIDATES = ("snp_id", "compound_id", "Compound", "compound", "ID", "id")
_SMILES_CANDIDATES = ("SMILES_std", "SMILES", "smiles", "canonical_smiles")
_MW_CANDIDATES = ("MW", "mw", "molecular_weight", "ExactMolWt")


def _find_column(df: pd.DataFrame, candidates: tuple[str, ...]) -> str | None:
    """Exact match first, then case-insensitive, then substring."""
    for cand in candidates:
        if cand in df.columns:
            return cand
    lowered = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lowered:
            return lowered[cand.lower()]
    for cand in candidates:
        for col in df.columns:
            if cand.lower() in col.lower():
                return col
    return None


def _molecular_weights(df: pd.DataFrame) -> tuple[np.ndarray | None, str]:
    """MW per row, from an MW column if present, else computed from SMILES."""
    mw_col = _find_column(df, _MW_CANDIDATES)
    if mw_col is not None:
        return pd.to_numeric(df[mw_col], errors="coerce").to_numpy(), f"column:{mw_col}"

    smiles_col = _find_column(df, _SMILES_CANDIDATES)
    if smiles_col is None:
        return None, "unavailable"
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
    except ImportError:
        return None, "rdkit_unavailable"

    values = []
    for smi in df[smiles_col]:
        mol = Chem.MolFromSmiles(str(smi)) if pd.notna(smi) else None
        values.append(Descriptors.ExactMolWt(mol) if mol is not None else np.nan)
    return np.asarray(values, dtype=float), f"computed_from:{smiles_col}"


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

    score_csv = Path(str(inputs.get("vina_score_csv", "")))
    df: pd.DataFrame | None = None
    read_error: str | None = None
    if score_csv and score_csv.is_file():
        try:
            df = pd.read_csv(score_csv)
        except Exception as exc:  # unparseable output is itself a finding
            read_error = f"{type(exc).__name__}: {exc}"
    else:
        read_error = f"File not found: {score_csv}"

    score_col = _find_column(df, _SCORE_CANDIDATES) if df is not None else None
    scores = (
        pd.to_numeric(df[score_col], errors="coerce")
        if df is not None and score_col is not None
        else None
    )

    for check in checks:
        cid = check["id"]
        sev = check.get("severity", "medium")

        # A critical check that cannot be evaluated must never report PASS.
        if df is None:
            fail(cid, sev, {"error": read_error or "score CSV unavailable"})
            continue

        if cid == "output_completeness":
            expected = check.get("expected_count", inputs.get("count_for_docking"))
            actual = len(df)
            if expected is None:
                warned.append(cid)
                evidence[cid] = {"actual": actual, "warning": "no expected count supplied"}
            elif actual == int(expected):
                passed.append(cid)
                evidence[cid] = {"expected": int(expected), "actual": actual}
            else:
                fail(cid, sev, {"expected": int(expected), "actual": actual})

        elif cid == "score_direction":
            if score_col is None:
                fail(cid, sev, {"error": "no score column found",
                                "columns": list(df.columns)})
                continue
            valid = scores.dropna()
            positive = valid[valid > 0]
            if len(valid) == 0:
                fail(cid, sev, {"error": "no numeric scores present"})
            elif len(positive) == 0:
                passed.append(cid)
                evidence[cid] = {"score_column": score_col, "n_scored": int(len(valid)),
                                 "n_positive": 0, "max_score": float(valid.max())}
            else:
                # Vina scores are binding free energies: more negative is better.
                # Any positive value means the sign convention was inverted or
                # the pose is non-physical; either way the ranking is unsafe.
                fail(cid, sev, {
                    "score_column": score_col,
                    "n_scored": int(len(valid)),
                    "n_positive": int(len(positive)),
                    "max_score": float(valid.max()),
                    "reason": "positive Vina scores indicate inverted sign convention",
                })

        elif cid == "seed_repeatability":
            tol = float(check.get("tolerance_kcal_mol", 0.01))
            rerun = inputs.get("seed_rerun_scores")
            if rerun is None:
                # Repeatability is critical; without a rerun it is unverified,
                # not satisfied.
                if sev == "critical":
                    fail(cid, sev, {"error": "no seed_rerun_scores supplied; "
                                             "repeatability unverified",
                                    "tolerance_kcal_mol": tol})
                else:
                    warned.append(cid)
                    evidence[cid] = {"warning": "no seed_rerun_scores supplied"}
                continue
            if score_col is None:
                fail(cid, sev, {"error": "no score column to compare against"})
                continue
            id_col = _find_column(df, _ID_CANDIDATES)
            original = (
                dict(zip(df[id_col], scores)) if id_col is not None
                else dict(enumerate(scores))
            )
            deviations = {
                k: abs(float(v) - float(original[k]))
                for k, v in dict(rerun).items()
                if k in original and pd.notna(original[k]) and pd.notna(v)
            }
            worst = max(deviations.values()) if deviations else 0.0
            offenders = {k: round(d, 6) for k, d in deviations.items() if d > tol}
            if not deviations:
                fail(cid, sev, {"error": "rerun scores matched no original compound"})
            elif offenders:
                fail(cid, sev, {"tolerance_kcal_mol": tol, "n_compared": len(deviations),
                                "max_deviation": round(worst, 6), "offenders": offenders})
            else:
                passed.append(cid)
                evidence[cid] = {"tolerance_kcal_mol": tol, "n_compared": len(deviations),
                                 "max_deviation": round(worst, 6)}

        elif cid == "missing_output_accounting":
            expected_ids = inputs.get("expected_compound_ids")
            ledger = inputs.get("failure_ledger") or {}
            if expected_ids is None:
                warned.append(cid)
                evidence[cid] = {"warning": "expected_compound_ids unavailable"}
                continue
            id_col = _find_column(df, _ID_CANDIDATES)
            if id_col is None:
                fail(cid, sev, {"error": "no ID column in score CSV"})
                continue
            scored = set(df[id_col].astype(str))
            explained = {str(k) for k in ledger}
            unaccounted = [
                str(c) for c in expected_ids
                if str(c) not in scored and str(c) not in explained
            ]
            if unaccounted:
                # Every compound must terminate in exactly one state: scored or
                # explained. Silent disappearance is the fault this check exists for.
                fail(cid, sev, {"n_unaccounted": len(unaccounted),
                                "unaccounted_sample": unaccounted[:20]})
            else:
                passed.append(cid)
                evidence[cid] = {"n_expected": len(list(expected_ids)),
                                 "n_scored": len(scored), "n_explained": len(explained)}

        elif cid == "score_distribution_sanity":
            lo, hi = check.get("expected_range_kcal_mol", [-12.0, -3.0])
            if score_col is None:
                fail(cid, sev, {"error": "no score column found"})
                continue
            valid = scores.dropna()
            median = float(valid.median()) if len(valid) else float("nan")
            in_range = bool(lo <= median <= hi) if len(valid) else False
            evidence[cid] = {"median": round(median, 4) if len(valid) else None,
                             "expected_range": [lo, hi], "in_range": in_range}
            if in_range:
                passed.append(cid)
            elif sev == "critical":
                fail(cid, sev, evidence[cid])
            else:
                warned.append(cid)

        elif cid == "mw_bias_check":
            if score_col is None:
                fail(cid, sev, {"error": "no score column found"})
                continue
            mw, mw_source = _molecular_weights(df)
            if mw is None:
                warned.append(cid)
                evidence[cid] = {"warning": f"MW unavailable ({mw_source})"}
                continue
            pair = pd.DataFrame({"mw": mw, "score": scores}).dropna()
            if len(pair) < 3:
                warned.append(cid)
                evidence[cid] = {"warning": "fewer than 3 comparable rows"}
                continue
            # Spearman is undefined when either series is constant; report that
            # explicitly rather than emitting a NaN through the threshold test.
            if pair["mw"].nunique() < 2 or pair["score"].nunique() < 2:
                warned.append(cid)
                evidence[cid] = {"warning": "correlation undefined (constant MW or score)",
                                 "mw_source": mw_source, "n": int(len(pair))}
                continue
            rho = float(pair["mw"].corr(pair["score"], method="spearman"))
            limit = float(check.get("abs_correlation_max", 0.5))
            # Two-tier. Moderate size correlation is a genuine property of
            # natural-product libraries and is reported as a ranking caveat.
            # A near-perfect correlation is a different claim entirely: it means
            # the score column is a monotone function of molecular weight, i.e.
            # it is not a docking score at all, so it fails outright.
            degenerate_limit = float(check.get("abs_correlation_degenerate", 0.95))
            abs_rho = abs(rho) if not np.isnan(rho) else float("nan")
            within = abs_rho <= limit if not np.isnan(abs_rho) else False
            degenerate = abs_rho >= degenerate_limit if not np.isnan(abs_rho) else False
            evidence[cid] = {"spearman_rho": round(rho, 4) if not np.isnan(rho) else None,
                             "abs_limit": limit, "degenerate_limit": degenerate_limit,
                             "within_limit": within, "degenerate": degenerate,
                             "mw_source": mw_source, "n": int(len(pair))}
            if within:
                passed.append(cid)
            elif degenerate:
                evidence[cid]["reason"] = (
                    "score is effectively a monotone function of molecular weight; "
                    "the column is not a docking score"
                )
                failed.append(cid)
                critical_failures.append(cid)
            elif sev == "critical":
                fail(cid, sev, evidence[cid])
            else:
                warned.append(cid)

        elif cid == "no_assay_columns_in_output":
            blocked = scan_dataframe(df)
            if blocked:
                fail(cid, sev, {"blocked_columns_found": blocked})
            else:
                passed.append(cid)
                evidence[cid] = {"blocked_columns_found": [], "n_columns": len(df.columns)}

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

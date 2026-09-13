"""
conformal_qsar.py — Mondrian inductive conformal prediction wrapped around
the frozen QSAR pipeline (E6, ADR-AIIM-073).

Rationale. G6's `applicability_domain` check is a hard cutoff: a compound is
either inside or outside domain, with no notion of how confidently. E2
(ADR-AIIM-073) showed this framing cannot in principle admit a structurally
novel library under any similarity-based descriptor space. Conformal
prediction reframes the question: instead of one boolean, every compound
gets a per-class p-value and a prediction set at a chosen significance level
(0.90 confidence, 0.80 confidence, ...) that is empirically validated on
held-out data rather than merely asserted, consistent with the OECD (Q)SAR
validation principles.

This script does NOT refit the model. It wraps the already-frozen, already-
reported pipeline (`qsar_pipeline_a2a.pkl`, CV AUC 0.907, held-out test AUC
0.922) exactly as trained. The calibration set is the SAME scaffold-disjoint
held-out split already used to report that AUC (n=1217) -- re-using labeled
held-out data to calibrate nonconformity scores is standard practice for
inductive conformal prediction and does not inflate any accuracy claim; it
answers a different question (how much can a single prediction be trusted)
than AUC does (how well does the model rank, in aggregate).

Method: Mondrian (label-conditional) inductive conformal prediction.
Nonconformity score for a candidate label y at compound x:
    alpha(x, y) = 1 - p_model(y | x)
p-value:
    p(x, y) = (1 + #{i in calibration, true label y_i = y : alpha_i >= alpha(x,y)})
              / (1 + n_y)
Prediction set at significance level `sig`: {y : p(x, y) > sig}.
Label-conditional calibration keeps validity per class despite the frozen
training set's severe class imbalance (5268 actives / 702 inactives) --
each class's p-values are calibrated against its own calibration subset,
not pooled.

Validity check. Before trusting the SNP-Lib p-values, the guarantee itself
is checked empirically: 5-fold scaffold-grouped cross-validation *within*
the calibration set (never touching the training set), reporting whether
the empirical miscoverage rate tracks the nominal significance level.

Caveat that must travel with every number this script produces. Conformal
validity assumes the calibration set and the query compounds are
exchangeable (drawn from the same distribution). E2 already established
that SNP-Lib is NOT exchangeable with the ChEMBL training/calibration
distribution in structural descriptor space (0% coverage across six
independent spaces). This script's SNP-Lib p-values should therefore be
read as a more informative, per-compound alternative to a hard AD cutoff --
not as a distribution-free statistical guarantee, which formally requires
exchangeability this library does not have. The validity check on held-out
ChEMBL data confirms the *mechanism* works correctly; it cannot certify
coverage under the specific distribution shift E2 documented.

Usage:
    python conformal_qsar.py --out <dir>
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import random
from pathlib import Path
from kdss_paths import project_root, resolve, under

import joblib
import numpy as np
import pandas as pd

ROOT = project_root()
QSAR_RUN = resolve("20_new_runs/qsar_lbdd_2026-08-07")
PIPELINE_PATH = under("20_new_runs/qsar_lbdd_2026-08-07", "qsar_pipeline_a2a.pkl")
TRAINING_SET_CSV = under("20_new_runs/qsar_lbdd_2026-08-07", "chembl_a2a_training_set.csv")
SNPLIB_LEDGER = resolve("20_new_runs/snp_screen_2026-08-04/G1_standardization/snp_source_ledger.csv")

SIGNIFICANCE_LEVELS = [0.05, 0.10, 0.20]   # 95% / 90% / 80% nominal confidence
N_VALIDITY_FOLDS = 5
RANDOM_STATE = 42


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load_qsar_module():
    """
    Import run_qsar_lbdd_2026-08-07.py by path (its filename is not a valid
    module identifier) to reuse its featurization and scaffold-split
    functions byte-for-byte -- guarantees the calibration set is exactly the
    same held-out split already reported (n=1217), not a re-derived
    approximation of it.
    """
    path = Path(__file__).parent / "run_qsar_lbdd_2026-08-07.py"
    spec = importlib.util.spec_from_file_location("_qsar_lbdd_source", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def mondrian_calibrate(p_active_calib: np.ndarray, y_calib: np.ndarray) -> dict[int, np.ndarray]:
    """Per-class nonconformity score arrays, sorted ascending."""
    scores = {}
    for cls in (0, 1):
        # nonconformity of the TRUE label at each calibration point of that class
        mask = y_calib == cls
        p_true = np.where(cls == 1, p_active_calib, 1 - p_active_calib)[mask]
        scores[cls] = np.sort(1.0 - p_true)
    return scores


def mondrian_pvalues(p_active_query: np.ndarray, calib_scores: dict[int, np.ndarray]) -> dict[int, np.ndarray]:
    """p(x, y) for y in {0,1}, vectorised over all query points."""
    out = {}
    for cls in (0, 1):
        p_y = np.where(cls == 1, p_active_query, 1 - p_active_query)
        alpha_xy = 1.0 - p_y
        n_y = len(calib_scores[cls])
        # rank via searchsorted: how many calibration scores are >= alpha_xy
        n_ge = n_y - np.searchsorted(calib_scores[cls], alpha_xy, side="left")
        out[cls] = (1.0 + n_ge) / (1.0 + n_y)
    return out


def prediction_set(pvals: dict[int, np.ndarray], sig: float) -> list[str]:
    n = len(pvals[0])
    labels = []
    for i in range(n):
        s = [str(c) for c in (0, 1) if pvals[c][i] > sig]
        labels.append("{" + ",".join(s) + "}" if s else "{}")
    return labels


def validity_check(smiles_calib: list[str], y_calib: np.ndarray, p_active_calib: np.ndarray,
                   scaffold_groups: np.ndarray, n_folds: int = N_VALIDITY_FOLDS,
                   seed: int = RANDOM_STATE) -> dict:
    """
    5-fold scaffold-grouped CV within the calibration set: for each fold,
    calibrate on the other folds and check the empirical miscoverage rate
    against each nominal significance level. This never touches the
    training set -- it tests the conformal MECHANISM, using data the model
    itself never saw during training either way.
    """
    rng = np.random.default_rng(seed)
    groups = np.unique(scaffold_groups)
    rng.shuffle(groups)
    fold_of_group = {g: i % n_folds for i, g in enumerate(groups)}
    fold_id = np.array([fold_of_group[g] for g in scaffold_groups])

    results = {sig: {"errors": 0, "n": 0} for sig in SIGNIFICANCE_LEVELS}
    empty_set_count, both_set_count = 0, 0
    total = 0
    for f in range(n_folds):
        is_test = fold_id == f
        is_cal = ~is_test
        if is_test.sum() == 0 or is_cal.sum() == 0:
            continue
        cs = mondrian_calibrate(p_active_calib[is_cal], y_calib[is_cal])
        pv = mondrian_pvalues(p_active_calib[is_test], cs)
        y_true_fold = y_calib[is_test]
        for sig in SIGNIFICANCE_LEVELS:
            p_true = np.where(y_true_fold == 1, pv[1], pv[0])
            results[sig]["errors"] += int((p_true <= sig).sum())
            results[sig]["n"] += len(p_true)
        sets = prediction_set(pv, 0.20)  # diagnostic at 80% confidence
        empty_set_count += sum(1 for s in sets if s == "{}")
        both_set_count += sum(1 for s in sets if s == "{0,1}")
        total += len(sets)

    report = {}
    for sig in SIGNIFICANCE_LEVELS:
        n = results[sig]["n"]
        err_rate = results[sig]["errors"] / n if n else float("nan")
        report[f"significance_{sig}"] = {
            "nominal_error_rate": sig, "empirical_error_rate": round(err_rate, 4),
            "n": n, "valid": err_rate <= sig + 0.03,  # small slack for finite-sample noise
        }
    report["diagnostic_at_80pct_confidence"] = {
        "empty_prediction_set_fraction": round(empty_set_count / total, 4) if total else None,
        "both_labels_fraction": round(both_set_count / total, 4) if total else None,
    }
    return report


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(QSAR_RUN))
    ap.add_argument("--query-csv", default=str(SNPLIB_LEDGER),
                    help="compound library to score (default: SNP-Lib ledger)")
    ap.add_argument("--query-smiles-col", default="smiles_canonical")
    ap.add_argument("--query-id-col", default="SNP_ID")
    ap.add_argument("--out-name", default="snplib_conformal_predictions.csv",
                    help="output filename for the per-compound predictions")
    args = ap.parse_args()
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    qsar = _load_qsar_module()
    pipe = joblib.load(PIPELINE_PATH)
    print(f"loaded frozen pipeline: {PIPELINE_PATH.name} "
          f"(sha256={sha256(PIPELINE_PATH)[:16]}...) -- not refit")

    combined = pd.read_csv(TRAINING_SET_CSV)
    smiles_all = combined["canonical_smiles"].tolist()
    y_all = combined["label"].to_numpy()
    train_idx, test_idx = qsar.scaffold_train_test_split(smiles_all, y_all)
    smiles_calib = [smiles_all[i] for i in test_idx]
    y_calib = y_all[test_idx]
    print(f"reproduced scaffold split: {len(train_idx)} train (unused here) / "
          f"{len(test_idx)} calibration (identical held-out set already reported "
          f"as external test AUC 0.922)")
    print(f"calibration set class balance: {int(y_calib.sum())} active / "
          f"{int((y_calib == 0).sum())} inactive\n")

    X_calib = qsar.compute_raw_features(smiles_calib)
    p_active_calib = pipe.predict_proba(X_calib)[:, 1]
    scaffold_groups_calib = qsar.get_scaffold_groups(smiles_calib)

    print("running validity check (5-fold scaffold-grouped CV within the "
          "calibration set)...")
    validity = validity_check(smiles_calib, y_calib, p_active_calib, scaffold_groups_calib)
    for sig in SIGNIFICANCE_LEVELS:
        v = validity[f"significance_{sig}"]
        flag = "OK" if v["valid"] else "MISCALIBRATED"
        print(f"  significance={sig:.2f} (nominal {100*(1-sig):.0f}% confidence): "
              f"empirical error rate={v['empirical_error_rate']:.4f}  n={v['n']}  [{flag}]")
    diag = validity["diagnostic_at_80pct_confidence"]
    print(f"  at 80% confidence: empty sets={diag['empty_prediction_set_fraction']:.1%}, "
          f"both-labels sets={diag['both_labels_fraction']:.1%}\n")

    # ---- final calibration on the FULL held-out set, applied to the query library ----
    calib_scores = mondrian_calibrate(p_active_calib, y_calib)

    query_path = Path(args.query_csv)
    query_lib = pd.read_csv(query_path)
    smi_col, id_col = args.query_smiles_col, args.query_id_col
    q_smiles = [s for s in query_lib[smi_col].tolist() if isinstance(s, str) and s]
    print(f"query library: {query_path.name} ({len(q_smiles)}/{len(query_lib)} valid SMILES)")
    X_q = qsar.compute_raw_features(q_smiles)
    p_active_q = pipe.predict_proba(X_q)[:, 1]
    pvals_q = mondrian_pvalues(p_active_q, calib_scores)

    out_df = pd.DataFrame({
        id_col: query_lib.loc[query_lib[smi_col].isin(q_smiles), id_col].values[: len(q_smiles)],
        "smiles": q_smiles,
        "p_active_calibrated": np.round(p_active_q, 4),
        "conformal_p_active": np.round(pvals_q[1], 4),
        "conformal_p_inactive": np.round(pvals_q[0], 4),
    })
    for sig in SIGNIFICANCE_LEVELS:
        out_df[f"prediction_set_sig{sig}"] = prediction_set(pvals_q, sig)

    print(f"{query_path.stem} prediction-set distribution by significance level:")
    for sig in SIGNIFICANCE_LEVELS:
        counts = out_df[f"prediction_set_sig{sig}"].value_counts().to_dict()
        print(f"  sig={sig:.2f} ({100*(1-sig):.0f}% confidence): {counts}")

    out_df.to_csv(outdir / args.out_name, index=False)
    meta = dict(
        pipeline_path=str(PIPELINE_PATH), pipeline_sha256=sha256(PIPELINE_PATH),
        calibration_set_n=len(smiles_calib), calibration_set_source=str(TRAINING_SET_CSV),
        calibration_set_note="identical scaffold-disjoint held-out split already used to "
                             "report external AUC=0.922; reused for conformal calibration, "
                             "not for re-training",
        method="Mondrian (label-conditional) inductive conformal prediction",
        significance_levels=SIGNIFICANCE_LEVELS,
        query_library=str(query_path), n_query=len(q_smiles),
        exchangeability_caveat=(
            f"E2 (ADR-AIIM-073) found SNP-Lib 0% in-domain across six independent "
            f"structural descriptor spaces vs. the ChEMBL training/calibration "
            f"distribution; the natural-product library expansion (ADR-AIIM-079) is "
            f"drawn from the same biosynthetic classes and is not independently "
            f"re-verified for exchangeability here. Conformal validity requires "
            f"exchangeability between calibration and query data, which this "
            f"divergence calls into question for any natural-product query library. "
            f"The validity check below confirms the conformal MECHANISM is correctly "
            f"implemented on held-out ChEMBL data; it does not and cannot certify "
            f"nominal coverage for {query_path.stem} specifically."
        ),
        validity_check=validity,
    )
    summary_name = f"conformal_qsar_summary_{query_path.stem}.json"
    (outdir / summary_name).write_text(json.dumps(meta, indent=2, default=str))
    print(f"\nwritten -> {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

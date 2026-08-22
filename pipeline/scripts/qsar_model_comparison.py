"""
qsar_model_comparison.py — is RandomForest actually the best algorithm for
this QSAR task, or was it never benchmarked against alternatives?

The frozen pipeline (`qsar_pipeline_a2a.pkl`) uses RandomForestClassifier,
inherited from a legacy script without a documented comparison against other
algorithm families. This script runs the identical scaffold-grouped nested
CV (5 outer StratifiedGroupKFold, 3 inner StratifiedKFold, ECFP4+FCFP4
features, VarianceThreshold inside the CV loop) against four algorithm
families, so the comparison is apples-to-apples with the existing frozen
model's own reported CV numbers (AUC 0.907 +/- 0.017).

Deep-learning / graph-neural-network approaches are deliberately excluded,
consistent with this project's own governance
(`11_AGENTIC_SKILLS_GOVERNANCE.md`: "the current number of independent fault
families does not justify deep learning") -- that reasoning applies equally
here (~6,000 compounds, 7.5:1 class imbalance).

SVM uses a linear kernel (LinearSVC, calibrated via CalibratedClassifierCV
for probability outputs) rather than RBF: at ~4,750 training compounds x
~4,096 binary fingerprint bits, an RBF-kernel SVC's O(n^2)-O(n^3) training
cost is not tractable for a 5-fold x hyperparameter-search comparison in a
reasonable wall time. This is a disclosed simplification, not a claim that
linear SVM represents SVM's full potential on this task.

Per-fold hyperparameter search is intentionally smaller than the original
RF-only search (n_iter reduced) to keep four algorithms' combined runtime
bounded; this is a screening comparison to check whether RF is competitive,
not a final from-scratch tuning of the winner.

Usage:
    python qsar_model_comparison.py --out <dir>
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import time
from pathlib import Path

# Cap thread usage globally (must be set before numpy/sklearn import): the
# SBDD full-protocol redock is running concurrently on the same physical CPU
# via WSL2 and is the priority workload. Without this, HistGradientBoosting's
# internal OpenMP threading and BLAS-backed linear algebra would each try to
# claim all available cores regardless of RandomForestClassifier's own
# n_jobs setting.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "3")

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.feature_selection import VarianceThreshold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, matthews_corrcoef, roc_auc_score
from sklearn.model_selection import RandomizedSearchCV, StratifiedGroupKFold, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.svm import LinearSVC
from sklearn.calibration import CalibratedClassifierCV

ROOT = Path(r"D:\Agentic VS Workflow")
QSAR_RUN = ROOT / "20_new_runs/qsar_lbdd_2026-08-07"
TRAINING_SET_CSV = QSAR_RUN / "chembl_a2a_training_set.csv"

VAR_THRESHOLD = 0.01
N_OUTER_FOLDS = 5
N_INNER_FOLDS = 3
N_HP_ITER = 6           # reduced from the frozen model's 20, for a 4-algorithm screen
RANDOM_STATE = 42

FROZEN_RF_CV_AUC = 0.907  # for reference in the printed comparison


def _load_qsar_module():
    path = Path(__file__).parent / "run_qsar_lbdd_2026-08-07.py"
    spec = importlib.util.spec_from_file_location("_qsar_lbdd_source", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def make_pipe(name: str, random_state: int = RANDOM_STATE) -> tuple[Pipeline, dict]:
    var = ("var", VarianceThreshold(threshold=VAR_THRESHOLD))
    if name == "random_forest":
        # n_jobs deliberately capped (not -1): the SBDD full-protocol redock is
        # running concurrently on the same physical CPU (via WSL2) and is the
        # priority workload -- this comparison should not compete for all cores.
        clf = RandomForestClassifier(class_weight="balanced", random_state=random_state, n_jobs=3)
        space = {"clf__n_estimators": [100, 200, 500], "clf__max_depth": [None, 15, 25],
                "clf__min_samples_leaf": [1, 5, 10]}
    elif name == "hist_gradient_boosting":
        clf = HistGradientBoostingClassifier(random_state=random_state, class_weight="balanced")
        space = {"clf__max_iter": [100, 200, 300], "clf__max_depth": [None, 6, 10],
                "clf__learning_rate": [0.05, 0.1, 0.2]}
    elif name == "linear_svm":
        base = LinearSVC(class_weight="balanced", random_state=random_state, max_iter=5000, dual="auto")
        clf = CalibratedClassifierCV(base, method="sigmoid", cv=3)
        space = {"clf__estimator__C": [0.01, 0.1, 1.0, 10.0]}
    elif name == "logistic_regression":
        clf = LogisticRegression(penalty="l2", class_weight="balanced", max_iter=2000,
                                 random_state=random_state, solver="lbfgs")
        space = {"clf__C": [0.01, 0.1, 1.0, 10.0, 100.0]}
    else:
        raise ValueError(name)
    return Pipeline([var, ("clf", clf)]), space


def nested_cv_evaluate(name: str, X: np.ndarray, y: np.ndarray, groups: np.ndarray) -> dict:
    outer_cv = StratifiedGroupKFold(n_splits=N_OUTER_FOLDS, shuffle=True, random_state=RANDOM_STATE)
    inner_cv = StratifiedKFold(n_splits=N_INNER_FOLDS, shuffle=True, random_state=RANDOM_STATE)
    fold_metrics = {k: [] for k in ("roc_auc", "mcc", "bal_acc")}
    t0 = time.time()
    for fold, (tr, val) in enumerate(outer_cv.split(X, y, groups=groups)):
        pipe, space = make_pipe(name)
        search = RandomizedSearchCV(pipe, space, n_iter=min(N_HP_ITER, _space_size(space)),
                                    cv=inner_cv, scoring="roc_auc", random_state=RANDOM_STATE,
                                    n_jobs=1, refit=True)
        search.fit(X[tr], y[tr])
        best = search.best_estimator_
        y_prob = best.predict_proba(X[val])[:, 1]
        y_pred = best.predict(X[val])
        fold_metrics["roc_auc"].append(roc_auc_score(y[val], y_prob))
        fold_metrics["mcc"].append(matthews_corrcoef(y[val], y_pred))
        fold_metrics["bal_acc"].append(balanced_accuracy_score(y[val], y_pred))
        print(f"  [{name}] fold {fold+1}/{N_OUTER_FOLDS}: AUC={fold_metrics['roc_auc'][-1]:.3f}", flush=True)
    elapsed = time.time() - t0
    out = {f"{k}_mean": float(np.mean(v)) for k, v in fold_metrics.items()}
    out.update({f"{k}_sd": float(np.std(v)) for k, v in fold_metrics.items()})
    out["elapsed_s"] = round(elapsed, 1)
    return out


def _space_size(space: dict) -> int:
    n = 1
    for v in space.values():
        n *= len(v)
    return n


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(QSAR_RUN))
    args = ap.parse_args()
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    qsar = _load_qsar_module()
    combined = pd.read_csv(TRAINING_SET_CSV)
    smiles_all = combined["canonical_smiles"].tolist()
    y_all = combined["label"].to_numpy()
    train_idx, _ = qsar.scaffold_train_test_split(smiles_all, y_all)
    smiles_train = [smiles_all[i] for i in train_idx]
    y_train = y_all[train_idx]
    print(f"training set (identical to the frozen model's own split): "
          f"{len(train_idx)} compounds, {int(y_train.sum())} active / {int((y_train==0).sum())} inactive")

    X_train = qsar.compute_raw_features(smiles_train)
    groups_train = qsar.get_scaffold_groups(smiles_train)

    algorithms = ["random_forest", "hist_gradient_boosting", "linear_svm", "logistic_regression"]
    results = {}
    for name in algorithms:
        print(f"\n=== {name} ===")
        results[name] = nested_cv_evaluate(name, X_train, y_train, groups_train)

    print(f"\n{'algorithm':>24} {'AUC':>16} {'MCC':>16} {'Bal.Acc':>16} {'time(s)':>8}")
    print("-" * 90)
    for name in algorithms:
        r = results[name]
        print(f"{name:>24} {r['roc_auc_mean']:.3f} +/- {r['roc_auc_sd']:.3f}   "
              f"{r['mcc_mean']:.3f} +/- {r['mcc_sd']:.3f}   "
              f"{r['bal_acc_mean']:.3f} +/- {r['bal_acc_sd']:.3f}   {r['elapsed_s']:>8.0f}")
    print(f"\nreference: frozen RandomForest model's own reported CV AUC = {FROZEN_RF_CV_AUC:.3f} "
          f"(20-iter search, not directly comparable to this {N_HP_ITER}-iter screen)")

    best = max(algorithms, key=lambda n: results[n]["roc_auc_mean"])
    print(f"\nbest by mean CV AUC in this screen: {best} "
          f"({results[best]['roc_auc_mean']:.3f})")

    (outdir / "qsar_model_comparison.json").write_text(json.dumps({
        "note": "Screening comparison with reduced hyperparameter search "
                f"(n_iter={N_HP_ITER}) for computational tractability across 4 algorithms; "
                "not a final from-scratch tuning of the winning algorithm.",
        "n_outer_folds": N_OUTER_FOLDS, "n_inner_folds": N_INNER_FOLDS, "n_hp_iter": N_HP_ITER,
        "training_set_n": len(train_idx),
        "results": results, "best_by_auc": best,
    }, indent=2, default=str))
    print(f"\nwritten -> {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

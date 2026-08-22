"""
qsar_y_randomization.py — y-randomisation of the frozen A2AR QSAR model.

Following Rücker, Rücker & Meringer (J Chem Inf Model 2007;47:2345-57,
10.1021/ci700157b): the real model's performance is compared against models
built on permuted response using the *same descriptor pool and the same
model-building procedure*. A model whose real performance does not clearly
exceed the permuted distribution is not learning structure-activity signal.

This reuses the frozen pipeline's own feature computation and scaffold-grouped
splitting (`run_qsar_lbdd_2026-08-07.py`), so the permuted runs differ from the
real one only in the label vector.

Thread use is capped: docking may be running concurrently on the same host.

Usage:
    python qsar_y_randomization.py --n-permutations 100 --out <dir>
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "4")

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import VarianceThreshold
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import Pipeline


def _root() -> Path:
    marker = "00_project_control/08_DECISION_LOG.md"
    for cand in (Path(r"D:\Agentic VS Workflow"), Path("/mnt/d/Agentic VS Workflow")):
        if (cand / marker).is_file():
            return cand
    raise SystemExit("project root not found")


ROOT = _root()
QSAR_RUN = ROOT / "20_new_runs/qsar_lbdd_2026-08-07"
TRAIN_CSV = QSAR_RUN / "chembl_a2a_training_set.csv"
VAR_THRESHOLD = 0.01
N_FOLDS = 5
RANDOM_STATE = 42


def _load_qsar_module():
    path = ROOT / "10_new_system/scripts/run_qsar_lbdd_2026-08-07.py"
    spec = importlib.util.spec_from_file_location("_qsar_src", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def cv_auc(X: np.ndarray, y: np.ndarray, groups: np.ndarray, seed: int) -> float:
    """One scaffold-grouped CV estimate under fixed hyperparameters.

    Hyperparameters are held fixed rather than re-searched: y-randomisation
    asks whether the *same procedure* extracts signal from a permuted
    response, and re-tuning on each permutation would let the search chase
    noise and inflate the permuted scores.
    """
    cv = StratifiedGroupKFold(n_splits=N_FOLDS, shuffle=True, random_state=RANDOM_STATE)
    aucs = []
    for tr, va in cv.split(X, y, groups=groups):
        pipe = Pipeline([
            ("var", VarianceThreshold(threshold=VAR_THRESHOLD)),
            ("clf", RandomForestClassifier(n_estimators=200, class_weight="balanced",
                                           random_state=seed, n_jobs=4)),
        ])
        pipe.fit(X[tr], y[tr])
        aucs.append(roc_auc_score(y[va], pipe.predict_proba(X[va])[:, 1]))
    return float(np.mean(aucs))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-permutations", type=int, default=100)
    ap.add_argument("--out", default=str(QSAR_RUN))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    qsar = _load_qsar_module()
    combined = pd.read_csv(TRAIN_CSV)
    smiles_all = combined["canonical_smiles"].tolist()
    y_all = combined["label"].to_numpy()
    train_idx, _ = qsar.scaffold_train_test_split(smiles_all, y_all)
    smiles = [smiles_all[i] for i in train_idx]
    y = y_all[train_idx]

    print(f"training set: {len(smiles)} compounds "
          f"({int(y.sum())} active / {int((y == 0).sum())} inactive)")
    X = qsar.compute_raw_features(smiles)
    groups = qsar.get_scaffold_groups(smiles)

    t0 = time.time()
    real = cv_auc(X, y, groups, seed=RANDOM_STATE)
    print(f"real model CV AUC = {real:.4f}  ({time.time()-t0:.0f}s per run)")

    rng = np.random.default_rng(RANDOM_STATE)
    perm: list[float] = []
    for i in range(args.n_permutations):
        y_perm = rng.permutation(y)
        perm.append(cv_auc(X, y_perm, groups, seed=RANDOM_STATE + i + 1))
        if (i + 1) % 10 == 0:
            print(f"  permutation {i+1}/{args.n_permutations}: "
                  f"mean so far {np.mean(perm):.4f}, max {np.max(perm):.4f}", flush=True)

    perm_arr = np.asarray(perm)
    # Empirical p: fraction of permuted models reaching the real model's score.
    # +1 in numerator and denominator is the standard non-zero-p correction.
    p_emp = float((np.sum(perm_arr >= real) + 1) / (len(perm_arr) + 1))

    result = {
        "reference": "Rucker C, Rucker G, Meringer M. J Chem Inf Model 2007;47:2345-57. "
                     "doi:10.1021/ci700157b",
        "procedure": "Response permuted; identical descriptor pool, identical "
                     "scaffold-grouped CV, hyperparameters held fixed so the "
                     "permuted runs cannot re-tune on noise.",
        "n_permutations": int(len(perm_arr)),
        "real_cv_auc": real,
        "permuted_mean": float(perm_arr.mean()),
        "permuted_sd": float(perm_arr.std()),
        "permuted_max": float(perm_arr.max()),
        "permuted_min": float(perm_arr.min()),
        "gap_real_minus_permuted_mean": float(real - perm_arr.mean()),
        "z_vs_permuted": float((real - perm_arr.mean()) / perm_arr.std())
        if perm_arr.std() > 0 else None,
        "empirical_p": p_emp,
        "permuted_aucs": [round(v, 5) for v in perm],
    }
    (out / "qsar_y_randomization.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8")

    print(f"\nreal            : {real:.4f}")
    print(f"permuted mean+-sd: {perm_arr.mean():.4f} +- {perm_arr.std():.4f}")
    print(f"permuted max    : {perm_arr.max():.4f}")
    print(f"gap             : {real - perm_arr.mean():.4f}")
    print(f"empirical p     : {p_emp:.4f}")
    print(f"\nwritten -> {out/'qsar_y_randomization.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

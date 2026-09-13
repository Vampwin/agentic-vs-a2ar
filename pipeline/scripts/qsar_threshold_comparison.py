"""
qsar_threshold_comparison.py — does a natural-product-appropriate hit
threshold change the QSAR model's behaviour vs. the frozen medicinal-
chemistry threshold (pChEMBL >=6.0/<=5.5, i.e. <=1 uM active / >=3.16 uM
inactive)?

Motivation (human author, this session): the frozen threshold labels
caffeine (all 21 records, pChEMBL 4.30-5.61) and galangin (all 3 records,
pChEMBL 4.48-4.78) -- both well-documented natural-product adenosine-
receptor ligands -- as inactive, because it was calibrated for optimised
synthetic hit-to-lead potency, not for natural products screened without
subsequent medicinal-chemistry optimisation. NP/HTS screening literature
commonly uses a looser "hit" convention, often <=10 uM (pChEMBL >=5.0).

This script builds a SECOND training set from the identical raw ChEMBL pull
using ACTIVE_CUTOFF=5.0 / INACTIVE_CUTOFF=4.5 (preserving the same 0.5-log-
unit grey-zone gap as the frozen threshold, shifted down by one log unit),
trains a RandomForest under the identical scaffold-grouped nested-CV
protocol used for the frozen model, and applies Mondrian conformal
prediction (identical method to ADR-AIIM-077/081) to SNP-Lib and the
natural-product library expansion, so the two threshold regimes' downstream
predictions are directly comparable.

The frozen model (`qsar_pipeline_a2a.pkl`) is NOT modified or replaced by
this script; this is a side-by-side comparison, not a swap.

Usage:
    python qsar_threshold_comparison.py --out <dir>
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import time
from pathlib import Path
from kdss_paths import project_root, resolve, under, dev_dir

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "3")

import joblib
import numpy as np
import pandas as pd
from rdkit import Chem
from sklearn.calibration import CalibratedClassifierCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import VarianceThreshold
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import RandomizedSearchCV, StratifiedGroupKFold, StratifiedKFold
from sklearn.pipeline import Pipeline

ROOT = project_root()
FROZEN = dev_dir("02_frozen_inputs")
CHEMBL_RAW_SOURCE = under("02_frozen_inputs", "public_benchmarks/chembl_a2ar/chembl_a2a_raw_v2.csv")
SNPLIB_LEDGER = resolve("20_new_runs/snp_screen_2026-08-04/G1_standardization/snp_source_ledger.csv")
NPLIB_LEDGER = resolve("20_new_runs/np_library_funnel_2026-08-15/np_library_dock_ledger.csv")

# Frozen (original, medicinal-chemistry-calibrated) regime
FROZEN_ACTIVE, FROZEN_INACTIVE = 6.0, 5.5
# Alternative (NP/HTS-screening-appropriate) regime -- same 0.5-log grey-zone
# gap, shifted down 1.0 log unit (10x looser), matching a commonly cited
# <=10 uM "hit" convention for unoptimised natural-product screening.
NP_ACTIVE, NP_INACTIVE = 5.0, 4.5

N_OUTER_FOLDS, N_INNER_FOLDS, N_HP_ITER = 5, 3, 8
VAR_THRESHOLD = 0.01
RANDOM_STATE = 42
SIGNIFICANCE_LEVELS = [0.05, 0.10, 0.20]

REF_COMPOUNDS = {  # for tracking specific NP actives across threshold regimes
    "caffeine": "Cn1c(=O)c2c(ncn2C)n(C)c1=O",
    "galangin": "O=c1c(O)c(-c2ccccc2)oc2cc(O)cc(O)c12",
}


def _load_qsar_module():
    path = Path(__file__).parent / "run_qsar_lbdd_2026-08-07.py"
    spec = importlib.util.spec_from_file_location("_qsar_lbdd_source", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def build_training_set(raw: pd.DataFrame, active_cut: float, inactive_cut: float,
                       qsar) -> pd.DataFrame:
    df = raw.copy()
    act_rel = df["standard_relation"].isin(["=", "<=", "<"])
    inact_rel = df["standard_relation"].isin(["=", ">=", ">"])
    actives = df[act_rel & (df["pchembl_value"] >= active_cut)].copy()
    inactives = df[inact_rel & (df["pchembl_value"] <= inactive_cut)].copy()
    actives["label"] = 1
    inactives["label"] = 0
    combined = pd.concat([actives, inactives], ignore_index=True)
    combined["inchikey"] = combined["canonical_smiles"].apply(qsar._smiles_to_inchikey)
    combined = combined.dropna(subset=["inchikey", "canonical_smiles"])
    combined = (combined.sort_values("pchembl_value", ascending=False)
               .drop_duplicates(subset=["inchikey"]).reset_index(drop=True))
    return combined


def make_rf_pipe(random_state: int = RANDOM_STATE) -> Pipeline:
    return Pipeline([
        ("var", VarianceThreshold(threshold=VAR_THRESHOLD)),
        ("rf", RandomForestClassifier(class_weight="balanced", random_state=random_state, n_jobs=3)),
    ])


def nested_cv_auc(X: np.ndarray, y: np.ndarray, groups: np.ndarray) -> dict:
    space = {"rf__n_estimators": [100, 200, 500], "rf__max_depth": [None, 15, 25],
            "rf__min_samples_leaf": [1, 5, 10]}
    outer = StratifiedGroupKFold(n_splits=N_OUTER_FOLDS, shuffle=True, random_state=RANDOM_STATE)
    inner = StratifiedKFold(n_splits=N_INNER_FOLDS, shuffle=True, random_state=RANDOM_STATE)
    aucs = []
    for tr, val in outer.split(X, y, groups=groups):
        search = RandomizedSearchCV(make_rf_pipe(), space, n_iter=N_HP_ITER, cv=inner,
                                    scoring="roc_auc", random_state=RANDOM_STATE, n_jobs=1)
        search.fit(X[tr], y[tr])
        y_prob = search.best_estimator_.predict_proba(X[val])[:, 1]
        aucs.append(roc_auc_score(y[val], y_prob))
    return {"cv_auc_mean": float(np.mean(aucs)), "cv_auc_sd": float(np.std(aucs))}


def mondrian_calibrate(p_active: np.ndarray, y: np.ndarray) -> dict[int, np.ndarray]:
    scores = {}
    for cls in (0, 1):
        mask = y == cls
        p_true = np.where(cls == 1, p_active, 1 - p_active)[mask]
        scores[cls] = np.sort(1.0 - p_true)
    return scores


def mondrian_pvalues(p_active: np.ndarray, calib: dict[int, np.ndarray]) -> dict[int, np.ndarray]:
    out = {}
    for cls in (0, 1):
        p_y = np.where(cls == 1, p_active, 1 - p_active)
        alpha_xy = 1.0 - p_y
        n_y = len(calib[cls])
        n_ge = n_y - np.searchsorted(calib[cls], alpha_xy, side="left")
        out[cls] = (1.0 + n_ge) / (1.0 + n_y)
    return out


def prediction_set(pvals: dict[int, np.ndarray], sig: float) -> list[str]:
    n = len(pvals[0])
    return ["{" + ",".join(str(c) for c in (0, 1) if pvals[c][i] > sig) + "}"
           if any(pvals[c][i] > sig for c in (0, 1)) else "{}"
           for i in range(n)]


def run_regime(name: str, active_cut: float, inactive_cut: float, raw: pd.DataFrame,
               qsar, query_libs: dict[str, tuple[str, str, str]]) -> dict:
    print(f"\n{'='*70}\nregime: {name}  (active>={active_cut}, inactive<={inactive_cut})\n{'='*70}")
    combined = build_training_set(raw, active_cut, inactive_cut, qsar)
    n_act, n_ina = (combined["label"] == 1).sum(), (combined["label"] == 0).sum()
    print(f"training set: {n_act} active / {n_ina} inactive (n={len(combined)})")

    for name_c, smi in REF_COMPOUNDS.items():
        ik = Chem.MolToInchiKey(Chem.MolFromSmiles(smi))
        hit = combined[combined["inchikey"] == ik]
        if len(hit):
            print(f"  {name_c}: {len(hit)} record(s), label={hit['label'].iloc[0]}, "
                  f"pChEMBL={hit['pchembl_value'].tolist()}")
        else:
            print(f"  {name_c}: not present under this threshold (grey zone or absent)")

    smiles_all = combined["canonical_smiles"].tolist()
    y_all = combined["label"].to_numpy()
    train_idx, test_idx = qsar.scaffold_train_test_split(smiles_all, y_all)
    smiles_train = [smiles_all[i] for i in train_idx]
    y_train = y_all[train_idx]
    smiles_calib = [smiles_all[i] for i in test_idx]
    y_calib = y_all[test_idx]
    print(f"scaffold split: {len(train_idx)} train / {len(test_idx)} calibration "
          f"(calib: {int(y_calib.sum())} active / {int((y_calib==0).sum())} inactive)")

    X_train = qsar.compute_raw_features(smiles_train)
    groups_train = qsar.get_scaffold_groups(smiles_train)
    t0 = time.time()
    cv = nested_cv_auc(X_train, y_train, groups_train)
    print(f"nested CV AUC: {cv['cv_auc_mean']:.3f} +/- {cv['cv_auc_sd']:.3f} "
          f"({time.time()-t0:.0f}s)")

    final_pipe = make_rf_pipe()
    final_pipe.fit(X_train, y_train)
    X_calib = qsar.compute_raw_features(smiles_calib)
    p_active_calib = final_pipe.predict_proba(X_calib)[:, 1]
    test_auc = roc_auc_score(y_calib, p_active_calib) if len(set(y_calib)) > 1 else float("nan")
    print(f"held-out scaffold-split test AUC: {test_auc:.3f} (n={len(y_calib)})")

    calib_scores = mondrian_calibrate(p_active_calib, y_calib)

    query_results = {}
    for qname, (path, smi_col, id_col) in query_libs.items():
        qdf = pd.read_csv(path)
        q_smiles = [s for s in qdf[smi_col].tolist() if isinstance(s, str) and s]
        X_q = qsar.compute_raw_features(q_smiles)
        p_active_q = final_pipe.predict_proba(X_q)[:, 1]
        pvals_q = mondrian_pvalues(p_active_q, calib_scores)
        sets_005 = prediction_set(pvals_q, 0.05)
        counts = pd.Series(sets_005).value_counts().to_dict()
        print(f"  [{qname}] n={len(q_smiles)}  sig=0.05 prediction sets: {counts}  "
              f"mean P(active)={p_active_q.mean():.3f}")
        query_results[qname] = {
            "n": len(q_smiles), "prediction_sets_sig0.05": counts,
            "mean_p_active": round(float(p_active_q.mean()), 4),
        }

    return {
        "active_cutoff": active_cut, "inactive_cutoff": inactive_cut,
        "n_active": int(n_act), "n_inactive": int(n_ina),
        "cv_auc_mean": cv["cv_auc_mean"], "cv_auc_sd": cv["cv_auc_sd"],
        "held_out_test_auc": float(test_auc), "n_calib": len(smiles_calib),
        "query_results": query_results,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(resolve("20_new_runs/qsar_lbdd_2026-08-07")))
    args = ap.parse_args()
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    qsar = _load_qsar_module()
    raw = pd.read_csv(CHEMBL_RAW_SOURCE)
    print(f"raw ChEMBL A2AR pull: {len(raw)} records")

    query_libs = {
        "SNP-Lib": (str(SNPLIB_LEDGER), "smiles_canonical", "SNP_ID"),
        "NP-library": (str(NPLIB_LEDGER), "SMILES", "compound_id"),
    }

    results = {}
    results["frozen_6.0_5.5"] = run_regime("frozen (6.0/5.5)", FROZEN_ACTIVE, FROZEN_INACTIVE,
                                           raw, qsar, query_libs)
    results["np_appropriate_5.0_4.5"] = run_regime("NP-appropriate (5.0/4.5)", NP_ACTIVE, NP_INACTIVE,
                                                    raw, qsar, query_libs)

    print(f"\n{'='*70}\nSUMMARY\n{'='*70}")
    print(f"{'regime':>24} {'n_active':>9} {'n_inactive':>11} {'CV AUC':>16} {'test AUC':>9}")
    for key, r in results.items():
        print(f"{key:>24} {r['n_active']:>9} {r['n_inactive']:>11} "
              f"{r['cv_auc_mean']:.3f} +/- {r['cv_auc_sd']:.3f}   {r['held_out_test_auc']:>9.3f}")

    (outdir / "qsar_threshold_comparison.json").write_text(json.dumps({
        "note": "Side-by-side comparison; the frozen model/pipeline "
                "(qsar_pipeline_a2a.pkl) is unchanged by this script.",
        "reference_compounds_checked": list(REF_COMPOUNDS),
        "results": results,
    }, indent=2, default=str))
    print(f"\nwritten -> {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""conformal_diagnostics.py -- C3. What the conformal predictor does off-domain.

The manuscript reports that at alpha = 0.05 the predictor assigns the singleton
{inactive} to 4,977 of 5,000 query compounds, and carries one caveat: nominal
coverage cannot be certified because exchangeability between calibration and
query is not established. That caveat is true and is weaker than the situation
warrants.

A Mondrian inductive conformal predictor has a way of saying "this object is
atypical of both classes": it returns the empty set, when both class p-values
fall below alpha. On held-out ChEMBL the mechanism does exactly that. On a
natural-product library that is 0% in-domain in four structural spaces, it never
does. That contrast is the diagnostic the manuscript is missing, and it is
computable from artifacts that already exist -- no re-docking, no re-training.

This script computes it. It reuses conformal_qsar.py's own calibration and
p-value functions rather than reimplementing them, so the numbers it reports
are produced by the same code that produced the released prediction sets.

Usage:
    python conformal_diagnostics.py [--out DIR]
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
from kdss_paths import project_root, resolve, under, dev_dir

import joblib
import numpy as np
import pandas as pd


def _root() -> Path:
    # Delegates to kdss_paths, which finds whichever layout this file is inside
    # -- the development tree or a deposit -- by walking up from its own
    # location. It previously tried two absolute paths and exited if neither
    # existed, which on the author's machine meant a copy of the deposit read
    # and wrote the development tree instead of itself.
    return project_root()


ROOT = _root()
RUNS = dev_dir("20_new_runs")
QSAR_RUN = under("20_new_runs", "qsar_lbdd_2026-08-07")
NP_RUN = under("20_new_runs", "np_library_funnel_2026-08-15")

SIGS = [0.05, 0.10, 0.20]
N_FOLDS = 5
SEED = 42


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def set_profile(p_inactive: np.ndarray, p_active: np.ndarray, sig: float) -> dict:
    """Prediction-set composition and efficiency at one significance level."""
    keep0 = p_inactive > sig
    keep1 = p_active > sig
    n = len(p_inactive)
    size = keep0.astype(int) + keep1.astype(int)
    return {
        "n": int(n),
        "empty_{}": round(float((~keep0 & ~keep1).mean()), 4),
        "singleton_inactive_{0}": round(float((keep0 & ~keep1).mean()), 4),
        "singleton_active_{1}": round(float((~keep0 & keep1).mean()), 4),
        "both_{0,1}": round(float((keep0 & keep1).mean()), 4),
        "mean_set_size": round(float(size.mean()), 4),
    }


def dist(x: np.ndarray) -> dict:
    q = np.percentile(x, [5, 25, 50, 75, 95])
    return {"min": round(float(x.min()), 4), "p05": round(float(q[0]), 4),
            "q1": round(float(q[1]), 4), "median": round(float(q[2]), 4),
            "q3": round(float(q[3]), 4), "p95": round(float(q[4]), 4),
            "max": round(float(x.max()), 4), "mean": round(float(x.mean()), 4)}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(under("20_new_runs", "conformal_diagnostics_2026-08-25")))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    scripts = resolve("10_new_system/scripts")
    cq = _load("conformal_qsar", scripts / "conformal_qsar.py")
    qsar = cq._load_qsar_module()

    pipe = joblib.load(cq.PIPELINE_PATH)
    print(f"frozen pipeline: {cq.PIPELINE_PATH.name} "
          f"(sha256={sha256(cq.PIPELINE_PATH)[:16]}...) -- not refit")

    combined = pd.read_csv(cq.TRAINING_SET_CSV)
    smiles_all = combined["canonical_smiles"].tolist()
    y_all = combined["label"].to_numpy()
    _, test_idx = qsar.scaffold_train_test_split(smiles_all, y_all)
    smiles_cal = [smiles_all[i] for i in test_idx]
    y_cal = y_all[test_idx]
    X_cal = qsar.compute_raw_features(smiles_cal)
    p_act_cal = pipe.predict_proba(X_cal)[:, 1]
    groups_cal = qsar.get_scaffold_groups(smiles_cal)
    print(f"calibration set reproduced: n={len(smiles_cal)} "
          f"({int(y_cal.sum())} active / {int((y_cal == 0).sum())} inactive)")

    # ---- honest calibration-set p-values, by the same 5-fold scaffold-grouped
    # cross-conformal split validity_check() uses. A point's p-values must come
    # from a calibration set that does not contain it.
    rng = np.random.default_rng(SEED)
    groups = np.unique(groups_cal)
    rng.shuffle(groups)
    fold_of = {g: i % N_FOLDS for i, g in enumerate(groups)}
    fold_id = np.array([fold_of[g] for g in groups_cal])

    pv_cal = {0: np.full(len(y_cal), np.nan), 1: np.full(len(y_cal), np.nan)}
    for f in range(N_FOLDS):
        te, ca = fold_id == f, fold_id != f
        if te.sum() == 0 or ca.sum() == 0:
            continue
        cs = cq.mondrian_calibrate(p_act_cal[ca], y_cal[ca])
        pv = cq.mondrian_pvalues(p_act_cal[te], cs)
        pv_cal[0][te], pv_cal[1][te] = pv[0], pv[1]
    assert not np.isnan(pv_cal[0]).any(), "some calibration point got no fold"

    pd.DataFrame({
        "smiles": smiles_cal, "label": y_cal, "fold": fold_id,
        "p_active_calibrated": np.round(p_act_cal, 4),
        "conformal_p_active": np.round(pv_cal[1], 4),
        "conformal_p_inactive": np.round(pv_cal[0], 4),
    }).to_csv(out / "calibration_pvalues.csv", index=False)

    # ---- the two query libraries, read from the released artifacts ----------
    np_lib = pd.read_csv(NP_RUN / "nplibrary_conformal_predictions.csv")
    snp = pd.read_csv(QSAR_RUN / "snplib_conformal_predictions.csv")

    panels = {
        "calibration_chembl_heldout": (pv_cal[0], pv_cal[1]),
        "snplib_236": (snp.conformal_p_inactive.to_numpy(),
                       snp.conformal_p_active.to_numpy()),
        "np_library_5000": (np_lib.conformal_p_inactive.to_numpy(),
                            np_lib.conformal_p_active.to_numpy()),
    }

    report = {
        "significance_levels": SIGS,
        # The Mondrian guarantee is per class, so the size of each calibration
        # class is what the corresponding p-value is estimated from. They are
        # very unequal here and that is the whole mechanism.
        "calibration_set": {
            "n": int(len(y_cal)),
            "n_active": int((y_cal == 1).sum()),
            "n_inactive": int((y_cal == 0).sum()),
        },
        "panels": {},
    }
    for name, (p0, p1) in panels.items():
        report["panels"][name] = {
            "p_inactive_distribution": dist(p0),
            "p_active_distribution": dist(p1),
            "by_significance": {str(s): set_profile(p0, p1, s) for s in SIGS},
        }

    # ---- the comparison the manuscript needs, stated as one number ---------
    e_cal = report["panels"]["calibration_chembl_heldout"]["by_significance"]["0.2"]["empty_{}"]
    e_np = report["panels"]["np_library_5000"]["by_significance"]["0.2"]["empty_{}"]
    e_snp = report["panels"]["snplib_236"]["by_significance"]["0.2"]["empty_{}"]
    report["headline"] = {
        "empty_set_rate_at_alpha_0.20": {
            "calibration_chembl_heldout": e_cal,
            "snplib_236": e_snp,
            "np_library_5000": e_np,
        },
        "reading":
            "The empty set is how a Mondrian conformal predictor says an object "
            "is atypical of BOTH classes. On held-out ChEMBL the mechanism "
            f"returns it for {e_cal:.1%} of compounds at alpha=0.20. On the two "
            "natural-product libraries -- which are 0% in-domain in four "
            "structural descriptor spaces against this same calibration set -- "
            f"it returns it for {e_snp:.1%} and {e_np:.1%}. The predictor is not "
            "abstaining on the query libraries; it is confidently assigning them "
            "to the inactive class. That is a property of the nonconformity "
            "measure off-domain, not a finding about the compounds.",
    }
    # ---- where the Tier-A criterion actually sits in this distribution -----
    # The tier rule reads "a ligand-based prediction set that does not
    # confidently exclude activity", which is p_active > alpha at alpha = 0.05.
    # Whether that is a signal depends entirely on the range p_active spans.
    e8 = pd.read_csv(NP_RUN / "e8_np_library_prioritization.csv")
    tierA = e8.loc[e8.tier == "A_converging_sbdd_lbdd", "compound_id"].tolist()
    pa = np_lib.conformal_p_active
    flagged = np_lib[pa > 0.05]
    tsub = np_lib[np_lib.compound_id.isin(tierA)]
    report["tier_a_criterion"] = {
        "alpha": 0.05,
        "n_flagged": int(len(flagged)),
        "library_p_active_max": round(float(pa.max()), 4),
        "library_p_active_median": round(float(pa.median()), 4),
        "flagged_p_active_range": [round(float(flagged.conformal_p_active.min()), 4),
                                   round(float(flagged.conformal_p_active.max()), 4)],
        "p_active_granularity": round(1.0 / (1 + int((y_cal == 1).sum())), 6),
        "tier_a": [
            {"compound_id": r.compound_id,
             "p_active_calibrated": round(float(r.p_active_calibrated), 4),
             "conformal_p_active": round(float(r.conformal_p_active), 4),
             "conformal_p_inactive": round(float(r.conformal_p_inactive), 4),
             "rank_by_p_active_of_5000": int(pa.rank(ascending=False,
                                                     method="min")[r.Index])}
            for r in tsub.itertuples()],
        "reading":
            "No compound in the 5,000 reaches a conformal p_active above "
            f"{float(pa.max()):.3f}; the median is {float(pa.median()):.3f}. The "
            "tier criterion cuts at 0.05, so the compounds it admits occupy a "
            f"band {float(flagged.conformal_p_active.min()):.3f}-"
            f"{float(flagged.conformal_p_active.max()):.3f} wide at the very top "
            "of a distribution that never leaves the bottom of its range. Two of "
            "the three Tier-A compounds sit two ranks above the cut. Reported "
            "beside this, the calibrated probability in the candidate table "
            "(0.66-0.71) is the more flattering of two numbers that point in "
            "opposite directions: 0.67 is a moderate probability in absolute "
            "terms and a 5th-percentile one relative to known actives.",
    }

    report["provenance"] = {
        "script": "conformal_diagnostics.py",
        "reuses": "conformal_qsar.py (mondrian_calibrate, mondrian_pvalues) -- "
                  "the same functions that produced the released prediction sets",
        "pipeline_sha256": sha256(cq.PIPELINE_PATH),
        "inputs": {
            "nplibrary_conformal_predictions.csv":
                sha256(NP_RUN / "nplibrary_conformal_predictions.csv"),
            "snplib_conformal_predictions.csv":
                sha256(QSAR_RUN / "snplib_conformal_predictions.csv"),
            "chembl_a2a_training_set.csv": sha256(cq.TRAINING_SET_CSV),
        },
        "cross_conformal": {"n_folds": N_FOLDS, "seed": SEED,
                            "grouping": "Bemis-Murcko scaffold"},
    }

    (out / "conformal_diagnostics.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")

    print("\nempty-set rate (the predictor abstaining):")
    for s in SIGS:
        row = {k: v["by_significance"][str(s)]["empty_{}"]
               for k, v in report["panels"].items()}
        print(f"  alpha={s:.2f}  " + "  ".join(f"{k}={v:.4f}" for k, v in row.items()))
    print("\nmean prediction-set size:")
    for s in SIGS:
        row = {k: v["by_significance"][str(s)]["mean_set_size"]
               for k, v in report["panels"].items()}
        print(f"  alpha={s:.2f}  " + "  ".join(f"{k}={v:.3f}" for k, v in row.items()))
    print("\np(inactive) median:  " + "  ".join(
        f"{k}={v['p_inactive_distribution']['median']:.4f}"
        for k, v in report["panels"].items()))
    print(f"\nwritten -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""ADR-AIIM-289, executed. Predictions frozen before this ran.

Pre-registration SHA-256
  ead611ad83c9bf5aab6d43d66d23902b5d0a8925cf41f50c249675befc0a9b7f
Decision log at freeze
  9b920ad1e2f0a4c0aea7beb7ffcd2bb12be4e04339ff53793f8173a6ed99fefa

ADR-AIIM-288 showed the Ki-only label (arm K) cuts convention sensitivity from
376x to 3.62x, and explicitly declared its discrimination UNTESTED. Arm K was
promoted to priority 2 on reproducibility evidence alone. This tests whether it
also discriminates better, or only more reproducibly.

`nested_cv_auc` is imported from the frozen tree and its fold structure and
hyperparameter space are used unchanged. The per-fold metrics beyond ROC-AUC
(MCC, balanced accuracy, PR-AUC and its lift over prevalence) are computed here
because the frozen function returns only mean and SD of ROC-AUC -- the fold
structure, search space, seeds and scoring are identical to it by construction,
and the reproduction of its ROC-AUC is the anchor (V3).

PR-AUC lift over prevalence is reported per ADR-AIIM-268's convention: at
76-88 % prevalence the raw PR-AUC flatters and the lift is the honest figure.

`10_new_system` is not modified.
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "3")

import numpy as np
import pandas as pd
from sklearn.metrics import (average_precision_score, balanced_accuracy_score,
                             matthews_corrcoef, roc_auc_score)
from sklearn.model_selection import (RandomizedSearchCV, StratifiedGroupKFold,
                                     StratifiedKFold)

ROOT = Path(r"D:\Agentic VS Workflow")
SCRIPTS = ROOT / "10_new_system" / "scripts"
REPULL = ROOT / "20_new_runs" / "chembl_censored_repull_2026-10-01" / "chembl_a2a_all_relations_v3.csv"
OUTDIR = ROOT / "20_new_runs" / "arm_k_auc_2026-10-05"
sys.path.insert(0, str(SCRIPTS))

PREREG = "ead611ad83c9bf5aab6d43d66d23902b5d0a8925cf41f50c249675befc0a9b7f"
LOG_AT_FREEZE = "9b920ad1e2f0a4c0aea7beb7ffcd2bb12be4e04339ff53793f8173a6ed99fefa"
RS, N_PERM = 42, 2
# ADR-AIIM-268's reported nested-CV ROC-AUC, the V3 anchor.
ANCHOR = {"A": 0.9181, "B": 0.9371}
ANCHOR_TOL = 0.02


def load(name, fname):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / fname)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def build_arms(qsar, cmp_mod):
    """Identical construction to ADR-AIIM-287 (label_balance_domain.py)."""
    raw = pd.read_csv(cmp_mod.CHEMBL_RAW_SOURCE)
    rep = pd.read_csv(REPULL).dropna(subset=["canonical_smiles"]).drop_duplicates()
    rep["ik"] = rep["canonical_smiles"].apply(qsar._smiles_to_inchikey)
    rep = rep.dropna(subset=["ik"])

    def recovery(have, ki_only):
        c = rep[rep["censoring_direction"] == "upper"].dropna(subset=["pchembl_upper_bound"])
        if ki_only:
            c = c[c["standard_type"] == "Ki"]
        cm = c[~c["ik"].isin(have)].groupby("ik").agg(
            b=("pchembl_upper_bound", "min"), smi=("canonical_smiles", "first"))
        return [r["smi"] for _, r in cm.iterrows() if r["b"] <= 5.5]

    out = {}
    for arm, src, ki in (("A", raw, None), ("B", raw, False),
                         ("K", raw[raw["standard_type"] == "Ki"], True)):
        d = cmp_mod.build_training_set(src, 6.0, 5.5, qsar)
        smi = d["canonical_smiles"].tolist()
        y = list(d["label"].to_numpy())
        if ki is None:
            out[arm] = (smi, np.array(y))
        else:
            add = recovery(set(d["inchikey"]), ki)
            out[arm] = (smi + add, np.array(y + [0] * len(add)))
        a = int(out[arm][1].sum())
        n = len(out[arm][1])
        print(f"  arm {arm}: {a} active / {n-a} inactive  n={n}  "
              f"prevalence {a/n*100:.1f} %", flush=True)
    return out


def nested_cv(X, y, groups, cmp_mod):
    """Fold structure, space, seeds and scoring exactly as
    qsar_threshold_comparison.nested_cv_auc; extra metrics added per fold."""
    space = {"rf__n_estimators": [100, 200, 500], "rf__max_depth": [None, 15, 25],
             "rf__min_samples_leaf": [1, 5, 10]}
    outer = StratifiedGroupKFold(n_splits=cmp_mod.N_OUTER_FOLDS, shuffle=True,
                                 random_state=cmp_mod.RANDOM_STATE)
    inner = StratifiedKFold(n_splits=cmp_mod.N_INNER_FOLDS, shuffle=True,
                            random_state=cmp_mod.RANDOM_STATE)
    folds = []
    for tr, val in outer.split(X, y, groups=groups):
        search = RandomizedSearchCV(cmp_mod.make_rf_pipe(), space,
                                    n_iter=cmp_mod.N_HP_ITER, cv=inner,
                                    scoring="roc_auc",
                                    random_state=cmp_mod.RANDOM_STATE, n_jobs=1)
        search.fit(X[tr], y[tr])
        est = search.best_estimator_
        p = est.predict_proba(X[val])[:, 1]
        yhat = est.predict(X[val])
        prev = float(y[val].mean())
        pr = float(average_precision_score(y[val], p))
        folds.append({
            "roc_auc": float(roc_auc_score(y[val], p)),
            "pr_auc": pr, "prevalence": prev, "pr_auc_lift": pr - prev,
            "mcc": float(matthews_corrcoef(y[val], yhat)),
            "balanced_accuracy": float(balanced_accuracy_score(y[val], yhat)),
            "n_val": int(len(val)), "best_params": {
                k: (v if v is not None else "None")
                for k, v in search.best_params_.items()}})
    agg = {}
    for m in ("roc_auc", "pr_auc", "pr_auc_lift", "mcc", "balanced_accuracy",
              "prevalence"):
        v = np.array([f[m] for f in folds])
        agg[m] = {"mean": round(float(v.mean()), 4), "sd": round(float(v.std()), 4)}
    return {"folds": folds, "aggregate": agg}


def main() -> int:
    (OUTDIR / "analysis").mkdir(parents=True, exist_ok=True)
    qsar = load("_qsar", "run_qsar_lbdd_2026-08-07.py")
    cmp_mod = load("_cmp", "qsar_threshold_comparison.py")
    print(f"frozen CV: outer={cmp_mod.N_OUTER_FOLDS} inner={cmp_mod.N_INNER_FOLDS} "
          f"n_iter={cmp_mod.N_HP_ITER} seed={cmp_mod.RANDOM_STATE}", flush=True)
    arms = build_arms(qsar, cmp_mod)

    out = {"adr": "ADR-AIIM-289", "prereg_sha256": PREREG,
           "log_sha256_at_freeze": LOG_AT_FREEZE,
           "adr_268_reported_roc_auc": ANCHOR,
           "adr_268_composition": {"A": "5244 active / 698 inactive",
                                   "B": "5244 active / 1429 inactive"},
           "note": "10_new_system not modified; nested_cv_auc's fold structure, "
                   "space, seeds and scoring used unchanged. Canonical InChIKey "
                   "row order throughout (ADR-AIIM-275).",
           "arms": {}}
    t0 = time.time()

    for arm, (smiles, y) in arms.items():
        X = qsar.compute_raw_features(smiles)
        ikeys = [qsar._smiles_to_inchikey(s) or f"__NONE_{i}__"
                 for i, s in enumerate(smiles)]
        groups_all = qsar.get_scaffold_groups(smiles)
        draws = []
        for k in range(N_PERM):
            idx = (np.arange(len(y)) if k == 0
                   else np.random.default_rng(10_000 + k).permutation(len(y)))
            canon = sorted(range(len(idx)), key=lambda j: ikeys[idx[j]])
            idx = idx[np.array(canon)]
            r = nested_cv(X[idx], y[idx], groups_all[idx], cmp_mod)
            draws.append(r)
            a = r["aggregate"]
            print(f"  arm {arm} draw {k}: ROC {a['roc_auc']['mean']:.4f} "
                  f"+/-{a['roc_auc']['sd']:.4f}  PR {a['pr_auc']['mean']:.4f} "
                  f"(lift {a['pr_auc_lift']['mean']:+.4f})  "
                  f"MCC {a['mcc']['mean']:.4f}  "
                  f"balacc {a['balanced_accuracy']['mean']:.4f}  "
                  f"[{time.time()-t0:.0f}s]", flush=True)
        roc = [d["aggregate"]["roc_auc"]["mean"] for d in draws]
        a0 = int(y.sum())
        out["arms"][arm] = {
            "composition": {"active": a0, "inactive": int(len(y) - a0),
                            "n": int(len(y)),
                            "prevalence": round(a0 / len(y), 4)},
            "draws": draws, "aggregate": draws[0]["aggregate"],
            "roc_auc_across_permutations": roc,
            "invariant_across_permutations": len(set(round(v, 6) for v in roc)) == 1}
        (OUTDIR / "arm_k_auc.json").write_text(json.dumps(out, indent=1),
                                               encoding="utf-8")

    # ---------------- grading ----------------
    g = lambda a, m, s="mean": out["arms"][a]["aggregate"][m][s]
    anchor_ok = {a: abs(g(a, "roc_auc") - ANCHOR[a]) <= ANCHOR_TOL
                 for a in ANCHOR}
    out["grading"] = {
        "V1_K_roc_exceeds_B": {
            "arm_K": g("K", "roc_auc"), "arm_B": g("B", "roc_auc"),
            "delta": round(g("K", "roc_auc") - g("B", "roc_auc"), 4),
            "verdict": "MET" if g("K", "roc_auc") > g("B", "roc_auc") else "MISSED"},
        "V2_K_sd_below_B": {
            "arm_K_sd": g("K", "roc_auc", "sd"), "arm_B_sd": g("B", "roc_auc", "sd"),
            "verdict": "MET" if g("K", "roc_auc", "sd") < g("B", "roc_auc", "sd")
            else "MISSED"},
        "V3_anchor_within_0.02": {
            "observed": {a: g(a, "roc_auc") for a in ANCHOR},
            "expected": ANCHOR, "within_tolerance": anchor_ok,
            "verdict": "MET" if all(anchor_ok.values()) else "MISSED"},
        "V4_K_pr_lift_exceeds_B": {
            "arm_K_lift": g("K", "pr_auc_lift"), "arm_B_lift": g("B", "pr_auc_lift"),
            "arm_K_pr_auc": g("K", "pr_auc"), "arm_B_pr_auc": g("B", "pr_auc"),
            "verdict": "MET" if g("K", "pr_auc_lift") > g("B", "pr_auc_lift")
            else "MISSED"},
    }
    out["seconds"] = round(time.time() - t0, 1)
    (OUTDIR / "arm_k_auc.json").write_text(json.dumps(out, indent=1), encoding="utf-8")

    print("\n" + "=" * 76)
    for t, gg in out["grading"].items():
        print(f"{t:28s} {gg['verdict']}")
    print("=" * 76)
    hdr = f"{'arm':4s} {'ROC-AUC':>16s} {'PR-AUC':>9s} {'prev':>7s} {'lift':>8s} {'MCC':>7s} {'bal-acc':>8s}"
    print(hdr)
    for a in out["arms"]:
        print(f"{a:4s} {g(a,'roc_auc'):8.4f} +/-{g(a,'roc_auc','sd'):.4f} "
              f"{g(a,'pr_auc'):9.4f} {g(a,'prevalence'):7.4f} "
              f"{g(a,'pr_auc_lift'):+8.4f} {g(a,'mcc'):7.4f} "
              f"{g(a,'balanced_accuracy'):8.4f}")
    print("\ninvariant across permutations: " + json.dumps(
        {a: out["arms"][a]["invariant_across_permutations"] for a in out["arms"]}))
    print(f"\nwritten -> {OUTDIR / 'arm_k_auc.json'} [{out['seconds']}s]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

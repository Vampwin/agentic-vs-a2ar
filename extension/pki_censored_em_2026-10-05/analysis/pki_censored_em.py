"""ADR-AIIM-298, executed. Predictions frozen before this ran.

Pre-registration SHA-256
  27c3752ba8105d0dd25bcca329bcf96690f71499c03a9c3513e93b45c79ee1d5
Decision log at freeze (text/LF, per the 2026-10-05 addendum)
  4014358c797ad0728cb0d36658931a537e52f334603b017cf967d1033e87a3f6

ADR-AIIM-297 held all 914 censored compounds out of training and tested the
model against their bounds. Not one passed, in any convention -- because that
design left a training set 86.2 % of which sits above pKi 6.0, and a forest on
such a set cannot reach down to bounds whose median is 5.00.

ONE VARIABLE MOVES HERE: how censored records enter. Same model family, same
features, same pipeline, same eight conventions, same canonical ordering. A
linear tobit would have changed the estimator too and is deliberately not used.

Censored-EM, the interval-respecting analogue of tobit for a tree model:
initialise each censored target at its bound, fit, predict, then PROJECT BACK
INTO THE INTERVAL the measurement actually asserts -- min(pred, bound) for
upper-censored, max(pred, bound) for lower -- and iterate. The target is never
allowed outside the interval, so a censored record contributes its information
without being assigned a value it does not have.

This is NOT tobit: no Gaussian assumption, no likelihood value, no standard
errors (ADR-AIIM-298 constraint 5).

Censored compounds are now split like everything else: those in the training
side contribute intervals, those in the held-out side are scored one-sided and
never fitted.

`10_new_system` is not modified; chembl251_raw.csv is read, not re-pulled.
"""
from __future__ import annotations

import importlib.util
import itertools
import json
import os
import random
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "3")

import numpy as np
import pandas as pd
import rdkit
from rdkit import RDLogger
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor
from sklearn.feature_selection import VarianceThreshold
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.pipeline import Pipeline

RDLogger.DisableLog("rdApp.*")

ROOT = Path(r"D:\Agentic VS Workflow")
SCRIPTS = ROOT / "10_new_system" / "scripts"
CUR_RUN = ROOT / "20_new_runs" / "chembl251_full_curation_2026-10-05"
REG_RUN = ROOT / "20_new_runs" / "pki_regression_2026-10-05"
NPLIB = ROOT / "20_new_runs" / "np_library_funnel_2026-08-15" / "np_library_dock_ledger.csv"
SNPLIB = ROOT / "20_new_runs" / "snp_screen_2026-08-04" / "G1_standardization" / "snp_source_ledger.csv"
OUTDIR = ROOT / "20_new_runs" / "pki_censored_em_2026-10-05"
sys.path.insert(0, str(SCRIPTS))

PREREG = "27c3752ba8105d0dd25bcca329bcf96690f71499c03a9c3513e93b45c79ee1d5"
LOG_AT_FREEZE = "4014358c797ad0728cb0d36658931a537e52f334603b017cf967d1033e87a3f6"
RS, TEST_FRAC, SIG = 42, 0.20, 0.05
THRESHOLDS = [5.0, 5.5, 6.0, 6.5, 7.0]
EM_MAX_ITER, EM_TOL = 10, 0.01
ADR297_SPREAD = 10.99          # ADR-AIIM-297, same compounds
ADR297_NP_MIN_PRED = 5.144     # ADR-AIIM-297 anchor, the T2 reference
CLASSIFIER_SPREAD = 38.88      # ADR-AIIM-293


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


CELLS = [{"d1_shuffle": a, "d2_train_first": b, "d3_divert": c}
         for a, b, c in itertools.product((True, False), (False, True), (False, True))]


def cell_name(c):
    return (f"D1{'shuffle' if c['d1_shuffle'] else 'sorted'}"
            f"_D2{'train' if c['d2_train_first'] else 'test'}first"
            f"_D3{'divert' if c['d3_divert'] else 'overflow'}")


def split(scaffolds, n_rows, d1_shuffle, d2_train_first, d3_divert):
    s2i: dict[str, list[int]] = {}
    for i, sc in enumerate(scaffolds):
        s2i.setdefault(sc, []).append(i)
    ordered = sorted(s2i.items(), key=lambda x: (-len(x[1]), x[0]))
    if d1_shuffle:
        random.Random(RS).shuffle(ordered)
    n_te = int(n_rows * TEST_FRAC)
    n_tr = n_rows - n_te
    te: set[int] = set()
    tr: set[int] = set()
    for _sc, idxs in ordered:
        primary, budget, other = ((tr, n_tr, te) if d2_train_first
                                  else (te, n_te, tr))
        if d3_divert:
            (primary if len(primary) + len(idxs) <= budget else other).update(idxs)
        else:
            (primary if len(primary) < budget else other).update(idxs)
    return np.array(sorted(tr)), np.array(sorted(te))


def pipe():
    return Pipeline([("var", VarianceThreshold(threshold=0.01)),
                     ("rf", RandomForestRegressor(random_state=RS, n_jobs=3))])


def censored_em(X, y0, kind, tr, max_iter=EM_MAX_ITER, tol=EM_TOL):
    """kind: 0 point, +1 upper-censored (y <= bound), -1 lower (y >= bound).
    Returns (fitted pipeline, n_iter, final_max_change, imputed targets)."""
    y = y0.copy()
    up = tr[kind[tr] == 1]
    lo = tr[kind[tr] == -1]
    n_iter, change = 0, float("nan")
    mdl = None
    for it in range(1, max_iter + 1):
        mdl = pipe().fit(X[tr], y[tr])
        n_iter = it
        if len(up) == 0 and len(lo) == 0:
            change = 0.0
            break
        new = y.copy()
        if len(up):
            new[up] = np.minimum(mdl.predict(X[up]), y0[up])
        if len(lo):
            new[lo] = np.maximum(mdl.predict(X[lo]), y0[lo])
        change = float(np.max(np.abs(new[tr] - y[tr]))) if len(tr) else 0.0
        y = new
        if change < tol:
            mdl = pipe().fit(X[tr], y[tr])
            break
    return mdl, n_iter, change, y


def main() -> int:
    (OUTDIR / "analysis").mkdir(parents=True, exist_ok=True)
    qsar = load("_qsar", SCRIPTS / "run_qsar_lbdd_2026-08-07.py")
    cur = load("_cur", CUR_RUN / "analysis" / "chembl251_curate.py")
    reg = load("_reg", REG_RUN / "analysis" / "pki_regression.py")
    t0 = time.time()

    tbl, n_raw = reg.rebuild_table(cur)
    b = tbl[tbl["endpoint_class"] == "binding"].reset_index(drop=True)
    keep = b["basis"].isin(["point_median", "point_median_gap",
                            "censored_upper", "censored_lower"])
    b = b[keep & np.isfinite(b["px"])].reset_index(drop=True)
    kind = np.where(b["basis"] == "censored_upper", 1,
                    np.where(b["basis"] == "censored_lower", -1, 0)).astype(int)
    smiles = b["std_smiles"].tolist()
    y0 = b["px"].to_numpy(dtype=float)

    ikeys = [qsar._smiles_to_inchikey(s) or f"__N{i}__" for i, s in enumerate(smiles)]
    order = sorted(range(len(smiles)), key=lambda j: ikeys[j])
    smiles = [smiles[i] for i in order]
    y0, kind = y0[order], kind[order]
    scaf = [qsar._murcko(s) for s in smiles]
    X = qsar.compute_raw_features(smiles)
    print(f"raw {n_raw} -> modelling set {len(y0)}: "
          f"point {(kind==0).sum()}  upper {(kind==1).sum()}  "
          f"lower {(kind==-1).sum()}", flush=True)
    ub = y0[kind == 1]
    print(f"  upper bounds: min {ub.min():.2f} median {np.median(ub):.2f} "
          f"max {ub.max():.2f}", flush=True)

    q_smiles, queries = {}, {}
    for lname, path, col in (("SNP-Lib", SNPLIB, "smiles_canonical"),
                             ("NP-library", NPLIB, "SMILES")):
        q = pd.read_csv(path)
        q = q[q[col].apply(lambda s: isinstance(s, str) and bool(s))]
        q_smiles[lname] = q[col].tolist()
        queries[lname] = qsar.compute_raw_features(q_smiles[lname])

    out = {"adr": "ADR-AIIM-298", "prereg_sha256": PREREG,
           "log_sha256_at_freeze": LOG_AT_FREEZE,
           "rdkit_version": rdkit.__version__,
           "method": "censored-EM (interval-projected target), NOT tobit",
           "em": {"max_iter": EM_MAX_ITER, "tol": EM_TOL},
           "composition": {
               "modelling_set": int(len(y0)), "point": int((kind == 0).sum()),
               "censored_upper": int((kind == 1).sum()),
               "censored_lower": int((kind == -1).sum()),
               "upper_bound_median": round(float(np.median(ub)), 3)},
           "thresholds": THRESHOLDS,
           "reference": {"adr297_spread": ADR297_SPREAD,
                         "adr297_np_min_pred": ADR297_NP_MIN_PRED,
                         "classifier_spread": CLASSIFIER_SPREAD},
           "conventions": {}}

    for cfg in CELLS:
        nm = cell_name(cfg)
        tr, te = split(scaf, len(y0), **cfg)
        mdl, n_it, chg, y_imp = censored_em(X, y0, kind, tr)

        te_pt = te[kind[te] == 0]
        te_up = te[kind[te] == 1]
        te_lo = te[kind[te] == -1]
        p_pt = mdl.predict(X[te_pt])
        resid = np.abs(y0[te_pt] - p_pt)
        k = int(np.ceil((1 - SIG) * (len(resid) + 1))) - 1
        qhat = float(np.sort(resid)[min(k, len(resid) - 1)])

        rec = {"em_iterations": n_it, "em_final_max_change": round(chg, 5),
               "n_train": int(len(tr)), "n_train_censored": int((kind[tr] != 0).sum()),
               "n_heldout_point": int(len(te_pt)),
               "n_heldout_upper": int(len(te_up)),
               "n_heldout_lower": int(len(te_lo)),
               "r2": round(float(r2_score(y0[te_pt], p_pt)), 4),
               "rmse": round(float(np.sqrt(mean_squared_error(y0[te_pt], p_pt))), 4),
               "mae": round(float(mean_absolute_error(y0[te_pt], p_pt)), 4),
               "spearman": round(float(spearmanr(y0[te_pt], p_pt).statistic), 4),
               "conformal_qhat": round(qhat, 4)}

        pu = mdl.predict(X[te_up]) if len(te_up) else np.array([])
        pl = mdl.predict(X[te_lo]) if len(te_lo) else np.array([])
        rec["censored_heldout_test"] = {
            "upper_n": int(len(te_up)),
            "upper_pass": int((pu <= y0[te_up]).sum()) if len(te_up) else None,
            "upper_pass_frac": round(float((pu <= y0[te_up]).mean()), 4) if len(te_up) else None,
            "upper_mean_violation": round(float(np.clip(pu - y0[te_up], 0, None).mean()), 4) if len(te_up) else None,
            "lower_n": int(len(te_lo)),
            "lower_pass_frac": round(float((pl >= y0[te_lo]).mean()), 4) if len(te_lo) else None}

        rec["counts"] = {}
        for lname, Xq in queries.items():
            pq = mdl.predict(Xq)
            rec["counts"][lname] = {
                "point": {str(t): int((pq >= t).sum()) for t in THRESHOLDS},
                "conformal_lower": {str(t): int((pq - qhat >= t).sum())
                                    for t in THRESHOLDS},
                "pred_min": round(float(pq.min()), 3),
                "pred_max": round(float(pq.max()), 3),
                "pred_mean": round(float(pq.mean()), 3)}
        out["conventions"][nm] = rec
        ct = rec["censored_heldout_test"]
        print(f"  {nm:40s} EM {n_it}it  rho {rec['spearman']:.3f}  "
              f"R2 {rec['r2']:.3f}  up-pass {ct['upper_pass_frac']}  "
              f"NPmin {rec['counts']['NP-library']['pred_min']:.3f}  "
              f"NP@6 {rec['counts']['NP-library']['point']['6.0']}  "
              f"conf@6 {rec['counts']['NP-library']['conformal_lower']['6.0']}  "
              f"[{time.time()-t0:.0f}s]", flush=True)

    anchor = cell_name({"d1_shuffle": False, "d2_train_first": True,
                        "d3_divert": True})
    A = out["conventions"][anchor]
    np6 = [out["conventions"][cell_name(c)]["counts"]["NP-library"]["point"]["6.0"]
           for c in CELLS]
    spread = max(np6) / max(min(np6), 1)
    conf6 = {cell_name(c): out["conventions"][cell_name(c)]["counts"]
             ["NP-library"]["conformal_lower"]["6.0"] for c in CELLS}

    out["grading"] = {
        "T1_upper_censored_pass_above_70pct": {
            "anchor": anchor,
            "pass_frac": A["censored_heldout_test"]["upper_pass_frac"],
            "n": A["censored_heldout_test"]["upper_n"],
            "adr297_value": 0.0,
            "all_conventions": {cell_name(c): out["conventions"][cell_name(c)]
                                ["censored_heldout_test"]["upper_pass_frac"]
                                for c in CELLS},
            "verdict": "MET" if (A["censored_heldout_test"]["upper_pass_frac"]
                                 or 0) > 0.70 else "MISSED"},
        "T2_np_min_pred_below_5.0": {
            "value": A["counts"]["NP-library"]["pred_min"],
            "adr297_value": ADR297_NP_MIN_PRED,
            "verdict": "MET" if A["counts"]["NP-library"]["pred_min"] < 5.0
            else "MISSED"},
        "T3_spread_below_adr297": {
            "np_counts_at_6.0": {cell_name(c): n for c, n in zip(CELLS, np6)},
            "min": min(np6), "max": max(np6), "fold_range": round(spread, 2),
            "adr297_fold_range": ADR297_SPREAD,
            "classifier_fold_range": CLASSIFIER_SPREAD,
            "verdict": "MET" if spread < ADR297_SPREAD else "MISSED"},
        "T4_conformal_count_still_zero": {
            "conformal_lower_at_6.0": conf6,
            "qhat_range": [min(out["conventions"][cell_name(c)]["conformal_qhat"]
                               for c in CELLS),
                           max(out["conventions"][cell_name(c)]["conformal_qhat"]
                               for c in CELLS)],
            "verdict": "MET" if all(v == 0 for v in conf6.values()) else "MISSED"},
    }
    out["seconds"] = round(time.time() - t0, 1)
    (OUTDIR / "pki_censored_em.json").write_text(json.dumps(out, indent=1),
                                                 encoding="utf-8")

    print("\n" + "=" * 80)
    for t, g in out["grading"].items():
        print(f"{t:40s} {g['verdict']}")
    print("=" * 80)
    print(f"{'convention':40s} {'rho':>6s} {'R2':>7s} {'up-pass':>8s} "
          f"{'NPmin':>6s} {'NP@6':>6s} {'conf@6':>7s} {'qhat':>6s}")
    for c in CELLS:
        nm = cell_name(c)
        r = out["conventions"][nm]
        print(f"{nm:40s} {r['spearman']:6.3f} {r['r2']:7.3f} "
              f"{r['censored_heldout_test']['upper_pass_frac']:8.3f} "
              f"{r['counts']['NP-library']['pred_min']:6.3f} "
              f"{r['counts']['NP-library']['point']['6.0']:6d} "
              f"{r['counts']['NP-library']['conformal_lower']['6.0']:7d} "
              f"{r['conformal_qhat']:6.3f}")
    print(f"\nNP@6.0 spread {min(np6)}-{max(np6)} = {spread:.2f}x   "
          f"(ADR-297 {ADR297_SPREAD}x, classifier {CLASSIFIER_SPREAD}x)")
    print("threshold curve, anchor, NP-library:")
    for t in THRESHOLDS:
        print(f"  >= {t}: point {A['counts']['NP-library']['point'][str(t)]:5d}   "
              f"conformal lower {A['counts']['NP-library']['conformal_lower'][str(t)]:5d}")
    print(f"\nwritten -> {OUTDIR / 'pki_censored_em.json'} [{out['seconds']}s]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

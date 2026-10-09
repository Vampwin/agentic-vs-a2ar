"""ADR-AIIM-295, executed. Band edges, metric and criterion frozen before this ran.

Pre-registration SHA-256
  30555270484983ce475618359c4f8a6c71859086266a97888cf40931bedf291e
Decision log at freeze (text/LF, per the 2026-10-05 addendum)
  a61c4017b61354af99084fd790bcc22ada08104d6e3608f5ad17a0cf74061c53

ADR-AIIM-294 measured two points on the distance axis -- the whole held-out set
(0.9217) and the natural-product band (0.6808) -- and nothing between them.
This measures the curve.

NO REFITTING. The model is reconstructed from the same seeded, canonical inputs
as ADR-AIIM-294, and the full held-out ROC-AUC of 0.9217 is the anchor: if it
does not reproduce, the curve is measured on a different model and nothing here
holds.

A1 (ECFP4 Tanimoto) and A8 (ECFP4 Dice) only. A2-A7 are excluded because their
bands covered 68-94 % of the held-out set -- natural products are not at a
distinguishable distance in those spaces -- and that exclusion is declared in
the pre-registration, not chosen after seeing results.

The 0.15-0.25 band holds 22 compounds. It is REPORTED AND NOT DECISION-BEARING,
declared in advance so that an inconvenient estimate cannot be dropped later.

`10_new_system` is not modified.
"""
from __future__ import annotations

import importlib.util
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
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import VarianceThreshold
from sklearn.metrics import matthews_corrcoef, roc_auc_score
from sklearn.pipeline import Pipeline

RDLogger.DisableLog("rdApp.*")

ROOT = Path(r"D:\Agentic VS Workflow")
SCRIPTS = ROOT / "10_new_system" / "scripts"
CURATED = ROOT / "20_new_runs" / "chembl251_full_curation_2026-10-05" / "curated_binding.csv"
OUTDIR = ROOT / "20_new_runs" / "ad_distance_curve_2026-10-05"
sys.path.insert(0, str(SCRIPTS))

PREREG = "30555270484983ce475618359c4f8a6c71859086266a97888cf40931bedf291e"
LOG_AT_FREEZE = "a61c4017b61354af99084fd790bcc22ada08104d6e3608f5ad17a0cf74061c53"
RS, TEST_FRAC, AD_K = 42, 0.20, 5
N_BOOT, TIER2_CUT = 2000, 0.70
ANCHOR_AUC = 0.9217                       # ADR-AIIM-294, must reproduce
BANDS = [(0.15, 0.25), (0.25, 0.35), (0.35, 0.45),
         (0.45, 0.55), (0.55, 0.70), (0.70, 1.01)]
NOT_DECISION_BEARING = (0.15, 0.25)       # declared in ADR-AIIM-295

_morgan = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def load(name, fname):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / fname)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def fps(smiles):
    out = []
    for s in smiles:
        mol = Chem.MolFromSmiles(s) if isinstance(s, str) else None
        out.append(None if mol is None else _morgan.GetFingerprint(mol))
    return out


def knn_scores(query_fps, train_fps, metric, k=AD_K):
    bulk = (DataStructs.BulkTanimotoSimilarity if metric == "tanimoto"
            else DataStructs.BulkDiceSimilarity)
    out = np.full(len(query_fps), np.nan)
    for i, fp in enumerate(query_fps):
        if fp is not None:
            out[i] = float(np.mean(sorted(bulk(fp, train_fps), reverse=True)[:k]))
    return out


def boot_auc(y, p, n_boot=N_BOOT, seed=RS):
    if len(y) == 0 or len(np.unique(y)) < 2:
        return None, None, None
    point = float(roc_auc_score(y, p))
    rng = np.random.default_rng(seed)
    v = []
    for _ in range(n_boot):
        ix = rng.integers(0, len(y), len(y))
        if len(np.unique(y[ix])) > 1:
            v.append(roc_auc_score(y[ix], p[ix]))
    if not v:
        return point, None, None
    lo, hi = np.percentile(v, [2.5, 97.5])
    return point, float(lo), float(hi)


def tier(point, lo):
    if point is None or lo is None:
        return "undetermined"
    if lo <= 0.50:
        return "Tier 0"
    return "Tier 2" if point > TIER2_CUT else "Tier 1"


def main() -> int:
    (OUTDIR / "analysis").mkdir(parents=True, exist_ok=True)
    qsar = load("_qsar", "run_qsar_lbdd_2026-08-07.py")
    t0 = time.time()

    b = pd.read_csv(CURATED)
    smiles = b["std_smiles"].tolist()
    y = b["label"].to_numpy()
    ikeys = [qsar._smiles_to_inchikey(s) or f"__N{i}__" for i, s in enumerate(smiles)]
    order = sorted(range(len(smiles)), key=lambda j: ikeys[j])
    smiles = [smiles[i] for i in order]
    y = y[order]
    scaf = [qsar._murcko(s) for s in smiles]

    s2i: dict[str, list[int]] = {}
    for i, sc in enumerate(scaf):
        s2i.setdefault(sc, []).append(i)
    ordered = sorted(s2i.items(), key=lambda x: (-len(x[1]), x[0]))
    n_tr = len(smiles) - int(len(smiles) * TEST_FRAC)
    tr_s: set[int] = set()
    te_s: set[int] = set()
    for _sc, ix in ordered:
        (tr_s if len(tr_s) + len(ix) <= n_tr else te_s).update(ix)
    tri, tei = sorted(tr_s), sorted(te_s)

    X = qsar.compute_raw_features(smiles)
    mdl = Pipeline([("var", VarianceThreshold(threshold=0.01)),
                    ("rf", RandomForestClassifier(class_weight="balanced",
                                                  random_state=RS, n_jobs=3))])
    mdl.fit(X[tri], y[tri])
    p_he = mdl.predict_proba(X[tei])[:, 1]
    yhat = mdl.predict(X[tei])
    y_he = y[tei]
    full_auc = float(roc_auc_score(y_he, p_he))
    anchor_ok = abs(full_auc - ANCHOR_AUC) < 5e-4
    print(f"train {len(tri)}  held-out {len(tei)} ({int(y_he.sum())} act)", flush=True)
    print(f"ANCHOR full held-out AUC {full_auc:.4f}  "
          f"(ADR-294: {ANCHOR_AUC})  reproduces: {anchor_ok}", flush=True)
    if not anchor_ok:
        print("  !! anchor failed -- the curve below is on a DIFFERENT model", flush=True)

    train_fps = [f for f in fps([smiles[i] for i in tri]) if f is not None]
    held_fps = fps([smiles[i] for i in tei])

    out = {"adr": "ADR-AIIM-295", "prereg_sha256": PREREG,
           "log_sha256_at_freeze": LOG_AT_FREEZE,
           "rdkit_version": rdkit.__version__,
           "anchor": {"observed": round(full_auc, 4), "expected": ANCHOR_AUC,
                      "reproduces": bool(anchor_ok)},
           "tier_criterion": {"tier0": "CI lower <= 0.50",
                              "tier1": "CI lower > 0.50 and point <= 0.70",
                              "tier2": "CI lower > 0.50 and point > 0.70",
                              "n_bootstrap": N_BOOT},
           "not_decision_bearing": list(NOT_DECISION_BEARING),
           "definitions": {}}

    for dname, metric in (("A1_ecfp4_tanimoto", "tanimoto"),
                          ("A8_ecfp4_dice", "dice")):
        s_he = knn_scores(held_fps, train_fps, metric)
        rows, dec_auc, dec_names = [], [], []
        for lo_e, hi_e in BANDS:
            m = np.isfinite(s_he) & (s_he >= lo_e) & (s_he < hi_e)
            n = int(m.sum())
            yb, pb = y_he[m], p_he[m]
            pt, lo, hi = boot_auc(yb, pb)
            mcc = (float(matthews_corrcoef(yb, yhat[m]))
                   if n and len(np.unique(yb)) > 1 else None)
            dec = (lo_e, hi_e) != NOT_DECISION_BEARING
            r = {"band": [lo_e, hi_e], "n": n,
                 "active": int(yb.sum()) if n else 0,
                 "inactive": int((yb == 0).sum()) if n else 0,
                 "roc_auc": None if pt is None else round(pt, 4),
                 "ci95": None if lo is None else [round(lo, 4), round(hi, 4)],
                 "mcc": None if mcc is None else round(mcc, 4),
                 "tier": tier(pt, lo), "decision_bearing": dec}
            rows.append(r)
            if dec and pt is not None:
                dec_auc.append(pt)
                dec_names.append(f"{lo_e:.2f}-{hi_e:.2f}")
            flag = "" if dec else "   (reported, not decision-bearing)"
            print(f"  {dname:18s} {lo_e:.2f}-{hi_e:.2f}  n={n:4d} "
                  f"({r['active']}/{r['inactive']})  AUC {r['roc_auc']}  "
                  f"CI {r['ci95']}  {r['tier']}{flag}", flush=True)

        # Y1: Spearman over the decision-bearing bands, bootstrapped jointly
        ranks = np.arange(len(dec_auc))
        rho = float(spearmanr(ranks, dec_auc).statistic) if len(dec_auc) > 2 else None
        rng = np.random.default_rng(RS)
        rhos = []
        for _ in range(N_BOOT):
            vals = []
            for (lo_e, hi_e) in BANDS:
                if (lo_e, hi_e) == NOT_DECISION_BEARING:
                    continue
                m = np.isfinite(s_he) & (s_he >= lo_e) & (s_he < hi_e)
                yb, pb = y_he[m], p_he[m]
                if len(yb) == 0:
                    continue
                ix = rng.integers(0, len(yb), len(yb))
                if len(np.unique(yb[ix])) < 2:
                    vals.append(np.nan)
                else:
                    vals.append(roc_auc_score(yb[ix], pb[ix]))
            vals = np.asarray(vals, dtype=float)
            if np.isnan(vals).any() or len(vals) < 3:
                continue
            rhos.append(spearmanr(np.arange(len(vals)), vals).statistic)
        rho_ci = ([round(float(np.percentile(rhos, 2.5)), 4),
                   round(float(np.percentile(rhos, 97.5)), 4)] if rhos else None)

        tier2 = [r for r in rows if r["decision_bearing"] and r["tier"] == "Tier 2"]
        crossing = min((r["band"][0] for r in tier2), default=None)
        out["definitions"][dname] = {
            "bands": rows, "decision_bearing_bands": dec_names,
            "spearman_rank_vs_auc": None if rho is None else round(rho, 4),
            "spearman_ci95": rho_ci,
            "lowest_tier2_band_start": crossing,
            "any_tier2": bool(tier2),
            "any_below_tier2": any(r["decision_bearing"] and r["tier"] != "Tier 2"
                                   for r in rows)}
        print(f"  {dname}: spearman {out['definitions'][dname]['spearman_rank_vs_auc']} "
              f"CI {rho_ci}  crossing at {crossing}", flush=True)

    A1 = out["definitions"]["A1_ecfp4_tanimoto"]
    A8 = out["definitions"]["A8_ecfp4_dice"]
    out["grading"] = {
        "Y1_monotone_degradation": {
            "spearman_A1": A1["spearman_rank_vs_auc"], "ci_A1": A1["spearman_ci95"],
            "spearman_A8": A8["spearman_rank_vs_auc"], "ci_A8": A8["spearman_ci95"],
            "verdict": "MET" if (A1["spearman_ci95"] and A1["spearman_rank_vs_auc"]
                                 and A1["spearman_rank_vs_auc"] < 0
                                 and A1["spearman_ci95"][1] < 0) else "MISSED"},
        "Y2_crossing_inside_range": {
            "any_tier2": A1["any_tier2"], "any_below": A1["any_below_tier2"],
            "verdict": "MET" if (A1["any_tier2"] and A1["any_below_tier2"])
            else "MISSED"},
        "Y3_crossing_between_0.35_and_0.55": {
            "lowest_tier2_band_start": A1["lowest_tier2_band_start"],
            "verdict": "MET" if A1["lowest_tier2_band_start"] in (0.35, 0.45)
            else "MISSED"},
        "Y4_A1_and_A8_agree": {
            "A1": A1["lowest_tier2_band_start"], "A8": A8["lowest_tier2_band_start"],
            "verdict": "MET" if A1["lowest_tier2_band_start"] ==
            A8["lowest_tier2_band_start"] else "MISSED"},
    }
    out["seconds"] = round(time.time() - t0, 1)
    (OUTDIR / "ad_distance_curve.json").write_text(json.dumps(out, indent=1),
                                                   encoding="utf-8")

    print("\n" + "=" * 70)
    print(f"{'anchor reproduces ADR-294':38s} {out['anchor']['reproduces']}")
    for t, g in out["grading"].items():
        print(f"{t:38s} {g['verdict']}")
    print("=" * 70)
    for dn, D in out["definitions"].items():
        print(f"\n{dn}:")
        for r in D["bands"]:
            mark = "" if r["decision_bearing"] else "  *"
            print(f"  {r['band'][0]:.2f}-{r['band'][1]:.2f}  n={r['n']:4d}  "
                  f"AUC {str(r['roc_auc']):>6s}  CI {str(r['ci95']):>18s}  "
                  f"{r['tier']}{mark}")
        print(f"  crossing at AD {D['lowest_tier2_band_start']}  "
              f"spearman {D['spearman_rank_vs_auc']} CI {D['spearman_ci95']}")
    print("\n  * reported, not decision-bearing (declared in ADR-AIIM-295)")
    print(f"\nwritten -> {OUTDIR / 'ad_distance_curve.json'} [{out['seconds']}s]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""ADR-AIIM-292, executed. Predictions and the tier criterion frozen before
this ran.

Pre-registration SHA-256
  56e831514329e9f67208f541009f3604c49ceb1a782507d96fce243da611bd26
Decision log at freeze
  a753eb24a9862fecfe405f084fe9e9a27ff88bd692440a64510f1bba053e7a09

Figure 4 varies the AD SPACE and holds the METHOD fixed. This adds the method
dimension (leverage, distance-to-centroid, Mahalanobis) and -- the part absent
from the manuscript entirely -- measures PERFORMANCE as a function of position
on the distance axis.

ONE model is fitted. The eight AD definitions re-slice the same held-out
predictions, so no difference between them can come from refitting.

Uniform threshold convention, declared: every definition's threshold admits
90 % of the training set, which is what "10th percentile of training
self-similarity" does. Leverage is the exception by pre-registration: it uses
the classical Williams h* = 3p/n, and the fraction of training it admits is
reported for comparability.

Every score is oriented so that HIGHER = MORE IN-DOMAIN.

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
from rdkit import DataStructs, RDLogger
from rdkit.Chem import Descriptors, MACCSkeys, rdFingerprintGenerator, rdMolDescriptors
from rdkit import Chem
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import VarianceThreshold
from sklearn.metrics import matthews_corrcoef, roc_auc_score
from sklearn.pipeline import Pipeline

RDLogger.DisableLog("rdApp.*")

ROOT = Path(r"D:\Agentic VS Workflow")
SCRIPTS = ROOT / "10_new_system" / "scripts"
CURATED = ROOT / "20_new_runs" / "chembl251_full_curation_2026-10-05" / "curated_binding.csv"
NPLIB = ROOT / "20_new_runs" / "np_library_funnel_2026-08-15" / "np_library_dock_ledger.csv"
SNPLIB = ROOT / "20_new_runs" / "snp_screen_2026-08-04" / "G1_standardization" / "snp_source_ledger.csv"
OUTDIR = ROOT / "20_new_runs" / "ad_methods_distance_performance_2026-10-05"
sys.path.insert(0, str(SCRIPTS))

PREREG = "56e831514329e9f67208f541009f3604c49ceb1a782507d96fce243da611bd26"
LOG_AT_FREEZE = "a753eb24a9862fecfe405f084fe9e9a27ff88bd692440a64510f1bba053e7a09"
RS, TEST_FRAC, AD_K, AD_PCT = 42, 0.20, 5, 10
N_BOOT, SELFSIM_SAMPLE = 2000, 500
FULL_HELDOUT_AUC = 0.9332          # ADR-AIIM-291, the X1 reference
TIER2_CUT = 0.70                   # declared QSAR convention


def load(name, fname):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / fname)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


# ------------------------------------------------------------- representations
_morgan = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
_atompair = rdFingerprintGenerator.GetAtomPairGenerator(fpSize=2048)
DESC_NAMES = ["mol_wt", "logp", "hbd", "hba", "tpsa", "rot", "rings", "arom"]


def _desc(mol):
    return [Descriptors.MolWt(mol), Descriptors.MolLogP(mol),
            rdMolDescriptors.CalcNumHBD(mol), rdMolDescriptors.CalcNumHBA(mol),
            rdMolDescriptors.CalcTPSA(mol),
            rdMolDescriptors.CalcNumRotatableBonds(mol),
            rdMolDescriptors.CalcNumRings(mol),
            rdMolDescriptors.CalcNumAromaticRings(mol)]


def represent(smiles):
    """Returns dict of representation -> list/array aligned to smiles, with a
    validity mask. Invalid molecules are carried as None / NaN."""
    out = {"ecfp4": [], "maccs": [], "atompair": [], "desc": []}
    ok = []
    for s in smiles:
        mol = Chem.MolFromSmiles(s) if isinstance(s, str) else None
        if mol is None:
            ok.append(False)
            out["ecfp4"].append(None)
            out["maccs"].append(None)
            out["atompair"].append(None)
            out["desc"].append([np.nan] * len(DESC_NAMES))
            continue
        ok.append(True)
        out["ecfp4"].append(_morgan.GetFingerprint(mol))
        out["maccs"].append(MACCSkeys.GenMACCSKeys(mol))
        out["atompair"].append(_atompair.GetFingerprint(mol))
        out["desc"].append(_desc(mol))
    out["desc"] = np.asarray(out["desc"], dtype=float)
    out["ok"] = np.asarray(ok)
    return out


# ------------------------------------------------------------------ AD scoring
def knn_fp_scores(query_fps, train_fps, metric="tanimoto", k=AD_K):
    """Mean of the top-k similarities to the training set. Higher = in-domain."""
    bulk = (DataStructs.BulkTanimotoSimilarity if metric == "tanimoto"
            else DataStructs.BulkDiceSimilarity)
    out = np.full(len(query_fps), np.nan)
    for i, fp in enumerate(query_fps):
        if fp is None:
            continue
        sims = bulk(fp, train_fps)
        out[i] = float(np.mean(sorted(sims, reverse=True)[:k]))
    return out


def knn_fp_selfscores(train_fps, metric="tanimoto", k=AD_K, n=SELFSIM_SAMPLE):
    """Training self-similarity: each sampled train compound's mean top-k to the
    OTHER train compounds. Matches run_qsar_lbdd's compute_ad_threshold."""
    bulk = (DataStructs.BulkTanimotoSimilarity if metric == "tanimoto"
            else DataStructs.BulkDiceSimilarity)
    idx = random.Random(RS).sample(range(len(train_fps)), min(n, len(train_fps)))
    out = []
    for i in idx:
        others = [f for j, f in enumerate(train_fps) if j != i]
        sims = bulk(train_fps[i], others)
        out.append(float(np.mean(sorted(sims, reverse=True)[:k])))
    return np.asarray(out)


def knn_euclid_scores(Q, T, k=AD_K):
    """Negative mean of the k smallest Euclidean distances. Higher = closer."""
    out = np.full(len(Q), np.nan)
    for i in range(len(Q)):
        if not np.all(np.isfinite(Q[i])):
            continue
        d = np.sqrt(((T - Q[i]) ** 2).sum(axis=1))
        out[i] = -float(np.mean(np.sort(d)[:k]))
    return out


def bootstrap_auc(y, p, n_boot=N_BOOT, seed=RS):
    if len(np.unique(y)) < 2:
        return None, None, None
    point = float(roc_auc_score(y, p))
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_boot):
        ix = rng.integers(0, len(y), len(y))
        if len(np.unique(y[ix])) < 2:
            continue
        vals.append(roc_auc_score(y[ix], p[ix]))
    if not vals:
        return point, None, None
    lo, hi = np.percentile(vals, [2.5, 97.5])
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

    # principled convention: sorted, train-first, divert (ADR-AIIM-286)
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
    train_smiles = [smiles[i] for i in tri]
    held_smiles = [smiles[i] for i in tei]
    y_tr, y_he = y[tri], y[tei]
    print(f"curated {len(smiles)}  train {len(tri)} ({int(y_tr.sum())} act)  "
          f"held-out {len(tei)} ({int(y_he.sum())} act)", flush=True)

    # ---- ONE model, fitted once ----
    X = qsar.compute_raw_features(smiles)
    mdl = Pipeline([("var", VarianceThreshold(threshold=0.01)),
                    ("rf", RandomForestClassifier(class_weight="balanced",
                                                  random_state=RS, n_jobs=3))])
    mdl.fit(X[tri], y_tr)
    p_he = mdl.predict_proba(X[tei])[:, 1]
    yhat_he = mdl.predict(X[tei])
    full_auc = float(roc_auc_score(y_he, p_he))
    print(f"full held-out ROC-AUC {full_auc:.4f} "
          f"(ADR-291 reported {FULL_HELDOUT_AUC})", flush=True)

    libs = {}
    for lname, path, col in (("SNP-Lib", SNPLIB, "smiles_canonical"),
                             ("NP-library", NPLIB, "SMILES")):
        q = pd.read_csv(path)
        q = q[q[col].apply(lambda s: isinstance(s, str) and bool(s))]
        libs[lname] = q[col].tolist()

    print("  computing representations ...", flush=True)
    R = {"train": represent(train_smiles), "held": represent(held_smiles),
         "NP-library": represent(libs["NP-library"]),
         "SNP-Lib": represent(libs["SNP-Lib"])}
    print(f"  done [{time.time()-t0:.0f}s]", flush=True)

    # standardise descriptors on TRAIN only
    Td = R["train"]["desc"]
    mu, sd = np.nanmean(Td, axis=0), np.nanstd(Td, axis=0)
    sd[sd == 0] = 1.0
    Z = {k: (R[k]["desc"] - mu) / sd for k in R}
    Zt = Z["train"]
    cov = np.cov(Zt, rowvar=False)
    cov_inv = np.linalg.pinv(cov)
    XtX_inv = np.linalg.pinv(Zt.T @ Zt)
    p_dim = Zt.shape[1]
    h_star = 3.0 * p_dim / len(Zt)

    def leverage(Zq):
        out = np.full(len(Zq), np.nan)
        for i in range(len(Zq)):
            if np.all(np.isfinite(Zq[i])):
                out[i] = float(Zq[i] @ XtX_inv @ Zq[i])
        return out

    def centroid_dist(Zq):
        c = Zt.mean(axis=0)
        return np.sqrt(((Zq - c) ** 2).sum(axis=1))

    def mahal(Zq):
        c = Zt.mean(axis=0)
        d = Zq - c
        return np.sqrt(np.einsum("ij,jk,ik->i", d, cov_inv, d))

    # ------------------------------------------------- the eight definitions
    defs = {}
    for name, rep, metric in (("A1_ecfp4_tanimoto_knn", "ecfp4", "tanimoto"),
                              ("A2_maccs_tanimoto_knn", "maccs", "tanimoto"),
                              ("A3_atompair_tanimoto_knn", "atompair", "tanimoto"),
                              ("A8_ecfp4_dice_knn", "ecfp4", "dice")):
        tfps = [f for f in R["train"][rep] if f is not None]
        self_sc = knn_fp_selfscores(tfps, metric)
        defs[name] = {
            "threshold": float(np.percentile(self_sc, AD_PCT)),
            "scores": {k: knn_fp_scores(R[k][rep], tfps, metric) for k in
                       ("held", "NP-library", "SNP-Lib")},
            "train_self": self_sc, "kind": "similarity"}
        print(f"  {name}: thr {defs[name]['threshold']:.4f} "
              f"[{time.time()-t0:.0f}s]", flush=True)

    self_e = knn_euclid_scores(Zt[random.Random(RS).sample(
        range(len(Zt)), min(SELFSIM_SAMPLE, len(Zt)))], Zt, k=AD_K + 1)
    defs["A4_desc_euclid_knn"] = {
        "threshold": float(np.percentile(self_e, AD_PCT)),
        "scores": {k: knn_euclid_scores(Z[k], Zt) for k in
                   ("held", "NP-library", "SNP-Lib")},
        "train_self": self_e, "kind": "neg_distance"}

    lev_t = leverage(Zt)
    defs["A5_desc_leverage_williams"] = {
        "threshold": float(-h_star),
        "scores": {k: -leverage(Z[k]) for k in ("held", "NP-library", "SNP-Lib")},
        "train_self": -lev_t, "kind": "neg_leverage",
        "h_star": float(h_star),
        "train_admitted_frac": float((lev_t <= h_star).mean())}

    cd_t = centroid_dist(Zt)
    defs["A6_desc_centroid_euclid"] = {
        "threshold": float(-np.percentile(cd_t, 100 - AD_PCT)),
        "scores": {k: -centroid_dist(Z[k]) for k in
                   ("held", "NP-library", "SNP-Lib")},
        "train_self": -cd_t, "kind": "neg_distance"}

    mh_t = mahal(Zt)
    defs["A7_desc_centroid_mahalanobis"] = {
        "threshold": float(-np.percentile(mh_t, 100 - AD_PCT)),
        "scores": {k: -mahal(Z[k]) for k in ("held", "NP-library", "SNP-Lib")},
        "train_self": -mh_t, "kind": "neg_distance"}
    print(f"  descriptor definitions done [{time.time()-t0:.0f}s]", flush=True)

    # ------------------------------------------------------------- measurement
    results = {}
    for name, D in defs.items():
        s_he = D["scores"]["held"]
        s_np = D["scores"]["NP-library"]
        s_snp = D["scores"]["SNP-Lib"]
        thr = D["threshold"]
        lo_np, hi_np = (float(np.nanpercentile(s_np, 1)),
                        float(np.nanpercentile(s_np, 99)))
        band = np.isfinite(s_he) & (s_he >= lo_np) & (s_he <= hi_np)
        n_b = int(band.sum())
        yb, pb = y_he[band], p_he[band]
        pt, lo, hi = bootstrap_auc(yb, pb) if n_b else (None, None, None)
        mcc = (float(matthews_corrcoef(yb, yhat_he[band]))
               if n_b and len(np.unique(yb)) > 1 else None)
        results[name] = {
            "threshold": thr, "kind": D["kind"],
            "np_occupied_range_p1_p99": [lo_np, hi_np],
            "np_full_range": [float(np.nanmin(s_np)), float(np.nanmax(s_np))],
            "heldout_range": [float(np.nanmin(s_he)), float(np.nanmax(s_he))],
            "in_domain": {
                "held": int(np.nansum(s_he >= thr)), "held_n": int(len(s_he)),
                "NP-library": int(np.nansum(s_np >= thr)), "NP_n": int(len(s_np)),
                "SNP-Lib": int(np.nansum(s_snp >= thr)), "SNP_n": int(len(s_snp))},
            "np_overlap_band": {
                "n": n_b, "active": int(yb.sum()) if n_b else 0,
                "inactive": int((yb == 0).sum()) if n_b else 0,
                "roc_auc": None if pt is None else round(pt, 4),
                "ci95": None if lo is None else [round(lo, 4), round(hi, 4)],
                "mcc": None if mcc is None else round(mcc, 4),
                "tier": tier(pt, lo),
                "delta_vs_full": None if pt is None else round(pt - full_auc, 4)},
        }
        if "h_star" in D:
            results[name]["h_star"] = D["h_star"]
            results[name]["train_admitted_frac"] = D["train_admitted_frac"]
        r = results[name]["np_overlap_band"]
        print(f"  {name:32s} band n={r['n']:4d} ({r['active']}/{r['inactive']})  "
              f"AUC {r['roc_auc']}  CI {r['ci95']}  {r['tier']}", flush=True)

    # ---------------------------------------------------------------- grading
    A1 = results["A1_ecfp4_tanimoto_knn"]["np_overlap_band"]
    tiers = {n: results[n]["np_overlap_band"]["tier"] for n in results}
    fp_defs = ["A1_ecfp4_tanimoto_knn", "A2_maccs_tanimoto_knn",
               "A3_atompair_tanimoto_knn", "A8_ecfp4_dice_knn"]
    desc_defs = ["A4_desc_euclid_knn", "A5_desc_leverage_williams",
                 "A6_desc_centroid_euclid", "A7_desc_centroid_mahalanobis"]
    np_in_fp = max(results[n]["in_domain"]["NP-library"] for n in fp_defs)
    np_in_desc = max(results[n]["in_domain"]["NP-library"] for n in desc_defs)

    grading = {
        "X1_degrades_by_gt_0.10": {
            "band_auc": A1["roc_auc"], "full_auc": round(full_auc, 4),
            "delta": A1["delta_vs_full"],
            "verdict": "MET" if (A1["roc_auc"] is not None
                                 and full_auc - A1["roc_auc"] > 0.10) else "MISSED"},
        "X2_at_least_tier1_under_A1": {
            "tier": A1["tier"], "ci95": A1["ci95"],
            "verdict": "MET" if A1["tier"] in ("Tier 1", "Tier 2") else "MISSED"},
        "X3_not_tier2_under_A1": {
            "tier": A1["tier"], "roc_auc": A1["roc_auc"],
            "verdict": "MET" if A1["tier"] != "Tier 2" else "MISSED"},
        "X4_definitions_disagree_on_tier": {
            "tiers": tiers, "n_distinct": len(set(tiers.values())),
            "verdict": "MET" if len(set(tiers.values())) >= 2 else "MISSED"},
        "X5_descriptor_admits_more_than_fingerprint": {
            "max_NP_in_domain_fingerprint": np_in_fp,
            "max_NP_in_domain_descriptor": np_in_desc,
            "per_definition": {n: results[n]["in_domain"]["NP-library"]
                               for n in results},
            "verdict": "MET" if np_in_desc > np_in_fp else "MISSED"},
        "X6_not_all_zero": {
            "any_nonzero": bool(max(results[n]["in_domain"]["NP-library"]
                                    for n in results) > 0),
            "verdict": "MET" if max(results[n]["in_domain"]["NP-library"]
                                    for n in results) > 0 else "MISSED"},
    }

    out = {"adr": "ADR-AIIM-292", "prereg_sha256": PREREG,
           "log_sha256_at_freeze": LOG_AT_FREEZE,
           "rdkit_version": rdkit.__version__,
           "curated_set": {"n": len(smiles), "train": len(tri),
                           "held_out": len(tei),
                           "held_out_active": int(y_he.sum()),
                           "held_out_inactive": int((y_he == 0).sum())},
           "full_heldout_auc": round(full_auc, 4),
           "tier_criterion": {"tier0": "bootstrap 95% CI lower <= 0.50",
                              "tier1": "CI lower > 0.50 and point <= 0.70",
                              "tier2": "CI lower > 0.50 and point > 0.70",
                              "n_bootstrap": N_BOOT},
           "threshold_convention": "each definition's threshold admits 90% of "
                                   "train; leverage uses Williams h*=3p/n and "
                                   "reports the fraction it admits",
           "definitions": results, "grading": grading,
           "seconds": round(time.time() - t0, 1)}
    (OUTDIR / "ad_methods_distance.json").write_text(json.dumps(out, indent=1),
                                                     encoding="utf-8")

    print("\n" + "=" * 76)
    for t, g in grading.items():
        print(f"{t:44s} {g['verdict']}")
    print("=" * 76)
    print(f"{'definition':32s} {'NP in-dom':>10s} {'held in-dom':>12s} "
          f"{'band n':>7s} {'AUC':>7s} {'tier':>8s}")
    for n, r in results.items():
        bb = r["np_overlap_band"]
        print(f"{n:32s} {r['in_domain']['NP-library']:10d} "
              f"{r['in_domain']['held']:12d} {bb['n']:7d} "
              f"{str(bb['roc_auc']):>7s} {bb['tier']:>8s}")
    print(f"\nwritten -> {OUTDIR / 'ad_methods_distance.json'} [{out['seconds']}s]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

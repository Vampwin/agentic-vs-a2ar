"""ADR-AIIM-310, executed. Pools, gates, windows and predictions frozen before
this ran.

Pre-registration SHA-256
  a0b8a137badd0cd37278b5751ede4156460f8556e7b5f969836551bd91c79250
Decision log at freeze (text/LF, per the 2026-10-05 addendum)
  d085b4a36203cdb915fde5a7b1b24bd3a7f77c059fc84191aa34f16144774e57

Three pools, two window conventions, eight split conventions.

THE MODEL IS FITTED ONCE PER SPLIT CONVENTION AND SCORES ALL THREE POOLS, so no
difference between pools can come from refitting.

GATES ARE IDENTICAL FOR BOTH LIBRARIES, which they never were before:
ADR-AIIM-300 applied only the AD gate to the NP library while ADR-AIIM-308
applied all four to SNP-Lib, so ADR-AIIM-302's six and ADR-AIIM-309's eight
have never been comparable. Ro5 is computed with RDKit for both rather than
read from a column that exists only for SNP-Lib -- anchored at 236/236 against
G1's own flag.

THE WINDOW IS A NINTH CONVENTION, not a nuisance parameter. ADR-AIIM-309
recorded that leaving it absolute confounds every cross-pool comparison; X3
measures that instead of recording it.

Pooling cannot reorder compounds within a library -- the model scores each
compound independently -- so the SNP-Lib subsequence of the pooled ranking is
necessarily the SNP-Lib-alone ranking. What pooling changes is how many slots
SNP-Lib wins against NP competition.

`10_new_system` is not modified.
"""
from __future__ import annotations

import importlib.util
import itertools
import json
import math
import os
import sys
import time
from collections import Counter
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "3")

import numpy as np
import pandas as pd
import rdkit
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Descriptors, rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

RDLogger.DisableLog("rdApp.*")

ROOT = Path(r"D:\Agentic VS Workflow")
SCRIPTS = ROOT / "10_new_system" / "scripts"
CUR_RUN = ROOT / "20_new_runs" / "chembl251_full_curation_2026-10-05"
REG_RUN = ROOT / "20_new_runs" / "pki_regression_2026-10-05"
EM_RUN = ROOT / "20_new_runs" / "pki_censored_em_2026-10-05"
SNPLIB = ROOT / "20_new_runs" / "snp_screen_2026-08-04" / "G1_standardization" / "snp_source_ledger.csv"
NPLIB = ROOT / "20_new_runs" / "np_library_funnel_2026-08-15" / "np_library_dock_ledger.csv"
OUTDIR = ROOT / "20_new_runs" / "pool_comparison_2026-10-06"
sys.path.insert(0, str(SCRIPTS))

PREREG = "a0b8a137badd0cd37278b5751ede4156460f8556e7b5f969836551bd91c79250"
LOG_AT_FREEZE = "d085b4a36203cdb915fde5a7b1b24bd3a7f77c059fc84191aa34f16144774e57"
RS, AD_FLOOR, MW_FLOOR = 42, 0.25, 200.0
TOP_ABS, WINDOW_FRAC = 20, 0.10
MAJORITY, SHORTLIST = 5, 10
CONFORMAL_INTERVAL = 1.5       # ADR-AIIM-299, the X4 reference
ADR309_SHORTLIST = ["SNP00163", "SNP00157", "SNP00010", "SNP00154",
                    "SNP00052", "SNP00144", "SNP00142", "SNP00006"]

_morgan = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
CATALOGS = {"PAINS": FilterCatalogParams.FilterCatalogs.PAINS,
            "BRENK": FilterCatalogParams.FilterCatalogs.BRENK,
            "NIH": FilterCatalogParams.FilterCatalogs.NIH,
            "CHEMBL": FilterCatalogParams.FilterCatalogs.CHEMBL}


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


CELLS = [{"d1_shuffle": a, "d2_train_first": b, "d3_divert": c}
         for a, b, c in itertools.product((True, False), (False, True), (False, True))]


def ro5(m):
    return (Descriptors.MolWt(m) <= 500 and Descriptors.MolLogP(m) <= 5
            and rdMolDescriptors.CalcNumHBD(m) <= 5
            and rdMolDescriptors.CalcNumHBA(m) <= 10)


def gate(df, id_col, smi_col, ref_fps, qsar, label):
    """dedupe -> Ro5 -> MW -> AD. Identical for both libraries."""
    wf = [{"pool": label, "step": 0, "name": "rows", "n": int(len(df))}]
    can, keep = [], []
    for s in df[smi_col]:
        m = Chem.MolFromSmiles(s) if isinstance(s, str) else None
        can.append(Chem.MolToSmiles(m) if m else None)
        keep.append(m is not None)
    df = df.assign(_can=can)[keep].reset_index(drop=True)
    n0 = len(df)
    df = df.sort_values(id_col).drop_duplicates("_can", keep="first").reset_index(drop=True)
    wf.append({"pool": label, "step": "0b", "name": "dedupe canonical SMILES",
               "n": int(len(df)), "rows_removed": int(n0 - len(df))})
    mols = [Chem.MolFromSmiles(s) for s in df["_can"]]
    m5 = np.array([ro5(m) for m in mols])
    df, mols = df[m5].reset_index(drop=True), [m for m, k in zip(mols, m5) if k]
    wf.append({"pool": label, "step": 1, "name": "gate Ro5 (RDKit)",
               "n": int(len(df))})
    mw = np.array([Descriptors.MolWt(m) >= MW_FLOOR for m in mols])
    df, mols = df[mw].reset_index(drop=True), [m for m, k in zip(mols, mw) if k]
    wf.append({"pool": label, "step": 2, "name": f"gate MW >= {MW_FLOOR:.0f}",
               "n": int(len(df))})
    ad = np.array([float(np.mean(sorted(
        DataStructs.BulkTanimotoSimilarity(_morgan.GetFingerprint(m), ref_fps),
        reverse=True)[:5])) for m in mols])
    ok = ad >= AD_FLOOR
    df = df[ok].reset_index(drop=True)
    df["ad_score"] = ad[ok]
    df["murcko"] = [qsar._murcko(s) for s in df["_can"]]
    df["source"] = label
    df["uid"] = df[id_col].astype(str)
    wf.append({"pool": label, "step": 3, "name": f"gate AD >= {AD_FLOOR}",
               "n": int(len(df)), "frameworks": int(df["murcko"].nunique())})
    for w in wf:
        print(f"  [{label}:{w['step']}] {w['name']:30s} n={w['n']:5d}"
              + (f" fw={w['frameworks']}" if "frameworks" in w else ""), flush=True)
    return df, wf


def flags(smiles, cats):
    out = []
    for s in smiles:
        m = Chem.MolFromSmiles(s)
        out.append({k: sorted({x.GetDescription() for x in fc.GetMatches(m)})
                    if m else [] for k, fc in cats.items()})
    return out


def select(pool_df, preds, window_k):
    """top-K per convention -> 5-of-8 majority -> scaffold uniqueness -> 10."""
    cnt = Counter()
    for nm, p in preds.items():
        top = np.argsort(-p)[:window_k]
        cnt.update(pool_df.loc[top, "uid"].tolist())
    surv = {k: v for k, v in cnt.items() if v >= MAJORITY}
    mean = dict(zip(pool_df["uid"], np.vstack(list(preds.values())).mean(axis=0)))
    ranked = sorted(surv, key=lambda c: (-surv[c], -mean[c]))
    info = pool_df.set_index("uid")
    sl, seen = [], set()
    for c in ranked:
        f = info.loc[c, "murcko"]
        if f in seen:
            continue
        seen.add(f)
        sl.append(c)
        if len(sl) == SHORTLIST:
            break
    return {"n_survivors": len(surv), "survivors": surv,
            "n_unanimous": sum(1 for v in cnt.values() if v == len(preds)),
            "shortlist": sl, "mean_pred": mean, "window_k": window_k}


def main() -> int:
    (OUTDIR / "analysis").mkdir(parents=True, exist_ok=True)
    qsar = load("_qsar", SCRIPTS / "run_qsar_lbdd_2026-08-07.py")
    cur = load("_cur", CUR_RUN / "analysis" / "chembl251_curate.py")
    rg = load("_reg", REG_RUN / "analysis" / "pki_regression.py")
    em = load("_em", EM_RUN / "analysis" / "pki_censored_em.py")
    t0 = time.time()

    tbl, _ = rg.rebuild_table(cur)
    b = tbl[tbl["endpoint_class"] == "binding"]
    b = b[b["basis"].isin(["point_median", "point_median_gap",
                           "censored_upper", "censored_lower"])]
    b = b[np.isfinite(b["px"])].reset_index(drop=True)
    kind = np.where(b["basis"] == "censored_upper", 1,
                    np.where(b["basis"] == "censored_lower", -1, 0)).astype(int)
    msm = b["std_smiles"].tolist()
    y0 = b["px"].to_numpy(float)
    iks = [qsar._smiles_to_inchikey(s) or f"__N{i}__" for i, s in enumerate(msm)]
    o = sorted(range(len(msm)), key=lambda j: iks[j])
    msm = [msm[i] for i in o]
    y0, kind = y0[o], kind[o]
    scaf = [qsar._murcko(s) for s in msm]
    X = qsar.compute_raw_features(msm)
    ref = [_morgan.GetFingerprint(Chem.MolFromSmiles(s)) for s in msm
           if Chem.MolFromSmiles(s)]
    print(f"modelling set {len(y0)}  AD reference {len(ref)}", flush=True)

    sl_raw = pd.read_csv(SNPLIB)
    sl_raw = sl_raw[sl_raw["smiles_canonical"].notna()].reset_index(drop=True)
    A, wfA = gate(sl_raw, "SNP_ID", "smiles_canonical", ref, qsar, "SNP-Lib")
    np_raw = pd.read_csv(NPLIB)
    np_raw = np_raw[np_raw["SMILES"].notna()].reset_index(drop=True)
    B, wfB = gate(np_raw, "compound_id", "SMILES", ref, qsar, "NP-library")

    # Ro5 anchor: RDKit against G1's own flag
    anchor = sum(1 for s, f in zip(sl_raw["smiles_canonical"],
                                   sl_raw["lipinski_ro5_pass"])
                 if (m := Chem.MolFromSmiles(s)) is not None and ro5(m) == bool(f))
    print(f"  Ro5 anchor vs G1 flag: {anchor}/{len(sl_raw)}", flush=True)

    keepc = ["uid", "_can", "ad_score", "murcko", "source"]
    C = pd.concat([A[keepc], B[keepc]], ignore_index=True)
    pools = {"A_SNPLib": A[keepc], "B_NPLib": B[keepc], "C_pooled": C}
    feats = {k: qsar.compute_raw_features(v["_can"].tolist()) for k, v in pools.items()}
    print(f"  pools: " + "  ".join(f"{k}={len(v)}" for k, v in pools.items())
          + f"  [{time.time()-t0:.0f}s]", flush=True)

    # one fit per convention, scoring all three pools
    preds = {k: {} for k in pools}
    for cfg in CELLS:
        nm = em.cell_name(cfg)
        tr, _te = em.split(scaf, len(y0), **cfg)
        mdl, n_it, chg, _ = em.censored_em(X, y0, kind, tr)
        for k in pools:
            preds[k][nm] = mdl.predict(feats[k])
        print(f"  {nm:40s} EM {n_it}it  [{time.time()-t0:.0f}s]", flush=True)

    cats = {k: FilterCatalog((lambda c: (lambda p: (p.AddCatalog(c), p)[1])(
        FilterCatalogParams()))(v)) for k, v in CATALOGS.items()}

    res, rows = {}, []
    for pk, pdf in pools.items():
        for wk, k in (("abs", TOP_ABS),
                      ("frac", max(int(math.ceil(WINDOW_FRAC * len(pdf))), 1))):
            s = select(pdf.reset_index(drop=True), preds[pk], k)
            info = pdf.set_index("uid")
            n_snp = sum(1 for c in s["shortlist"] if info.loc[c, "source"] == "SNP-Lib")
            surv_snp = sum(1 for c in s["survivors"]
                           if info.loc[c, "source"] == "SNP-Lib")
            fl = flags([info.loc[c, "_can"] for c in s["shortlist"]], cats)
            res[f"{pk}|{wk}"] = {
                "pool": pk, "window": wk, "window_k": k,
                "pool_n": int(len(pdf)),
                "window_frac": round(k / len(pdf), 4),
                "n_survivors": s["n_survivors"],
                "n_unanimous": s["n_unanimous"],
                "survivors_snplib": surv_snp,
                "shortlist": s["shortlist"],
                "shortlist_snplib": n_snp,
                "shortlist_snplib_frac": round(n_snp / max(len(s["shortlist"]), 1), 4),
                "top_pred_pki": round(max(s["mean_pred"][c]
                                          for c in s["shortlist"]), 3)
                if s["shortlist"] else None,
                "survivor_set": sorted(s["survivors"])}
            for i, c in enumerate(s["shortlist"], 1):
                rows.append({"pool": pk, "window": wk, "rank": i, "uid": c,
                             "source": info.loc[c, "source"],
                             "mean_pred_pki": round(s["mean_pred"][c], 3),
                             "ad_score": round(float(info.loc[c, "ad_score"]), 4),
                             "murcko": info.loc[c, "murcko"],
                             "PAINS": "; ".join(fl[i-1]["PAINS"]) or "-",
                             "BRENK": "; ".join(fl[i-1]["BRENK"]) or "-",
                             "NIH": "; ".join(fl[i-1]["NIH"]) or "-",
                             "smiles": info.loc[c, "_can"]})
            r = res[f"{pk}|{wk}"]
            print(f"  {pk:10s} {wk:4s} k={k:4d} ({r['window_frac']:.1%})  "
                  f"surv={r['n_survivors']:4d} (SNP {surv_snp:3d})  "
                  f"unan={r['n_unanimous']:3d}  SL={len(s['shortlist'])} "
                  f"(SNP {n_snp})  top={r['top_pred_pki']}", flush=True)
    pd.DataFrame(rows).to_csv(OUTDIR / "shortlists.csv", index=False)

    # grading
    a_abs = res["A_SNPLib|abs"]
    x1 = a_abs["shortlist"] == ADR309_SHORTLIST
    c_abs, c_frac = res["C_pooled|abs"], res["C_pooled|frac"]
    s1, s2 = set(c_abs["survivor_set"]), set(c_frac["survivor_set"])
    ov = len(s1 & s2) / max(min(len(s1), len(s2)), 1)
    topA = max(res["A_SNPLib|abs"]["top_pred_pki"] or 0,
               res["A_SNPLib|frac"]["top_pred_pki"] or 0)
    topB = max(res["B_NPLib|abs"]["top_pred_pki"] or 0,
               res["B_NPLib|frac"]["top_pred_pki"] or 0)

    out = {"adr": "ADR-AIIM-310", "prereg_sha256": PREREG,
           "log_sha256_at_freeze": LOG_AT_FREEZE,
           "rdkit_version": rdkit.__version__,
           "ro5_anchor_vs_g1": f"{anchor}/{len(sl_raw)}",
           "gates": ["dedupe canonical SMILES", "Ro5 (RDKit, both libraries)",
                     f"MW >= {MW_FLOOR:.0f}", f"AD >= {AD_FLOOR}"],
           "waterfall": wfA + wfB,
           "pool_sizes": {k: int(len(v)) for k, v in pools.items()},
           "snplib_share_of_pooled": round(len(A) / len(C), 4),
           "cells": res,
           "grading": {
               "X1_anchor_reproduces_adr309": {
                   "observed": a_abs["shortlist"], "expected": ADR309_SHORTLIST,
                   "verdict": "MET" if x1 else "MISSED"},
               "X2_snplib_survives_pooled": {
                   "survivors_snplib": c_abs["survivors_snplib"],
                   "shortlist_snplib": c_abs["shortlist_snplib"],
                   "verdict": "MET" if c_abs["survivors_snplib"] >= 1 else "MISSED"},
               "X3_window_overlap_below_50pct": {
                   "overlap_frac": round(ov, 4),
                   "n_abs": len(s1), "n_frac": len(s2),
                   "verdict": "MET" if ov < 0.50 else "MISSED"},
               "X4_libraries_indistinguishable": {
                   "top_A": topA, "top_B": topB,
                   "gap": round(abs(topA - topB), 3),
                   "interval": CONFORMAL_INTERVAL,
                   "verdict": "MET" if abs(topA - topB) < CONFORMAL_INTERVAL
                   else "MISSED"},
               "X5_nplib_majority_of_pooled_shortlist": {
                   "snplib_frac": c_abs["shortlist_snplib_frac"],
                   "verdict": "MET" if c_abs["shortlist_snplib_frac"] < 0.5
                   else "MISSED"}},
           "caveats": [
               "Pooling is measured, not recommended: a pooled shortlist mixes "
               "assayable with non-assayable compounds (ADR-AIIM-310 c.3).",
               "Pooling cannot reorder within a library; only slot competition "
               "changes.",
               "No shortlisted compound is a predicted active. ADR-AIIM-299: "
               "conformal count 0 at pKi >= 6.0 under every framing.",
               "X3 measures the window's consequence, not its correctness.",
               "Alerts are flags; nothing was removed on their strength.",
               "Nothing in this project has been assayed."],
           "seconds": round(time.time() - t0, 1)}
    (OUTDIR / "pool_comparison.json").write_text(json.dumps(out, indent=1),
                                                 encoding="utf-8")
    (OUTDIR / "waterfall.json").write_text(json.dumps(
        {"waterfall": wfA + wfB, "pool_sizes": out["pool_sizes"],
         "ro5_anchor": out["ro5_anchor_vs_g1"]}, indent=1), encoding="utf-8")

    print("\n" + "=" * 80)
    for t, g in out["grading"].items():
        print(f"{t:42s} {g['verdict']}")
    print("=" * 80)
    print(f"pools: A={len(A)}  B={len(B)}  C={len(C)}   "
          f"SNP share of pooled {out['snplib_share_of_pooled']:.1%}")
    print(f"\n{'cell':16s} {'k':>5s} {'win%':>6s} {'surv':>5s} {'SNPsurv':>8s} "
          f"{'unan':>5s} {'SL':>3s} {'SNPinSL':>8s} {'top pKi':>8s}")
    for k, r in res.items():
        print(f"{k:16s} {r['window_k']:5d} {r['window_frac']:6.1%} "
              f"{r['n_survivors']:5d} {r['survivors_snplib']:8d} "
              f"{r['n_unanimous']:5d} {len(r['shortlist']):3d} "
              f"{r['shortlist_snplib']:8d} {str(r['top_pred_pki']):>8s}")
    print(f"\nwritten -> {OUTDIR} [{out['seconds']}s]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""ADR-AIIM-311, executed. Predictions frozen before this ran.

Pre-registration SHA-256
  bde7df91688656656406d6d7a64470a100bc043d5d6d80405b94bfd93757742d
Decision log at freeze (text/LF, per the 2026-10-05 addendum)
  2c78bd0b186f2230d39a26502b62413d2ac17b83eaaf6fc3361a2503ef3cfa99

The censored-EM results of ADR-AIIM-298 to 310 are not reproducible. The cause
is a random forest at n_jobs=3 inside a ten-iteration feedback loop: a single
fit carries float noise of order 1e-15 from thread accumulation order, forest
splits are threshold comparisons, and ten rounds of feeding predictions back
into targets amplify that to 0.127 -- enough to change a top-20.

Nothing here re-docks, re-pulls or refits the classifier work, and
`10_new_system` is not modified. The fix is n_jobs=1 and lives in this file.

Three measurements, all CPU:
  M1  determinism at n_jobs=1, three repeats on the anchor convention
  M2  replicate spread at n_jobs=3, the setting every prior entry used,
      five repeats per convention -- the error bar those entries should
      have carried
  M3  the corrected eight-convention run at n_jobs=1, with the SNP-Lib
      shortlist rebuilt from it through ADR-AIIM-310 gate and majority-rule
      code, imported rather than copied
"""
from __future__ import annotations

import importlib.util
import json
import os
import statistics
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "3")

import numpy as np
import pandas as pd
import rdkit
import sklearn
from rdkit import Chem, RDLogger
from sklearn.ensemble import RandomForestRegressor
from sklearn.feature_selection import VarianceThreshold
from sklearn.pipeline import Pipeline

RDLogger.DisableLog("rdApp.*")

ROOT = Path(r"D:\Agentic VS Workflow")
SCRIPTS = ROOT / "10_new_system" / "scripts"
CUR_RUN = ROOT / "20_new_runs" / "chembl251_full_curation_2026-10-05"
REG_RUN = ROOT / "20_new_runs" / "pki_regression_2026-10-05"
EM_RUN = ROOT / "20_new_runs" / "pki_censored_em_2026-10-05"
PC_RUN = ROOT / "20_new_runs" / "pool_comparison_2026-10-06"
NPLIB = ROOT / "20_new_runs" / "np_library_funnel_2026-08-15" / "np_library_dock_ledger.csv"
SNPLIB = ROOT / "20_new_runs" / "snp_screen_2026-08-04" / "G1_standardization" / "snp_source_ledger.csv"
OUTDIR = ROOT / "20_new_runs" / "em_reproducibility_2026-10-06"
sys.path.insert(0, str(SCRIPTS))

PREREG = "bde7df91688656656406d6d7a64470a100bc043d5d6d80405b94bfd93757742d"
LOG_AT_FREEZE = "2c78bd0b186f2230d39a26502b62413d2ac17b83eaaf6fc3361a2503ef3cfa99"

RS, SIG = 42, 0.05
THRESHOLDS = [5.0, 5.5, 6.0, 6.5, 7.0]
M1_REPEATS, M2_REPEATS = 3, 5

# the numbers these measurements are checked against, all produced at n_jobs=3
ADR298_SPREAD = 4.10
ADR298_NP6 = {"D1shuffle_D2testfirst_D3overflow": 29,
              "D1shuffle_D2testfirst_D3divert": 29,
              "D1shuffle_D2trainfirst_D3overflow": 71,
              "D1shuffle_D2trainfirst_D3divert": 52,
              "D1sorted_D2testfirst_D3overflow": 39,
              "D1sorted_D2testfirst_D3divert": 60,
              "D1sorted_D2trainfirst_D3overflow": 78,
              "D1sorted_D2trainfirst_D3divert": 119}
# the ADR-AIIM-309 eight-structure priority order, SNP-Lib alone, absolute top-20
ADR309_TOP_PKI = 5.88


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def make_pipe(n_jobs):
    def _p():
        return Pipeline([("var", VarianceThreshold(threshold=0.01)),
                         ("rf", RandomForestRegressor(random_state=RS,
                                                      n_jobs=n_jobs))])
    return _p


def adr309_shortlist():
    """Read the reported shortlist from the ADR-AIIM-309 artifact rather than
    hard-coding it, so the comparison cannot drift from what was published."""
    p = ROOT / "20_new_runs" / "snplib_dedup_2026-10-05" / "shortlist.csv"
    if not p.exists():
        return []
    df = pd.read_csv(p)
    for col in ("SNP_ID", "compound_id", "uid"):
        if col in df.columns:
            return df[col].astype(str).tolist()
    return df.iloc[:, 0].astype(str).tolist()


def main() -> int:
    (OUTDIR / "analysis").mkdir(parents=True, exist_ok=True)
    qsar = load("_qsar", SCRIPTS / "run_qsar_lbdd_2026-08-07.py")
    cur = load("_cur", CUR_RUN / "analysis" / "chembl251_curate.py")
    rg = load("_reg", REG_RUN / "analysis" / "pki_regression.py")
    em = load("_em", EM_RUN / "analysis" / "pki_censored_em.py")
    pc = load("_pc", PC_RUN / "analysis" / "pool_comparison.py")
    t0 = time.time()

    reported = adr309_shortlist()
    print(f"ADR-309 reported shortlist ({len(reported)}): {reported}", flush=True)

    # ---- data, exactly as ADR-AIIM-298 and 310 build it -------------------
    tbl, n_raw = rg.rebuild_table(cur)
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
    print(f"modelling set {len(y0)}  point {(kind==0).sum()}  "
          f"upper {(kind==1).sum()}  lower {(kind==-1).sum()}", flush=True)

    npl = pd.read_csv(NPLIB)
    npl = npl[npl["SMILES"].notna()].reset_index(drop=True)
    X_np = qsar.compute_raw_features(npl["SMILES"].tolist())

    anchor = em.cell_name({"d1_shuffle": False, "d2_train_first": True,
                           "d3_divert": True})
    cells = [em.cell_name(c) for c in em.CELLS]
    print(f"features ready, anchor {anchor}  [{time.time()-t0:.0f}s]", flush=True)

    out = {"adr": "ADR-AIIM-311", "prereg_sha256": PREREG,
           "log_sha256_at_freeze": LOG_AT_FREEZE,
           "rdkit_version": rdkit.__version__,
           "sklearn_version": sklearn.__version__,
           "defect": ("RandomForestRegressor at n_jobs>1 is non-deterministic "
                      "at 1e-15; the ten-iteration censored-EM feedback loop "
                      "amplifies that to 1.3e-1"),
           "modelling_set": int(len(y0)), "raw_rows": int(n_raw),
           "np_library_rows": int(len(npl)),
           "reference_n_jobs3": {"adr298_spread": ADR298_SPREAD,
                                 "adr298_np6": ADR298_NP6,
                                 "adr309_shortlist": reported,
                                 "adr309_top_pki": ADR309_TOP_PKI}}

    def run_cell(cell, n_jobs):
        em.pipe = make_pipe(n_jobs)
        cfg = next(c for c in em.CELLS if em.cell_name(c) == cell)
        tr, te = em.split(scaf, len(y0), **cfg)
        mdl, n_it, chg, _ = em.censored_em(X, y0, kind, tr)
        te_pt = te[kind[te] == 0]
        p_pt = mdl.predict(X[te_pt])
        resid = np.abs(y0[te_pt] - p_pt)
        k = int(np.ceil((1 - SIG) * (len(resid) + 1))) - 1
        qhat = float(np.sort(resid)[min(k, len(resid) - 1)])
        p_np = mdl.predict(X_np)
        return mdl, {"em_iterations": n_it,
                     "em_final_max_change": round(float(chg), 5),
                     "conformal_qhat": round(qhat, 4),
                     "np_point": {str(t): int((p_np >= t).sum())
                                  for t in THRESHOLDS},
                     "np_conformal_lower": {str(t): int((p_np - qhat >= t).sum())
                                            for t in THRESHOLDS},
                     "np_pred_min": round(float(p_np.min()), 4),
                     "np_pred_max": round(float(p_np.max()), 4)}, p_np

    # ---- M1: determinism at n_jobs=1 -------------------------------------
    print(f"\n== M1  n_jobs=1, {M1_REPEATS} repeats on {anchor}", flush=True)
    m1 = {"n_jobs": 1, "convention": anchor, "repeats": M1_REPEATS, "runs": []}
    ref_np = None
    for r in range(M1_REPEATS):
        _m, rec, p_np = run_cell(anchor, 1)
        if ref_np is None:
            ref_np, d = p_np, 0.0
        else:
            d = float(np.max(np.abs(p_np - ref_np)))
        rec["max_abs_diff_vs_repeat1"] = d
        m1["runs"].append(rec)
        print(f"  repeat {r+1}: EM {rec['em_iterations']}it  "
              f"change {rec['em_final_max_change']:.5f}  "
              f"NP@6 {rec['np_point']['6.0']}  "
              f"max|diff| {d:.3e}  [{time.time()-t0:.0f}s]", flush=True)
    m1["all_identical"] = all(r["max_abs_diff_vs_repeat1"] == 0.0
                              for r in m1["runs"])
    m1["max_abs_diff"] = max(r["max_abs_diff_vs_repeat1"] for r in m1["runs"])
    out["M1_determinism_n_jobs_1"] = m1

    # ---- M2: replicate spread at n_jobs=3 --------------------------------
    print(f"\n== M2  n_jobs=3, {M2_REPEATS} repeats x {len(cells)} conventions",
          flush=True)
    m2 = {"n_jobs": 3, "repeats": M2_REPEATS, "by_convention": {}}
    for cell in cells:
        np6, chgs, qh = [], [], []
        base, maxd = None, 0.0
        for r in range(M2_REPEATS):
            _m, rec, p_np = run_cell(cell, 3)
            np6.append(rec["np_point"]["6.0"])
            chgs.append(rec["em_final_max_change"])
            qh.append(rec["conformal_qhat"])
            if base is None:
                base = p_np
            else:
                maxd = max(maxd, float(np.max(np.abs(p_np - base))))
        m2["by_convention"][cell] = {
            "np6_per_repeat": np6, "np6_min": min(np6), "np6_max": max(np6),
            "np6_range": max(np6) - min(np6),
            "np6_mean": round(statistics.mean(np6), 2),
            "np6_sd": round(statistics.stdev(np6), 2) if len(np6) > 1 else 0.0,
            "em_change_per_repeat": chgs,
            "qhat_per_repeat": qh,
            "max_abs_pred_diff": maxd,
            "adr298_reported": ADR298_NP6[cell],
            "adr298_inside_replicate_range": bool(
                min(np6) <= ADR298_NP6[cell] <= max(np6))}
        c = m2["by_convention"][cell]
        print(f"  {cell:40s} NP@6 {np6}  range {c['np6_range']:3d}  "
              f"sd {c['np6_sd']:5.2f}  max|diff| {maxd:.3e}  "
              f"ADR298={c['adr298_reported']} "
              f"{'IN' if c['adr298_inside_replicate_range'] else 'OUT'}  "
              f"[{time.time()-t0:.0f}s]", flush=True)
    within = [v["np6_range"] for v in m2["by_convention"].values()]
    between = max(ADR298_NP6.values()) - min(ADR298_NP6.values())
    m2["median_within_convention_range"] = statistics.median(within)
    m2["between_convention_range_adr298"] = between
    m2["ratio_within_over_between"] = round(
        statistics.median(within) / between, 4) if between else None
    m2["n_conventions_whose_replicate_range_covers_adr298"] = sum(
        v["adr298_inside_replicate_range"] for v in m2["by_convention"].values())
    out["M2_replicate_spread_n_jobs_3"] = m2

    # ---- M3: corrected run at n_jobs=1 -----------------------------------
    print(f"\n== M3  n_jobs=1, {len(cells)} conventions (corrected)", flush=True)
    m3 = {"n_jobs": 1, "by_convention": {}}
    models = {}
    for cell in cells:
        mdl, rec, _p = run_cell(cell, 1)
        models[cell] = mdl
        rec["adr298_reported_np6"] = ADR298_NP6[cell]
        rec["delta_vs_adr298"] = rec["np_point"]["6.0"] - ADR298_NP6[cell]
        m3["by_convention"][cell] = rec
        print(f"  {cell:40s} EM {rec['em_iterations']}it  "
              f"NP@6 {rec['np_point']['6.0']:4d} "
              f"(ADR298 {ADR298_NP6[cell]:4d}, d {rec['delta_vs_adr298']:+4d})  "
              f"conf@6 {rec['np_conformal_lower']['6.0']}  "
              f"qhat {rec['conformal_qhat']:.3f}  [{time.time()-t0:.0f}s]",
              flush=True)
    np6c = [m3["by_convention"][c]["np_point"]["6.0"] for c in cells]
    m3["np6_min"], m3["np6_max"] = min(np6c), max(np6c)
    m3["spread_fold"] = round(max(np6c) / max(min(np6c), 1), 2)
    m3["conformal_at_6_all_zero"] = all(
        m3["by_convention"][c]["np_conformal_lower"]["6.0"] == 0 for c in cells)
    out["M3_corrected_n_jobs_1"] = m3

    # ---- M3b: the SNP-Lib shortlist, rebuilt through ADR-310 own code ----
    print("\n== M3b  SNP-Lib shortlist at n_jobs=1, "
          "gates and majority rule imported from ADR-AIIM-310", flush=True)
    sl_raw = pd.read_csv(SNPLIB)
    sl_raw = sl_raw[sl_raw["smiles_canonical"].notna()].reset_index(drop=True)
    ref_fps = [pc._morgan.GetFingerprint(Chem.MolFromSmiles(s)) for s in msm
               if Chem.MolFromSmiles(s)]
    A, wfA = pc.gate(sl_raw, "SNP_ID", "smiles_canonical", ref_fps, qsar, "SNP-Lib")
    keepc = ["uid", "_can", "ad_score", "murcko", "source"]
    A = A[keepc].reset_index(drop=True)
    X_A = qsar.compute_raw_features(A["_can"].tolist())
    predsA = {c: models[c].predict(X_A) for c in cells}
    sel = pc.select(A, predsA, pc.TOP_ABS)
    sl = sel["shortlist"]
    mean = sel["mean_pred"]
    m3b = {"pool_size": int(len(A)), "window_k": pc.TOP_ABS,
           "majority": pc.MAJORITY,
           "waterfall": wfA,
           "n_survivors": sel["n_survivors"],
           "n_unanimous": sel["n_unanimous"],
           "shortlist": sl,
           "shortlist_mean_pki": {c: round(float(mean[c]), 4) for c in sl},
           "shortlist_conventions": {c: int(sel["survivors"][c]) for c in sl},
           "top_mean_pki": round(float(max(mean[c] for c in sl)), 4) if sl else None,
           "any_crosses_6.0": bool(any(mean[c] >= 6.0 for c in sl)),
           "adr309_shortlist": reported,
           "adr309_top_pki": ADR309_TOP_PKI}
    m3b["overlap_with_adr309"] = sorted(set(sl) & set(reported))
    m3b["n_overlap"] = len(m3b["overlap_with_adr309"])
    m3b["identical_to_adr309"] = (sl == reported)
    m3b["jaccard_vs_adr309"] = round(
        len(set(sl) & set(reported)) / max(len(set(sl) | set(reported)), 1), 4)
    out["M3b_snplib_shortlist_corrected"] = m3b
    print(f"  pool {len(A)}  survivors {sel['n_survivors']}  "
          f"unanimous {sel['n_unanimous']}  shortlist {len(sl)}", flush=True)
    for c in sl:
        print(f"    {c:12s} {sel['survivors'][c]}/8 conv  "
              f"pKi {mean[c]:.4f}", flush=True)
    print(f"  overlap with ADR-309: {m3b['n_overlap']}  "
          f"Jaccard {m3b['jaccard_vs_adr309']}", flush=True)

    # ---- grading ---------------------------------------------------------
    y2_delta = abs(m3["spread_fold"] - ADR298_SPREAD) / ADR298_SPREAD
    ratio = m2["ratio_within_over_between"]
    out["grading"] = {
        "Y1_n_jobs_1_bit_identical": {
            "repeats": M1_REPEATS, "max_abs_diff": m1["max_abs_diff"],
            "verdict": "MET" if m1["all_identical"] else "MISSED"},
        "Y2_corrected_spread_differs_by_over_20pct": {
            "corrected": m3["spread_fold"], "reported": ADR298_SPREAD,
            "relative_change": round(y2_delta, 4),
            "verdict": "MET" if y2_delta > 0.20 else "MISSED"},
        "Y3_run_noise_at_least_half_the_convention_effect": {
            "median_within_convention_range": m2["median_within_convention_range"],
            "between_convention_range": between,
            "ratio": ratio,
            "n_conventions_whose_replicate_range_covers_the_reported_value":
                m2["n_conventions_whose_replicate_range_covers_adr298"],
            "verdict": "MET" if (ratio is not None and ratio >= 0.5) else "MISSED"},
        "Y4_qualitative_conclusions_survive": {
            "conformal_count_at_6.0_still_zero": m3["conformal_at_6_all_zero"],
            "no_shortlist_compound_crosses_6.0": not m3b["any_crosses_6.0"],
            "top_shortlist_pki": m3b["top_mean_pki"],
            "verdict": "MET" if (m3["conformal_at_6_all_zero"]
                                 and not m3b["any_crosses_6.0"]) else "MISSED"},
        "Y5_shortlist_changes": {
            "corrected": sl, "reported": reported,
            "n_overlap": m3b["n_overlap"],
            "jaccard": m3b["jaccard_vs_adr309"],
            "verdict": "MET" if not m3b["identical_to_adr309"] else "MISSED"},
    }
    out["seconds"] = round(time.time() - t0, 1)
    (OUTDIR / "em_reproducibility.json").write_text(
        json.dumps(out, indent=1), encoding="utf-8")

    print("\n" + "=" * 78)
    for k, g in out["grading"].items():
        print(f"{k:56s} {g['verdict']}")
    print("=" * 78)
    print(f"M1 n_jobs=1 max|diff| over {M1_REPEATS} repeats: "
          f"{m1['max_abs_diff']:.3e}")
    print(f"M2 median within-convention NP@6 range "
          f"{m2['median_within_convention_range']}"
          f"  vs between-convention {between}  ratio {ratio}")
    print(f"M3 corrected spread {m3['spread_fold']}x  "
          f"(reported {ADR298_SPREAD}x, {y2_delta*100:.1f}% change)")
    print(f"M3 NP@6 range {m3['np6_min']}-{m3['np6_max']}  "
          f"conformal@6 all zero: {m3['conformal_at_6_all_zero']}")
    print(f"M3b shortlist {len(sl)} cmpds, top pKi {m3b['top_mean_pki']}, "
          f"overlap with ADR-309 {m3b['n_overlap']}")
    print(f"\nwritten -> {OUTDIR / 'em_reproducibility.json'} [{out['seconds']}s]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

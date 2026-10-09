"""ADR-AIIM-312 addendum: recover the full survivor set the ADR-AIIM-311
harness computed but did not save.

That harness stored only the eight shortlisted compounds, so the identity of
the fifteen survivors -- and of the unanimous compound the scaffold-uniqueness
rule dropped -- could not be read back from the artifact. This re-derives them.

Decision log at launch (text/LF)
  35f117716b5dd093a4bf4787c323c43a3523f3edd4785050af52fad524b27472

No new method and no new prediction, with one anchor stated before the run:
n_jobs=1 was shown bit-identical in ADR-AIIM-312, so this must reproduce the
eight-compound shortlist in the same order, 15 survivors and 5 unanimous. If it
does not, an unidentified source of non-determinism remains, and that matters
more than the compound names this run was launched to recover.

`10_new_system` is not modified. Nothing is re-docked or re-pulled.
"""
from __future__ import annotations

import importlib.util
import json
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
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, QED, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
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
REPRO = ROOT / "20_new_runs" / "em_reproducibility_2026-10-06"
SNPLIB = (ROOT / "20_new_runs" / "snp_screen_2026-08-04" / "G1_standardization"
          / "snp_source_ledger.csv")
OUTDIR = ROOT / "20_new_runs" / "em_survivors_2026-10-06"
sys.path.insert(0, str(SCRIPTS))

LOG_AT_LAUNCH = "35f117716b5dd093a4bf4787c323c43a3523f3edd4785050af52fad524b27472"
RS = 42

CATS = {"PAINS": FilterCatalogParams.FilterCatalogs.PAINS,
        "BRENK": FilterCatalogParams.FilterCatalogs.BRENK,
        "NIH": FilterCatalogParams.FilterCatalogs.NIH,
        "CHEMBL": FilterCatalogParams.FilterCatalogs.CHEMBL}


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def pipe1():
    return Pipeline([("var", VarianceThreshold(threshold=0.01)),
                     ("rf", RandomForestRegressor(random_state=RS, n_jobs=1))])


def catalog(c):
    p = FilterCatalogParams()
    p.AddCatalog(c)
    return FilterCatalog(p)


def main() -> int:
    (OUTDIR / "analysis").mkdir(parents=True, exist_ok=True)
    qsar = load("_qsar", SCRIPTS / "run_qsar_lbdd_2026-08-07.py")
    cur = load("_cur", CUR_RUN / "analysis" / "chembl251_curate.py")
    rg = load("_reg", REG_RUN / "analysis" / "pki_regression.py")
    em = load("_em", EM_RUN / "analysis" / "pki_censored_em.py")
    pc = load("_pc", PC_RUN / "analysis" / "pool_comparison.py")
    em.pipe = pipe1
    t0 = time.time()

    prev = json.loads((REPRO / "em_reproducibility.json").read_text("utf-8"))
    prev_sl = prev["M3b_snplib_shortlist_corrected"]["shortlist"]
    prev_nsurv = prev["M3b_snplib_shortlist_corrected"]["n_survivors"]
    prev_nunan = prev["M3b_snplib_shortlist_corrected"]["n_unanimous"]

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
    print(f"modelling set {len(y0)}  [{time.time()-t0:.0f}s]", flush=True)

    sl_raw = pd.read_csv(SNPLIB)
    sl_raw = sl_raw[sl_raw["smiles_canonical"].notna()].reset_index(drop=True)
    ref = [pc._morgan.GetFingerprint(Chem.MolFromSmiles(s)) for s in msm
           if Chem.MolFromSmiles(s)]
    A, wf = pc.gate(sl_raw, "SNP_ID", "smiles_canonical", ref, qsar, "SNP-Lib")
    keepc = ["uid", "_can", "ad_score", "murcko", "source"]
    A = A[keepc].reset_index(drop=True)
    X_A = qsar.compute_raw_features(A["_can"].tolist())

    cells = [em.cell_name(c) for c in em.CELLS]
    preds = {}
    for cfg in em.CELLS:
        nm = em.cell_name(cfg)
        tr, _te = em.split(scaf, len(y0), **cfg)
        mdl, n_it, _chg, _ = em.censored_em(X, y0, kind, tr)
        preds[nm] = mdl.predict(X_A)
        print(f"  {nm:40s} EM {n_it}it  [{time.time()-t0:.0f}s]", flush=True)

    # the selection, reproduced step by step so every stage is inspectable
    cnt = Counter()
    per_conv_top = {}
    for nm, p in preds.items():
        top = np.argsort(-p)[:pc.TOP_ABS]
        ids = A.loc[top, "uid"].tolist()
        per_conv_top[nm] = ids
        cnt.update(ids)
    mean = dict(zip(A["uid"], np.vstack([preds[c] for c in cells]).mean(axis=0)))
    surv = {k: v for k, v in cnt.items() if v >= pc.MAJORITY}
    ranked = sorted(surv, key=lambda c: (-surv[c], -mean[c]))
    info = A.set_index("uid")

    # scaffold uniqueness, recording WHY each drop happened
    shortlist, seen, dropped = [], {}, []
    for c in ranked:
        fw = info.loc[c, "murcko"]
        if fw in seen:
            dropped.append({"uid": c, "conventions": int(surv[c]),
                            "mean_pred_pki": round(float(mean[c]), 4),
                            "murcko": fw,
                            "blocked_by": seen[fw],
                            "blocked_by_conventions": int(surv[seen[fw]]),
                            "blocked_by_mean_pki": round(float(mean[seen[fw]]), 4)})
            continue
        seen[fw] = c
        shortlist.append(c)
        if len(shortlist) == pc.SHORTLIST:
            break

    # names, flags, properties for every survivor
    led = pd.read_csv(SNPLIB)
    namecol = next((c for c in ("name", "compound_name", "preferred_name")
                    if c in led.columns), None)
    names = dict(zip(led["SNP_ID"].astype(str), led[namecol])) if namecol else {}
    cats = {k: catalog(v) for k, v in CATS.items()}

    rows = []
    for i, c in enumerate(ranked, 1):
        s = info.loc[c, "_can"]
        m = Chem.MolFromSmiles(s)
        fl = {k: sorted({x.GetDescription() for x in fc.GetMatches(m)})
              for k, fc in cats.items()} if m else {k: [] for k in CATS}
        rows.append({
            "survivor_rank": i, "SNP_ID": c, "name": names.get(c, ""),
            "conventions_of_8": int(surv[c]),
            "distinct_splits_max": None,
            "mean_pred_pki": round(float(mean[c]), 4),
            "ad_score": round(float(info.loc[c, "ad_score"]), 4),
            "murcko_framework": info.loc[c, "murcko"],
            "in_shortlist": c in shortlist,
            "shortlist_rank": (shortlist.index(c) + 1) if c in shortlist else "",
            "dropped_by_scaffold_uniqueness": any(d["uid"] == c for d in dropped),
            "mol_wt": round(Descriptors.MolWt(m), 1) if m else None,
            "logp": round(Descriptors.MolLogP(m), 2) if m else None,
            "tpsa": round(rdMolDescriptors.CalcTPSA(m), 1) if m else None,
            "qed": round(QED.qed(m), 3) if m else None,
            "n_flags": sum(len(v) for v in fl.values()),
            "PAINS": "; ".join(fl["PAINS"]) or "-",
            "BRENK": "; ".join(fl["BRENK"]) or "-",
            "NIH": "; ".join(fl["NIH"]) or "-",
            "CHEMBL": "; ".join(fl["CHEMBL"]) or "-",
            "smiles": s})
    df = pd.DataFrame(rows)
    df.to_csv(OUTDIR / "survivors.csv", index=False, encoding="utf-8")

    # how many DISTINCT splits back each survivor, given the D3 degeneracy
    sig = {}
    for nm in cells:
        key = tuple(np.round(preds[nm], 10))
        sig.setdefault(key, []).append(nm)
    groups = list(sig.values())
    distinct = {}
    for c in ranked:
        n = sum(1 for g in groups if any(c in per_conv_top[nm] for nm in g))
        distinct[c] = n
    df["distinct_splits_max"] = df["SNP_ID"].map(distinct)
    df.to_csv(OUTDIR / "survivors.csv", index=False, encoding="utf-8")
    for r in rows:                      # the rows dicts back the console table
        r["distinct_splits_max"] = distinct[r["SNP_ID"]]

    unan = [c for c in ranked if surv[c] == len(cells)]
    unan_dropped = [d for d in dropped if d["conventions"] == len(cells)]

    out = {"adr": "ADR-AIIM-312 addendum",
           "log_sha256_at_launch": LOG_AT_LAUNCH,
           "rdkit_version": rdkit.__version__,
           "purpose": ("recover the survivor set the ADR-AIIM-311 harness "
                       "computed but did not save"),
           "anchor": {
               "shortlist_reproduced_in_order": shortlist == prev_sl,
               "n_survivors_reproduced": len(surv) == prev_nsurv,
               "n_unanimous_reproduced": len(unan) == prev_nunan,
               "previous_shortlist": prev_sl, "this_shortlist": shortlist,
               "previous_n_survivors": prev_nsurv, "this_n_survivors": len(surv),
               "previous_n_unanimous": prev_nunan, "this_n_unanimous": len(unan)},
           "pool_size": int(len(A)), "window_k": pc.TOP_ABS,
           "majority": pc.MAJORITY, "shortlist_cap": pc.SHORTLIST,
           "n_survivors": len(surv), "n_unanimous": len(unan),
           "unanimous": unan,
           "shortlist": shortlist,
           "dropped_by_scaffold_uniqueness": dropped,
           "unanimous_dropped_by_scaffold_uniqueness": unan_dropped,
           "distinct_split_groups": [g for g in sig.values()],
           "n_distinct_splits": len(groups),
           "distinct_splits_per_survivor": distinct,
           "seconds": round(time.time() - t0, 1)}
    (OUTDIR / "em_survivors.json").write_text(json.dumps(out, indent=1),
                                              encoding="utf-8")

    print("\n" + "=" * 92)
    a = out["anchor"]
    print(f"ANCHOR  shortlist in order {a['shortlist_reproduced_in_order']}   "
          f"survivors {a['this_n_survivors']}/{a['previous_n_survivors']}   "
          f"unanimous {a['this_n_unanimous']}/{a['previous_n_unanimous']}")
    print("=" * 92)
    print(f"{'#':>2s} {'SNP_ID':10s} {'conv':>4s} {'dist':>4s} {'pKi':>6s} "
          f"{'AD':>6s} {'fl':>3s} {'SL':>3s}  name")
    for r in rows:
        mark = (f"#{r['shortlist_rank']}" if r["in_shortlist"]
                else ("CUT" if r["dropped_by_scaffold_uniqueness"] else "-"))
        print(f"{r['survivor_rank']:2d} {r['SNP_ID']:10s} "
              f"{r['conventions_of_8']:3d}/8 {r['distinct_splits_max']:3d}/"
              f"{len(groups)} {r['mean_pred_pki']:6.3f} {r['ad_score']:6.3f} "
              f"{r['n_flags']:3d} {mark:>3s}  {r['name'][:44]}")
    print(f"\ndistinct split groups: {len(groups)}")
    for g in groups:
        print("   " + " == ".join(g))
    print(f"\ndropped by scaffold uniqueness: {len(dropped)}")
    for d in dropped:
        print(f"   {d['uid']:10s} {d['conventions']}/8  pKi "
              f"{d['mean_pred_pki']:.4f}  framework {d['murcko']}")
        print(f"       blocked by {d['blocked_by']} "
              f"({d['blocked_by_conventions']}/8, "
              f"pKi {d['blocked_by_mean_pki']:.4f}) sharing that framework")
    print(f"\nunanimous ({len(unan)}): {unan}")
    print(f"unanimous dropped: {[d['uid'] for d in unan_dropped]}")
    print(f"\nwritten -> {OUTDIR / 'em_survivors.json'}")
    print(f"           {OUTDIR / 'survivors.csv'}  [{out['seconds']}s]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

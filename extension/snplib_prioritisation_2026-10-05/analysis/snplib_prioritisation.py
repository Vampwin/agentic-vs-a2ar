"""ADR-AIIM-306, executed. Pool, gates, flags and shortlist length frozen
before this ran.

Pre-registration SHA-256
  830e195a1e42f2eeb4d093e1eb167f6591e1007538ea46e7e8ce5d97aa5cb0cc
Decision log at freeze (text/LF, per the 2026-10-05 addendum)
  b2d269a04cd0d4a52194f06e902abc8c6aae8eb8ade96e4ff937a3f0ffc646d1

ADR-AIIM-300 drew its pool from one file -- np_library_dock_ledger.csv -- and
never declared why SNP-Lib was excluded. That was my undeclared harness
decision. The author's ruling settles its cost: SNP-Lib is obtainable and the
NP-library is not, so a LOTUS shortlist cannot be assayed.

This is a SCOPE correction, not a result-driven one: the pool changes because
of what is obtainable, which is independent of what any ranking returns.

GATES:  Lipinski Ro5 (from G1's own ledger column), then AD >= 0.25.
FLAGS:  PAINS, Brenk, NIH, ChEMBL -- reported per compound, REMOVING NOTHING.
        SNP-Lib is 52 % phenanthrenes and 20 % phenolics, catechol matches both
        PAINS and Brenk, and PAINS was derived from HTS of synthetic libraries
        and is known for false positives on natural products. A hard gate would
        remove compounds for a reason that does not apply here; omitting the
        flags would withhold what a reviewer will ask for.

DECLARED DEFECT CARRIED FORWARD: the 0.25 boundary was measured by
ADR-AIIM-296 against a 4,566-compound TRAINING side; this scores against the
5,707-compound curated binding set, as ADR-AIIM-300 did. ADR-AIIM-302 recorded
this reference-size dependence and it is repeated here for comparability,
declared rather than discovered afterwards.

`medchem` SKILL.md was read in full and medchem is NOT invoked -- not installed,
no install requested, RDKit's catalogs used instead (ADR-AIIM-306 constraint 1).

`10_new_system` is not modified.
"""
from __future__ import annotations

import importlib.util
import itertools
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
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

RDLogger.DisableLog("rdApp.*")

ROOT = Path(r"D:\Agentic VS Workflow")
SCRIPTS = ROOT / "10_new_system" / "scripts"
CUR_RUN = ROOT / "20_new_runs" / "chembl251_full_curation_2026-10-05"
REG_RUN = ROOT / "20_new_runs" / "pki_regression_2026-10-05"
EM_RUN = ROOT / "20_new_runs" / "pki_censored_em_2026-10-05"
SNPLIB = ROOT / "20_new_runs" / "snp_screen_2026-08-04" / "G1_standardization" / "snp_source_ledger.csv"
OUTDIR = ROOT / "20_new_runs" / "snplib_prioritisation_2026-10-05"
sys.path.insert(0, str(SCRIPTS))

PREREG = "830e195a1e42f2eeb4d093e1eb167f6591e1007538ea46e7e8ce5d97aa5cb0cc"
LOG_AT_FREEZE = "b2d269a04cd0d4a52194f06e902abc8c6aae8eb8ade96e4ff937a3f0ffc646d1"
RS, TEST_FRAC = 42, 0.20
AD_PRIMARY, AD_SENSITIVITY = 0.25, 0.20
TOP_N, MAJORITY, SHORTLIST = 20, 5, 10
ADR302_TOP_PKI = 6.191          # the S3 reference
S2_FLOOR = 0.30
NPLIB_THREE = ["NPLIB01185", "NPLIB01643", "NPLIB03056"]

_morgan = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
CATALOGS = {
    "PAINS": FilterCatalogParams.FilterCatalogs.PAINS,
    "BRENK": FilterCatalogParams.FilterCatalogs.BRENK,
    "NIH": FilterCatalogParams.FilterCatalogs.NIH,
    "CHEMBL": FilterCatalogParams.FilterCatalogs.CHEMBL,
}


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


CELLS = [{"d1_shuffle": a, "d2_train_first": b, "d3_divert": c}
         for a, b, c in itertools.product((True, False), (False, True), (False, True))]


def build_catalogs():
    out = {}
    for name, cat in CATALOGS.items():
        p = FilterCatalogParams()
        p.AddCatalog(cat)
        out[name] = FilterCatalog(p)
    return out


def alert_flags(mol, cats):
    """Every matching alert description per catalog. Flags only -- nothing is
    removed on the strength of these (ADR-AIIM-306 constraint 5)."""
    res = {}
    for name, fc in cats.items():
        if mol is None:
            res[name] = None
            continue
        ms = fc.GetMatches(mol)
        res[name] = sorted({m.GetDescription() for m in ms})
    return res


def main() -> int:
    (OUTDIR / "analysis").mkdir(parents=True, exist_ok=True)
    qsar = load("_qsar", SCRIPTS / "run_qsar_lbdd_2026-08-07.py")
    cur = load("_cur", CUR_RUN / "analysis" / "chembl251_curate.py")
    rg = load("_reg", REG_RUN / "analysis" / "pki_regression.py")
    em = load("_em", EM_RUN / "analysis" / "pki_censored_em.py")
    t0 = time.time()

    # ---- modelling set, identical to ADR-AIIM-298 / 300 / 306 -------------
    tbl, _ = rg.rebuild_table(cur)
    b = tbl[tbl["endpoint_class"] == "binding"].reset_index(drop=True)
    keep = b["basis"].isin(["point_median", "point_median_gap",
                            "censored_upper", "censored_lower"])
    b = b[keep & np.isfinite(b["px"])].reset_index(drop=True)
    kind = np.where(b["basis"] == "censored_upper", 1,
                    np.where(b["basis"] == "censored_lower", -1, 0)).astype(int)
    msm = b["std_smiles"].tolist()
    y0 = b["px"].to_numpy(dtype=float)
    iks = [qsar._smiles_to_inchikey(s) or f"__N{i}__" for i, s in enumerate(msm)]
    order = sorted(range(len(msm)), key=lambda j: iks[j])
    msm = [msm[i] for i in order]
    y0, kind = y0[order], kind[order]
    scaf = [qsar._murcko(s) for s in msm]
    X = qsar.compute_raw_features(msm)
    ref_fps = [f for f in (_morgan.GetFingerprint(Chem.MolFromSmiles(s))
                           if Chem.MolFromSmiles(s) else None for s in msm)
               if f is not None]
    print(f"modelling set {len(y0)}  AD reference fps {len(ref_fps)}", flush=True)

    # ---- the pool: SNP-Lib, G1's canonical SMILES unchanged --------------
    d = pd.read_csv(SNPLIB)
    d = d[d["smiles_canonical"].notna()].reset_index(drop=True)
    wf = [{"step": 0, "name": "SNP-Lib as standardised by G1", "n": int(len(d))}]

    r5 = d[d["lipinski_ro5_pass"] == True].reset_index(drop=True)  # noqa: E712
    wf.append({"step": 1, "name": "gate: Lipinski Ro5 (G1 ledger column)",
               "n": int(len(r5))})

    ad = np.full(len(r5), np.nan)
    fw = [None] * len(r5)
    mols = []
    for i, s in enumerate(r5["smiles_canonical"]):
        m = Chem.MolFromSmiles(s)
        mols.append(m)
        if m is None:
            continue
        sims = DataStructs.BulkTanimotoSimilarity(_morgan.GetFingerprint(m), ref_fps)
        ad[i] = float(np.mean(sorted(sims, reverse=True)[:5]))
        fw[i] = qsar._murcko(s)
    r5["ad_score"] = ad
    r5["murcko_framework"] = fw

    elig = r5[np.isfinite(r5["ad_score"]) & (r5["ad_score"] >= AD_PRIMARY)].reset_index(drop=True)
    sens = r5[np.isfinite(r5["ad_score"]) & (r5["ad_score"] >= AD_SENSITIVITY)]
    wf.append({"step": 2, "name": f"gate: AD >= {AD_PRIMARY}", "n": int(len(elig)),
               "frameworks": int(elig["murcko_framework"].nunique())})
    wf.append({"step": "2s", "name": f"sensitivity only: AD >= {AD_SENSITIVITY}",
               "n": int(len(sens)),
               "frameworks": int(sens["murcko_framework"].nunique()),
               "note": "reported, not used for selection"})
    for w in wf:
        print(f"  [{w['step']}] {w['name']:44s} n={w['n']:4d}"
              + (f"  frameworks={w['frameworks']}" if "frameworks" in w else ""),
              flush=True)

    # ---- flags, removing nothing -----------------------------------------
    cats = build_catalogs()
    flags = [alert_flags(Chem.MolFromSmiles(s), cats)
             for s in elig["smiles_canonical"]]
    any_alert = np.array([bool(f["PAINS"] or f["BRENK"]) for f in flags])
    print(f"\n  flags (PAINS or Brenk): {int(any_alert.sum())}/{len(elig)} "
          f"({any_alert.mean()*100:.1f} %)  -- nothing removed", flush=True)
    for nm in CATALOGS:
        n = sum(1 for f in flags if f[nm])
        print(f"    {nm:7s} {n:4d}/{len(elig)}", flush=True)

    # ---- rank under the eight conventions --------------------------------
    Xq = qsar.compute_raw_features(elig["smiles_canonical"].tolist())
    per_conv, preds = {}, {}
    for cfg in CELLS:
        nm = em.cell_name(cfg)
        tr, _te = em.split(scaf, len(y0), **cfg)
        mdl, n_it, chg, _ = em.censored_em(X, y0, kind, tr)
        pq = mdl.predict(Xq)
        preds[nm] = pq
        top = np.argsort(-pq)[:TOP_N]
        per_conv[nm] = {"em_iterations": n_it,
                        "em_final_max_change": round(float(chg), 5),
                        "top_ids": elig.loc[top, "SNP_ID"].tolist(),
                        "top_preds": [round(float(pq[i]), 3) for i in top]}
        print(f"  {nm:40s} EM {n_it}it  top1 {elig.loc[top[0],'SNP_ID']} "
              f"{pq[top[0]]:.3f}  [{time.time()-t0:.0f}s]", flush=True)

    cnt = Counter()
    for r in per_conv.values():
        cnt.update(r["top_ids"])
    survivors = {k: v for k, v in cnt.items() if v >= MAJORITY}
    P = np.vstack([preds[em.cell_name(c)] for c in CELLS])
    mean_pred = {elig.loc[i, "SNP_ID"]: float(P[:, i].mean()) for i in range(len(elig))}
    ranked = sorted(survivors, key=lambda c: (-survivors[c], -mean_pred[c]))

    info = elig.set_index("SNP_ID")
    fl = {elig.loc[i, "SNP_ID"]: flags[i] for i in range(len(elig))}
    n_fw_first10 = len({info.loc[c, "murcko_framework"] for c in ranked[:10]})

    shortlist, seen = [], set()
    for c in ranked:
        f = info.loc[c, "murcko_framework"]
        if f in seen:
            continue
        seen.add(f)
        shortlist.append(c)
        if len(shortlist) == SHORTLIST:
            break

    rows = []
    for i, c in enumerate(shortlist, 1):
        rows.append({
            "rank": i, "SNP_ID": c, "name": info.loc[c, "name"],
            "conventions_in_top20": int(survivors[c]),
            "mean_pred_pki": round(mean_pred[c], 3),
            "ad_score": round(float(info.loc[c, "ad_score"]), 4),
            "groups": info.loc[c, "groups"],
            "murcko_framework": info.loc[c, "murcko_framework"],
            "mol_wt": round(float(info.loc[c, "mol_wt"]), 1),
            "logp": round(float(info.loc[c, "logp"]), 2),
            "tpsa": round(float(info.loc[c, "tpsa"]), 1),
            "qed": round(float(info.loc[c, "qed"]), 3),
            "n_stereocenters": int(info.loc[c, "n_stereocenters"]),
            "PAINS": "; ".join(fl[c]["PAINS"]) or "-",
            "BRENK": "; ".join(fl[c]["BRENK"]) or "-",
            "NIH": "; ".join(fl[c]["NIH"]) or "-",
            "CHEMBL": "; ".join(fl[c]["CHEMBL"]) or "-",
            "smiles": info.loc[c, "smiles_canonical"]})
    pd.DataFrame(rows).to_csv(OUTDIR / "shortlist.csv", index=False)

    # reported observation: do the three NPLIB compounds' frameworks occur here?
    nplib = pd.read_csv(ROOT / "20_new_runs" / "np_library_funnel_2026-08-15" /
                        "np_library_dock_ledger.csv").set_index("compound_id")
    three_fw = {c: qsar._murcko(nplib.loc[c, "SMILES"]) for c in NPLIB_THREE
                if c in nplib.index}
    shared = {c: sorted(info.index[info["murcko_framework"] == f].tolist())
              for c, f in three_fw.items()}

    top_pki = rows[0]["mean_pred_pki"] if rows else None
    out = {"adr": "ADR-AIIM-306", "prereg_sha256": PREREG,
           "log_sha256_at_freeze": LOG_AT_FREEZE,
           "rdkit_version": rdkit.__version__,
           "pool": "SNP-Lib (obtainable focus library), replacing the "
                   "NP-library expansion on availability grounds",
           "gates": ["Lipinski Ro5 (G1 ledger)", f"AD >= {AD_PRIMARY}"],
           "flags_not_gates": list(CATALOGS),
           "ad_reference": {
               "n": len(ref_fps),
               "declared_defect": "the 0.25 boundary was measured by "
                                  "ADR-AIIM-296 against a 4,566-compound "
                                  "training side; this scores against the "
                                  "5,707-compound curated binding set, as "
                                  "ADR-AIIM-300 did (ADR-AIIM-302)."},
           "medchem": "SKILL.md read in full; package not installed; not "
                      "invoked; no install requested. Forgone: rule_of_cns, "
                      "NIBR severity, Lilly demerits, ZINC-15 complexity.",
           "waterfall": wf,
           "alerts_on_eligible": {
               "any_pains_or_brenk": int(any_alert.sum()),
               "of": int(len(elig)),
               "frac": round(float(any_alert.mean()), 4),
               "per_catalog": {nm: sum(1 for f in flags if f[nm]) for nm in CATALOGS}},
           "per_convention": per_conv,
           "n_survivors_majority": len(survivors),
           "n_unanimous_all8": sum(1 for v in cnt.values() if v == 8),
           "frameworks_in_first_10_survivors": n_fw_first10,
           "shortlist": rows,
           "nplib_three_framework_overlap": {
               "frameworks": three_fw, "snplib_compounds_sharing": shared},
           "caveats": [
               "A shortlisted compound is NOT a predicted active. ADR-AIIM-299: "
               "conformal count 0 at pKi >= 6.0 under every framing, interval "
               "+/-1.5 log against a 0.5 decision gap.",
               "Alerts are flags. Nothing was removed on their strength.",
               "The top-20 window and the scaffold-uniqueness rule remain "
               "unvalidated against labels; ADR-AIIM-304 closed only the "
               "majority rule.",
               "ADR-AIIM-299's limitations carry over: the censored-EM E-step "
               "is biased upward at censored points and the EM does not reach "
               "its 0.01 tolerance in 10 iterations.",
               "Nothing in this project has been assayed."]}

    out["grading"] = {
        "S1_at_least_10_survivors": {
            "n_survivors": len(survivors), "majority": MAJORITY,
            "verdict": "MET" if len(survivors) >= 10 else "MISSED"},
        "S2_alerts_above_30pct": {
            "frac": round(float(any_alert.mean()), 4), "floor": S2_FLOOR,
            "verdict": "MET" if any_alert.mean() > S2_FLOOR else "MISSED"},
        "S3_top_pki_below_adr302": {
            "top_pki": top_pki, "adr302_reference": ADR302_TOP_PKI,
            "verdict": ("MET" if (top_pki is not None and top_pki < ADR302_TOP_PKI)
                        else "MISSED")},
        "S4_uniqueness_binds": {
            "shortlist_len": len(shortlist), "target": SHORTLIST,
            "frameworks_in_first_10": n_fw_first10,
            "verdict": "MET" if len(shortlist) < SHORTLIST else "MISSED"},
    }
    out["seconds"] = round(time.time() - t0, 1)
    (OUTDIR / "snplib_prioritisation.json").write_text(json.dumps(out, indent=1),
                                                       encoding="utf-8")
    (OUTDIR / "waterfall.json").write_text(json.dumps(
        {"waterfall": wf, "alerts": out["alerts_on_eligible"]}, indent=1),
        encoding="utf-8")

    print("\n" + "=" * 92)
    for t, g in out["grading"].items():
        print(f"{t:34s} {g['verdict']}")
    print("=" * 92)
    print(f"eligible {len(elig)}  survivors(>={MAJORITY}/8) {len(survivors)}  "
          f"unanimous 8/8 {out['n_unanimous_all8']}  shortlist {len(shortlist)}")
    print(f"\nPRIORITY ORDER FOR ASSAY -- obtainable focus library, not predicted actives\n")
    print(f"{'#':>2s} {'SNP_ID':9s} {'conv':>4s} {'pKi':>6s} {'AD':>6s} "
          f"{'MW':>6s} {'logP':>5s} {'TPSA':>6s} {'QED':>5s} {'group':20s} alerts")
    for r in rows:
        al = []
        if r["PAINS"] != "-":
            al.append("PAINS")
        if r["BRENK"] != "-":
            al.append("BRENK")
        if r["NIH"] != "-":
            al.append("NIH")
        print(f"{r['rank']:2d} {r['SNP_ID']:9s} {r['conventions_in_top20']:4d} "
              f"{r['mean_pred_pki']:6.3f} {r['ad_score']:6.4f} {r['mol_wt']:6.1f} "
              f"{r['logp']:5.2f} {r['tpsa']:6.1f} {r['qed']:5.3f} "
              f"{str(r['groups'])[:20]:20s} {','.join(al) or '-'}")
    print(f"\nwritten -> {OUTDIR} [{out['seconds']}s]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

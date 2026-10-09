"""ADR-AIIM-300, executed. The rule, every threshold in it, and the shortlist
length were fixed before any compound identity was read.

Pre-registration SHA-256
  2e269726f50f8bcabb4d8f43e07fe7d06d2d987d0569c8f6bb81ae48997bc030
Decision log at freeze (text/LF, per the 2026-10-05 addendum)
  c5c7c4732e4b33192cbbd718a4b074bea5fa047f539d291e4b31986c4df32901

ADR-AIIM-299 established what this model can and cannot do: it RANKS
(Spearman 0.73-0.84) and it CANNOT assert absolute activity (interval +/-1.5
log against a 0.5 decision gap; conformal count 0 at every threshold under
every convention).

So this produces a PRIORITY ORDER FOR EXPERIMENT, not a set of predicted
actives. The word "active" is applied to nothing here (ADR-AIIM-300
constraint 3).

The rule, fixed in advance:
  model        censored-EM exactly as ADR-AIIM-298 (10 iter, tol 0.01)
  eligible     ECFP4 mean-top-5 similarity to training >= 0.25, the boundary
               ADR-AIIM-296 MEASURED -- used as a competence boundary read off
               performance, not as a similarity admission gate (ADR-AIIM-073)
  rank         predicted pKi descending, within the eligible set
  robustness   in the top 20 under AT LEAST 5 OF THE 8 conventions. Instead of
               choosing a convention -- which ADR-AIIM-286 showed is choosing
               an answer -- a compound must survive a declared majority.
  diversity    at most one compound per Bemis-Murcko framework, higher rank kept
  length       10

NPLIB01185 / NPLIB01643 / NPLIB03056 are reported wherever they fall. Their
presence or absence is an observation, not a criterion (ADR-AIIM-300
constraint 7).

`10_new_system` is not modified.
"""
from __future__ import annotations

import importlib.util
import itertools
import json
import os
import random
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

RDLogger.DisableLog("rdApp.*")

ROOT = Path(r"D:\Agentic VS Workflow")
SCRIPTS = ROOT / "10_new_system" / "scripts"
CUR_RUN = ROOT / "20_new_runs" / "chembl251_full_curation_2026-10-05"
REG_RUN = ROOT / "20_new_runs" / "pki_regression_2026-10-05"
EM_RUN = ROOT / "20_new_runs" / "pki_censored_em_2026-10-05"
NPLIB = ROOT / "20_new_runs" / "np_library_funnel_2026-08-15" / "np_library_dock_ledger.csv"
OUTDIR = ROOT / "20_new_runs" / "rank_prioritisation_2026-10-05"
sys.path.insert(0, str(SCRIPTS))

PREREG = "2e269726f50f8bcabb4d8f43e07fe7d06d2d987d0569c8f6bb81ae48997bc030"
LOG_AT_FREEZE = "c5c7c4732e4b33192cbbd718a4b074bea5fa047f539d291e4b31986c4df32901"
RS, TEST_FRAC = 42, 0.20
AD_FLOOR = 0.25          # ADR-AIIM-296's measured competence boundary
TOP_N = 20               # per-convention rank window
MAJORITY = 5             # of 8 conventions
SHORTLIST = 10
NOVELTY_PROBE = 0.45     # Q3: below the band where the model scores above 0.80
RELEASED_THREE = ["NPLIB01185", "NPLIB01643", "NPLIB03056"]

_morgan = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


CELLS = [{"d1_shuffle": a, "d2_train_first": b, "d3_divert": c}
         for a, b, c in itertools.product((True, False), (False, True), (False, True))]


def main() -> int:
    (OUTDIR / "analysis").mkdir(parents=True, exist_ok=True)
    qsar = load("_qsar", SCRIPTS / "run_qsar_lbdd_2026-08-07.py")
    cur = load("_cur", CUR_RUN / "analysis" / "chembl251_curate.py")
    reg = load("_reg", REG_RUN / "analysis" / "pki_regression.py")
    em = load("_em", EM_RUN / "analysis" / "pki_censored_em.py")
    t0 = time.time()

    # ---- modelling set, identical construction to ADR-AIIM-298 -------------
    tbl, _ = reg.rebuild_table(cur)
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
    print(f"modelling set {len(y0)}: point {(kind==0).sum()} "
          f"upper {(kind==1).sum()} lower {(kind==-1).sum()}", flush=True)

    # ---- library ----------------------------------------------------------
    lib = pd.read_csv(NPLIB)
    lib = lib[lib["SMILES"].apply(lambda s: isinstance(s, str) and bool(s))]
    lib = lib.reset_index(drop=True)
    Xq = qsar.compute_raw_features(lib["SMILES"].tolist())
    print(f"library {len(lib)} compounds", flush=True)

    # ---- competence score: ECFP4 mean top-5 to the FULL modelling set ------
    tfps = [f for f in (_morgan.GetFingerprint(Chem.MolFromSmiles(s))
                        if Chem.MolFromSmiles(s) else None for s in smiles)
            if f is not None]
    ad = np.full(len(lib), np.nan)
    for i, s in enumerate(lib["SMILES"]):
        m = Chem.MolFromSmiles(s)
        if m is None:
            continue
        sims = DataStructs.BulkTanimotoSimilarity(_morgan.GetFingerprint(m), tfps)
        ad[i] = float(np.mean(sorted(sims, reverse=True)[:5]))
    lib["ad_score"] = ad
    eligible = np.isfinite(ad) & (ad >= AD_FLOOR)
    print(f"eligible (AD >= {AD_FLOOR}): {int(eligible.sum())} of {len(lib)} "
          f"({eligible.mean()*100:.1f} %)  [{time.time()-t0:.0f}s]", flush=True)

    # ---- rank under each of the eight conventions -------------------------
    per_conv = {}
    preds = {}
    for cfg in CELLS:
        nm = em.cell_name(cfg)
        tr, _te = em.split(scaf, len(y0), **cfg)
        mdl, n_it, chg, _y = em.censored_em(X, y0, kind, tr)
        pq = mdl.predict(Xq)
        preds[nm] = pq
        idx = np.where(eligible)[0]
        ordered = idx[np.argsort(-pq[idx])]
        top = ordered[:TOP_N]
        per_conv[nm] = {"em_iterations": n_it,
                        "em_final_max_change": round(float(chg), 5),
                        "top_ids": lib.loc[top, "compound_id"].tolist(),
                        "top_preds": [round(float(pq[i]), 3) for i in top]}
        print(f"  {nm:40s} EM {n_it}it  top1 {lib.loc[top[0],'compound_id']} "
              f"{pq[top[0]]:.3f}  [{time.time()-t0:.0f}s]", flush=True)

    # ---- the majority rule ------------------------------------------------
    cnt = Counter()
    for nm, r in per_conv.items():
        cnt.update(r["top_ids"])
    survivors = {cid: n for cid, n in cnt.items() if n >= MAJORITY}
    P = np.vstack([preds[em.cell_name(c)] for c in CELLS])
    mean_pred = {lib.loc[i, "compound_id"]: float(P[:, i].mean())
                 for i in range(len(lib))}
    # Ordering of survivors is not fixed by ADR-AIIM-300; declared here as
    # (conventions desc, mean predicted pKi desc) and reported as an
    # implementation choice the pre-registration did not cover.
    ranked = sorted(survivors, key=lambda c: (-survivors[c], -mean_pred[c]))

    info = lib.set_index("compound_id")
    frameworks = {}
    for cid in ranked:
        frameworks[cid] = qsar._murcko(info.loc[cid, "SMILES"])
    n_frameworks_first10 = len({frameworks[c] for c in ranked[:10]})

    shortlist, seen = [], set()
    for cid in ranked:
        fw = frameworks[cid]
        if fw in seen:
            continue
        seen.add(fw)
        shortlist.append(cid)
        if len(shortlist) == SHORTLIST:
            break

    # intersection across all eight, for Q4
    sets = [set(r["top_ids"]) for r in per_conv.values()]
    inter = set.intersection(*sets) if sets else set()

    rows = []
    for rank, cid in enumerate(shortlist, 1):
        rows.append({"rank": rank, "compound_id": cid,
                     "conventions_in_top20": int(survivors[cid]),
                     "mean_pred_pki": round(mean_pred[cid], 3),
                     "ad_score": round(float(info.loc[cid, "ad_score"]), 4),
                     "murcko_framework": frameworks[cid],
                     "npclassifier_class": info.loc[cid, "structure_taxonomy_npclassifier_03class"],
                     "organism": info.loc[cid, "organism_name"],
                     "mw": round(float(info.loc[cid, "mw"]), 2),
                     "funnel_rank": int(info.loc[cid, "funnel_rank"])})
    sl = pd.DataFrame(rows)
    sl.to_csv(OUTDIR / "shortlist.csv", index=False)

    three = {}
    for cid in RELEASED_THREE:
        if cid not in info.index:
            three[cid] = {"present_in_library": False}
            continue
        pos = {}
        for nm, r in per_conv.items():
            pos[nm] = (r["top_ids"].index(cid) + 1) if cid in r["top_ids"] else None
        three[cid] = {"present_in_library": True,
                      "ad_score": round(float(info.loc[cid, "ad_score"]), 4),
                      "eligible": bool(info.loc[cid, "ad_score"] >= AD_FLOOR),
                      "mean_pred_pki": round(mean_pred[cid], 3),
                      "conventions_in_top20": int(cnt.get(cid, 0)),
                      "in_shortlist": cid in shortlist,
                      "rank_per_convention": pos}

    out = {"adr": "ADR-AIIM-300", "prereg_sha256": PREREG,
           "log_sha256_at_freeze": LOG_AT_FREEZE,
           "rdkit_version": rdkit.__version__,
           "rule": {"ad_floor": AD_FLOOR, "top_n": TOP_N,
                    "majority_of_8": MAJORITY, "shortlist": SHORTLIST,
                    "survivor_ordering": "conventions desc, then mean predicted "
                                         "pKi desc -- implementation choice, not "
                                         "fixed by ADR-AIIM-300"},
           "eligible": {"n": int(eligible.sum()), "of": int(len(lib)),
                        "frac": round(float(eligible.mean()), 4)},
           "per_convention": per_conv,
           "n_survivors_majority": len(survivors),
           "n_intersection_all8": len(inter),
           "frameworks_in_first_10_survivors": n_frameworks_first10,
           "shortlist": rows,
           "released_three": three,
           "caveats": [
               "A shortlisted compound is NOT a predicted active. ADR-AIIM-299: "
               "no compound in this library can be asserted at pKi >= 6.0 under "
               "any framing (conformal count 0, interval +/-1.5 log).",
               "The competence band was measured on ChEMBL compounds; carrying "
               "it to natural products assumes degradation tracks distance "
               "rather than chemotype (ADR-AIIM-296 constraint 4).",
               "ADR-AIIM-299's limitations carry over: the censored-EM E-step is "
               "biased upward at censored points and the EM did not reach its "
               "0.01 tolerance inside 10 iterations.",
               "Exchangeability is not established, so no conformal statement is "
               "made about any shortlisted compound.",
               "Nothing in this project has been assayed."]}

    out["grading"] = {
        "Q1_at_least_10_survivors": {
            "n_survivors": len(survivors), "majority": MAJORITY,
            "verdict": "MET" if len(survivors) >= 10 else "MISSED"},
        "Q2_at_least_5_frameworks_in_first_10": {
            "n_frameworks": n_frameworks_first10,
            "verdict": "MET" if n_frameworks_first10 >= 5 else "MISSED"},
        "Q3_at_least_one_below_AD_0.45": {
            "min_ad_in_shortlist": round(float(sl["ad_score"].min()), 4) if len(sl) else None,
            "n_below": int((sl["ad_score"] < NOVELTY_PROBE).sum()) if len(sl) else 0,
            "verdict": "MET" if len(sl) and (sl["ad_score"] < NOVELTY_PROBE).any()
            else "MISSED"},
        "Q4_majority_beats_intersection": {
            "n_majority": len(survivors), "n_intersection_all8": len(inter),
            "verdict": "MET" if len(survivors) > len(inter) else "MISSED"},
    }
    out["seconds"] = round(time.time() - t0, 1)
    (OUTDIR / "rank_prioritisation.json").write_text(json.dumps(out, indent=1),
                                                     encoding="utf-8")

    print("\n" + "=" * 86)
    for t, g in out["grading"].items():
        print(f"{t:40s} {g['verdict']}")
    print("=" * 86)
    print(f"eligible {int(eligible.sum())}/{len(lib)}   survivors (>={MAJORITY}/8) "
          f"{len(survivors)}   intersection all 8: {len(inter)}")
    print(f"\nPRIORITY ORDER FOR EXPERIMENT -- not predicted actives\n")
    print(f"{'#':>2s} {'compound':11s} {'conv':>4s} {'pKi':>6s} {'AD':>6s} "
          f"{'MW':>7s} {'class':34s} organism")
    for r in rows:
        print(f"{r['rank']:2d} {r['compound_id']:11s} {r['conventions_in_top20']:4d} "
              f"{r['mean_pred_pki']:6.3f} {r['ad_score']:6.4f} {r['mw']:7.1f} "
              f"{str(r['npclassifier_class'])[:34]:34s} {str(r['organism'])[:28]}")
    print("\nthe three originally prioritised compounds, reported not scored:")
    for cid, v in three.items():
        if not v.get("present_in_library"):
            print(f"  {cid}: not in library")
            continue
        print(f"  {cid}: AD {v['ad_score']:.4f}  eligible {v['eligible']}  "
              f"mean pKi {v['mean_pred_pki']:.3f}  in top20 of "
              f"{v['conventions_in_top20']}/8  shortlisted {v['in_shortlist']}")
    print(f"\nwritten -> {OUTDIR} [{out['seconds']}s]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

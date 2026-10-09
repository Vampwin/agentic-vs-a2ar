"""ADR-AIIM-301, executed. Rule and comparison frozen before the run, and
before ADR-AIIM-300's shortlist was read.

Pre-registration SHA-256
  50ecb09367a619bfc14b84022397bda00800810bd0f8d47d09fd39c3c2b133ef
Decision log at freeze (text/LF, per the 2026-10-05 addendum)
  b8bf35ef5f90db0c2f2cb76d52addfd1fc76917d4400104a9061280d7f96aed7

ADR-AIIM-299: no compound in this library can be asserted active -- interval
+/-1.5 log against a 0.5 decision gap, conformal count 0 everywhere. Confidence
cannot come from tightening one model. It can only come from requiring
agreement between channels that do not share their evidence.

The project already measured the AVERAGING version and it failed: equal-weight
z-score ensemble AUC 0.530, worse than Asn253 H-bond + Vina at 0.577. SET
INTERSECTION is a different operation with a different operating characteristic
and has never been measured here. R1 is where it can fail.

TWO CHANNELS ONLY:
  L  censored-EM predicted pKi      -- ligand structures (ADR-AIIM-298)
  S  Vina docking score             -- receptor structure and pose
The funnel's 2D/3D pharmacophore scores are LIGAND SIMILARITY and would
double-count L; including them would manufacture agreement rather than test it.
p_active_calibrated / prediction_set_sig0.05 are excluded: ADR-AIIM-286
measured that column spanning 376x.

VALIDATED ON LABELS FIRST. DUD-E carries 482 actives and 1,919
property-matched decoys with production-protocol Vina scores. The NP library
has no labels and never will inside this project.

THE COMPARISON IS AT MATCHED SELECTION SIZE. An intersection selects fewer
compounds than either channel, so comparing it against a channel's full top-X%
would flatter it for free.

APPLICATION IS CONDITIONAL: if R1 fails, NO intersection shortlist is produced
(ADR-AIIM-301, "Application to the library").

`10_new_system` is not modified; no re-docking; Vina scores read from frozen
artifacts.
"""
from __future__ import annotations

import importlib.util
import itertools
import json
import math
import os
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "3")

import numpy as np
import pandas as pd
import rdkit
from rdkit import RDLogger
from sklearn.metrics import roc_auc_score

RDLogger.DisableLog("rdApp.*")

ROOT = Path(r"D:\Agentic VS Workflow")
SCRIPTS = ROOT / "10_new_system" / "scripts"
CUR_RUN = ROOT / "20_new_runs" / "chembl251_full_curation_2026-10-05"
REG_RUN = ROOT / "20_new_runs" / "pki_regression_2026-10-05"
EM_RUN = ROOT / "20_new_runs" / "pki_censored_em_2026-10-05"
RANK_RUN = ROOT / "20_new_runs" / "rank_prioritisation_2026-10-05"
DUDE = ROOT / "20_new_runs" / "dude_enrichment_2026-09-22" / "analysis" / "dude_enrichment_per_compound.csv"
E8 = ROOT / "20_new_runs" / "np_library_funnel_2026-08-15" / "e8_np_library_prioritization.csv"
OUTDIR = ROOT / "20_new_runs" / "orthogonal_intersection_2026-10-05"
sys.path.insert(0, str(SCRIPTS))

PREREG = "50ecb09367a619bfc14b84022397bda00800810bd0f8d47d09fd39c3c2b133ef"
LOG_AT_FREEZE = "b8bf35ef5f90db0c2f2cb76d52addfd1fc76917d4400104a9061280d7f96aed7"
FRACTIONS = [0.01, 0.05, 0.10, 0.20]
PRIMARY_X = 0.10
ENSEMBLE_REFERENCE = 0.530      # section 3.4, equal-weight z-score ensemble
VINA_ANCHOR_MAX = 0.70          # R2


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


CELLS = [{"d1_shuffle": a, "d2_train_first": b, "d3_divert": c}
         for a, b, c in itertools.product((True, False), (False, True), (False, True))]


def sel_metrics(mask, y):
    n = int(mask.sum())
    if n == 0:
        return {"n": 0, "precision": None, "recall": None, "enrichment": None}
    hits = int(y[mask].sum())
    base = float(y.mean())
    prec = hits / n
    return {"n": n, "hits": hits, "precision": round(prec, 4),
            "recall": round(hits / max(int(y.sum()), 1), 4),
            "enrichment": round(prec / base, 3) if base else None}


def topk_mask(score, k):
    """True for the k highest scores. Ties broken by index, deterministically."""
    m = np.zeros(len(score), dtype=bool)
    if k > 0:
        m[np.argsort(-score, kind="stable")[:k]] = True
    return m


def main() -> int:
    (OUTDIR / "analysis").mkdir(parents=True, exist_ok=True)
    qsar = load("_qsar", SCRIPTS / "run_qsar_lbdd_2026-08-07.py")
    cur = load("_cur", CUR_RUN / "analysis" / "chembl251_curate.py")
    rg = load("_reg", REG_RUN / "analysis" / "pki_regression.py")
    em = load("_em", EM_RUN / "analysis" / "pki_censored_em.py")
    t0 = time.time()

    # ---- modelling set, identical to ADR-AIIM-298 / 300 --------------------
    tbl, _ = rg.rebuild_table(cur)
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
    model_ik = {ikeys[i] for i in order}
    scaf = [qsar._murcko(s) for s in smiles]
    X = qsar.compute_raw_features(smiles)
    print(f"modelling set {len(y0)}  ({len(model_ik)} distinct InChIKeys)",
          flush=True)

    # ---- DUD-E ------------------------------------------------------------
    du = pd.read_csv(DUDE)
    du = du[du["SMILES"].notna() & du["vina_score"].notna()].reset_index(drop=True)
    n_in = len(du)
    # leakage control: standardise the same way the curation does, then InChIKey
    std = [cur.standardize(s) for s in du["SMILES"]]
    du["std_smiles"] = std
    du = du[du["std_smiles"].notna()].reset_index(drop=True)
    du["ik"] = [qsar._smiles_to_inchikey(s) for s in du["std_smiles"]]
    leak = du["ik"].isin(model_ik)
    removed = int(leak.sum())
    du_before = du["label"].value_counts().to_dict()
    du = du[~leak].reset_index(drop=True)
    y = du["label"].to_numpy(dtype=int)
    print(f"DUD-E {n_in} -> parsed {n_in - (n_in - len(std))}, "
          f"leakage removed {removed}, surviving {len(du)}: "
          f"{int(y.sum())} active / {int((y==0).sum())} decoy "
          f"({y.mean()*100:.1f} % active)", flush=True)

    Xd = qsar.compute_raw_features(du["std_smiles"].tolist())

    # ---- channel L: mean censored-EM prediction across the eight conventions
    # ADR-AIIM-301 leaves this open between one convention and the mean; the
    # mean is used because it is the quantity ADR-AIIM-300 ranks on. Declared.
    P = []
    for cfg in CELLS:
        tr, _te = em.split(scaf, len(y0), **cfg)
        mdl, n_it, chg, _ = em.censored_em(X, y0, kind, tr)
        P.append(mdl.predict(Xd))
        print(f"  {em.cell_name(cfg):40s} EM {n_it}it  [{time.time()-t0:.0f}s]",
              flush=True)
    L = np.vstack(P).mean(axis=0)
    S = -du["vina_score"].to_numpy(dtype=float)      # more negative = better

    auc_L = float(roc_auc_score(y, L))
    auc_S = float(roc_auc_score(y, S))
    print(f"\n  AUC  L {auc_L:.4f}   S {auc_S:.4f}   "
          f"(section 3.4 ensemble reference {ENSEMBLE_REFERENCE})", flush=True)

    frac_out = {}
    for X_frac in FRACTIONS:
        k = max(int(math.ceil(X_frac * len(du))), 1)
        mL, mS = topk_mask(L, k), topk_mask(S, k)
        mI = mL & mS
        nI = int(mI.sum())
        row = {"fraction": X_frac, "k_per_channel": k,
               "L_at_k": sel_metrics(mL, y), "S_at_k": sel_metrics(mS, y),
               "intersection": sel_metrics(mI, y),
               "matched_L_at_nI": sel_metrics(topk_mask(L, nI), y) if nI else None,
               "matched_S_at_nI": sel_metrics(topk_mask(S, nI), y) if nI else None}
        frac_out[str(X_frac)] = row
        i, ml, ms = row["intersection"], row["matched_L_at_nI"], row["matched_S_at_nI"]
        print(f"  X={X_frac:.0%}  k={k:4d}  I n={i['n']:3d} prec={i['precision']}  "
              f"| matched L prec={ml['precision'] if ml else None}  "
              f"S prec={ms['precision'] if ms else None}", flush=True)

    pr = frac_out[str(PRIMARY_X)]
    I, ML, MS = pr["intersection"], pr["matched_L_at_nI"], pr["matched_S_at_nI"]
    r1 = (I["precision"] is not None and ML and MS
          and I["precision"] > ML["precision"] and I["precision"] > MS["precision"])

    out = {"adr": "ADR-AIIM-301", "prereg_sha256": PREREG,
           "log_sha256_at_freeze": LOG_AT_FREEZE,
           "rdkit_version": rdkit.__version__,
           "channels": {"L": "censored-EM predicted pKi, mean of 8 conventions "
                             "(implementation choice, ADR-AIIM-301 left it open)",
                        "S": "Vina score from frozen artifacts, sign flipped",
                        "excluded": ["2D pharmacophore (ligand similarity, "
                                     "double-counts L)",
                                     "3D shape/pharmacophore (same)",
                                     "p_active_calibrated / prediction_set "
                                     "(ADR-AIIM-286: spans 376x)"]},
           "validation_set": {
               "source": str(DUDE.relative_to(ROOT)),
               "rows_in": n_in, "leakage_removed": removed,
               "composition_before_removal": {str(k): int(v)
                                              for k, v in du_before.items()},
               "surviving": int(len(du)), "actives": int(y.sum()),
               "decoys": int((y == 0).sum()),
               "active_rate": round(float(y.mean()), 4),
               "note": "surviving actives are those this project's curation did "
                       "not capture, and may be a biased subset (ADR-AIIM-301)"},
           "auc": {"L": round(auc_L, 4), "S": round(auc_S, 4),
                   "section_3_4_ensemble_reference": ENSEMBLE_REFERENCE},
           "by_fraction": frac_out,
           "primary_fraction": PRIMARY_X}

    out["grading"] = {
        "R1_intersection_beats_both_at_matched_n": {
            "intersection_precision": I["precision"], "n": I["n"],
            "matched_L_precision": ML["precision"] if ML else None,
            "matched_S_precision": MS["precision"] if MS else None,
            "verdict": "MET" if r1 else "MISSED"},
        "R2_vina_auc_below_0.70": {
            "auc_S": round(auc_S, 4), "threshold": VINA_ANCHOR_MAX,
            "verdict": "MET" if auc_S < VINA_ANCHOR_MAX else "MISSED"},
        "R3_L_auc_exceeds_S": {
            "auc_L": round(auc_L, 4), "auc_S": round(auc_S, 4),
            "verdict": "MET" if auc_L > auc_S else "MISSED"},
        "R4_intersection_recall_below_50pct": {
            "recall": I["recall"],
            "verdict": "MET" if (I["recall"] is not None and I["recall"] < 0.50)
            else "MISSED"},
    }

    # ---- application, CONDITIONAL on R1 -----------------------------------
    app = {"applied": False,
           "reason": "R1 MISSED -- ADR-AIIM-301 forbids producing an "
                     "intersection shortlist in that case"}
    if r1:
        rank_json = RANK_RUN / "rank_prioritisation.json"
        if not rank_json.exists():
            app = {"applied": False,
                   "reason": "R1 MET but ADR-AIIM-300 output not yet on disk; "
                             "application deferred rather than improvised"}
        else:
            rj = json.loads(rank_json.read_text(encoding="utf-8"))
            robust = {c for r in rj["per_convention"].values()
                      for c in r["top_ids"]}
            e8 = pd.read_csv(E8)
            e8 = e8[e8["vina_score"].notna()].reset_index(drop=True)
            k8 = max(int(math.ceil(PRIMARY_X * len(e8))), 1)
            s8 = topk_mask(-e8["vina_score"].to_numpy(dtype=float), k8)
            vina_top = set(e8.loc[s8, "compound_id"])
            survivors = {c: n for c, n in
                         ((c, sum(1 for r in rj["per_convention"].values()
                                  if c in r["top_ids"])) for c in robust)
                         if n >= rj["rule"]["majority_of_8"]}
            inter = sorted(set(survivors) & vina_top,
                           key=lambda c: -survivors[c])
            rows = []
            info = e8.set_index("compound_id")
            for i, cid in enumerate(inter, 1):
                rows.append({"rank": i, "compound_id": cid,
                             "conventions_in_top20": survivors[cid],
                             "vina_score": float(info.loc[cid, "vina_score"]),
                             "on_pocket": bool(info.loc[cid, "on_pocket"]),
                             "npclassifier_class": info.loc[cid, "structure_taxonomy_npclassifier_03class"],
                             "organism": info.loc[cid, "organism_name"]})
            if rows:
                pd.DataFrame(rows).to_csv(OUTDIR / "intersection_shortlist.csv",
                                          index=False)
            app = {"applied": True, "e8_docked": int(len(e8)),
                   "vina_top_k": k8, "rank_survivors": len(survivors),
                   "intersection": rows,
                   "scope_limit": "covers only the docked compounds, never the "
                                  "full 5,000 (ADR-AIIM-301 constraint 4)"}
    out["application"] = app

    out["caveats"] = [
        "Nothing selected here is a predicted active. ADR-AIIM-299: this "
        "instrument cannot resolve pKi >= 6.0 for any compound in this library.",
        "Precision measured on property-matched decoys transfers to natural "
        "products only if both channels behave the same way on both, which is "
        "not tested (ADR-AIIM-301 constraint 5).",
        "ADR-AIIM-299's limitations carry over: the censored-EM E-step is "
        "biased upward at censored points and the EM did not reach tolerance.",
        "Exchangeability is not established; no conformal statement is made.",
        "Nothing in this project has been assayed."]
    out["seconds"] = round(time.time() - t0, 1)
    (OUTDIR / "intersection_validation.json").write_text(
        json.dumps(out, indent=1), encoding="utf-8")

    print("\n" + "=" * 80)
    for t, g in out["grading"].items():
        print(f"{t:46s} {g['verdict']}")
    print("=" * 80)
    print(f"surviving labelled set: {len(du)}  ({int(y.sum())} active, "
          f"{y.mean()*100:.1f} %)   leakage removed {removed}")
    print(f"AUC   L {auc_L:.4f}   S {auc_S:.4f}")
    print(f"\nat X = {PRIMARY_X:.0%}, k = {pr['k_per_channel']} per channel:")
    print(f"  intersection   n={I['n']:3d}  precision {I['precision']}  "
          f"recall {I['recall']}  enrichment {I['enrichment']}")
    if ML and MS:
        print(f"  matched L      n={ML['n']:3d}  precision {ML['precision']}")
        print(f"  matched S      n={MS['n']:3d}  precision {MS['precision']}")
    print(f"\napplication: {'YES' if app.get('applied') else 'NO'} -- "
          f"{app.get('reason', app.get('scope_limit',''))}")
    if app.get("applied"):
        for r in app["intersection"]:
            print(f"  {r['rank']:2d} {r['compound_id']:11s} conv "
                  f"{r['conventions_in_top20']}/8  vina {r['vina_score']:7.3f}  "
                  f"{str(r['npclassifier_class'])[:30]}")
    print(f"\nwritten -> {OUTDIR} [{out['seconds']}s]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

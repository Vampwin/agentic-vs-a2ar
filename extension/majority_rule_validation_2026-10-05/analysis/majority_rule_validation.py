"""ADR-AIIM-303, executed. Predictions frozen before this ran.

Pre-registration SHA-256
  81011897d6976ee0f17623812675e92082b83858d9378a9dbe31d9b370ab660f
Decision log at freeze (text/LF, per the 2026-10-05 addendum)
  561f84e9945fe4019be4be333b03872a7407bb2a21cc68a3f7201a81ba23e8e8

ADR-AIIM-302's audit table: of the five criteria that produced the
six-compound priority order, TWO were validated against labelled actives and
inactives and THREE were not. The 5-of-8 majority rule is one of the three.

ADR-AIIM-300 defended it on DETERMINACY -- ADR-AIIM-286 had shown that
choosing a convention is choosing an answer -- and Q4 confirmed it yields
something where unanimity yields nothing (10 survivors against 0). But YIELD
IS NOT ENRICHMENT. This tests whether majority survivors are more likely to be
active than the compounds a single convention would have picked.

M and A are two aggregations of the SAME eight models, not independent
channels: M counts top-K memberships and discards magnitude, A averages
magnitudes and discards rank. This tests which aggregation is better, not
whether aggregation helps.

EVERY COMPARISON IS AT MATCHED SELECTION SIZE n = |M|.

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
DUDE = ROOT / "20_new_runs" / "dude_enrichment_2026-09-22" / "analysis" / "dude_enrichment_per_compound.csv"
OUTDIR = ROOT / "20_new_runs" / "majority_rule_validation_2026-10-05"
sys.path.insert(0, str(SCRIPTS))

PREREG = "81011897d6976ee0f17623812675e92082b83858d9378a9dbe31d9b370ab660f"
LOG_AT_FREEZE = "561f84e9945fe4019be4be333b03872a7407bb2a21cc68a3f7201a81ba23e8e8"
WINDOWS = [0.01, 0.05, 0.10, 0.20, 0.30]
PRIMARY_K = 0.10
MAJORITY = 5          # of 8, exactly as ADR-AIIM-300
V4_FLOOR = 0.50


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


CELLS = [{"d1_shuffle": a, "d2_train_first": b, "d3_divert": c}
         for a, b, c in itertools.product((True, False), (False, True), (False, True))]


def topk_mask(score, k):
    m = np.zeros(len(score), dtype=bool)
    if k > 0:
        m[np.argsort(-score, kind="stable")[:k]] = True
    return m


def metrics(mask, y):
    n = int(mask.sum())
    if n == 0:
        return {"n": 0, "hits": 0, "precision": None, "recall": None,
                "enrichment": None}
    hits = int(y[mask].sum())
    base = float(y.mean())
    prec = hits / n
    return {"n": n, "hits": hits, "precision": round(prec, 4),
            "recall": round(hits / max(int(y.sum()), 1), 4),
            "enrichment": round(prec / base, 3) if base else None}


def main() -> int:
    (OUTDIR / "analysis").mkdir(parents=True, exist_ok=True)
    qsar = load("_qsar", SCRIPTS / "run_qsar_lbdd_2026-08-07.py")
    cur = load("_cur", CUR_RUN / "analysis" / "chembl251_curate.py")
    rg = load("_reg", REG_RUN / "analysis" / "pki_regression.py")
    em = load("_em", EM_RUN / "analysis" / "pki_censored_em.py")
    t0 = time.time()

    # ---- modelling set, identical to ADR-AIIM-298 / 300 / 301 -------------
    tbl, _ = rg.rebuild_table(cur)
    b = tbl[tbl["endpoint_class"] == "binding"].reset_index(drop=True)
    keep = b["basis"].isin(["point_median", "point_median_gap",
                            "censored_upper", "censored_lower"])
    b = b[keep & np.isfinite(b["px"])].reset_index(drop=True)
    kind = np.where(b["basis"] == "censored_upper", 1,
                    np.where(b["basis"] == "censored_lower", -1, 0)).astype(int)
    smiles = b["std_smiles"].tolist()
    y0 = b["px"].to_numpy(dtype=float)
    iks = [qsar._smiles_to_inchikey(s) or f"__N{i}__" for i, s in enumerate(smiles)]
    order = sorted(range(len(smiles)), key=lambda j: iks[j])
    smiles = [smiles[i] for i in order]
    y0, kind = y0[order], kind[order]
    model_ik = {iks[i] for i in order}
    scaf = [qsar._murcko(s) for s in smiles]
    X = qsar.compute_raw_features(smiles)

    # ---- DUD-E with the ADR-AIIM-301 leakage control ----------------------
    du = pd.read_csv(DUDE)
    du = du[du["SMILES"].notna() & du["vina_score"].notna()].reset_index(drop=True)
    n_in = len(du)
    du["std_smiles"] = [cur.standardize(s) for s in du["SMILES"]]
    du = du[du["std_smiles"].notna()].reset_index(drop=True)
    du["ik"] = [qsar._smiles_to_inchikey(s) for s in du["std_smiles"]]
    leak = du["ik"].isin(model_ik)
    removed = int(leak.sum())
    before = du["label"].value_counts().to_dict()
    du = du[~leak].reset_index(drop=True)
    y = du["label"].to_numpy(dtype=int)
    print(f"DUD-E {n_in} -> leakage removed {removed} -> surviving {len(du)}: "
          f"{int(y.sum())} active / {int((y==0).sum())} decoy "
          f"({y.mean()*100:.1f} %)", flush=True)
    Xd = qsar.compute_raw_features(du["std_smiles"].tolist())

    # ---- the eight conventions -------------------------------------------
    preds, names = [], []
    for cfg in CELLS:
        nm = em.cell_name(cfg)
        tr, _te = em.split(scaf, len(y0), **cfg)
        mdl, n_it, _chg, _ = em.censored_em(X, y0, kind, tr)
        preds.append(mdl.predict(Xd))
        names.append(nm)
        print(f"  {nm:40s} EM {n_it}it  [{time.time()-t0:.0f}s]", flush=True)
    P = np.vstack(preds)
    A_score = P.mean(axis=0)

    auc = {"mean_A": round(float(roc_auc_score(y, A_score)), 4)}
    for i, nm in enumerate(names):
        auc[nm] = round(float(roc_auc_score(y, P[i])), 4)

    out = {"adr": "ADR-AIIM-303", "prereg_sha256": PREREG,
           "log_sha256_at_freeze": LOG_AT_FREEZE,
           "rdkit_version": rdkit.__version__,
           "majority_of_8": MAJORITY, "windows": WINDOWS,
           "primary_window": PRIMARY_K,
           "validation_set": {
               "rows_in": n_in, "leakage_removed": removed,
               "composition_before_removal": {str(k): int(v) for k, v in before.items()},
               "surviving": int(len(du)), "actives": int(y.sum()),
               "decoys": int((y == 0).sum()),
               "active_rate": round(float(y.mean()), 4),
               "note": "all removals were actives; decoys remain matched to the "
                       "original 482, which loosens the matching "
                       "(ADR-AIIM-303 constraint 4)"},
           "auc": auc, "by_window": {}}

    for K in WINDOWS:
        k = max(int(math.ceil(K * len(du))), 1)
        masks = [topk_mask(P[i], k) for i in range(len(names))]
        votes = np.sum(masks, axis=0)
        M = votes >= MAJORITY
        U = votes == len(names)
        nM = int(M.sum())
        mM = metrics(M, y)
        mA = metrics(topk_mask(A_score, nM), y) if nM else None
        singles = {names[i]: metrics(topk_mask(P[i], nM), y) for i in range(len(names))} if nM else {}
        sp = sorted(v["precision"] for v in singles.values()) if singles else []
        row = {"window": K, "k_per_convention": k,
               "n_majority": nM, "n_unanimous": int(U.sum()),
               "majority": mM, "mean_rank_A_at_nM": mA,
               "singles_at_nM": singles,
               "single_precision_min": sp[0] if sp else None,
               "single_precision_median": round(float(np.median(sp)), 4) if sp else None,
               "single_precision_max": sp[-1] if sp else None,
               "unanimous": metrics(U, y)}
        out["by_window"][str(K)] = row
        print(f"  K={K:.0%} k={k:4d}  |M|={nM:4d} |U|={int(U.sum()):4d}  "
              f"M prec={mM['precision']}  A prec={mA['precision'] if mA else None}  "
              f"single med={row['single_precision_median']}  "
              f"(min {row['single_precision_min']} max {row['single_precision_max']})",
              flush=True)

    pr = out["by_window"][str(PRIMARY_K)]
    Mp = pr["majority"]["precision"]
    Ap = pr["mean_rank_A_at_nM"]["precision"] if pr["mean_rank_A_at_nM"] else None
    med = pr["single_precision_median"]
    nests = all(out["by_window"][str(K)]["n_unanimous"]
                <= out["by_window"][str(K)]["n_majority"]
                <= out["by_window"][str(K)]["k_per_convention"] for K in WINDOWS)

    out["grading"] = {
        "V1_majority_beats_median_single": {
            "majority_precision": Mp, "single_median": med,
            "verdict": "MET" if (Mp is not None and med is not None and Mp > med)
            else "MISSED"},
        "V2_majority_does_not_beat_mean_rank": {
            "majority_precision": Mp, "mean_rank_precision": Ap,
            "verdict": "MET" if (Mp is not None and Ap is not None and Mp <= Ap)
            else "MISSED"},
        "V3_sets_nest": {"nests_at_every_window": bool(nests),
                         "verdict": "MET" if nests else "MISSED"},
        "V4_majority_precision_above_0.50": {
            "majority_precision": Mp, "base_rate": round(float(y.mean()), 4),
            "verdict": "MET" if (Mp is not None and Mp > V4_FLOOR) else "MISSED"},
    }
    out["caveats"] = [
        "This validates a selection rule, not the compounds. No compound is "
        "reclassified and nothing is called active (ADR-AIIM-303 constraint 3).",
        "M and A aggregate the SAME eight models; neither is an independent "
        "channel, so this compares aggregations rather than testing whether "
        "aggregation helps.",
        "Only the majority rule is tested. The top-20 window and the "
        "scaffold-uniqueness rule remain unvalidated.",
        "ADR-AIIM-299's limitations carry over: the censored-EM E-step is "
        "biased upward at censored points and the EM does not reach tolerance.",
        "Nothing in this project has been assayed."]
    out["seconds"] = round(time.time() - t0, 1)
    (OUTDIR / "majority_rule_validation.json").write_text(
        json.dumps(out, indent=1), encoding="utf-8")

    print("\n" + "=" * 82)
    for t, g in out["grading"].items():
        print(f"{t:44s} {g['verdict']}")
    print("=" * 82)
    print(f"surviving {len(du)}  ({int(y.sum())} active, {y.mean()*100:.1f} %)  "
          f"leakage removed {removed}")
    print(f"AUC  mean-A {auc['mean_A']}   single range "
          f"{min(v for k,v in auc.items() if k!='mean_A')}-"
          f"{max(v for k,v in auc.items() if k!='mean_A')}")
    print(f"\n{'K':>5s} {'k':>5s} {'|M|':>5s} {'|U|':>5s} {'M prec':>8s} "
          f"{'A prec':>8s} {'single med':>11s} {'min':>7s} {'max':>7s}")
    for K in WINDOWS:
        r = out["by_window"][str(K)]
        a = r["mean_rank_A_at_nM"]
        print(f"{K:5.0%} {r['k_per_convention']:5d} {r['n_majority']:5d} "
              f"{r['n_unanimous']:5d} {str(r['majority']['precision']):>8s} "
              f"{str(a['precision'] if a else None):>8s} "
              f"{str(r['single_precision_median']):>11s} "
              f"{str(r['single_precision_min']):>7s} "
              f"{str(r['single_precision_max']):>7s}")
    print(f"\nwritten -> {OUTDIR} [{out['seconds']}s]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

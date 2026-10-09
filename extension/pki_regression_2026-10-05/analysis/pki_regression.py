"""ADR-AIIM-297, executed. Predictions frozen before this ran.

Pre-registration SHA-256
  b46f97a981c575287f5faffc414571e1b96a976895d40deb9a41b25934ca569c
Decision log at freeze (text/LF, per the 2026-10-05 addendum)
  a36db9f9f3bcd8df25ce239e93e3aa163d02bd7d03978d1074837cde1e702e58

Models pKi as a continuous quantity and thresholds afterwards, which pays none
of the three costs the binary framing does: the 1,252 rows discarded to the
5.5-6.0 gap, the 914 censored compounds flattened to a hard 0/1, and a cut
that has to be chosen before training.

It also removes the Mondrian conformal machinery whose calibration-actives
composition ADR-AIIM-286 showed swinging the count by two orders of magnitude.

THE CENSORED RECORDS ARE A TEST SET, NOT TRAINING DATA. 617 upper- and 297
lower-censored compounds are excluded from fitting and from the held-out
regression metrics, then checked against the single inequality each one
actually asserts. No ground-truth value is needed and no fitting step has seen
them.

Counts on the NP library are reported as a CURVE over five thresholds, never a
single number: picking one after seeing the counts is the failure mode this
discipline exists to prevent.

Curation steps 1-8 are replayed by importing ADR-AIIM-291's harness, not
rewritten. `10_new_system` is not modified; chembl251_raw.csv is read, not
re-pulled.
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
NPLIB = ROOT / "20_new_runs" / "np_library_funnel_2026-08-15" / "np_library_dock_ledger.csv"
SNPLIB = ROOT / "20_new_runs" / "snp_screen_2026-08-04" / "G1_standardization" / "snp_source_ledger.csv"
OUTDIR = ROOT / "20_new_runs" / "pki_regression_2026-10-05"
sys.path.insert(0, str(SCRIPTS))

PREREG = "b46f97a981c575287f5faffc414571e1b96a976895d40deb9a41b25934ca569c"
LOG_AT_FREEZE = "a36db9f9f3bcd8df25ce239e93e3aa163d02bd7d03978d1074837cde1e702e58"
RS, TEST_FRAC, SIG = 42, 0.20, 0.05
THRESHOLDS = [5.0, 5.5, 6.0, 6.5, 7.0]
CLASSIFIER_SPREAD = 38.88      # ADR-AIIM-293, same curated set
RELEASED = 72                  # the figure section 3.6 reports


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


def rebuild_table(cur):
    """Replay ADR-AIIM-291 steps 1-8 on the saved raw pull, keeping the rows
    step 9 discarded. Returns one row per (compound, endpoint_class)."""
    d = pd.read_csv(CUR_RUN / "chembl251_raw.csv", low_memory=False)
    n_raw = len(d)
    d = d[d["target_organism"] == "Homo sapiens"]
    d = d[d["assay_type"].isin(["B", "F"])]
    bad = d["data_validity_comment"].notna() & (
        d["data_validity_comment"] != "Manually validated")
    d = d[~bad]
    d = d[d["standard_type"].isin(cur.BINDING | cur.FUNCTIONAL)].copy()
    d["endpoint_class"] = np.where(d["standard_type"].isin(cur.BINDING),
                                   "binding", "functional")
    px_kind = d.apply(cur.to_plog, axis=1, result_type="expand")
    d["px"], d["kind"] = px_kind[0], px_kind[1]
    d = d[d["kind"] != "unusable"]
    uniq = d[["molecule_chembl_id", "canonical_smiles"]].drop_duplicates(
        "molecule_chembl_id")
    std = {r.molecule_chembl_id: cur.standardize(r.canonical_smiles)
           for r in uniq.itertuples()}
    d["std_smiles"] = d["molecule_chembl_id"].map(std)
    d = d[d["std_smiles"].notna()]
    recs = []
    for (mid, ec), g in d.groupby(["molecule_chembl_id", "endpoint_class"]):
        px, lab, basis = cur.aggregate_and_label(g)
        recs.append({"molecule_chembl_id": mid, "endpoint_class": ec,
                     "std_smiles": g["std_smiles"].iloc[0], "px": px,
                     "label": lab, "basis": basis})
    return pd.DataFrame(recs), n_raw


def main() -> int:
    (OUTDIR / "analysis").mkdir(parents=True, exist_ok=True)
    qsar = load("_qsar", SCRIPTS / "run_qsar_lbdd_2026-08-07.py")
    cur = load("_cur", CUR_RUN / "analysis" / "chembl251_curate.py")
    t0 = time.time()

    tbl, n_raw = rebuild_table(cur)
    b = tbl[tbl["endpoint_class"] == "binding"].reset_index(drop=True)
    pts = b[b["basis"] == "point_median"].reset_index(drop=True)
    up = b[b["basis"] == "censored_upper"].reset_index(drop=True)
    lo = b[b["basis"] == "censored_lower"].reset_index(drop=True)
    gap = b[b["basis"] == "point_median_gap"].reset_index(drop=True)
    print(f"raw {n_raw} -> binding {len(b)} compound-endpoint rows", flush=True)
    print(f"  point {len(pts)}  (of which gap-only, recovered from step 9: "
          f"{len(gap)})", flush=True)
    print(f"  censored upper {len(up)}  lower {len(lo)}  -- TEST ONLY", flush=True)

    # Point measurements, gap included: this is the regression set.
    reg = pd.concat([pts, gap], ignore_index=True)
    reg = reg[np.isfinite(reg["px"])].reset_index(drop=True)
    smiles = reg["std_smiles"].tolist()
    y = reg["px"].to_numpy(dtype=float)
    ikeys = [qsar._smiles_to_inchikey(s) or f"__N{i}__" for i, s in enumerate(smiles)]
    order = sorted(range(len(smiles)), key=lambda j: ikeys[j])
    smiles = [smiles[i] for i in order]
    y = y[order]
    scaf = [qsar._murcko(s) for s in smiles]
    X = qsar.compute_raw_features(smiles)
    print(f"  regression set n={len(y)}  pKi {y.min():.2f}-{y.max():.2f} "
          f"mean {y.mean():.2f}", flush=True)

    Xup = qsar.compute_raw_features(up["std_smiles"].tolist())
    yup = up["px"].to_numpy(dtype=float)
    Xlo = qsar.compute_raw_features(lo["std_smiles"].tolist())
    ylo = lo["px"].to_numpy(dtype=float)

    q_smiles, queries = {}, {}
    for lname, path, col in (("SNP-Lib", SNPLIB, "smiles_canonical"),
                             ("NP-library", NPLIB, "SMILES")):
        q = pd.read_csv(path)
        q = q[q[col].apply(lambda s: isinstance(s, str) and bool(s))]
        q_smiles[lname] = q[col].tolist()
        queries[lname] = qsar.compute_raw_features(q_smiles[lname])

    out = {"adr": "ADR-AIIM-297", "prereg_sha256": PREREG,
           "log_sha256_at_freeze": LOG_AT_FREEZE,
           "rdkit_version": rdkit.__version__,
           "composition": {
               "raw_rows": int(n_raw),
               "binding_compound_endpoint_rows": int(len(b)),
               "point": int(len(pts)), "gap_recovered": int(len(gap)),
               "regression_set": int(len(y)),
               "censored_upper_test_only": int(len(up)),
               "censored_lower_test_only": int(len(lo))},
           "thresholds": THRESHOLDS,
           "classifier_spread_same_set": CLASSIFIER_SPREAD,
           "conventions": {}}

    for cfg in CELLS:
        nm = cell_name(cfg)
        tr, te = split(scaf, len(y), **cfg)
        mdl = pipe().fit(X[tr], y[tr])
        p_te = mdl.predict(X[te])
        resid = np.abs(y[te] - p_te)
        # conformal regression: absolute-residual quantile at 1-SIG
        k = int(np.ceil((1 - SIG) * (len(resid) + 1))) - 1
        qhat = float(np.sort(resid)[min(k, len(resid) - 1)])

        rho = float(spearmanr(y[te], p_te).statistic)
        rec = {"n_train": int(len(tr)), "n_heldout": int(len(te)),
               "r2": round(float(r2_score(y[te], p_te)), 4),
               "rmse": round(float(np.sqrt(mean_squared_error(y[te], p_te))), 4),
               "mae": round(float(mean_absolute_error(y[te], p_te)), 4),
               "spearman": round(rho, 4),
               "conformal_qhat": round(qhat, 4)}

        # censored one-sided test: never trained on, no ground truth needed
        pu = mdl.predict(Xup)
        pl = mdl.predict(Xlo)
        rec["censored_test"] = {
            "upper_n": int(len(yup)),
            "upper_pass": int((pu <= yup).sum()),
            "upper_pass_frac": round(float((pu <= yup).mean()), 4),
            "upper_mean_violation": round(float(np.clip(pu - yup, 0, None).mean()), 4),
            "lower_n": int(len(ylo)),
            "lower_pass": int((pl >= ylo).sum()),
            "lower_pass_frac": round(float((pl >= ylo).mean()), 4)}

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
        c6 = rec["counts"]["NP-library"]["point"]["6.0"]
        print(f"  {nm:40s} rho {rho:.3f}  R2 {rec['r2']:.3f}  "
              f"RMSE {rec['rmse']:.3f}  cens-up {rec['censored_test']['upper_pass_frac']:.3f}  "
              f"NP@6.0 {c6}  [{time.time()-t0:.0f}s]", flush=True)

    # ------------------------------------------------------------- grading
    anchor = cell_name({"d1_shuffle": False, "d2_train_first": True,
                        "d3_divert": True})
    A = out["conventions"][anchor]
    np6 = [out["conventions"][cell_name(c)]["counts"]["NP-library"]["point"]["6.0"]
           for c in CELLS]
    spread = max(np6) / max(min(np6), 1)
    near72 = [n for n in np6 if 36 <= n <= 108]

    out["grading"] = {
        "Z1_spearman_above_0.60": {
            "anchor": anchor, "spearman": A["spearman"],
            "verdict": "MET" if A["spearman"] > 0.60 else "MISSED"},
        "Z2_upper_censored_pass_above_70pct": {
            "anchor_pass_frac": A["censored_test"]["upper_pass_frac"],
            "n": A["censored_test"]["upper_n"],
            "all_conventions": {cell_name(c): out["conventions"][cell_name(c)]
                                ["censored_test"]["upper_pass_frac"] for c in CELLS},
            "verdict": "MET" if A["censored_test"]["upper_pass_frac"] > 0.70
            else "MISSED"},
        "Z3_spread_under_10x": {
            "np_counts_at_6.0": {cell_name(c): n for c, n in zip(CELLS, np6)},
            "min": min(np6), "max": max(np6),
            "fold_range": round(spread, 2),
            "classifier_fold_range": CLASSIFIER_SPREAD,
            "verdict": "MET" if spread < 10.0 else "MISSED"},
        "Z4_none_within_50pct_of_72": {
            "band": [36, 108], "counts_in_band": near72,
            "verdict": "MET" if not near72 else "MISSED"},
    }
    out["seconds"] = round(time.time() - t0, 1)
    (OUTDIR / "pki_regression.json").write_text(json.dumps(out, indent=1),
                                                encoding="utf-8")

    print("\n" + "=" * 78)
    for t, g in out["grading"].items():
        print(f"{t:40s} {g['verdict']}")
    print("=" * 78)
    print(f"{'convention':40s} {'rho':>6s} {'R2':>7s} {'RMSE':>6s} "
          f"{'cens':>6s} {'NP@6':>6s} {'conf@6':>7s}")
    for c in CELLS:
        nm = cell_name(c)
        r = out["conventions"][nm]
        print(f"{nm:40s} {r['spearman']:6.3f} {r['r2']:7.3f} {r['rmse']:6.3f} "
              f"{r['censored_test']['upper_pass_frac']:6.3f} "
              f"{r['counts']['NP-library']['point']['6.0']:6d} "
              f"{r['counts']['NP-library']['conformal_lower']['6.0']:7d}")
    print(f"\nNP@6.0 spread {min(np6)}-{max(np6)} = {spread:.2f}x  "
          f"(classifier on the same set: {CLASSIFIER_SPREAD}x)")
    print(f"threshold curve, anchor convention, NP-library point:")
    for t in THRESHOLDS:
        print(f"  >= {t}: {A['counts']['NP-library']['point'][str(t)]:5d}   "
              f"conformal lower: {A['counts']['NP-library']['conformal_lower'][str(t)]:5d}")
    print(f"\nwritten -> {OUTDIR / 'pki_regression.json'} [{out['seconds']}s]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

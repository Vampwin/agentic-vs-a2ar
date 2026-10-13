"""ADR-AIIM-276, executed. Predictions frozen before this ran.

Pre-registration SHA-256
  00eeda4d6074caa619e207e53e98cb34319e36e48fe2373d521a880b1b18a118
Decision log at freeze
  3284c5a280f816fc4d36035db85eae00ce17c0355c527f410e56910ff3b38410

Three questions about the library the project set aside:

  R1  does SNP-Lib still return zero at the stricter significance levels,
      0.10 and 0.20, where ADR-AIIM-077's claim is also pitched?
  R2  do the released-split and fixed-split canonical variants flag the same
      SNP-Lib compounds, or only the same number of them?
  R3  among the 500 full-protocol compounds, does the count whose conformal
      set does not exclude activity exceed the 3 that ADR-AIIM-264 calls the
      binding constraint on the whole screening arm?

Canonical row order throughout, per ADR-AIIM-275: the labelled table is
sorted by InChIKey before anything reads it, so split, forest and CV see one
fixed table. Four configurations -- two arms by two split rules -- and three
significance levels. Compound identities are recorded, because R2 and R3 are
about which molecules, not how many.

Nothing is re-docked. No model promoted, no contract amended, no published
number changed. A compound flagged here is one the model declines to
exclude; it is not an active compound and nothing here has been assayed.
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
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import VarianceThreshold
from sklearn.pipeline import Pipeline

ROOT = Path(r"D:\Agentic VS Workflow")
SCRIPTS = ROOT / "10_new_system" / "scripts"
REPULL = ROOT / "20_new_runs" / "chembl_censored_repull_2026-10-01" / "chembl_a2a_all_relations_v3.csv"
SNPLIB = ROOT / "20_new_runs" / "snp_screen_2026-08-04" / "G1_standardization" / "snp_source_ledger.csv"
TOP500 = ROOT / "20_new_runs" / "np_library_funnel_2026-08-15" / "top500_for_full_protocol.csv"
OUTDIR = ROOT / "20_new_runs" / "snplib_canonical_2026-10-02"
sys.path.insert(0, str(SCRIPTS))

PREREG = "00eeda4d6074caa619e207e53e98cb34319e36e48fe2373d521a880b1b18a118"
RS, TEST_FRAC = 42, 0.20
SIGS = (0.05, 0.10, 0.20)


def load(name, fname):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / fname)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def split(scaffolds, n_rows, fixed):
    rng = random.Random(RS)
    m: dict[str, list[int]] = {}
    for i, sc in enumerate(scaffolds):
        m.setdefault(sc, []).append(i)
    ordered = (sorted(m.items(), key=lambda x: (-len(x[1]), x[0])) if fixed
               else sorted(m.items(), key=lambda x: len(x[1]), reverse=True))
    rng.shuffle(ordered)
    target = int(n_rows * TEST_FRAC)
    test, train = set(), set()
    for _sc, idxs in ordered:
        (test if len(test) < target else train).update(idxs)
    return np.array(sorted(train)), np.array(sorted(test))


def pipe():
    return Pipeline([("var", VarianceThreshold(threshold=0.01)),
                     ("rf", RandomForestClassifier(class_weight="balanced",
                                                   random_state=RS, n_jobs=3))])


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    qsar = load("_qsar", "run_qsar_lbdd_2026-08-07.py")
    cmp_mod = load("_cmp", "qsar_threshold_comparison.py")
    raw = pd.read_csv(cmp_mod.CHEMBL_RAW_SOURCE)

    d = cmp_mod.build_training_set(raw, 6.0, 5.5, qsar)
    smiles_A, y_A = d["canonical_smiles"].tolist(), list(d["label"].to_numpy())
    have = set(d["inchikey"])
    rep = pd.read_csv(REPULL).dropna(subset=["canonical_smiles"]).drop_duplicates()
    rep["ik"] = rep["canonical_smiles"].apply(qsar._smiles_to_inchikey)
    rep = rep.dropna(subset=["ik"])
    c = rep[rep["censoring_direction"] == "upper"].dropna(subset=["pchembl_upper_bound"])
    cm = c[~c["ik"].isin(have)].groupby("ik").agg(
        b=("pchembl_upper_bound", "min"), smi=("canonical_smiles", "first"))
    add = [r["smi"] for _, r in cm.iterrows() if r["b"] <= 5.5]
    arms = {"A": (smiles_A, np.array(y_A)),
            "B": (smiles_A + add, np.array(y_A + [0] * len(add)))}

    libs = {}
    sn = pd.read_csv(SNPLIB)
    sn = sn[sn["smiles_canonical"].apply(lambda s: isinstance(s, str) and bool(s))]
    libs["SNP-Lib"] = (sn["SNP_ID"].astype(str).tolist(),
                       qsar.compute_raw_features(sn["smiles_canonical"].tolist()))
    t5 = pd.read_csv(TOP500)
    t5 = t5[t5["SMILES"].apply(lambda s: isinstance(s, str) and bool(s))]
    libs["top500"] = (t5["compound_id"].astype(str).tolist(),
                      qsar.compute_raw_features(t5["SMILES"].tolist()))
    print(f"SNP-Lib {len(libs['SNP-Lib'][0])} | top500 {len(libs['top500'][0])}",
          flush=True)

    out = {"adr": "ADR-AIIM-276", "prereg_sha256": PREREG,
           "canonical_row_order": "labelled table sorted by InChIKey before use",
           "significance_levels": list(SIGS), "configurations": {}}
    t0 = time.time()

    for arm, (smiles, y) in arms.items():
        ikeys = [qsar._smiles_to_inchikey(s) or f"__NONE_{i}__"
                 for i, s in enumerate(smiles)]
        order = np.array(sorted(range(len(smiles)), key=lambda j: ikeys[j]))
        sm = [smiles[i] for i in order]
        yy = y[order]
        X = qsar.compute_raw_features(sm)
        scaf = [qsar._murcko(s) for s in sm]

        for fixed in (False, True):
            tag = f"{arm}_{'fixed' if fixed else 'released'}_split"
            tr, cal = split(scaf, len(yy), fixed)
            mdl = pipe().fit(X[tr], yy[tr])
            p_cal = mdl.predict_proba(X[cal])[:, 1]
            y_cal = yy[cal]
            cal1 = np.sort(1.0 - p_cal[y_cal == 1])
            n1 = len(cal1)
            cfg = {"calib_active": int(y_cal.sum()),
                   "calib_inactive": int((y_cal == 0).sum()), "libraries": {}}
            for lname, (ids, Xq) in libs.items():
                pa = mdl.predict_proba(Xq)[:, 1]
                p1 = (1.0 + (n1 - np.searchsorted(cal1, 1.0 - pa, side="left"))) / (1.0 + n1)
                rec = {"n": len(ids), "max_conformal_p": round(float(p1.max()), 4),
                       "by_significance": {}}
                for s in SIGS:
                    hits = [ids[j] for j in range(len(ids)) if p1[j] > s]
                    rec["by_significance"][str(s)] = {
                        "n_possibly_active": len(hits), "ids": sorted(hits)}
                cfg["libraries"][lname] = rec
                print(f"  {tag} [{lname}] " + "  ".join(
                    f"a={s}:{rec['by_significance'][str(s)]['n_possibly_active']}"
                    for s in SIGS) + f"  max_p={p1.max():.4f}", flush=True)
            out["configurations"][tag] = cfg

    # R2: do the two split rules flag the same SNP-Lib compounds?
    out["R2_identity"] = {}
    for arm in ("A", "B"):
        a = set(out["configurations"][f"{arm}_released_split"]["libraries"]["SNP-Lib"]["by_significance"]["0.05"]["ids"])
        b = set(out["configurations"][f"{arm}_fixed_split"]["libraries"]["SNP-Lib"]["by_significance"]["0.05"]["ids"])
        out["R2_identity"][arm] = {
            "released_only": sorted(a - b), "fixed_only": sorted(b - a),
            "shared": sorted(a & b), "identical": a == b}
    out["seconds"] = round(time.time() - t0, 1)
    (OUTDIR / "snplib_canonical.json").write_text(json.dumps(out, indent=1),
                                                  encoding="utf-8")
    print(f"\nwritten -> {OUTDIR / 'snplib_canonical.json'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

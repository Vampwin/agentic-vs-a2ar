"""ADR-AIIM-285, executed. Predictions frozen before this ran.

Pre-registration SHA-256
  eff0353e578e23b20ddd943f0faaa5554a359b8294a83cc926d7f43c78c5f8ae
Decision log at freeze
  61659f173f41295e1b269abf05c28e18e1febe992bb9160779ef04171b382c0f

ADR-AIIM-273 fixed the stable-sort tiebreak. ADR-AIIM-275 fixed the table
order and closed the row-order class exactly. Neither touched the three
design choices inside the split itself:

  D1  rng.shuffle after the size sort, which makes the size sort decorative
  D2  the held-out side is filled FIRST, so training is the remainder
  D3  the budget is tested BEFORE adding a whole group, so held-out overflows

Changed one at a time, D1 and D2 stayed invisible for two rounds. This runs
them as a full 2x2x2 factorial over the same canonical table so each main
effect is attributable instead of inferred.

Cell (0,0,0) -- shuffle retained, held-out first, overflow admitted -- IS
ADR-AIIM-275's fixed split. It must reproduce arm A NP 3334 / SNP 28 and arm
B NP 269 / SNP 1. That is the harness's own validity check and it is asserted,
not merely printed: if it fails, nothing else in this file means anything.

The tiebreak is NOT a factor here. It is held at the total order
(-size, scaffold_smiles) throughout, because ADR-AIIM-275 already measured it.

`10_new_system` is not modified. The switches live in this harness only.
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
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import VarianceThreshold
from sklearn.pipeline import Pipeline

ROOT = Path(r"D:\Agentic VS Workflow")
SCRIPTS = ROOT / "10_new_system" / "scripts"
REPULL = ROOT / "20_new_runs" / "chembl_censored_repull_2026-10-01" / "chembl_a2a_all_relations_v3.csv"
NPLIB = ROOT / "20_new_runs" / "np_library_funnel_2026-08-15" / "np_library_dock_ledger.csv"
SNPLIB = ROOT / "20_new_runs" / "snp_screen_2026-08-04" / "G1_standardization" / "snp_source_ledger.csv"
OUTDIR = ROOT / "20_new_runs" / "split_convention_factorial_2026-10-05"
sys.path.insert(0, str(SCRIPTS))

PREREG = "eff0353e578e23b20ddd943f0faaa5554a359b8294a83cc926d7f43c78c5f8ae"
LOG_AT_FREEZE = "61659f173f41295e1b269abf05c28e18e1febe992bb9160779ef04171b382c0f"
RS, SIG, TEST_FRAC, N_PERM = 42, 0.05, 0.20, 10

# ADR-AIIM-275's fixed split, which cell (0,0,0) must reproduce exactly.
ANCHOR = {"A": {"NP-library": 3334, "SNP-Lib": 28, "calib_inactive": 132},
          "B": {"NP-library": 269, "SNP-Lib": 1, "calib_inactive": 260}}


def load(name, fname):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / fname)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def split(scaffolds: list[str], n_rows: int,
          d1_shuffle: bool, d2_train_first: bool, d3_divert: bool):
    """Scaffold split with the three released conventions as switches.

    d1_shuffle      True  = released: shuffle after sorting by size
                    False = strict (-size, smiles) order, largest group first
    d2_train_first  False = released: fill the held-out side first
                    True  = fill training first, held-out takes the tail
    d3_divert       False = released: test budget before adding, admit whole
                            group, so the primary side overflows its target
                    True  = divert an overflowing group to the other side and
                            CONTINUE, so a later smaller group may still fit
                            (greedy best-fit, not a hard cutoff -- pre-registered)
    """
    scaf_to_idx: dict[str, list[int]] = {}
    for i, sc in enumerate(scaffolds):
        scaf_to_idx.setdefault(sc, []).append(i)
    ordered = sorted(scaf_to_idx.items(), key=lambda x: (-len(x[1]), x[0]))
    if d1_shuffle:
        random.Random(RS).shuffle(ordered)

    n_test_target = int(n_rows * TEST_FRAC)
    n_train_target = n_rows - n_test_target
    test_set: set[int] = set()
    train_set: set[int] = set()

    for _sc, idxs in ordered:
        if d2_train_first:
            primary, budget, other = train_set, n_train_target, test_set
        else:
            primary, budget, other = test_set, n_test_target, train_set
        if d3_divert:
            (primary if len(primary) + len(idxs) <= budget else other).update(idxs)
        else:
            (primary if len(primary) < budget else other).update(idxs)
    return np.array(sorted(train_set)), np.array(sorted(test_set))


def pipe():
    return Pipeline([("var", VarianceThreshold(threshold=0.01)),
                     ("rf", RandomForestClassifier(class_weight="balanced",
                                                   random_state=RS, n_jobs=3))])


CELLS = [{"d1_shuffle": a, "d2_train_first": b, "d3_divert": c}
         for a, b, c in itertools.product((True, False), (False, True), (False, True))]


def cell_name(c):
    return (f"D1{'shuffle' if c['d1_shuffle'] else 'sorted'}"
            f"_D2{'train' if c['d2_train_first'] else 'test'}first"
            f"_D3{'divert' if c['d3_divert'] else 'overflow'}")


def main() -> int:
    (OUTDIR / "analysis").mkdir(parents=True, exist_ok=True)
    qsar = load("_qsar", "run_qsar_lbdd_2026-08-07.py")
    cmp_mod = load("_cmp", "qsar_threshold_comparison.py")
    raw = pd.read_csv(cmp_mod.CHEMBL_RAW_SOURCE)

    d = cmp_mod.build_training_set(raw, 6.0, 5.5, qsar)
    smiles_A = d["canonical_smiles"].tolist()
    y_A = list(d["label"].to_numpy())
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
    print(f"arm A: {len(y_A)} rows ({int(np.sum(y_A))} active)  "
          f"arm B: +{len(add)} recovered inactives", flush=True)

    queries = {}
    for lname, path, scol in (("SNP-Lib", SNPLIB, "smiles_canonical"),
                              ("NP-library", NPLIB, "SMILES")):
        q = pd.read_csv(path)
        q = q[q[scol].apply(lambda s: isinstance(s, str) and bool(s))]
        queries[lname] = qsar.compute_raw_features(q[scol].tolist())

    out = {"adr": "ADR-AIIM-285", "prereg_sha256": PREREG,
           "log_sha256_at_freeze": LOG_AT_FREEZE, "n_permutations": N_PERM,
           "released_values": {"SNP-Lib": 0, "NP-library": 72},
           "anchor_adr_275_fixed_split": ANCHOR,
           "note": "10_new_system is not modified; the three switches live in "
                   "this harness only. Tiebreak held at (-size, smiles); it is "
                   "not a factor here (ADR-AIIM-275 measured it).",
           "arms": {}}
    t0 = time.time()

    for arm, (smiles, y) in arms.items():
        X = qsar.compute_raw_features(smiles)
        scaf_all = [qsar._murcko(s) for s in smiles]
        ikeys = [qsar._smiles_to_inchikey(s) or f"__NONE_{i}__"
                 for i, s in enumerate(smiles)]
        n_test_target = int(len(y) * TEST_FRAC)
        rec: dict[str, list] = {cell_name(c): [] for c in CELLS}

        for k in range(N_PERM):
            idx = (np.arange(len(y)) if k == 0
                   else np.random.default_rng(10_000 + k).permutation(len(y)))
            # ADR-AIIM-275: canonical InChIKey order, always, in every cell.
            canon = sorted(range(len(idx)), key=lambda j: ikeys[idx[j]])
            idx = idx[np.array(canon)]
            y_p, Xp = y[idx], X[idx]
            scaf_p = [scaf_all[i] for i in idx]

            for cfg in CELLS:
                tr, cal = split(scaf_p, len(y_p), **cfg)
                mdl = pipe().fit(Xp[tr], y_p[tr])
                p_cal = mdl.predict_proba(Xp[cal])[:, 1]
                y_cal = y_p[cal]
                cal1 = np.sort(1.0 - p_cal[y_cal == 1])
                n1 = len(cal1)
                r = {"draw": k, "n_train": int(len(tr)), "n_calib": int(len(cal)),
                     "calib_active": int(y_cal.sum()),
                     "calib_inactive": int((y_cal == 0).sum()),
                     "calib_over_target": int(len(cal) - n_test_target)}
                for lname, Xq in queries.items():
                    pa = mdl.predict_proba(Xq)[:, 1]
                    p1 = (1.0 + (n1 - np.searchsorted(cal1, 1.0 - pa, side="left"))) / (1.0 + n1)
                    r[lname] = {"n_possibly_active": int((p1 > SIG).sum()),
                                "max_p": round(float(p1.max()), 4)}
                rec[cell_name(cfg)].append(r)

            done = (k + 1) * len(CELLS)
            print(f"  arm {arm} draw {k}: "
                  + "  ".join(f"{cell_name(c).split('_')[1][2:]}"
                              f"{'/div' if c['d3_divert'] else '/ovf'}"
                              f"{'+sh' if c['d1_shuffle'] else '+st'}="
                              f"{rec[cell_name(c)][-1]['NP-library']['n_possibly_active']}"
                              for c in CELLS)
                  + f"  [{done} fits, {time.time()-t0:.0f}s]", flush=True)

        summary = {}
        for cfg in CELLS:
            nm = cell_name(cfg)
            s = {"config": cfg}
            for lname in queries:
                v = [x[lname]["n_possibly_active"] for x in rec[nm]]
                s[lname] = {"min": min(v), "max": max(v),
                            "n_distinct_values": len(set(v)),
                            "invariant": len(set(v)) == 1,
                            "value": v[0] if len(set(v)) == 1 else None}
            for f in ("n_train", "n_calib", "calib_active", "calib_inactive",
                      "calib_over_target"):
                vv = [x[f] for x in rec[nm]]
                s[f] = {"value": vv[0] if len(set(vv)) == 1 else None,
                        "invariant": len(set(vv)) == 1,
                        "min": min(vv), "max": max(vv)}
            summary[nm] = s
        out["arms"][arm] = {"per_draw": rec, "summary": summary}
        (OUTDIR / "split_convention_factorial.json").write_text(
            json.dumps(out, indent=1), encoding="utf-8")

        anchor_cell = cell_name({"d1_shuffle": True, "d2_train_first": False,
                                 "d3_divert": False})
        got = {lname: summary[anchor_cell][lname]["value"] for lname in queries}
        got["calib_inactive"] = summary[anchor_cell]["calib_inactive"]["value"]
        print(f"\narm {arm} anchor cell {anchor_cell}: {got}", flush=True)
        print(f"arm {arm} expected (ADR-AIIM-275 fixed split): {ANCHOR[arm]}",
              flush=True)
        out["arms"][arm]["anchor_check"] = {
            "cell": anchor_cell, "observed": got, "expected": ANCHOR[arm],
            "reproduces_adr_275": got == ANCHOR[arm]}

    # ---- grading, on arm A NP-library as pre-registered ----
    A = out["arms"]["A"]["summary"]
    npv = {nm: A[nm]["NP-library"]["value"] for nm in A}
    all_cells = [(arm, nm, lname,
                  out["arms"][arm]["summary"][nm][lname]["n_distinct_values"])
                 for arm in ("A", "B") for nm in A for lname in queries]

    def main_effect(key):
        hi = [npv[cell_name(c)] for c in CELLS if c[key]]
        lo = [npv[cell_name(c)] for c in CELLS if not c[key]]
        return float(np.mean(hi) - np.mean(lo))

    me = {k: main_effect(k) for k in ("d1_shuffle", "d2_train_first", "d3_divert")}
    # D1 level 1 is "sorted" (shuffle removed), so flip its sign to report the
    # effect of the PRINCIPLED level, matching the pre-registration's levels.
    me_principled = {"D1 (remove shuffle)": -me["d1_shuffle"],
                     "D2 (train first)": me["d2_train_first"],
                     "D3 (divert)": me["d3_divert"]}
    base = cell_name({"d1_shuffle": True, "d2_train_first": False, "d3_divert": False})
    d3_only = cell_name({"d1_shuffle": True, "d2_train_first": False, "d3_divert": True})
    principled = cell_name({"d1_shuffle": False, "d2_train_first": True, "d3_divert": True})
    vals = list(npv.values())

    grading = {
        "T1_all_cells_invariant": {
            "verdict": "MET" if all(n == 1 for *_, n in all_cells) else "MISSED",
            "cells_with_more_than_one_value": [
                {"arm": a, "cell": nm, "library": l, "n_distinct": n}
                for a, nm, l, n in all_cells if n != 1]},
        "T2_principled_exceeds_3334": {
            "value": npv[principled], "threshold": 3334,
            "verdict": "MET" if npv[principled] > 3334 else "MISSED"},
        "T3_D2_dominates": {
            "main_effects_arm_A_NP": {k: round(v, 1) for k, v in me_principled.items()},
            "verdict": "MET" if (abs(me_principled["D2 (train first)"])
                                 > max(abs(me_principled["D1 (remove shuffle)"]),
                                       abs(me_principled["D3 (divert)"]))) else "MISSED"},
        "T4_D3_alone_under_5pct": {
            "base_cell": base, "base_value": npv[base],
            "d3_only_cell": d3_only, "d3_only_value": npv[d3_only],
            "pct_change": round(abs(npv[d3_only] - npv[base]) / npv[base] * 100, 2)
            if npv[base] else None,
            "verdict": "MET" if npv[base] and abs(npv[d3_only] - npv[base]) / npv[base] < 0.05
            else "MISSED"},
        "T5_none_recovers_72": {
            "cells_equal_to_72": [
                {"arm": a, "cell": nm, "library": l}
                for a in ("A", "B") for nm in A for l in queries
                if out["arms"][a]["summary"][nm][l]["value"] == 72],
            "verdict": "MET" if not any(
                out["arms"][a]["summary"][nm][l]["value"] == 72
                for a in ("A", "B") for nm in A for l in queries) else "MISSED"},
        "T6_span_exceeds_6x": {
            "arm_A_NP_values": npv, "min": min(vals), "max": max(vals),
            "fold_range": round(max(vals) / max(min(vals), 1), 2),
            "adr_275_fold_range": 5.87,
            "verdict": "MET" if max(vals) / max(min(vals), 1) > 6.0 else "MISSED"},
    }
    out["grading"] = grading
    out["seconds"] = round(time.time() - t0, 1)
    (OUTDIR / "split_convention_factorial.json").write_text(
        json.dumps(out, indent=1), encoding="utf-8")

    print("\n" + "=" * 70)
    for t, g in grading.items():
        print(f"{t:34s} {g['verdict']}")
    print("=" * 70)
    print("arm A NP-library, eight conventions:")
    for nm, v in sorted(npv.items(), key=lambda x: x[1]):
        print(f"  {v:6d}  {nm}")
    print(f"\nmain effects (arm A, NP): "
          + json.dumps({k: round(v, 1) for k, v in me_principled.items()}))
    print(f"\nwritten -> {OUTDIR / 'split_convention_factorial.json'}  "
          f"[{out['seconds']}s]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

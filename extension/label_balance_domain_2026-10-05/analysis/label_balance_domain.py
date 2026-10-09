"""ADR-AIIM-287, executed. Predictions frozen before this ran.

Pre-registration SHA-256
  3a0d096a1ba9ff3aa031d00c81285d334d286e846723dcb7db8dcdd311af969a
Decision log at freeze
  1ed078b4a8b48be37e507cc01737e70d50351935da81b82e56ca7099304e16bb

Three claims were asserted after ADR-AIIM-286 and none was tested:
  U1  Ki-only labels would not collapse the convention spread
  U2  calibration balance is not the driver (argued from an OBSERVATIONAL
      correlation across eight cells where balance was never manipulated)
  U3  restricting to the applicability domain cannot help

Here balance is MANIPULATED, the label source is switched, and the query set
is restricted by the project's own AD rule. 8 conventions x 3 label arms x
2 calibration modes, each reporting full-library and in-domain counts.

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
OUTDIR = ROOT / "20_new_runs" / "label_balance_domain_2026-10-05"
sys.path.insert(0, str(SCRIPTS))

PREREG = "3a0d096a1ba9ff3aa031d00c81285d334d286e846723dcb7db8dcdd311af969a"
LOG_AT_FREEZE = "1ed078b4a8b48be37e507cc01737e70d50351935da81b82e56ca7099304e16bb"
RS, SIG, TEST_FRAC, N_PERM = 42, 0.05, 0.20, 2
AD_SAMPLE = 500


def load(name, fname):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / fname)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def split(scaffolds, n_rows, d1_shuffle, d2_train_first, d3_divert):
    """Identical to ADR-AIIM-286's split. Tiebreak held at (-size, smiles)."""
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


def build_arms(qsar, cmp_mod):
    """A = released. B = A + censored recovery. K = Ki-only, same recovery."""
    raw = pd.read_csv(cmp_mod.CHEMBL_RAW_SOURCE)
    rep = pd.read_csv(REPULL).dropna(subset=["canonical_smiles"]).drop_duplicates()
    rep["ik"] = rep["canonical_smiles"].apply(qsar._smiles_to_inchikey)
    rep = rep.dropna(subset=["ik"])

    def recovery(have, ki_only):
        c = rep[rep["censoring_direction"] == "upper"].dropna(subset=["pchembl_upper_bound"])
        if ki_only:
            c = c[c["standard_type"] == "Ki"]
        cm = c[~c["ik"].isin(have)].groupby("ik").agg(
            b=("pchembl_upper_bound", "min"), smi=("canonical_smiles", "first"))
        return [r["smi"] for _, r in cm.iterrows() if r["b"] <= 5.5]

    out = {}
    for arm, src, ki in (("A", raw, False), ("B", raw, False),
                         ("K", raw[raw["standard_type"] == "Ki"], True)):
        d = cmp_mod.build_training_set(src, 6.0, 5.5, qsar)
        smi = d["canonical_smiles"].tolist()
        y = list(d["label"].to_numpy())
        if arm == "A":
            out[arm] = (smi, np.array(y))
            continue
        add = recovery(set(d["inchikey"]), ki)
        out[arm] = (smi + add, np.array(y + [0] * len(add)))
        print(f"  arm {arm}: {len(smi)} base + {len(add)} recovered inactives",
              flush=True)
    return out


def main() -> int:
    (OUTDIR / "analysis").mkdir(parents=True, exist_ok=True)
    qsar = load("_qsar", "run_qsar_lbdd_2026-08-07.py")
    cmp_mod = load("_cmp", "qsar_threshold_comparison.py")
    arms = build_arms(qsar, cmp_mod)

    q_smiles, queries = {}, {}
    for lname, path, scol in (("SNP-Lib", SNPLIB, "smiles_canonical"),
                              ("NP-library", NPLIB, "SMILES")):
        q = pd.read_csv(path)
        q = q[q[scol].apply(lambda s: isinstance(s, str) and bool(s))]
        q_smiles[lname] = q[scol].tolist()
        queries[lname] = qsar.compute_raw_features(q_smiles[lname])

    out = {"adr": "ADR-AIIM-287", "prereg_sha256": PREREG,
           "log_sha256_at_freeze": LOG_AT_FREEZE, "n_permutations": N_PERM,
           "released_values": {"SNP-Lib": 0, "NP-library": 72},
           "ad_note": "AD computed once per arm against a seeded 500-compound "
                      "sample of that arm's labelled set, not per cell. "
                      "Declared simplification (ADR-AIIM-287).",
           "arms": {}}
    t0 = time.time()

    for arm, (smiles, y) in arms.items():
        X = qsar.compute_raw_features(smiles)
        scaf_all = [qsar._murcko(s) for s in smiles]
        ikeys = [qsar._smiles_to_inchikey(s) or f"__NONE_{i}__"
                 for i, s in enumerate(smiles)]

        # ---- applicability domain, once per arm (declared simplification) ----
        ref = random.Random(RS).sample(smiles, min(AD_SAMPLE, len(smiles)))
        ad_thr = qsar.compute_ad_threshold(ref)
        in_dom = {}
        for lname in queries:
            sc = qsar.compute_ad_scores(q_smiles[lname], ref)
            in_dom[lname] = sc >= ad_thr
            print(f"  arm {arm} AD: {lname} {int(in_dom[lname].sum())}/"
                  f"{len(sc)} in-domain (thr {ad_thr:.3f})", flush=True)

        rec: dict[str, list] = {}
        for k in range(N_PERM):
            idx = (np.arange(len(y)) if k == 0
                   else np.random.default_rng(10_000 + k).permutation(len(y)))
            canon = sorted(range(len(idx)), key=lambda j: ikeys[idx[j]])
            idx = idx[np.array(canon)]
            y_p, Xp = y[idx], X[idx]
            scaf_p = [scaf_all[i] for i in idx]

            for cfg in CELLS:
                tr, cal = split(scaf_p, len(y_p), **cfg)
                mdl = pipe().fit(Xp[tr], y_p[tr])
                p_cal = mdl.predict_proba(Xp[cal])[:, 1]
                y_cal = y_p[cal]
                act_pos = np.where(y_cal == 1)[0]
                n_inact = int((y_cal == 0).sum())

                for mode in ("raw", "balanced"):
                    if mode == "balanced":
                        # U2: MANIPULATE the balance the p-value actually reads.
                        # cal1 is built from actives alone, so subsample actives
                        # down to the inactive count. Seeded, order-independent.
                        keep = np.sort(np.random.default_rng(RS).choice(
                            act_pos, size=min(n_inact, len(act_pos)),
                            replace=False))
                    else:
                        keep = act_pos
                    cal1 = np.sort(1.0 - p_cal[keep])
                    n1 = len(cal1)
                    key = f"{cell_name(cfg)}|{mode}"
                    r = {"draw": k, "n_calib": int(len(cal)),
                         "calib_active_used": n1,
                         "calib_inactive": n_inact}
                    for lname, Xq in queries.items():
                        pa = mdl.predict_proba(Xq)[:, 1]
                        p1 = (1.0 + (n1 - np.searchsorted(
                            cal1, 1.0 - pa, side="left"))) / (1.0 + n1)
                        flag = p1 > SIG
                        r[lname] = {
                            "n_possibly_active": int(flag.sum()),
                            "n_in_domain": int(in_dom[lname].sum()),
                            "n_possibly_active_in_domain":
                                int((flag & in_dom[lname]).sum())}
                    rec.setdefault(key, []).append(r)
            print(f"  arm {arm} draw {k} done [{time.time()-t0:.0f}s]", flush=True)

        summary = {}
        for key, rows in rec.items():
            s = {}
            for lname in queries:
                for f in ("n_possibly_active", "n_possibly_active_in_domain"):
                    v = [x[lname][f] for x in rows]
                    s[f"{lname}.{f}"] = {"value": v[0] if len(set(v)) == 1 else None,
                                         "min": min(v), "max": max(v),
                                         "invariant": len(set(v)) == 1}
                s[f"{lname}.n_in_domain"] = rows[0][lname]["n_in_domain"]
            for f in ("n_calib", "calib_active_used", "calib_inactive"):
                s[f] = rows[0][f]
            summary[key] = s
        out["arms"][arm] = {"per_draw": rec, "summary": summary}
        (OUTDIR / "label_balance_domain.json").write_text(
            json.dumps(out, indent=1), encoding="utf-8")

    # ---------------- grading ----------------
    def spans(arm, mode, field):
        S = out["arms"][arm]["summary"]
        v = [S[f"{cell_name(c)}|{mode}"][field]["value"] for c in CELLS]
        v = [x for x in v if x is not None]
        return (min(v), max(v), round(max(v) / max(min(v), 1), 2)) if v else None

    full = "NP-library.n_possibly_active"
    dom = "NP-library.n_possibly_active_in_domain"
    sp_raw = {a: spans(a, "raw", full) for a in arms}
    sp_bal = {a: spans(a, "balanced", full) for a in arms}
    sp_dom = {a: spans(a, "raw", dom) for a in arms}

    princ = cell_name({"d1_shuffle": False, "d2_train_first": True, "d3_divert": True})
    pB = out["arms"]["B"]["summary"][f"{princ}|raw"][full]["value"]
    pK = out["arms"]["K"]["summary"][f"{princ}|raw"][full]["value"]

    n_bal_collapsed = sum(1 for a in arms if sp_bal[a] and sp_bal[a][2] <= 6.0)
    n_dom_collapsed = sum(1 for a in arms if sp_dom[a] and sp_dom[a][2] <= 6.0)

    out["grading"] = {
        "U1_Ki_only_does_not_collapse": {
            "arm_K_span": sp_raw["K"],
            "verdict": "MET" if sp_raw["K"] and sp_raw["K"][2] > 6.0 else "MISSED"},
        "U2_balance_does_not_collapse": {
            "spans_balanced": sp_bal, "spans_raw": sp_raw,
            "arms_collapsed": n_bal_collapsed,
            "verdict": "MET" if n_bal_collapsed < 2 else "MISSED"},
        "U3_in_domain_does_not_collapse": {
            "spans_in_domain": sp_dom,
            "arms_collapsed": n_dom_collapsed,
            "verdict": "MET" if n_dom_collapsed < 2 else "MISSED"},
        "U4_label_source_moves_answer": {
            "arm_B_principled": pB, "arm_K_principled": pK,
            "pct_diff": round(abs(pK - pB) / max(pB, 1) * 100, 1),
            "verdict": "MET" if pB and abs(pK - pB) / pB > 0.20 else "MISSED"},
    }
    out["seconds"] = round(time.time() - t0, 1)
    (OUTDIR / "label_balance_domain.json").write_text(
        json.dumps(out, indent=1), encoding="utf-8")

    print("\n" + "=" * 72)
    for t, g in out["grading"].items():
        print(f"{t:36s} {g['verdict']}")
    print("=" * 72)
    print(f"{'arm':4s} {'raw span (min,max,fold)':30s} {'balanced':30s} {'in-domain':30s}")
    for a in arms:
        print(f"{a:4s} {str(sp_raw[a]):30s} {str(sp_bal[a]):30s} {str(sp_dom[a]):30s}")
    print(f"\nin-domain sizes: " + json.dumps(
        {a: out['arms'][a]['summary'][f'{princ}|raw']['NP-library.n_in_domain']
         for a in arms}))
    print(f"principled cell, raw: arm B {pB}  arm K {pK}")
    print(f"\nwritten -> {OUTDIR / 'label_balance_domain.json'} [{out['seconds']}s]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""A1 and A2 of ADR-AIIM-399: the focus-library selection with one vote per distinct split, and its sensitivity to the
applicability threshold. Nothing is refitted: the predictions of the 69 pool compounds under the eight conventions are read
from snplib_dock_2026-10-06/pool69_predictions.csv (one-thread fits), and a prediction does not depend on the pool.

The applicability score of the expansion library against the 6,284-compound reference is recomputed with the registered code
path (pool_comparison.py: rebuild_table, InChIKey order, ECFP4 radius 2 / 2,048 bits, mean of the five highest Tanimoto).
"""
import importlib.util
import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger

RDLogger.DisableLog("rdApp.*")
ROOT = Path(r"D:\Agentic VS Workflow")
RUNS = ROOT / "20_new_runs"
OUT = RUNS / "v9_sensitivity_2026-10-11"
SCRIPTS = ROOT / "10_new_system" / "scripts"
sys.path.insert(0, str(SCRIPTS))
TOP, MAJ8, MAJ5, CAP = 20, 5, 3, 10
N_MATCH = 1560            # expansion compounds at or above 0.25 against the 4,566-compound training side (31.2%)


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


t0 = time.time()
P = pd.read_csv(RUNS / "snplib_dock_2026-10-06" / "pool69_predictions.csv")
cells = [c for c in P.columns if c.startswith("pred_")]
ES = json.loads((RUNS / "em_survivors_2026-10-06" / "em_survivors.json").read_text("utf-8"))
groups = [["pred_" + n for n in g] for g in ES["distinct_split_groups"]]
assert sorted(c for g in groups for c in g) == sorted(cells) and len(groups) == 5
SURV = pd.read_csv(RUNS / "em_survivors_2026-10-06" / "survivors.csv")
FIFTEEN, EIGHT = SURV["SNP_ID"].tolist(), ES["shortlist"]
P["mean8"] = P[cells].mean(axis=1)


def select(pool, vote):
    """pool: rows of P. vote 'conv': at least 5 of 8 conventions; vote 'split': at least 3 of 5 distinct splits."""
    pool = pool.reset_index(drop=True)
    cols = cells if vote == "conv" else [g[0] for g in groups]
    need = MAJ8 if vote == "conv" else MAJ5
    cnt = Counter()
    for c in cols:
        top = np.argsort(-pool[c].to_numpy(), kind="stable")[:TOP]
        cnt.update(pool.loc[top, "SNP_ID"].tolist())
    mean = dict(zip(pool["SNP_ID"], pool["mean8"]))
    fw = dict(zip(pool["SNP_ID"], pool["murcko"]))
    surv = {k: v for k, v in cnt.items() if v >= need}
    ranked = sorted(surv, key=lambda c: (-surv[c], -mean[c]))
    short, seen = [], set()
    for c in ranked:
        if fw[c] in seen:
            continue
        seen.add(fw[c]); short.append(c)
        if len(short) == CAP:
            break
    return ranked, short, surv


def jac(a, b):
    a, b = set(a), set(b)
    return round(len(a & b) / len(a | b), 3) if a | b else None


res = {"adr": "ADR-AIIM-399", "inputs": {"pool69_predictions.csv": "snplib_dock_2026-10-06", "distinct_split_groups": ES["distinct_split_groups"]}}

# ------------------------------------------------------------------ A1
r8, s8, _ = select(P, "conv")
res["A1"] = {"anchor": {"fifteen_reproduced": set(r8) == set(FIFTEEN) and len(r8) == 15,
                        "eight_reproduced_in_order": s8 == EIGHT}}
r5, s5, c5 = select(P, "split")
res["A1"].update({
    "vote_3_of_5": {"selected": r5, "n_selected": len(r5), "shortlist": s5, "votes": {k: int(c5[k]) for k in r5}},
    "added_to_fifteen": [c for c in r5 if c not in FIFTEEN], "removed_from_fifteen": [c for c in FIFTEEN if c not in r5],
    "jaccard_selected_vs_fifteen": jac(r5, FIFTEEN), "jaccard_shortlist_vs_eight": jac(s5, EIGHT),
    "eight_kept_in_shortlist": [c for c in EIGHT if c in s5], "shortlist_new": [c for c in s5 if c not in EIGHT],
    "P1_anchor": set(r8) == set(FIFTEEN) and s8 == EIGHT,
    "P2_all_fifteen_selected_and_at_most_three_added": set(FIFTEEN) <= set(r5) and len(set(r5) - set(FIFTEEN)) <= 3,
    "P3_at_least_six_of_eight_in_shortlist": len(set(EIGHT) & set(s5)) >= 6})
print("A1", json.dumps(res["A1"], indent=1), flush=True)

# ------------------------------------------------------------------ A2
qsar = load("_qsar", SCRIPTS / "run_qsar_lbdd_2026-08-07.py")
cur = load("_cur", RUNS / "chembl251_full_curation_2026-10-05" / "analysis" / "chembl251_curate.py")
rg = load("_reg", RUNS / "pki_regression_2026-10-05" / "analysis" / "pki_regression.py")
pc = load("_pc", RUNS / "pool_comparison_2026-10-06" / "analysis" / "pool_comparison.py")
tbl, _ = rg.rebuild_table(cur)
b = tbl[tbl["endpoint_class"] == "binding"]
b = b[b["basis"].isin(["point_median", "point_median_gap", "censored_upper", "censored_lower"])]
b = b[np.isfinite(b["px"])].reset_index(drop=True)
msm = b["std_smiles"].tolist()
ref = [pc._morgan.GetFingerprint(Chem.MolFromSmiles(s)) for s in msm if Chem.MolFromSmiles(s)]
print(f"reference {len(ref)}  [{time.time()-t0:.0f}s]", flush=True)


def ad(smiles):
    out = []
    for s in smiles:
        m = Chem.MolFromSmiles(s) if isinstance(s, str) else None
        out.append(float(np.mean(sorted(DataStructs.BulkTanimotoSimilarity(pc._morgan.GetFingerprint(m), ref), reverse=True)[:5]))
                   if m else np.nan)
    return np.array(out)


npl = pd.read_csv(pc.NPLIB)
npl = npl[npl["SMILES"].notna()].reset_index(drop=True)
ad_np = ad(npl["SMILES"])
n_025 = int(np.nansum(ad_np >= 0.25))
srt = np.sort(ad_np[np.isfinite(ad_np)])[::-1]
t_star = float(srt[N_MATCH - 1])
# the SNP-Lib pool: scores recomputed for the stored 69 from their canonical SMILES in the gate
sl_raw = pd.read_csv(pc.SNPLIB)
sl_raw = sl_raw[sl_raw["smiles_canonical"].notna()].reset_index(drop=True)
pc.AD_FLOOR = -1.0                                         # this run's copy of the module only: keep every gated compound, with its score
A_all, _wf = pc.gate(sl_raw, "SNP_ID", "smiles_canonical", ref, qsar, "SNP-Lib")
A_all = A_all.set_index("uid")
pool69 = A_all[A_all["ad_score"] >= 0.25]
stored = P.set_index("SNP_ID")["ad_score"]
anchor_pool = (len(pool69) == 69 and set(pool69.index) == set(stored.index)
               and float(np.max(np.abs(pool69.loc[stored.index, "ad_score"].to_numpy() - stored.to_numpy()))) < 1e-6)
res["A2"] = {"reference_n": len(ref), "expansion_n": int(np.isfinite(ad_np).sum()),
             "anchor": {"expansion_at_or_above_0.25": n_025, "expected": 2645, "holds": n_025 == 2645,
                        "pool69_scores_reproduced": bool(anchor_pool)},
             "n_gated_before_AD": int(len(A_all)),
             "t_star": round(t_star, 4), "expansion_at_or_above_t_star": int(np.nansum(ad_np >= t_star)),
             "thresholds": {}}
for name, t in (("0.25", 0.25), ("0.30", 0.30), ("0.35", 0.35), ("t_star", t_star)):
    pool = P[P["ad_score"] >= t]
    rc, sc, _ = select(pool, "conv") if len(pool) else ([], [], {})
    rs, ss, _ = select(pool, "split") if len(pool) else ([], [], {})
    res["A2"]["thresholds"][name] = {
        "threshold": round(float(t), 4), "pool": int(len(pool)),
        "expansion_fraction_admitted": round(float(np.nanmean(ad_np >= t)), 4),
        "fifteen_in_pool": [c for c in FIFTEEN if c in set(pool["SNP_ID"])],
        "eight_in_pool": [c for c in EIGHT if c in set(pool["SNP_ID"])],
        "vote_5_of_8": {"n_selected": len(rc), "selected": rc, "shortlist": sc,
                        "fifteen_still_selected": [c for c in FIFTEEN if c in rc], "new": [c for c in rc if c not in FIFTEEN]},
        "vote_3_of_5": {"n_selected": len(rs), "selected": rs, "shortlist": ss,
                        "fifteen_still_selected": [c for c in FIFTEEN if c in rs], "new": [c for c in rs if c not in FIFTEEN]}}
T = res["A2"]["thresholds"]["t_star"]
res["A2"].update({"P4_anchors": n_025 == 2645 and bool(anchor_pool),
                  "P5_t_star_between_0.28_and_0.33": 0.28 <= t_star <= 0.33,
                  "P6_pool_at_t_star_between_15_and_40": 15 <= T["pool"] <= 40,
                  "P7_at_most_7_of_fifteen_in_pool_at_t_star": len(T["fifteen_in_pool"]) <= 7})
res["seconds"] = round(time.time() - t0, 1)
(OUT / "a1_a2_vote_and_threshold.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
npl.assign(ad_score_vs_6284=ad_np)[["compound_id", "ad_score_vs_6284"]].to_csv(OUT / "a2_expansion_ad_scores.csv", index=False)
A_all.reset_index()[["uid", "ad_score"]].to_csv(OUT / "a2_snplib_gated_ad_scores.csv", index=False)
print("A2", json.dumps({k: v for k, v in res["A2"].items() if k != "thresholds"}, indent=1))
for k, v in res["A2"]["thresholds"].items():
    print(k, v["threshold"], "pool", v["pool"], "| 5of8:", v["vote_5_of_8"]["n_selected"], "kept", len(v["vote_5_of_8"]["fifteen_still_selected"]),
          "| 3of5:", v["vote_3_of_5"]["n_selected"], "kept", len(v["vote_3_of_5"]["fifteen_still_selected"]),
          "| fifteen in pool", len(v["fifteen_in_pool"]), "eight in pool", len(v["eight_in_pool"]))

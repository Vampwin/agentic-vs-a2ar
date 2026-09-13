"""analyse_r2_smiles_redock.py -- R2, evaluated against ADR-AIIM-133.

The pre-registration sets no correlation threshold, deliberately: the median
adjacent score gap across the 497 redocked compounds is 0.0020 kcal/mol, so that
ordering is not resolved by the score and a Spearman cut-off would test almost
nothing. The endpoint is instead binary and needs no invented number -- which
released claims change.

Claims under test, transcribed from the pre-registration:
  * the identity of the three Tier-A candidates
  * 497 scored, 482 on-pocket, 3 failed to converge
  * the tier counts 3 / 119 / 360 / 15 / 3
  * that the strongest raw docking score belongs to the second tier

Everything else -- deviations, Spearman, rank movement -- is reported
descriptively with no verdict attached.

Usage:
    python analyse_r2_smiles_redock.py [--run DIR]
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
from pathlib import Path
from kdss_paths import project_root, resolve, under, dev_dir

import numpy as np
import pandas as pd
from scipy import stats


def _root() -> Path:
    # Delegates to kdss_paths, which finds whichever layout this file is inside
    # -- the development tree or a deposit -- by walking up from its own
    # location. It previously tried two absolute paths and exited if neither
    # existed, which on the author's machine meant a copy of the deposit read
    # and wrote the development tree instead of itself.
    return project_root()


ROOT = _root()
RUNS = dev_dir("20_new_runs")
PROD = under("20_new_runs", "g4_verification_2026-08-21/full_scores_MERGED.csv")
E8 = under("20_new_runs", "np_library_funnel_2026-08-15/e8_np_library_prioritization.csv")
ON_POCKET_CUTOFF = 8.0

# released claims, ADR-AIIM-133
CLAIM_N_OK, CLAIM_N_ONPOCKET, CLAIM_N_FAILED = 497, 482, 3
CLAIM_TIERS = {"A": 3, "B": 119, "C": 360, "off_pocket": 15, "failed": 3}


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default=str(under("20_new_runs", "smiles_redock_500_2026-08-26")))
    args = ap.parse_args()
    run = Path(args.run)

    files = sorted(glob.glob(str(run / "scores_R2_chunk*.csv")))
    if not files:
        raise SystemExit(f"no R2 score files in {run}")
    r2 = pd.concat([pd.read_csv(f) for f in files]).drop_duplicates("compound_id")
    prod = pd.read_csv(PROD)
    e8 = pd.read_csv(E8)

    tierA_released = sorted(e8.loc[e8.tier == "A_converging_sbdd_lbdd",
                                   "compound_id"].tolist())
    # the ligand-based flag is graph-derived and cannot move when only the
    # SMILES string changes; it is carried across rather than recomputed
    flagged = set(e8.loc[e8.tier.isin(["A_converging_sbdd_lbdd"]), "compound_id"])

    m = prod.merge(r2, on="compound_id", suffixes=("_prod", "_r2"))
    both = m[(m.status_prod == "ok") & (m.status_r2 == "ok")].copy()
    both["delta"] = both.vina_score_r2 - both.vina_score_prod

    n_ok_r2 = int((r2.status == "ok").sum())
    n_failed_r2 = int((r2.status != "ok").sum())
    onpocket_r2 = r2[(r2.status == "ok") & (r2.centroid_dist_A <= ON_POCKET_CUTOFF)]
    n_onpocket_r2 = int(len(onpocket_r2))

    # tiers rebuilt under the R2 scores by the same rule as the released run
    q25 = float(onpocket_r2.vina_score.quantile(0.25))
    tiers_r2 = {"A": 0, "B": 0, "C": 0,
                "off_pocket": int(((r2.status == "ok") &
                                   (r2.centroid_dist_A > ON_POCKET_CUTOFF)).sum()),
                "failed": n_failed_r2}
    tierA_r2 = []
    for _, r in onpocket_r2.iterrows():
        if r.compound_id in flagged:
            tiers_r2["A"] += 1
            tierA_r2.append(r.compound_id)
        elif r.vina_score <= q25:
            tiers_r2["B"] += 1
        else:
            tiers_r2["C"] += 1

    rank_prod = both.sort_values("vina_score_prod").reset_index(drop=True)
    rank_r2 = both.sort_values("vina_score_r2").reset_index(drop=True)
    pos_prod = {c: i + 1 for i, c in enumerate(rank_prod.compound_id)}
    pos_r2 = {c: i + 1 for i, c in enumerate(rank_r2.compound_id)}
    movement = {c: {"rank_production": pos_prod.get(c), "rank_r2": pos_r2.get(c),
                    "moved": (None if c not in pos_prod or c not in pos_r2
                              else pos_r2[c] - pos_prod[c])}
                for c in tierA_released}

    strongest_r2 = rank_r2.compound_id.iloc[0]
    strongest_in_tierA = bool(strongest_r2 in set(tierA_r2))

    claims = {
        "tier_a_identity": {
            "released": tierA_released, "r2": sorted(tierA_r2),
            "unchanged": sorted(tierA_r2) == tierA_released},
        "n_scored": {"released": CLAIM_N_OK, "r2": n_ok_r2,
                     "unchanged": n_ok_r2 == CLAIM_N_OK},
        "n_on_pocket": {"released": CLAIM_N_ONPOCKET, "r2": n_onpocket_r2,
                        "unchanged": n_onpocket_r2 == CLAIM_N_ONPOCKET},
        "n_failed": {"released": CLAIM_N_FAILED, "r2": n_failed_r2,
                     "unchanged": n_failed_r2 == CLAIM_N_FAILED},
        "tier_counts": {"released": CLAIM_TIERS, "r2": tiers_r2,
                        "unchanged": tiers_r2 == CLAIM_TIERS},
        "strongest_raw_score_not_in_tier_a": {
            "released": True, "r2": not strongest_in_tierA,
            "unchanged": not strongest_in_tierA},
    }
    all_unchanged = all(c["unchanged"] for c in claims.values())

    d = both.delta.to_numpy(float)
    descriptive = {
        "n_paired": int(len(both)),
        "mean_abs_deviation": round(float(np.abs(d).mean()), 4),
        "max_abs_deviation": round(float(np.abs(d).max()), 4),
        "n_identical": int((d == 0).sum()),
        "spearman": round(float(stats.spearmanr(both.vina_score_prod,
                                                both.vina_score_r2).statistic), 4),
        "on_pocket_flips": int(((prod.set_index("compound_id").centroid_dist_A
                                 <= ON_POCKET_CUTOFF).reindex(r2.compound_id).fillna(False).to_numpy()
                                != (r2.set_index("compound_id").centroid_dist_A
                                    <= ON_POCKET_CUTOFF).to_numpy()).sum()),
        "tier_a_rank_movement": movement,
        "comparison_on_25_compounds": {"mean": 0.166, "max": 0.974},
        "note": "Reported without a verdict attached, per ADR-AIIM-133.",
    }

    report = {"pre_registration": "ADR-AIIM-133",
              "released_claims": claims,
              "all_released_claims_unchanged": all_unchanged,
              "verdict": ("STABLE -- no released claim changes under "
                          "re-representation" if all_unchanged else
                          "CHANGED -- at least one released claim differs; "
                          "§3.6 must report it"),
              "descriptive": descriptive,
              "inputs": {"production": sha256(PROD),
                         "r2_scores": [sha256(Path(f)) for f in files]}}
    (run / "r2_results.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    print("released claims under re-representation:")
    for k, v in claims.items():
        print(f"  {'OK ' if v['unchanged'] else '!! '} {k}: "
              f"released={v['released']}  r2={v['r2']}")
    print(f"\n{report['verdict']}\n")
    print(f"descriptive: mean |delta| {descriptive['mean_abs_deviation']:.4f}  "
          f"max {descriptive['max_abs_deviation']:.4f}  "
          f"identical {descriptive['n_identical']}/{descriptive['n_paired']}  "
          f"rho {descriptive['spearman']:.4f}")
    for c, v in movement.items():
        print(f"  {c}: rank {v['rank_production']} -> {v['rank_r2']} "
              f"({v['moved']:+d})" if v["moved"] is not None else f"  {c}: not paired")
    print(f"\nwritten -> {run / 'r2_results.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

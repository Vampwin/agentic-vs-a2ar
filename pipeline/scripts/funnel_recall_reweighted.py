"""funnel_recall_reweighted.py -- what the funnel discards, weighted correctly.

The exclusion control sampled 250 compounds from each of the two ways the funnel
rejects and reported one pooled figure: 7.0 % of rejected compounds would have
reached the top-500 cutoff. Review comment M9 objects that pooling an equal
sample from unequal strata is not an unbiased estimate of the population rate,
and that the absolute number was never stated. Both are arithmetic on data
already collected; neither needs a new docking run (ADR-AIIM-150).

**The strata.** The funnel rejects twice, and the two rejections are very
different sizes:

    arm_2d   21,882 -> 8,000, so 13,882 rejected by the 2D pharmacophore/ErG screen
    arm_3d    8,000 -> 5,000, so  3,000 rejected by the 3D shape refinement

250 were sampled from each. A pooled rate over those 500 weights the smaller
stratum 4.6 times too heavily. The population estimate weights each arm's rate
by its own size.

**A correction to the review comment itself.** M9 describes "strata of ~16,882
and ~3,000". 16,882 is the *total* rejected, 13,882 + 3,000; it is not a
stratum. Its arithmetic -- 7 % of 16,882 -- therefore applies a pooled rate to a
population that is the sum of both strata, which is two errors that partly
cancel. The figures below do it once, correctly.

**What this does not do.** The control docked its sample under the *fast pass*
-- two conformers, exhaustiveness 8 -- the same protocol that chose the top 500,
which is the right comparison for asking whether the funnel's own cut was sound.
It is not the full protocol the 500 finally received. R4 answers that separate
question and is not superseded by this.

Usage:
    python funnel_recall_reweighted.py [--out DIR]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from kdss_paths import project_root, resolve, under, dev_dir

import pandas as pd

# --- pre-registered before the run, from the funnel's own definition --------
# Stratum sizes are the funnel's stage boundaries, not estimates.
STRATA = {
    "arm_2d": {"population": 21882 - 8000,
               "rejected_by": "2D pharmacophore / ErG screen, 21,882 -> 8,000"},
    "arm_3d": {"population": 8000 - 5000,
               "rejected_by": "3D shape / pharmacophore refinement, 8,000 -> 5,000"},
}
N_SELECTED = 500          # the compounds that received the full protocol
Z = 1.959963984540054     # two-sided normal quantile at 0.95


def _root() -> Path:
    # Delegates to kdss_paths, which finds whichever layout this file is inside
    # -- the development tree or a deposit -- by walking up from its own
    # location. It previously tried two absolute paths and exited if neither
    # existed, which on the author's machine meant a copy of the deposit read
    # and wrote the development tree instead of itself.
    return project_root()


ROOT = _root()
RUNS = dev_dir("20_new_runs")
CTRL = under("20_new_runs", "funnel_exclusion_control_2026-08-21")


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def wilson(k: int, n: int) -> tuple[float, float]:
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + Z * Z / n
    c = (p + Z * Z / (2 * n)) / d
    h = Z * math.sqrt(p * (1 - p) / n + Z * Z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(under("20_new_runs", "funnel_recall_reweighted_2026-08-31")))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    ledger = pd.read_csv(CTRL / "exclusion_control_ledger.csv", low_memory=False)
    scores = pd.read_csv(CTRL / "exclusion_scores.csv", low_memory=False)
    summary = json.loads((CTRL / "exclusion_control_summary.json").read_text(
        encoding="utf-8"))
    cutoff = summary["top500_cutoff"]

    d = ledger[["compound_id", "arm"]].merge(
        scores[["compound_id", "status", "vina_score"]], on="compound_id")
    ok = d[(d.status == "ok") & d.vina_score.notna()].copy()
    ok["reaches"] = ok.vina_score <= cutoff       # more negative is better

    by_arm, weighted_rate, weighted_abs = {}, 0.0, 0.0
    total_population = sum(s["population"] for s in STRATA.values())
    for arm, meta in STRATA.items():
        sub = ok[ok.arm == arm]
        n, k = len(sub), int(sub.reaches.sum())
        rate = k / n if n else float("nan")
        lo, hi = wilson(k, n)
        pop = meta["population"]
        by_arm[arm] = {
            "rejected_by": meta["rejected_by"],
            "population": pop,
            "n_sampled": n,
            "n_docked_ok": n,
            "n_reaching_cutoff": k,
            "rate": round(rate, 4),
            "wilson_95_ci": [round(lo, 4), round(hi, 4)],
            "implied_population_count": round(rate * pop),
            "weight_in_population": round(pop / total_population, 4),
        }
        weighted_rate += rate * pop / total_population
        weighted_abs += rate * pop

    pooled_n = len(ok)
    pooled_k = int(ok.reaches.sum())
    pooled_rate = pooled_k / pooled_n

    report = {
        "question": "Of the compounds the expansion funnel discarded, what "
                    "fraction would have scored well enough to enter the "
                    "top-500 -- estimated over the population rather than over "
                    "an equal sample from unequal strata?",
        "authority": "ADR-AIIM-150",
        "closes": "M9, the reweighting and the absolute number. R4 remains open "
                  "and answers the different question of recall at full protocol.",
        "sources": {
            "ledger": {"file": str((CTRL / "exclusion_control_ledger.csv")
                                   .relative_to(ROOT).as_posix()),
                       "sha256": sha256(CTRL / "exclusion_control_ledger.csv")},
            "scores": {"file": str((CTRL / "exclusion_scores.csv")
                                   .relative_to(ROOT).as_posix()),
                       "sha256": sha256(CTRL / "exclusion_scores.csv")},
        },
        "cutoff_kcal_per_mol": cutoff,
        "protocol": "fast pass -- two conformers, exhaustiveness 8 -- identical "
                    "to the screen that chose the top 500. Not the full protocol "
                    "those 500 later received.",
        "as_published_pooled": {
            "n": pooled_n,
            "n_reaching_cutoff": pooled_k,
            "rate": round(pooled_rate, 4),
            "why_it_is_biased": "250 were drawn from each arm, but the arms hold "
                                "13,882 and 3,000 compounds. The pooled figure "
                                "weights the smaller stratum 4.6 times too "
                                "heavily.",
        },
        "by_stratum": by_arm,
        "population_estimate": {
            "total_rejected": total_population,
            "weighted_rate": round(weighted_rate, 4),
            "implied_count": round(weighted_abs),
            "against_n_selected": N_SELECTED,
            "ratio_to_selected": round(weighted_abs / N_SELECTED, 2),
        },
        "review_comment_correction": (
            "M9 describes strata of about 16,882 and 3,000. 16,882 is the total "
            "rejected -- 13,882 plus 3,000 -- and is not a stratum. Applying the "
            "pooled 7.0 % to it compounds two errors that partly cancel."),
    }

    dest = out / "funnel_recall_reweighted.json"
    dest.write_text(json.dumps(report, indent=2), encoding="utf-8")

    print(f"cutoff {cutoff:.3f} kcal/mol   (fast-pass protocol)\n")
    print(f"{'stratum':<10}{'population':>12}{'sampled':>9}{'reach':>7}"
          f"{'rate':>8}{'95% CI':>18}{'implied':>10}")
    for arm, a in by_arm.items():
        ci = f"{a['wilson_95_ci'][0]:.3f}-{a['wilson_95_ci'][1]:.3f}"
        print(f"{arm:<10}{a['population']:>12,}{a['n_sampled']:>9}"
              f"{a['n_reaching_cutoff']:>7}{a['rate']:>8.4f}{ci:>18}"
              f"{a['implied_population_count']:>10,}")
    pe = report["population_estimate"]
    print(f"\npooled as published : {pooled_rate:.4f}  ({pooled_k}/{pooled_n})")
    print(f"population-weighted : {pe['weighted_rate']:.4f}")
    print(f"implied count       : {pe['implied_count']:,} of "
          f"{pe['total_rejected']:,} rejected, against {N_SELECTED} selected "
          f"({pe['ratio_to_selected']}x)")
    print(f"\nwritten -> {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

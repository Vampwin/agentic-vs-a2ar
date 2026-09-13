"""analyse_r4_funnel_recall.py -- funnel recall at the full protocol.

Criteria are module constants, transcribed from ADR-AIIM-132 and ADR-AIIM-150
**before the docking they judge**, which is the convention this project uses for
R2 and R3 and the reason their outcomes could be reported without argument.

ADR-AIIM-150 measured the discarded compounds at the fast pass and weighted them
by stratum: 7.41 %, about 1,251 of 16,882, against the 500 selected. The fast
pass is the right comparison for whether the funnel's own cut was sound, because
the fast pass made the cut. It cannot say what those compounds would have scored
had they been given the treatment the kept compounds received.

This measures that. Same 500 sampled compounds, same seeded draw, no re-draw;
8 conformers and exhaustiveness 32 instead of the fast pass. The stratum weights
and the cutoff are unchanged.

**What would falsify the paper's claim.** The manuscript argues the funnel is a
compute-allocation device with modest discriminative power, not a filter that
finds where the good scores are. If the full protocol collapsed the discarded
arm's success rate toward zero, that claim would be wrong and the funnel would
be doing real work. The rate that would say so is fixed below at 1 %, an order
of magnitude under the fast-pass figure, so the outcome cannot be argued into
agreement afterwards.

Usage:
    python analyse_r4_funnel_recall.py [--out DIR]
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
import math
from pathlib import Path
from kdss_paths import project_root, resolve, under, dev_dir

import pandas as pd

# --- pre-registered, ADR-AIIM-132 and ADR-AIIM-150 --------------------------
STRATA = {"arm_2d": 21882 - 8000, "arm_3d": 8000 - 5000}
N_SELECTED = 500
CUTOFF_SOURCE = "top500_cutoff in exclusion_control_summary.json (fast pass)"
# The funnel's discriminative power as released, kept rate / rejected rate.
RELEASED_ENRICHMENT = 1.43
# Below this weighted rate the funnel would be doing real work and the
# manuscript's characterisation of it would be wrong.
RATE_FUNNEL_VINDICATED = 0.01
# Above this the fast-pass and full-protocol pictures agree and nothing changes.
RATE_UNCHANGED_FLOOR = 0.05
Z = 1.959963984540054


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
R4 = under("20_new_runs", "funnel_recall_fullprotocol_2026-09-01")
FASTPASS = under("20_new_runs", "funnel_recall_reweighted_2026-08-31/funnel_recall_reweighted.json")


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
    ap.add_argument("--out", default=str(R4))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    files = sorted(glob.glob(str(R4 / "scores_R4_chunk*.csv")))
    if not files:
        raise SystemExit("no R4 score files; the docking has not run")
    scores = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    ledger = pd.read_csv(CTRL / "exclusion_control_ledger.csv")
    summary = json.loads((CTRL / "exclusion_control_summary.json").read_text(
        encoding="utf-8"))
    cutoff = summary["top500_cutoff"]
    fast = json.loads(FASTPASS.read_text(encoding="utf-8"))

    d = ledger[["compound_id", "arm"]].merge(
        scores[["compound_id", "status", "vina_score"]], on="compound_id", how="left")
    n_missing = int(d.vina_score.isna().sum())
    ok = d[(d.status == "ok") & d.vina_score.notna()].copy()
    ok["reaches"] = ok.vina_score <= cutoff

    by_arm, weighted_rate, weighted_abs = {}, 0.0, 0.0
    total_pop = sum(STRATA.values())
    for arm, pop in STRATA.items():
        sub = ok[ok.arm == arm]
        n, k = len(sub), int(sub.reaches.sum())
        rate = k / n if n else float("nan")
        lo, hi = wilson(k, n)
        by_arm[arm] = {
            "population": pop, "n_docked": n, "n_reaching_cutoff": k,
            "rate": round(rate, 4), "wilson_95_ci": [round(lo, 4), round(hi, 4)],
            "implied_population_count": round(rate * pop),
            "fastpass_rate": fast["by_stratum"][arm]["rate"],
        }
        weighted_rate += rate * pop / total_pop
        weighted_abs += rate * pop

    if weighted_rate <= RATE_FUNNEL_VINDICATED:
        outcome = "FUNNEL_DISCRIMINATES"
    elif weighted_rate >= RATE_UNCHANGED_FLOOR:
        outcome = "UNCHANGED"
    else:
        outcome = "REDUCED"

    report = {
        "pre_registration": {
            "criteria_fixed": "before the R4 docking, from ADR-AIIM-132 and "
                              "ADR-AIIM-150",
            "rate_funnel_vindicated_at_or_below": RATE_FUNNEL_VINDICATED,
            "rate_unchanged_at_or_above": RATE_UNCHANGED_FLOOR,
            "strata": STRATA,
            "cutoff_source": CUTOFF_SOURCE,
        },
        "protocol": "8 conformers, exhaustiveness 32, 10 poses — the protocol "
                    "the 500 selected compounds received. The fast-pass figures "
                    "are ADR-AIIM-150's and are not superseded; they answer a "
                    "different question.",
        "cutoff_kcal_per_mol": cutoff,
        "inputs": {
            "ledger": {"file": str((CTRL / "exclusion_control_ledger.csv")
                                   .relative_to(ROOT).as_posix()),
                       "sha256": sha256(CTRL / "exclusion_control_ledger.csv")},
            "n_docked": len(ok), "n_missing_or_failed": n_missing + int((d.status != "ok").sum()),
        },
        "by_stratum": by_arm,
        "population_estimate": {
            "total_rejected": total_pop,
            "weighted_rate": round(weighted_rate, 4),
            "implied_count": round(weighted_abs),
            "against_n_selected": N_SELECTED,
            "ratio_to_selected": round(weighted_abs / N_SELECTED, 2),
        },
        "fastpass_comparison": {
            "weighted_rate": fast["population_estimate"]["weighted_rate"],
            "implied_count": fast["population_estimate"]["implied_count"],
        },
        "outcome": outcome,
    }
    dest = out / "r4_results.json"
    dest.write_text(json.dumps(report, indent=2), encoding="utf-8")

    print(f"cutoff {cutoff:.3f} kcal/mol   full protocol   "
          f"{len(ok)} docked, {report['inputs']['n_missing_or_failed']} not\n")
    print(f"{'stratum':<10}{'pop':>10}{'docked':>8}{'reach':>7}{'rate':>8}"
          f"{'fastpass':>10}{'implied':>10}")
    for arm, a in by_arm.items():
        print(f"{arm:<10}{a['population']:>10,}{a['n_docked']:>8}"
              f"{a['n_reaching_cutoff']:>7}{a['rate']:>8.4f}"
              f"{a['fastpass_rate']:>10.4f}{a['implied_population_count']:>10,}")
    pe = report["population_estimate"]
    print(f"\nweighted rate  full protocol {pe['weighted_rate']:.4f}   "
          f"fast pass {report['fastpass_comparison']['weighted_rate']:.4f}")
    print(f"implied count  {pe['implied_count']:,} of {pe['total_rejected']:,} "
          f"against {N_SELECTED} selected ({pe['ratio_to_selected']}x)")
    print(f"\noutcome: {outcome}")
    print(f"written -> {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

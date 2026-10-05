"""run_selection_scorecard.py -- was each screening arm licensed by the project's own rules?

Executes M1.1. The question is not whether the SBVS and LBVS arms produced good
numbers, but whether each one **met the entry and validation conditions the
project froze for it before the campaign**, and therefore whether combining
them in a dual-criterion rule was licensed.

The verdict is not this script's opinion. The frozen evaluator
`evaluators/_confidence_tier.py` is imported and run unchanged; this file only
supplies the measured values and prints what that evaluator returns. The
thresholds it applies are the ones in the module:

    pose RMSD pass            2.0 A
    similarity coverage       0.50        scaffold coverage   0.10
    AD coverage  high 0.80    low 0.40
    CV AUC       high 0.75    low 0.60
    external LBDD AUC min     0.60
    training positives        >= 30 recommended, >= 10 floor

Two cases are scored separately, because the decision was taken twice:

  A. the original 236-compound SNP-Lib
  B. the expanded library, as represented by the 497 redocked compounds

Case B is additionally scored under each descriptor space that the project's
own applicability-domain sweep reports, because **the choice of gating space is
the decision under examination** and the spaces disagree from 0 % to 86 %.

Sources, all frozen artifacts:
  redocking control RMSD            3.2 / ADR-AIIM-037 chain
  control representativeness        20_new_runs/control_representativeness_validation_2026-08-04
  AD coverage, 236                  ADR-AIIM-045/047 (0 of 236)
  AD coverage by space, expansion   20_new_runs/g6_selection_2026-09-02
  model CV / external AUC, n_pos    20_new_runs/qsar_lbdd_2026-08-07

Output -> 20_new_runs/selection_scorecard_2026-09-22/
"""
from __future__ import annotations

import argparse
import io
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kdss_paths import under  # noqa: E402

ROOT = Path(under("10_new_system", ""))
sys.path.insert(0, str(ROOT))
from evaluators._confidence_tier import (assess_lbdd, assess_sbdd,  # noqa: E402
                                         combine)

G6 = under("20_new_runs", "g6_selection_2026-09-02/G6_contract_outcome_tierA.json")

# --- measured values, each traceable to a frozen artifact --------------------
REDOCK_RMSD = 0.379          # the redocking control that passed
MODEL = {"cv_auc": 0.907, "external_auc": 0.922, "n_positive": 5268}

# The only control-representativeness measurement on record is the production
# case: similarity 0.000, scaffold 0.000 against the 236-compound library.
CTRL_236 = {"available": True, "similarity_coverage": 0.0, "scaffold_coverage": 0.0}
# No control-representativeness was recomputed after the library was expanded.
# Passing None is not a convenience: the frozen rater treats an unmeasured
# control as very_serious, which is exactly the right reading of "never
# re-evaluated".
CTRL_EXPANSION = None

AD_236 = 0.0                 # 0 of 236 in domain (ADR-AIIM-045/047)


def tier_allows_ranking(tier: str) -> bool:
    """The decision flow: LBDD output drives a ranking only at HIGH/MODERATE."""
    return tier in ("HIGH", "MODERATE")


def score(case: str, ad_coverage, ctrl, note: str) -> dict:
    sb = assess_sbdd(rmsd_angstrom=REDOCK_RMSD,
                     control_representativeness_report=ctrl)
    lb = assess_lbdd(ad_coverage_fraction=ad_coverage,
                     cv_auc=MODEL["cv_auc"], n_positive=MODEL["n_positive"],
                     external_auc=MODEL["external_auc"])
    both = combine(sb, lb)
    sb_ok = tier_allows_ranking(sb["final_tier"])
    lb_ok = tier_allows_ranking(lb["final_tier"])
    return {
        "case": case, "note": note,
        "ad_coverage_used": ad_coverage,
        "sbdd_tier": sb["final_tier"], "sbdd_licensed_to_rank": sb_ok,
        "lbdd_tier": lb["final_tier"], "lbdd_licensed_to_rank": lb_ok,
        "combined_tier": both["overall_tier"],
        "consensus_licensed": bool(sb_ok and lb_ok),
        "sbdd_downgrades": sb["downgrade_reasons"],
        "lbdd_downgrades": lb["downgrade_reasons"],
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(
        under("20_new_runs", "selection_scorecard_2026-09-22")))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    g6 = json.load(io.open(G6, encoding="utf-8"))
    by_space = g6["evidence"]["applicability_domain_multi_space"]["coverage_by_space"]

    rows = [score("A. original 236-compound SNP-Lib", AD_236, CTRL_236,
                  "the case the decision flow walked through in 2026-08-08")]

    for space, cov in sorted(by_space.items(), key=lambda kv: -kv[1]):
        rows.append(score(f"B. expansion, gated on {space}", cov, CTRL_EXPANSION,
                          "control representativeness never recomputed after "
                          "expansion"))

    licensed = [r for r in rows if r["consensus_licensed"]]
    summary = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "executes": "M1.1",
        "verdict_produced_by": "evaluators/_confidence_tier.py, imported and run "
                               "unchanged; this script supplies measured values only",
        "ranking_rule": "the decision flow licenses an arm to drive a ranking "
                        "only at HIGH or MODERATE",
        "n_configurations_scored": len(rows),
        "n_in_which_consensus_was_licensed": len(licensed),
        "control_representativeness_recomputed_after_expansion": False,
        "rows": rows,
    }
    with open(out / "selection_scorecard.json", "w", encoding="utf-8",
              newline="\n") as f:
        json.dump(summary, f, indent=1, ensure_ascii=False)
        f.write("\n")

    print("=" * 86)
    print("WAS EACH ARM LICENSED BY THE PROJECT'S OWN FROZEN RULES?  (M1.1)")
    print("=" * 86)
    print("verdicts produced by evaluators/_confidence_tier.py, run unchanged\n")
    print(f"{'configuration':40} {'AD':>7} {'SBDD':>9} {'LBDD':>9} {'consensus':>10}")
    print("-" * 86)
    for r in rows:
        ad = "n/a" if r["ad_coverage_used"] is None else f"{r['ad_coverage_used']:.3f}"
        print(f"{r['case'][:40]:40} {ad:>7} {r['sbdd_tier']:>9} "
              f"{r['lbdd_tier']:>9} {('YES' if r['consensus_licensed'] else 'no'):>10}")
    print("-" * 86)
    print(f"consensus licensed in {len(licensed)} of {len(rows)} configurations")
    print()
    print("why each arm was downgraded, case A:")
    for k in ("sbdd_downgrades", "lbdd_downgrades"):
        for d in rows[0][k]:
            print(f"   {k.split('_')[0].upper():5} {d}")
    print()
    print("NOTE: control representativeness was never recomputed after the")
    print("library was expanded, so every case B row scores it as unmeasured,")
    print("which the frozen rater treats as very_serious.")
    print()
    print(f"written -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

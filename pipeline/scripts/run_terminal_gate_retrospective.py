"""run_terminal_gate_retrospective.py -- execute the release gate that was never executed.

Executes M1.3.

`00_project_control/18_SBDD_LBDD_DECISION_FLOW.md` names G8's
`confidence_tier_assessment` as the final release decision: the check that
decides whether a ranking may be reported at all. It is defined in
`contracts/G8_reporting.yaml`, implemented in `evaluators/G8_evaluator.py`, and
unit-tested.

Across every run directory in this project, `contract_outcome` records exist
for G0 (3), G1 (3), G2 (3), G3 (3), G4 (2) and G6 (3) — and **none for G8**
(ADR-AIIM-211). The gate that authorises release has never run.

This runs it, on the expansion-library selection, using the contract's own
frozen thresholds rather than the Python module defaults, exactly as the
evaluator does in production. Two things are reported and neither is softened:

  1. the tier the gate returns, and
  2. which of G8's other nine checks could not be evaluated because their
     inputs were never assembled — which is the mechanical reason the gate was
     never run.

Measured values, each from a frozen artifact:

  redocking RMSD 0.379 A                      3.2, ADR-AIIM-037 chain
  control representativeness                  never recomputed after expansion,
                                              so passed as unmeasured
  AD coverage by descriptor space             g6_selection_2026-09-02
  CV AUC, external AUC, n positives           qsar_lbdd_2026-08-07

Output -> 20_new_runs/terminal_gate_retrospective_2026-09-22/
"""
from __future__ import annotations

import argparse
import io
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kdss_paths import under  # noqa: E402

ROOT = Path(under("10_new_system", ""))
sys.path.insert(0, str(ROOT))
from evaluators import G8_evaluator  # noqa: E402

CONTRACT = ROOT / "contracts" / "G8_reporting.yaml"
G6 = under("20_new_runs", "g6_selection_2026-09-02/G6_contract_outcome_tierA.json")

REDOCK_RMSD = 0.379
MODEL = {"cv_auc": 0.907, "external_auc": 0.922, "n_positive": 5268}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(
        under("20_new_runs", "terminal_gate_retrospective_2026-09-22")))
    ap.add_argument("--ad-space", default="ECFP4",
                    help="descriptor space used as the gating domain")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    contract = yaml.safe_load(io.open(CONTRACT, encoding="utf-8"))
    g6 = json.load(io.open(G6, encoding="utf-8"))
    by_space = g6["evidence"]["applicability_domain_multi_space"]["coverage_by_space"]

    results = {}
    for space, cov in sorted(by_space.items(), key=lambda kv: -kv[1]):
        inputs = {
            "confidence_tier_inputs": {
                "paths_attempted": ["SBDD", "LBDD"],
                "sbdd": {
                    "rmsd_angstrom": REDOCK_RMSD,
                    # Never recomputed after the library was expanded. The
                    # frozen rater treats an unmeasured control as
                    # very_serious, which is the correct reading.
                    "control_representativeness_report": None,
                },
                "lbdd": {
                    "ad_coverage_fraction": cov,
                    "cv_auc": MODEL["cv_auc"],
                    "n_positive": MODEL["n_positive"],
                    "external_auc": MODEL["external_auc"],
                },
            }
        }
        res = G8_evaluator.evaluate(contract=contract, inputs=inputs,
                                    registry=None)
        ev = res.get("evidence", {}).get("confidence_tier_assessment", {})
        results[space] = {
            "ad_coverage": cov,
            "state": res.get("state"),
            "overall_tier": ev.get("overall_tier"),
            "sbdd_tier": (ev.get("sbdd") or {}).get("final_tier"),
            "lbdd_tier": (ev.get("lbdd") or {}).get("final_tier"),
            "checks_passed": res.get("passed", []),
            "checks_warned": res.get("warned", []),
            "checks_failed": res.get("failed", []),
        }

    first = results[args.ad_space] if args.ad_space in results else \
        next(iter(results.values()))
    # What does this gate return when handed nothing at all? Only three of
    # its ten checks turned out to be input-dependent, which prompted the
    # question. Recorded because a manuscript figure depends on it.
    empty = G8_evaluator.evaluate(contract=contract, inputs={}, registry=None)
    empty_probe = {
        "inputs_supplied": "none",
        "state": empty.get("state"),
        "n_passed_on_no_evidence": len(empty.get("passed", [])),
        "passed": empty.get("passed", []),
        "warned": empty.get("warned", []),
        "failed": empty.get("failed", []),
        "evidence_of_vacuous_passes": {
            k: empty.get("evidence", {}).get(k)
            for k in empty.get("passed", [])},
        "on_pass_action": contract.get("on_pass"),
        "reading": "Each passes by vacuous truth: all zero claims are "
                   "traceable, all zero figures verified. The one check "
                   "that warns rather than passes on absent input is "
                   "confidence_tier_assessment.",
    }

    all_checks = [c.get("id") for c in contract.get("checks", [])]
    uneval = first["checks_warned"]

    summary = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "executes": "M1.3",
        "question": "What does the release gate say, and why was it never run?",
        "gate": "G8_reporting :: confidence_tier_assessment",
        "thresholds_source": "the contract's own frozen values, not the Python "
                             "module defaults",
        "contract_outcomes_on_record_before_this": {
            "G0": 3, "G1": 3, "G2": 3, "G3": 3, "G4": 2, "G6": 3, "G8": 0},
        "control_representativeness_recomputed_after_expansion": False,
        "n_g8_checks": len(all_checks),
        "checks_not_evaluable_for_want_of_inputs": uneval,
        "by_descriptor_space": results,
        "empty_input_probe": empty_probe,
        "reading": "The tier is the same in every space in which the domain "
                   "coverage falls below the frozen low cut point, and the "
                   "structure-based arm is capped regardless of space because "
                   "its precondition was never re-measured.",
    }
    with open(out / "terminal_gate_retrospective.json", "w", encoding="utf-8",
              newline="\n") as f:
        json.dump(summary, f, indent=1, ensure_ascii=False)
        f.write("\n")

    print("=" * 78)
    print("THE RELEASE GATE, RUN RETROSPECTIVELY  (M1.3)")
    print("=" * 78)
    print(f"  contract checks: {len(all_checks)}")
    print(f"  not evaluable for want of inputs: {len(uneval)} -> {uneval}")
    print()
    print(f"  {'gating space':22} {'AD':>7} {'SBDD':>9} {'LBDD':>9} {'OVERALL':>9}")
    print("  " + "-" * 60)
    for space, r in results.items():
        print(f"  {space:22} {r['ad_coverage']:>7.3f} {str(r['sbdd_tier']):>9} "
              f"{str(r['lbdd_tier']):>9} {str(r['overall_tier']):>9}")
    print()
    print(f"  contract state returned: {first['state']}")
    print(f"  (the check passes when the assessment is computed and disclosed;")
    print(f"   a low tier is a reportable outcome, not a check failure)")
    print()
    print(f"  handed NO inputs at all this contract returns "
          f"{empty_probe['state']} with "
          f"{empty_probe['n_passed_on_no_evidence']} checks passing,")
    print(f"  and its on_pass action is {contract.get('on_pass')!r}")
    print()
    print(f"written -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

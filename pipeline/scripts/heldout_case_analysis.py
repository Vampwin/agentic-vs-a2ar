"""heldout_case_analysis.py -- what the five blind cases can and cannot support.

Five fault cases were authored by a project-blind session from a fixed brief and
evaluated once against the unmodified evaluators. One was detected. The
manuscript has been reporting that as a rate beside the 96.9 % measured on the
development families, which invites a comparison the design does not support,
for two reasons this script measures rather than asserts (M4, ADR-AIIM-135).

**The interval.** One in five is 20 % with a 95 % Wilson interval running from
roughly 1 % to 62 %. An interval that wide is not a measurement of anything; it
is compatible with a system that catches almost nothing and with one that
catches most things.

**The confound.** The five cases are not distributed over the gates the
development benchmark exercised. Two of them target gates that carry no
implemented development fault family at all, so for those two the comparison is
not "known faults versus novel faults" but "a tested gate versus an untested
one". The other three target gates the development set does exercise -- and two
of those three were missed at gates that score 100 % on the development
families, which is the more uncomfortable observation and the one the case study
is actually about.

This script emits that decomposition. It reports no rate as a headline, by
design: the point of ADR-AIIM-135 is that a rate is the wrong summary of five
cases.

**Nothing here discloses a fault mechanism.** It reads the evaluation's own
per-case outcome record -- gate, contract check name, expected state, observed
state -- which is what the evaluation measured and what Section 3.1 already
reports. The faulted workdirs stay sealed under ADR-AIIM-113, and H06-H12 have
never been exposed and are not touched.

Usage:
    python heldout_case_analysis.py [--out DIR]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from kdss_paths import project_root, resolve, under, dev_dir

CONFIDENCE = 0.95
Z = 1.959963984540054   # two-sided normal quantile at 0.95


def _root() -> Path:
    # Delegates to kdss_paths, which finds whichever layout this file is inside
    # -- the development tree or a deposit -- by walking up from its own
    # location. It previously tried two absolute paths and exited if neither
    # existed, which on the author's machine meant a copy of the deposit read
    # and wrote the development tree instead of itself.
    return project_root()


ROOT = _root()
RUNS = dev_dir("20_new_runs")
HELDOUT = under("20_new_runs", "heldout_evaluation_2026-08-08/heldout_evaluation_summary.json")
RECOUNT = RUNS / ("seeded_defect_recount_2026-08-29/"
                  "mutation_score_by_gate_recounted.json")


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def wilson(k: int, n: int, z: float = Z) -> tuple[float, float]:
    """Wilson score interval. Chosen over Wald because Wald on 1/5 gives a
    lower bound below zero, which is the sort of number that gets quoted."""
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, centre - half), min(1.0, centre + half))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(under("20_new_runs", "heldout_case_analysis_2026-08-29")))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    ho = json.loads(HELDOUT.read_text(encoding="utf-8"))
    rc = json.loads(RECOUNT.read_text(encoding="utf-8"))

    dev_gates = {g: v for g, v in rc["by_gate"].items()}
    untested = set(rc["coverage_limitation"]["gates_with_no_implemented_fault_family"])

    cases = []
    for r in ho["results"]:
        gate = r["gate"]
        dev = dev_gates.get(gate)
        cases.append({
            "case": r["case"],
            "gate": gate,
            "contract_check": r["check"],
            "expected": r["expected"],
            "observed": r["actual_state"],
            "detected": bool(r.get("detected")),
            "evaluator_ran": r.get("evaluated_at_all", True),
            "gate_exercised_by_development_benchmark": gate not in untested,
            "development_detection_rate_at_this_gate":
                round(dev["detection_rate"], 4) if dev else None,
            "development_cases_at_this_gate": dev["n_fault_cases"] if dev else 0,
        })

    on_untested = [c for c in cases if not c["gate_exercised_by_development_benchmark"]]
    on_tested = [c for c in cases if c["gate_exercised_by_development_benchmark"]]
    missed_at_perfect_gate = [
        c for c in on_tested
        if not c["detected"] and c["development_detection_rate_at_this_gate"] == 1.0]

    lo, hi = wilson(ho["n_detected"], ho["n_cases"])

    report = {
        "purpose": "M4 / ADR-AIIM-135. The five blind cases are reported as a "
                   "qualitative case study. This artifact records why a rate is "
                   "the wrong summary of them, in numbers, so that the decision "
                   "not to quote one is itself checkable.",
        "authority": "ADR-AIIM-135",
        "sources": {
            "heldout_evaluation": {"file": str(HELDOUT.relative_to(ROOT).as_posix()),
                                   "sha256": sha256(HELDOUT)},
            "development_benchmark": {"file": str(RECOUNT.relative_to(ROOT).as_posix()),
                                      "sha256": sha256(RECOUNT)},
        },
        "disclosure": "Per-case gate, contract check name and observed state "
                      "only -- what the evaluation measured. No fault mechanism "
                      "is described here. The faulted workdirs remain sealed "
                      "(ADR-AIIM-113); H06-H12 have never been exposed.",
        "n_cases": ho["n_cases"],
        "n_detected": ho["n_detected"],
        "point_estimate": round(ho["n_detected"] / ho["n_cases"], 4),
        "wilson_95_ci": [round(lo, 3), round(hi, 3)],
        "ci_width": round(hi - lo, 3),
        "why_no_rate_is_reported":
            f"A 95 % Wilson interval of [{lo:.3f}, {hi:.3f}] spans "
            f"{(hi - lo) * 100:.0f} percentage points. It excludes neither a "
            f"system that detects almost nothing nor one that detects most "
            f"things, so quoting 20 % beside a development figure measured on "
            f"128 cases compares an estimate with a number.",
        "gate_coverage_confound": {
            "gates_with_no_implemented_development_family": sorted(untested),
            "n_cases_on_such_gates": len(on_untested),
            "cases_on_such_gates": [c["case"] for c in on_untested],
            "gates_of_those_cases": sorted({c["gate"] for c in on_untested}),
            "reading": "For these cases the contrast is not between anticipated "
                       "and unanticipated fault mechanisms but between a gate "
                       "the benchmark exercised and one it never did. They "
                       "cannot speak to generalisation across mechanisms.",
        },
        "the_uncomfortable_part": {
            "n_missed_at_gates_scoring_100pc_on_development":
                len(missed_at_perfect_gate),
            "cases": [{"case": c["case"], "gate": c["gate"],
                       "development_cases_at_this_gate":
                           c["development_cases_at_this_gate"]}
                      for c in missed_at_perfect_gate],
            "reading": "These are the cases the confound does not explain. The "
                       "gate was exercised by the development benchmark, scored "
                       "perfectly there, and still missed a case authored blind "
                       "to it. That is the finding the case study carries.",
        },
        "cases": cases,
    }

    dest = out / "heldout_case_analysis.json"
    dest.write_text(json.dumps(report, indent=2), encoding="utf-8")

    print(f"cases: {ho['n_cases']}, detected: {ho['n_detected']}, "
          f"point {report['point_estimate']:.2f}")
    print(f"Wilson 95% CI: [{lo:.3f}, {hi:.3f}]  width {hi - lo:.3f}")
    print()
    print(f"{'case':<6}{'gate':<6}{'det':<5}{'ran':<5}{'gate tested?':<14}"
          f"{'dev rate':>9}{'dev n':>7}")
    for c in cases:
        rate = ("--" if c["development_detection_rate_at_this_gate"] is None
                else f"{c['development_detection_rate_at_this_gate'] * 100:.1f}%")
        print(f"{c['case']:<6}{c['gate']:<6}"
              f"{'yes' if c['detected'] else 'no':<5}"
              f"{'yes' if c['evaluator_ran'] else 'NO':<5}"
              f"{'yes' if c['gate_exercised_by_development_benchmark'] else 'no':<14}"
              f"{rate:>9}{c['development_cases_at_this_gate']:>7}")
    print()
    print(f"on gates with no development family : {len(on_untested)} "
          f"({', '.join(c['case'] for c in on_untested)})")
    print(f"missed at gates scoring 100 % on dev: {len(missed_at_perfect_gate)} "
          f"({', '.join(c['case'] for c in missed_at_perfect_gate)})")
    print(f"\nwritten -> {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""analyse_r3_conformer_sweep.py -- R3, evaluated against ADR-AIIM-133.

The criteria are transcribed from the pre-registration and are constants at the
top of this file. They were fixed before any arm was docked. Nothing below
chooses a threshold after looking at a number.

Order of evaluation is the order the ADR sets:
  1. control -- the 8-conformer arms must reproduce the runs they reuse the
     ledgers of, to 0.001 kcal/mol. If either fails, nothing else is reported.
  2. primary -- d(128)/d(8), where d(n) is the mean absolute canonical-versus-
     rewritten deviation at n conformers.
  3. secondary -- TOST equivalence on the paired differences at 128 conformers,
     margin 0.096 kcal/mol, the smallest gap the reported Tier-A ordering
     depends on.
  4. monotonicity -- a separate declared prediction, reported either way.

Usage:
    python analyse_r3_conformer_sweep.py [--run DIR]
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# ---- pre-registered constants, ADR-AIIM-133 -------------------------------
CONFORMERS = [8, 32, 128]
CONTROL_TOL = 0.001          # kcal/mol, control reproduction tolerance
RATIO_ARTEFACT = 0.10        # d(128)/d(8) at or below this: sampling explains it
RATIO_NOT_EXPLAINED = 0.50   # at or above this: it does not
TOST_MARGIN = 0.096          # smallest adjacent Tier-A gap, kcal/mol
D8_PRIOR = 0.166             # mean deviation already measured at 8 conformers


def _root() -> Path:
    marker = "00_project_control/08_DECISION_LOG.md"
    for cand in (Path(r"D:\Agentic VS Workflow"), Path("/mnt/d/Agentic VS Workflow")):
        if (cand / marker).is_file():
            return cand
    raise SystemExit("project root not found")


ROOT = _root()
RUNS = ROOT / "20_new_runs"
CTRL_CANON = RUNS / "g4_verification_2026-08-21/rerun_passA_cpu6.csv"
CTRL_REWRITTEN = RUNS / "metamorphic_2026-08-21/scores_M2.csv"


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def scores(p: Path) -> pd.Series:
    d = pd.read_csv(p)
    d = d[d.get("status", "ok").astype(str) == "ok"] if "status" in d else d
    return d.set_index("compound_id")["vina_score"].astype(float)


def tost(diff: np.ndarray, margin: float, alpha: float = 0.05) -> dict:
    """Two one-sided tests. Equivalent iff both one-sided p-values < alpha."""
    n = len(diff)
    m, se = float(diff.mean()), float(diff.std(ddof=1) / np.sqrt(n))
    if se == 0.0:
        equiv = abs(m) < margin
        return {"n": n, "mean_difference": round(m, 6), "se": 0.0,
                "p_lower": 0.0 if equiv else 1.0, "p_upper": 0.0 if equiv else 1.0,
                "equivalent_at_margin": bool(equiv),
                "note": "zero variance: every paired difference identical"}
    t_lo = (m + margin) / se
    t_hi = (m - margin) / se
    p_lo = float(stats.t.sf(t_lo, n - 1))
    p_hi = float(stats.t.cdf(t_hi, n - 1))
    return {"n": n, "mean_difference": round(m, 6), "se": round(se, 6),
            "p_lower": round(p_lo, 6), "p_upper": round(p_hi, 6),
            "equivalent_at_margin": bool(p_lo < alpha and p_hi < alpha)}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default=str(RUNS / "conformer_budget_sweep_2026-08-27"))
    # The control gate decides whether the remaining arms are worth waiting
    # for, so it has to be runnable before they exist.
    ap.add_argument("--control-only", action="store_true",
                    help="evaluate the 8-conformer control gate and stop")
    args = ap.parse_args()
    run = Path(args.run)

    needed = [8] if args.control_only else CONFORMERS
    arms = {}
    for n in needed:
        for rep in ("canonical", "rewritten"):
            p = run / f"scores_{rep}_c{n}.csv"
            if not p.exists():
                raise SystemExit(f"missing arm: {p}")
            arms[(rep, n)] = scores(p)

    report = {"pre_registration": "ADR-AIIM-133",
              "criteria": {"control_tolerance_kcal": CONTROL_TOL,
                           "ratio_sampling_artefact_at_or_below": RATIO_ARTEFACT,
                           "ratio_not_explained_at_or_above": RATIO_NOT_EXPLAINED,
                           "tost_margin_kcal": TOST_MARGIN,
                           "d8_previously_measured": D8_PRIOR}}

    # ---- 1. control, evaluated first --------------------------------------
    control = {}
    for rep, ref in (("canonical", CTRL_CANON), ("rewritten", CTRL_REWRITTEN)):
        a, b = arms[(rep, 8)], scores(ref)
        common = a.index.intersection(b.index)
        dev = float(np.abs(a[common] - b[common]).max())
        control[rep] = {"reference": str(ref.relative_to(ROOT).as_posix()),
                        "reference_sha256": sha256(ref),
                        "n_compared": int(len(common)),
                        "max_abs_deviation": round(dev, 6),
                        "passes": bool(dev <= CONTROL_TOL)}
    report["control"] = control
    report["control_passed"] = all(v["passes"] for v in control.values())

    if args.control_only:
        (run / "r3_control_gate.json").write_text(
            json.dumps({k: report[k] for k in
                        ("pre_registration", "criteria", "control",
                         "control_passed")}, indent=2), encoding="utf-8")
        for rep, v in control.items():
            print(f"  {rep:>10} n={v['n_compared']}  max dev "
                  f"{v['max_abs_deviation']:.6f}  "
                  f"{'PASS' if v['passes'] else 'FAIL'}")
        print("\ncontrol gate: "
              + ("PASSED -- the remaining arms are interpretable"
                 if report["control_passed"] else
                 "FAILED -- per ADR-AIIM-133, stop; no other arm is interpretable"))
        return 0 if report["control_passed"] else 1

    if not report["control_passed"]:
        report["verdict"] = ("CONTROL FAILED -- the 8-conformer arms do not "
                             "reproduce the runs whose ledgers they reuse. Per "
                             "ADR-AIIM-133 no conclusion is drawn from any "
                             "other arm.")
        (run / "r3_results.json").write_text(json.dumps(report, indent=2),
                                             encoding="utf-8")
        print(json.dumps(report["control"], indent=2))
        print("\n" + report["verdict"])
        return 1

    # ---- 2. deviation by conformer budget ---------------------------------
    per_n, paired = {}, {}
    for n in CONFORMERS:
        a, b = arms[("canonical", n)], arms[("rewritten", n)]
        common = a.index.intersection(b.index)
        d = (a[common] - b[common]).to_numpy(float)
        paired[n] = d
        per_n[n] = {
            "n_compared": int(len(common)),
            "mean_abs_deviation": round(float(np.abs(d).mean()), 6),
            "max_abs_deviation": round(float(np.abs(d).max()), 6),
            "n_identical": int((d == 0).sum()),
            "spearman_vs_canonical": round(
                float(stats.spearmanr(a[common], b[common]).statistic), 4),
        }
    report["by_conformer_budget"] = {str(k): v for k, v in per_n.items()}

    d8, d128 = per_n[8]["mean_abs_deviation"], per_n[128]["mean_abs_deviation"]
    ratio = d128 / d8 if d8 else float("nan")
    if ratio <= RATIO_ARTEFACT:
        outcome = "SAMPLING_ARTEFACT"
    elif ratio >= RATIO_NOT_EXPLAINED:
        outcome = "NOT_EXPLAINED_BY_SAMPLING"
    else:
        outcome = "PARTIAL"
    report["primary"] = {"d8": d8, "d128": d128,
                         "ratio_d128_over_d8": round(ratio, 4),
                         "outcome": outcome}

    # ---- 3. equivalence at the decision-relevant margin -------------------
    report["secondary_tost_at_128"] = tost(paired[128], TOST_MARGIN)

    # ---- 4. monotonicity ---------------------------------------------------
    means = [per_n[n]["mean_abs_deviation"] for n in CONFORMERS]
    maxes = [per_n[n]["max_abs_deviation"] for n in CONFORMERS]
    report["monotonicity"] = {
        "mean_abs_deviation_by_budget": dict(zip(map(str, CONFORMERS), means)),
        "max_abs_deviation_by_budget": dict(zip(map(str, CONFORMERS), maxes)),
        "mean_monotonic_decreasing": bool(all(x >= y for x, y in zip(means, means[1:]))),
        "max_monotonic_decreasing": bool(all(x >= y for x, y in zip(maxes, maxes[1:]))),
    }

    report["inputs"] = {f"{rep}_c{n}": sha256(run / f"scores_{rep}_c{n}.csv")
                        for n in CONFORMERS for rep in ("canonical", "rewritten")}
    (run / "r3_results.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    print("control:")
    for rep, v in control.items():
        print(f"  {rep:>10} max dev {v['max_abs_deviation']:.6f}  "
              f"{'PASS' if v['passes'] else 'FAIL'}")
    print("\ncanonical vs rewritten, by conformer budget:")
    for n in CONFORMERS:
        v = per_n[n]
        print(f"  c{n:<4} mean {v['mean_abs_deviation']:.4f}  max {v['max_abs_deviation']:.4f}  "
              f"identical {v['n_identical']}/{v['n_compared']}  rho {v['spearman_vs_canonical']:.4f}")
    print(f"\nprimary: d(128)/d(8) = {ratio:.4f}  ->  {outcome}")
    t = report["secondary_tost_at_128"]
    print(f"secondary: TOST margin {TOST_MARGIN} -> "
          f"{'EQUIVALENT' if t['equivalent_at_margin'] else 'NOT equivalent'} "
          f"(mean diff {t['mean_difference']:+.5f})")
    m = report["monotonicity"]
    print(f"monotonic (mean): {m['mean_monotonic_decreasing']}   "
          f"(max): {m['max_monotonic_decreasing']}")
    print(f"\nwritten -> {run / 'r3_results.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

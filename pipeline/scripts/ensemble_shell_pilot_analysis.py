"""
Pilot analysis: does ensemble docking across water shells recover the
discriminative power that single-reference IFP rescoring could not
(ADR-AIIM-075, ADR-AIIM-076)?

Pre-registered before the pilot dock was run, in
20_new_runs/ensemble_shell_pilot_2026-08-15/PREREGISTRATION.json:

  subset        all 50 successfully-docked benchmark actives + 50 inactives
                sampled at random_state=20260815 (n=100, ledger SHA-256 recorded)
  protocol      the same subset re-docked at the 3 A water shell (4 waters),
                identical box, seed and exhaustiveness as the 5 A run
  comparison    paired on the identical subset -- AUC(cut5) vs AUC(cut3) vs
                AUC(ensemble = best Vina energy across both shells)
  decision      pursue full ensemble docking only if
                AUC(ensemble) - AUC(cut5) >= +0.05; otherwise stop the SBDD
                rescue and report as null

The decision rule is applied mechanically here rather than judged after the
fact, which is the whole point of having declared it first.

Usage:
    python ensemble_shell_pilot_analysis.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit.ML.Scoring import Scoring

ROOT = Path(r"D:\Agentic VS Workflow")
PILOT = ROOT / "20_new_runs/ensemble_shell_pilot_2026-08-15"
CUT5_SCORES = ROOT / "20_new_runs/chembl245_ifp_rescoring_2026-08-15/scores.csv"

BOOTSTRAP_N, BOOTSTRAP_SEED = 1000, 42
EF_FRACTIONS = [0.01, 0.05, 0.10]
BEDROC_ALPHA = 20.0
DECISION_THRESHOLD = 0.05


def _rows(score_label):
    return [[int(l)] for _, l in sorted(score_label, key=lambda t: -t[0])]


def metrics(score_label):
    rows = _rows(score_label)
    n_pos = sum(r[0] for r in rows)
    out = {"n": len(rows), "n_active": n_pos, "n_inactive": len(rows) - n_pos}
    if n_pos in (0, len(rows)):
        out["error"] = "degenerate labels"
        return out
    out["auc_roc"] = round(float(Scoring.CalcAUC(rows, 0)), 4)
    out["bedroc20"] = round(float(Scoring.CalcBEDROC(rows, 0, BEDROC_ALPHA)), 4)
    for f, v in zip(EF_FRACTIONS, Scoring.CalcEnrichment(rows, 0, EF_FRACTIONS)):
        out[f"ef_{int(f*100)}pct"] = round(float(v), 4)
    return out


def boot_ci(score_label, key="auc_roc", n_boot=BOOTSTRAP_N, seed=BOOTSTRAP_SEED):
    rng = np.random.default_rng(seed)
    arr = np.array(score_label, dtype=object)
    n, vals = len(arr), []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        s = [tuple(arr[i]) for i in idx]
        if len({l for _, l in s}) < 2:
            continue
        m = metrics(s)
        if key in m:
            vals.append(m[key])
    if not vals:
        return None
    return [round(float(np.percentile(vals, 2.5)), 4),
            round(float(np.percentile(vals, 97.5)), 4)]


def paired_delta_ci(sl_a, sl_b, n_boot=BOOTSTRAP_N, seed=BOOTSTRAP_SEED):
    """Bootstrap CI on AUC(b) - AUC(a), resampling compounds jointly."""
    rng = np.random.default_rng(seed)
    a = np.array(sl_a, dtype=object)
    b = np.array(sl_b, dtype=object)
    n, deltas = len(a), []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        sa = [tuple(a[i]) for i in idx]
        sb = [tuple(b[i]) for i in idx]
        if len({l for _, l in sa}) < 2:
            continue
        ma, mb = metrics(sa), metrics(sb)
        if "auc_roc" in ma and "auc_roc" in mb:
            deltas.append(mb["auc_roc"] - ma["auc_roc"])
    if not deltas:
        return None
    return [round(float(np.percentile(deltas, 2.5)), 4),
            round(float(np.percentile(deltas, 97.5)), 4)]


def main() -> int:
    prereg = json.loads((PILOT / "PREREGISTRATION.json").read_text())
    ledger = pd.read_csv(PILOT / "pilot_ledger.csv")
    labels = dict(zip(ledger["compound_id"].astype(str), ledger["label"].astype(int)))

    cut3 = pd.read_csv(PILOT / "scores_cut3.csv")
    cut3 = cut3[cut3["status"] == "ok"][["compound_id", "vina_score"]]
    cut3.columns = ["compound_id", "vina_cut3"]

    cut5 = pd.read_csv(CUT5_SCORES)
    cut5 = cut5[cut5["status"] == "ok"][["compound_id", "vina_score"]]
    cut5.columns = ["compound_id", "vina_cut5"]

    df = cut3.merge(cut5, on="compound_id", how="inner")
    df["compound_id"] = df["compound_id"].astype(str)
    df = df[df["compound_id"].isin(labels)]
    df["label"] = df["compound_id"].map(labels)
    # Ensemble: best (most negative) energy across the two receptor states.
    df["vina_ensemble"] = df[["vina_cut3", "vina_cut5"]].min(axis=1)

    print(f"pilot subset docked at both shells: n={len(df)} "
          f"({int(df['label'].sum())} active / {int((df['label'] == 0).sum())} inactive)")
    print(f"pre-registered decision rule: "
          f"{prereg['decision_rule_declared_before_running']}\n")
    if df.empty:
        print("no overlapping compounds yet -- rerun once the pilot dock finishes")
        return 0

    methods = {
        "cut5_only (5A, 15 waters)": -df["vina_cut5"],
        "cut3_only (3A, 4 waters)": -df["vina_cut3"],
        "ensemble_best_of_both": -df["vina_ensemble"],
    }

    print(f"{'method':>28} {'AUC':>7} {'95% CI':>16} {'BEDROC20':>9} {'EF@5%':>7}")
    print("-" * 74)
    summary, sls = {}, {}
    for name, col in methods.items():
        sl = list(zip(col.tolist(), df["label"].tolist()))
        sls[name] = sl
        m = metrics(sl)
        m["auc_ci95"] = boot_ci(sl)
        summary[name] = m
        ci = m["auc_ci95"]
        cis = f"[{ci[0]:.3f}, {ci[1]:.3f}]" if ci else "n/a"
        print(f"{name:>28} {m['auc_roc']:>7.3f} {cis:>16} {m['bedroc20']:>9.3f} "
              f"{m.get('ef_5pct', float('nan')):>7.2f}")

    base = summary["cut5_only (5A, 15 waters)"]["auc_roc"]
    ens = summary["ensemble_best_of_both"]["auc_roc"]
    delta = ens - base
    dci = paired_delta_ci(sls["cut5_only (5A, 15 waters)"], sls["ensemble_best_of_both"])

    print(f"\npaired delta AUC (ensemble - cut5): {delta:+.4f}")
    if dci:
        print(f"  95% CI on the paired difference: [{dci[0]:+.4f}, {dci[1]:+.4f}]")
    print(f"  pre-registered threshold: >= +{DECISION_THRESHOLD}")

    meets = delta >= DECISION_THRESHOLD
    verdict = ("PURSUE full ensemble docking" if meets
               else "STOP the SBDD rescue -- report as null and move to the LBDD/QSAR arm")
    print(f"\n  DECISION (rule applied mechanically): {verdict}")

    out = {
        "preregistration": prereg,
        "n": len(df),
        "methods": summary,
        "paired_delta_auc_ensemble_minus_cut5": round(delta, 4),
        "paired_delta_ci95": dci,
        "decision_threshold": DECISION_THRESHOLD,
        "threshold_met": bool(meets),
        "decision": verdict,
    }
    (PILOT / "pilot_result.json").write_text(json.dumps(out, indent=2, default=str))
    df.to_csv(PILOT / "pilot_per_compound.csv", index=False)
    print(f"\nwritten -> {PILOT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

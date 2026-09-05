"""qsar_calibration_check.py -- m4 and m5.

Two loose ends from the internal review, both answerable from artifacts that
already exist.

m4. The manuscript quotes three different AUCs for the ligand-based model --
0.909, 0.911 and 0.922 -- without saying that they measure three different
things. A fourth, 0.907, exists in the run summary and is not quoted at all.
None is wrong; the reader is simply given no way to tell them apart. This
script assembles them from their own artifacts with what each one is.

m5. Table 3 reports a calibrated probability for each candidate and the paper
never says how well calibrated it is. Brier score and a reliability curve are
computed here on the same scaffold-disjoint held-out set the external AUC comes
from, using the per-compound probabilities recovered in the C3 run
(ADR-AIIM-137) so the model does not have to be reloaded.

Usage:
    python qsar_calibration_check.py [--out DIR]
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import brier_score_loss, roc_auc_score


def _root() -> Path:
    marker = "00_project_control/08_DECISION_LOG.md"
    for cand in (Path(r"D:\Agentic VS Workflow"), Path("/mnt/d/Agentic VS Workflow")):
        if (cand / marker).is_file():
            return cand
    raise SystemExit("project root not found")


ROOT = _root()
RUNS = ROOT / "20_new_runs"
QSAR = RUNS / "qsar_lbdd_2026-08-07"
CDIAG = RUNS / "conformal_diagnostics_2026-08-25"
N_BINS = 10


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def reliability(y: np.ndarray, p: np.ndarray, n_bins: int = N_BINS) -> tuple[list, float, float]:
    """Equal-width reliability curve, plus ECE and MCE over occupied bins."""
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    idx = np.clip(np.digitize(p, edges[1:-1], right=False), 0, n_bins - 1)
    rows, ece, mce = [], 0.0, 0.0
    for b in range(n_bins):
        m = idx == b
        if not m.any():
            continue
        conf, acc, n = float(p[m].mean()), float(y[m].mean()), int(m.sum())
        gap = abs(conf - acc)
        ece += (n / len(p)) * gap
        mce = max(mce, gap)
        rows.append({"bin": f"[{edges[b]:.1f},{edges[b + 1]:.1f})",
                     "n": n, "mean_predicted": round(conf, 4),
                     "observed_active_fraction": round(acc, 4),
                     "gap": round(gap, 4)})
    return rows, round(ece, 4), round(mce, 4)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(RUNS / "qsar_calibration_2026-08-26"))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    # ---------------- m5: calibration on the held-out set -------------------
    cal = pd.read_csv(CDIAG / "calibration_pvalues.csv")
    y = cal["label"].to_numpy(int)
    p = cal["p_active_calibrated"].to_numpy(float)

    brier = float(brier_score_loss(y, p))
    auc = float(roc_auc_score(y, p))
    base = float(y.mean())
    # Brier of the no-skill model that always predicts the base rate. Without
    # it, 0.0x on a 90%-active set reads as excellent when it may be no better
    # than a constant.
    brier_base = float(np.mean((base - y) ** 2))
    rows, ece, mce = reliability(y, p)

    calibration = {
        "held_out_n": int(len(y)),
        "n_active": int(y.sum()),
        "n_inactive": int((y == 0).sum()),
        "base_rate_active": round(base, 4),
        "roc_auc": round(auc, 4),
        "brier_score": round(brier, 4),
        "brier_score_of_base_rate_predictor": round(brier_base, 4),
        "brier_skill_score_vs_base_rate": round(1.0 - brier / brier_base, 4),
        "expected_calibration_error": ece,
        "maximum_calibration_error": mce,
        "reliability_curve": rows,
        "caveat":
            "This is the calibration of the classifier on data drawn from its "
            "own distribution. It says nothing about the natural-product "
            "libraries, which sit outside it -- see the conformal diagnostics "
            "(ADR-AIIM-137) for what happens there. The held-out set is also "
            "90.1 % active, so the Brier score must be read against the "
            "base-rate predictor rather than against zero.",
    }

    # ---------------- m4: the AUCs, each with what it measures --------------
    s_sum = json.loads((QSAR / "qsar_lbdd_summary.json").read_text(encoding="utf-8"))
    s_cmp = json.loads((QSAR / "qsar_model_comparison.json").read_text(encoding="utf-8"))
    s_yr = json.loads((QSAR / "qsar_y_randomization.json").read_text(encoding="utf-8"))

    aucs = [
        {"value": round(s_cmp["results"]["random_forest"]["roc_auc_mean"], 4),
         "sd": round(s_cmp["results"]["random_forest"]["roc_auc_sd"], 4),
         "quoted_in_manuscript": True,
         "what_it_is": "Random forest in the four-algorithm screening comparison: "
                       "nested CV, scaffold-grouped, on the 4,753-compound training "
                       "portion, with a reduced hyperparameter search (n_iter=6) "
                       "chosen for tractability across four algorithms.",
         "source": "qsar_model_comparison.json"},
        {"value": round(s_sum["cv_results"]["roc_auc_mean"], 4),
         "sd": round(s_sum["cv_results"]["roc_auc_sd"], 4),
         "quoted_in_manuscript": False,
         "what_it_is": "The retained model's own nested CV on the same 4,753 "
                       "training compounds, from the run that produced the frozen "
                       "pipeline. Differs from the line above because it is a "
                       "separate fit with a fuller search, not a different dataset.",
         "source": "qsar_lbdd_summary.json"},
        {"value": round(s_yr["real_cv_auc"], 4),
         "sd": None,
         "quoted_in_manuscript": True,
         "what_it_is": "The same frozen model re-evaluated inside the "
                       "y-randomisation harness, so that the real model and the 100 "
                       "permuted models are scored by identical code. It is the "
                       "correct comparator for the permuted mean of 0.500 and is "
                       "not a fourth estimate of performance.",
         "source": "qsar_y_randomization.json"},
        {"value": round(s_sum["held_out_test"]["roc_auc"], 4),
         "sd": None,
         "quoted_in_manuscript": True,
         "what_it_is": "External performance on the scaffold-disjoint held-out set "
                       "(n = 1,217). This is the split that doubles as the conformal "
                       "calibration set, so it is an optimistic estimate -- already "
                       "declared in Section 2.4.",
         "source": "qsar_lbdd_summary.json"},
    ]

    report = {
        "calibration_of_the_qsar_classifier": calibration,
        "auc_reconciliation": {
            "note": "Four numbers, four measurements, none of them wrong. Three "
                    "appear in the manuscript without being distinguished; the "
                    "second appears nowhere.",
            "values": aucs,
        },
        "provenance": {
            "script": "qsar_calibration_check.py",
            "n_bins": N_BINS,
            "inputs": {
                "calibration_pvalues.csv": sha256(CDIAG / "calibration_pvalues.csv"),
                "qsar_lbdd_summary.json": sha256(QSAR / "qsar_lbdd_summary.json"),
                "qsar_model_comparison.json": sha256(QSAR / "qsar_model_comparison.json"),
                "qsar_y_randomization.json": sha256(QSAR / "qsar_y_randomization.json"),
            },
        },
    }
    (out / "qsar_calibration.json").write_text(json.dumps(report, indent=2),
                                               encoding="utf-8")
    pd.DataFrame(rows).to_csv(out / "reliability_curve.csv", index=False)

    print(f"held-out n={len(y)}  ({int(y.sum())} active, {int((y == 0).sum())} inactive, "
          f"base rate {base:.3f})")
    print(f"ROC-AUC        {auc:.4f}")
    print(f"Brier          {brier:.4f}   (base-rate predictor {brier_base:.4f}, "
          f"skill score {1 - brier / brier_base:+.4f})")
    print(f"ECE {ece:.4f}   MCE {mce:.4f}\n")
    print("reliability curve:")
    for r in rows:
        print(f"  {r['bin']:>12}  n={r['n']:>4}  predicted={r['mean_predicted']:.3f}  "
              f"observed={r['observed_active_fraction']:.3f}  gap={r['gap']:.3f}")
    print("\nAUCs:")
    for a in aucs:
        q = "quoted" if a["quoted_in_manuscript"] else "NOT quoted"
        sd = f" ± {a['sd']}" if a["sd"] else ""
        print(f"  {a['value']}{sd}  [{q}]  {a['source']}")
    print(f"\nwritten -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

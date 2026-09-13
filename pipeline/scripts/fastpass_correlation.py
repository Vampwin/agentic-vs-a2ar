"""fastpass_correlation.py -- C4. Was the top-500 selection reliable?

Every reported docking score comes from the full protocol (8 conformers,
exhaustiveness 32). But the 500 compounds that received it were chosen by a
reduced fast-pass protocol (2 conformers, exhaustiveness 8) run over 5,000. The
manuscript says the fast pass was "never used to produce a reported score",
which is true and beside the point: it chose which compounds would have one.

This script measures the agreement between the two protocols on the 497
compounds that carry both, and is explicit about what that measurement cannot
settle. The 500 are the top 500 by fast-pass score, so the correlation is
computed on a truncated range and is attenuated by construction. The
unrestricted estimate needs full-protocol scores on compounds the funnel did
not select, which is R4.

What can be settled here is how tightly the fast-pass score resolves the
selection boundary: given the residual spread of full-protocol score around
fast-pass score, how many unselected compounds sit close enough to the cutoff
that the fast pass cannot separate them from the ones that were taken.

Usage:
    python fastpass_correlation.py [--out DIR]
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
NP_RUN = under("20_new_runs", "np_library_funnel_2026-08-15")
G4_RUN = under("20_new_runs", "g4_verification_2026-08-21")
SEED = 20260825
N_BOOT = 2000


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def boot_ci(x, y, fn, n=N_BOOT, seed=SEED):
    rng = np.random.default_rng(seed)
    idx = np.arange(len(x))
    vals = [fn(x[s], y[s]) for s in (rng.choice(idx, len(idx), replace=True)
                                     for _ in range(n))]
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return round(float(lo), 4), round(float(hi), 4)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(under("20_new_runs", "fastpass_correlation_2026-08-25")))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    ledger = pd.concat([pd.read_csv(f) for f in
                        sorted(glob.glob(str(NP_RUN / "full_ledger_chunk*.csv")))])
    full = pd.read_csv(G4_RUN / "full_scores_MERGED.csv")
    fp_all = pd.read_csv(NP_RUN / "fastpass_scores_MERGED.csv")

    m = full.merge(ledger[["compound_id", "fastpass_vina_score", "funnel_rank"]],
                   on="compound_id", how="left")
    ok = m[(m.status == "ok") & m.vina_score.notna()
           & m.fastpass_vina_score.notna()].copy()
    x = ok.fastpass_vina_score.to_numpy(float)
    y = ok.vina_score.to_numpy(float)

    rho = float(stats.spearmanr(x, y).statistic)
    r = float(stats.pearsonr(x, y).statistic)
    rho_ci = boot_ci(x, y, lambda a, b: stats.spearmanr(a, b).statistic)
    r_ci = boot_ci(x, y, lambda a, b: stats.pearsonr(a, b).statistic)

    # residual spread of the full-protocol score around the fast-pass score
    slope, icept, *_ = stats.linregress(x, y)
    resid = y - (icept + slope * x)
    s_resid = float(resid.std(ddof=2))

    # the selection boundary, from the 4,999 fast-pass scores that exist
    fps = fp_all[fp_all.vina_score.notna()].vina_score.to_numpy(float)
    fps_sorted = np.sort(fps)
    cutoff = float(fps_sorted[499])           # 500th best (most negative) score
    unselected = fps_sorted[500:]
    within_1s = int((unselected <= cutoff + s_resid).sum())
    within_2s = int((unselected <= cutoff + 2 * s_resid).sum())

    report = {
        "n_with_both_scores": int(len(ok)),
        "agreement_within_the_selected_500": {
            "spearman_rho": round(rho, 4), "spearman_95CI": list(rho_ci),
            "pearson_r": round(r, 4), "pearson_95CI": list(r_ci),
            "caveat": "Range-restricted by construction: these 497 are the top "
                      "500 by fast-pass score, so both correlations are "
                      "attenuated and are lower bounds on the agreement over "
                      "the full 5,000. The unrestricted estimate needs "
                      "full-protocol scores outside the selection (R4).",
        },
        "residual_spread": {
            "regression": {"slope": round(float(slope), 4),
                           "intercept": round(float(icept), 4)},
            "residual_sd_kcal_per_mol": round(s_resid, 4),
            "meaning": "Spread of the full-protocol score around what the "
                       "fast-pass score predicts, on the selected compounds.",
        },
        "selection_boundary": {
            "fastpass_cutoff_rank500_kcal_per_mol": round(cutoff, 4),
            "fastpass_range_selected": [round(float(fps_sorted[0]), 3),
                                        round(cutoff, 3)],
            "fastpass_range_all_scored": [round(float(fps_sorted[0]), 3),
                                          round(float(fps_sorted[-1]), 3)],
            "n_scored_total": int(len(fps)),
            "n_unselected_within_1_residual_sd_of_cutoff": within_1s,
            "n_unselected_within_2_residual_sd_of_cutoff": within_2s,
            "reading":
                "The fast-pass score at the cutoff is separated from the "
                f"compounds just below it by far less than the {s_resid:.2f} "
                "kcal/mol of residual spread between the two protocols. "
                f"{within_1s} unselected compounds lie within one residual SD "
                f"of the cutoff and {within_2s} within two, so at the boundary "
                "the fast pass is not resolving which compounds deserved the "
                "full protocol. This bounds the reliability of the selection, "
                "not of any reported score.",
        },
        "provenance": {
            "script": "fastpass_correlation.py",
            "bootstrap": {"n": N_BOOT, "seed": SEED},
            "inputs": {
                "full_scores_MERGED.csv": sha256(G4_RUN / "full_scores_MERGED.csv"),
                "fastpass_scores_MERGED.csv":
                    sha256(NP_RUN / "fastpass_scores_MERGED.csv"),
            },
        },
    }

    ok[["compound_id", "funnel_rank", "fastpass_vina_score", "vina_score",
        "on_pocket"]].assign(residual=np.round(resid, 4)).to_csv(
        out / "fastpass_vs_fullprotocol.csv", index=False)
    (out / "fastpass_correlation.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")

    print(f"n with both scores: {len(ok)}")
    print(f"Spearman rho = {rho:.3f}  95% CI {rho_ci}")
    print(f"Pearson  r   = {r:.3f}  95% CI {r_ci}")
    print(f"residual SD  = {s_resid:.3f} kcal/mol")
    print(f"fast-pass cutoff at rank 500 = {cutoff:.3f}")
    print(f"unselected within 1 residual SD of the cutoff: {within_1s}")
    print(f"unselected within 2 residual SD of the cutoff: {within_2s}")
    print(f"\nwritten -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""
figure_lbdd_sbdd_performance.py -- the screening-performance figure.

Every figure in this manuscript so far depicts the governance layer. None
shows what the screening itself did, and the three results that most need
seeing are all currently numbers in a table:

  A  The pre-registered rescue attempts (Table 1) as ROC curves. The claim
     is a null one, and a null is easier to overstate in prose than in a
     picture: the IFP-alone curve sits *below* the diagonal along most of
     its length, which is the observation the text turns into its me-too
     argument. Drawing all five on one axis also shows that the two Vina
     variants are indistinguishable, which "0.603 vs 0.601" states but does
     not convey.

  B  Y-randomisation. The ligand-based model's headline AUC means nothing
     without the permuted null beside it. 100 permutations, the real model
     marked, and the gap visible rather than asserted.

  C  Early recognition. IFP-alone is the one method whose BEDROC20 beats
     baseline while its AUC is below chance. That combination is the whole
     of Section 3.3's argument and is invisible in a ROC plot, because it
     lives entirely in the first few percent of the ranking.

A tie-handling note, because it changes a reported number. The benchmark's
AUCs were computed with RDKit's `Scoring.CalcAUC`, which walks a ranked list;
a ROC curve's area is the tie-averaged quantity sklearn computes. These agree
to <0.001 for every continuous score here, but `max_ifp_tanimoto` takes only
90 distinct values across 245 compounds, and on that score the two conventions
differ: 0.4487 (RDKit, as reported in Table 1) against 0.4422 (tie-averaged).
Panel A is a ROC curve, so it is drawn and labelled with the tie-averaged
value, and the caption says so. Both are below chance and no conclusion moves.

Usage:
    python figure_lbdd_sbdd_performance.py [--out DIR]
"""
from __future__ import annotations
import sys
sys.path.insert(0, r'D:\Agentic VS Workflow\10_new_system\scripts')

import argparse
import hashlib
import json
from pathlib import Path
from kdss_paths import project_root, resolve, under, dev_dir

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import figstyle as fs
import numpy as np
import pandas as pd
from rdkit.ML.Scoring import Scoring
from sklearn.metrics import roc_auc_score, roc_curve


def _root() -> Path:
    # Delegates to kdss_paths, which finds whichever layout this file is inside
    # -- the development tree or a deposit -- by walking up from its own
    # location. It previously tried two absolute paths and exited if neither
    # existed, which on the author's machine meant a copy of the deposit read
    # and wrote the development tree instead of itself.
    return project_root()


ROOT = _root()
RUNS = dev_dir("20_new_runs")
IFP = under("20_new_runs", "chembl245_ifp_rescoring_2026-08-15/ifp_rescoring_per_compound.csv")
SUMMARY = under("20_new_runs", "chembl245_ifp_rescoring_2026-08-15/ifp_rescoring_summary.json")
YRAND = under("20_new_runs", "qsar_lbdd_2026-08-07/qsar_y_randomization.json")

# (label, column, negate) -- negate where a more negative score is better
METHODS = [
    ("Vina alone (baseline)",        "vina_score",            True,  fs.COL["vina"], fs.LW["main"]),
    ("Vina at best-IFP pose",        "vina_at_best_ifp_pose", True,  "#4a7fb5", fs.LW["sec"]),
    ("Asn253 H-bond filter + Vina",  "asn253_filtered_score", False, "#7aa9d6", fs.LW["sec"]),
    ("Equal-weight z-score ensemble", "combined_zscore",      False, "#b0c4de", fs.LW["sec"]),
    ("IFP similarity alone",         "max_ifp_tanimoto",      False, fs.COL["alert"], fs.LW["main"]),
]

fs.apply()


def sha256(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def rdkit_auc(score: np.ndarray, y: np.ndarray) -> float:
    """The benchmark's own convention, kept so the figure can state both."""
    rows = sorted(zip(score, y), key=lambda t: -t[0])
    return float(Scoring.CalcAUC([[lab] for _, lab in rows], 0))


def enrichment_curve(score: np.ndarray, y: np.ndarray):
    """Fraction of all actives recovered against fraction of the list screened."""
    order = np.argsort(-score, kind="mergesort")
    hits = np.cumsum(y[order])
    frac_screened = np.arange(1, len(y) + 1) / len(y)
    return frac_screened, hits / y.sum()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(resolve("30_manuscript/figures")))
    args = ap.parse_args()

    d = pd.read_csv(IFP)
    y = d["label"].values.astype(int)
    summary = json.loads(SUMMARY.read_text(encoding="utf-8"))["methods"]
    yr = json.loads(YRAND.read_text(encoding="utf-8"))

    fig = plt.figure(figsize=(7.48, 3.05))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.06, 1.0, 1.0], wspace=0.40,
                          left=0.062, right=0.988, top=0.80, bottom=0.17)

    # ---- A: ROC ----------------------------------------------------------
    axA = fig.add_subplot(gs[0, 0])
    axA.plot([0, 1], [0, 1], color="#999999", lw=fs.LW["ref"], ls=fs.REF_DASH, zorder=1)

    # A pointwise bootstrap band on the baseline, and only the baseline.
    # Five bands would be mush, and the claim the band has to support is not
    # "each curve has an interval" but "these curves are not distinguishable":
    # every other method lies inside the baseline's band over almost its whole
    # range. Drawing it once makes that visible; drawing it five times hides it.
    # Resampling matches the benchmark's own bootstrap (n=1000, seed 42) as
    # recorded in ifp_rescoring_summary.json.
    grid = np.linspace(0, 1, 201)
    rng = np.random.default_rng(42)
    base_s = (-d[METHODS[0][1]]).values.astype(float)
    curves = []
    for _ in range(1000):
        idx = rng.integers(0, len(y), len(y))
        if len(np.unique(y[idx])) < 2:
            continue
        f, t, _ = roc_curve(y[idx], base_s[idx])
        curves.append(np.interp(grid, f, t))
    band = np.percentile(np.vstack(curves), [2.5, 97.5], axis=0)
    axA.fill_between(grid, band[0], band[1], color=METHODS[0][3], alpha=0.13,
                     linewidth=0, zorder=1.5)

    provenance = {}
    for label, col, neg, colour, lw in METHODS:
        s = (-d[col] if neg else d[col]).values.astype(float)
        fpr, tpr, _ = roc_curve(y, s)
        a_tie = roc_auc_score(y, s)
        a_rdk = rdkit_auc(s, y)
        provenance[label] = {"column": col, "negated": neg,
                             "auc_tie_averaged": round(a_tie, 4),
                             "auc_rdkit_calcauc": round(a_rdk, 4)}
        axA.plot(fpr, tpr, color=colour, lw=lw,
                 label=f"{label} ({a_tie:.3f})", zorder=3 if lw > 1.5 else 2)
    axA.set_xlim(0, 1); axA.set_ylim(0, 1)
    axA.set_xlabel("false positive rate", fontsize=fs.FS["label"])
    axA.set_ylabel("true positive rate", fontsize=fs.FS["label"])
    axA.tick_params(labelsize=fs.FS["tick"])
    fs.panel_label(axA, "A")

    # ---- B: y-randomisation (broken axis) ---------------------------------
    sub = gs[0, 1].subgridspec(1, 2, width_ratios=[6.0, 1.0], wspace=0.10)
    axB = fig.add_subplot(sub[0, 0]); axB2 = fig.add_subplot(sub[0, 1])
    perm = np.asarray(yr["permuted_aucs"], dtype=float)
    assert perm.min() > 0.44 and perm.max() < 0.56 and 0.90 < yr["real_cv_auc"] < 0.92
    axB.hist(perm, bins=np.arange(0.44, 0.5601, 0.005), color="#b0b8c4", edgecolor="white", linewidth=fs.LW["bar_edge"])
    axB.axvline(0.5, color="#999999", lw=fs.LW["ref"], ls=fs.REF_DASH)
    axB.set_xlim(0.44, 0.56); axB.set_xticks([0.45, 0.50, 0.55])
    axB.set_ylim(0, axB.get_ylim()[1] * 1.06)
    axB.set_ylabel("permutations", fontsize=fs.FS["label"])
    axB.tick_params(labelsize=fs.FS["tick"])
    axB2.axvline(yr["real_cv_auc"], color=fs.COL["alert"], lw=fs.LW["main"])
    axB2.set_xlim(0.90, 0.92); axB2.set_xticks([round(yr["real_cv_auc"], 3)]); axB2.set_ylim(axB.get_ylim())
    axB2.spines["left"].set_visible(False); axB2.set_yticks([])
    axB2.tick_params(labelsize=fs.FS["tick"])
    axB2.annotate("real\nlabels", xy=(yr["real_cv_auc"], axB.get_ylim()[1] * 0.55), xytext=(-5, 0), textcoords="offset points",
                  ha="right", va="center", fontsize=fs.FS["tick"], color=fs.COL["alert"])
    for ax_, x_ in ((axB, 1.0), (axB2, 0.0)):
        ax_.plot([x_ - 0.04, x_ + 0.04], [-0.02, 0.02], transform=ax_.transAxes, color=fs.COL["axis"], lw=fs.LW["axis"], clip_on=False)
    axB.set_xlabel("cross-validated ROC-AUC", fontsize=fs.FS["label"])
    axB.xaxis.set_label_coords(0.58, -0.135)
    fs.panel_label(axB, "B")

    # ---- C: early recognition --------------------------------------------
    axC = fig.add_subplot(gs[0, 2])
    _x = np.logspace(np.log10(1/len(y)), 0, 200)
    axC.plot(_x, _x, color="#999999", lw=fs.LW["ref"], ls=fs.REF_DASH)
    for label, col, neg, colour, lw in METHODS:
        if label not in ("Vina alone (baseline)", "IFP similarity alone"):
            continue
        s = (-d[col] if neg else d[col]).values.astype(float)
        fx, fy = enrichment_curve(s, y)
        key = "vina_only" if "baseline" in label else "ifp_only"
        axC.plot(fx, fy, color=colour, lw=fs.LW["main"],
                 label=label)
    axC.axvline(0.01, color=fs.COL["alert"], lw=fs.LW["thr"], ls=fs.THR_DASH)
    axC.set_xscale("log")
    axC.set_xlim(1/len(y), 1); axC.set_ylim(0, 1)
    axC.set_xticks([0.01, 0.05, 0.1, 1.0])
    axC.set_xticklabels(["1%", "5%", "10%", "100%"])
    axC.set_xlabel("fraction of ranked list screened (log)", fontsize=fs.FS["label"])
    axC.set_ylabel("fraction of actives recovered", fontsize=fs.FS["label"])
    axC.tick_params(labelsize=fs.FS["tick"])
    fs.panel_label(axC, "C")

    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    hA, lA = axA.get_legend_handles_labels()
    hA = hA + [Patch(facecolor=METHODS[0][3], alpha=0.13, linewidth=0), Line2D([], [], color="#999999", lw=fs.LW["ref"], ls=fs.REF_DASH)]
    lA = lA + ["95% bootstrap band, Vina alone", "chance"]
    fig.legend(hA, lA, loc="upper center", ncol=4, frameon=False, fontsize=fs.FS["legend"], bbox_to_anchor=(0.5, 1.0),
               handlelength=1.6, columnspacing=1.5, handletextpad=0.5, labelspacing=0.35)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    stem = out / "figure3_screening_performance_harmonised"
    fs.save_figure(fig, stem)

    prov = {
        # Figure 3 in the manuscript. The file stem says 4 for historical reasons
        # and the applicability-domain figure is the real Figure 4, so name both.
        "figure": "Figure 3 in the manuscript -- screening performance "
                  "(file stem figure4_screening_performance)",
        "generated": "2026-08-24",
        "inputs": {
            "ifp_per_compound": {"path": str(IFP.relative_to(ROOT)),
                                 "sha256": sha256(IFP), "n_rows": int(len(d)),
                                 "n_active": int(y.sum()),
                                 "n_inactive": int((y == 0).sum())},
            "ifp_summary": {"path": str(SUMMARY.relative_to(ROOT)),
                            "sha256": sha256(SUMMARY)},
            "y_randomisation": {"path": str(YRAND.relative_to(ROOT)),
                                "sha256": sha256(YRAND),
                                "n_permutations": yr["n_permutations"]},
        },
        "panelA_auc_by_convention": provenance,
        "tie_handling_note": (
            "Table 1 quotes RDKit Scoring.CalcAUC. Panel A is a ROC curve and "
            "is labelled with the tie-averaged area. The two agree to <0.001 "
            "except for max_ifp_tanimoto, which takes only "
            f"{int(d['max_ifp_tanimoto'].nunique())} distinct values across "
            f"{len(d)} compounds: 0.4487 (RDKit) vs 0.4422 (tie-averaged). "
            "Both are below chance."),
        "outputs": {p.name: sha256(p) for p in
                    (stem.with_suffix(".png"), stem.with_suffix(".svg"))},
    }
    (out / "figure4_screening_performance_provenance.json").write_text(
        json.dumps(prov, indent=2), encoding="utf-8")

    print(f"written -> {stem}.png / .svg")
    print(f"  panel A: {len(METHODS)} ROC curves, n={len(d)} "
          f"({int(y.sum())} active / {int((y==0).sum())} inactive)")
    print(f"  panel B: {yr['n_permutations']} permutations, "
          f"real {yr['real_cv_auc']:.3f} vs permuted {yr['permuted_mean']:.3f}")
    print(f"  panel C: enrichment, baseline vs IFP-alone")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

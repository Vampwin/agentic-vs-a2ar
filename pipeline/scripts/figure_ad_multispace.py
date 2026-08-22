"""
figure_ad_multispace.py -- the applicability-domain figure.

The claim this figure supports is narrow and should stay narrow: an
applicability-domain verdict is not a property of a library, it is a property
of a library *and a descriptor space*, and reporting one without the other is
uninterpretable. Nothing here says which space is correct.

Three panels, each answering a different objection:

  A  "Maybe the library is simply exotic."     -> Deliberate 2D-pharmacophore
     and 3D-shape selection toward 52 known A2AR ligands (236 -> 497
     compounds) moved physicochemical coverage 64.4% -> 86.3% and moved
     structural coverage not at all. The two are decoupled.

  B  "Maybe they are just below an arbitrary threshold."  -> The whole
     distribution is shown against the threshold, so how far outside is
     visible, not merely that it is outside.

  C  "Maybe one space is unlucky."  -> Per-compound count of spaces inside.

Usage:
    python figure_ad_multispace.py [--out DIR]
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch


def _root() -> Path:
    marker = "00_project_control/08_DECISION_LOG.md"
    for cand in (Path(r"D:\Agentic VS Workflow"), Path("/mnt/d/Agentic VS Workflow")):
        if (cand / marker).is_file():
            return cand
    raise SystemExit("project root not found")


ROOT = _root()
SWEEP = ROOT / "20_new_runs/ad_sweep_expanded_2026-08-21"

# Ordered so the four zero-coverage structural spaces lead and PhysChem --
# the space that disagrees with all of them -- sits last. The ordering is
# presentational; no space was preferred when the thresholds were set.
ORDER = ["ECFP4", "FCFP4", "TopTorsion", "Gobbi_Pharm2D",
         "AtomPair", "ErG", "MACCS", "PhysChem"]
LABEL = {"ECFP4": "ECFP4", "FCFP4": "FCFP4", "TopTorsion": "Topol.\ntorsion",
         "Gobbi_Pharm2D": "Gobbi 2D\npharmacophore", "AtomPair": "Atom\npair",
         "ErG": "ErG", "MACCS": "MACCS", "PhysChem": "Physico-\nchemical"}
KIND = {s: ("structural" if s in ("ECFP4", "FCFP4", "TopTorsion", "AtomPair", "MACCS")
            else "pharmacophoric" if s in ("Gobbi_Pharm2D", "ErG")
            else "physicochemical") for s in ORDER}
CKIND = {"structural": "#2c5f8a", "pharmacophoric": "#7a9a3c",
         "physicochemical": "#b5651d"}

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.2,
    "axes.linewidth": 0.7, "axes.edgecolor": "#333333",
    "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "savefig.dpi": 600, "figure.dpi": 150,
})


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "30_manuscript/figures"))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    cov = pd.read_csv(SWEEP / "ad_coverage_by_descriptor.csv").set_index("descriptor_space")
    per = pd.read_csv(SWEEP / "ad_per_compound.csv")
    comp = pd.DataFrame(
        json.loads((SWEEP / "comparison_vs_snplib.json").read_text(encoding="utf-8"))["comparison"]
    ).set_index("descriptor_space")

    # 7.48 in = 190 mm, Elsevier double-column: draw at final size so no
    # rescaling shrinks the type below the 7 pt floor.
    fig = plt.figure(figsize=(7.48, 6.1))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.88, 1.15],
                          width_ratios=[1.0, 0.52],
                          hspace=0.42, wspace=0.30,
                          left=0.085, right=0.985, top=0.925, bottom=0.085)

    # ---- Panel A : coverage before vs after similarity selection -----------
    axA = fig.add_subplot(gs[0, :])
    x = np.arange(len(ORDER))
    w = 0.38
    before = [comp.loc[s, "coverage_fraction_snplib"] * 100 for s in ORDER]
    after = [comp.loc[s, "coverage_fraction_expanded"] * 100 for s in ORDER]
    axA.bar(x - w / 2, before, w, color="#c9c9c9", edgecolor="#555555", linewidth=0.6)
    axA.bar(x + w / 2, after, w, color=[CKIND[KIND[s]] for s in ORDER],
            edgecolor="#333333", linewidth=0.6)
    for xi, (b, a) in enumerate(zip(before, after)):
        axA.text(xi - w / 2, b + 1.6, f"{b:.1f}", ha="center", va="bottom",
                 fontsize=6.8, color="#555555")
        axA.text(xi + w / 2, a + 1.6, f"{a:.1f}", ha="center", va="bottom",
                 fontsize=6.9, weight="bold" if a == 0 else "normal")
    axA.set_xticks(x)
    axA.set_xticklabels([LABEL[s] for s in ORDER], fontsize=7)
    axA.set_ylabel("compounds inside the\napplicability domain (%)")
    axA.set_ylim(0, 100)
    axA.spines[["top", "right"]].set_visible(False)
    handles = [Patch(facecolor="#c9c9c9", edgecolor="#555555",
                     label="before selection (n = 236)"),
               Patch(facecolor="#2c5f8a", edgecolor="#333333", label="structural"),
               Patch(facecolor="#7a9a3c", edgecolor="#333333", label="pharmacophoric"),
               Patch(facecolor="#b5651d", edgecolor="#333333", label="physicochemical")]
    axA.legend(handles=handles, ncol=4, fontsize=6.9, frameon=False,
               loc="upper left", bbox_to_anchor=(0.0, 1.19), columnspacing=1.4)
    axA.set_title("A   Similarity selection moved physicochemical coverage by 21.9 points\n"
                  "     and structural coverage by none",
                  loc="left", fontsize=8.4, weight="bold", pad=22)

    # ---- Panel B : full similarity distribution against each threshold -----
    axB = fig.add_subplot(gs[1, 0])
    for i, s in enumerate(ORDER):
        v = per[f"{s}_sim"].to_numpy()
        thr = float(cov.loc[s, "threshold"])
        parts = axB.violinplot([v], positions=[i], widths=0.72,
                               showextrema=False, showmedians=False)
        for b in parts["bodies"]:
            b.set_facecolor(CKIND[KIND[s]])
            b.set_alpha(0.34)
            b.set_edgecolor(CKIND[KIND[s]])
            b.set_linewidth(0.6)
        q1, med, q3 = np.percentile(v, [25, 50, 75])
        axB.vlines(i, q1, q3, color=CKIND[KIND[s]], linewidth=2.6)
        axB.plot(i, med, "o", ms=2.6, color="white", zorder=4,
                 markeredgecolor=CKIND[KIND[s]], markeredgewidth=0.6)
        axB.hlines(thr, i - 0.40, i + 0.40, color="#b00020", linewidth=1.3, zorder=5)
    axB.set_xticks(range(len(ORDER)))
    axB.set_xticklabels([LABEL[s].replace("\n", " ") for s in ORDER],
                        fontsize=6.8, rotation=38, ha="right")
    axB.set_ylabel("mean similarity to 5 nearest\nChEMBL A2AR training neighbours")
    axB.set_ylim(0, 1.02)
    axB.spines[["top", "right"]].set_visible(False)
    axB.plot([], [], color="#b00020", linewidth=1.3,
             label="domain threshold (10th percentile\nof training self-similarity)")
    axB.legend(fontsize=6.8, frameon=False, loc="lower left")
    axB.set_title("B   Distance from the domain, not merely in or out",
                  loc="left", fontsize=8.4, weight="bold", pad=6)

    # ---- Panel C : how many spaces call a compound in-domain ---------------
    axC = fig.add_subplot(gs[1, 1])
    counts = per["n_spaces_inside"].value_counts().reindex(range(9), fill_value=0)
    axC.bar(counts.index, counts.values, color="#4a4a4a", edgecolor="#222222",
            linewidth=0.6, width=0.72)
    for i, c in counts.items():
        if c:
            axC.text(i, c + 6, str(int(c)), ha="center", va="bottom", fontsize=6.8)
    axC.set_xlabel("descriptor spaces calling the\ncompound in-domain (of 8)", fontsize=7)
    axC.set_ylabel("compounds", fontsize=7.4)
    axC.set_xticks(range(9))
    axC.tick_params(labelsize=6.6)
    axC.set_ylim(0, max(counts.values) * 1.32)
    axC.spines[["top", "right"]].set_visible(False)
    n_max = int(per["n_spaces_inside"].max())

    # Which space does the in-domain verdict actually come from? Nearly always
    # the physicochemical one -- so a study reporting a single AD number is,
    # in practice, reporting that one space and calling it the domain.
    inside_cols = [f"{s}_inside" for s in ORDER]
    n_calls = int(per[inside_cols].sum().sum())
    n_phys = int(per["PhysChem_inside"].sum())
    axC.text(0.97, 0.80, f"{n_phys} of the {n_calls}\nin-domain calls come\n"
                         f"from the physico-\nchemical space",
             transform=axC.transAxes, ha="right", va="top", fontsize=6.8,
             color="#b5651d", linespacing=1.35)
    axC.set_title(f"C   No compound is inside\n     more than {n_max} of the 8 spaces",
                  loc="left", fontsize=8.4, weight="bold", pad=6)

    stem = out / "figure_ad_multispace_2026-08-21"
    fig.savefig(stem.with_suffix(".png"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)

    prov = {
        "figure": stem.name,
        "generated": "2026-08-21",
        "script": "figure_ad_multispace.py",
        "inputs": {
            "ad_coverage_by_descriptor.csv": sha256(SWEEP / "ad_coverage_by_descriptor.csv"),
            "ad_per_compound.csv": sha256(SWEEP / "ad_per_compound.csv"),
            "comparison_vs_snplib.json": sha256(SWEEP / "comparison_vs_snplib.json"),
        },
        "outputs": {p.name: sha256(p) for p in
                    (stem.with_suffix(".png"), stem.with_suffix(".svg"))},
        "panel_A_numbers": {
            s: {"before_pct": round(float(comp.loc[s, "coverage_fraction_snplib"]) * 100, 2),
                "after_pct": round(float(comp.loc[s, "coverage_fraction_expanded"]) * 100, 2)}
            for s in ORDER},
        "panel_B_thresholds": {s: float(cov.loc[s, "threshold"]) for s in ORDER},
        "panel_C_counts": {int(k): int(v) for k, v in counts.items()},
        "max_spaces_inside": n_max,
    }
    (SWEEP / "figure_ad_multispace_provenance.json").write_text(
        json.dumps(prov, indent=2), encoding="utf-8")

    print("panel C counts:", {int(k): int(v) for k, v in counts.items() if v})
    print("max spaces inside for any compound:", n_max)
    print(f"written -> {stem}.png / .svg")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

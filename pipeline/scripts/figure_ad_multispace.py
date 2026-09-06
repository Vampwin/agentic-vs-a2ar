"""
figure_ad_multispace.py -- the applicability-domain figure.

The claim this figure supports is narrow and should stay narrow: an
applicability-domain verdict is interpretable only when the *kind* of
descriptor is named. Structural exclusion and physicochemical admission are
both true of this library and are not in conflict; neither is about activity.
Nothing here says which space is correct.

Revised 2026-08-25 under ADR-AIIM-136. The earlier version argued from eight
admission rates spanning 0% to 86.3%, which reads as eight spaces disagreeing.
ADR-AIIM-126 showed that is wrong: seven of them agree closely once the
quantity the rule thresholds is shown directly. Panel B was rebuilt to show it.

Three panels, each answering a different objection:

  A  "Maybe the library is simply exotic."     -> A second library, anchored to
     the same biosynthetic superclasses and then selected toward 52 known A2AR
     ligands, raised physicochemical coverage 64.4% -> 86.3% and structural
     coverage not at all. These are two libraries, not one before and after a
     filter; if similarity selection could move a structural applicability
     domain, it would have moved it here.

  B  "Then the spaces contradict each other."  -> They do not. A Tanimoto of
     0.24 in ECFP4 and 0.91 in ErG are not the same quantity and cannot share
     an axis. A percentile inside each space's own training self-similarity
     distribution can, and on that axis seven of the eight spaces put the
     median compound between 0.23 and 3.82, all below one shared threshold
     line. The spread in admission rates is where each threshold happens to
     fall in a differently shaped distribution.

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
# Repointed 2026-08-25 (ADR-AIIM-136). The 2026-08-25 sweep reproduces the
# 2026-08-21 outputs byte-for-byte after normalising line endings (ADR-AIIM-126)
# and additionally carries the per-space training self-similarity distribution
# that panel B now needs. G6 and Supplementary S10 continue to read the older
# directory, which is correct: the contents are identical and re-pointing them
# would force a re-freeze for no change in any number.
SWEEP = ROOT / "20_new_runs/ad_sweep_selfsim_2026-08-25"

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
    pct = pd.read_csv(SWEEP / "ad_percentile_position.csv")

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
                     label="SNP-Lib (n = 236)"),
               Patch(facecolor="#2c5f8a", edgecolor="#333333",
                     label="LOTUS expansion, redocked (n = 497):  structural"),
               Patch(facecolor="#7a9a3c", edgecolor="#333333", label="pharmacophoric"),
               Patch(facecolor="#b5651d", edgecolor="#333333", label="physicochemical")]
    axA.legend(handles=handles, ncol=4, fontsize=6.6, frameon=False,
               loc="upper left", bbox_to_anchor=(0.0, 1.19), columnspacing=1.1)
    # Not "before and after". These are two libraries: the second was anchored
    # to the biosynthetic superclasses that dominate the first, then selected
    # toward 52 known A2AR ligands. If similarity selection could move a
    # structural applicability domain, it would have moved it here.
    axA.set_title("A", loc="left", fontsize=8.4, weight="bold", pad=8)

    # ---- Panel B : percentile position inside each space's own distribution -
    # A Tanimoto of 0.24 in ECFP4 and one of 0.91 in ErG are not the same
    # quantity and cannot share an axis, which is why the previous version of
    # this panel could not support the comparison the text was drawing from it.
    # A percentile can: it says where a compound sits inside the training set's
    # own self-similarity distribution, in that same space. The rule admits at
    # the 10th percentile everywhere, so one horizontal line serves all eight.
    axB = fig.add_subplot(gs[1, 0])
    # Log axis cannot show zero. Three compounds sit at percentile 0 exactly
    # (two in ECFP4, one physicochemical): less similar to the training set
    # than any training compound is to its own neighbours. The floor is set
    # below the smallest non-zero value the data can take, 1/5970 = 0.0168,
    # so only those three are moved, and the caption says so.
    FLOOR = 0.01
    med_pct = {}
    for i, s in enumerate(ORDER):
        v = np.clip(pct[f"{s}_pctile"].to_numpy(), FLOOR, 100.0)
        med_pct[s] = float(np.median(v))
        lv = np.log10(v)
        parts = axB.violinplot([lv], positions=[i], widths=0.72,
                               showextrema=False, showmedians=False)
        for b in parts["bodies"]:
            b.set_facecolor(CKIND[KIND[s]])
            b.set_alpha(0.34)
            b.set_edgecolor(CKIND[KIND[s]])
            b.set_linewidth(0.6)
        q1, m, q3 = np.percentile(lv, [25, 50, 75])
        axB.vlines(i, q1, q3, color=CKIND[KIND[s]], linewidth=2.6)
        axB.plot(i, m, "o", ms=2.6, color="white", zorder=4,
                 markeredgecolor=CKIND[KIND[s]], markeredgewidth=0.6)
    axB.axhline(np.log10(10.0), color="#b00020", linewidth=1.3, zorder=5)
    axB.set_xticks(range(len(ORDER)))
    axB.set_xticklabels([LABEL[s].replace("\n", " ") for s in ORDER],
                        fontsize=6.8, rotation=38, ha="right")
    ticks = [0.1, 0.3, 1, 3, 10, 30, 100]
    axB.set_yticks([np.log10(t) for t in ticks])
    axB.set_yticklabels([("%g" % t) for t in ticks], fontsize=6.8)
    axB.set_ylim(np.log10(FLOOR * 0.8), np.log10(140))
    axB.set_ylabel("percentile of the compound inside the\n"
                   "training set's own self-similarity (log)")
    axB.spines[["top", "right"]].set_visible(False)
    axB.plot([], [], color="#b00020", linewidth=1.3,
             label="the domain threshold, now one line\nfor all eight spaces (10th percentile)")
    axB.legend(fontsize=6.8, frameon=False, loc="upper left")
    axB.set_title("B", loc="left", fontsize=8.4, weight="bold", pad=6)

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
    axC.set_title("C", loc="left", fontsize=8.4, weight="bold", pad=6)

    stem = out / "figure_ad_multispace_2026-08-25"
    fig.savefig(stem.with_suffix(".png"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)

    prov = {
        "figure": stem.name,
        "generated": "2026-08-25",
        "script": "figure_ad_multispace.py",
        "supersedes": "figure_ad_multispace_2026-08-21",
        "adr": "ADR-AIIM-136",
        "sweep_dir": str(SWEEP.relative_to(ROOT).as_posix()),
        "inputs": {
            "ad_coverage_by_descriptor.csv": sha256(SWEEP / "ad_coverage_by_descriptor.csv"),
            "ad_per_compound.csv": sha256(SWEEP / "ad_per_compound.csv"),
            "ad_percentile_position.csv": sha256(SWEEP / "ad_percentile_position.csv"),
            "comparison_vs_snplib.json": sha256(SWEEP / "comparison_vs_snplib.json"),
        },
        "outputs": {p.name: sha256(p) for p in
                    (stem.with_suffix(".png"), stem.with_suffix(".svg"))},
        "panel_A_numbers": {
            s: {"snplib_pct": round(float(comp.loc[s, "coverage_fraction_snplib"]) * 100, 2),
                "expanded_pct": round(float(comp.loc[s, "coverage_fraction_expanded"]) * 100, 2)}
            for s in ORDER},
        "panel_B_median_percentile": {s: round(med_pct[s], 2) for s in ORDER},
        "panel_B_threshold_percentile": 10,
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

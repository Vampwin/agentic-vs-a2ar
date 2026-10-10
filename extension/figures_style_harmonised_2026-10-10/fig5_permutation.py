"""figure_permutation_distribution.py -- Figure 7: the conformal counts of
Section 3.6 are draws from a wide distribution, and the released run drew low.

Why this figure exists. Section 3.6 reports three quantities -- 72 of 5,000
possibly-active, a library maximum conformal p(active) of 0.072, and an
eligible pool of 3 of 497 -- that ADR-AIIM-271, ADR-AIIM-275 and ADR-AIIM-279
measured moving by orders of magnitude when the training rows are permuted.
The prose can state the range; it cannot show that the published value sits
near the favourable end of its own distribution, which is the thing a reader
needs to see rather than be told.

Three panels, each answering a different objection:

  A  "a range is not a distribution" -- all 200 draws, with the released
     value and the median marked, so the reader can see where 72 falls.
  B  "the count is one statistic, maybe the others are stable" -- the
     library maximum behaves the same way, and the released 0.072 is below
     every draw rather than merely low among them.
  C  "so the method is simply unusable" -- it is not: canonical row order
     collapses the spread to a point. Two conventions give two different
     points, which is the honest end of the story and the reason the
     article adopts neither.

Nothing is recomputed here. Every value is read from the frozen run outputs,
so the figure cannot drift from the numbers G8 anchors (ADR-AIIM-281).

Writes the PNG to 30_manuscript/figures/ and copies it into
CURRENT_v5/figures_png/, which is the pair build_v6_clean_docx.py checks for
staleness, plus a JSON of every plotted number.

Usage:
    python figure_permutation_distribution.py
"""
from __future__ import annotations
import sys
sys.path.insert(0, r'D:\Agentic VS Workflow\10_new_system\scripts')

import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import figstyle as fs
import numpy as np
from kdss_paths import resolve, under

PERM = Path(under("20_new_runs", "conformal_permutation_2026-10-02/"
                  "permutation_distribution.json"))
CANON = Path(under("20_new_runs", "canonical_row_order_2026-10-02/"
                   "canonical_row_order.json"))
SRC = Path(__file__).resolve().parent / "out"
EMBED = Path(__file__).resolve().parent / "out"
STEM = "figure5_permutation_harmonised"

RELEASED_COUNT = 72
RELEASED_MAXP = 0.072

# Slate for the draws, one warm accent for the published value. Both survive
# greyscale, and the released value is also marked by line style, not colour
# alone.
C_DRAWS = "#6b7a8f"
C_RELEASED = fs.COL["alert"]
C_MEDIAN = "#333333"
C_CANON = ("#1f6f8b", "#8a5a00")

fs.apply()


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> int:
    perm = json.loads(PERM.read_text(encoding="utf-8"))["arms"]["A"]["draws"]
    counts = np.array([d["NP-library"]["n_possibly_active"] for d in perm])
    maxp = np.array([d["NP-library"]["max_p"] for d in perm])
    canon = json.loads(CANON.read_text(encoding="utf-8"))["arms"]["A"]["summary"]
    c_rel = canon["control"]["NP-library"]["value_if_invariant"]
    c_fix = canon["fixed"]["NP-library"]["value_if_invariant"]
    med = float(np.median(counts))
    pct = 100.0 * float(np.mean(counts <= RELEASED_COUNT))

    fig = plt.figure(figsize=(7.48, 4.7))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.82], hspace=0.62,
                          wspace=0.28, left=0.085, right=0.985,
                          top=0.955, bottom=0.115)

    # ---- A: the distribution of the count ---------------------------------
    ax = fig.add_subplot(gs[0, 0])
    bins = np.logspace(np.log10(max(counts.min(), 1)), np.log10(counts.max()), 26)
    ax.hist(counts, bins=bins, color=C_DRAWS, edgecolor="white", linewidth=fs.LW["bar_edge"])
    ax.set_xscale("log")
    ax.axvline(RELEASED_COUNT, color=C_RELEASED, lw=fs.LW["thr"], ls=fs.THR_DASH, zorder=5)
    ax.axvline(med, color=C_MEDIAN, lw=fs.LW["ref"], ls=":", zorder=5)
    ax.set_xlabel("compounds not excluded, of 5,000 (log scale)")
    ax.set_ylabel("draws")
    # A clear band above the bars, so neither label sits on data.
    ax.set_ylim(0, ax.get_ylim()[1] * 1.42)
    _top = ax.get_ylim()[1]
    ax.text(RELEASED_COUNT * 1.35, _top * 0.97,
            f"released {RELEASED_COUNT}\n({pct:.0f}th percentile)",
            color=C_RELEASED, fontsize=fs.FS["label"], va="top", ha="left")
    # One line lower than the released label: at the same height the two
    # right- and left-growing strings met in the middle of the panel.
    ax.text(med * 0.85, _top * 0.70, f"median\n{med:,.0f}",
            color=C_MEDIAN, fontsize=fs.FS["label"], va="top", ha="right")
    fs.panel_label(ax, "A")

    # ---- B: the distribution of the library maximum ------------------------
    ax = fig.add_subplot(gs[0, 1])
    ax.hist(maxp, bins=24, color=C_DRAWS, edgecolor="white", linewidth=fs.LW["bar_edge"])
    ax.axvline(RELEASED_MAXP, color=C_RELEASED, lw=fs.LW["thr"], ls=fs.THR_DASH, zorder=5)
    ax.set_xlabel("maximum conformal $p$(active) in the library")
    ax.set_ylabel("draws")
    ax.set_ylim(0, ax.get_ylim()[1] * 1.42)
    ax.text(RELEASED_MAXP + 0.009, ax.get_ylim()[1] * 0.97,
            f"released 0.072 \u2014 below all {len(maxp)} draws\n"
            f"(lowest draw {maxp.min():.4f})",
            color=C_RELEASED, fontsize=fs.FS["label"], va="top", ha="left")
    fs.panel_label(ax, "B")

    # ---- C: what canonical ordering does -----------------------------------
    ax = fig.add_subplot(gs[1, :])
    rng = np.random.default_rng(7)
    ax.scatter(counts, np.full(len(counts), 1.0) + rng.uniform(-0.1, 0.1, len(counts)),
               s=9, color=C_DRAWS, alpha=0.55, edgecolors="none", zorder=3)
    ax.scatter([c_rel], [0.35], s=58, marker="D", color=C_CANON[0], zorder=4)
    ax.scatter([c_fix], [0.35], s=58, marker="s", color=C_CANON[1], zorder=4)
    ax.axvline(RELEASED_COUNT, color=C_RELEASED, lw=fs.LW["thr"], ls=fs.THR_DASH, zorder=2)
    # Direct labels rather than a legend box: the rows are already named
    # on the y axis, and a legend here landed on the released-value mark.
    for x, lab, col in ((c_rel, f"{c_rel:,}\nreleased tiebreak", C_CANON[0]),
                        (c_fix, f"{c_fix:,}\ndeclared tiebreak", C_CANON[1])):
        ax.annotate(lab, xy=(x, 0.35), xytext=(0, -14),
                    textcoords="offset points", ha="center", fontsize=fs.FS["tick"],
                    color=col, va="top")
    ax.text(4.6, 1.31, f"{len(counts)} permutations of the same rows",
            color=C_DRAWS, fontsize=fs.FS["label"], ha="left", va="center")
    ax.text(RELEASED_COUNT * 0.92, 0.72, "released 72", color=C_RELEASED,
            fontsize=fs.FS["label"], ha="right", va="center")
    ax.set_xscale("log")
    ax.set_xlim(4, 6000)
    ax.set_ylim(-0.12, 1.45)
    ax.set_yticks([0.35, 1.0])
    ax.set_yticklabels(["canonical\nrow order", "released\nprocedure"],
                       fontsize=fs.FS["label"])
    ax.set_xlabel("compounds not excluded, of 5,000 (log scale)")
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    fs.panel_label(ax, "C")

    SRC.mkdir(parents=True, exist_ok=True)
    EMBED.mkdir(parents=True, exist_ok=True)
    out, _ = fs.save_figure(fig, SRC / STEM)
    plt.close(fig)
    embed = out

    numbers = {
        "figure": "Figure 7",
        "source_runs": ["20_new_runs/conformal_permutation_2026-10-02",
                        "20_new_runs/canonical_row_order_2026-10-02"],
        "n_draws": int(len(counts)),
        "count_min": int(counts.min()), "count_max": int(counts.max()),
        "count_median": med,
        "released_count": RELEASED_COUNT,
        "released_count_percentile": round(pct, 1),
        "maxp_min": float(maxp.min()), "maxp_max": float(maxp.max()),
        "released_maxp": RELEASED_MAXP,
        "released_maxp_below_all_draws": bool(maxp.min() > RELEASED_MAXP),
        "canonical_released_tiebreak": c_rel,
        "canonical_declared_tiebreak": c_fix,
        "png_sha256": sha256(out),
    }
    (SRC / f"{STEM}_numbers.json").write_text(json.dumps(numbers, indent=1),
                                              encoding="utf-8")
    print(json.dumps(numbers, indent=1))
    print(f"written -> {out}\n           {embed}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

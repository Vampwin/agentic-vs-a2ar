"""figure_m6_panel_canonical.py -- Figure 6: the completed gate replaces the
panel rather than shortening it.

Why this figure exists. ADR-AIIM-264 concluded from an eligible pool of three
that enforcing scaffold uniqueness "does not reorder this panel; it shortens
it", and Section 3.6 reports three of 497 as the binding constraint on the
screening arm. ADR-AIIM-279 re-ran that selection with one column recomputed
under canonical row order and found the pool at 9 to 298, the panel at 6 to
67, and -- in both configurations trained on the recovered inactives -- none
of the three released compounds even eligible.

A count can be written in a sentence. The disappearance of three named
compounds cannot: the reader has to see the same three rows go from filled
to absent as the training data improves, which is panel B.

Panel A carries the magnitudes with the released values as reference lines.
Panel B carries the identities, because identity is the finding.

Nothing is recomputed here; every value is read from the frozen run output,
so the figure cannot drift from the numbers G8 anchors.

Usage:
    python figure_m6_panel_canonical.py
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from kdss_paths import resolve, under

M6C = Path(under("20_new_runs", "m6_canonical_2026-10-02/m6_canonical.json"))
SRC = resolve("30_manuscript/figures")
EMBED = resolve("30_manuscript/CURRENT_v5/figures_png")
STEM = "figure8_m6_panel_canonical_2026-10-02"

REL_ELIGIBLE, REL_PANEL = 3, 2
C_ELIG = "#9bb0c4"
C_PANEL = "#42566b"
C_RELEASED = "#b03a2e"
C_ABSENT = "#c8ccd0"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.2,
    "axes.linewidth": 0.7, "axes.edgecolor": "#333333",
    "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "axes.spines.top": False, "axes.spines.right": False,
    "savefig.dpi": 600, "figure.dpi": 150, "savefig.facecolor": "white",
})

ORDER = [("A_released_split", "A\nreleased tiebreak"),
         ("A_fixed_split", "A\ndeclared tiebreak"),
         ("B_released_split", "B\nreleased tiebreak"),
         ("B_fixed_split", "B\ndeclared tiebreak")]


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> int:
    d = json.loads(M6C.read_text(encoding="utf-8"))
    cfg = d["configurations"]
    rel3 = d["released_panel"]
    elig = [cfg[k]["n_eligible"] for k, _ in ORDER]
    panel = [cfg[k]["panel_size"] for k, _ in ORDER]
    entered = [cfg[k]["n_entered"] for k, _ in ORDER]
    labels = [lab for _, lab in ORDER]
    x = np.arange(len(ORDER))

    fig = plt.figure(figsize=(7.48, 4.45))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 0.86], hspace=0.30,
                          left=0.145, right=0.985, top=0.955, bottom=0.135)

    # ---- A: how many -------------------------------------------------------
    ax = fig.add_subplot(gs[0])
    w = 0.36
    ax.bar(x - w / 2, elig, w, color=C_ELIG, edgecolor="white", linewidth=0.5,
           label="eligible for Tier A, of 497")
    ax.bar(x + w / 2, panel, w, color=C_PANEL, edgecolor="white", linewidth=0.5,
           label="panel after scaffold uniqueness")
    for xi, (e, p) in enumerate(zip(elig, panel)):
        ax.text(xi - w / 2, e * 1.09, f"{e}", ha="center", fontsize=7.4,
                color=C_PANEL)
        ax.text(xi + w / 2, p * 1.09, f"{p}", ha="center", fontsize=7.4,
                color=C_PANEL)
    ax.axhline(REL_ELIGIBLE, color=C_RELEASED, lw=1.2, ls="--", zorder=1)
    ax.axhline(REL_PANEL, color=C_RELEASED, lw=1.2, ls=":", zorder=1)
    ax.set_yscale("log")
    # Headroom for the two-column legend: at 900 the 298 bar's value label
    # met the legend's second row.
    ax.set_ylim(1, 6000)
    ax.set_xticks(x)
    ax.set_xticklabels([])
    ax.set_ylabel("compounds (log scale)")
    from matplotlib.lines import Line2D
    handles, labs = ax.get_legend_handles_labels()
    handles += [Line2D([], [], color=C_RELEASED, lw=1.2, ls="--"),
                Line2D([], [], color=C_RELEASED, lw=1.2, ls=":")]
    labs += ["released run: 3 eligible", "released run: panel of 2"]
    ax.legend(handles, labs, loc="upper left", frameon=False,
              fontsize=7.4, ncol=2, columnspacing=1.6,
              handlelength=1.4, handletextpad=0.5, borderaxespad=0.2)
    ax.text(-0.115, 1.06, "(A)", transform=ax.transAxes, fontweight="bold",
            fontsize=9.2, va="top")

    # ---- B: which compounds ------------------------------------------------
    ax = fig.add_subplot(gs[1])
    rows = list(reversed(rel3))          # top row = first compound
    for yi, cid in enumerate(rows):
        for xi, (k, _) in enumerate(ORDER):
            c = cfg[k]
            eligible = c["released_three_still_eligible"][cid]
            in_panel = cid in c["panel"]
            if in_panel:
                ax.scatter(xi, yi, s=92, marker="o", color=C_PANEL, zorder=3)
            elif eligible:
                ax.scatter(xi, yi, s=92, marker="o", facecolors="none",
                           edgecolors=C_PANEL, linewidths=1.2, zorder=3)
            else:
                ax.scatter(xi, yi, s=58, marker="x", color=C_ABSENT,
                           linewidths=1.6, zorder=3)
    for xi, n in enumerate(entered):
        ax.text(xi, -1.0, f"+{n}", ha="center", va="center", fontsize=8.0,
                color=C_PANEL, fontweight="bold")
    ax.set_yticks(list(range(len(rows))) + [-1.0])
    ax.set_yticklabels([c.replace("NPLIB", "NPLIB ") for c in rows]
                       + ["new in panel"], fontsize=7.6)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7.6)
    ax.set_xlim(-0.55, len(ORDER) - 0.45)
    ax.set_ylim(-1.6, len(rows) + 0.55)
    ax.axhline(-0.5, color="#dddddd", lw=0.8)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(length=0)
    ax.set_xlabel("training set and tiebreak convention", labelpad=6)
    # A compact key, placed in the empty upper right of the grid.
    _ky = len(rows) + 0.18
    ax.scatter(-0.42, _ky, s=72, marker="o", color=C_PANEL, zorder=3)
    ax.text(-0.33, _ky, "selected into the panel", fontsize=7.2,
            va="center")
    ax.scatter(1.05, _ky, s=72, marker="o", facecolors="none",
               edgecolors=C_PANEL, linewidths=1.2, zorder=3)
    ax.text(1.14, _ky, "eligible, displaced by a better-scoring\n"
            "compound of the same scaffold", fontsize=7.2, va="center")
    ax.scatter(2.92, _ky, s=48, marker="x", color=C_ABSENT,
               linewidths=1.6, zorder=3)
    ax.text(3.01, _ky, "not eligible", fontsize=7.2, va="center")
    ax.text(-0.115, 1.06, "(B)", transform=ax.transAxes, fontweight="bold",
            fontsize=9.2, va="top")

    SRC.mkdir(parents=True, exist_ok=True)
    EMBED.mkdir(parents=True, exist_ok=True)
    out = SRC / f"{STEM}.png"
    fig.savefig(out)
    plt.close(fig)
    (EMBED / f"{STEM}.png").write_bytes(out.read_bytes())

    numbers = {
        "figure": "Figure 6",
        "source_run": "20_new_runs/m6_canonical_2026-10-02",
        "released_eligible": REL_ELIGIBLE, "released_panel": REL_PANEL,
        "configurations": [k for k, _ in ORDER],
        "eligible": elig, "panel": panel, "entered": entered,
        "released_three": rel3,
        "released_three_eligible": {
            k: cfg[k]["released_three_still_eligible"] for k, _ in ORDER},
        "released_three_in_panel": {
            k: [c for c in rel3 if c in cfg[k]["panel"]] for k, _ in ORDER},
        "png_sha256": sha256(out),
    }
    (SRC / f"{STEM}_numbers.json").write_text(json.dumps(numbers, indent=1),
                                              encoding="utf-8")
    print(json.dumps(numbers, indent=1))
    print(f"written -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

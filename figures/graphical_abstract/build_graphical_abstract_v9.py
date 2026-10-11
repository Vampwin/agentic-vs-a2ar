"""Graphical abstract, v9 (ADR-AIIM-401): the v8 funnel with three changes asked for by the internal review.

  1. Proportion. The canvas is 590 x 236 pt, 2.5 to 1, the proportion the publisher asks for (v8: 511 x 232, 2.2 to 1).
  2. What the numbers belong to. The left panel is labelled as the frozen screening run; the third layer is labelled as a
     re-analysis of the ligand-based step and carries its own denominator (of 5,000), so that "3 of 497" and "72 of 5,000"
     cannot be read as successive stages of one funnel.
  3. The seeded-fault result under both contracts. Layer 2 gives 124 / 128 for the benchmark's harness contract and 101 / 128
     for the contracts as deployed (Section 3.1; Supplementary Table S6).
Secondary text is reduced so that the smallest type is larger than in v8.

Every number is a literal copied from the manuscript text (Sections 2.7, 3.1, 3.5, 3.6 and the abstract); nothing here
recomputes them. Pictograms, palette, fonts and stroke weights are those of v7/v8 and of Figures 1 and 2. Wording is neutral:
no compound is described as a candidate or as prioritised.

    python build_graphical_abstract_v9.py
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from PIL import Image

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from build_graphical_abstract_v7 import (  # noqa: E402
    AMBER, AMBER_F, BLUE, BLUE_F, FILL, INK, MOSS, MOSS_F, RULE, SOFT, VIOLET, VIOLET_F, SCRIPTS,
    Svg, bug, contract_page, fits, path, shuffle_glyph)

GW, GH = 590.0, 236.0
STEM = "figure_graphical_abstract_v9_2026-10-11"

ROWS = [("25,604", "biosynthetic-superclass matches"),
        ("21,882", "after standardisation and filters"),
        ("8,000", "2D pharmacophore screen"),
        ("5,000", "3D refinement"),
        ("500", "best reduced-docking scores"),
        ("497", "completed full-protocol docking"),
        ("482", "with an on-pocket pose"),
        ("3", "released by the original rule")]

# colour, fill, title, glyph, (big, label) x 2, note
LAYERS = [
    (MOSS, MOSS_F, "2 Contract coverage: 128 seeded faults, 11 of 18 fault families", "contract",
     ("124 / 128", "detected, harness contract"), ("101 / 128", "detected, contracts as deployed"),
     "snapshot comparator: 128 / 128"),
    (AMBER, AMBER_F, "3 Beyond the written contracts", "bug",
     ("1 / 5", "held-out faults detected"), ("25 / 25", "scores changed by equivalent SMILES"),
     "no contract encoded this relation"),
    (VIOLET, VIOLET_F, "4 Decision stability: re-analysis of the ligand-based step", "shuffle",
     ("72", "of 5,000 eligible, released run"), ("6 to 4,119", "of 5,000 over 200 row orders"),
     "median 1,004; same molecules and model"),
]


def build(out: Path) -> Path:
    s = Svg(GW, GH, "Auditing agent-generated biomedical code with executable contracts",
            "Left, the frozen screening run, an agent-written A2A receptor virtual screen drawn as a funnel: 25,604 compounds, "
            "21,882 after filters, 8,000 after a 2D pharmacophore screen, 5,000 after 3D refinement, 500 selected for full "
            "docking, 497 completed, 482 on pocket, 3 released by the original rule. Right, three layers of evaluation of that "
            "workflow. Contract coverage: of 128 seeded faults in 11 of 18 fault families, 124 were detected under the "
            "benchmark's harness contract and 101 under the contracts as deployed; a snapshot comparator detected all 128. "
            "Beyond the written contracts: 1 of 5 held-out faults was detected and equivalent SMILES changed 25 of 25 docking "
            "scores. Decision stability, a re-analysis of the ligand-based step: 72 of 5,000 compounds were eligible in the "
            "released run against 6 to 4,119 over 200 orders of the training rows.")
    s.text(GW / 2, 14.5, fits("Auditing agent-generated biomedical code with executable contracts", 11.5, GW - 12, bold=True),
           size=11.5, weight="bold")

    top, bot = 23.0, 206.0
    # ------------------------------------------------------------ left: the frozen screening run
    px, pw = 6.0, 216.0
    s.box(px, top, pw, bot - top, fill="#FFFFFF", stroke=BLUE, sw=0.9, r=3)
    s.box(px, top, pw, 14.5, fill=BLUE_F, stroke=BLUE, sw=0.9, r=3)
    s.text(px + 8, top + 9.9, fits("1 Frozen run: an agent-written A2A receptor screen", 8.2, pw - 12, bold=True),
           size=8.2, weight="bold", fill=BLUE, anchor="start")
    rx, rw, rh, rgap, ry0 = px + 6, pw - 20.0, 15.2, 4.4, 42.0
    ys = [ry0 + i * (rh + rgap) for i in range(len(ROWS))]
    for i, (n, d) in enumerate(ROWS):
        last = i == len(ROWS) - 1
        s.box(rx, ys[i], rw, rh, fill=BLUE_F if last else FILL, stroke=BLUE if last else RULE, sw=0.9 if last else 0.6, r=2.5)
        s.text(rx + 42, ys[i] + 10.8, n, size=8.8 if not last else 9.8, weight="bold", fill=BLUE if last else INK, anchor="end")
        s.text(rx + 47, ys[i] + 10.6, fits(d, 7.8, rw - 50, bold=last), size=7.8, fill=INK if last else SOFT,
               anchor="start", weight="bold" if last else "normal")
        if i < len(ROWS) - 1:
            s.arrow(rx + 22, ys[i] + rh + 0.4, rx + 22, ys[i + 1] - 0.4, sw=0.6)
    s.text(px + pw / 2, bot - 4.2, fits("the released decision; any follow-up needs experiment", 7.6, pw - 10), size=7.6,
           fill=SOFT, style="italic")

    # ------------------------------------------------------------ right: three layers of evaluation
    lx, gap = 236.0, 6.5
    lw = GW - 6.0 - lx
    lh = (bot - top - 2 * gap) / 3
    for k, (col, fil, title, glyph, a, b, note) in enumerate(LAYERS):
        y = top + k * (lh + gap)
        s.box(lx, y, lw, lh, fill="#FFFFFF", stroke=col, sw=0.9, r=3)
        s.box(lx, y, lw, 14.5, fill=fil, stroke=col, sw=0.9, r=3)
        s.text(lx + 8, y + 9.9, fits(title, 8.2, lw - 12, bold=True), size=8.2, weight="bold", fill=col, anchor="start")
        cy = y + 14.5 + (lh - 14.5) / 2
        gx = lx + 20
        if glyph == "contract":
            contract_page(s, gx, cy, col)
        elif glyph == "bug":
            bug(s, gx, cy - 1, col, fil)
        else:
            shuffle_glyph(s, gx, cy, col, fil)
        colw = (lw - 44) / 2
        for j, (big, label) in enumerate((a, b)):
            ox = lx + 40 + colw * (j + 0.5)
            s.text(ox, cy - 5.5, fits(big, 15, colw - 4, bold=True), size=15, weight="bold", fill=col)
            s.text(ox, cy + 5.2, fits(label, 7.8, colw - 2, bold=True), size=7.8, weight="bold")
        s.text(lx + 40 + colw, cy + 15.6, fits(note, 7.6, lw - 48), size=7.6, fill=SOFT)
        yc = y + lh / 2
        if k < 2:     # layers 2 and 3 examine the workflow as a whole
            s.arrow(lx - 0.5, yc, px + pw + 1.2, yc, stroke=col, sw=0.9)
        else:         # layer 4 re-ran the ligand-based step applied to the 5,000 compounds
            bx = px + pw - 6.0
            y0, y1 = ys[3] + 1.0, ys[7] + rh - 1.0
            path(s, f"M{bx - 3:g},{y0:g} L{bx:g},{y0:g} L{bx:g},{y1:g} L{bx - 3:g},{y1:g}", stroke=col, sw=1.1, cap="butt",
                 join="miter")
            s.arrow(lx - 0.5, yc, bx + 1.0, yc, stroke=col, sw=0.9)

    s.box(6, 213, GW - 12, 18.5, fill=FILL, stroke=RULE, sw=0.7)
    s.text(GW / 2, 225.2, fits("Contract coverage, decision stability and biomedical validation must be evaluated separately.",
                               8.6, GW - 28, bold=True), size=8.6, weight="bold")
    p = out / f"{STEM}.svg"
    s.save(p)
    return p


def main() -> int:
    svg = build(HERE)
    png = svg.with_suffix(".png")
    r = subprocess.run([sys.executable, str(SCRIPTS / "render_svg_to_png.py"), str(svg), str(png)])
    if r.returncode:
        return r.returncode
    with Image.open(png) as im:
        im.convert("RGB").save(png, dpi=(600, 600), optimize=True)
        print(png.name, im.size, "ratio %.3f" % (im.size[0] / im.size[1]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

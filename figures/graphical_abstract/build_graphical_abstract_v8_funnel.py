"""Graphical abstract, v8 (ADR-AIIM-371): the audited workflow as a funnel, the three checks as layers of control on it.

Why it was redrawn. The v7 graphic read as four sequential stages, 1 -> 2 -> 3 -> 4. That is not the structure of the study.
Only stage 1 is a pipeline that compounds pass through; stages 2 to 4 are three evaluations of the *control layer* that audits
that pipeline (Figure 1 draws the same thing: nine gates and, beneath them, the append-only registry). So here:

  left   the screen the agent wrote, as a vertical funnel with the number of compounds that passed each step, down to the three
         the workflow released under its original dual-criterion rule (Sections 2.7 and 3.6). The released decision is the
         object that the checks examine; the manuscript puts none forward as a candidate, which the line under the funnel leaves to experiment.
  right  three layers of control, each with the one result that answers its own question:
           2  seeded faults   124 / 128 detected (a golden-file snapshot comparator detected 128 / 128)
           3  blind faults    1 / 5 detected; the SMILES-invariance relation, which no contract encoded, identified a defect (25 of 25 scores changed)
           4  reanalysis      eligible compounds of 5,000: 72 in the released run, median 1,004 over 200 row orders (range 6 to 4,119)
         Layers 2 and 3 point at the workflow as a whole; layer 4 points at the part of the funnel it re-ran (5,000 -> 3).

Every number is a literal copied from the manuscript text (Sections 2.7, 3.1, 3.2, 3.5, 3.6, abstract); nothing here recomputes
them. Pictograms, palette, fonts and stroke weights are those of v7 and of Figures 1 and 2. The screen is the original LOTUS
campaign; the later SNP-Lib / DiffDock / MD work (Sections 3.8 onward) is not drawn.

    python build_graphical_abstract_v8_funnel.py
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from PIL import Image

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from build_graphical_abstract_v7 import (  # noqa: E402
    AMBER, AMBER_F, BLUE, BLUE_F, FILL, INK, LINE, MOSS, MOSS_F, RULE, SOFT, VIOLET, VIOLET_F, SCRIPTS,
    Svg, bug, contract_page, fits, path, shuffle_glyph)

GW, GH = 511.0, 232.0
STEM = "figure_graphical_abstract_v8_funnel_2026-10-10"

ROWS = [("25,604", "biosynthetic-superclass matches"),
        ("21,882", "after standardisation and filters"),
        ("8,000", "2D pharmacophore screen"),
        ("5,000", "3D refinement"),
        ("500", "best reduced-docking scores"),
        ("497", "completed full-protocol docking"),
        ("482", "with an on-pocket pose"),
        ("3", "released by the original rule")]

LAYERS = [
    (MOSS, MOSS_F, "2 Specify and test contracts",
     ["10 contracts, 94 checks,", "gates G0 to G8;", "128 seeded faults"],
     "124 / 128", "seeded faults detected", "snapshots: 128 / 128", "contract"),
    (AMBER, AMBER_F, "3 Challenge the contracts",
     ["5 blind-authored faults;", "4 metamorphic relations"],
     "1 / 5", "blind faults detected", "SMILES test adds coverage", "bug"),
    (VIOLET, VIOLET_F, "4 Reanalyse the decision",
     ["discrimination,", "applicability domain,", "stability to row order"],
     "72 vs 1,004", "eligible compounds of 5,000", "released vs 200-order median", "shuffle"),
]


def build(out: Path) -> Path:
    s = Svg(GW, GH, "Auditing agent-generated biomedical code with executable contracts",
            "Left, the agent-written A2A receptor virtual screen as a funnel: 25,604 candidates, 21,882 after filters, 8,000 after "
            "a 2D pharmacophore screen, 5,000 after 3D refinement, 500 selected for full docking, 497 completed, 482 on pocket, 3 "
            "released by the original rule. Right, three layers of control applied to that workflow: ten contracts detected 124 of "
            "128 seeded faults, where deterministic snapshots detected all 128; 1 of 5 blind-authored faults was detected and the "
            "SMILES-invariance relation identified a defect; and the number of eligible compounds of 5,000 was 72 in the released run against a "
            "median of 1,004 over 200 reorderings of the training rows.")
    s.text(GW / 2, 14, fits("Auditing agent-generated biomedical code with executable contracts", 11, GW - 12, bold=True),
           size=11, weight="bold")

    top, bot = 22.0, 202.5
    # ------------------------------------------------------------ left: the workflow
    px, pw = 6.0, 230.0
    s.box(px, top, pw, bot - top, fill="#FFFFFF", stroke=BLUE, sw=0.9, r=3)
    s.box(px, top, pw, 14, fill=BLUE_F, stroke=BLUE, sw=0.9, r=3)
    s.text(px + 8, top + 9.2, fits("1 The workflow: an agent-written A2A receptor screen", 7.8, pw - 12, bold=True),
           size=7.8, weight="bold", fill=BLUE, anchor="start")
    rx, rw, rh, rgap, ry0 = px + 6, 212.0, 15.0, 4.5, 40.0
    ys = [ry0 + i * (rh + rgap) for i in range(len(ROWS))]
    for i, (n, d) in enumerate(ROWS):
        last = i == len(ROWS) - 1
        s.box(rx, ys[i], rw, rh, fill=BLUE_F if last else FILL, stroke=BLUE if last else RULE, sw=0.9 if last else 0.6, r=2.5)
        s.text(rx + 44, ys[i] + 10.6, n, size=8.6 if not last else 9.6, weight="bold", fill=BLUE if last else INK, anchor="end")
        s.text(rx + 50, ys[i] + 10.4, fits(d, 7.4, rw - 56, bold=last), size=7.4, fill=INK if last else SOFT,
               anchor="start", weight="bold" if last else "normal")
        if i < len(ROWS) - 1:
            s.arrow(rx + 22, ys[i] + rh + 0.4, rx + 22, ys[i + 1] - 0.4, sw=0.6)
    s.text(px + pw / 2, bot - 4.0, fits("the decision examined by layers 2 to 4; follow-up needs experiment", 7.2, pw - 10), size=7.2, fill=SOFT,
           style="italic")

    # ------------------------------------------------------------ right: three layers of control
    lx, lw, gap = 250.0, 255.0, 7.0
    lh = (bot - top - 2 * gap) / 3
    for k, (col, fil, title, lines, big, label, note, glyph) in enumerate(LAYERS):
        y = top + k * (lh + gap)
        s.box(lx, y, lw, lh, fill="#FFFFFF", stroke=col, sw=0.9, r=3)
        s.box(lx, y, lw, 13, fill=fil, stroke=col, sw=0.9, r=3)
        s.text(lx + 8, y + 9.4, fits(title, 7.8, lw - 12, bold=True), size=7.8, weight="bold", fill=col, anchor="start")
        cy = y + 13 + (lh - 13) / 2
        if glyph == "contract":
            contract_page(s, lx + 21, cy, col)
        elif glyph == "bug":
            bug(s, lx + 21, cy - 1, col, fil)
        else:
            shuffle_glyph(s, lx + 21, cy, col, fil)
        ty = cy - (len(lines) - 1) * 4.7 + 2.4
        for j, t in enumerate(lines):
            s.text(lx + 42, ty + j * 9.4, fits(t, 7.2, 108), size=7.2, fill=SOFT, anchor="start")
        ox = lx + 203
        s.text(ox, cy - 4.5, fits(big, 15, 98, bold=True), size=15, weight="bold", fill=col)
        s.text(ox, cy + 6.2, fits(label, 7.2, 98, bold=True), size=7.2, weight="bold")
        s.text(ox, cy + 15.2, fits(note, 7.0, 98), size=7.0, fill=SOFT)
        yc = y + lh / 2
        if k < 2:     # layers 2 and 3 audit the workflow as a whole
            s.arrow(lx - 0.5, yc, px + pw + 1.2, yc, stroke=col, sw=0.9)
        else:         # layer 4 re-ran the part of the funnel from 5,000 compounds to the released three
            bx = px + pw - 6.0
            y0, y1 = ys[3] + 1.0, ys[7] + rh - 1.0
            path(s, f"M{bx - 3:g},{y0:g} L{bx:g},{y0:g} L{bx:g},{y1:g} L{bx - 3:g},{y1:g}", stroke=col, sw=1.1, cap="butt",
                 join="miter")
            s.arrow(lx - 0.5, yc, bx + 1.0, yc, stroke=col, sw=0.9)

    s.box(6, 210, GW - 12, 18, fill=FILL, stroke=RULE, sw=0.7)
    s.text(GW / 2, 221.8, fits("Contract coverage, decision stability and biomedical validation must be evaluated separately.",
                               8, GW - 28, bold=True), size=8, weight="bold")
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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

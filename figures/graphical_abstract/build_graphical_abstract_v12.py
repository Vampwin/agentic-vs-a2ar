"""Graphical abstract, v12 (ADR-AIIM-454): v11 with card 3 showing both rounds of errors written blind to the project (1 of 5, then 3 of 6 evaluable).
Original docstring of v11 follows.

v11 (ADR-AIIM-447): v10 with the left panel changed from the funnel to the sequence of the work in four steps (one audit-and-correction
cycle, not a loop). Original docstring of v10 follows.

v10 (ADR-AIIM-440): the right-hand side redrawn so that it can be read without the paper.

Why. v9 put three cards of jargon on the right ("harness contract", "contracts as deployed", "snapshot comparator", "held-out faults",
"scores changed by equivalent SMILES", "re-analysis of the ligand-based step") with a number and two clipped labels in each, and arrows
that tied the cards to the funnel without saying how. A reader had to know the paper to read the numbers.

What changed.
  1. Each card on the right opens with the question it answers, in plain words.
  2. Each number is drawn as what it counts: bars out of 128 planted errors (card 2), dots out of 5 and out of 25 (card 3), and a
     logarithmic strip of the 5,000 compounds that shows the released value, the median and the range of the 200 reorderings (card 4).
  3. The arrows are gone; the card that concerns the ligand-based step says so in its title and counts out of the same 5,000 as the
     funnel row. The footer states what is reported, not what "must" be done.
Panel 1 is the v9 funnel with the same numbers and a plainer title. The canvas is 590 x 236 pt (2.5 to 1), as in v9.

Every number is a literal copied from the manuscript text (Sections 2.7, 3.1, 3.5, 3.6 and the abstract); nothing here recomputes
them, apart from the positions on the logarithmic strip. Wording is neutral: no compound is described as a candidate or as prioritised.

    python build_graphical_abstract_v12.py
"""
from __future__ import annotations

import math
import subprocess
import sys
from pathlib import Path

from PIL import Image

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from build_graphical_abstract_v7 import (  # noqa: E402
    AMBER, AMBER_F, BLUE, BLUE_F, FILL, INK, MOSS, MOSS_F, RULE, SOFT, VIOLET, VIOLET_F, SCRIPTS,
    Svg, circle, fits, path)

GW, GH = 590.0, 236.0
STEM = "figure_graphical_abstract_v12_2026-10-14"

ROWS = [("25,604", "biosynthetic-superclass matches"),
        ("21,882", "after standardisation and filters"),
        ("8,000", "2D pharmacophore screen"),
        ("5,000", "3D refinement"),
        ("500", "best reduced-docking scores"),
        ("497", "completed full-protocol docking"),
        ("482", "with an on-pocket pose"),
        ("3", "met the original two-part rule")]


def build(out: Path) -> Path:
    s = Svg(GW, GH, "Auditing agent-generated biomedical code with executable contracts",
            "Left, the work in four steps, once: (1) build and freeze, an agent-written A2A receptor screen of an expanded library, 497 docked, 3 met the "
            "original two-part rule; (2) audit; (3) correct, once, human-approved and pre-registered, returning to the SNP-Lib, 236, 69, 15, then 8 compounds; "
            "(4) report as persistence, 3 of the 15 selected compounds persist under every choice tested, which is not a prediction of activity. Right, the "
            "audit as three questions. Do the checks catch planted errors? Of 128 planted errors, the contract checks caught 124 with the benchmark's check "
            "set and 101 as deployed, and a simple file-comparison test caught 128. What lies outside the written checks? Errors written blind to "
            "the project were caught in 1 of 5 cases in a first round and in 3 of 6 evaluable cases in a second round, and rewriting the same molecules as equivalent SMILES strings changed 25 of 25 docking scores. Is the final selection "
            "stable? Of 5,000 compounds, 72 passed one filter in the released run; over 200 reorderings of the same data the number ranged from 6 to 4,119 "
            "(median 1,004).")
    s.text(GW / 2, 14.5, fits("Auditing agent-generated biomedical code with executable contracts", 11.5, GW - 12, bold=True),
           size=11.5, weight="bold")

    top, bot = 23.0, 206.0
    # ------------------------------------------------------------ left: the sequence of the work, once
    px, pw = 6.0, 216.0
    s.box(px, top, pw, bot - top, fill="#FFFFFF", stroke=BLUE, sw=0.9, r=3)
    s.box(px, top, pw, 14.5, fill=BLUE_F, stroke=BLUE, sw=0.9, r=3)
    s.text(px + 8, top + 9.9, fits("The work in four steps, once", 8.6, pw - 12, bold=True), size=8.6, weight="bold", fill=BLUE, anchor="start")
    steps = [("1  Build and freeze", ["an agent-written screen of an expanded library:", "497 docked, 3 met the original two-part rule"]),
             ("2  Audit", ["checks against planted and blind errors, rewritten", "molecules and reordered data (the three panels at right)"]),
             ("3  Correct, once", ["human-approved, each pre-registered and tested;", "back to the SNP-Lib: 236, 69, 15, then 8 compounds"]),
             ("4  Report as persistence", ["3 of the 15 selected compounds persist under every", "choice tested; persistence is not a prediction of activity"])]
    sx, sw_, sh, sgap, sy0 = px + 7, pw - 14.0, 34.0, 8.0, top + 22.0
    for k, (head, lines) in enumerate(steps):
        yy = sy0 + k * (sh + sgap)
        last = k == len(steps) - 1
        s.box(sx, yy, sw_, sh, fill=BLUE_F if last else FILL, stroke=BLUE if last else RULE, sw=0.9 if last else 0.6, r=2.5)
        s.text(sx + 7, yy + 10.8, fits(head, 8.8, sw_ - 12, bold=True), size=8.8, weight="bold", fill=BLUE if last else INK, anchor="start")
        for m, ln in enumerate(lines):
            s.text(sx + 7, yy + 20.6 + m * 8.8, fits(ln, 7.4, sw_ - 12), size=7.4, fill=SOFT if not last else INK, anchor="start")
        if k < len(steps) - 1:
            s.arrow(sx + sw_ / 2, yy + sh + 0.4, sx + sw_ / 2, yy + sh + sgap - 0.4, sw=0.7)

    # ------------------------------------------------------------ right: three questions
    lx, gap = 232.0, 6.5
    lw = GW - 6.0 - lx
    lh = (bot - top - 2 * gap) / 3
    cards = []
    for k, (col, fil, title) in enumerate([
            (MOSS, MOSS_F, "Audit: do the checks catch planted errors?"),
            (AMBER, AMBER_F, "Audit: what lies outside the written checks?"),
            (VIOLET, VIOLET_F, "Audit: is the final selection stable? (5,000 compounds)")]):
        y = top + k * (lh + gap)
        s.box(lx, y, lw, lh, fill="#FFFFFF", stroke=col, sw=0.9, r=3)
        s.box(lx, y, lw, 14.5, fill=fil, stroke=col, sw=0.9, r=3)
        s.text(lx + 8, y + 10.0, fits(title, 8.6, lw - 12, bold=True), size=8.6, weight="bold", fill=col, anchor="start")
        cards.append((y, col, fil))

    # ---- card 2: bars out of 128 planted errors
    y, col, fil = cards[0]
    bx, bw, bh = lx + 150.0, 128.0, 6.4
    s.text(lx + 8, y + 14.5 + 7.8, "128 errors deliberately planted in the workflow", size=7.4, fill=SOFT, anchor="start")
    for i, (label, n, c) in enumerate([("contract checks, benchmark set", 124, col),
                                       ("contract checks, as deployed", 101, col),
                                       ("simple file comparison", 128, "#7C8A82")]):
        yy = y + 14.5 + 13.4 + i * 9.4
        s.text(lx + 8, yy + 5.5, fits(label, 7.8, 138, bold=True), size=7.8, weight="bold", anchor="start")
        s.box(bx, yy, bw, bh, fill="#FFFFFF", stroke=RULE, sw=0.5, r=1.2)
        s.box(bx, yy, bw * n / 128, bh, fill=c, stroke=c, sw=0.5, r=1.2)
        s.text(bx + bw + 6, yy + 5.9, f"{n} / 128", size=8.6, weight="bold", fill=c, anchor="start")

    # ---- card 3: dots out of 5 and out of 25
    y, col, fil = cards[1]
    r0 = y + 14.5 + (lh - 14.5) / 2
    tx = lx + 70
    for i in range(11):                       # first round: 5 cases, 1 caught; second round: 6 evaluable cases, 3 caught
        gx = lx + 9 + (i + (0.7 if i >= 5 else 0)) * 4.8
        caught = i == 0 or (i >= 5 and i - 5 < 3)
        circle(s, gx, r0 - 9.6, 2.0, stroke=col, sw=0.8, fill=col if caught else "#FFFFFF")
    s.text(tx, r0 - 6.6, fits("errors written blind caught: 1 of 5, then 3 of 6", 7.8, lw - 74, bold=True),
           size=7.8, weight="bold", anchor="start")
    for i in range(25):
        s.box(lx + 8 + i * 2.2, r0 + 1.4, 1.8, 5.4, fill=col, stroke=col, sw=0.2, r=0.4)
    s.text(tx, r0 + 6.6, fits("25 of 25 docking scores changed when molecules were rewritten", 7.8, lw - 74, bold=True),
           size=7.8, weight="bold", anchor="start")
    s.text(tx, r0 + 15.4, fits("(same structures, written as different SMILES strings)", 7.2, lw - 74), size=7.2, fill=SOFT,
           anchor="start")

    # ---- card 4: logarithmic strip of the number of compounds passing the filter
    y, col, fil = cards[2]
    ax0, ax1 = lx + 20.0, lx + lw - 20.0
    sy = y + 14.5 + 17.6

    def X(v):
        return ax0 + (ax1 - ax0) * math.log10(v) / math.log10(5000)

    s.box(ax0, sy, ax1 - ax0, 6.0, fill="#FFFFFF", stroke=RULE, sw=0.5, r=1.2)
    s.box(X(6), sy, X(4119) - X(6), 6.0, fill=fil, stroke=col, sw=0.7, r=1.2)
    s.line(X(1004), sy - 1.5, X(1004), sy + 7.5, stroke=col, sw=1.1)
    path(s, f"M{X(72):g},{sy - 1.0:g} L{X(72) - 3.0:g},{sy - 6.5:g} L{X(72) + 3.0:g},{sy - 6.5:g} z", stroke=col, sw=0.6, fill=col)
    s.text(X(72), sy - 8.6, "72 released run", size=7.6, weight="bold", fill=col)
    s.text(X(1004) + 2.0, sy - 8.6, "1,004 median", size=7.6, weight="bold", fill=col, anchor="start")
    for v, lab in ((1, "1"), (10, "10"), (100, "100"), (1000, "1,000"), (5000, "5,000")):
        s.line(X(v), sy + 6.0, X(v), sy + 8.2, stroke=RULE, sw=0.5)
        s.text(X(v), sy + 14.6, lab, size=7.2, fill=SOFT)
    s.text(X(6) + 1.0, sy + 3.8, "6", size=7.0, fill=col, anchor="end") if False else None
    s.text(lx + lw / 2, y + lh - 1.9,
           fits("compounds passing one filter, over 200 reorderings of the same data", 7.4, lw - 12), size=7.4, fill=SOFT)

    s.box(6, 213, GW - 12, 18.5, fill=FILL, stroke=RULE, sw=0.7)
    s.text(GW / 2, 225.2, fits("Contract coverage, decision stability and biomedical validation are reported separately.",
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

"""Graphical abstract, v14 (ADR-AIIM-455): v13 with plain-language captions (a reader without the article should be able to say what was done and found), a stated
basis for each outcome, and the focus-library re-analysis named. Original docstring of v13 follows.

v13: the method row of v7 (four stages, each drawn as the act it is) with an outcome drawn under each stage, and a strip for the
one correction cycle and the persistence of the selected compounds; numbers are kept small and the outcomes are drawn as bars, dots, a funnel and a strip.

Why. v12 told the order of the work (build, audit, correct, report) and dropped what the framework is (the contracts, the gates, the seeded and the blind faults, the
reanalyses). The author asked for an infographic that is read as a picture, keeps the key words, and shows less text and fewer numbers. v13 takes the pictograms of v7
(imported from build_graphical_abstract_v7.py, so the glyph rules of that file apply: gates numbered and not ticked, the ROC glyph at 0.60, a run record and not source
code, a generic molecule), widens them to the 590 x 236 pt frame of the article (2.5 to 1), and replaces each outcome tile by a drawing of its result.

Every number is a literal copied from the manuscript text (Sections 2.2, 3.1, 3.5, 3.6, 3.8.4 and the abstract); nothing here recomputes them apart from positions.
No compound is described as a candidate; persistence is described as not a prediction of activity.

    python build_graphical_abstract_v14.py
"""
from __future__ import annotations

import math
import subprocess
import sys
from pathlib import Path

from PIL import Image

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import build_graphical_abstract_v7 as V7  # noqa: E402
from build_graphical_abstract_v7 import (  # noqa: E402
    AMBER, AMBER_F, BLUE, BLUE_F, FILL, INK, MOSS, MOSS_F, RULE, SOFT, VIOLET, VIOLET_F, Svg, SCRIPTS,
    bars_glyph, bug, chip, circle, code_page, contract_page, fits, molecule, path, receptor, record_page, roc_glyph, shuffle_glyph)

GW, GH = 590.0, 236.0
STEM = "figure_graphical_abstract_v14_2026-10-14"
GREY = "#7C8A82"


def build(out: Path) -> Path:
    s = Svg(GW, GH, "Layered verification of LLM-agent-generated drug-discovery workflows",
            "Four stages drawn left to right, each with its outcome below. 1, an LLM agent writes the code of an A2A receptor virtual screen; the decision under audit "
            "released 3 of 497 docked compounds. 2, ten executable contracts over nine gates, 94 checks, tested on 128 seeded faults: 124 detected under the benchmark's "
            "harness contract, 101 as deployed, and a simple snapshot comparison 128. 3, the contracts are challenged with blind-authored faults and metamorphic relations: "
            "1 of 5 blind faults caught in a first round and 3 of 6 evaluable in a second, and rewriting molecules as equivalent SMILES changed all 25 docking scores. "
            "4, the decision is reanalysed over the benchmark, eight descriptor spaces and the order of the rows: 72 of 5,000 compounds passed one filter in the released "
            "run, against a median of 1,004 over 200 reorderings. After one human-approved correction cycle, 3 of the 15 selected compounds persist under every analytic "
            "choice tested; persistence is not a prediction of activity. Contract coverage, decision stability and biomedical validation are reported separately.")
    s.text(GW / 2, 14.5, fits("Layered verification of LLM-agent-generated drug-discovery workflows", 11.5, GW - 12, bold=True), size=11.5, weight="bold")

    colw, gap = 137.0, 12.0
    total = 4 * colw + 3 * gap
    x0 = (GW - total) / 2
    k = colw / 117.25                     # the glyph offsets of v7 were laid out for a column of 117.25 pt
    inner = colw - 10
    top, ch, band, R = 22.0, 86.0, 14.0, 3.0
    cy = top + 33.0
    lab = top + 58.0
    CLEAR = 2.0

    def panel(i, col, fil, head, caption):
        x = x0 + i * (colw + gap)
        s.parts.append(f'  <path d="M{x:g},{top + band:g} V{top + R:g} A{R:g},{R:g} 0 0 1 {x + R:g},{top:g} H{x + colw - R:g} A{R:g},{R:g} 0 0 1 '
                       f'{x + colw:g},{top + R:g} V{top + band:g} Z" fill="{fil}"/>\n')
        s.line(x, top + band, x + colw, top + band, stroke=col, sw=0.7)
        s.box(x, top, colw, ch, fill="none", stroke=col, sw=0.9, r=R)
        s.text(x + 8, top + 10.0, fits(f"{i + 1}  {head}", 7.8, inner, bold=True), size=7.8, weight="bold", fill=col, anchor="start")
        s.lines(x + colw / 2, top + 70, [fits(t, 7.2, inner) for t in caption], size=7.2, lead=1.4)
        if i:
            s.arrow(x - gap + 1.6, cy, x - CLEAR - 0.5, cy, sw=0.65)
        return x

    def label(x, txt, col=SOFT):
        s.text(x, lab, txt, size=7, fill=col)

    # ------------------------------------------------------------------ row A: what was done
    x = panel(0, BLUE, BLUE_F, "Generate the workflow", ["An LLM agent writes a virtual screen", "of natural products for A2A receptor"])
    molecule(s, x + 22 * k, cy - 10.5, 4.4, BLUE)
    receptor(s, x + 22 * k, cy + 9, BLUE, BLUE_F)
    s.arrow(x + 42 * k, cy, x + 49 * k, cy, sw=0.65)
    chip(s, x + 66 * k, cy, BLUE, BLUE_F, "LLM")
    s.arrow(x + 82 * k, cy, x + 89 * k, cy, sw=0.65)
    code_page(s, x + 104 * k, cy, BLUE)
    label(x + 22 * k, "inputs")
    label(x + 66 * k, "agent")
    label(x + 104 * k, "code")

    x = panel(1, MOSS, MOSS_F, "Specify and test contracts", ["Automated checks (contracts) for nine", "gates, tested on 128 planted errors"])
    sq, sg = 9.6, 2.6
    gx = x + (colw - (9 * sq + 8 * sg)) / 2
    gy, rail2 = cy + 6.5, cy + 3.0
    mid = [gx + j * (sq + sg) + sq / 2 for j in range(9)]
    contract_page(s, x + 34, cy - 9, MOSS)
    s.text(x + 48, cy - 11, "executable", size=7, fill=SOFT, anchor="start")
    s.text(x + 48, cy - 3, "contracts", size=7, fill=SOFT, anchor="start")
    path(s, f"M{mid[0]:.2f},{gy:g} V{rail2:g} H{mid[8]:.2f} V{gy:g} " + " ".join(f"M{m:.2f},{rail2:g} V{gy:g}" for m in mid[1:8]) + f" M{x + 34:g},{cy + 0.5:g} V{rail2:g}",
         stroke=MOSS, sw=0.6)
    for j, m in enumerate(mid):
        s.box(m - sq / 2, gy, sq, sq, fill=MOSS_F, stroke=MOSS, sw=0.7, r=1.6)
        s.text(m, gy + 7.3, str(j), size=7, weight="bold", fill=MOSS)
    s.text(x + colw / 2, lab, "gates G0 to G8", size=7, fill=SOFT)

    x = panel(2, AMBER, AMBER_F, "Challenge the contracts", ["Try to fool the checks: blind-written", "errors and metamorphic tests"])
    bug(s, x + 20 * k, cy, AMBER, AMBER_F)
    s.arrow(x + 32 * k, cy, x + 45 * k, cy, sw=0.65)
    record_page(s, x + 58 * k, cy, AMBER, hot=AMBER)
    s.arrow(x + 71 * k, cy, x + 90 * k - 5, cy, sw=0.65)
    s.box(x + 90 * k - 3, cy - 10.5, 21, 21, fill=MOSS_F, stroke=MOSS, sw=0.9)
    s.text(x + 90 * k + 7.5, cy + 4.6, "?", size=13, weight="bold", fill=AMBER)
    label(x + 20 * k, "fault")
    label(x + 58 * k, "faulted run")
    label(x + 90 * k + 7.5, "detected?")

    x = panel(3, VIOLET, VIOLET_F, "Reanalyse the decision", ["Does the decision hold under changes", "that should not matter?"])
    roc_glyph(s, x + 23 * k, cy, VIOLET)
    bars_glyph(s, x + colw / 2, cy, VIOLET)
    shuffle_glyph(s, x + colw - 24 * k, cy, VIOLET, VIOLET_F)
    label(x + 23 * k, "benchmark")
    label(x + colw / 2, "8 spaces")
    label(x + colw - 24 * k, "row order")

    # ------------------------------------------------------------------ row B: what was found, drawn
    top2, ch2 = 116.0, 62.0
    start, tip = top + ch + 0.45, top2 - CLEAR - 0.5
    for i in range(4):
        xc = x0 + i * (colw + gap) + colw / 2
        path(s, f"M{xc:.2f},{start:g} V{tip:g}", sw=0.65, arrow=True, cap="butt")
    for i, col in enumerate((BLUE, MOSS, AMBER, VIOLET)):
        s.box(x0 + i * (colw + gap), top2, colw, ch2, fill="#FFFFFF", stroke=col, sw=0.9)

    # B1: a funnel that ends in the decision under audit (the dots illustrate a narrowing; they are not a count)
    x = x0
    fx, fy = x + colw / 2, top2 + 8
    import random
    rnd = random.Random(7)
    path(s, f"M{fx - 34:g},{fy:g} L{fx + 34:g},{fy:g} L{fx + 7:g},{fy + 24:g} L{fx - 7:g},{fy + 24:g} Z", stroke=BLUE, sw=0.8, fill=BLUE_F)
    for _ in range(30):
        px = fx + rnd.uniform(-30, 30)
        py = fy + rnd.uniform(1.5, 8.5)
        circle(s, px, py, 0.9, stroke=BLUE, sw=0.3, fill=BLUE)
    for j in range(3):
        circle(s, fx + (j - 1) * 5.2, fy + 29.5, 2.3, stroke=BLUE, sw=0.8, fill=BLUE)
    s.text(fx, top2 + 49.0, fits("497 docked, 3 passed", 7.8, colw - 10, bold=True), size=7.8, weight="bold", fill=BLUE)
    s.text(fx, top2 + 57.0, fits("the original selection rule (audited)", 7.2, colw - 10), size=7.2, fill=SOFT)

    # B2: planted errors caught, out of 128
    x = x0 + (colw + gap)
    bx, bw, bh = x + 82, 32.0, 6.0
    s.text(x + 8, top2 + 10.5, "planted errors caught, of 128", size=7.2, fill=SOFT, anchor="start")
    for j, (name, n, c) in enumerate((("checks, benchmark", 124, MOSS), ("checks, as deployed", 101, MOSS), ("file comparison", 128, GREY))):
        yy = top2 + 15.5 + j * 10.5
        s.text(x + 8, yy + 5.4, fits(name, 7.2, 72, bold=True), size=7.2, weight="bold", anchor="start")
        s.box(bx, yy, bw, bh, fill="#FFFFFF", stroke=RULE, sw=0.5, r=1.2)
        s.box(bx, yy, bw * n / 128, bh, fill=c, stroke=c, sw=0.5, r=1.2)
        s.text(bx + bw + 4, yy + 5.5, str(n), size=7.8, weight="bold", fill=c, anchor="start")
    s.text(x + colw / 2, top2 + 50.0, fits("contracts missed some; a plain file", 7.0, colw - 8), size=7.0, fill=SOFT)
    s.text(x + colw / 2, top2 + 57.0, fits("comparison caught all 128 (benchmark)", 7.0, colw - 8), size=7.0, fill=SOFT)

    # B3: blind-written errors caught in two rounds, and the SMILES relation
    x = x0 + 2 * (colw + gap)
    s.text(x + 8, top2 + 11.5, "first set", size=7.4, weight="bold", anchor="start")
    for j in range(5):
        circle(s, x + 52 + j * 8.6, top2 + 9.2, 3.0, stroke=AMBER, sw=0.8, fill=AMBER if j == 0 else "#FFFFFF")
    s.text(x + 52 + 5 * 8.6 - 2, top2 + 11.9, "1 of 5", size=7.0, fill=SOFT, anchor="start")
    s.text(x + 8, top2 + 23.5, "second set", size=7.4, weight="bold", anchor="start")
    for j in range(6):
        circle(s, x + 52 + j * 8.6, top2 + 21.2, 3.0, stroke=AMBER, sw=0.8, fill=AMBER if j < 3 else "#FFFFFF")
    s.text(x + 52 + 6 * 8.6 - 2, top2 + 23.9, "3 of 6", size=7.0, fill=SOFT, anchor="start")
    s.text(x + colw / 2, top2 + 35.0, fits("blind-written errors caught", 7.8, colw - 10, bold=True), size=7.8, weight="bold", fill=AMBER)
    for j in range(25):
        s.box(x + (colw - 25 * 2.3) / 2 + j * 2.3, top2 + 39.5, 1.9, 5.0, fill=AMBER, stroke=AMBER, sw=0.2, r=0.4)
    s.text(x + colw / 2, top2 + 51.0, fits("second set: 6 of 7 cases evaluable", 7.0, colw - 8), size=7.0, fill=SOFT)
    s.text(x + colw / 2, top2 + 58.0, fits("rewritten SMILES: all 25 scores changed", 7.0, colw - 8), size=7.0, fill=SOFT)

    # B4: a logarithmic strip of the compounds kept by one filter
    x = x0 + 3 * (colw + gap)
    ax0, ax1 = x + 14.0, x + colw - 14.0
    sy = top2 + 22.0

    def X(v):
        return ax0 + (ax1 - ax0) * math.log10(v) / math.log10(5000)

    s.box(ax0, sy, ax1 - ax0, 6.0, fill="#FFFFFF", stroke=RULE, sw=0.5, r=1.2)
    s.box(X(6), sy, X(4119) - X(6), 6.0, fill=VIOLET_F, stroke=VIOLET, sw=0.7, r=1.2)
    s.line(X(1004), sy - 1.5, X(1004), sy + 7.5, stroke=VIOLET, sw=1.1)
    path(s, f"M{X(72):g},{sy - 1.0:g} L{X(72) - 3.0:g},{sy - 6.5:g} L{X(72) + 3.0:g},{sy - 6.5:g} z", stroke=VIOLET, sw=0.6, fill=VIOLET)
    s.text(X(72), sy - 9.0, "original run", size=7.2, weight="bold", fill=VIOLET)
    s.text(X(1004) + 2.0, sy - 9.0, "median", size=7.2, weight="bold", fill=VIOLET, anchor="start")
    for v, lab_ in ((1, "1"), (100, "100"), (5000, "5,000")):
        s.line(X(v), sy + 6.0, X(v), sy + 8.2, stroke=RULE, sw=0.5)
        s.text(X(v), sy + 15.0, lab_, size=7.0, fill=SOFT)
    s.text(x + colw / 2, top2 + 47.0, fits("of 5,000 kept by one filter", 7.6, colw - 8, bold=True), size=7.4, weight="bold", fill=VIOLET)
    s.text(x + colw / 2, top2 + 57.0, fits("only the order of training rows changed", 7.0, colw - 8), size=7.0, fill=SOFT)

    # ------------------------------------------------------------------ row C: one correction cycle, persistence, and the three layers
    ctop, chh = 184.0, 48.0
    lw = total - 160.0 - gap               # left card; the right card is 160 pt wide
    s.box(x0, ctop, lw, chh, fill=FILL, stroke=RULE, sw=0.7)
    s.text(x0 + 8, ctop + 11.5, "After one human-approved correction cycle, a focus-library re-analysis selected 15 compounds", size=7.8, weight="bold", anchor="start")
    s.text(x0 + 8, ctop + 21.0, "persistence: how many of the 15 are still selected under every analytic choice tested", size=7.2, fill=SOFT, anchor="start")
    d0 = x0 + 12
    for j in range(15):
        circle(s, d0 + j * 11.0, ctop + 33.5, 3.3, stroke=BLUE, sw=0.8, fill=BLUE if j < 3 else "#FFFFFF")
    s.text(d0 + 15 * 11.0 + 2, ctop + 33.0, "3 of 15 persist", size=7.8, weight="bold", fill=BLUE, anchor="start")
    s.text(d0 + 15 * 11.0 + 2, ctop + 41.5, "not a prediction of activity", size=7.2, fill=SOFT, anchor="start")

    rx0 = x0 + lw + gap
    rw = total - lw - gap
    s.box(rx0, ctop, rw, chh, fill=FILL, stroke=RULE, sw=0.7)
    s.text(rx0 + rw / 2, ctop + 11.5, "Three answers, reported separately", size=7.8, weight="bold")
    for j, (txt, col, fil) in enumerate((("contract coverage", MOSS, MOSS_F), ("decision stability", VIOLET, VIOLET_F), ("biomedical validation", AMBER, AMBER_F))):
        yy = ctop + 16.5 + j * 10.2
        s.box(rx0 + 18, yy, rw - 36, 8.6, fill=fil, stroke=col, sw=0.6, r=2.0)
        s.text(rx0 + rw / 2, yy + 6.3, fits(txt, 7.2, rw - 40, bold=True), size=7.2, weight="bold", fill=col)

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

"""
build_graphical_abstract_v7.py -- the v7 graphical abstract as a workflow
picture, in the line style of Figures 1 and 2.

Why a redraw. The graphical abstract shipped in
AIIM_submission_package_v7_2026-10-03.zip was drawn with Pillow on its own
terms: 2.5 pt outlines, filled number discs, pill headers, block arrows and a
plum that is in no other figure. Beside Figures 1 and 2 it read as a figure
from a different paper. Its fourth outcome box also ran off the right edge of
the canvas, and the DOCX placed the 2.33:1 image in a 2.48:1 frame, so it was
stretched as well as clipped.

Why pictograms. A first redraw kept the v7 layout, which was eight boxes of
text: it told a reader what was counted but not what was done. Each stage now
carries a small line drawing of the act itself -- inputs to agent to code;
a contract over the nine gates; a fault pushed into a run; the three
reanalyses -- with one outcome under each stage, in that stage's colour.

The v7 graphic also had four tiles under four stages, but they did not
correspond: the two seeded-fault results sat under "Generate" and "Specify".
Here the snapshot result is the note on the contract tile it qualifies, and
the freed first tile carries what the generated workflow produced. Stage 2 is
"specify and test" so that the seeded-fault result under it is its own.

What each pictogram is allowed to say. A glyph in a published figure is a
claim, so each one is held to what the manuscript reports:
  - the gates are numbered, not ticked. A tick in every gate reads as "every
    gate passed", which is a result the specification stage does not report;
  - the ROC glyph encloses 0.60 of its square, the value reported, and is not
    bowed toward the corner;
  - the faulted item is a run record with one altered line, not source code,
    because the faults were seeded into benchmark workspaces;
  - the molecule is generic and depicts no compound in the paper.

This version does not restyle by hand. It imports the `Svg` class and the
palette from `build_figures_1_2.py`, so rule weights, corner radius, arrowhead,
font stack and colours are the same objects Figures 1 and 2 are drawn with.
The canvas is 511 x 206 pt, the frame the manuscript already reserves.

Every number is taken from the v7 graphic or the v7 manuscript text. "3 of
497" and "6 to 4,119 of 5,000" are from Section 3.6 and "4 metamorphic
relations" from Section 2.2; the rest were in the v7 graphic. "0 of 55 false
escalations", "236" and "0 / 3" were in the v7 graphic and are not in this
one, so the manuscript's sentence listing the graphic's four outcomes no
longer matches it. All are hand-written literals; nothing in this file
recomputes them.

It lives outside `30_manuscript/figures` and `10_new_system/scripts` on
purpose: `build_deposit.py` copies both, and the deposit's file count is cited
in the manuscript.

Usage:
    python build_graphical_abstract_v7.py
"""
from __future__ import annotations

import math
import subprocess
import sys
from pathlib import Path

from PIL import Image, ImageFont

HERE = Path(__file__).resolve().parent
SCRIPTS = HERE.parents[1] / "10_new_system" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from build_figures_1_2 import (  # noqa: E402
    AMBER, AMBER_F, BLUE, BLUE_F, FILL, INK, LINE, MOSS, MOSS_F, RULE, SOFT,
    VIOLET, VIOLET_F, Svg)

GW, GH = 511.0, 206.0
STEM = "figure_graphical_abstract_v7_2026-10-03"

# The font stack resolves to Arial on the machine that rasterises the figure,
# so widths are measured in Arial rather than estimated per character.
_FONTS = {False: ImageFont.truetype("arial.ttf", 100),
          True: ImageFont.truetype("arialbd.ttf", 100)}


def fits(txt: str, size: float, width: float, bold: bool = False) -> str:
    need = _FONTS[bold].getlength(txt) * size / 100
    if need > width:
        raise SystemExit(f"graphical abstract: {txt!r} needs {need:.0f} pt at "
                         f"{size} pt but has {width:.0f} pt")
    return txt


# ------------------------------------------------------------- pictograms
# `Svg` draws boxes, text and straight lines. The glyphs need paths, so they
# are appended as raw elements with the same stroke weights the class uses.
def path(s: Svg, d: str, stroke=LINE, sw=0.8, fill="none", dash=None,
         arrow=False, join="round", cap="round"):
    da = f' stroke-dasharray="{dash}"' if dash else ""
    mk = ' marker-end="url(#a)"' if arrow else ""
    s.parts.append(
        f'  <path d="{d}" fill="{fill}" stroke="{stroke}" stroke-width="{sw:g}" '
        f'stroke-linecap="{cap}" stroke-linejoin="{join}"{da}{mk}/>\n')


def circle(s: Svg, cx, cy, r, stroke=LINE, sw=0.8, fill="none"):
    s.parts.append(
        f'  <circle cx="{cx:g}" cy="{cy:g}" r="{r:g}" fill="{fill}" '
        f'stroke="{stroke}" stroke-width="{sw:g}"/>\n')


def _hex(cx, cy, r):
    return [(cx + r * math.cos(math.radians(90 + 60 * k)),
             cy - r * math.sin(math.radians(90 + 60 * k))) for k in range(6)]


def _poly(pts):
    return "M" + " L".join(f"{x:.2f},{y:.2f}" for x, y in pts) + " Z"


def molecule(s: Svg, cx, cy, r, col):
    """An aromatic ring fused to a saturated one, with substituents.

    Deliberately generic. It says "small organic molecule" and nothing about
    which scaffolds the library holds; no compound in the paper is drawn.
    """
    dx = r * math.sqrt(3) / 2
    lc, rc = (cx - dx, cy), (cx + dx, cy)
    left, right = _hex(*lc, r), _hex(*rc, r)
    path(s, _poly(left), stroke=col, fill="#FFFFFF")
    path(s, _poly(right), stroke=col, fill="#FFFFFF")
    # Kekule bonds on alternate edges of the left ring, inset toward its centre
    for k in (0, 2, 4):
        (ax, ay), (bx, by) = left[k], left[(k + 1) % 6]
        f, g = 0.24, 0.16
        px = ax + (lc[0] - ax) * f + (bx - ax) * g
        py = ay + (lc[1] - ay) * f + (by - ay) * g
        qx = bx + (lc[0] - bx) * f + (ax - bx) * g
        qy = by + (lc[1] - by) * f + (ay - by) * g
        path(s, f"M{px:.2f},{py:.2f} L{qx:.2f},{qy:.2f}", stroke=col, sw=0.6)
    for (vx, vy), (ox, oy) in ((left[1], lc), (right[4], rc), (right[0], rc)):
        ux, uy = vx - ox, vy - oy
        n = math.hypot(ux, uy)
        path(s, f"M{vx:.2f},{vy:.2f} L{vx + 3.2 * ux / n:.2f},"
                f"{vy + 3.2 * uy / n:.2f}", stroke=col)


def receptor(s: Svg, cx, cy, col, fil):
    """Seven helices joined by alternating loops across a bilayer: a GPCR."""
    pitch, hw, hh = 3.9, 1.35, 5.5
    xs = [cx + (i - 3) * pitch for i in range(7)]
    for y in (cy - 3.4, cy + 3.4):
        s.line(cx - 15.5, y, cx + 15.5, y, stroke=RULE, sw=0.6)
    for i in range(6):
        y, sweep = (cy - hh, 1) if i % 2 == 0 else (cy + hh, 0)
        path(s, f"M{xs[i]:.2f},{y:.2f} A{pitch / 2:g},{pitch / 2:g} 0 0 {sweep} "
                f"{xs[i + 1]:.2f},{y:.2f}", stroke=col, sw=0.7)
    for x in xs:
        s.box(x - hw, cy - hh, 2 * hw, 2 * hh, fill=fil, stroke=col, sw=0.7, r=hw)


def chip(s: Svg, cx, cy, col, fil, label):
    h = 10.5
    for k in (-5.5, 0, 5.5):
        path(s, f"M{cx - h - 2.6:g},{cy + k:g} L{cx - h:g},{cy + k:g} "
                f"M{cx + h:g},{cy + k:g} L{cx + h + 2.6:g},{cy + k:g} "
                f"M{cx + k:g},{cy - h - 2.6:g} L{cx + k:g},{cy - h:g} "
                f"M{cx + k:g},{cy + h:g} L{cx + k:g},{cy + h + 2.6:g}",
             stroke=col)
    s.box(cx - h, cy - h, 2 * h, 2 * h, fill=fil, stroke=col, sw=0.9)
    s.text(cx, cy + 2.5, label, size=7, weight="bold", fill=col)


def page(s: Svg, cx, cy, col, w=17.0, h=22.0):
    x, y, c = cx - w / 2, cy - h / 2, 5.0
    path(s, f"M{x:g},{y:g} L{x + w - c:g},{y:g} L{x + w:g},{y + c:g} "
            f"L{x + w:g},{y + h:g} L{x:g},{y + h:g} Z", stroke=col, sw=0.9,
         fill="#FFFFFF")
    path(s, f"M{x + w - c:g},{y:g} L{x + w - c:g},{y + c:g} L{x + w:g},{y + c:g}",
         stroke=col, sw=0.7)
    return x, y


def code_page(s: Svg, cx, cy, col):
    """A page carrying the angle-bracket mark: source code, not a document."""
    x, y = page(s, cx, cy, col, w=19.0, h=23.0)
    s.text(cx, y + 13.4, "</>", size=7, weight="bold", fill=col)
    for k, (ind, ln) in enumerate([(0, 9), (2.5, 8)]):
        yy = y + 16.9 + k * 3.0
        path(s, f"M{x + 3.5 + ind:g},{yy:g} L{x + 3.5 + ind + ln:g},{yy:g}",
             stroke=RULE)


def record_page(s: Svg, cx, cy, col, hot):
    """A run's output with one record altered: where the faults were seeded."""
    x, y = page(s, cx, cy, col)
    for k, ln in enumerate([7, 11, 9, 11, 8]):
        yy = y + 7.5 + k * 3.1
        path(s, f"M{x + 3:g},{yy:g} L{x + 3 + ln:g},{yy:g}",
             stroke=hot if k == 2 else RULE, sw=1.6 if k == 2 else 0.8)


def tick(s: Svg, x, y, col, k=1.0, sw=0.9):
    path(s, f"M{x:g},{y:g} l{1.3 * k:g},{1.5 * k:g} l{2.5 * k:g},{-3.1 * k:g}",
         stroke=col, sw=sw)


def contract_page(s: Svg, cx, cy, col):
    """A page of checked items: a contract is a list of checks."""
    x, y = page(s, cx, cy, col, h=19.0)
    for k in range(4):
        yy = y + 6.9 + k * 3.25
        tick(s, x + 2.8, yy - 0.4, col, k=0.62, sw=0.8)
        path(s, f"M{x + 7.6:g},{yy:g} L{x + 14:g},{yy:g}", stroke=RULE)


def bug(s: Svg, cx, cy, col, fil):
    for sx in (-1, 1):
        for y0, y1 in ((-1.2, -3.6), (1.8, 1.8), (4.8, 7.2)):
            path(s, f"M{cx + sx * 3.8:g},{cy + y0:g} L{cx + sx * 7.8:g},{cy + y1:g}",
                 stroke=col)
        path(s, f"M{cx + sx * 1.2:g},{cy - 6.6:g} L{cx + sx * 3.4:g},{cy - 10:g}",
             stroke=col)
    s.parts.append(
        f'  <ellipse cx="{cx:g}" cy="{cy + 1.8:g}" rx="4.4" ry="5.8" fill="{fil}" '
        f'stroke="{col}" stroke-width="0.9"/>\n')
    circle(s, cx, cy - 5.2, 2.5, stroke=col, sw=0.9, fill=fil)
    path(s, f"M{cx:g},{cy - 2.6:g} L{cx:g},{cy + 7.4:g}", stroke=col, sw=0.7)


def roc_glyph(s: Svg, cx, cy, col):
    """A ROC curve enclosing 0.60 of the unit square.

    The area between a quadratic Bezier and its chord is two thirds of the
    control triangle, so a control point at (0.35, 0.65) gives 0.5 + 0.1. The
    manuscript reports 0.601 for docking on the ChEMBL benchmark; a glyph
    bowed toward the corner would have drawn a better result than was found.
    """
    a = 10.0
    path(s, f"M{cx - a:g},{cy + a:g} L{cx + a:g},{cy - a:g}", stroke=RULE,
         sw=0.6, dash="1.6 1.6")
    path(s, f"M{cx - a:g},{cy + a:g} Q{cx - a + 0.35 * 2 * a:g},"
            f"{cy + a - 0.65 * 2 * a:g} {cx + a:g},{cy - a:g}", stroke=col, sw=1.0)
    path(s, f"M{cx - a:g},{cy - a:g} L{cx - a:g},{cy + a:g} L{cx + a:g},{cy + a:g}",
         stroke=LINE, join="miter", cap="square")


def bars_glyph(s: Svg, cx, cy, col):
    """Eight descriptor spaces, one of which admits most of the library."""
    base, pitch, bw = cy + 10.0, 2.8, 1.9
    for k, hgt in enumerate([0.7, 0.7, 0.7, 0.7, 1.1, 1.4, 1.8, 17.0]):
        x = cx - 11 + k * pitch + 0.2
        s.parts.append(
            f'  <rect x="{x:g}" y="{base - hgt:g}" width="{bw:g}" '
            f'height="{hgt:g}" fill="{col}"/>\n')
    path(s, f"M{cx - 12:g},{base:g} L{cx + 12:g},{base:g}", stroke=LINE)


def shuffle_glyph(s: Svg, cx, cy, col, fil):
    """The same rows in another order."""
    ys = (cy - 7.5, cy, cy + 7.5)
    for a, b in zip(ys, (ys[2], ys[0], ys[1])):
        path(s, f"M{cx - 5.2:g},{a:g} L{cx + 5.2:g},{b:g}", stroke=col, sw=0.7)
    for y in ys:
        for x in (cx - 12, cx + 5.2):
            s.box(x, y - 1.7, 6.8, 3.4, fill=fil, stroke=col, sw=0.7, r=1.2)


# ------------------------------------------------------- graphical abstract
def graphical_abstract(out: Path) -> Path:
    s = Svg(GW, GH,
            "Auditing agent-generated biomedical code with executable contracts",
            "A four-stage workflow with one outcome under each stage. An LLM "
            "agent writes an A2A receptor virtual screen of natural-product "
            "libraries, which released 3 of 497 docked compounds; ten "
            "executable contracts over nine gates detected 124 of 128 seeded "
            "faults, where deterministic snapshots detected all 128; 1 of 5 "
            "blind-authored faults was detected and equivalent SMILES changed "
            "all 25 docking scores; and the released run kept 72 of 5,000 "
            "compounds as eligible, against a median of 1,004 over 200 "
            "reorderings of the training rows.")

    s.text(GW / 2, 14,
           fits("Auditing agent-generated biomedical code with executable "
                "contracts", 11, GW - 12, bold=True), size=11, weight="bold")

    # An arrowhead is 3.9 pt long at this stroke weight. The 7 pt gap of
    # Figure 1's gate strip leaves a head with no shaft once it has to stop
    # short of a 0.9 pt coloured border, so the gap here is wider.
    colw, gap = 117.25, 11.0
    total = 4 * colw + 3 * gap
    x0 = (GW - total) / 2
    inner = colw - 10
    top, ch, band, R = 22.0, 92.0, 15.0, 3.0
    cy = top + 37.0                      # centre line of every pictogram
    lab = top + 62.0                     # baseline of the glyph labels
    CLEAR = 2.0                          # arrow tip to the border it points at

    def panel(i, col, fil, head, caption):
        x = x0 + i * (colw + gap)
        # One outline, drawn last. The header is a fill clipped to the rounded
        # top and closed by a single rule; a second stroked box here doubled
        # the side borders and left a notch where its lower corners curved in.
        s.parts.append(
            f'  <path d="M{x:g},{top + band:g} V{top + R:g} A{R:g},{R:g} 0 0 1 '
            f'{x + R:g},{top:g} H{x + colw - R:g} A{R:g},{R:g} 0 0 1 '
            f'{x + colw:g},{top + R:g} V{top + band:g} Z" fill="{fil}"/>\n')
        s.line(x, top + band, x + colw, top + band, stroke=col, sw=0.7)
        s.box(x, top, colw, ch, fill="none", stroke=col, sw=0.9, r=R)
        s.text(x + 8, top + 10.5, fits(f"{i + 1}  {head}", 7.6, inner, bold=True),
               size=7.6, weight="bold", fill=col, anchor="start")
        s.lines(x + colw / 2, top + 75, [fits(t, 7.2, inner) for t in caption],
                size=7.2, lead=1.4)
        if i:
            # On the pictogram centre line, so the flow reads straight through.
            s.arrow(x - gap + 1.6, cy, x - CLEAR - 0.5, cy, sw=0.65)
        return x

    def label(x, txt, col=SOFT):
        s.text(x, lab, txt, size=7, fill=col)

    # 1 -- what was built
    # No library is named here, and no library size. Checked against the run
    # outputs on 2026-10-03: all 500 compounds that reached full-protocol
    # docking are from the LOTUS expansion, none is from Orchidaceae
    # (`np_library_funnel_2026-08-15/e8_np_library_prioritization.csv`), and
    # no SNP-Lib docking result is reported in v6 or v7 -- the 236-compound
    # docking run is formally withdrawn and must not be read. Naming the
    # orchid library beside "3 compounds prioritised" would have said the
    # three came from it; they are Pongamia pinnata flavonoids.
    x = panel(0, BLUE, BLUE_F, "Generate the workflow",
              ["A2A receptor virtual screen of", "natural-product libraries"])
    molecule(s, x + 20, cy - 10.5, 4.4, BLUE)
    receptor(s, x + 20, cy + 9, BLUE, BLUE_F)
    s.arrow(x + 38.5, cy, x + 45.5, cy, sw=0.65)
    chip(s, x + 60.5, cy, BLUE, BLUE_F, "LLM")
    s.arrow(x + 75.5, cy, x + 82.5, cy, sw=0.65)
    code_page(s, x + 96.5, cy, BLUE)
    label(x + 20, "inputs")
    label(x + 60.5, "agent")
    label(x + 96.5, "code")

    # 2 -- how it was specified
    # The seeded-fault test is named in the stage because its result is the
    # tile underneath; a stage that only specified could not have produced it.
    x = panel(1, MOSS, MOSS_F, "Specify and test contracts",
              ["94 deterministic checks,", "tested on 128 seeded faults"])
    sq, sg = 9.6, 2.2
    gx = x + (colw - (9 * sq + 8 * sg)) / 2
    gy, rail2 = cy + 6.5, cy + 3.0
    mid = [gx + k * (sq + sg) + sq / 2 for k in range(9)]
    contract_page(s, x + 30, cy - 9, MOSS)
    s.text(x + 43, cy - 11, "10 executable", size=7, fill=SOFT, anchor="start")
    s.text(x + 43, cy - 3, "contracts", size=7, fill=SOFT, anchor="start")
    # One path, so the comb's corners are joined rather than abutted.
    path(s, f"M{mid[0]:.2f},{gy:g} V{rail2:g} H{mid[8]:.2f} V{gy:g} "
            + " ".join(f"M{m:.2f},{rail2:g} V{gy:g}" for m in mid[1:8])
            + f" M{x + 30:g},{cy + 0.5:g} V{rail2:g}", stroke=MOSS, sw=0.6)
    for k, m in enumerate(mid):
        s.box(m - sq / 2, gy, sq, sq, fill=MOSS_F, stroke=MOSS, sw=0.7, r=1.6)
        s.text(m, gy + 7.3, str(k), size=7, weight="bold", fill=MOSS)
    s.text(x + colw / 2, lab, "gates G0 to G8", size=7, fill=SOFT)

    # 3 -- how it was challenged
    # The seeded faults are reported under stage 2, whose contracts they
    # measure; this stage is the part the contracts were not written for.
    x = panel(2, AMBER, AMBER_F, "Challenge the contracts",
              ["5 blind-authored faults,", "4 metamorphic relations"])
    bug(s, x + 19, cy, AMBER, AMBER_F)
    s.arrow(x + 30.5, cy, x + 43, cy, sw=0.65)
    record_page(s, x + 56, cy, AMBER, hot=AMBER)
    s.arrow(x + 68.5, cy, x + 81, cy, sw=0.65)
    # The gate is drawn as panel 2 draws it, so it reads as the same object.
    s.box(x + 85.5, cy - 10.5, 21, 21, fill=MOSS_F, stroke=MOSS, sw=0.9)
    s.text(x + 96, cy + 4.6, "?", size=13, weight="bold", fill=AMBER)
    label(x + 19, "fault")
    label(x + 56, "faulted run")
    label(x + 96, "detected?")

    # 4 -- how the decision was re-examined
    x = panel(3, VIOLET, VIOLET_F, "Reanalyse the decision",
              ["Discrimination, applicability", "domain, selection stability"])
    roc_glyph(s, x + 20.5, cy, VIOLET)
    bars_glyph(s, x + colw / 2, cy, VIOLET)
    shuffle_glyph(s, x + colw - 21, cy, VIOLET, VIOLET_F)
    label(x + 20.5, "benchmark")
    label(x + colw / 2, "8 spaces")
    label(x + colw - 21, "row order")

    # One outcome under each stage, in that stage's colour, joined by one
    # arrow. A reader takes a tile to be the result of the box above it, so
    # the layout has to make that true: an earlier draft hung three tiles from
    # stage 3 and read as a mistake.
    #
    # The snapshot result is the note on the contract tile, not a tile of its
    # own. It is the comparator for the same 128 faults, and the manuscript
    # reads it as a limit on the 96.9 %, not as a property of the workflow.
    #
    # The last tile is the row-order range, not "0 / 3". A bare zero under
    # the final stage read as the study failing, and "survived the
    # better-informed analysis" did not say what had been done. The range is
    # the instability itself: same molecules, labels, model, seed and split
    # code, 200 orderings of the training rows.
    top2, ch2 = 128.0, 43.0
    start, tip = top + ch + 0.45, top2 - CLEAR - 0.5
    for i in range(4):
        xc = x0 + i * (colw + gap) + colw / 2
        path(s, f"M{xc:.2f},{start:g} V{tip:g}", sw=0.65, arrow=True, cap="butt")

    outcomes = [
        (BLUE, "3 / 497", "compounds released",
         "the decision under audit"),
        (MOSS, "124 / 128", "seeded faults detected",
         "snapshots detected 128 / 128"),
        (AMBER, "1 / 5", "blind faults detected",
         "SMILES invariance also failed"),
        (VIOLET, "72 vs 1,004", "eligible compounds of 5,000",
         "released vs 200-order median"),
    ]
    for i, (col, big, text, note) in enumerate(outcomes):
        x = x0 + i * (colw + gap)
        s.box(x, top2, colw, ch2, fill="#FFFFFF", stroke=col, sw=0.9)
        s.text(x + colw / 2, top2 + 19.5, fits(big, 17, inner, bold=True),
               size=17, weight="bold", fill=col)
        s.text(x + colw / 2, top2 + 29.5, fits(text, 7.6, inner, bold=True),
               size=7.6, weight="bold")
        s.text(x + colw / 2, top2 + 38.5, fits(note, 7.2, inner), size=7.2,
               fill=SOFT)

    s.box(x0, 178, total, 18, fill=FILL, stroke=RULE, sw=0.7)
    s.text(GW / 2, 189.8,
           fits("Contract coverage, decision stability and biomedical "
                "validation must be evaluated separately.", 8, total - 16,
                bold=True), size=8, weight="bold")

    p = out / f"{STEM}.svg"
    s.save(p)
    return p


def main() -> int:
    svg = graphical_abstract(HERE)
    png = svg.with_suffix(".png")
    r = subprocess.run([sys.executable, str(SCRIPTS / "render_svg_to_png.py"),
                        str(svg), str(png)])
    if r.returncode:
        return r.returncode
    # Figures 3-8 carry their resolution in the file; Chrome's screenshot does
    # not write one, and a submission system reads a PNG without it as 72 dpi.
    with Image.open(png) as im:
        im.convert("RGB").save(png, dpi=(600, 600), optimize=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""
build_figures_1_2.py -- redraw Figures 1 and 2 at journal type size.

Why a redraw and not a re-export. The previous figures were authored on a
1500-unit canvas with 8.5-13 unit labels. Placed at Elsevier's 190 mm
double-column width those render at 3.1-3.8 pt against a 7 pt floor, and no
choice of dpi changes that -- 600-dpi 3.4 pt type is still 3.4 pt type. The
type has to grow relative to the canvas, which means each figure has to carry
less.

The fix here splits the old Figure 1, which was two diagrams stacked: a
pipeline overview on top and a single gate's internals below. They become
Figure 1 and Figure 2. Halving the content is what buys the type size.

The old Figure 2 (SBDD vs LBDD decision flow) is NOT redrawn, and the reason
is not its size. It ends in "D. Neither validates -> HALT / Refuse ranking"
and states "G4-G8 were not executed for this application". Version 6 reports
the opposite on both counts: G4 was executed and passed 7 of 7 checks, G6 was
executed, and a confidence-tiered prioritisation was released. The figure
depicts a superseded conclusion and is retired rather than enlarged.

Everything is laid out in POINTS. The viewBox is 539 pt wide, which is 190 mm,
so a font-size of 8 in this file is 8 pt on the printed page -- no mental
arithmetic, and no way to accidentally drop below the floor.

Usage:
    python build_figures_1_2.py [--out DIR]
"""
from __future__ import annotations

import argparse
from pathlib import Path

W = 539.0          # 190 mm in points: Elsevier double-column width
MIN_PT = 7.0       # nothing below this

INK = "#1B2420"
SOFT = "#55625B"
RULE = "#9AA79F"
FILL = "#EEF2EE"
LINE = "#6E7C74"
BLUE = "#2C5F8A"
BLUE_F = "#E3EBF3"
MOSS = "#4F7130"
MOSS_F = "#E7EEDC"
AMBER = "#A85C18"
AMBER_F = "#F7E9D8"
ALERT = "#9B2226"
ALERT_F = "#F7E2E2"
VIOLET = "#6B4E8C"
VIOLET_F = "#EFE8F5"

FONT = "Helvetica, Arial, sans-serif"


def esc(t: str) -> str:
    return (t.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


class Svg:
    def __init__(self, w: float, h: float, title: str, desc: str):
        self.h = h
        self.parts: list[str] = []
        self.head = (
            f'<svg viewBox="0 0 {w:g} {h:g}" width="{w/72*25.4:.1f}mm" '
            f'height="{h/72*25.4:.1f}mm" xmlns="http://www.w3.org/2000/svg" '
            f'role="img" font-family="{FONT}">\n'
            f"  <title>{esc(title)}</title>\n  <desc>{esc(desc)}</desc>\n"
            f'  <rect width="{w:g}" height="{h:g}" fill="#FFFFFF"/>\n'
            '  <defs><marker id="a" viewBox="0 0 8 8" refX="7" refY="4" '
            'markerWidth="6" markerHeight="6" orient="auto-start-reverse">'
            f'<path d="M0,0.5 L8,4 L0,7.5 z" fill="{LINE}"/></marker></defs>\n'
        )

    def box(self, x, y, w, h, fill=FILL, stroke=RULE, sw=0.7, r=3):
        self.parts.append(
            f'  <rect x="{x:g}" y="{y:g}" width="{w:g}" height="{h:g}" rx="{r:g}" '
            f'fill="{fill}" stroke="{stroke}" stroke-width="{sw:g}"/>\n')

    def text(self, x, y, s, size=8, fill=INK, weight="normal", anchor="middle",
             style="normal", spacing=None):
        if size < MIN_PT:
            raise SystemExit(f"font size {size} pt is below the {MIN_PT} pt floor: {s!r}")
        sp = f' letter-spacing="{spacing:g}"' if spacing else ""
        self.parts.append(
            f'  <text x="{x:g}" y="{y:g}" font-size="{size:g}" fill="{fill}" '
            f'font-weight="{weight}" font-style="{style}" text-anchor="{anchor}"'
            f'{sp}>{esc(s)}</text>\n')

    def lines(self, x, y, rows, size=8, fill=SOFT, lead=1.28, **kw):
        for i, r in enumerate(rows):
            self.text(x, y + i * size * lead, r, size=size, fill=fill, **kw)

    def arrow(self, x1, y1, x2, y2, stroke=LINE, sw=0.7, dash=None):
        d = f' stroke-dasharray="{dash}"' if dash else ""
        self.parts.append(
            f'  <line x1="{x1:g}" y1="{y1:g}" x2="{x2:g}" y2="{y2:g}" '
            f'stroke="{stroke}" stroke-width="{sw:g}"{d} marker-end="url(#a)"/>\n')

    def line(self, x1, y1, x2, y2, stroke=RULE, sw=0.6, dash=None):
        d = f' stroke-dasharray="{dash}"' if dash else ""
        self.parts.append(
            f'  <line x1="{x1:g}" y1="{y1:g}" x2="{x2:g}" y2="{y2:g}" '
            f'stroke="{stroke}" stroke-width="{sw:g}"{d}/>\n')

    def save(self, p: Path):
        p.write_text(self.head + "".join(self.parts) + "</svg>\n", encoding="utf-8")


# ---------------------------------------------------------------- Figure 1
def figure1(out: Path) -> Path:
    s = Svg(W, 250,
            "Architecture of the audited agentic pipeline",
            "Frozen hash-anchored inputs feed nine sequential gates, G0 input "
            "provenance through G8 reporting integrity. Every gate writes to an "
            "append-only artifact registry.")

    s.text(W / 2, 16, "Nine gates, ten contracts, one append-only record", size=11.5,
           weight="bold")
    s.text(W / 2, 29, "Every stage runs under a frozen contract the acting agent cannot modify",
           size=8, fill=SOFT)

    # frozen inputs
    s.box(18, 42, 154, 50, fill="#FFFFFF", stroke=INK, sw=0.9)
    s.text(95, 57, "Frozen inputs", size=9, weight="bold")
    s.lines(95, 69, ["SNP-Lib natural products (n = 236)",
                     "A2AR co-crystals · ChEMBL benchmark"], size=7.4)
    s.text(95, 86, "hash-anchored, read-only", size=7.2, fill=BLUE, style="italic")
    s.arrow(95, 93, 95, 104)

    # gate strip
    gates = [("G0", "Input", "provenance"), ("G1", "Standard-", "isation"),
             ("G2", "Prepar-", "ation"), ("G3", "Pose", "validation"),
             ("G4", "Scoring", ""), ("G5", "External", "validity"),
             ("G6", "Selection", ""), ("G7", "Dynamics", ""), ("G8", "Reporting", "")]
    gw, gap = 51, 7.0
    x0 = (W - (len(gates) * gw + (len(gates) - 1) * gap)) / 2
    gy, gh = 106, 48
    for i, (gid, l1, l2) in enumerate(gates):
        x = x0 + i * (gw + gap)
        hot = gid in ("G4", "G6")
        s.box(x, gy, gw, gh, fill=BLUE_F if hot else FILL,
              stroke=BLUE if hot else RULE, sw=1.0 if hot else 0.7)
        s.text(x + gw / 2, gy + 15, gid, size=9.5, weight="bold",
               fill=BLUE if hot else INK)
        s.text(x + gw / 2, gy + 27, l1, size=7.4, fill=SOFT)
        if l2:
            s.text(x + gw / 2, gy + 36, l2, size=7.4, fill=SOFT)
        if i:
            s.arrow(x - gap + 0.5, gy + gh / 2, x - 1.5, gy + gh / 2, sw=0.65)
        s.line(x + gw / 2, gy + gh, x + gw / 2, 176, dash="2 2.5")

    s.text(x0 + gw / 2 + 5 * (gw + gap), gy + gh + 11,
           "G4 and G6: the gates the natural-product arm passed in this study",
           size=7.2, fill=BLUE, style="italic")

    # registry band
    s.box(18, 178, W - 36, 34, fill=VIOLET_F, stroke=VIOLET, sw=0.8)
    s.text(W / 2, 192, "Artifact registry", size=9, weight="bold", fill=VIOLET)
    s.text(W / 2, 204,
           "append-only: every file read, every artifact written with its SHA-256, "
           "every contract outcome, every human intervention", size=7.4, fill=SOFT)

    s.line(18, 226, W - 18, 226)
    s.text(W / 2, 238,
           "A check that cannot be evaluated reports failure, never PASS. "
           "Figure 2 shows what happens inside one gate.",
           size=7.6, fill=SOFT, style="italic")

    p = out / "figure1_architecture_2026-08-21.svg"
    s.save(p)
    return p


# ---------------------------------------------------------------- Figure 2
def figure2(out: Path) -> Path:
    s = Svg(W, 268,
            "Inside one gate: contract, evaluator, monitor, policy",
            "A frozen YAML contract and a learned risk monitor feed a policy that "
            "resolves to one of four states. A hard-contract failure cannot be "
            "overridden by a low risk score.")

    s.text(W / 2, 16, "Inside one gate", size=11.5, weight="bold")
    s.text(W / 2, 29,
           "The deterministic contract decides first; the learned monitor can raise the "
           "state but never lower it", size=8, fill=SOFT)

    # inputs
    s.box(14, 48, 146, 46, fill="#FFFFFF", stroke=INK, sw=0.9)
    s.text(87, 63, "Frozen contract", size=9, weight="bold")
    s.lines(87, 75, ["YAML, version-controlled;", "checks and thresholds pre-registered"],
            size=7.4)

    s.box(14, 106, 146, 46, fill="#FFFFFF", stroke=RULE, sw=0.7)
    s.text(87, 121, "Learned risk monitor", size=9, weight="bold")
    s.lines(87, 133, ["calibrated on pre-release", "workflow telemetry only"], size=7.4)

    s.box(14, 164, 146, 40, fill=FILL, stroke=RULE, sw=0.7)
    s.text(87, 178, "Independent evaluator", size=9, weight="bold")
    s.text(87, 190, "deterministic; the agent cannot override it", size=7.4, fill=SOFT)

    s.arrow(160, 71, 196, 108)
    s.arrow(160, 129, 196, 124)
    s.arrow(160, 184, 196, 140)

    # policy
    s.box(198, 96, 100, 60, fill=BLUE_F, stroke=BLUE, sw=1.0)
    s.text(248, 114, "Policy", size=9.5, weight="bold", fill=BLUE)
    s.lines(248, 127, ["t_escalate = 0.2857", "t_halt = 0.8958", "frozen before evaluation"],
            size=7.3)

    # four states
    states = [("PASS", MOSS, MOSS_F, "continue to the next gate"),
              ("AUTO_REPAIR", AMBER, AMBER_F, "pre-specified fix, continue"),
              ("ESCALATE", AMBER, AMBER_F, "human review before release"),
              ("HALT", ALERT, ALERT_F, "workflow stops; reason logged")]
    sx, sw_, sh = 328, 78, 30
    for i, (name, col, fil, note) in enumerate(states):
        y = 62 + i * 42
        s.box(sx, y, sw_, sh, fill=fil, stroke=col, sw=0.9, r=13)
        s.text(sx + sw_ / 2, y + 19, name, size=8.6, weight="bold", fill=col)
        s.arrow(298, 126, sx - 2, y + sh / 2, sw=0.65)
        s.text(sx + sw_ + 7, y + 19, note, size=7.2, fill=SOFT, anchor="start")

    # fail-closed rule
    s.box(18, 218, W - 36, 34, fill="#FFFFFF", stroke=ALERT, sw=0.9)
    s.text(W / 2, 232, "The rule that makes the architecture fail closed", size=8.6,
           weight="bold", fill=ALERT)
    s.text(W / 2, 244,
           "A hard-contract failure overrides any risk score, and a check that cannot be "
           "evaluated returns failure rather than passing unverified.",
           size=7.5, fill=SOFT)

    p = out / "figure2_inside_a_gate_2026-08-21.svg"
    s.save(p)
    return p



# ------------------------------------------------------- Graphical abstract
def graphical_abstract(out: Path) -> Path:
    """Elsevier graphical abstract: landscape, and readable at small size.

    Not routed through an image-generation model, deliberately. This abstract
    carries measured results -- coverage fractions, a detection rate, compound
    counts -- and generated imagery must not depict real results. It is a
    data-bearing schematic, so it is drawn the way Figures 1 and 2 are: in
    points, with a hard floor on type size, from numbers that trace to
    artifacts.

    The previous version was 1560 units wide with 8.5-unit labels, which is
    2.0 pt at the size a graphical abstract is actually displayed, and its
    1.56:1 aspect was closer to square than the landscape format expected.

    Every string here is length-checked against its column before drawing, so
    a later edit cannot silently push text outside its box the way the first
    draft of this figure did.
    """
    GW, GH = 511.0, 206.0
    s = Svg(GW, GH,
            "Agentic virtual screening against the adenosine A2A receptor",
            "Three stages: a natural-product library that standard "
            "applicability-domain gates reject wholesale; an agentic pipeline "
            "whose code is verified under executable contracts; and a "
            "confidence-tiered candidate set.")

    # Short form: the full title does not fit 511 pt at 11 pt. The complete
    # title is carried by the SVG <title> element and by the caption.
    s.text(GW / 2, 15, "Verifying agent-generated code in A2A virtual screening",
           size=11, weight="bold")

    colw, gap = 155.0, 20.0
    x0 = (GW - (3 * colw + 2 * gap)) / 2
    top, ch = 26.0, 146.0
    inner = colw - 18

    def fits(txt: str, size: float, width: float = inner) -> str:
        # Helvetica averages ~0.5 em per character across mixed case; 0.52 is a
        # deliberately conservative estimate so the check errs toward catching
        # overflow rather than missing it.
        if len(txt) * size * 0.52 > width:
            raise SystemExit(
                f"graphical abstract: {txt!r} needs "
                f"~{len(txt) * size * 0.52:.0f} pt at {size} pt but the column "
                f"gives {width:.0f} pt")
        return txt

    cols = [
        (BLUE, BLUE_F, "1  Standard QC rejects the library",
         "236", ["Thai orchid natural products", "screened against A2AR"],
         ["0 % in-domain in ECFP4, FCFP4,", "topological torsion, Gobbi 2D",
          "86.3 % in physicochemical space"],
         "The same library, the same rule.",
         ["A similarity gate cannot tell", "unfamiliar from unreliable."]),
        (MOSS, MOSS_F, "2  The code itself is verified",
         "93.2 %", ["of 133 seeded defects caught;", "0 % when run ungoverned"],
         ["9 gates · 10 contracts · 94 checks", "false alarms on clean runs: 0.0",
          "learned monitor alone: 0.8 %"],
         "All 94 checks passed.",
         ["One metamorphic invariant still", "found a defect none encode."]),
        (AMBER, AMBER_F, "3  Corrected gates, tiered output",
         "3", ["convergent flavones, all from", "Pongamia pinnata"],
         ["25,604 → 500 full redocks", "structure- and ligand-based",
          "arms agree; p = 0.002"],
         "Computational hypotheses.",
         ["No wet-laboratory validation", "was performed."]),
    ]

    for i, (col, fil, head, big, note, bullets, kicker, caveat) in enumerate(cols):
        x = x0 + i * (colw + gap)
        s.box(x, top, colw, ch, fill="#FFFFFF", stroke=col, sw=0.9)
        s.box(x, top, colw, 15, fill=fil, stroke=col, sw=0.9)
        s.text(x + 8, top + 10.5, fits(head, 7.4), size=7.4, weight="bold",
               fill=col, anchor="start")

        s.text(x + colw / 2, top + 42, big, size=19, weight="bold", fill=col)
        for k, t in enumerate(note):
            # Species names are set in italic; nothing else in the note is.
            ital = t.strip() in ("Pongamia pinnata",)
            s.text(x + colw / 2, top + 53 + k * 9.2, fits(t, 7.2), size=7.2,
                   fill=SOFT, style="italic" if ital else "normal")

        s.line(x + 12, top + 70, x + colw - 12, top + 70)
        s.lines(x + colw / 2, top + 81, [fits(t, 7.2) for t in bullets],
                size=7.2, lead=1.42)

        s.line(x + 12, top + 108, x + colw - 12, top + 108)
        s.text(x + colw / 2, top + 119, fits(kicker, 7.4), size=7.4,
               weight="bold", fill=INK)
        # Two lines of caveat have to clear the box floor, descenders included.
        s.lines(x + colw / 2, top + 128, [fits(t, 7.0) for t in caveat],
                size=7.0, fill=SOFT, style="italic")

        if i < 2:
            s.arrow(x + colw + 3, top + ch / 2, x + colw + gap - 3,
                    top + ch / 2, sw=0.9)

    s.text(GW / 2, 186,
           "An applicability-domain verdict is a property of the library and the "
           "descriptor space, not of the library alone.", size=7.6, fill=INK)
    s.text(GW / 2, 197,
           "Contracts encode what their authors thought to check — and here the "
           "authors were the agent whose work was being checked.",
           size=7.2, fill=SOFT, style="italic")

    p = out / "figure_graphical_abstract_2026-08-21.svg"
    s.save(p)
    return p


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    root = next((c for c in (Path(r"D:\Agentic VS Workflow"),
                             Path("/mnt/d/Agentic VS Workflow"))
                 if (c / "00_project_control/08_DECISION_LOG.md").is_file()), None)
    if root is None:
        raise SystemExit("project root not found")
    out = Path(args.out) if args.out else root / "30_manuscript/figures"
    out.mkdir(parents=True, exist_ok=True)

    for p in (figure1(out), figure2(out), graphical_abstract(out)):
        txt = p.read_text(encoding="utf-8")
        sizes = sorted({float(v) for v in
                        __import__("re").findall(r'font-size="([\d.]+)"', txt)})
        print(f"{p.name}")
        print(f"   canvas 190 mm wide; font sizes in pt: {sizes}")
        print(f"   smallest {min(sizes)} pt  (floor {MIN_PT} pt)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

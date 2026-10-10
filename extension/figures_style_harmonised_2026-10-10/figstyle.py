"""One style for every figure of the manuscript (ADR-AIIM-367).

The figures were made at different times by different scripts. Before this module they differed in font family (DejaVu Sans
against Helvetica/Arial), in the size of the same kind of text (6.0 to 9.2 pt), in the panel label (A, (A), position and
size), in line widths, in the colour of the same kind of mark (three reds, two blues for Vina) and in width on the page
(7.5 to 8.0 inches). Everything below is the single place those choices are made. Scripts import it, call apply(), take
their colours from COL and their widths from LW, label panels with panel_label() and save with save_figure().

Rules (all figures are 190 mm = 7.48 inch wide, the journal's full-width size, so a size in points here is the size on the page):
  font        Arial (the family of Figures 1 and 2); nothing below 7.0 pt
  text        panel labels 9.0 pt bold, uppercase letter without brackets, 4 pt above the top-left corner of the panel;
              axis labels 7.5; tick labels, legends, annotations 7.0
  lines       axes and ticks 0.7; reference lines (diagonal, chance, median) 0.8 dashed grey; threshold lines 1.2 dashed;
              main data 1.6; secondary data 1.0; time series 0.8; bars edged 0.5
  colour      one hue per kind of object, used for that object in every figure
"""
from __future__ import annotations

import io
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
from PIL import Image             # noqa: E402

WIDTH_IN = 7.48          # 190 mm
DPI = 600
FS = {"panel": 9.0, "label": 7.5, "tick": 7.0, "legend": 7.0, "note": 7.0, "title": 7.5}
LW = {"axis": 0.7, "ref": 0.8, "thr": 1.2, "main": 1.6, "sec": 1.0, "series": 0.8, "bar_edge": 0.5}
COL = {
    "ink": "#222222", "axis": "#333333", "grey": "#999999", "light": "#a9b3c4", "dark": "#2b3a55",
    "alert": "#9b2226",                       # thresholds, "released" values, and the second engine (DiffDock, IFP) wherever a red is needed
    "vina": "#2c5f8a",                        # Vina and structural spaces (the BLUE of Figures 1 and 2)
    "green": "#7a9a3c", "amber": "#b5651d",   # pharmacophoric and physicochemical spaces
    "ref_bar": "#c9c9c9",                     # the reference (unchanged) bar in paired bar charts
    "hist": "#6b7a8f",                        # histogram bodies
    "control": "#111111",
    # compounds that appear in more than one figure keep their colour
    "NPLIB01185": "#b5651d", "NPLIB01643": "#2c5f8a", "NPLIB03056": "#7a9a3c",
    "SNP00157": "#8c4f9f", "SNP00163": "#c2487a",
}
REF_DASH = (0, (3, 3))
THR_DASH = (0, (4, 2))


def apply():
    plt.rcParams.update({
        "font.family": "Arial", "font.size": FS["tick"],
        "axes.labelsize": FS["label"], "xtick.labelsize": FS["tick"], "ytick.labelsize": FS["tick"],
        "legend.fontsize": FS["legend"], "axes.titlesize": FS["title"],
        "axes.linewidth": LW["axis"], "axes.edgecolor": COL["axis"],
        "xtick.major.width": LW["axis"], "ytick.major.width": LW["axis"],
        "xtick.major.size": 3.0, "ytick.major.size": 3.0,
        "axes.spines.top": False, "axes.spines.right": False,
        "lines.linewidth": LW["main"], "patch.linewidth": LW["bar_edge"],
        "savefig.dpi": DPI, "figure.dpi": 150, "savefig.facecolor": "white",
        "svg.fonttype": "path", "pdf.fonttype": 42,
    })


def panel_label(ax, letter, offset=4):
    """Panel label: bold uppercase letter, `offset` pt (4 as standard) above the top-left corner of the axes."""
    ax.annotate(letter, xy=(0, 1), xycoords="axes fraction", xytext=(0, offset), textcoords="offset points",
                fontsize=FS["panel"], weight="bold", ha="left", va="bottom", annotation_clip=False)


def panel_label_fig(fig, letter, x_in, y_in_from_top, fig_h_in):
    """Panel label for a panel that is not an axes (an image): position in inches from the top-left of the figure."""
    fig.text(x_in / fig.get_figwidth(), 1 - y_in_from_top / fig_h_in, letter, fontsize=FS["panel"], weight="bold", ha="left", va="top")


def save_figure(fig, stem, **kw):
    """Save PNG (600 dpi) and SVG at exactly 190 mm wide: render tight, then scale to the standard width."""
    stem = Path(stem)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=DPI, bbox_inches="tight", pad_inches=0.04, **kw)
    buf.seek(0)
    im = Image.open(buf).convert("RGB")
    h = round(im.height * WIDTH_IN * DPI / im.width)
    im = im.resize((round(WIDTH_IN * DPI), h), Image.LANCZOS)
    im.save(stem.with_suffix(".png"), dpi=(DPI, DPI))
    sbuf = io.BytesIO()
    fig.savefig(sbuf, format="svg", bbox_inches="tight", pad_inches=0.04, **kw)
    svg = sbuf.getvalue().decode("utf-8")
    m = re.search(r'<svg[^>]*?width="([\d.]+)pt"[^>]*?height="([\d.]+)pt"', svg)
    if m:
        w_pt, h_pt = float(m.group(1)), float(m.group(2))
        k = WIDTH_IN * 72 / w_pt
        svg = svg.replace(f'width="{m.group(1)}pt"', f'width="{w_pt * k:.2f}pt"', 1).replace(f'height="{m.group(2)}pt"', f'height="{h_pt * k:.2f}pt"', 1)
    stem.with_suffix(".svg").write_text(svg, encoding="utf-8")
    return stem.with_suffix(".png"), (WIDTH_IN, h / DPI)

"""Builds the harmonised copies of the figure scripts (ADR-AIIM-367) from their originals, which are not modified.

The originals live in 10_new_system/scripts (frozen tree) and in 20_new_runs/figures_docking_md_2026-10-09/analysis.
Each copy is the original with mechanical substitutions (the style module, panel labels, saving, fonts, the one red and the one
blue) and a few per-figure edits listed below. Run once; the copies are then run by hand.

    python make_harmonised.py
"""
from __future__ import annotations

import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCRIPTS = Path(r"D:\Agentic VS Workflow\10_new_system\scripts")
FIGRUN = Path(r"D:\Agentic VS Workflow\20_new_runs\figures_docking_md_2026-10-09\analysis")


def fs_size(m):
    v = float(m.group(1))
    key = "tick" if v < 7.25 else "label" if v < 8.8 else "panel"
    return f'fontsize=fs.FS["{key}"]'


def common(src: str) -> str:
    path_line = "import sys\nsys.path.insert(0, r'D:\\Agentic VS Workflow\\10_new_system\\scripts')"
    src = src.replace("from __future__ import annotations", "from __future__ import annotations\n" + path_line, 1)
    # style module
    src = re.sub(r"plt\.rcParams\.update\(\{.*?\}\)\n", "fs.apply()\n", src, count=1, flags=re.S)
    src = src.replace("import matplotlib.pyplot as plt", "import matplotlib.pyplot as plt\nimport figstyle as fs", 1)
    # panel labels written with set_title
    src = re.sub(r'(\w+)\.set_title\(\s*"([A-E])",\s*loc="left",\s*fontsize=[\d.]+,\s*weight="bold",\s*pad=[^)]*\)', r'fs.panel_label(\1, "\2")', src)
    # one red and one blue
    for red in ("#c1272d", "#b03a2e", "#b00020", "#c0392b", "#9b2226", "#d6403f"):
        src = src.replace(f'"{red}"', 'fs.COL["alert"]')
    for blue in ("#1f3b73", "#1f5fa8"):
        src = src.replace(f'"{blue}"', 'fs.COL["vina"]')
    src = re.sub(r"labelsize=(\d+(?:\.\d+)?)", 'labelsize=fs.FS["tick"]', src)
    src = src.replace('edgecolor="white", linewidth=0.5', 'edgecolor="white", linewidth=fs.LW["bar_edge"]').replace('edgecolor="white", linewidth=0.4', 'edgecolor="white", linewidth=fs.LW["bar_edge"]')
    # numeric font sizes -> the standard sizes
    src = re.sub(r"fontsize=(\d+(?:\.\d+)?)", fs_size, src)
    return src


def write(name: str, src: str):
    (HERE / name).write_text(src, encoding="utf-8")
    print("written", name)


# ---------------------------------------------------------------- Figure 3 (screening performance)
s = (SCRIPTS / "figure_lbdd_sbdd_performance.py").read_text(encoding="utf-8")
s = common(s)
s = s.replace('("Vina alone (baseline)",        "vina_score",            True,  fs.COL["vina"], 2.0)', '("Vina alone (baseline)",        "vina_score",            True,  fs.COL["vina"], fs.LW["main"])')
s = s.replace('True,  "#4a7fb5", 1.3)', 'True,  "#4a7fb5", fs.LW["sec"])').replace('False, "#7aa9d6", 1.3)', 'False, "#7aa9d6", fs.LW["sec"])').replace('False, "#b0c4de", 1.3)', 'False, "#b0c4de", fs.LW["sec"])')
s = s.replace('False, fs.COL["alert"], 2.0)', 'False, fs.COL["alert"], fs.LW["main"])')
s = s.replace("lw=0.8, ls=(0, (3, 3))", 'lw=fs.LW["ref"], ls=fs.REF_DASH')
s = s.replace('color=fs.COL["alert"], lw=1.8)', 'color=fs.COL["alert"], lw=fs.LW["main"])')
s = s.replace('color=fs.COL["alert"], lw=0.7, ls=(0, (2, 2)), alpha=0.7)', 'color=fs.COL["alert"], lw=fs.LW["thr"], ls=fs.THR_DASH)')
s = s.replace("colour, lw=2.0,", 'colour, lw=fs.LW["main"],')
s = re.sub(r'stem = out / "figure4_screening_performance_2026-08-24"\n\s*fig\.savefig\(stem\.with_suffix\("\.png"\), bbox_inches="tight"\)\n\s*fig\.savefig\(stem\.with_suffix\("\.svg"\), bbox_inches="tight"\)',
           'stem = out / "figure3_screening_performance_harmonised"\n    fs.save_figure(fig, stem)', s)
s = re.sub(r"(axC\.legend\([^)]*?bbox_to_anchor=)\([0-9.]+, ?[0-9.]+\)", lambda m: m.group(1) + "(0.15, 0.97)", s, flags=re.S)
import fig3_redesign   # ADR-AIIM-368
s = fig3_redesign.apply(s)
write("fig3_screening.py", s)


# ---------------------------------------------------------------- Figure 4 (applicability domain across spaces)
s = (SCRIPTS / "figure_ad_multispace.py").read_text(encoding="utf-8")
s = common(s)
s = s.replace("labelsize=6.6", 'labelsize=fs.FS["tick"]')
s = s.replace("linewidth=0.6", 'linewidth=fs.LW["bar_edge"]').replace("b.set_linewidth(0.6)", 'b.set_linewidth(fs.LW["bar_edge"])')
s = s.replace('linewidth=1.3, zorder=5)', 'linewidth=fs.LW["thr"], linestyle=fs.THR_DASH, zorder=5)')
s = s.replace('axB.plot([], [], color=fs.COL["alert"], linewidth=1.3,', 'axB.plot([], [], color=fs.COL["alert"], linewidth=fs.LW["thr"], linestyle=fs.THR_DASH,')
s = s.replace("bbox_to_anchor=(0.0, 1.19)", "bbox_to_anchor=(0.045, 1.19)")
s = re.sub(r'stem = out / "figure_ad_multispace_2026-08-25"\n\s*fig\.savefig\(stem\.with_suffix\("\.png"\), bbox_inches="tight"\)\n\s*fig\.savefig\(stem\.with_suffix\("\.svg"\), bbox_inches="tight"\)',
           'stem = out / "figure4_ad_multispace_harmonised"\n    fs.save_figure(fig, stem)', s)
# never overwrite the original run's provenance file (G8 reads it): the harmonised script writes its own record into out/
s = s.replace('(under("20_new_runs/ad_sweep_selfsim_2026-08-25", "figure_ad_multispace_provenance.json")).write_text(',
              '(out / "figure4_ad_multispace_harmonised_provenance.json").write_text(')
write("fig4_ad.py", s)


# ---------------------------------------------------------------- Figure 5 (permutation distribution)
s = (SCRIPTS / "figure_permutation_distribution.py").read_text(encoding="utf-8")
s = common(s)
s = re.sub(r'ax\.text\(-?[\d.]+, [\d.]+, "\(([A-C])\)", transform=ax\.transAxes, fontweight="bold",\s*fontsize=[^)]*\)', r'fs.panel_label(ax, "\1")', s)
s = re.sub(r'ax\.text\(-?[\d.]+, [\d.]+, "\(([A-C])\)", transform=ax\.transAxes, fontweight="bold",\s*fontsize=fs\.FS\["panel"\], va="top"\)', r'fs.panel_label(ax, "\1")', s)
s = s.replace('lw=1.5, ls="--"', 'lw=fs.LW["thr"], ls=fs.THR_DASH').replace('lw=1.0, ls=":"', 'lw=fs.LW["ref"], ls=":"')
s = s.replace('SRC = resolve("30_manuscript/figures")', 'SRC = Path(__file__).resolve().parent / "out"').replace('EMBED = resolve("30_manuscript/CURRENT_v5/figures_png")', 'EMBED = Path(__file__).resolve().parent / "out"')
s = s.replace('STEM = "figure7_permutation_distribution_2026-10-02"', 'STEM = "figure5_permutation_harmonised"')
s = s.replace('    out = SRC / f"{STEM}.png"\n    fig.savefig(out)', '    out, _ = fs.save_figure(fig, SRC / STEM)')
s = s.replace('    embed = EMBED / f"{STEM}.png"\n    embed.write_bytes(out.read_bytes())', '    embed = out')
write("fig5_permutation.py", s)

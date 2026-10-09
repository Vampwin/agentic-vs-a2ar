"""Figure 6 (orthogonal pose comparison) with the fifteen focus-library survivors added.

The three structure panels are PyMOL renderings made together by render_snplib_pymol.py with one camera orientation
(set on the crystallographic ZM241385 pose, orthographic, zoom per panel): A the control, B SNP00157, C SNP00163, each with
the crystallographic pose as a thin grey reference. Panels D and E are drawn from pose_comparison.csv of ADR-AIIM-327,
so that the figure carries the result that Section 3.8 reports in prose.

Legend (manuscript-figures-tables skill, phase 2):

  The figure legend was written before the artwork. It is the legend of the manuscript's figure and is not
  reproduced here, because the manuscript's text is not deposited.

Anchors: every value plotted equals diffdock_snplib.json / pose_comparison.csv, and the counts the legend
rests on (none of fifteen within 2 A, eleven within 8 A for DiffDock against seven for Vina, five agreeing
classifications) are recomputed from the plotted points. The structure panels are recorded by hash.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg   # noqa: E402
import matplotlib.pyplot as plt    # noqa: E402
import numpy as np                 # noqa: E402
import pandas as pd                # noqa: E402

ROOT = Path(r"D:\Agentic VS Workflow")
RUNS = ROOT / "20_new_runs"
DD = RUNS / "diffdock_snplib_2026-10-07"
PKG = ROOT / "30_manuscript" / "AIIM_submission_package_v8_working_2026-10-03"
OUT = RUNS / "figures_docking_md_2026-10-09"
PA = OUT / "source" / "_panel_control.png"
PB = OUT / "source" / "_panel_SNP00157.png"
PC = OUT / "source" / "_panel_SNP00163.png"
CTRL = RUNS / "diffdock_orthogonal_2026-08-23" / "orthogonal_check_result.json"
AGREE_A, SITE_A = 2.0, 8.0
DARK, LIGHT = "#2b3a55", "#a9b3c4"

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.2, "axes.linewidth": 0.7,
                     "axes.edgecolor": "#333333", "xtick.major.width": 0.7, "ytick.major.width": 0.7,
                     "savefig.dpi": 600, "figure.dpi": 150})


def crop_content(path, pad=40):
    """Trim the transparent margin of a PyMOL render (alpha > 0 or non-white) so each panel uses its box."""
    im = mpimg.imread(path)
    mask = (im[..., 3] > 0.02) if im.shape[2] == 4 else (im[..., :3].min(axis=2) < 0.98)
    ys, xs = np.where(mask)
    y0, y1, x0, x1 = max(ys.min() - pad, 0), min(ys.max() + pad, im.shape[0]), max(xs.min() - pad, 0), min(xs.max() + pad, im.shape[1])
    return im[y0:y1, x0:x1]


def main():
    j = json.loads((DD / "diffdock_snplib.json").read_text("utf-8"))
    df = pd.read_csv(DD / "pose_comparison.csv")
    assert len(df) == 15
    # ---- anchors from the plotted points
    n_agree = int((df["rmsd_dd_vs_vina_A"] <= AGREE_A).sum())
    dd_site = int((df["dd_centroid_dist_A"] <= SITE_A).sum())
    vina_site = int((df["vina_centroid_dist_recomputed_A"] <= SITE_A).sum())
    same = int(((df["dd_centroid_dist_A"] <= SITE_A) == (df["vina_centroid_dist_recomputed_A"] <= SITE_A)).sum())
    g = j["grading"]; d = j["descriptive_not_graded"]
    checks = {"none_within_2A": n_agree == 0 == g["E1_fewer_than_8_of_15_within_2A"]["n_within_2A"],
              "dd_on_pocket_11": dd_site == d["dd_on_pocket_n"] == 11,
              "vina_on_pocket_7": vina_site == d["vina_on_pocket_n"] == 7,
              "classifications_agree_5": same == d["on_off_pocket_classification_agrees_n"] == 5,
              "min_max_rmsd": [round(float(df["rmsd_dd_vs_vina_A"].min()), 3), round(float(df["rmsd_dd_vs_vina_A"].max()), 3)] == d["min_max_rmsd_A"],
              "panel_hashes": {q.name: hashlib.sha256(q.read_bytes()).hexdigest() for q in (PA, PB, PC)}}
    assert all(v for k, v in checks.items() if k != "panel_hashes") and j["anchors"]["1_frame"]["verdict"] == "MET", checks
    print("anchors met:", {k: v for k, v in checks.items() if k != "panel_hashes"})

    ctrl = json.loads(CTRL.read_text("utf-8"))["control"]["rmsd_to_crystal_A"]
    rows = {q: df[df["SNP_ID"] == q].iloc[0] for q in ("SNP00157", "SNP00163")}
    W = 7.48
    top_h = 2.85
    bot_h = 3.05
    H = top_h + bot_h + 0.15
    fig = plt.figure(figsize=(W, H))
    titles = [("A", "ZM241385 (control)", f"DiffDock {ctrl:.2f} Å from the crystal pose"),
              ("B", "SNP00157", f"DiffDock {rows['SNP00157']['rmsd_dd_vs_vina_A']:.1f} Å from Vina"),
              ("C", "SNP00163", f"DiffDock {rows['SNP00163']['rmsd_dd_vs_vina_A']:.1f} Å from Vina")]
    for i, (path, (lab, name, sub)) in enumerate(zip((PA, PB, PC), titles)):
        x0 = 0.005 + i * 0.3317
        ax = fig.add_axes([x0, (bot_h + 0.15) / H + 0.30 / H, 0.325, (top_h - 0.62) / H])
        ax.imshow(crop_content(path)); ax.axis("off"); ax.set_anchor("C")
        fig.text(x0 + 0.004, (H - 0.04) / H, lab, fontsize=8.4, weight="bold", va="top")
        fig.text(x0 + 0.03, (H - 0.04) / H, name, fontsize=7.4, weight="bold", va="top")
        fig.text(x0 + 0.03, (H - 0.17) / H, sub, fontsize=6.6, va="top", color="#444444")
    fig.legend(handles=[plt.Line2D([0], [0], color="#7a7a7a", lw=1.3, label="crystallographic ZM241385 (thin = reference)"),
                        plt.Line2D([0], [0], color="#1f5fa8", lw=2.6, label="Vina rank-1 pose"),
                        plt.Line2D([0], [0], color="#d6403f", lw=2.6, label="DiffDock rank-1 pose")],
               loc="lower center", ncol=3, frameon=False, fontsize=7.0, bbox_to_anchor=(0.5, (bot_h + 0.15) / H),
               handlelength=1.8, columnspacing=1.4, handletextpad=0.6)
    hb = bot_h / H
    axC = fig.add_axes([0.14, 0.085, 0.36, hb - 0.075])
    axD = fig.add_axes([0.62, 0.085, 0.355, hb - 0.075])

    s = df.sort_values("rmsd_dd_vs_vina_A", ascending=True).reset_index(drop=True)
    y = np.arange(len(s))
    cols = [DARK if b else LIGHT for b in s["in_shortlist"]]
    axC.barh(y, s["rmsd_dd_vs_vina_A"], color=cols, height=0.72, lw=0)
    axC.axvline(AGREE_A, color="#c0392b", lw=1.0, dashes=(4, 2), zorder=3)
    axC.set_yticks(y); axC.set_yticklabels(s["SNP_ID"], fontsize=6.4)
    axC.set_xlim(0, 36); axC.set_ylim(-0.7, len(s) - 0.3); axC.tick_params(axis="x", labelsize=7)
    axC.set_xlabel("DiffDock to Vina pose RMSD (\u00c5)", fontsize=7.6)
    axC.invert_yaxis()
    for sp in ("top", "right"):
        axC.spines[sp].set_visible(False)
    axC.annotate("2.0 Å", xy=(AGREE_A, 1.0), xycoords=("data", "axes fraction"), xytext=(2, 1), textcoords="offset points",
                 color="#c0392b", fontsize=6.4, ha="left", va="bottom")
    axC.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=DARK, label="shortlisted (8)"),
                        plt.Rectangle((0, 0), 1, 1, color=LIGHT, label="other survivors (7)")],
               fontsize=6.4, frameon=False, loc="upper right", handlelength=1.0)

    for _, r in df.iterrows():
        axD.plot(r["vina_centroid_dist_recomputed_A"], r["dd_centroid_dist_A"], "o", ms=4.8,
                 mfc=DARK if r["in_shortlist"] else "white", mec=DARK, mew=0.9, zorder=3)
    for sid in ("SNP00157", "SNP00163", "SNP00006"):
        r = df[df["SNP_ID"] == sid].iloc[0]
        axD.annotate(sid[-3:], (r["vina_centroid_dist_recomputed_A"], r["dd_centroid_dist_A"]), xytext=(4, 3),
                     textcoords="offset points", fontsize=6.4, color=DARK)
    axD.plot([0, 32], [0, 32], color="#999999", lw=0.7, zorder=1)
    axD.axvline(SITE_A, color="#666666", lw=0.8, dashes=(1.5, 2), zorder=1)
    axD.axhline(SITE_A, color="#666666", lw=0.8, dashes=(1.5, 2), zorder=1)
    axD.set_xlim(0, 32); axD.set_ylim(0, 32); axD.set_aspect("equal", adjustable="box")
    axD.tick_params(labelsize=7)
    axD.set_xlabel("Vina pose: distance from the ZM241385 centroid (\u00c5)", fontsize=7.2)
    axD.set_ylabel("DiffDock pose: distance (\u00c5)", fontsize=7.2)
    for sp in ("top", "right"):
        axD.spines[sp].set_visible(False)
    axD.legend(handles=[plt.Line2D([], [], marker="o", ls="", mfc=DARK, mec=DARK, ms=4.6, label="shortlisted"),
                        plt.Line2D([], [], marker="o", ls="", mfc="white", mec=DARK, ms=4.6, label="other survivors")],
               fontsize=6.4, frameon=False, loc="upper right", borderaxespad=0.1)
    fig.text(0.012, hb - 0.012, "D", fontsize=8.4, weight="bold", va="top")
    fig.text(0.545, hb - 0.012, "E", fontsize=8.4, weight="bold", va="top")
    stem = OUT / "Figure_6_docking_snplib"
    fig.savefig(stem.with_suffix(".png"), bbox_inches="tight"); fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    prov = {"figure": "Figure 6", "script": "figure6_docking_snplib.py", "adr": "ADR-AIIM-327",
            "structure_panels": {"A": "30_manuscript/figures/_panelA_control.png, placed unchanged", "B_C": "source/_panel_SNP00157.png and _panel_SNP00163.png, from render_snplib_pymol.py"}, "checks": checks,
            "inputs_sha256": {str((DD / f).relative_to(ROOT)): hashlib.sha256((DD / f).read_bytes()).hexdigest()
                              for f in ("diffdock_snplib.json", "pose_comparison.csv")}}
    (OUT / "Figure_6_docking_snplib_provenance.json").write_text(json.dumps(prov, indent=1), "utf-8")
    print("written", stem.with_suffix(".png"))


if __name__ == "__main__":
    main()

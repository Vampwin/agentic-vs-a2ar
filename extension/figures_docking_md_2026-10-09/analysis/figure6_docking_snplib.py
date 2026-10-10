"""Supplementary Figure S1 (formerly Figure 6): DiffDock-L run blind against the protein alone, compared with Vina in the 22 A
site box, for the ZM241385 control, SNP00157 and SNP00163 and for the fifteen focus-library survivors. Harmonised with the other
figures (ADR-AIIM-367); moved to the supplement in ADR-AIIM-366.

A-C  PyMOL renders from render_snplib_pymol.py: one camera orientation, a zoom per panel that fits the two poses, whole ribbon
     segments, the Asn253 label drawn here at the calibrated position; crystal pose thin grey, Vina blue, DiffDock-L red.
D    heavy-atom RMSD between the DiffDock-L and Vina rank-1 poses of the fifteen survivors (pose_comparison.csv, ADR-AIIM-327).
E    distance of each pose's centroid from the ZM241385 centroid, Vina horizontal, DiffDock-L vertical.
SNP00157 and SNP00163 keep, in D and E, the colours they have in Figure 7.

The legend of the manuscript's figure was written before the artwork and is not reproduced here (the manuscript's text is
not deposited). Anchors, all recomputed from the plotted points: none of fifteen within 2 A, eleven within 8 A for DiffDock-L
against seven for Vina, five agreeing classifications, the min and max RMSD of ADR-AIIM-327.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg   # noqa: E402
import matplotlib.pyplot as plt    # noqa: E402
import numpy as np                 # noqa: E402
import pandas as pd                # noqa: E402

sys.path.insert(0, "D:/Agentic VS Workflow/20_new_runs/figures_style_harmonised_2026-10-10")
import figstyle as fs              # noqa: E402

ROOT = Path(r"D:\Agentic VS Workflow")
RUNS = ROOT / "20_new_runs"
DD = RUNS / "diffdock_snplib_2026-10-07"
OUT = RUNS / "figures_docking_md_2026-10-09"
PANELS = [("A", "control", "ZM241385 (control)"), ("B", "SNP00157", "SNP00157"), ("C", "SNP00163", "SNP00163")]
AGREE_A, SITE_A = 2.0, 8.0
DARK, LIGHT = fs.COL["dark"], fs.COL["light"]
fs.apply()


def crop_content(path, pad=30):
    im = mpimg.imread(path)
    mask = (im[..., 3] > 0.02) if im.shape[2] == 4 else (im[..., :3].min(axis=2) < 0.98)
    ys, xs = np.where(mask)
    y0, y1 = max(ys.min() - pad, 0), min(ys.max() + pad, im.shape[0])
    x0, x1 = max(xs.min() - pad, 0), min(xs.max() + pad, im.shape[1])
    return im[y0:y1, x0:x1], (x0, y0)


def point_colour(c, shortlisted):
    if c in fs.COL:
        return fs.COL[c]
    return DARK if shortlisted else LIGHT


def main():
    j = json.loads((DD / "diffdock_snplib.json").read_text("utf-8"))
    df = pd.read_csv(DD / "pose_comparison.csv")
    assert len(df) == 15
    n_agree = int((df["rmsd_dd_vs_vina_A"] <= AGREE_A).sum())
    dd_site = int((df["dd_centroid_dist_A"] <= SITE_A).sum())
    vina_site = int((df["vina_centroid_dist_recomputed_A"] <= SITE_A).sum())
    same = int(((df["dd_centroid_dist_A"] <= SITE_A) == (df["vina_centroid_dist_recomputed_A"] <= SITE_A)).sum())
    g = j["grading"]; d = j["descriptive_not_graded"]
    checks = {"none_within_2A": n_agree == 0 == g["E1_fewer_than_8_of_15_within_2A"]["n_within_2A"],
              "dd_on_pocket_11": dd_site == d["dd_on_pocket_n"] == 11,
              "vina_on_pocket_7": vina_site == d["vina_on_pocket_n"] == 7,
              "classifications_agree_5": same == d["on_off_pocket_classification_agrees_n"] == 5,
              "min_max_rmsd": [round(float(df["rmsd_dd_vs_vina_A"].min()), 3), round(float(df["rmsd_dd_vs_vina_A"].max()), 3)] == d["min_max_rmsd_A"]}
    assert all(checks.values()) and j["anchors"]["1_frame"]["verdict"] == "MET", checks
    print("anchors met:", checks)

    W, top_h, bot_h = fs.WIDTH_IN, 2.75, 3.05
    H = top_h + bot_h + 0.15
    fig = plt.figure(figsize=(W, H))
    sub = {"control": "DiffDock-L 1.92 \u00c5 from the crystal pose",
           "SNP00157": f"DiffDock-L {df.set_index('SNP_ID').loc['SNP00157', 'rmsd_dd_vs_vina_A']:.1f} \u00c5 from Vina",
           "SNP00163": f"DiffDock-L {df.set_index('SNP_ID').loc['SNP00163', 'rmsd_dd_vs_vina_A']:.1f} \u00c5 from Vina"}
    hashes = {}
    for i, (lab, cid, name) in enumerate(PANELS):
        f = OUT / "source" / f"_panel_{cid}.png"
        hashes[f.name] = hashlib.sha256(f.read_bytes()).hexdigest()
        im, (cx0, cy0) = crop_content(f)
        x0 = 0.005 + i * 0.3317
        ax = fig.add_axes([x0, (bot_h + 0.15) / H + 0.30 / H, 0.325, (top_h - 0.62) / H])
        ax.imshow(im); ax.axis("off")
        cal = mpimg.imread(OUT / "source" / f"_calib_{cid}.png")
        yy, xx = np.where(cal[..., 3] > 0.5) if cal.shape[2] == 4 else np.where(cal[..., :3].min(axis=2) < 0.9)
        ax.annotate("Asn253", xy=(xx.mean() - cx0, yy.mean() - cy0), xytext=(10, 12), textcoords="offset points",
                    fontsize=fs.FS["note"], color="#44484f", ha="left", va="bottom")
        fs.panel_label_fig(fig, lab, x0 * W + 0.02, 0.05, H)
        fig.text(x0 + 0.032, (H - 0.05) / H, name, fontsize=fs.FS["label"], weight="bold", va="top")
        fig.text(x0 + 0.032, (H - 0.21) / H, sub[cid], fontsize=fs.FS["note"], va="top", color="#444444")
    fig.legend(handles=[plt.Line2D([0], [0], color="#7a7a7a", lw=fs.LW["sec"], label="crystallographic ZM241385 (thin = reference)"),
                        plt.Line2D([0], [0], color=fs.COL["vina"], lw=fs.LW["main"] * 1.6, label="Vina rank-1 pose, 22 \u00c5 site box"),
                        plt.Line2D([0], [0], color=fs.COL["alert"], lw=fs.LW["main"] * 1.6, label="DiffDock-L rank-1 pose, blind")],
               loc="lower center", ncol=3, frameon=False, fontsize=fs.FS["legend"], bbox_to_anchor=(0.5, (bot_h + 0.15) / H),
               handlelength=1.8, columnspacing=1.4, handletextpad=0.6)
    hb = bot_h / H
    axD = fig.add_axes([0.145, 0.085, 0.345, hb - 0.095])
    axE = fig.add_axes([0.625, 0.085, 0.350, hb - 0.095])

    s = df.sort_values("rmsd_dd_vs_vina_A", ascending=True).reset_index(drop=True)
    y = np.arange(len(s))
    axD.barh(y, s["rmsd_dd_vs_vina_A"], color=[point_colour(c, b) for c, b in zip(s["SNP_ID"], s["in_shortlist"])], height=0.72, lw=0)
    axD.axvline(AGREE_A, color=fs.COL["alert"], lw=fs.LW["thr"], ls=fs.THR_DASH, zorder=3)
    axD.set_yticks(y); axD.set_yticklabels(s["SNP_ID"])
    axD.set_xlim(0, 36); axD.set_ylim(-0.7, len(s) - 0.3); axD.invert_yaxis()
    axD.set_xlabel("DiffDock-L to Vina pose RMSD (\u00c5)")
    axD.annotate("2.0 \u00c5", xy=(AGREE_A, 1.0), xycoords=("data", "axes fraction"), xytext=(3, 1), textcoords="offset points",
                 color=fs.COL["alert"], fontsize=fs.FS["note"], ha="left", va="bottom")
    axD.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=DARK, label="shortlisted (8)"),
                        plt.Rectangle((0, 0), 1, 1, color=LIGHT, label="other survivors (7)")],
               frameon=False, loc="lower right", handlelength=1.0, bbox_to_anchor=(1.0, 0.12))
    fs.panel_label_fig(fig, "D", 0.03, top_h + 0.15 + 0.10, H)

    for _, r in df.iterrows():
        c = r["SNP_ID"]; col = point_colour(c, bool(r["in_shortlist"])); solid = bool(r["in_shortlist"]) or c in fs.COL
        axE.plot(r["vina_centroid_dist_recomputed_A"], r["dd_centroid_dist_A"], "o", ms=5.0,
                 mfc=col if solid else "white", mec=col if solid else DARK, mew=0.9, zorder=3)
    for sid in ("SNP00157", "SNP00163", "SNP00006"):
        r = df[df["SNP_ID"] == sid].iloc[0]
        axE.annotate(sid[-3:], (r["vina_centroid_dist_recomputed_A"], r["dd_centroid_dist_A"]), xytext=(4, 3),
                     textcoords="offset points", fontsize=fs.FS["note"], color=point_colour(sid, True))
    axE.plot([0, 32], [0, 32], color=fs.COL["grey"], lw=fs.LW["ref"], zorder=1)
    axE.axvline(SITE_A, color=fs.COL["grey"], lw=fs.LW["ref"], ls=fs.REF_DASH, zorder=1)
    axE.axhline(SITE_A, color=fs.COL["grey"], lw=fs.LW["ref"], ls=fs.REF_DASH, zorder=1)
    axE.set_xlim(0, 32); axE.set_ylim(0, 32); axE.set_aspect("equal", adjustable="box")
    axE.set_xlabel("Vina pose: distance from the\nZM241385 centroid (\u00c5)")
    axE.set_ylabel("DiffDock-L pose: distance (\u00c5)")
    axE.legend(handles=[plt.Line2D([], [], marker="o", ls="", mfc=DARK, mec=DARK, ms=5.0, label="shortlisted"),
                        plt.Line2D([], [], marker="o", ls="", mfc="white", mec=DARK, ms=5.0, label="other survivors")],
               frameon=False, loc="upper right", borderaxespad=0.1)
    fs.panel_label_fig(fig, "E", 0.545 * W, top_h + 0.15 + 0.10, H)
    stem = OUT / "Figure_6_docking_snplib"
    fs.save_figure(fig, stem)
    prov = {"figure": "Supplementary Figure S1 (formerly Figure 6), harmonised style", "script": "figure6_docking_snplib.py", "adr": "ADR-AIIM-327, ADR-AIIM-367",
            "checks": checks, "panel_png_sha256": hashes,
            "inputs_sha256": {str((DD / f).relative_to(ROOT)): hashlib.sha256((DD / f).read_bytes()).hexdigest()
                              for f in ("diffdock_snplib.json", "pose_comparison.csv")}}
    (OUT / "Figure_6_docking_snplib_provenance.json").write_text(json.dumps(prov, indent=1), "utf-8")
    print("written", stem.with_suffix(".png"))


if __name__ == "__main__":
    main()

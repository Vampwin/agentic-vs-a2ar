"""Figure 6, redesigned (ADR-AIIM-366) and harmonised with the other figures (ADR-AIIM-367).

A-C  PyMOL renders from render_pocket_pymol.py, one camera and one scale: the ZM241385 control, SNP00157, SNP00163;
     crystal pose thin grey, Vina rank-1 pose from the 22 A site box blue, DiffDock-Pocket rank-1 pose red. The ribbons are whole
     secondary-structure segments, so every ribbon is continuous; the Asn253 label is drawn here in the figure font at the
     position calibrated by the renderer.
D    heavy-atom RMSD between the DiffDock-Pocket and Vina (box) rank-1 poses of the fifteen survivors (filter_pocket_poses.csv).
E    distance of the Vina rank-1 pose from the ZM241385 centroid, site box (x) against the whole receptor (y) (matched_poses.csv).
SNP00157 and SNP00163 keep, in D and E, the colours they have in Figure 7.

The legend of the manuscript's figure was written before the artwork and is not reproduced here (the manuscript's text is
not deposited). Anchors, all asserted before the figure is written and recomputed from the plotted points: none of fifteen
below 2.0 A in D and its median; 7 of 15 on pocket in the box against 2 of 15 over the whole receptor; 13 of 15 whole-receptor
poses 11.9-13.3 A from the site; the control 0.35 A; the two compounds named in the panel titles.
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

sys.path.insert(0, r"D:\Agentic VS Workflow\20_new_runs\figures_style_harmonised_2026-10-10")
import figstyle as fs              # noqa: E402

ROOT = Path(r"D:\Agentic VS Workflow")
RUNS = ROOT / "20_new_runs"
OUT = RUNS / "figures_docking_md_2026-10-09"
FP = pd.read_csv(RUNS / "docking_matched_filter_pocket_2026-10-09" / "filter_pocket_poses.csv").set_index("compound")
MP = pd.read_csv(RUNS / "docking_matched_search_space_2026-10-09" / "matched_poses.csv").set_index("compound")
PANELS = [("A", "ZM241385", "ZM241385 (control)"), ("B", "SNP00157", "SNP00157"), ("C", "SNP00163", "SNP00163")]
AGREE_A, SITE_A = 2.0, 8.0
DARK, LIGHT = fs.COL["dark"], fs.COL["light"]
fs.apply()


def union_crop(files, pad=30):
    """Crop all panels with the same box (the union of their non-transparent pixels): the shared scale and alignment are kept."""
    ims = [mpimg.imread(f) for f in files]
    boxes = []
    for im in ims:
        mask = (im[..., 3] > 0.02) if im.shape[2] == 4 else (im[..., :3].min(axis=2) < 0.98)
        ys, xs = np.where(mask)
        boxes.append((ys.min(), ys.max(), xs.min(), xs.max()))
    y0 = max(min(b[0] for b in boxes) - pad, 0); y1 = min(max(b[1] for b in boxes) + pad, ims[0].shape[0])
    x0 = max(min(b[2] for b in boxes) - pad, 0); x1 = min(max(b[3] for b in boxes) + pad, ims[0].shape[1])
    return [im[y0:y1, x0:x1] for im in ims], (x0, y0)


def point_colour(c, shortlisted):
    """Compounds that appear in Figure 7 keep their colour; other shortlisted compounds are dark, the rest light."""
    if c in fs.COL:
        return fs.COL[c]
    return DARK if shortlisted else LIGHT


def main():
    surv = [c for c in FP.index if c != "ZM241385"]
    assert len(surv) == 15 and set(surv) == set(MP.index) - {"ZM241385"}
    rm = FP.loc[surv, "pocket_C_vs_vina_box_rmsd_A"]
    xb = MP.loc[surv, "vina_pocket_centroid_dist_A"]
    yb = MP.loc[surv, "vina_blind_centroid_dist_A"]
    ctrl = float(FP.loc["ZM241385", "pocket_C_rmsd_to_crystal_A"])
    shortlist = FP.loc[surv, "in_shortlist"].astype(bool)
    checks = {"none_below_2A": int((rm <= AGREE_A).sum()) == 0, "median_12.8": round(float(rm.median()), 1) == 12.8,
              "box_on_pocket_7": int((xb <= SITE_A).sum()) == 7, "blind_on_pocket_2": int((yb <= SITE_A).sum()) == 2,
              "blind_13_of_15_11.9_to_13.3": int(((yb >= 11.85) & (yb <= 13.3)).sum()) == 13,
              "control_0.35": round(ctrl, 2) == 0.35, "shortlist_8": int(shortlist.sum()) == 8,
              "titles": [round(float(FP.loc[c, "pocket_C_vs_vina_box_rmsd_A"]), 1) for c in ("SNP00157", "SNP00163")] == [9.9, 3.4]}
    assert all(checks.values()), checks
    print("anchors met:", checks)

    W, top_h, bot_h = fs.WIDTH_IN, 2.75, 3.05
    H = top_h + bot_h + 0.15
    fig = plt.figure(figsize=(W, H))
    sub = {"ZM241385": f"DiffDock-Pocket {ctrl:.2f} \u00c5 from the crystal pose",
           "SNP00157": f"DiffDock-Pocket {FP.loc['SNP00157', 'pocket_C_vs_vina_box_rmsd_A']:.1f} \u00c5 from Vina",
           "SNP00163": f"DiffDock-Pocket {FP.loc['SNP00163', 'pocket_C_vs_vina_box_rmsd_A']:.1f} \u00c5 from Vina"}
    hashes = {}
    files = [OUT / "source" / f"_pocket_panel_{cid}.png" for _, cid, _ in PANELS]
    cropped, (cx0, cy0) = union_crop(files)
    for i, ((lab, cid, name), f) in enumerate(zip(PANELS, files)):
        hashes[f.name] = hashlib.sha256(f.read_bytes()).hexdigest()
        x0 = 0.005 + i * 0.3317
        ax = fig.add_axes([x0, (bot_h + 0.15) / H + 0.30 / H, 0.325, (top_h - 0.62) / H])
        ax.imshow(cropped[i]); ax.axis("off")
        cal = mpimg.imread(OUT / "source" / f"_pocket_calib_{cid}.png")
        yy, xx = np.where(cal[..., 3] > 0.5) if cal.shape[2] == 4 else np.where(cal[..., :3].min(axis=2) < 0.9)
        ax.annotate("Asn253", xy=(xx.mean() - cx0, yy.mean() - cy0), xytext=(10, 12), textcoords="offset points",
                    fontsize=fs.FS["note"], color="#44484f", ha="left", va="bottom")
        fs.panel_label_fig(fig, lab, x0 * W + 0.02, 0.05, H)
        fig.text(x0 + 0.032, (H - 0.05) / H, name, fontsize=fs.FS["label"], weight="bold", va="top")
        fig.text(x0 + 0.032, (H - 0.21) / H, sub[cid], fontsize=fs.FS["note"], va="top", color="#444444")
    fig.legend(handles=[plt.Line2D([0], [0], color="#7a7a7a", lw=fs.LW["sec"], label="crystallographic ZM241385 (thin = reference)"),
                        plt.Line2D([0], [0], color=fs.COL["vina"], lw=fs.LW["main"] * 1.6, label="Vina rank-1 pose, 22 \u00c5 site box"),
                        plt.Line2D([0], [0], color=fs.COL["alert"], lw=fs.LW["main"] * 1.6, label="DiffDock-Pocket rank-1 pose")],
               loc="lower center", ncol=3, frameon=False, fontsize=fs.FS["legend"], bbox_to_anchor=(0.5, (bot_h + 0.15) / H),
               handlelength=1.8, columnspacing=1.4, handletextpad=0.6)
    hb = bot_h / H
    axD = fig.add_axes([0.145, 0.085, 0.345, hb - 0.095])
    axE = fig.add_axes([0.625, 0.085, 0.350, hb - 0.095])

    s = rm.sort_values().reset_index()
    y = np.arange(len(s))
    axD.barh(y, s["pocket_C_vs_vina_box_rmsd_A"], color=[point_colour(c, shortlist[c]) for c in s["compound"]],
             height=0.72, lw=0)
    axD.axvline(AGREE_A, color=fs.COL["alert"], lw=fs.LW["thr"], ls=fs.THR_DASH, zorder=3)
    axD.set_yticks(y); axD.set_yticklabels(s["compound"])
    axD.set_xlim(0, 16); axD.set_ylim(-0.7, len(s) - 0.3); axD.invert_yaxis()
    axD.set_xlabel("DiffDock-Pocket to Vina pose RMSD (\u00c5)")
    axD.annotate("2.0 \u00c5", xy=(AGREE_A, 1.0), xycoords=("data", "axes fraction"), xytext=(3, 1), textcoords="offset points",
                 color=fs.COL["alert"], fontsize=fs.FS["note"], ha="left", va="bottom")
    axD.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=DARK, label="shortlisted (8)"),
                        plt.Rectangle((0, 0), 1, 1, color=LIGHT, label="other survivors (7)")],
               frameon=False, loc="upper right", handlelength=1.0)
    fs.panel_label_fig(fig, "D", 0.03, top_h + 0.15 + 0.10, H)

    for c in surv:
        col = point_colour(c, shortlist[c])
        solid = shortlist[c] or c in fs.COL
        axE.plot(xb[c], yb[c], "o", ms=5.0, mfc=col if solid else "white", mec=col if solid else DARK, mew=0.9, zorder=3)
    for c in ("SNP00157", "SNP00163", "SNP00006"):
        axE.annotate(c[-3:], (xb[c], yb[c]), xytext=(4, 3), textcoords="offset points", fontsize=fs.FS["note"], color=point_colour(c, True))
    axE.plot([0, 16], [0, 16], color=fs.COL["grey"], lw=fs.LW["ref"], zorder=1)
    axE.axvline(SITE_A, color=fs.COL["grey"], lw=fs.LW["ref"], ls=fs.REF_DASH, zorder=1)
    axE.axhline(SITE_A, color=fs.COL["grey"], lw=fs.LW["ref"], ls=fs.REF_DASH, zorder=1)
    axE.set_xlim(0, 16); axE.set_ylim(0, 16); axE.set_aspect("equal", adjustable="box")
    axE.set_xlabel("Vina in the 22 \u00c5 site box:\ndistance from the site (\u00c5)")
    axE.set_ylabel("Vina over the whole receptor (\u00c5)")
    axE.legend(handles=[plt.Line2D([], [], marker="o", ls="", mfc=DARK, mec=DARK, ms=5.0, label="shortlisted"),
                        plt.Line2D([], [], marker="o", ls="", mfc="white", mec=DARK, ms=5.0, label="other survivors")],
               frameon=False, loc="upper left", borderaxespad=0.1)
    fs.panel_label_fig(fig, "E", 0.545 * W, top_h + 0.15 + 0.10, H)
    stem = OUT / "Figure_6_pocket"
    fs.save_figure(fig, stem)
    prov = {"figure": "Figure 6 (pocket version, harmonised style)", "script": "figure6_pocket.py", "adr": "ADR-AIIM-366, ADR-AIIM-367", "checks": checks,
            "panel_png_sha256": hashes,
            "inputs_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in (
                RUNS / "docking_matched_filter_pocket_2026-10-09" / "filter_pocket_poses.csv",
                RUNS / "docking_matched_search_space_2026-10-09" / "matched_poses.csv")}}
    (OUT / "Figure_6_pocket_provenance.json").write_text(json.dumps(prov, indent=1), "utf-8")
    print("written", stem.with_suffix(".png"))


if __name__ == "__main__":
    main()

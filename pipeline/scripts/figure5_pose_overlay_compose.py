"""
figure5_pose_overlay_compose.py -- compose the two PyMOL renders into
Figure 5, with panel labels, a colour legend and the RMSD figures that
license the reading, matching this manuscript's other figures (points-based
layout, 7 pt type floor, source digests recorded).

Run render_pose_overlay_pymol.py first to produce the two source PNGs.

Usage:
    python figure5_pose_overlay_compose.py
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from PIL import Image


def _root() -> Path:
    marker = "00_project_control/08_DECISION_LOG.md"
    for cand in (Path(r"D:\Agentic VS Workflow"), Path("/mnt/d/Agentic VS Workflow")):
        if (cand / marker).is_file():
            return cand
    raise SystemExit("project root not found")


ROOT = _root()
FIGDIR = ROOT / "30_manuscript/figures"
PANEL_A = FIGDIR / "_panelA_control.png"
PANEL_B = FIGDIR / "_panelB_candidate.png"

REF_GREY = "#7a7a7a"
VINA_BLUE = "#1f5fa8"
DIFFDOCK_RED = "#d6403f"

MIN_PT = 7.0


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> int:
    for p in (PANEL_A, PANEL_B):
        if not p.is_file():
            raise SystemExit(f"missing source render: {p} -- run "
                             "render_pose_overlay_pymol.py first")

    imgA = Image.open(PANEL_A)
    imgB = Image.open(PANEL_B)

    fig = plt.figure(figsize=(7.48, 3.55))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1], wspace=0.03,
                          left=0.005, right=0.995, top=0.90, bottom=0.14)

    axA = fig.add_subplot(gs[0, 0])
    axA.imshow(imgA)
    axA.axis("off")
    axA.set_title("A  Positive control: ZM241385", loc="left",
                  fontsize=8.6, weight="bold", pad=5)
    axA.text(0.5, -0.045,
             "Vina 0.379 \u00c5  \u00b7  DiffDock 1.92 \u00c5  (both vs crystal)",
             transform=axA.transAxes, ha="center", va="top", fontsize=MIN_PT)

    axB = fig.add_subplot(gs[0, 1])
    axB.imshow(imgB)
    axB.axis("off")
    axB.set_title("B  Tier-A candidate: NPLIB01643", loc="left",
                  fontsize=8.6, weight="bold", pad=5)
    axB.text(0.5, -0.045, "Vina vs DiffDock: 8.80 \u00c5 apart, both on-pocket",
             transform=axB.transAxes, ha="center", va="top", fontsize=MIN_PT)

    legend_elems = [
        Line2D([0], [0], color=REF_GREY, lw=1.5, label="crystallographic pose (thin = reference)"),
        Line2D([0], [0], color=VINA_BLUE, lw=3, label="Vina rank-1 pose"),
        Line2D([0], [0], color=DIFFDOCK_RED, lw=3, label="DiffDock rank-1 pose"),
    ]
    fig.legend(handles=legend_elems, loc="lower center", ncol=3, frameon=False,
              fontsize=MIN_PT, bbox_to_anchor=(0.5, 0.0),
              handlelength=1.8, columnspacing=1.4, handletextpad=0.6)

    out_stem = FIGDIR / "figure5_pose_overlay_2026-08-24"
    fig.savefig(out_stem.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(out_stem.with_suffix(".svg"), bbox_inches="tight")

    prov = {
        "figure": "Figure 5 -- pose overlay: agreement (control) vs disagreement (candidate)",
        "generated": "2026-08-24",
        "renderer": "PyMOL 3.1.0 (open-source, conda-forge), "
                    "10_new_system/scripts/render_pose_overlay_pymol.py",
        "rmsd_values_A_vs_crystal": {"vina_redock": 0.379, "diffdock": 1.917,
                                     "verified_against": "manuscript SS3.5/2.5, "
                                     "recomputed RDKit rdMolAlign.CalcRMS, in-place, "
                                     "symmetry-corrected"},
        "rmsd_value_B_vina_vs_diffdock": 8.800,
        "cartoon_scope_note": ("Cartoon restricted to residues within 13 A of the "
                              "displayed ligand poses. The full receptor carries "
                              "three real chain discontinuities (two unresolved "
                              "loops, one BRIL-excision boundary; residue "
                              "numbering gaps 208-219, 307-1001, 1042-1058), "
                              "none near the orthosteric site; excluded by the "
                              "pocket-local view rather than hidden by any "
                              "rendering override."),
        "inputs": {
            "receptor": "20_new_runs/snp_screen_2026-08-04/G2_preparation/"
                       "4EIY_protein_plus_waters.pdb",
            "zma_crystal": "20_new_runs/snp_screen_2026-08-04/G2_preparation/"
                          "ZMA_reference_4EIY.sdf",
            "zma_vina_redock": "20_new_runs/snp_screen_2026-08-04/G3_redocking/"
                              "ZMA_redocked_4EIY.sdf",
            "zma_diffdock": "20_new_runs/diffdock_orthogonal_2026-08-23/results/"
                           "ZM241385_control/rank1.sdf",
            "n1643_vina": "20_new_runs/np_library_funnel_2026-08-15/"
                         "full_poses_chunk1/NPLIB01643.sdf",
            "n1643_diffdock": "20_new_runs/diffdock_orthogonal_2026-08-23/results/"
                             "NPLIB01643/rank1.sdf",
        },
        "input_sha256": {
            name: sha256(ROOT / relpath) for name, relpath in {
                "receptor": "20_new_runs/snp_screen_2026-08-04/G2_preparation/4EIY_protein_plus_waters.pdb",
                "zma_crystal": "20_new_runs/snp_screen_2026-08-04/G2_preparation/ZMA_reference_4EIY.sdf",
                "zma_vina_redock": "20_new_runs/snp_screen_2026-08-04/G3_redocking/ZMA_redocked_4EIY.sdf",
                "zma_diffdock": "20_new_runs/diffdock_orthogonal_2026-08-23/results/ZM241385_control/rank1.sdf",
                "n1643_vina": "20_new_runs/np_library_funnel_2026-08-15/full_poses_chunk1/NPLIB01643.sdf",
                "n1643_diffdock": "20_new_runs/diffdock_orthogonal_2026-08-23/results/NPLIB01643/rank1.sdf",
            }.items()
        },
        "outputs": {p.name: sha256(p) for p in
                    (out_stem.with_suffix(".png"), out_stem.with_suffix(".svg"))},
    }
    (FIGDIR / "figure5_pose_overlay_provenance.json").write_text(
        json.dumps(prov, indent=2), encoding="utf-8")

    print(f"written -> {out_stem}.png / .svg")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

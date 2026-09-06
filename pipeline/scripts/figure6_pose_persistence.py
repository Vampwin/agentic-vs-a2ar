"""figure6_pose_persistence.py -- do the predicted poses stay where they were put?

The pose-persistence result (ADR-AIIM-149) is the one follow-up in this study
that argues against a candidate the study itself selected, and it reached the
manuscript as prose only. The numbers separate cleanly enough that a reader
should be able to see the separation rather than reconstruct it from three
means: one compound leaves the pocket in all three replicates, the other two
stay and drift to roughly twice the control's deviation.

Two panels, because the two facts are on incompatible scales:

  A  Every replicate against time, log y. NPLIB01185 runs to sixty ångström and
     the other two sit near five; on a linear axis the second group is a flat
     line at the bottom and the comparison is unreadable. The control is drawn
     as a band, not a line, because it is one trajectory and drawing it like the
     nine would imply it carries the same weight.

  B  The threshold ladder, per replicate. The 2.0 Å redocking criterion is far
     too strict for a 310 K trajectory -- the control itself keeps only 17 % of
     frames below it -- so the analysis reports a ladder rather than a verdict,
     and the figure shows the ladder rather than choosing a rung.

**No threshold line is drawn as a pass mark.** This project pre-registered no
pose-persistence criterion, and a horizontal line at 2 or 4 Å would read as one.
The control band is a reference, and it is labelled as a single trajectory.

Usage:
    python figure6_pose_persistence.py [--out DIR]
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402


def _root() -> Path:
    marker = "00_project_control/08_DECISION_LOG.md"
    for cand in (Path(r"D:\Agentic VS Workflow"), Path("/mnt/d/Agentic VS Workflow")):
        if (cand / marker).is_file():
            return cand
    raise SystemExit("project root not found")


ROOT = _root()
RUNS = ROOT / "20_new_runs"
CAND = RUNS / "md_candidates_2026-08-22"
SUMMARY = RUNS / "md_pose_persistence_2026-08-31/md_pose_persistence.json"
OCCUPANCY = RUNS / "md_pose_persistence_2026-08-31/pocket_occupancy.json"
MOTIF = RUNS / "md_recognition_motif_2026-09-02"
# The hydrogen-bond cutoff is the one fixed in ADR-AIIM-161 before the
# measurement ran, and matches the polar-contact definition the §2.5
# interaction fingerprint already used. It is drawn in panel C because it
# is a pre-registered criterion, unlike the RMSD thresholds in panel B,
# where no line is drawn for exactly that reason.
HBOND_CUTOFF_A = 3.5
# NPLIB01185's RMSD is corrupted by a periodic image (ADR-AIIM-156) and is not
# drawn. Plotting it beside two trustworthy traces would put an artefact in the
# figure with nothing marking it as one.
RMSD_TRUSTWORTHY = ["NPLIB01643", "NPLIB03056"]

COMPOUNDS = ["NPLIB01185", "NPLIB01643", "NPLIB03056"]
REPS = [1, 2, 3]
LADDER = ["2.0", "3.0", "4.0", "5.0"]
# One colour per compound, shared across panels; replicates differ by line style
# so a reader tracks the compound first and the replicate second.
CCOL = {"NPLIB01185": "#b5651d", "NPLIB01643": "#2c5f8a", "NPLIB03056": "#7a9a3c"}
DASH = {1: "solid", 2: (0, (3.2, 1.6)), 3: (0, (1.1, 1.2))}
CONTROL_MEAN, CONTROL_SD = 2.59, 0.51

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.2,
    "axes.linewidth": 0.7, "axes.edgecolor": "#333333",
    "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "savefig.dpi": 600, "figure.dpi": 150,
})


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read_xvg(cid: str, rep: int, name: str) -> tuple[np.ndarray, np.ndarray]:
    xvg = CAND / cid / "04_run" / f"analysis_prod_r{rep}" / name
    t, v = [], []
    for line in xvg.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = line.strip()
        if not line or line[0] in "#@":
            continue
        p = line.split()
        if len(p) >= 2:
            t.append(float(p[0]))               # already ns; the xvg says so
            v.append(float(p[1]) * 10.0)        # nm -> angstrom
    return np.asarray(t), np.asarray(v)


def read_occ(cid: str, rep: int):
    return read_xvg(cid, rep, "pocket_start_dist.xvg")


def read_rmsd(cid: str, rep: int):
    return read_xvg(cid, rep, "ligand_rmsd.xvg")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "30_manuscript/figures"))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    summ = json.loads(SUMMARY.read_text(encoding="utf-8"))["candidates"]

    fig = plt.figure(figsize=(7.48, 5.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.75, 1.0],
                          height_ratios=[1.0, 0.92], wspace=0.24, hspace=0.28,
                          left=0.075, right=0.985, top=0.925, bottom=0.085)
    axA, axB = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, :])

    # ---- A: pocket occupancy, the measurement PBC cannot corrupt -----------
    for cid in COMPOUNDS:
        for r in REPS:
            t, v = read_occ(cid, r)
            axA.plot(t, v, color=CCOL[cid], lw=0.55, alpha=0.85, zorder=3,
                     linestyle=DASH[r])
    axA.axhline(10.0, color="#666666", lw=0.8, dashes=(4, 2), zorder=1)
    axA.set_ylim(0, 12)
    axA.set_xlim(0, 100)
    axA.set_yticks([0, 2, 4, 6, 8, 10, 12])
    axA.set_yticklabels(["0", "2", "4", "6", "8", "10", "12"], fontsize=7)
    axA.set_xticks([0, 20, 40, 60, 80, 100])
    axA.set_xticklabels(["0", "20", "40", "60", "80", "100"], fontsize=7)
    axA.set_xlabel("time (ns)", fontsize=7.6)
    axA.set_ylabel("distance to the docked binding site (Å)", fontsize=7.6)
    for sp in ("top", "right"):
        axA.spines[sp].set_visible(False)
    axA.set_title("A", loc="left", fontsize=8.4, weight="bold", pad=6)
    handles = [plt.Line2D([], [], color=CCOL[c], lw=1.4, label=c)
               for c in COMPOUNDS]
    axA.legend(handles=handles, fontsize=6.6, frameon=False, loc="lower right",
               ncol=3, handlelength=1.8, borderaxespad=0.3)

    # ---- B: the threshold ladder -------------------------------------------
    xs = np.arange(len(LADDER))
    width = 0.26
    for i, cid in enumerate(RMSD_TRUSTWORTHY):
        fr = [[summ[cid]["replicates"][f"r{r}"]["fraction_below"][L] * 100
               for r in REPS] for L in LADDER]
        means = [np.mean(f) for f in fr]
        axB.bar(xs + (i - 0.5) * width, means, width * 0.92,
                color=CCOL[cid], alpha=0.85, lw=0, zorder=2)
        for j, f in enumerate(fr):
            axB.plot([xs[j] + (i - 0.5) * width] * 3, f, ".", ms=2.6,
                     color="#222222", zorder=3)
    axB.axhline(100, color="#666666", lw=0.8, zorder=1)
    axB.set_xticks(xs)
    axB.set_xticklabels([f"< {L} Å" for L in LADDER], fontsize=7)
    axB.set_ylim(0, 104)
    axB.set_yticks([0, 25, 50, 75, 100])
    axB.set_yticklabels(["0", "25", "50", "75", "100"], fontsize=7)
    axB.set_ylabel("frames below the threshold (%)", fontsize=7.6)
    for sp in ("top", "right"):
        axB.spines[sp].set_visible(False)
    axB.set_title("B", loc="left", fontsize=8.4, weight="bold", pad=6)

    # ---- C: the interaction that defines antagonism, over time -------------
    def read_motif(system):
        xvg = MOTIF / system / "N6_55.xvg"
        tt, vv = [], []
        for line in xvg.read_text(encoding="utf-8", errors="ignore").splitlines():
            line = line.strip()
            if not line or line[0] in "#@":
                continue
            parts = line.split()
            if len(parts) >= 2:
                tt.append(float(parts[0]))          # ns; -tu ns was passed
                vv.append(float(parts[1]) * 10.0)   # nm -> angstrom
        return np.asarray(tt), np.asarray(vv)

    for cid in COMPOUNDS:
        for r in REPS:
            t_, v_ = read_motif(f"{cid}_r{r}")
            axC.plot(t_, v_, color=CCOL[cid], lw=0.45, alpha=0.5, zorder=2,
                     linestyle=DASH[r])
    tc, vc = read_motif("control_ZM241385")
    axC.axhline(HBOND_CUTOFF_A, color="#c0392b", lw=1.0, dashes=(5, 2), zorder=5)
    axC.plot(tc, vc, color="#111111", lw=1.25, alpha=1.0, zorder=6)
    axC.set_xlim(0, 100)
    axC.set_ylim(1.8, 10.2)
    axC.set_xticks([0, 20, 40, 60, 80, 100])
    axC.set_xticklabels(["0", "20", "40", "60", "80", "100"], fontsize=7)
    axC.set_yticks([2, 4, 6, 8, 10])
    axC.set_yticklabels(["2", "4", "6", "8", "10"], fontsize=7)
    axC.set_xlabel("time (ns)", fontsize=7.6)
    axC.set_ylabel("ligand to Asn253 sidechain (\u00c5)", fontsize=7.6)
    for sp in ("top", "right"):
        axC.spines[sp].set_visible(False)
    axC.set_title("C", loc="left", fontsize=8.4, weight="bold", pad=6)
    hC = ([plt.Line2D([], [], color="#333333", lw=1.4,
                      label="ZM241385 control")]
          + [plt.Line2D([], [], color=CCOL[c], lw=1.4, label=c) for c in COMPOUNDS])
    axC.legend(handles=hC, fontsize=6.6, frameon=False, loc="upper left",
               ncol=4, handlelength=1.8, borderaxespad=0.3)

    stem = out / "figure6_pose_persistence_2026-08-31"
    fig.savefig(stem.with_suffix(".png"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")

    prov = {
        "figure": "Figure 6, pose persistence",
        "authority": "ADR-AIIM-149",
        "summary_artifact": {"file": str(SUMMARY.relative_to(ROOT).as_posix()),
                             "sha256": sha256(SUMMARY)},
        "trajectory_analyses": {
            f"{c}_r{r}": sha256(CAND / c / "04_run" /
                                f"analysis_prod_r{r}" / "ligand_rmsd.xvg")
            for c in COMPOUNDS for r in REPS},
        "control": {"mean_A": CONTROL_MEAN, "sd_A": CONTROL_SD,
                    "n_replicates": 1,
                    "note": "drawn as a band and labelled as one trajectory, "
                            "because it is not three"},
        "panel_C": {
            "source": str(MOTIF.relative_to(ROOT).as_posix()),
            "authority": "ADR-AIIM-161",
            "measured_by": "gmx mindist, minimum-image convention applied by the "
                           "tool itself",
            "cutoff_A": HBOND_CUTOFF_A,
            "why_a_line_is_drawn_here_and_not_in_B": "3.5 A was fixed before the "
                "measurement ran and matches the §2.5 polar-contact definition. "
                "No pose-persistence RMSD threshold was ever pre-registered, so "
                "panel B draws none.",
            "covers_all_three": "unlike panel B, which cannot show NPLIB01185",
        },
        "no_threshold_drawn": "This project pre-registered no pose-persistence "
                              "criterion. A horizontal line at 2 or 4 A would "
                              "read as a pass mark and none exists.",
    }
    (out / "figure6_pose_persistence_provenance.json").write_text(
        json.dumps(prov, indent=2), encoding="utf-8")

    for c in COMPOUNDS:
        s = summ[c]
        print(f"{c}: mean {s['mean_of_replicate_means_A']:.2f} A, "
              f"spread {s['spread_of_replicate_means_A']:.2f}, "
              f"finals {s['final_rmsd_range_A'][0]:.1f}-{s['final_rmsd_range_A'][1]:.1f}")
    print(f"written -> {stem}.png / .svg")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

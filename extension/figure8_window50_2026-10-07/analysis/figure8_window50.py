"""Figure 8 over a chosen window (default the first 50 ns), from the time series.

The plotting code is figure6_pose_persistence.py with ONE change: every series is
cut at the window. Panel layout, colours, line styles, axis limits and fonts are
unchanged, so the 100 ns figure and the 50 ns figure can be compared like for like.

Two anchors, both run before the 50 ns figure is trusted:

  1. The threshold ladder of panel B is computed here from ligand_rmsd.xvg and
     must equal the values ADR-AIIM-321 wrote to md_window50.json, to 4 places.
  2. Run with --window 100 the same code must reproduce the figure published in
     v8 (Figure_8.png) pixel for pixel. If it does not, the 50 ns figure differs
     from the published one in more than the window.

Legend, written before the artwork (manuscript-figures-tables skill, phase 2):

  Figure 8. Over the first 50 ns the three originally prioritised compounds stayed
  in the pocket but did not hold the control's Asn253 contact. Membrane simulations
  of NPLIB01185, NPLIB01643 and NPLIB03056 (three replicates each, differing only in
  velocity seed; line style identifies the replicate) and of the ZM241385 control
  (one trajectory), over the first 50 ns, the window shared with the control
  replicates; the 100 ns summaries are in Supplementary Table S9. (A) Minimum
  distance between the ligand and the binding site, defined as the receptor C-alpha
  atoms within 10 A of the ligand in the energy-minimised starting structure; the
  dashed line marks that 10 A radius. (B) Fraction of frames with ligand heavy-atom
  RMSD, after receptor C-alpha fitting and against the pre-equilibration structure,
  below each threshold, for NPLIB01643 and NPLIB03056: bars are the mean of three
  replicates, dots the individual replicates; no threshold is a pass mark.
  NPLIB01185 RMSD is omitted because of an unresolved periodic-image artifact.
  (C) Minimum distance from the ligand's polar atoms to the Asn253 side chain; the
  red dashed line marks the 3.5 A hydrogen-bond cutoff fixed before the measurement.
  The black trace is the control, a single trajectory.

`10_new_system` is not modified; no published artifact is overwritten.
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

ROOT = Path(r"D:\Agentic VS Workflow")
RUNS = ROOT / "20_new_runs"
CAND = RUNS / "md_candidates_2026-08-22"
MOTIF = RUNS / "md_recognition_motif_2026-09-02"
WIN50 = RUNS / "md_window50_2026-10-07" / "md_window50.json"
OUTDIR = RUNS / "figure8_window50_2026-10-07"

HBOND_CUTOFF_A = 3.5
RMSD_TRUSTWORTHY = ["NPLIB01643", "NPLIB03056"]       # 01185: periodic-image artefact
COMPOUNDS = ["NPLIB01185", "NPLIB01643", "NPLIB03056"]
REPS = [1, 2, 3]
LADDER = ["2.0", "3.0", "4.0", "5.0"]
CCOL = {"NPLIB01185": "#b5651d", "NPLIB01643": "#2c5f8a", "NPLIB03056": "#7a9a3c"}
DASH = {1: "solid", 2: (0, (3.2, 1.6)), 3: (0, (1.1, 1.2))}
EPS = 1e-3

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.2,
    "axes.linewidth": 0.7, "axes.edgecolor": "#333333",
    "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "savefig.dpi": 600, "figure.dpi": 150,
})


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read_series(p: Path, full_ns: float = 100.0):
    """Time (ns) and value (nm -> angstrom). Asserts the time axis really is ns:
    the first execution of the 50 ns analysis cut a picoseconds file as ns."""
    t, v = [], []
    for line in p.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = line.strip()
        if not line or line[0] in "#@":
            continue
        parts = line.split()
        if len(parts) >= 2:
            t.append(float(parts[0]))
            v.append(float(parts[1]) * 10.0)
    t, v = np.asarray(t), np.asarray(v)
    assert abs(t[-1] - full_ns) < 0.5, f"{p.name}: last time {t[-1]}, not ~{full_ns} ns"
    return t, v


def cut(t, v, tmax):
    m = t <= tmax + EPS
    return t[m], v[m]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--window", type=float, default=50.0)
    ap.add_argument("--out", default=str(OUTDIR))
    ap.add_argument("--stem", default="")
    ap.add_argument("--control-runs", type=int, default=1, choices=(1, 3),
                    help="1 = the original control only (reproduces the published figure); "
                         "3 = also control replicates 2 and 3 (ADR-AIIM-324, 325)")
    args = ap.parse_args()
    tmax = args.window
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    stem = out / (args.stem or f"Figure_8_{int(tmax)}ns")

    # ---- data, read once -------------------------------------------------
    inputs = {}

    def need(p):
        inputs[str(p.relative_to(ROOT)).replace("\\", "/")] = sha(p)
        return p

    occ, rmsd, n655 = {}, {}, {}
    for cid in COMPOUNDS:
        for r in REPS:
            d = CAND / cid / "04_run" / f"analysis_prod_r{r}"
            occ[(cid, r)] = cut(*read_series(need(d / "pocket_start_dist.xvg")), tmax)
            rmsd[(cid, r)] = cut(*read_series(need(d / "ligand_rmsd.xvg")), tmax)
            n655[(cid, r)] = cut(*read_series(need(MOTIF / f"{cid}_r{r}" / "N6_55.xvg")), tmax)
    ctrl = cut(*read_series(need(MOTIF / "control_ZM241385" / "N6_55.xvg")), tmax)
    ctrl_runs = [("original", ctrl)]
    if args.control_runs == 3:
        REPD = RUNS / "md_control_replicates_2026-10-06"
        for nm in ("r2", "r3"):
            ctrl_runs.append((nm, cut(*read_series(need(REPD / f"analysis_prod_{nm}" / "N6_55.xvg"),
                                                    full_ns=53.2 if nm == "r2" else 50.0), tmax)))   # r2 ran to 53.2 ns, cut here

    # ---- anchor 1: the ladder equals the stored 50 ns analysis ------------
    ladder = {cid: {f"r{r}": {L: float(np.mean(rmsd[(cid, r)][1] < float(L)))
                               for L in LADDER} for r in REPS}
              for cid in RMSD_TRUSTWORTHY}
    stored = json.loads(need(WIN50).read_text("utf-8"))["windows"][str(int(tmax))]["candidates"] \
        if str(int(tmax)) in ("50", "100") else None
    ladder_ok = None
    if stored is not None:
        ladder_ok = all(
            abs(ladder[cid][f"r{r}"][L]
                - stored[cid]["replicates"][f"r{r}"]["rmsd"]["fraction_below"][L]) < 5e-5
            for cid in RMSD_TRUSTWORTHY for r in REPS for L in LADDER)
        print(f"anchor 1, ladder equals md_window50.json at {tmax:g} ns: {ladder_ok}")
        assert ladder_ok, "panel B disagrees with the stored analysis"

    # ---- figure ----------------------------------------------------------
    fig = plt.figure(figsize=(7.48, 5.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.75, 1.0], height_ratios=[1.0, 0.92],
                          wspace=0.24, hspace=0.28,
                          left=0.075, right=0.985, top=0.925, bottom=0.085)
    axA, axB = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, :])
    xt = [0, 20, 40, 60, 80, 100] if tmax >= 100 else [0, 10, 20, 30, 40, 50]
    xl = [str(x) for x in xt]

    for cid in COMPOUNDS:
        for r in REPS:
            t, v = occ[(cid, r)]
            axA.plot(t, v, color=CCOL[cid], lw=0.55, alpha=0.85, zorder=3, linestyle=DASH[r])
    axA.axhline(10.0, color="#666666", lw=0.8, dashes=(4, 2), zorder=1)
    axA.set_ylim(0, 12); axA.set_xlim(0, tmax)
    axA.set_yticks([0, 2, 4, 6, 8, 10, 12])
    axA.set_yticklabels(["0", "2", "4", "6", "8", "10", "12"], fontsize=7)
    axA.set_xticks(xt); axA.set_xticklabels(xl, fontsize=7)
    axA.set_xlabel("time (ns)", fontsize=7.6)
    axA.set_ylabel("distance to the docked binding site (\u00c5)", fontsize=7.6)
    for sp in ("top", "right"):
        axA.spines[sp].set_visible(False)
    axA.set_title("A", loc="left", fontsize=8.4, weight="bold", pad=6)
    axA.legend(handles=[plt.Line2D([], [], color=CCOL[c], lw=1.4, label=c) for c in COMPOUNDS],
               fontsize=6.6, frameon=False, loc="lower right", ncol=3,
               handlelength=1.8, borderaxespad=0.3)

    xs = np.arange(len(LADDER)); width = 0.26
    for i, cid in enumerate(RMSD_TRUSTWORTHY):
        fr = [[ladder[cid][f"r{r}"][L] * 100 for r in REPS] for L in LADDER]
        axB.bar(xs + (i - 0.5) * width, [np.mean(f) for f in fr], width * 0.92,
                color=CCOL[cid], alpha=0.85, lw=0, zorder=2)
        for j, f in enumerate(fr):
            axB.plot([xs[j] + (i - 0.5) * width] * 3, f, ".", ms=2.6, color="#222222", zorder=3)
    axB.axhline(100, color="#666666", lw=0.8, zorder=1)
    axB.set_xticks(xs)
    axB.set_xticklabels([f"< {L} \u00c5" for L in LADDER], fontsize=7)
    axB.set_ylim(0, 104); axB.set_yticks([0, 25, 50, 75, 100])
    axB.set_yticklabels(["0", "25", "50", "75", "100"], fontsize=7)
    axB.set_ylabel("frames below the threshold (%)", fontsize=7.6)
    for sp in ("top", "right"):
        axB.spines[sp].set_visible(False)
    axB.set_title("B", loc="left", fontsize=8.4, weight="bold", pad=6)

    for cid in COMPOUNDS:
        for r in REPS:
            t_, v_ = n655[(cid, r)]
            axC.plot(t_, v_, color=CCOL[cid], lw=0.45, alpha=0.5, zorder=2, linestyle=DASH[r])
    axC.axhline(HBOND_CUTOFF_A, color="#c0392b", lw=1.0, dashes=(5, 2), zorder=5)
    CSTY = {"original": "solid", "r2": (0, (3.2, 1.6)), "r3": (0, (1.1, 1.2))}
    for nm, (tc_, vc_) in ctrl_runs:
        axC.plot(tc_, vc_, color="#111111", lw=1.25, alpha=1.0, zorder=6,
                 linestyle=CSTY[nm] if args.control_runs == 3 else "solid")
    axC.set_xlim(0, tmax); axC.set_ylim(1.8, 10.2)
    axC.set_xticks(xt); axC.set_xticklabels(xl, fontsize=7)
    axC.set_yticks([2, 4, 6, 8, 10]); axC.set_yticklabels(["2", "4", "6", "8", "10"], fontsize=7)
    axC.set_xlabel("time (ns)", fontsize=7.6)
    axC.set_ylabel("ligand to Asn253 sidechain (\u00c5)", fontsize=7.6)
    for sp in ("top", "right"):
        axC.spines[sp].set_visible(False)
    axC.set_title("C", loc="left", fontsize=8.4, weight="bold", pad=6)
    axC.legend(handles=[plt.Line2D([], [], color="#333333", lw=1.4, label="ZM241385 control" + (" (3 runs)" if args.control_runs == 3 else ""))]
               + [plt.Line2D([], [], color=CCOL[c], lw=1.4, label=c) for c in COMPOUNDS],
               fontsize=6.6, frameon=False, loc="upper left", ncol=4,
               handlelength=1.8, borderaxespad=0.3)

    fig.savefig(stem.with_suffix(".png"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")

    # ---- facts the legend asserts, computed from the plotted data ---------
    pocket_max = {f"{c}_r{r}": round(float(occ[(c, r)][1].max()), 2) for c in COMPOUNDS for r in REPS}
    claims = {
        "window_ns": tmax,
        "all_replicates_stay_within_10A_of_site": bool(all(v < 10.0 for v in pocket_max.values())),
        "max_pocket_distance_A": pocket_max,
        "control_fraction_n655_below_3.5A": round(float(np.mean(ctrl[1] <= HBOND_CUTOFF_A)), 4),
        "control_runs_n655_below_3.5A": {nm: round(float(np.mean(v[1] <= HBOND_CUTOFF_A)), 4)
                                        for nm, v in ctrl_runs},
        "candidate_fraction_n655_below_3.5A": {
            c: [round(float(np.mean(n655[(c, r)][1] <= HBOND_CUTOFF_A)), 4) for r in REPS]
            for c in COMPOUNDS},
        "frames": {"control_n655": int(len(ctrl[0])),
                   "candidate_rmsd": int(len(rmsd[(COMPOUNDS[1], 1)][0]))},
        "ladder_percent_by_replicate": {c: {r: {L: round(v * 100, 2) for L, v in d.items()}
                                            for r, d in ladder[c].items()}
                                        for c in RMSD_TRUSTWORTHY},
        "anchor_ladder_equals_md_window50_json": ladder_ok}
    prov = {"figure": "Figure 8", "adr": "ADR-AIIM-320, ADR-AIIM-321",
            "script": "figure8_window50.py", "window_ns": tmax,
            "changed_from_published": "the window only; layout, colours, line styles, "
                                      "axis limits and fonts are those of figure6_pose_persistence.py",
            "claims_checked_against_plotted_data": claims, "inputs_sha256": inputs}
    (out / f"{stem.name}_provenance.json").write_text(json.dumps(prov, indent=2), encoding="utf-8")
    print(json.dumps({k: claims[k] for k in ("window_ns", "all_replicates_stay_within_10A_of_site",
                                             "control_fraction_n655_below_3.5A",
                                             "candidate_fraction_n655_below_3.5A")}, indent=1))
    print("written ->", stem.with_suffix(".png"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

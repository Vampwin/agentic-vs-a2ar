"""Figure 7 (membrane simulations) with SNP00157 and SNP00163 added, over the first 50 ns.

Replaces Figure 7 of the v8 working article (formerly numbered 8; script figure8_window50.py), which showed
the three originally prioritised compounds and the control only. Layout, fonts and the colours of the
earlier compounds are those of that script; what is new is the two SNP-Lib compounds in every panel, the
control as a reference bar in panel B, and a panel D that shows the quantity the text reports for every run.

Legend (manuscript-figures-tables skill, phase 2):

  The figure legend was written before the artwork. It is the legend of the manuscript's figure and is not
  reproduced here, because the manuscript's text is not deposited.

Anchors, all required before the figure is trusted:
  1. panel B's ladder equals the stored values (md_window50.json for the NP compounds,
     md_control_replicates_window50.json for the control) to 4 places, and for SNP00157 and SNP00163 the
     fraction below 4 A equals adr331_results.json
  2. panel D's persistence equals the stored values for every run: control (ADR-AIIM-325), NP compounds
     (md_window50.json), SNP compounds (adr331_results.json), to 4 places
  3. every pocket distance is below 10 A (the claim in the legend), computed from the plotted series

Every number in the figure comes from a time series read here; nothing is typed. `10_new_system` is not modified.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import sys                        # noqa: E402
sys.path.insert(0, "D:/Agentic VS Workflow/20_new_runs/figures_style_harmonised_2026-10-10")
import figstyle as fs             # noqa: E402
import numpy as np                # noqa: E402

ROOT = Path(r"D:\Agentic VS Workflow")
RUNS = ROOT / "20_new_runs"
CAND = RUNS / "md_candidates_2026-08-22"
MOTIF = RUNS / "md_recognition_motif_2026-09-02"
CTRL_ANCHOR = RUNS / "md_control_anchor_2026-10-07" / "analysis_prod"
CTRL_REP = RUNS / "md_control_replicates_2026-10-06"
WIN50 = RUNS / "md_window50_2026-10-07" / "md_window50.json"
CTRL50 = RUNS / "md_control_replicates_window50_2026-10-07" / "md_control_replicates_window50.json"
ADR331 = RUNS / "md_snplib_duplicate_2026-10-08" / "adr331_results.json"
OUT = RUNS / "figures_docking_md_2026-10-09"

TMAX, EPS = 50.0, 1e-3
HB = 3.5
NP = ["NPLIB01185", "NPLIB01643", "NPLIB03056"]
SNP = ["SNP00157", "SNP00163"]
ALL = "--all" in sys.argv            # --all: also draw the three originally prioritised compounds (Supplementary Figure S2)
NPD = NP if ALL else []              # compounds of the earlier set that are drawn
RMSD_NP_OK = ["NPLIB01643", "NPLIB03056"]
LADDER = ["2.0", "3.0", "4.0", "5.0"]
CCOL = {"NPLIB01185": "#b5651d", "NPLIB01643": "#2c5f8a", "NPLIB03056": "#7a9a3c",
        "SNP00157": "#8c4f9f", "SNP00163": "#c2487a"}
CTRL_COL = "#111111"
DASH = {1: "solid", 2: (0, (3.2, 1.6)), 3: (0, (1.1, 1.2))}
MARK = {1: "o", 2: "s", 3: "^"}

fs.apply()
inputs = {}


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read(p: Path):
    inputs[str(p.relative_to(ROOT)).replace("\\", "/")] = sha(p)
    t, v = [], []
    for ln in p.read_text(errors="ignore").splitlines():
        ln = ln.strip()
        if ln and ln[0] not in "#@":
            a = ln.split(); t.append(float(a[0])); v.append(float(a[1]) * 10.0)
    t, v = np.array(t), np.array(v)
    assert t[-1] >= TMAX - 0.5, f"{p}: ends at {t[-1]} ns"
    m = t <= TMAX + EPS
    return t[m], v[m]


def main():
    A = {}      # pocket-from-start
    R = {}      # ligand rmsd
    N = {}      # N6.55 distance
    for c in NP:
        for r in (1, 2, 3):
            d = CAND / c / "04_run" / f"analysis_prod_r{r}"
            A[(c, r)] = read(d / "pocket_start_dist.xvg"); R[(c, r)] = read(d / "ligand_rmsd.xvg")
            N[(c, r)] = read(MOTIF / f"{c}_r{r}" / "N6_55.xvg")
    for c in SNP:
        for r in (1, 2):
            d = CAND / c / "04_run" / f"analysis_prod_r{r}"
            A[(c, r)] = read(d / "pocket_start_dist.xvg"); R[(c, r)] = read(d / "ligand_rmsd.xvg"); N[(c, r)] = read(d / "N6_55.xvg")
    CT = {"r1": CTRL_ANCHOR, "r2": CTRL_REP / "analysis_prod_r2", "r3": CTRL_REP / "analysis_prod_r3"}
    CN = {"r1": read(MOTIF / "control_ZM241385" / "N6_55.xvg"),
          "r2": read(CT["r2"] / "N6_55.xvg"), "r3": read(CT["r3"] / "N6_55.xvg")}
    CR = {k: read(v / "ligand_rmsd.xvg") for k, v in CT.items()}

    # ---- stored values for the anchors
    w50 = json.loads(WIN50.read_text("utf-8"))["windows"]["50"]["candidates"]
    c50 = json.loads(CTRL50.read_text("utf-8"))
    a331 = json.loads(ADR331.read_text("utf-8"))
    ctrl_stored = {"r1": c50["existing_control"]["50"], "r2": c50["replicates"]["r2"], "r3": c50["replicates"]["r3"]}
    persist = lambda v: float(np.mean(v <= HB + 1e-9))
    P = {("control", k): persist(CN[k][1]) for k in CN}
    P.update({(c, r): persist(N[(c, r)][1]) for (c, r) in N})
    ladder = lambda v: {L: float(np.mean(v < float(L))) for L in LADDER}

    anchors = {}
    ok = True
    for k, s in ctrl_stored.items():
        okp = abs(P[("control", k)] - s["n655"]["persistence"]) < 5e-5
        lad = ladder(CR[k][1]); okl = all(abs(lad[L] - s["rmsd"]["fraction_below"][L]) < 5e-5 for L in LADDER)
        anchors[f"control_{k}"] = {"persistence_equals_stored": okp, "ladder_equals_stored": okl}; ok &= okp and okl
    for c in NP:
        for r in (1, 2, 3):
            s = w50[c]["replicates"][f"r{r}"]
            okp = abs(P[(c, r)] - s["n655"]["persistence"]) < 5e-5
            okl = True
            if c in RMSD_NP_OK:
                lad = ladder(R[(c, r)][1]); okl = all(abs(lad[L] - s["rmsd"]["fraction_below"][L]) < 5e-5 for L in LADDER)
            anchors[f"{c}_r{r}"] = {"persistence_equals_stored": okp, "ladder_equals_stored": okl}; ok &= okp and okl
    for c in SNP:
        for r in (1, 2):
            s = a331["candidate_runs"][c][f"r{r}"]
            okp = abs(P[(c, r)] - s["n655_persistence"]) < 5e-5
            okl = abs(ladder(R[(c, r)][1])["4.0"] - s["rmsd_frac_lt_4A"]) < 5e-5
            anchors[f"{c}_r{r}"] = {"persistence_equals_stored": okp, "ladder_4A_equals_stored": okl}; ok &= okp and okl
    pocket_max = {f"{c}_r{r}": round(float(A[(c, r)][1].max()), 2) for (c, r) in A}
    all10 = all(v < 10.0 for v in pocket_max.values()); ok &= all10
    anchors["all_pocket_distances_below_10A"] = all10
    print("anchors:", "ALL MET" if ok else "FAILED"); assert ok, anchors

    # ---- figure
    fig = plt.figure(figsize=(7.48, 7.9 if ALL else 6.5))
    gs = fig.add_gridspec(3, 2, width_ratios=[1.75, 1.0], height_ratios=[1.0, 0.92, 0.62 if ALL else 0.36],
                          wspace=0.24, hspace=0.36, left=0.075, right=0.985, top=0.955, bottom=0.07)
    axA, axB = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    axC, axD = fig.add_subplot(gs[1, :]), fig.add_subplot(gs[2, :])
    xt = [0, 10, 20, 30, 40, 50]
    for ax in (axA, axB, axC, axD):
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)

    for (c, r), (t, v) in A.items():
        if c in NP and not ALL:
            continue
        axA.plot(t, v, color=CCOL[c], lw=fs.LW["series"], alpha=0.85, zorder=3, linestyle=DASH[r])
    axA.axhline(10.0, color=fs.COL["grey"], lw=fs.LW["ref"], ls=fs.REF_DASH, zorder=1)
    axA.set_ylim(0, 12); axA.set_xlim(0, TMAX); axA.set_yticks([0, 2, 4, 6, 8, 10, 12]); axA.set_xticks(xt)
    axA.tick_params(labelsize=fs.FS["tick"])
    axA.set_xlabel("time (ns)", fontsize=fs.FS["label"]); axA.set_ylabel("distance to the docked binding site (\u00c5)", fontsize=fs.FS["label"])
    axA.legend(handles=[plt.Line2D([], [], color=CCOL[c], lw=fs.LW["main"], label=c) for c in NPD + SNP], fontsize=fs.FS["tick"],
               frameon=False, loc="lower right", ncol=3 if ALL else 2, handlelength=1.6, borderaxespad=0.3, columnspacing=1.0)

    groups = [("ZM241385 control", "#555555", [ladder(CR[k][1]) for k in ("r1", "r2", "r3")])]
    groups += [(c, CCOL[c], [ladder(R[(c, r)][1]) for r in (1, 2, 3)]) for c in RMSD_NP_OK if ALL]
    groups += [(c, CCOL[c], [ladder(R[(c, r)][1]) for r in (1, 2)]) for c in SNP]
    xs = np.arange(len(LADDER)); width = 0.15
    for i, (name, col, lads) in enumerate(groups):
        off = (i - (len(groups) - 1) / 2) * width
        fr = [[l[L] * 100 for l in lads] for L in LADDER]
        axB.bar(xs + off, [np.mean(f) for f in fr], width * 0.92, color=col, alpha=0.85, lw=0, zorder=2)
        for j, f in enumerate(fr):
            axB.plot([xs[j] + off] * len(f), f, ".", ms=2.4, color="#222222", zorder=3)
    axB.axhline(100, color=fs.COL["grey"], lw=fs.LW["ref"], zorder=1)
    axB.set_xticks(xs); axB.set_xticklabels([f"< {L} \u00c5" for L in LADDER], fontsize=fs.FS["tick"])
    axB.set_ylim(0, 104); axB.set_yticks([0, 25, 50, 75, 100]); axB.tick_params(axis="y", labelsize=fs.FS["tick"])
    axB.set_ylabel("frames below the threshold (%)", fontsize=fs.FS["label"])
    axB.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=col, alpha=0.85, label=n) for n, col, _ in groups[:1]],   # the control only; the compounds carry the colours of A and C
               fontsize=fs.FS["tick"], frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.93), handlelength=1.0, borderaxespad=0.1)

    for (c, r), (t, v) in N.items():
        if c in NP and not ALL:
            continue
        axC.plot(t, v, color=CCOL[c], lw=fs.LW["series"], alpha=0.5, zorder=2, linestyle=DASH[r])
    axC.axhline(HB, color=fs.COL["alert"], lw=fs.LW["thr"], ls=fs.THR_DASH, zorder=5)
    for k in ("r1", "r2", "r3"):
        axC.plot(*CN[k], color=CTRL_COL, lw=fs.LW["main"], zorder=6, linestyle={"r1": "solid", "r2": DASH[2], "r3": DASH[3]}[k])
    axC.set_xlim(0, TMAX); axC.set_ylim(1.8, 10.2); axC.set_xticks(xt); axC.set_yticks([2, 4, 6, 8, 10]); axC.tick_params(labelsize=fs.FS["tick"])
    axC.set_xlabel("time (ns)", fontsize=fs.FS["label"]); axC.set_ylabel("ligand to Asn253 sidechain (\u00c5)", fontsize=fs.FS["label"])
    axC.legend(handles=[plt.Line2D([], [], color="#333333", lw=fs.LW["main"], label="ZM241385 control (3 runs)")]
               + [plt.Line2D([], [], color=CCOL[c], lw=fs.LW["main"], label=c) for c in NPD + SNP], fontsize=fs.FS["tick"], frameon=False,
               loc="upper left", ncol=6 if ALL else 3, handlelength=1.6, borderaxespad=0.3, columnspacing=1.0)

    rows = [("ZM241385 control", CTRL_COL, [("control", k) for k in ("r1", "r2", "r3")])] \
        + [(c, CCOL[c], [(c, r) for r in (1, 2, 3)]) for c in NPD] + [(c, CCOL[c], [(c, r) for r in (1, 2)]) for c in SNP]
    cvals = [P[("control", k)] for k in ("r1", "r2", "r3")]
    axD.axvspan(min(cvals), max(cvals), color="#bbbbbb", alpha=0.45, lw=0, zorder=1)
    for yi, (name, col, keys) in enumerate(rows):
        y = len(rows) - 1 - yi
        for j, key in enumerate(keys):
            rep = j + 1
            axD.plot(P[key], y, MARK[rep], color=col, ms=4.6, mec="#222222", mew=0.4, zorder=3)
    axD.set_yticks(range(len(rows))); axD.set_yticklabels([r[0] for r in rows][::-1], fontsize=fs.FS["tick"])
    axD.set_xlim(0, 1.04); axD.set_ylim(-0.6, len(rows) - 0.4); axD.tick_params(axis="x", labelsize=fs.FS["tick"])
    axD.set_xlabel("fraction of frames with the ligand to Asn253 distance \u2264 3.5 \u00c5", fontsize=fs.FS["label"])
    axD.legend(handles=[plt.Line2D([], [], marker=MARK[r], ls="", color="#888888", mec="#222222", ms=4.2, label=f"run {r}") for r in (1, 2, 3)]
               + [plt.Rectangle((0, 0), 1, 1, color="#bbbbbb", alpha=0.45, label="control range")],
               fontsize=fs.FS["tick"], frameon=False, loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=4, handlelength=1.0, borderaxespad=0.1)

    for ax, lab in ((axA, "A"), (axB, "B"), (axC, "C"), (axD, "D")):
        fs.panel_label(ax, lab, offset=4 if lab != "D" else 18)
    stem = OUT / ("Figure_S2_md_all_compounds" if ALL else "Figure_7_md_snplib")
    fs.save_figure(fig, stem)

    claims = {"window_ns": TMAX, "all_runs_within_10A_of_site": all10, "max_pocket_distance_A": pocket_max,
              "persistence": {f"{a}_{b}": round(v, 4) for (a, b), v in P.items()},
              "control_range": [round(min(cvals), 4), round(max(cvals), 4)],
              "runs_inside_control_range": [f"{a}_{b}" for (a, b), v in P.items() if a != "control" and min(cvals) <= v <= max(cvals)],
              "runs_below_control_range": [f"{a}_{b}" for (a, b), v in P.items() if a != "control" and v < min(cvals)],
              "anchors": anchors}
    (OUT / (stem.name + "_provenance.json")).write_text(json.dumps({"figure": "Figure 7", "script": "figure7_md_snplib.py",
                                                                         "claims_checked_against_plotted_data": claims,
                                                                         "inputs_sha256": inputs}, indent=1), "utf-8")
    print("inside control range:", claims["runs_inside_control_range"])
    print("below control range :", claims["runs_below_control_range"])
    print("written", stem.with_suffix(".png"))


if __name__ == "__main__":
    main()

"""ADR-AIIM-320, executed. Predictions frozen before this ran.

Pre-registration SHA-256
  7ce8f72966ab63b276c943d1d60cf1d3193a146d16fd2d20dc5fe59b9b448254
Decision log at freeze (text/LF)
  392061d7ba3b5d45b5b0e799e0251ebfd219d4ae6aa47118ff8a46036ab5cd69

Recomputes the published molecular-dynamics figures over a common window,
because the new control replicates run 50 ns and the candidates ran 100 ns.
Nothing is re-simulated: the inputs are the time series the published analyses
already wrote. The summary functions are those of md_14_summarise_persistence.py
(RMSD, threshold ladder, bilayer area) and md_18_recognition_motif.py (N6.55
persistence), with the window as the only change. The same code run at a 100 ns
window is the anchor (W1) and must reproduce the published numbers first.

`10_new_system` is not modified and no published artifact is overwritten.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(r"D:\Agentic VS Workflow")
RUNS = ROOT / "20_new_runs"
CAND = RUNS / "md_candidates_2026-08-22"
CTRL = RUNS / "md_production_2026-08-21" / "04_run" / "analysis_prod"
MOTIF = RUNS / "md_recognition_motif_2026-09-02"
OUT = RUNS / "md_window50_2026-10-07"

PREREG = "7ce8f72966ab63b276c943d1d60cf1d3193a146d16fd2d20dc5fe59b9b448254"
LOG_AT_FREEZE = "392061d7ba3b5d45b5b0e799e0251ebfd219d4ae6aa47118ff8a46036ab5cd69"

CANDIDATES = ["NPLIB01185", "NPLIB01643", "NPLIB03056"]
REPLICATES = [1, 2, 3]
LADDER = [2.0, 3.0, 4.0, 5.0]
HBOND_CUTOFF_A = 3.5
RHO_REPRODUCES, RHO_DOES_NOT = 0.5, 0.1          # ADR-AIIM-161, unchanged
WINDOWS = {"100": 100.0, "50": 50.0}
EPS = 1e-3                                       # xtc times such as 100.0000076

# published values, read from the artifacts that carry them or the article
PUB_CONTROL = {"mean_rmsd_A": 2.59, "sd_rmsd_A": 0.51, "box_percent": -1.05,
               "persistence": 0.974, "n_frames_persistence": 1001}


def read_xvg(p: Path, ncols: int = 2):
    rows = []
    for line in p.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = line.strip()
        if not line or line[0] in "#@":
            continue
        parts = line.split()
        if len(parts) >= ncols:
            try:
                rows.append(tuple(float(x) for x in parts[:ncols]))
            except ValueError:
                pass
    return rows


def window(rows, tmax):
    return [r for r in rows if r[0] <= tmax + EPS]


def summarise_rmsd(rows):
    vals = [v * 10.0 for _, v in rows]            # gmx rms writes nm
    ts = [t for t, _ in rows]
    n = len(vals)
    mean = sum(vals) / n
    sd = (sum((v - mean) ** 2 for v in vals) / (n - 1)) ** 0.5 if n > 1 else 0.0
    imax = max(range(n), key=lambda i: vals[i])
    return {"n_frames": n, "length_ns": round(ts[-1], 1),
            "mean_rmsd_A": round(mean, 2), "sd_rmsd_A": round(sd, 2),
            "max_rmsd_A": round(vals[imax], 2), "max_at_ns": round(ts[imax], 1),
            "final_rmsd_A": round(vals[-1], 2),
            "fraction_below": {f"{c:.1f}": round(sum(1 for v in vals if v < c) / n, 4)
                               for c in LADDER}}


def box_trend(rows):
    if len(rows) < 8:
        return None
    q = len(rows) // 4
    areas = [x * y for _, x, y in rows]           # xy area, as md_14 does
    first = sum(areas[:q]) / q
    last = sum(areas[-q:]) / q
    return {"first_quarter_nm2": round(first, 2), "last_quarter_nm2": round(last, 2),
            "percent_change": round((last - first) / first * 100.0, 2)}


def persistence(rows, cutoff=HBOND_CUTOFF_A):
    sel = [(t, v * 10.0) for t, v in rows]        # nm -> angstrom
    hits = sum(1 for _, v in sel if v <= cutoff)
    vals = sorted(v for _, v in sel)
    return {"n_frames": len(sel), "n_satisfying": hits,
            "persistence": round(hits / len(sel), 4),
            "median_A": round(vals[len(vals) // 2], 2),
            "min_A": round(vals[0], 2)}


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def outcome(rho):
    return ("REPRODUCES" if rho >= RHO_REPRODUCES else
            "DOES_NOT_REPRODUCE" if rho < RHO_DOES_NOT else "PARTIAL")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    inputs = {}

    def need(p):
        if not p.is_file():
            raise SystemExit(f"missing input: {p}")
        inputs[str(p.relative_to(ROOT)).replace("\\", "/")] = sha(p)
        return p

    # ---- raw series, read once -------------------------------------------
    raw = {"control": {
        "rmsd": read_xvg(need(CTRL / "ligand_rmsd.xvg")),
        "box": read_xvg(need(CTRL / "box.xvg"), ncols=3),
        "n655": read_xvg(need(MOTIF / "control_ZM241385" / "N6_55.xvg"))}}
    for cid in CANDIDATES:
        for r in REPLICATES:
            d = CAND / cid / "04_run" / f"analysis_prod_r{r}"
            raw[f"{cid}_r{r}"] = {
                "rmsd": read_xvg(need(d / "ligand_rmsd.xvg")),
                "box": read_xvg(need(d / "box.xvg"), ncols=3),
                "n655": read_xvg(need(MOTIF / f"{cid}_r{r}" / "N6_55.xvg"))}

    res = {}
    for wname, tmax in WINDOWS.items():
        w = {}
        for tag, s in raw.items():
            # RMSD and N6.55 series are in ns; box.xvg (gmx energy) is in ps.
            # The first run of this script cut box.xvg at "100" and "50" as if
            # they were ns, kept 11 and 6 rows, and the anchor (W1) caught it:
            # control bilayer area came out +0.15 % against the published -1.05 %.
            rm, bx, nn = (window(s["rmsd"], tmax), window(s["box"], tmax * 1000.0),
                          window(s["n655"], tmax))
            w[tag] = {"rmsd": summarise_rmsd(rm), "box": box_trend(bx),
                      "n655": persistence(nn)}
        ctrl_p = w["control"]["n655"]["persistence"]
        per = {}
        for cid in CANDIDATES:
            reps = {f"r{r}": w[f"{cid}_r{r}"] for r in REPLICATES}
            means = [reps[f"r{r}"]["rmsd"]["mean_rmsd_A"] for r in REPLICATES]
            ps = [reps[f"r{r}"]["n655"]["persistence"] for r in REPLICATES]
            mp = sum(ps) / len(ps)
            rho = mp / ctrl_p
            per[cid] = {"replicates": reps,
                        "mean_of_replicate_means_A": round(sum(means) / 3, 2),
                        "spread_of_replicate_means_A": round(max(means) - min(means), 2),
                        "replicate_persistence": ps,
                        "mean_persistence": round(mp, 4),
                        "rho_vs_control": round(rho, 3), "outcome": outcome(rho),
                        "rmsd_reportable": cid != "NPLIB01185"}
        res[wname] = {"tmax_ns": tmax, "control": w["control"], "candidates": per}

    a100, a50 = res["100"], res["50"]

    # ---- the published values the anchor is checked against --------------
    pub_p = json.loads((RUNS / "md_pose_persistence_2026-08-31" /
                        "md_pose_persistence.json").read_text("utf-8"))
    pub_m = json.loads((MOTIF / "recognition_motif.json").read_text("utf-8"))
    anchor = {}
    c100 = a100["control"]
    anchor["control_rmsd_mean"] = (c100["rmsd"]["mean_rmsd_A"], PUB_CONTROL["mean_rmsd_A"])
    anchor["control_rmsd_sd"] = (c100["rmsd"]["sd_rmsd_A"], PUB_CONTROL["sd_rmsd_A"])
    anchor["control_box_percent"] = (c100["box"]["percent_change"], PUB_CONTROL["box_percent"])
    anchor["control_persistence"] = (c100["n655"]["persistence"], PUB_CONTROL["persistence"])
    anchor["control_persistence_frames"] = (c100["n655"]["n_frames"], PUB_CONTROL["n_frames_persistence"])
    for cid in CANDIDATES:
        anchor[f"{cid}_rmsd_mean_of_means"] = (
            a100["candidates"][cid]["mean_of_replicate_means_A"],
            pub_p["candidates"][cid]["mean_of_replicate_means_A"])
        anchor[f"{cid}_mean_persistence"] = (
            a100["candidates"][cid]["mean_persistence"],
            pub_m["outcomes"][cid]["mean_persistence"])
    anchor_ok = {k: (abs(a - b) < 1e-9) for k, (a, b) in anchor.items()}
    w1 = all(anchor_ok.values())

    # ---- predictions W2-W5 ------------------------------------------------
    d_ctrl = abs(a50["control"]["rmsd"]["mean_rmsd_A"] - c100["rmsd"]["mean_rmsd_A"])
    w3_rows = {cid: (a50["candidates"][cid]["mean_of_replicate_means_A"],
                     a100["candidates"][cid]["mean_of_replicate_means_A"])
               for cid in ("NPLIB01643", "NPLIB03056")}
    w3 = all(a <= b for a, b in w3_rows.values())
    outs50 = {cid: a50["candidates"][cid]["outcome"] for cid in CANDIDATES}
    w4 = (all(o == "PARTIAL" for o in outs50.values())
          and a50["control"]["n655"]["persistence"] >= 0.95)
    order = sorted(CANDIDATES, key=lambda c: a50["candidates"][c]["mean_persistence"])
    w5 = order == ["NPLIB01185", "NPLIB01643", "NPLIB03056"]

    grading = {
        "W1_anchor_reproduces_published_at_100ns": {
            "checks": {k: {"computed": a, "published": b, "match": anchor_ok[k]}
                       for k, (a, b) in anchor.items()},
            "verdict": "MET" if w1 else "MISSED"},
        "W2_control_stationary_within_0.3A": {
            "mean_50_A": a50["control"]["rmsd"]["mean_rmsd_A"],
            "mean_100_A": c100["rmsd"]["mean_rmsd_A"], "abs_difference_A": round(d_ctrl, 3),
            "verdict": "MET" if d_ctrl < 0.3 else "MISSED"},
        "W3_drift_accumulates": {
            "mean_50_vs_100_A": {k: list(v) for k, v in w3_rows.items()},
            "verdict": "MET" if w3 else "MISSED"},
        "W4_verdict_unchanged_PARTIAL": {
            "outcomes_50ns": outs50,
            "rho_50ns": {c: a50["candidates"][c]["rho_vs_control"] for c in CANDIDATES},
            "control_persistence_50ns": a50["control"]["n655"]["persistence"],
            "verdict": "MET" if w4 else "MISSED"},
        "W5_ordering_unchanged": {
            "order_50ns_ascending": order, "verdict": "MET" if w5 else "MISSED"}}

    out = {"adr": "ADR-AIIM-320", "prereg_sha256": PREREG,
           "log_sha256_at_freeze": LOG_AT_FREEZE,
           "window_rule": f"frames with t <= tmax + {EPS} ns kept; nothing re-simulated",
           "frame_spacing_ns": {"control": 0.1, "candidates": 0.05},
           "windows": res, "grading": grading, "inputs_sha256": inputs}
    (OUT / "md_window50.json").write_text(json.dumps(out, indent=1), encoding="utf-8")

    # ---- console ----------------------------------------------------------
    print("=" * 84)
    for k, g in grading.items():
        print(f"{k:50s} {g['verdict']}")
    print("=" * 84)
    print("\nanchor (same code, 100 ns window) vs published:")
    for k, (a, b) in anchor.items():
        print(f"  {k:36s} computed {a!s:>10}  published {b!s:>10}  {'OK' if anchor_ok[k] else 'DIFFERS'}")
    print(f"\n{'':22s}{'100 ns':>12s}{'50 ns':>12s}")
    print(f"{'control mean RMSD (A)':22s}{c100['rmsd']['mean_rmsd_A']:>12.2f}{a50['control']['rmsd']['mean_rmsd_A']:>12.2f}")
    print(f"{'control sd (A)':22s}{c100['rmsd']['sd_rmsd_A']:>12.2f}{a50['control']['rmsd']['sd_rmsd_A']:>12.2f}")
    for lad in ("3.0", "4.0"):
        print(f"{'control frac < '+lad+' A':22s}{c100['rmsd']['fraction_below'][lad]:>12.3f}{a50['control']['rmsd']['fraction_below'][lad]:>12.3f}")
    print(f"{'control bilayer area %':22s}{c100['box']['percent_change']:>12.2f}{a50['control']['box']['percent_change']:>12.2f}")
    print(f"{'control N6.55 persist.':22s}{c100['n655']['persistence']:>12.3f}{a50['control']['n655']['persistence']:>12.3f}"
          f"   frames {c100['n655']['n_frames']} / {a50['control']['n655']['n_frames']}")
    print(f"\n{'compound':12s}{'RMSD100':>9s}{'RMSD50':>9s}{'P100':>8s}{'P50':>8s}{'rho50':>7s}  outcome50   per-replicate persistence at 50 ns")
    for cid in CANDIDATES:
        c1, c5 = a100["candidates"][cid], a50["candidates"][cid]
        flag = "" if c5["rmsd_reportable"] else "  (RMSD not reportable)"
        print(f"{cid:12s}{c1['mean_of_replicate_means_A']:>9.2f}{c5['mean_of_replicate_means_A']:>9.2f}"
              f"{c1['mean_persistence']:>8.3f}{c5['mean_persistence']:>8.3f}{c5['rho_vs_control']:>7.2f}"
              f"  {c5['outcome']:10s} {c5['replicate_persistence']}{flag}")
    print(f"\nwritten -> {OUT / 'md_window50.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

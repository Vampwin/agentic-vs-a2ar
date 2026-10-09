"""ADR-AIIM-324, executed. Predictions frozen before any replicate was analysed.

Pre-registration SHA-256
  57f0cf7cb0b8eef7316b3337a5cd77f5eba27c8cec7b5fe2be40dd46a69f871e
Decision log at freeze (text/LF)
  5b36b654a07fe9fa8425632a3aa78b7ab9c317e666b8c194f457adc04ecddb80

Summaries of control replicates 2 and 3 over the first 50 ns, from the time series
the unchanged original pipeline wrote (run_control_analysis.sh). The summary
functions are those of md_window50.py (ADR-AIIM-320), which are md_14's and md_18's
with the window as the only change. The anchor is the same pipeline on a COPY of
the existing control; if it does not reproduce the published numbers nothing about
the replicates may be reported.

Replicate 3 is graded only if its trajectory is complete (prod_r3.gro exists and the
series reach 50 ns). A partial trajectory is never summarised.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

ROOT = Path(r"D:\Agentic VS Workflow")
RUNS = ROOT / "20_new_runs"
REP = RUNS / "md_control_replicates_2026-10-06"
ANCH = RUNS / "md_control_anchor_2026-10-07"
WIN = RUNS / "md_window50_2026-10-07"
OUT = RUNS / "md_control_replicates_window50_2026-10-07"
PREREG = "57f0cf7cb0b8eef7316b3337a5cd77f5eba27c8cec7b5fe2be40dd46a69f871e"
LOG_AT_FREEZE = "5b36b654a07fe9fa8425632a3aa78b7ab9c317e666b8c194f457adc04ecddb80"
TMAX = 50.0

# Published values, read from the artifact ADR-AIIM-321 wrote rather than typed:
# the first execution of this script hard-coded 0.966 while the stored value is
# 0.9661 and reported the anchor as missed on that rounding alone.
_pub = json.loads((RUNS / "md_window50_2026-10-07" / "md_window50.json").read_text("utf-8"))["windows"]
PUB = {w: {"rmsd_mean": _pub[w]["control"]["rmsd"]["mean_rmsd_A"],
           "rmsd_sd": _pub[w]["control"]["rmsd"]["sd_rmsd_A"],
           "box": _pub[w]["control"]["box"]["percent_change"],
           "persist": _pub[w]["control"]["n655"]["persistence"]} for w in ("100", "50")}

spec = importlib.util.spec_from_file_location("w50", WIN / "analysis" / "md_window50.py")
w50 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(w50)            # defines functions only; main() is guarded


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def summarise(adir: Path, tmax: float, inputs: dict):
    def need(p):
        inputs[str(p.relative_to(ROOT)).replace("\\", "/")] = sha(p)
        return p
    rm = w50.window(w50.read_xvg(need(adir / "ligand_rmsd.xvg")), tmax)
    bx = w50.window(w50.read_xvg(need(adir / "box.xvg"), ncols=3), tmax * 1000.0)   # ps
    nn = w50.window(w50.read_xvg(need(adir / "N6_55.xvg")), tmax)
    return {"rmsd": w50.summarise_rmsd(rm), "box": w50.box_trend(bx),
            "n655": w50.persistence(nn)}


def complete(adir: Path, gro: Path) -> bool:
    if not gro.is_file():
        return False
    try:
        t = w50.read_xvg(adir / "ligand_rmsd.xvg")
        return bool(t) and t[-1][0] >= TMAX - 0.2
    except FileNotFoundError:
        return False


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    inputs: dict = {}
    existing = {"100": summarise(ANCH / "analysis_prod", 100.0, inputs),
                "50": summarise(ANCH / "analysis_prod", 50.0, inputs)}

    # ---- anchor ------------------------------------------------------------
    chk = {}
    for w in ("100", "50"):
        s, p = existing[w], PUB[w]
        chk[f"{w}ns_rmsd_mean"] = (s["rmsd"]["mean_rmsd_A"], p["rmsd_mean"])
        chk[f"{w}ns_rmsd_sd"] = (s["rmsd"]["sd_rmsd_A"], p["rmsd_sd"])
        chk[f"{w}ns_box_percent"] = (s["box"]["percent_change"], p["box"])
        chk[f"{w}ns_persistence"] = (s["n655"]["persistence"], p["persist"])
    anchor_ok = {k: abs(a - b) < 1e-9 for k, (a, b) in chk.items()}
    A = all(anchor_ok.values())

    # ---- replicates --------------------------------------------------------
    reps, status = {}, {}
    for name, deffnm in (("r2", "prod_r2"), ("r3", "prod_r3")):
        adir = REP / f"analysis_{deffnm}"
        done = complete(adir, REP / f"{deffnm}.gro")
        status[name] = "complete" if done else "not complete: not summarised"
        if done:
            reps[name] = summarise(adir, TMAX, inputs)

    ex50 = existing["50"]
    ids = ["existing"] + list(reps)
    means = {"existing": ex50["rmsd"]["mean_rmsd_A"],
             **{k: v["rmsd"]["mean_rmsd_A"] for k, v in reps.items()}}
    persist = {"existing": ex50["n655"]["persistence"],
               **{k: v["n655"]["persistence"] for k, v in reps.items()}}

    def verdict(cond, need_r3=False):
        if need_r3 and "r3" not in reps:
            return "PENDING (replicate 3 incomplete)"
        return "MET" if cond else "MISSED"

    r1 = {k: round(v["rmsd"]["fraction_below"]["4.0"], 4) for k, v in reps.items()}
    r5 = {k: v["box"]["percent_change"] for k, v in reps.items()}
    cand = json.loads((WIN / "md_window50.json").read_text("utf-8"))["windows"]["50"]["candidates"]
    ctrl_mean_p = sum(persist.values()) / len(persist)
    ratios = {c: round(cand[c]["mean_persistence"] / ctrl_mean_p, 3) for c in cand}
    spread = max(means.values()) - min(means.values())
    lo_hi = [means[k] for k in reps]
    grading = {
        "A_anchor_reproduces_published": {
            "checks": {k: {"computed": a, "published": b, "match": anchor_ok[k]}
                       for k, (a, b) in chk.items()},
            "verdict": "MET" if A else "MISSED"},
        "R1_each_replicate_below_4A_in_95pct_of_frames": {
            "fraction_below_4A": r1,
            "verdict": verdict(all(v >= 0.95 for v in r1.values()) and len(r1) == 2, True)},
        "R2_means_span_under_1A": {
            "means_A": means, "span_A": round(spread, 3),
            "verdict": verdict(spread < 1.0, True)},
        "R3_existing_run_inside_replicate_range": {
            "existing_A": means["existing"],
            "replicate_range_A": [min(lo_hi), max(lo_hi)] if lo_hi else None,
            "verdict": verdict(bool(lo_hi) and min(lo_hi) <= means["existing"] <= max(lo_hi), True)},
        "R4_persistence_at_least_0.90_each": {
            "persistence": {k: v for k, v in persist.items() if k != "existing"},
            "verdict": verdict(all(persist[k] >= 0.90 for k in reps) and len(reps) == 2, True)},
        "R5_bilayer_area_within_3pct_each": {
            "percent_change": r5,
            "verdict": verdict(all(abs(v) <= 3.0 for v in r5.values()) and len(r5) == 2, True)},
        "R6_candidates_stay_PARTIAL_against_mean_control": {
            "control_mean_persistence": round(ctrl_mean_p, 4), "runs_averaged": ids,
            "candidate_over_control": ratios,
            "verdict": verdict(all(0.1 <= r < 0.5 for r in ratios.values()), True)}}

    out = {"adr": "ADR-AIIM-324", "prereg_sha256": PREREG,
           "log_sha256_at_freeze": LOG_AT_FREEZE, "window_ns": TMAX,
           "status": status, "existing_control": existing, "replicates": reps,
           "grading": grading, "inputs_sha256": inputs}
    (OUT / "md_control_replicates_window50.json").write_text(json.dumps(out, indent=1), encoding="utf-8")

    print("=" * 88)
    for k, g in grading.items():
        print(f"{k:56s} {g['verdict']}")
    print("=" * 88)
    print("status:", status)
    print("\nanchor (existing control through the same pipeline) vs published:")
    for k, (a, b) in chk.items():
        print(f"  {k:22s} computed {a:>8}  published {b:>8}  {'OK' if anchor_ok[k] else 'DIFFERS'}")
    print(f"\n{'run (first 50 ns)':20s}{'frames':>7s}{'RMSD':>8s}{'sd':>6s}{'<3A%':>6s}{'<4A%':>6s}{'area%':>7s}{'N6.55':>8s}")
    rows = [("existing (continuation)", ex50)] + [(k, v) for k, v in reps.items()]
    for k, v in rows:
        print(f"{k:20s}{v['rmsd']['n_frames']:>7d}{v['rmsd']['mean_rmsd_A']:>8.2f}{v['rmsd']['sd_rmsd_A']:>6.2f}"
              f"{v['rmsd']['fraction_below']['3.0']*100:>6.0f}{v['rmsd']['fraction_below']['4.0']*100:>6.0f}"
              f"{v['box']['percent_change']:>7.2f}{v['n655']['persistence']:>8.3f}")
    print(f"\ncontrol persistence (mean of {ids}): {ctrl_mean_p:.4f}")
    print("candidate persistence / control:", ratios)
    print(f"\nwritten -> {OUT / 'md_control_replicates_window50.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

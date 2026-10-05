"""analyse_spikein_recovery.py -- the M2.1 recovery funnel, written before the data exists.

Executes the analysis pre-registered in ADR-AIIM-229. **Written while the
docking was still running and before any `vina_score` of it had been read**, for
the reason ADR-AIIM-228 gives: an analysis written after the results are visible
is not the analysis that was pre-registered, whatever its header claims. The
M5.1 precedent held under that discipline; the M2.2 arm got the weaker guarantee
and ADR-AIIM-241 says so rather than glossing it.

`wilson` is imported from `analyse_r4_funnel_recall.py` rather than restated, so
the interval convention is this project's own.

### What the pre-registration fixed

| | |
|---|---|
| Metric | count surviving each stage, Wilson 95 % interval, beside the background rate |
| Arms | primary **50**, pre-registered replication **55**, combined **105** |
| P1 | the AD gate and the ligand-based arm admit most surviving spike-ins |
| P2 | **the structure-based rank step loses the majority** |
| P3 | if P1 and P2 hold, the loss is localised to docking rank and 3.5's undifferentiated negative becomes a diagnosis |

Falsification, from ADR-AIIM-229: if most spike-ins survive every stage and reach
the shortlist, **P2 is wrong** and the funnel recovers known actives end to end.
If the AD gate rather than docking is the largest loss, **P1 is wrong** and
ADR-AIIM-210's 9-of-10 does not generalise past clinical-stage chemistry.

### The one thing ADR-AIIM-229 left underspecified, fixed here before any score

"Count surviving each stage" does not define what surviving **stage 4** means:
docking produces a score for every compound, not a pass or a fail. A cut-off is
needed and it must not be chosen after seeing where the spike-ins land.

**Primary cut-off: the top 500 by docking score**, because that is the cut the
production funnel used -- 3.6 reports "of the 500 highest-ranked compounds by
fast-pass score, 497 completed full-protocol redocking" -- so recovery at 500 is
the number that speaks to the campaign this project actually ran.

The full recovery curve at 100 / 250 / 500 / 1,000 / 2,500 is reported beside it,
so no reader depends on the cut-off having been well chosen. **All of these were
written down before a single score of this run was read.**

### Stages 5 and 6 are not analysed here

The conformal predictor, the AD gate and the dual-criterion rule need
`run_g6_on_np_library` executed against the spiked shortlist, which is a further
run rather than an analysis. **P1 therefore cannot be graded by this script** and
is reported as pending, not as met or missed. P2 is answerable now and is graded.

Output -> 20_new_runs/spike_in_dock_2026-09-23/analysis/
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kdss_paths import under  # noqa: E402

# This project's own Wilson convention, imported rather than restated.
from analyse_r4_funnel_recall import wilson  # noqa: E402

RUN = Path(under("20_new_runs", "spike_in_dock_2026-09-23"))
SCORES = RUN / "scores.csv"
FUNNEL = Path(under("20_new_runs", "spike_in_recovery_2026-09-23"))
SHORTLIST = FUNNEL / "np_library_docking_shortlist.csv"
STAGE2 = FUNNEL / "np_library_stage2_2d_ranked.csv"
POOL = Path(under("20_new_runs", "spikein_panel_2026-09-22",
                  "spiked_pool_25709.csv"))
PRIMARY50 = Path(under("20_new_runs", "spikein_panel_2026-09-22",
                       "spikein_panel_primary_50.csv"))

N_SHORTLIST = 5000
N_SPIKED = 105
CUTOFFS = (100, 250, 500, 1000, 2500)
PRIMARY_CUTOFF = 500

# ADR-AIIM-229's derived table. A run that disagrees has a fault, not a finding.
DERIVED = {"stage1": {"panel": 102, "primary": 50},
           "stage2": {"panel": 94, "primary": 47}}


def rate(k: int, n: int) -> dict:
    lo, hi = wilson(k, n)
    return {"k": int(k), "n": int(n),
            "rate": round(k / n, 4) if n else None,
            "wilson95": [round(lo, 4), round(hi, 4)]}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(RUN / "analysis"))
    ap.add_argument("--scores", default=str(SCORES))
    ap.add_argument("--allow-partial", action="store_true",
                    help="analyse an unfinished run; the output is marked "
                         "PROVISIONAL and must not be reported")
    args = ap.parse_args()
    out = Path(args.out)

    pool = pd.read_csv(POOL)
    spikes = set(pool.loc[pool.spike_in.astype(bool),
                          "structure_inchikey"].astype(str))
    if len(spikes) != N_SPIKED:
        print(f"ABORT: the frozen pool carries {len(spikes)} spike-ins, "
              f"expected {N_SPIKED}.")
        return 4
    try:
        prim = pd.read_csv(PRIMARY50)
        pcol = next((c for c in ("structure_inchikey", "skeleton_inchikey",
                                 "compound_id") if c in prim.columns), None)
        primary = set(prim[pcol].astype(str)) if pcol else set()
    except Exception:                                        # noqa: BLE001
        primary = set()

    # ---- stages 1-3, the implementation check ---------------------------
    st2 = pd.read_csv(STAGE2)
    sl = pd.read_csv(SHORTLIST)
    s2 = int(st2.structure_inchikey.astype(str).isin(spikes).sum())
    s3 = int(sl.structure_inchikey.astype(str).isin(spikes).sum())
    impl_ok = s2 == DERIVED["stage2"]["panel"]

    # ---- stage 4 ---------------------------------------------------------
    d = pd.read_csv(args.scores)
    scored = d[d.vina_score.notna() & (d.vina_score != "")].copy()
    if len(d) < N_SHORTLIST and not args.allow_partial:
        print(f"REFUSING: {len(d)} of {N_SHORTLIST} shortlist compounds are "
              f"scored. Re-run when the docking completes, or pass "
              f"--allow-partial and accept the label.")
        return 1
    scored["vina_score"] = scored.vina_score.astype(float)
    bad = scored[scored.n_conformers_docked != 1]
    if len(bad):
        print(f"ABORT: {len(bad)} rows docked with a conformer count other "
              f"than 1; not comparable with the other benchmarks.")
        return 2

    key = "structure_inchikey"
    scored[key] = scored[key].astype(str)
    scored["is_spike"] = scored[key].isin(spikes)
    scored["is_primary"] = scored[key].isin(primary) if primary else False

    # more negative is better
    scored = scored.sort_values("vina_score", kind="mergesort").reset_index(drop=True)
    scored["dock_rank"] = range(1, len(scored) + 1)

    n_spikes_docked = int(scored.is_spike.sum())
    curve = {}
    for c in CUTOFFS:
        top = scored.head(c)
        k = int(top.is_spike.sum())
        curve[str(c)] = {
            **rate(k, n_spikes_docked),
            "background_rate": round((c - k) / (len(scored) - n_spikes_docked), 4)
            if len(scored) > n_spikes_docked else None,
            "enrichment_over_background": (
                round((k / n_spikes_docked) / (c / len(scored)), 3)
                if n_spikes_docked and len(scored) else None),
        }

    p = curve[str(PRIMARY_CUTOFF)]
    # P2 asks whether the rank step loses the MAJORITY: fewer than half the
    # docked spike-ins reach the production cut.
    p2_met = bool(p["k"] < n_spikes_docked / 2)

    summary = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "executes": "M2.1 step 4",
        "pre_registration": "ADR-AIIM-229; step 4 authorised in ADR-AIIM-245",
        "analysis_written": "while the docking was still running, before any "
                            "vina_score of this run had been read",
        "provisional": bool(len(d) < N_SHORTLIST),
        "n_shortlist": int(len(d)),
        "n_scored": int(len(scored)),
        "n_failed": int(len(d) - len(scored)),
        "recovery_funnel": {
            "spiked_into_pool": rate(N_SPIKED, N_SPIKED),
            "stage2_top8000": rate(s2, N_SPIKED),
            "stage3_shortlist5000": rate(s3, N_SPIKED),
            "stage4_docked": rate(n_spikes_docked, N_SPIKED),
        },
        "implementation_check": {
            "stage2_expected": DERIVED["stage2"]["panel"],
            "stage2_observed": s2,
            "agrees": impl_ok,
            "note": "ADR-AIIM-229: a run that disagrees with the derived "
                    "stage-1/2 table has a fault, not a finding",
        },
        "stage4_cutoff": {
            "primary": PRIMARY_CUTOFF,
            "why": "the cut the production funnel used (section 3.6: the 500 "
                   "highest-ranked compounds went to full-protocol redocking)",
            "fixed_before_any_score_was_read": True,
            "curve": curve,
        },
        "predictions": {
            "P1": {"statement": "the AD gate and the ligand-based arm admit "
                                "most surviving spike-ins",
                   "status": "PENDING",
                   "why": "stages 5 and 6 need run_g6_on_np_library executed "
                          "against the spiked shortlist; that is a further run, "
                          "not an analysis, and this script does not guess it"},
            "P2": {"statement": "the structure-based rank step loses the majority",
                   "met": p2_met,
                   "observed": {"spike_ins_in_top_500": p["k"],
                                "of_docked": n_spikes_docked,
                                "rate": p["rate"], "wilson95": p["wilson95"]}},
            "P3": {"statement": "if P1 and P2 hold the loss is localised to "
                                "docking rank",
                   "status": "PENDING on P1"},
        },
        "not_claimed": [
            "that any panel compound is active in any assay run by this "
            "project; no compound reported anywhere here has been assayed",
            "that the LOTUS background is a property-matched decoy set -- this "
            "is a recovery control, not an enrichment benchmark (ADR-AIIM-229)",
            "that recovery through stages 1 and 2 is evidence the pipeline "
            "works: the stage-2 screen ranks by similarity to known A2A "
            "ligands and the panel are known A2A ligands",
        ],
    }

    out.mkdir(parents=True, exist_ok=True)
    scored.to_csv(out / "spikein_docked_ranked.csv", index=False)
    with open(out / "spikein_recovery.json", "w", encoding="utf-8",
              newline="\n") as f:
        json.dump(summary, f, indent=1, ensure_ascii=False)
        f.write("\n")

    print("=" * 78)
    print("M2.1 RECOVERY FUNNEL"
          + ("   ** PROVISIONAL **" if summary["provisional"] else ""))
    print("=" * 78)
    print(f"  {'stage':34}{'k / n':>12}{'rate':>9}{'Wilson 95 %':>20}")
    print("  " + "-" * 74)
    for k, v in summary["recovery_funnel"].items():
        ci = f"[{v['wilson95'][0]:.3f}, {v['wilson95'][1]:.3f}]"
        print(f"  {k:34}{str(v['k']) + ' / ' + str(v['n']):>12}"
              f"{v['rate']:>9.3f}{ci:>20}")
    print()
    print(f"  implementation check, stage 2: observed {s2}, "
          f"expected {DERIVED['stage2']['panel']} -> "
          f"{'AGREES' if impl_ok else 'DISAGREES -- this is a fault'}")
    print()
    print(f"  {'docking cut-off':>18}{'spike-ins':>12}{'rate':>9}"
          f"{'enrichment':>13}")
    for c in CUTOFFS:
        v = curve[str(c)]
        star = "  <- primary" if c == PRIMARY_CUTOFF else ""
        print(f"  {'top ' + str(c):>18}{v['k']:>12}{v['rate']:>9.3f}"
              f"{str(v['enrichment_over_background']):>13}{star}")
    print()
    print(f"  P2 (rank step loses the majority): "
          f"{'MET' if p2_met else 'MISSED'}")
    print("  P1 and P3: PENDING -- stages 5 and 6 are a further run")
    print(f"\nwritten -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

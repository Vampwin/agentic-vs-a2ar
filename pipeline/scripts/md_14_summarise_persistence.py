"""md_14_summarise_persistence.py -- the nine candidate replicates, in one place.

`md_11_analyse_pose_persistence.sh` writes one directory per replicate. This
collects them, adds the control for comparison, and emits the artifact Section 3
reports and G8 audits (ADR-AIIM-146, decision 2).

The three analysis choices were fixed on the control before any candidate ran and
are not re-decided here: fit on receptor C-alpha and measure the ligand;
reference the pre-equilibration structure, which holds the docked pose for a
candidate and the crystallographic pose for the control; and report a ladder of
thresholds rather than one pass/fail, because the 2.0 A redocking criterion
compares a prediction to a cryogenic crystal structure and is far too strict for
a 310 K trajectory -- on the control only 17 % of frames sit below 2.0 A while
100 % are below 4.0 A.

No verdict is computed. The question "did the pose persist" has no threshold in
this project's pre-registration, and inventing one after seeing three compounds
behave differently is exactly the move the rest of this revision has been undoing.
The numbers are reported and the reader draws the line.

Usage:
    python md_14_summarise_persistence.py [--out DIR]
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from kdss_paths import project_root, resolve, under, dev_dir


def _root() -> Path:
    # Delegates to kdss_paths, which finds whichever layout this file is inside
    # -- the development tree or a deposit -- by walking up from its own
    # location. It previously tried two absolute paths and exited if neither
    # existed, which on the author's machine meant a copy of the deposit read
    # and wrote the development tree instead of itself.
    return project_root()


ROOT = _root()
RUNS = dev_dir("20_new_runs")
CAND = under("20_new_runs", "md_candidates_2026-08-22")
CANDIDATES = ["NPLIB01185", "NPLIB01643", "NPLIB03056"]
REPLICATES = [1, 2, 3]
LADDER = [2.0, 3.0, 4.0, 5.0]


def read_xvg(p: Path, ncols: int = 2) -> list[tuple[float, ...]]:
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


def summarise(xvg: Path) -> dict:
    rows = read_xvg(xvg)
    if not rows:
        raise SystemExit(f"no data in {xvg}")
    # gmx rms writes nm; every number this project reports is in angstrom.
    # The TIME column is already nanoseconds -- the analysis script passes
    # `-tu ns` and the xvg header says so. Dividing it by 1000 recorded a
    # 100 ns run as 0.1 ns in the first version of this artifact. The RMSD
    # values were unaffected and no manuscript number read the time fields,
    # but a deposited artifact said a run was a thousandth of its length.
    # Found only when the trajectories were plotted and every line collapsed
    # onto x = 0.
    vals = [v * 10.0 for _, v in rows]
    ts = [t for t, _ in rows]
    n = len(vals)
    mean = sum(vals) / n
    sd = (sum((v - mean) ** 2 for v in vals) / (n - 1)) ** 0.5 if n > 1 else 0.0
    imax = max(range(n), key=lambda i: vals[i])
    return {
        "n_frames": n,
        "length_ns": round(ts[-1], 1),
        "mean_rmsd_A": round(mean, 2),
        "sd_rmsd_A": round(sd, 2),
        "max_rmsd_A": round(vals[imax], 2),
        "max_at_ns": round(ts[imax], 1),
        "final_rmsd_A": round(vals[-1], 2),
        "fraction_below": {f"{c:.1f}": round(sum(1 for v in vals if v < c) / n, 4)
                           for c in LADDER},
    }


def box_trend(analysis_dir: Path) -> dict | None:
    box = analysis_dir / "box.xvg"
    if not box.is_file():
        return None
    rows = read_xvg(box, ncols=3)
    if len(rows) < 8:
        return None
    q = len(rows) // 4
    # box.xvg holds time, Box-X, Box-Y. The bilayer quantity is the xy AREA,
    # X * Y, which is what md_box_trend.awk reports and what §3.6 quotes. The
    # first version of this function averaged column 2 alone and produced 1.69 %
    # where the shell check said 3.35 %, under a comment asserting the two
    # agreed by construction. They measured different things; G8 caught it.
    areas = [x * y for _, x, y in rows]
    first = sum(areas[:q]) / q
    last = sum(areas[-q:]) / q
    return {"first_quarter_nm2": round(first, 2), "last_quarter_nm2": round(last, 2),
            "percent_change": round((last - first) / first * 100.0, 2)}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(under("20_new_runs", "md_pose_persistence_2026-08-31")))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    per_compound = {}
    for cid in CANDIDATES:
        reps = {}
        for r in REPLICATES:
            d = CAND / cid / "04_run" / f"analysis_prod_r{r}"
            xvg = d / "ligand_rmsd.xvg"
            if not xvg.is_file():
                raise SystemExit(f"missing analysis: {xvg}")
            s = summarise(xvg)
            bt = box_trend(d)
            if bt:
                s["bilayer_box_trend"] = bt
            reps[f"r{r}"] = s
        means = [reps[f"r{r}"]["mean_rmsd_A"] for r in REPLICATES]
        finals = [reps[f"r{r}"]["final_rmsd_A"] for r in REPLICATES]
        per_compound[cid] = {
            "replicates": reps,
            "mean_of_replicate_means_A": round(sum(means) / len(means), 2),
            "spread_of_replicate_means_A": round(max(means) - min(means), 2),
            "final_rmsd_range_A": [min(finals), max(finals)],
            "n_replicates": len(REPLICATES),
        }

    report = {
        "question": "Does each candidate's docked pose stay where it was placed, "
                    "over three 100 ns replicates in an explicit membrane?",
        "authority": "ADR-AIIM-146 (decision 2), ADR-AIIM-148",
        "method": {
            "fit": "receptor C-alpha",
            "measured": "ligand heavy atoms",
            "reference": "the pre-equilibration structure (em.tpr), which holds "
                         "the docked pose -- not the production start, which "
                         "would measure drift since equilibration and understate "
                         "deviation by about 0.3 A",
            "replicates": "three per compound, 100 ns each, differing only in "
                          "the velocity seed",
            "thresholds": "a ladder, not a verdict. Fixed on the control before "
                          "any candidate ran.",
            "note": "All three choices are md_11_analyse_pose_persistence.sh's "
                    "and predate every candidate trajectory.",
        },
        "control_for_comparison": {
            "compound": "ZM241385, crystallographic pose, 4EIY",
            "n_replicates": 1,
            "mean_rmsd_A": 2.59,
            "sd_rmsd_A": 0.51,
            "fraction_below_4A": 1.00,
            "source": "20_new_runs/md_production_2026-08-21, reported in §2.9",
        },
        "candidates": per_compound,
        "no_verdict": "This project pre-registered no pose-persistence threshold. "
                      "Choosing one now, after seeing that the three candidates "
                      "behave differently, would be the post-hoc move the rest of "
                      "this revision has been removing. The numbers stand as "
                      "measured.",
    }

    dest = out / "md_pose_persistence.json"
    dest.write_text(json.dumps(report, indent=2), encoding="utf-8")

    print(f"{'compound':<12}{'rep':<5}{'mean':>8}{'sd':>7}{'max':>8}{'final':>8}"
          f"{'<4A':>7}{'<5A':>7}")
    for cid, c in per_compound.items():
        for r in REPLICATES:
            s = c["replicates"][f"r{r}"]
            print(f"{cid:<12}r{r:<4}{s['mean_rmsd_A']:>8.2f}{s['sd_rmsd_A']:>7.2f}"
                  f"{s['max_rmsd_A']:>8.2f}{s['final_rmsd_A']:>8.2f}"
                  f"{s['fraction_below']['4.0'] * 100:>6.0f}%"
                  f"{s['fraction_below']['5.0'] * 100:>6.0f}%")
        print(f"{'':<12}{'mean':<5}{c['mean_of_replicate_means_A']:>8.2f}"
              f"   spread {c['spread_of_replicate_means_A']:.2f} A"
              f"   final {c['final_rmsd_range_A'][0]:.1f}-"
              f"{c['final_rmsd_range_A'][1]:.1f} A")
        print()
    print(f"control (crystallographic ZM241385, 1 replicate): "
          f"2.59 +/- 0.51 A, 100 % of frames below 4 A")
    print(f"\nwritten -> {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

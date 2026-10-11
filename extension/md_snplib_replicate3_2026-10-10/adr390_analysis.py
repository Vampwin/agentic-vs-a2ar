"""ADR-AIIM-390 analysis: replicate 3 of SNP00157 and SNP00163, against the three-run control.

Written before the runs end and not changed afterwards (PREREG_ADR390.md). It imports the ADR-AIIM-331 machinery
unchanged (series readers, index groups, anchors) and only adds replicate 3, the grading of R1 to R5 and the
three-against-three summary. Nothing is typed as a threshold except the ones of PREREG_ADR390.md (R2 to R5).

Order of use, after `prod_r3.gro` exists for both compounds:
  1. bash 20_new_runs/md_snplib_duplicate_2026-10-08/analysis/run_candidate_md11.sh <ID> prod_r3     (in WSL, per compound)
  2. python adr390_analysis.py --analyse
The anchors of ADR-AIIM-331 are re-run first; if one fails no grade is reported.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "md_snplib_duplicate_2026-10-08" / "analysis"))
import adr331_analysis as a331  # noqa: E402  (unchanged)

OUT = HERE
ROOT, CAND, COMPOUNDS = a331.ROOT, a331.CAND, a331.COMPOUNDS
# thresholds of PREREG_ADR390.md
R2_MAX_CALPHA_A = 10.0
R3_MAX_PERSIST_157 = 0.35
R4_RANGE_163 = (0.40, 0.90)
R5_MAX_FRAC = {"SNP00157": 0.50, "SNP00163": 0.60}


def main() -> int:
    anchors = a331.anchors()
    (OUT / "adr390_anchors.json").write_text(json.dumps(anchors, indent=1), "utf-8")
    if not anchors["all_anchors_met"]:
        print("ANCHOR FAILED -- no grade reported")
        return 1
    ctrl = {k: a331.run_metrics(d) for k, d in a331.CONTROL_SERIES.items()}
    per_run, inputs = {}, {}
    for c in COMPOUNDS:
        per_run[c] = {}
        for r in (1, 2, 3):
            deffnm = f"prod_r{r}"
            adir = CAND / c / "04_run" / f"analysis_{deffnm}"
            if not (adir / "ligand_rmsd.xvg").is_file():
                continue
            a331.n655_series(c, deffnm)
            per_run[c][f"r{r}"] = a331.run_metrics(adir)
            for f in ("ligand_rmsd.xvg", "box.xvg", "N6_55.xvg", "mindist.xvg"):
                inputs[str((adir / f).relative_to(ROOT))] = a331.sha(adir / f)

    grades = {}
    for c in COMPOUNDS:
        v = per_run[c].get("r3")
        if v is None:
            grades[c] = {"R1": "MISSED: no replicate 3 analysis"}
            continue
        g = {"R1_reached_50ns": "MET" if v["length_ns"] >= 49.95 else f"MISSED ({v['length_ns']} ns)",
             "R2_calpha_max_le_10A": "MET" if v["lig_calpha_mindist_max_A"] <= R2_MAX_CALPHA_A else "MISSED",
             "R5_frac_lt4A_max": R5_MAX_FRAC[c]}
        g["R5_rmsd_frac"] = "MET" if v["rmsd_frac_lt_4A"] <= R5_MAX_FRAC[c] else "MISSED"
        if c == "SNP00157":
            g["R3_persistence_le_0.35"] = "MET" if v["n655_persistence"] <= R3_MAX_PERSIST_157 else "MISSED"
        else:
            lo, hi = R4_RANGE_163
            g["R4_persistence_0.40_to_0.90"] = "MET" if lo <= v["n655_persistence"] <= hi else "MISSED"
        g["observed"] = {k: v[k] for k in ("length_ns", "n655_persistence", "rmsd_frac_lt_4A", "rmsd_mean_A",
                                           "lig_calpha_mindist_max_A", "n655_by_10ns", "area_change_pct")}
        grades[c] = g

    def rng(vals):
        return [round(min(vals), 3), round(max(vals), 3)]
    summary = {"control_r1-r3_persistence": rng([ctrl[k]["n655_persistence"] for k in ("r1", "r2", "r3")]),
               "control_r1-r3_rmsd_frac_lt4A": rng([ctrl[k]["rmsd_frac_lt_4A"] for k in ("r1", "r2", "r3")]),
               "control_r1-r3_rmsd_mean_A": rng([ctrl[k]["rmsd_mean_A"] for k in ("r1", "r2", "r3")])}
    for c in COMPOUNDS:
        rs = [per_run[c][k] for k in ("r1", "r2", "r3") if k in per_run[c]]
        summary[f"{c}_n_runs"] = len(rs)
        summary[f"{c}_persistence"] = rng([x["n655_persistence"] for x in rs])
        summary[f"{c}_rmsd_frac_lt4A"] = rng([x["rmsd_frac_lt_4A"] for x in rs])
        summary[f"{c}_rmsd_mean_A"] = rng([x["rmsd_mean_A"] for x in rs])
        # non-overlap of the per-run persistence ranges against the control (no test on n = 3)
        summary[f"{c}_persistence_below_every_control_run"] = bool(
            max(x["n655_persistence"] for x in rs) < min(ctrl[k]["n655_persistence"] for k in ("r1", "r2", "r3")))
    out = {"adr": "ADR-AIIM-390", "anchors_met": True, "control_runs": ctrl, "candidate_runs": per_run,
           "grading": grades, "three_against_three": summary, "inputs_sha256": inputs,
           "note": "per-run values and ranges only, no test on n = 3; a stable trajectory is not evidence the pose is "
                   "correct; no compound is called active"}
    (OUT / "adr390_results.json").write_text(json.dumps(out, indent=1), "utf-8")
    with open(OUT / "adr390_per_run.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["system", "run", "length_ns", "rmsd_mean_A", "rmsd_frac_lt_4A", "n655_persistence", "n655_by_10ns",
                    "area_change_pct", "lig_calpha_mindist_max_A"])
        for k, v in ctrl.items():
            w.writerow(["ZM241385 control", k, v["length_ns"], v["rmsd_mean_A"], v["rmsd_frac_lt_4A"], v["n655_persistence"],
                        v["n655_by_10ns"], v["area_change_pct"], v["lig_calpha_mindist_max_A"]])
        for c in COMPOUNDS:
            for k, v in per_run[c].items():
                w.writerow([c, k, v["length_ns"], v["rmsd_mean_A"], v["rmsd_frac_lt_4A"], v["n655_persistence"],
                            v["n655_by_10ns"], v["area_change_pct"], v["lig_calpha_mindist_max_A"]])
    print(json.dumps({"grading": grades, "three_against_three": summary}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

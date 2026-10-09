"""ADR-AIIM-331 analysis: SNP00157 and SNP00163, per-run, against the three-run control.

Pre-registration entry hash and log hash are in the ADR itself; the predictions are S1 to S4
(plus D3 of ADR-AIIM-316 for SNP00157). Nothing here changes a threshold: every number the
predictions compare against is read from the control artifact (ADR-AIIM-324/325), not typed.

Modes
  --anchors   no candidate trajectory is read. Checks that the machinery is right:
              A1  the N6.55 index groups this script builds from the control's own structure
                  equal the borrowed index file the control analysis used (so the same
                  series would result)
              A2  the series readers reproduce the stored control numbers (RMSD mean,
                  fraction at or below 4 A, N6.55 persistence, bilayer area change) for all
                  three control runs
              A3  for each candidate: the Asn is located by sequence position, the TM6
                  spacing check passes in the candidate's own structure, and the ligand
                  polar group equals the number of N and O heavy atoms of its SMILES
              A4  the nearest Anchor-to-Ligand-polar distance in each candidate's pre-
                  equilibration structure (em.gro) is within 0.5 A of the 2.87 / 2.67 A that
                  ADR-AIIM-331 read from the Vina pose
  --analyse   needs `run_candidate_md11.sh` to have written analysis_<deffnm>/ for each
              replicate. Builds the compound's own index groups, runs the N6.55 mindist, grades.
              Refuses to report a grade if an anchor failed.

A run that has not reached 50 ns is analysed over the length it has reached and labelled
partial; one under 25 ns is excluded and said to be (ADR-AIIM-331, Deadline).
The N6.55 contact is a heavy-atom N or O of the ligand within 3.5 A of OD1 or ND2 of the Asn
side chain, as in md_18_recognition_motif.py, whose functions are imported unchanged.
A stable trajectory is not evidence that the pose is correct, and no compound is called active.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Agentic VS Workflow")
SCRIPTS = ROOT / "10_new_system" / "scripts"
sys.path.insert(0, str(SCRIPTS))
import md_18_recognition_motif as m18  # noqa: E402  (frozen; imported, not modified)

RUNS = ROOT / "20_new_runs"
CAND = RUNS / "md_candidates_2026-08-22"
OUT = RUNS / "md_snplib_duplicate_2026-10-08"
CTRL_JSON = RUNS / "md_control_replicates_window50_2026-10-07" / "md_control_replicates_window50.json"
BORROWED_NDX = RUNS / "md_recognition_motif_2026-09-02" / "control_ZM241385" / "N6_55.ndx"
CONTROL_SERIES = {"r1": RUNS / "md_control_anchor_2026-10-07" / "analysis_prod",
                  "r2": RUNS / "md_control_replicates_2026-10-06" / "analysis_prod_r2",
                  "r3": RUNS / "md_control_replicates_2026-10-06" / "analysis_prod_r3"}
COMPOUNDS = ["SNP00157", "SNP00163"]
START_POLAR_A = {"SNP00157": 2.87, "SNP00163": 2.67}      # ADR-AIIM-331 table, from the Vina pose
WINDOW_NS, MIN_NS = 50.0, 25.0
CONTACT_A = m18.HBOND_CUTOFF_A                            # 3.5
RMSD_SITE_A, SITE_FRACTION = 4.0, 0.90      # S2; strictly below 4.0 A, the measure the quoted control figures use


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


# ------------------------------------------------------------------ series readers
def xvg(path: Path, scale=1.0, tmax=None, tscale=1.0):
    a = np.array([[float(x) for x in ln.split()[:3]] for ln in path.read_text().splitlines()
                  if ln.strip() and ln[0] not in "#@"])
    t = a[:, 0] * tscale
    keep = np.ones(len(t), bool) if tmax is None else t <= tmax + 1e-3
    return t[keep], a[keep, 1:] * scale


def run_metrics(adir: Path):
    """RMSD, box and N6.55 series of one run, over at most the first 50 ns."""
    t_r, r = xvg(adir / "ligand_rmsd.xvg", scale=10.0, tmax=WINDOW_NS)             # A
    t_n, n = xvg(adir / "N6_55.xvg", scale=10.0, tmax=WINDOW_NS)                   # A
    t_b, b = xvg(adir / "box.xvg", tmax=WINDOW_NS, tscale=1e-3)                    # box.xvg is in ps
    t_m, mm = xvg(adir / "mindist.xvg", scale=10.0, tmax=WINDOW_NS)                # ligand to receptor C-alpha, A
    rr, nn = r[:, 0], n[:, 0]
    area = b[:, 0] * b[:, 1]
    q = max(1, len(area) // 4)
    return {"length_ns": round(float(min(t_r[-1], t_n[-1])), 2),
            "n_frames_rmsd": int(len(rr)), "n_frames_n655": int(len(nn)),
            "rmsd_mean_A": round(float(rr.mean()), 2),
            "rmsd_frac_lt_4A": round(float((rr < RMSD_SITE_A).mean()), 4),
            "n655_persistence": round(float((nn <= CONTACT_A).mean()), 4),
            "lig_calpha_mindist_median_A": round(float(np.median(mm[:, 0])), 2),
            "lig_calpha_mindist_max_A": round(float(mm[:, 0].max()), 2),
            "n655_by_10ns": [round(float((nn[(t_n > a0) & (t_n <= a0 + 10)] <= CONTACT_A).mean()), 3)
                             for a0 in range(0, 50, 10) if ((t_n > a0) & (t_n <= a0 + 10)).any()],
            "area_change_pct": round(float((area[-q:].mean() / area[:q].mean() - 1) * 100), 2)}


# ------------------------------------------------------------------ index groups
def build_groups(gro: Path, anchor_pos: int, expect_n_polar=None):
    atoms = m18.read_gro(gro)
    inst = m18.residue_instances(atoms)
    r = inst[anchor_pos]
    assert r["resname"] == "ASN", f"{gro}: position {anchor_pos} is {r['resname']}"
    for off, want in m18.TM6_SPACING.items():
        got = inst[anchor_pos + off]["resname"]
        want_t = (want,) if isinstance(want, str) else want
        assert got in want_t, f"TM6 spacing failed at {off:+d}: {got}"
    a_idx = [a[1] for a in r["atoms"] if a[0] in m18.ANCHOR_ATOMS["N6.55"][1]]
    _, lig_polar = m18.ligand_atoms(atoms, polar_only=True)
    l_idx = [a[1] for a in lig_polar]
    if expect_n_polar is not None:
        assert len(l_idx) == expect_n_polar, f"polar atoms {len(l_idx)} != {expect_n_polar}"
    return atoms, a_idx, l_idx


def control_anchor_position():
    atoms = m18.read_gro(m18.CONTROL / "em.gro")
    inst = m18.residue_instances(atoms)
    _, lig_all = m18.ligand_atoms(atoms)
    return m18.locate_anchors(inst, lig_all)["N6.55"]["position"]


def read_ndx(p: Path):
    g, cur = {}, None
    for ln in p.read_text().splitlines():
        if ln.startswith("["):
            cur = ln.strip("[] ").strip(); g[cur] = []
        elif cur and ln.strip():
            g[cur] += [int(x) for x in ln.split()]
    return g


def min_image_distance(atoms, a_idx, l_idx, gro: Path):
    box = [float(x) * 10 for x in gro.read_text().splitlines()[-1].split()[:3]]
    xyz = {a[3]: np.array(a[4:7]) for a in atoms}
    best = 1e9
    for i in a_idx:
        for j in l_idx:
            d = xyz[i] - xyz[j]
            d -= np.round(d / box) * box
            best = min(best, float(np.linalg.norm(d)))
    return best


# ------------------------------------------------------------------ anchors
def anchors():
    res, ok = {}, True
    pos = control_anchor_position()
    # A1
    atoms, a_idx, l_idx = build_groups(m18.CONTROL / "em.gro", pos, expect_n_polar=9)
    bor = read_ndx(BORROWED_NDX)
    a1 = (sorted(a_idx) == sorted(bor["Anchor"]) and sorted(l_idx) == sorted(bor["Ligand"]))
    res["A1_groups_equal_borrowed_control_ndx"] = {"anchor_atoms": len(a_idx), "ligand_polar_atoms": len(l_idx),
                                                   "verdict": "MET" if a1 else "MISSED"}
    ok &= a1
    # A2
    stored = json.loads(CTRL_JSON.read_text("utf-8"))
    ref = {"r1": stored["existing_control"]["50"], **stored["replicates"]}
    a2, rows = True, {}
    for k, d in CONTROL_SERIES.items():
        m = run_metrics(d)
        s = ref[k]
        good = (abs(m["rmsd_mean_A"] - s["rmsd"]["mean_rmsd_A"]) < 0.006
                and abs(m["rmsd_frac_lt_4A"] - s["rmsd"]["fraction_below"]["4.0"]) < 0.0006
                and abs(m["n655_persistence"] - s["n655"]["persistence"]) < 0.0006
                and abs(m["area_change_pct"] - s["box"]["percent_change"]) < 0.006)
        rows[k] = {"recomputed": {x: m[x] for x in ("rmsd_mean_A", "rmsd_frac_lt_4A", "n655_persistence", "area_change_pct")},
                   "stored": {"rmsd_mean_A": s["rmsd"]["mean_rmsd_A"], "rmsd_frac_lt_4A": s["rmsd"]["fraction_below"]["4.0"],
                              "n655_persistence": s["n655"]["persistence"], "area_change_pct": s["box"]["percent_change"]},
                   "verdict": "MET" if good else "MISSED"}
        a2 &= good
    res["A2_readers_reproduce_stored_control_numbers"] = {"runs": rows, "verdict": "MET" if a2 else "MISSED"}
    ok &= a2
    # A3, A4
    cand = {}
    import pandas as pd
    led = pd.read_csv(RUNS / "snp_screen_2026-08-04/G1_standardization/snp_source_ledger.csv").set_index("SNP_ID")
    from rdkit import Chem
    for c in COMPOUNDS:
        mol = Chem.MolFromSmiles(led.loc[c, "smiles_canonical"])
        n_no = sum(1 for a in mol.GetAtoms() if a.GetSymbol() in ("N", "O"))
        gro = CAND / c / "04_run" / "em.gro"
        try:
            atoms, a_idx, l_idx = build_groups(gro, pos, expect_n_polar=n_no)
            d = min_image_distance(atoms, a_idx, l_idx, gro)
            a3 = True
            a4 = abs(d - START_POLAR_A[c]) <= 0.5
            cand[c] = {"A3_anchor_and_polar_group": {"anchor_atoms": len(a_idx), "ligand_polar_atoms": len(l_idx),
                                                     "smiles_N_plus_O": n_no, "tm6_spacing": "passed", "verdict": "MET"},
                       "A4_start_polar_distance": {"em_gro_A": round(d, 2), "vina_pose_A": START_POLAR_A[c],
                                                   "tolerance_A": 0.5, "verdict": "MET" if a4 else "MISSED"}}
            ok &= a3 and a4
        except AssertionError as e:
            cand[c] = {"A3_anchor_and_polar_group": {"verdict": "MISSED", "reason": str(e)[:120]}}
            ok = False
    res["candidates"] = cand
    res["all_anchors_met"] = bool(ok)
    return res


# ------------------------------------------------------------------ analyse
def n655_series(c: str, deffnm: str):
    run = CAND / c / "04_run"
    adir = run / f"analysis_{deffnm}"
    gro = run / "em.gro"
    atoms, a_idx, l_idx = build_groups(gro, control_anchor_position())
    m18.write_index(adir / "N6_55.ndx", {"Anchor": a_idx, "Ligand": l_idx})
    ok = m18.run_mindist(run / f"{deffnm}.tpr", run / f"{deffnm}.xtc", adir / "N6_55.ndx", "Anchor",
                         adir / "N6_55.xvg", adir / "N6_55.log")
    if not ok:
        raise SystemExit(f"{c} {deffnm}: gmx mindist produced nothing")


def grade(per_run, ctrl):
    cmin = min(r["n655_persistence"] for r in ctrl.values())
    cmax_rmsd = max(r["rmsd_mean_A"] for r in ctrl.values())
    cmean = float(np.mean([r["n655_persistence"] for r in ctrl.values()]))
    out = {"thresholds_read_from_control": {"min_persistence": cmin, "max_rmsd_mean_A": cmax_rmsd,
                                            "mean_persistence": round(cmean, 4),
                                            "min_frac_le_4A": min(r["rmsd_frac_lt_4A"] for r in ctrl.values())}}
    for c in COMPOUNDS:
        runs = {k: v for k, v in per_run.get(c, {}).items() if v["length_ns"] >= MIN_NS}
        excl = [k for k, v in per_run.get(c, {}).items() if v["length_ns"] < MIN_NS]
        if not runs:
            out[c] = {"status": "NOT EVALUABLE: no replicate reached 25 ns"}
            continue
        partial = [k for k, v in runs.items() if v["length_ns"] < WINDOW_NS - 0.05]
        s2 = all(v["rmsd_frac_lt_4A"] >= SITE_FRACTION for v in runs.values())
        s3 = all(v["n655_persistence"] < cmin for v in runs.values())
        s4 = all(v["rmsd_mean_A"] > cmax_rmsd for v in runs.values())
        d3 = float(np.mean([v["n655_persistence"] for v in runs.values()])) < cmean
        out[c] = {"replicates_analysed": sorted(runs), "partial": partial, "excluded_under_25ns": excl,
                  "S2_in_site_ge_90pct_each": "MET" if s2 else "MISSED",
                  "S3_persistence_each_below_control_min": "MET" if s3 else "MISSED",
                  "S4_rmsd_each_above_control_max": "MET" if s4 else "MISSED",
                  "mean_persistence": round(float(np.mean([v["n655_persistence"] for v in runs.values()])), 4)}
        if c == "SNP00157":
            out[c]["D3_mean_persistence_below_control_mean"] = "MET" if d3 else "MISSED"
    return out


def analyse():
    a = anchors()
    print(json.dumps({k: v for k, v in a.items() if k != "candidates"}, indent=1)[:1500])
    if not a["all_anchors_met"]:
        print("ANCHOR FAILED -- no grade reported")
        (OUT / "adr331_anchors.json").write_text(json.dumps(a, indent=1), "utf-8")
        return 1
    stored = json.loads(CTRL_JSON.read_text("utf-8"))
    ctrl = {k: run_metrics(d) for k, d in CONTROL_SERIES.items()}
    per_run, inputs = {}, {}
    for c in COMPOUNDS:
        per_run[c] = {}
        for r in (1, 2):                                    # r3 is not run (ADR-AIIM-331 design 2)
            deffnm = f"prod_r{r}"
            adir = CAND / c / "04_run" / f"analysis_{deffnm}"
            if not (adir / "ligand_rmsd.xvg").is_file():
                continue
            n655_series(c, deffnm)
            per_run[c][f"r{r}"] = run_metrics(adir)
            for f in ("ligand_rmsd.xvg", "box.xvg", "N6_55.xvg"):
                inputs[str((adir / f).relative_to(ROOT))] = sha(adir / f)
    attn = {c: json.loads((CAND / c / "build_metadata.json").read_text())["parmchk2_attn_count"] for c in COMPOUNDS}
    g = grade(per_run, ctrl)
    for c in COMPOUNDS:
        g.setdefault(c, {})["S1_zero_ATTN"] = "MET" if attn[c] == 0 else "MISSED"
    pair = {"control_r1+r2": round(float(np.mean([ctrl[k]["n655_persistence"] for k in ("r1", "r2")])), 3)}
    for c in COMPOUNDS:
        rr = [per_run[c][k]["n655_persistence"] for k in ("r1", "r2") if k in per_run[c]]
        pair[f"{c}_r1+r2"] = round(float(np.mean(rr)), 3) if len(rr) == 2 else None
    out = {"adr": "ADR-AIIM-331", "anchors": a, "control_runs": ctrl, "candidate_runs": per_run,
           "matched_duplicate_persistence": pair,
           "all_three_control_mean_persistence": g["thresholds_read_from_control"]["mean_persistence"],
           "grading": g, "parmchk2_attn": attn, "inputs_sha256": inputs,
           "note": "stable trajectory is not evidence the pose is correct; no compound is called active"}
    (OUT / "adr331_results.json").write_text(json.dumps(out, indent=1), "utf-8")
    with open(OUT / "adr331_per_run.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["system", "run", "length_ns", "rmsd_mean_A", "rmsd_frac_lt_4A", "n655_persistence",
                    "n655_by_10ns", "area_change_pct", "lig_calpha_mindist_max_A (post hoc)"])
        for k, v in ctrl.items():
            w.writerow(["ZM241385 control", k, v["length_ns"], v["rmsd_mean_A"], v["rmsd_frac_lt_4A"],
                        v["n655_persistence"], v["n655_by_10ns"], v["area_change_pct"], v["lig_calpha_mindist_max_A"]])
        for c in COMPOUNDS:
            for k, v in per_run[c].items():
                w.writerow([c, k, v["length_ns"], v["rmsd_mean_A"], v["rmsd_frac_lt_4A"],
                            v["n655_persistence"], v["n655_by_10ns"], v["area_change_pct"], v["lig_calpha_mindist_max_A"]])
    print(json.dumps(g, indent=1))
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--anchors", action="store_true")
    ap.add_argument("--analyse", action="store_true")
    args = ap.parse_args()
    if args.anchors:
        a = anchors()
        (OUT / "adr331_anchors.json").write_text(json.dumps(a, indent=1), "utf-8")
        print(json.dumps(a, indent=1))
        return 0 if a["all_anchors_met"] else 1
    if args.analyse:
        return analyse()
    ap.print_help()
    return 2


if __name__ == "__main__":
    raise SystemExit(main())

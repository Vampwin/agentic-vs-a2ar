"""POST HOC (after ADR-AIIM-392; not pre-registered): the same analysis as posthoc_pose_clusters.py, for SNP00157, SNP00163 and the ZM241385 control.
posthoc_pose_clusters.py is kept unchanged (it produced ADR-AIIM-392); this file only generalises where the three runs of a system are read from.

Original text of the question: which poses does the ligand visit in each replicate, and which one holds N6.55?

Why. Replicate 3 of SNP00157 formed the N6.55 contact late and held it (persistence 0.773) where replicates 1 and 2 did not (0.148, 0.087).
This asks whether that is a different pose or the same pose reached at a different time.

Method (fixed here, before it was run; the cutoff was not tuned on the outcome).
  - frames: every 5th of the 1,001 frames of `whole.xtc` (0.25 ns), the trajectory the frozen md_11 wrote (molecules made whole, centred);
  - each frame is superposed on the receptor C-alpha atoms of the pre-equilibration structure (em.gro); the ligand is NOT fitted, so the
    RMSD between two frames is the movement of the ligand heavy atoms in the receptor frame;
  - the frames of the 3 replicates of one compound are pooled and clustered, average linkage on the pairwise heavy-atom RMSD,
    cut at 3.0 A (2.0 and 4.0 A are reported as sensitivity);
  - the N6.55 contact state of every frame is read from the N6_55.xvg that adr390_analysis.py wrote (<= 3.5 A), at the same frames.
Everything reported is descriptive. No test, no ranking; a cluster is a pose family seen in the trajectories, not evidence that it is correct.

    python posthoc_pose_clusters_all.py SNP00157 SNP00163 CONTROL
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import align
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

ROOT = Path(r"D:\Agentic VS Workflow")
CAND = ROOT / "20_new_runs" / "md_candidates_2026-08-22"
OUT = Path(__file__).resolve().parent
RUNS = ROOT / "20_new_runs"
CTRL_TOP = RUNS / "md_production_2026-08-21" / "04_run" / "em_prod.gro"
CTRL_DIRS = {1: RUNS / "md_control_anchor_2026-10-07" / "analysis_prod",
             2: RUNS / "md_control_replicates_2026-10-06" / "analysis_prod_r2",
             3: RUNS / "md_control_replicates_2026-10-06" / "analysis_prod_r3"}


def where(c: str, r: int):
    """(topology, analysis directory) of run r of a system."""
    if c == "CONTROL":
        return CTRL_TOP, CTRL_DIRS[r]
    run = CAND / c / "04_run"
    return run / "em.gro", run / f"analysis_prod_r{r}"
STEP = 5
CUTS = (3.0, 2.0, 4.0)
CONTACT_A = 3.5


def read_ndx(p: Path):
    g, cur = {}, None
    for ln in p.read_text().splitlines():
        if ln.startswith("["):
            cur = ln.strip("[] ").strip(); g[cur] = []
        elif cur and ln.strip():
            g[cur] += [int(x) for x in ln.split()]
    return g


def xvg(p: Path):
    a = np.array([[float(x) for x in ln.split()[:2]] for ln in p.read_text().splitlines() if ln.strip() and ln[0] not in "#@"])
    return a[:, 0], a[:, 1] * 10.0


def frames_of(c: str, r: int, ref: mda.Universe, ca_idx, lig_idx):
    top, adir = where(c, r)
    u = mda.Universe(str(top), str(adir / "whole.xtc"))
    assert u.atoms.n_atoms == ref.atoms.n_atoms
    ca, lig = u.atoms[ca_idx], u.atoms[lig_idx]
    ref_ca = ref.atoms[ca_idx].positions - ref.atoms[ca_idx].positions.mean(0)
    xyz, t = [], []
    for ts in u.trajectory[::STEP]:
        mob = ca.positions - ca.positions.mean(0)
        rot, _ = align.rotation_matrix(mob, ref_ca)
        xyz.append((lig.positions - ca.positions.mean(0)) @ rot.T)
        t.append(ts.time / 1000.0)
    tn, d = xvg(adir / "N6_55.xvg")
    contact = np.array([d[np.argmin(abs(tn - x))] <= CONTACT_A for x in t])
    return np.array(xyz), np.array(t), contact


def pairwise_rmsd(x):
    n = len(x)
    f = x.reshape(n, -1)
    sq = (f ** 2).sum(1)
    d2 = (sq[:, None] + sq[None, :] - 2 * f @ f.T) / x.shape[1]
    return np.sqrt(np.clip(d2, 0, None))


def analyse(c: str):
    top, adir1 = where(c, 1)
    ndx = read_ndx(adir1 / "analysis.ndx")
    ca_idx = np.array(ndx["Calpha"]) - 1
    ref = mda.Universe(str(top))
    lig_all = np.array(ndx["Ligand"]) - 1
    names = ref.atoms.names[lig_all]
    lig_idx = lig_all[[not n.upper().startswith("H") for n in names]]
    ref_ca = ref.atoms[ca_idx].positions
    ref_lig = ref.atoms[lig_idx].positions - ref_ca.mean(0)       # the pose as placed, in the receptor frame
    X, T, C, R = [], [], [], []
    for r in (1, 2, 3):
        x, t, k = frames_of(c, r, ref, ca_idx, lig_idx)
        X.append(x); T.append(t); C.append(k); R.append(np.full(len(t), r))
    X, T, C, R = map(np.concatenate, (X, T, C, R))
    allx = np.concatenate([X, ref_lig[None]])                       # last row = the starting pose
    D = pairwise_rmsd(allx)
    to_start = D[-1, :-1]
    D = D[:-1, :-1]
    Z = linkage(squareform(D, checks=False), method="average")
    res = {"compound": c, "n_heavy_atoms": int(len(lig_idx)), "frames_per_run": int((R == 1).sum()), "step_ns": STEP * 0.05}
    for cut in CUTS:
        lab = fcluster(Z, t=cut, criterion="distance")
        ids, cnt = np.unique(lab, return_counts=True)
        order = ids[np.argsort(-cnt)]
        rows = []
        for k in order:
            m = lab == k
            if m.sum() < 5:
                continue
            sub = D[np.ix_(m, m)]
            med = np.where(m)[0][np.argmin(sub.sum(0))]
            rows.append({"cluster": int(k), "frames": int(m.sum()),
                         "by_run": {f"r{r}": int((m & (R == r)).sum()) for r in (1, 2, 3)},
                         "contact_fraction": round(float(C[m].mean()), 3),
                         "median_rmsd_to_start_pose_A": round(float(np.median(to_start[m])), 2),
                         "medoid_run_time_ns": [int(R[med]), round(float(T[med]), 2)],
                         "mean_within_rmsd_A": round(float(sub[np.triu_indices(m.sum(), 1)].mean()), 2) if m.sum() > 1 else 0.0})
        res[f"cut_{cut}"] = {"n_clusters_with_>=5_frames": len(rows), "clusters": rows}
        if cut == CUTS[0]:
            main_lab = lab
            main_rows = rows
    # time course of the main clusters per replicate, 5 ns bins
    top = [r["cluster"] for r in main_rows[:5]]
    tc = {}
    for r in (1, 2, 3):
        m = R == r
        bins = []
        for b0 in range(0, 50, 5):
            sel = m & (T > b0) & (T <= b0 + 5)
            if not sel.any():
                continue
            occ = {f"c{k}": round(float((main_lab[sel] == k).mean()), 2) for k in top}
            bins.append({"ns": f"{b0}-{b0 + 5}", **occ, "contact": round(float(C[sel].mean()), 2)})
        tc[f"r{r}"] = bins
    res["time_course_cut_3.0"] = tc
    # the contact pose of replicate 3 against everything else
    for r in (1, 2, 3):
        m = (R == r) & C
        res[f"contact_frames_r{r}"] = {"n": int(m.sum()),
                                      "clusters": {int(k): int((main_lab[m] == k).sum()) for k in np.unique(main_lab[m])}}
    # pose before contact (r3, frames up to the first 3 ns of the contact being held) versus pose while in contact
    m3 = R == 3
    t3, c3 = T[m3], C[m3]
    run_len = 0; onset = None
    for t_, k in zip(t3, c3):
        run_len = run_len + 1 if k else 0
        if run_len >= 8:                       # 2 ns of contact in a row
            onset = float(t_ - 7 * STEP * 0.05); break
    res["r3_contact_onset_ns_(first 2 ns of contact in a row)"] = onset
    if onset is not None:
        pre = m3 & (T < onset - 0.01)
        post = m3 & (T >= onset) & C
        res["r3_pre_vs_post_contact"] = {
            "pre_frames": int(pre.sum()), "post_contact_frames": int(post.sum()),
            "median_rmsd_pre_to_start_A": round(float(np.median(to_start[pre])), 2) if pre.any() else None,
            "median_rmsd_post_to_start_A": round(float(np.median(to_start[post])), 2),
            "median_rmsd_pre_to_post_A": round(float(np.median(D[np.ix_(np.where(pre)[0], np.where(post)[0])])), 2) if pre.any() else None}
    # r1 and r2 frames against the r3 contact pose: how close does each ever get?
    post3 = np.where((R == 3) & C & (T >= (onset or 0)))[0]
    res["distance_to_r3_contact_pose_A"] = {
        f"r{r}": {"min": round(float(D[np.ix_(np.where(R == r)[0], post3)].min()), 2),
                  "median_of_per_frame_min": round(float(np.median(D[np.ix_(np.where(R == r)[0], post3)].min(1))), 2)}
        for r in (1, 2, 3)}
    return res


if __name__ == "__main__":
    out = {c: analyse(c) for c in sys.argv[1:]}
    (OUT / "posthoc_pose_clusters_all.json").write_text(json.dumps(out, indent=1), "utf-8")
    for c, r in out.items():
        print("=====", c, "heavy atoms", r["n_heavy_atoms"], "frames/run", r["frames_per_run"])
        for cut in CUTS:
            print(f"-- cut {cut} A: {r[f'cut_{cut}']['n_clusters_with_>=5_frames']} clusters (>=5 frames)")
            for k in r[f"cut_{cut}"]["clusters"][:6]:
                print("  ", k)
        print("-- r3 contact onset", r["r3_contact_onset_ns_(first 2 ns of contact in a row)"], r.get("r3_pre_vs_post_contact"))
        print("-- contact frames per run, by cluster (cut 3.0)", {k: r[k] for k in ("contact_frames_r1", "contact_frames_r2", "contact_frames_r3")})
        print("-- distance of each run to the r3 contact pose", r["distance_to_r3_contact_pose_A"])
        for run, bins in r["time_course_cut_3.0"].items():
            print("-- time course", run)
            for b in bins:
                print("   ", b)

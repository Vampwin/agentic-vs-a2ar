"""POST HOC (ADR-AIIM-394; not pre-registered): receptor C-alpha RMSF and ligand-protein hydrogen bonds, per run, for SNP00157, SNP00163 and the ZM241385 control.

Definitions, fixed here before the numbers were read.
  Window       the first 50 ns of every run (frames with t <= 50 ns), the window of every other MD number of the study. The control
               runs 1 and 2 are longer (100 ns, 53.1 ns) and are cut at 50 ns; the control frames are 100 ps apart, the compounds' 50 ps.
  Chains       the receptor is two chains (212 + 91 residues, each capped and each numbered from 1), so a residue is named by chain and number
               ("A:ASN36"), never by number alone. The numbers are those of the simulation structure, not of the receptor sequence; only the N6.55
               asparagine, located by sequence position in the earlier analyses, is named as such.
  Image jumps  in one run (SNP00163 replicate 2) 20 of the 1,001 frames of whole.xtc place the second chain in another periodic image (C-alpha
               superposition RMSD about 36 A against 1 to 2 A in every other frame; ligand RMSD 58 A while gmx mindist keeps the ligand 4.7 A
               from the receptor). Frames whose C-alpha superposition RMSD exceeds 5 A are excluded from the RMSF and from the hydrogen-bond
               occupancy and are counted (check_image_jumps.py).
  RMSF         C-alpha atoms of the receptor (the 'Calpha' group of the md_11 index), each frame superposed on the receptor C-alpha atoms
               of the pre-equilibration structure, RMSF about the run's own mean position, in angstrom, per residue and per run.
               "Site" residues = receptor residues with a heavy atom within 4.5 A of a ligand heavy atom in the pre-equilibration structure.
  Hydrogen     gmx hbond (GROMACS 2026.3, the 2024 implementation), geometric definition at its recommended defaults: donor-acceptor
  bonds        distance <= 0.35 nm and hydrogen-donor-acceptor angle <= 30 degrees; ligand against protein only (no water, no lipid),
               frame by frame (run_hbond.sh). Occupancy = share of frames in the window with the bond, a donor-acceptor pair counted once
               per frame whichever hydrogen is used.
Everything is descriptive. No test; no free energy; a hydrogen bond present in a frame is not evidence that a pose is correct.

    python posthoc_rmsf_hbond.py
"""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path

import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import align

ROOT = Path(r"D:\Agentic VS Workflow")
RUNS = ROOT / "20_new_runs"
CAND = RUNS / "md_candidates_2026-08-22"
OUT = Path(__file__).resolve().parent
HB = OUT / "hbond"
WINDOW_PS = 50000.0
FIT_BAD_A = 5.0
SYSTEMS = {
    "CONTROL": {"label": "ZM241385 control", "resname": "ZMA", "top": RUNS / "md_production_2026-08-21/04_run/em_prod.gro",
                "dirs": {1: RUNS / "md_control_anchor_2026-10-07/analysis_prod",
                         2: RUNS / "md_control_replicates_2026-10-06/analysis_prod_r2",
                         3: RUNS / "md_control_replicates_2026-10-06/analysis_prod_r3"}},
}
for c in ("SNP00157", "SNP00163"):
    SYSTEMS[c] = {"label": c, "resname": "LIG", "top": CAND / c / "04_run/em.gro",
                  "dirs": {r: CAND / c / "04_run" / f"analysis_prod_r{r}" for r in (1, 2, 3)}}
N655_ANCHOR_ATOM = 3814        # atom of the N6.55 asparagine side chain in the N6_55.ndx of the earlier analyses (1-based)


def read_ndx(p: Path):
    g, cur = {}, None
    for ln in p.read_text().splitlines():
        if ln.startswith("["):
            cur = ln.strip("[] ").strip(); g[cur] = []
        elif cur and ln.strip():
            g[cur] += [int(x) for x in ln.split()]
    return g


def chain_labels(u):
    """Chain index of every residue: a new protein chain begins where the residue number falls."""
    rid = u.residues.resids
    prot = np.isin(u.residues.resindices, u.select_atoms("protein").residues.resindices)
    ch = np.zeros(len(rid), int)
    k = 0
    for i in range(1, len(rid)):
        if prot[i] and prot[i - 1] and rid[i] < rid[i - 1]:
            k += 1
        ch[i] = k
    return ch


def rlabel(u, ch, resindex):
    return f"{'ABCD'[ch[resindex]]}:{u.residues.resnames[resindex]}{u.residues.resids[resindex]}"


def rmsf_run(system: str, r: int):
    s = SYSTEMS[system]
    adir = s["dirs"][r]
    ref = mda.Universe(str(s["top"]))
    ca_idx = np.array(read_ndx(adir / "analysis.ndx")["Calpha"]) - 1
    u = mda.Universe(str(s["top"]), str(adir / "whole.xtc"))
    assert u.atoms.n_atoms == ref.atoms.n_atoms
    ca = u.atoms[ca_idx]
    ref_ca = ref.atoms[ca_idx].positions
    ref_c = ref_ca.mean(0)
    X, fitr = [], []
    for ts in u.trajectory:
        if ts.time > WINDOW_PS + 1e-3:
            break
        mob = ca.positions - ca.positions.mean(0)
        rot, v = align.rotation_matrix(mob, ref_ca - ref_c)
        X.append(mob @ rot.T)
        fitr.append(v)
    X, fitr = np.array(X), np.array(fitr)
    bad = fitr > FIT_BAD_A
    Xg = X[~bad]
    rmsf = np.sqrt(((Xg - Xg.mean(0)) ** 2).sum(2).mean(0))
    return rmsf, int((~bad).sum()), int(bad.sum()), np.where(bad)[0], ref.atoms[ca_idx].resindices


def site_residues(system: str):
    s = SYSTEMS[system]
    u = mda.Universe(str(s["top"]))
    lig = u.select_atoms(f"resname {s['resname']} and not name H*")
    prot = u.select_atoms("protein and not name H* and around 4.5 group lig", lig=lig)
    return sorted(set(int(x) for x in prot.resindices))


def hbond_run(system: str, r: int, top: mda.Universe, ch, exclude=()):
    lab = f"{system}_r{r}"
    t = np.array([float(ln.split()[0]) for ln in (HB / f"{lab}.num.xvg").read_text().splitlines() if ln.strip() and ln[0] not in "#@"])
    frames, cur = {}, None
    for ln in (HB / f"{lab}.ndx").read_text().splitlines():
        m = re.match(r"\[ hbonds_.*_frame_(\d+) \]", ln)
        if m:
            cur = int(m.group(1)); frames[cur] = []
        elif cur is not None and ln.strip() and not ln.startswith("["):
            frames[cur] += [int(x) for x in ln.split()]
    lig_res = SYSTEMS[system]["resname"]
    names, resname, ridx = top.atoms.names, top.atoms.resnames, top.atoms.resindices
    ex = set(int(i) for i in exclude)
    keep = [i for i, x in enumerate(t) if x <= WINDOW_PS + 1e-3 and i not in ex]
    pair_frames, per_frame = {}, []
    for i in keep:
        trip = np.array(frames.get(i, []), dtype=int).reshape(-1, 3)
        seen = set()
        for d, h, a in trip:
            d0, a0 = d - 1, a - 1
            lig_is_donor = resname[d0] == lig_res
            la, pa = (d0, a0) if lig_is_donor else (a0, d0)
            seen.add(("ligand donor" if lig_is_donor else "ligand acceptor", str(names[la]), rlabel(top, ch, int(ridx[pa])), str(names[pa])))
        per_frame.append(len(seen))
        for k in seen:
            pair_frames[k] = pair_frames.get(k, 0) + 1
    n = len(keep)
    pf = np.array(per_frame)
    pairs = sorted(({"role": k[0], "ligand_atom": k[1], "residue": k[2], "protein_atom": k[3], "occupancy": round(v / n, 3)}
                    for k, v in pair_frames.items()), key=lambda x: -x["occupancy"])
    return {"frames": n, "time_range_ns": [round(float(t[keep[0]]) / 1000, 2), round(float(t[keep[-1]]) / 1000, 2)],
            "mean_bonds_per_frame": round(float(pf.mean()), 3), "frames_with_ge1": round(float((pf >= 1).mean()), 3),
            "pairs_ge_0.05": [p for p in pairs if p["occupancy"] >= 0.05]}


def main():
    res, rows, ref_labels = {}, [], None
    for system in SYSTEMS:
        top = mda.Universe(str(SYSTEMS[system]["top"]))
        ch = chain_labels(top)
        site = site_residues(system)
        n655 = int(top.atoms.resindices[N655_ANCHOR_ATOM - 1])
        out = {"site_residues": [rlabel(top, ch, i) for i in site], "n655_residue_in_structure": rlabel(top, ch, n655), "rmsf": {}, "hbond": {}}
        for r in (1, 2, 3):
            f, n, nbad, bad_idx, rix = rmsf_run(system, r)
            labels = [rlabel(top, ch, int(i)) for i in rix]
            if ref_labels is None:
                ref_labels = labels
            assert labels == ref_labels, f"{system}: receptor differs from the first system"
            ins = np.isin(rix, site)
            out["rmsf"][f"r{r}"] = {"frames_used": n, "frames_excluded_image_jump": nbad, "all_mean_A": round(float(f.mean()), 3),
                                    "site_mean_A": round(float(f[ins].mean()), 3), "non_site_mean_A": round(float(f[~ins].mean()), 3),
                                    "site_max_A": round(float(f[ins].max()), 3),
                                    "most_mobile_5": [[labels[i], round(float(f[i]), 2)] for i in np.argsort(-f)[:5]],
                                    "profile_A": [round(float(x), 3) for x in f]}
            for i in range(len(f)):
                rows.append([system, r, labels[i], round(float(f[i]), 3), int(ins[i])])
            out["hbond"][f"r{r}"] = hbond_run(system, r, top, ch, exclude=bad_idx)
            out["hbond"][f"r{r}"]["frames_excluded_image_jump"] = nbad
        res[system] = out
        print(f"{system}: site residues {len(site)}, N6.55 = {out['n655_residue_in_structure']}")
    res["_residues"] = ref_labels
    (OUT / "posthoc_rmsf_hbond.json").write_text(json.dumps(res, indent=1), "utf-8")
    with open(OUT / "posthoc_rmsf_per_residue.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["system", "run", "residue_chain:name+number_in_structure", "calpha_rmsf_A", "site_residue"])
        w.writerows(rows)
    for system in SYSTEMS:
        o = res[system]
        print("=====", system, "| N6.55 asparagine =", o["n655_residue_in_structure"], "| site residues:", o["site_residues"])
        for r in (1, 2, 3):
            f, h = o["rmsf"][f"r{r}"], o["hbond"][f"r{r}"]
            print(f"  r{r}: RMSF all {f['all_mean_A']} site {f['site_mean_A']} non-site {f['non_site_mean_A']} (frames {f['frames_used']}, excluded {f['frames_excluded_image_jump']}) | "
                  f"H-bonds/frame {h['mean_bonds_per_frame']}, frames with >=1: {h['frames_with_ge1']}")
            for p in h["pairs_ge_0.05"][:6]:
                print("       ", p)
            print("       most mobile:", f["most_mobile_5"])


if __name__ == "__main__":
    main()

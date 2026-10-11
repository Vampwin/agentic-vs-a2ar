"""POST HOC (ADR-AIIM-395): how far is the ICL3 gap (where BRIL was removed and the two receptor segments were capped) from the ligand?

Gap residues: the last real residue of segment A and the first real residue of segment B (caps excluded), plus the caps themselves.
Distances are minimum-image (box-aware), so they do not depend on how whole.xtc was written.
  - starting structure (em.gro / em_prod.gro): ligand heavy atoms to (a) the gap C-alpha atoms, (b) any heavy atom of the two gap residues and their caps,
    and the vertical (z, membrane normal) separation of the ligand centroid from the gap;
  - trajectory: the same minimum heavy-atom distance and the centroid distance to the gap C-alpha atoms over the first 50 ns, every 5th frame.
For orientation the same distances are given for the N6.55 asparagine side chain, and the length of the gap in the starting structure.
"""
import json
from pathlib import Path
import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import distance_array
import posthoc_rmsf_hbond as P

out = {}
for system, s in P.SYSTEMS.items():
    top = mda.Universe(str(s["top"]))
    ch = P.chain_labels(top)
    prot = top.select_atoms("protein")
    res = prot.residues
    ra = [r for r in res if ch[r.resindex] == 0 and r.resname not in ("ACE", "NME")]
    rb = [r for r in res if ch[r.resindex] == 1 and r.resname not in ("ACE", "NME")]
    lastA, firstB = ra[-1], rb[0]
    gap_res = [r for r in res if r.resindex in (lastA.resindex - 0, lastA.resindex + 1, firstB.resindex - 1, firstB.resindex)]
    gap_heavy = sum([list(r.atoms.select_atoms("not name H*")) for r in gap_res], [])
    gap_heavy = mda.AtomGroup(gap_heavy)
    gap_ca = mda.AtomGroup([lastA.atoms.select_atoms("name CA")[0], firstB.atoms.select_atoms("name CA")[0]])
    lig = top.select_atoms(f"resname {s['resname']} and not name H*")
    n655 = top.atoms[P.N655_ANCHOR_ATOM - 1].residue.atoms.select_atoms("name OD1 ND2")
    box = top.dimensions
    def dists(u_lig, u_gap_h, u_gap_ca, u_n655, box):
        return {"min_heavy_to_gap_residues_A": float(distance_array(u_lig.positions, u_gap_h.positions, box=box).min()),
                "centroid_to_gap_CA_A": [float(distance_array(u_lig.positions.mean(0)[None], u_gap_ca.positions[i][None], box=box)[0, 0]) for i in range(2)],
                "min_heavy_to_N655_side_chain_A": float(distance_array(u_lig.positions, u_n655.positions, box=box).min()),
                "z_centroid_minus_gap_CA_A": [float(u_lig.positions.mean(0)[2] - u_gap_ca.positions[i][2]) for i in range(2)]}
    start = dists(lig, gap_heavy, gap_ca, n655, box)
    start["gap_CA_to_CA_A"] = float(np.linalg.norm(gap_ca.positions[0] - gap_ca.positions[1]))
    start["gap_residues"] = [P.rlabel(top, ch, r.resindex) for r in gap_res]
    traj = {}
    for r in (1, 2, 3):
        u = mda.Universe(str(s["top"]), str(s["dirs"][r] / "whole.xtc"))
        l = u.select_atoms(f"resname {s['resname']} and not name H*")
        gh = mda.AtomGroup([u.atoms[a.index] for a in gap_heavy])
        gc = mda.AtomGroup([u.atoms[a.index] for a in gap_ca])
        nn = mda.AtomGroup([u.atoms[a.index] for a in n655])
        mh, c0, c1 = [], [], []
        for ts in u.trajectory[::5]:
            if ts.time > P.WINDOW_PS + 1e-3:
                break
            d = dists(l, gh, gc, nn, ts.dimensions)
            mh.append(d["min_heavy_to_gap_residues_A"]); c0.append(d["centroid_to_gap_CA_A"][0]); c1.append(d["centroid_to_gap_CA_A"][1])
        traj[f"r{r}"] = {"frames": len(mh), "min_heavy_to_gap_min_A": round(min(mh), 1), "min_heavy_to_gap_median_A": round(float(np.median(mh)), 1),
                         "centroid_to_gap_CA_min_A": round(min(min(c0), min(c1)), 1),
                         "centroid_to_gap_CA_median_A": [round(float(np.median(c0)), 1), round(float(np.median(c1)), 1)]}
    out[system] = {"start": {k: (round(v, 1) if isinstance(v, float) else [round(x, 1) for x in v] if isinstance(v, list) and v and isinstance(v[0], float) else v) for k, v in start.items()}, "trajectory": traj}
    print(system, json.dumps(out[system], indent=1))
json.dump(out, open(Path(__file__).resolve().parent / "measure_icl3_distance.json", "w"), indent=1)

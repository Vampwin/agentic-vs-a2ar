"""POST HOC check (ADR-AIIM-394): frames of whole.xtc in which the receptor C-alpha atoms do not superpose on the reference (image jumps of one chain).
Per frame: RMSD of the C-alpha atoms after least-squares superposition on the pre-equilibration structure. A bound receptor gives 1 to 2 A; a
chain placed in another periodic image gives tens of A. Frames above 5 A are flagged; for each flagged frame the ligand RMSD stored by md_11 is read."""
import json, sys
from pathlib import Path
import MDAnalysis as mda, numpy as np
from MDAnalysis.analysis import align
sys.path.insert(0, str(Path(__file__).resolve().parent))
import posthoc_rmsf_hbond as P

out = {}
for system, s in P.SYSTEMS.items():
    ref = mda.Universe(str(s["top"]))
    for r in (1, 2, 3):
        adir = s["dirs"][r]
        ca = np.array(P.read_ndx(adir / "analysis.ndx")["Calpha"]) - 1
        u = mda.Universe(str(s["top"]), str(adir / "whole.xtc"))
        rp = ref.atoms[ca].positions; rc = rp.mean(0)
        t, rm = [], []
        for ts in u.trajectory:
            if ts.time > P.WINDOW_PS + 1e-3:
                break
            m = u.atoms[ca].positions
            _, v = align.rotation_matrix(m - m.mean(0), rp - rc)
            t.append(ts.time / 1000); rm.append(v)
        t, rm = np.array(t), np.array(rm)
        rt = np.array([[float(x) for x in ln.split()[:2]] for ln in (adir / "ligand_rmsd.xvg").read_text().splitlines() if ln.strip() and ln[0] not in "#@"])
        bad = rm > 5.0
        rr = rt[:len(t), 1] * 10.0                     # ligand RMSD in A, same frames as t (window of 50 ns)
        keep = ~bad
        corrected = {"rmsd_mean_A_all_frames": round(float(rr.mean()), 2), "rmsd_mean_A_without_flagged": round(float(rr[keep].mean()), 2),
                     "frac_lt_4A_all_frames": round(float((rr < 4.0).mean()), 4), "frac_lt_4A_without_flagged": round(float((rr[keep] < 4.0).mean()), 4)}
        lig = [round(float(rt[np.argmin(abs(rt[:, 0] - tt)), 1] * 10), 1) for tt in t[bad]]
        out[f"{system}_r{r}"] = {"frames": int(len(t)), "flagged_gt_5A": int(bad.sum()), "flagged_times_ns": [round(float(x), 2) for x in t[bad]][:12], "flagged_times_ns_all": [round(float(x), 3) for x in t[bad]],
                                 "calpha_fit_rmsd_flagged_A": [round(float(x), 1) for x in rm[bad]][:12],
                                 "ligand_rmsd_at_flagged_A": lig[:12], "median_calpha_fit_rmsd_A": round(float(np.median(rm)), 2),
                                 "max_unflagged_A": round(float(rm[~bad].max()), 2), "ligand_rmsd_effect": corrected}
        print(f"{system}_r{r}", out[f"{system}_r{r}"])
json.dump(out, open(Path(__file__).resolve().parent / "check_image_jumps.json", "w"), indent=1)

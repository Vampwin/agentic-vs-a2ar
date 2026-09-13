"""md_17_pocket_from_start.py -- is the ligand still in the pocket it was docked into?

The RMSD to the docked pose could not be made trustworthy for NPLIB01185. Its
time series jumps by up to 56 Å between frames 50 ps apart — physically
impossible — and in every one of the 400 frames reading over 40 Å, `gmx mindist`
puts the ligand 3.0 to 4.7 Å from the receptor, in contact. Three trjconv
treatments were tried and none removed it: `-pbc mol -center`, adding
`-ur compact`, and `-pbc nojump`; `-pbc cluster` refuses the C-alpha group
because it is not whole molecules.

So the RMSD is abandoned for that compound rather than filtered, because
filtering a bimodal distribution down to its low mode is choosing the answer.

This measures something the periodic image cannot corrupt and that answers the
question directly: the distance from the ligand to the **binding site as the
docked pose itself defines it** — every receptor C-alpha within 10 Å of the
ligand in the energy-minimised starting structure. `gmx mindist` applies the
minimum-image convention itself, so no trajectory representation enters.

A ligand that stays in its pocket keeps this distance at a few ångström. A
ligand that migrates elsewhere on the protein keeps `mindist`-to-receptor small
while this grows. The two measurements together separate the cases that the
RMSD conflated.

Usage:
    python md_17_pocket_from_start.py [--cutoff 10.0]
"""
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path
from kdss_paths import project_root, resolve, under, wsl_root


def _root() -> Path:
    # Delegates to kdss_paths, which finds whichever layout this file is inside
    # -- the development tree or a deposit -- by walking up from its own
    # location. It previously tried two absolute paths and exited if neither
    # existed, which on the author's machine meant a copy of the deposit read
    # and wrote the development tree instead of itself.
    return project_root()


ROOT = _root()
BASE = resolve("20_new_runs/md_candidates_2026-08-22")
COMPOUNDS = ["NPLIB01185", "NPLIB01643", "NPLIB03056"]
REPS = [1, 2, 3]
GMX = "/opt/gromacs/2026.3-cuda13.1-sm120/bin/gmx"


def read_gro(path: Path):
    """(atom_index_1based, resname, atomname, x, y, z) in nm."""
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    n = int(lines[1])
    out = []
    for i in range(n):
        L = lines[2 + i]
        out.append((i + 1, L[5:10].strip(), L[10:15].strip(),
                    float(L[20:28]), float(L[28:36]), float(L[36:44])))
    return out


def pocket_indices(atoms, cutoff_nm: float) -> list[int]:
    lig = [(a[3], a[4], a[5]) for a in atoms if a[1] == "LIG"]
    if not lig:
        raise SystemExit("no LIG atoms in the structure")
    ca = [a for a in atoms if a[2] == "CA"]
    keep = []
    c2 = cutoff_nm * cutoff_nm
    for idx, _rn, _an, x, y, z in ca:
        for lx, ly, lz in lig:
            if (x - lx) ** 2 + (y - ly) ** 2 + (z - lz) ** 2 <= c2:
                keep.append(idx)
                break
    return keep


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cutoff", type=float, default=10.0, help="angstrom")
    args = ap.parse_args()
    cutoff_nm = args.cutoff / 10.0
    result = {}

    for cid in COMPOUNDS:
        run = BASE / cid / "04_run"
        gro = run / "em.gro"
        if not gro.is_file():
            gro = BASE / cid / "03_build" / "system_gmx.gro"
        atoms = read_gro(gro)
        pocket = pocket_indices(atoms, cutoff_nm)
        lig = [a[0] for a in atoms if a[1] == "LIG"]
        print(f"{cid}: pocket = {len(pocket)} C-alpha within {args.cutoff:.0f} Å "
              f"of the docked ligand ({len(lig)} atoms), from {gro.name}")

        for r in REPS:
            d = run / f"analysis_prod_r{r}"
            d.mkdir(parents=True, exist_ok=True)
            ndx = d / "pocket_start.ndx"
            ndx.write_text(
                "[ Pocket ]\n"
                + "\n".join(" ".join(str(i) for i in pocket[k:k + 15])
                            for k in range(0, len(pocket), 15))
                + "\n[ Ligand ]\n"
                + "\n".join(" ".join(str(i) for i in lig[k:k + 15])
                            for k in range(0, len(lig), 15)) + "\n",
                encoding="utf-8")
            xvg = d / "pocket_start_dist.xvg"
            if not xvg.is_file():
                cmd = ("cd {wsl_root()!r} && printf 'Pocket\\nLigand\\n' | "
                       f"{GMX} mindist -s '{(run / f'prod_r{r}.tpr').relative_to(ROOT).as_posix()}' "
                       f"-f '{(run / f'prod_r{r}.xtc').relative_to(ROOT).as_posix()}' "
                       f"-n '{ndx.relative_to(ROOT).as_posix()}' "
                       f"-od '{xvg.relative_to(ROOT).as_posix()}' -tu ns "
                       f"> '{(d / 'pocket_start.log').relative_to(ROOT).as_posix()}' 2>&1")
                subprocess.run(["wsl", "-d", "A2A-MD", "--", "bash", "-c", cmd],
                               capture_output=True, text=True, timeout=3600)
            if not xvg.is_file():
                print(f"  {cid} r{r}: mindist FAILED")
                continue
            v = [float(l.split()[1]) * 10.0
                 for l in xvg.read_text(encoding="utf-8").splitlines()
                 if l and l[0] not in "@#"]
            v_sorted = sorted(v)
            med = v_sorted[len(v) // 2]
            result.setdefault(cid, {})[f"r{r}"] = {
                "n_frames": len(v), "min_A": round(min(v), 2),
                "median_A": round(med, 2), "max_A": round(max(v), 2),
                "final_A": round(v[-1], 2),
                "frames_beyond_10A": sum(1 for x in v if x > 10.0),
                "fraction_beyond_10A": round(sum(1 for x in v if x > 10.0) / len(v), 4),
            }
            s = result[cid][f"r{r}"]
            print(f"  r{r}: min {s['min_A']:.2f}  median {s['median_A']:.2f}  "
                  f"max {s['max_A']:.2f}  final {s['final_A']:.2f} Å   "
                  f"beyond 10 Å: {s['frames_beyond_10A']} of {s['n_frames']}")

    out = resolve("20_new_runs/md_pose_persistence_2026-08-31")
    out.mkdir(parents=True, exist_ok=True)
    (out / "pocket_occupancy.json").write_text(json.dumps({
        "question": "Is the ligand still in the site it was docked into?",
        "pocket_definition": f"receptor C-alpha atoms within {args.cutoff:.0f} Å of "
                             "the ligand in the energy-minimised starting structure "
                             "— the binding site as the docked pose itself defines it",
        "why_not_rmsd": "The RMSD to the docked pose is unusable for NPLIB01185: "
                        "its time series jumps up to 56 Å between frames 50 ps "
                        "apart, and the high-RMSD frames are in van der Waals "
                        "contact with the receptor. Four trjconv treatments failed "
                        "to remove it. Filtering to the low mode would be choosing "
                        "the answer.",
        "measured_by": "gmx mindist, which applies the minimum-image convention "
                       "itself, on the raw trajectory",
        "by_compound": result,
    }, indent=2), encoding="utf-8")
    print(f"\nwritten -> {out / 'pocket_occupancy.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

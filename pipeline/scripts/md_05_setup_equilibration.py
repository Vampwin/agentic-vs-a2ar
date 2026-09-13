"""
md_05_setup_equilibration.py -- prepare the staged equilibration for the A2AR
membrane system.

Three things have to be built before GROMACS can run a restrained membrane
equilibration, and none of them survives being guessed at:

1. **Position restraints per molecule type.** The topology parmed wrote has
   eight moleculetypes; restraints belong only on the receptor's two protein
   segments and on the ligand. GROMACS needs a separate `#ifdef POSRES`
   include inside each of those blocks, with atom numbers local to that
   molecule -- not global. Written here from the topology's own `[ atoms ]`
   masses, so an atom counts as heavy if the force field says it is.

2. **Three temperature-coupling groups.** The benchmark run coupled `System`
   as one group and said in its own comment that production must not. A
   bilayer, its solvent and the solute have different heat capacities and
   relaxation times; one thermostat over all three lets the membrane and water
   exchange energy through the coupling rather than through physics. Groups
   here are solute (protein + ligand + the crystallographic Na+), membrane
   (PA + PC + OL, the three Lipid21 residues of POPC), and solvent (water +
   added ions).

3. **Semi-isotropic pressure coupling.** Isotropic coupling on a bilayer
   forces the membrane area and the box height to scale together, which
   squeezes or stretches the lipids away from their equilibrium area per
   lipid. xy and z are coupled separately.

Restraints are released in stages rather than at once: 1000 -> 200 -> 50 -> 0
kJ/mol/nm^2. Dropping straight to zero after minimisation is how a membrane
system gets a lipid tail threaded through an aromatic ring.

`refcoord-scaling = com` is set wherever restraints and pressure coupling run
together. Without it the reference coordinates stay fixed while the box
scales, and the restraints slowly pull the solute against the barostat.

Usage:
    python md_05_setup_equilibration.py
"""
from __future__ import annotations

import re
from pathlib import Path
from kdss_paths import project_root, resolve, under

RESTRAINED = ("system1", "system2", "ZMA")   # receptor segments and ligand
FC = 1000                                     # kJ/mol/nm^2, stage 1


def _root() -> Path:
    # Delegates to kdss_paths, which finds whichever layout this file is inside
    # -- the development tree or a deposit -- by walking up from its own
    # location. It previously tried two absolute paths and exited if neither
    # existed, which on the author's machine meant a copy of the deposit read
    # and wrote the development tree instead of itself.
    return project_root()


ROOT = _root()
RUN = resolve("20_new_runs/md_production_2026-08-21")
BUILD = under("20_new_runs/md_production_2026-08-21", "03_build")
MDP = under("20_new_runs/md_production_2026-08-21", "04_run")


def split_moleculetypes(lines: list[str]) -> list[tuple[str, int, int]]:
    """(name, start_index, end_index) for every [ moleculetype ] block."""
    heads = [i for i, l in enumerate(lines) if re.match(r"^\s*\[\s*moleculetype\s*\]", l)]
    out = []
    for k, h in enumerate(heads):
        name = None
        for l in lines[h + 1:]:
            t = l.split(";")[0].strip()
            if t:
                name = t.split()[0]
                break
        end = heads[k + 1] if k + 1 < len(heads) else next(
            (i for i, l in enumerate(lines[h:], h)
             if re.match(r"^\s*\[\s*system\s*\]", l)), len(lines))
        out.append((name, h, end))
    return out


def heavy_atoms(lines: list[str], start: int, end: int) -> list[int]:
    """Local atom numbers with mass > 2 amu inside one moleculetype block."""
    block = lines[start:end]
    ai = next((i for i, l in enumerate(block)
               if re.match(r"^\s*\[\s*atoms\s*\]", l)), None)
    if ai is None:
        return []
    out = []
    for l in block[ai + 1:]:
        if re.match(r"^\s*\[", l):
            break
        t = l.split(";")[0].split()
        if len(t) >= 8:
            try:
                if float(t[7]) > 2.0:
                    out.append(int(t[0]))
            except ValueError:
                continue
    return out


def main() -> int:
    top = BUILD / "system_gmx.top"
    lines = top.read_text(encoding="utf-8").split("\n")

    if "#ifdef POSRES" in "\n".join(lines):
        print("topology already carries POSRES includes; leaving it alone")
    else:
        blocks = split_moleculetypes(lines)
        found = {n for n, _, _ in blocks}
        missing = [m for m in RESTRAINED if m not in found]
        if missing:
            raise SystemExit(f"moleculetypes not found in the topology: {missing}")

        # Insert from the bottom up so earlier indices stay valid.
        for name, h, end in sorted(blocks, key=lambda b: -b[1]):
            if name not in RESTRAINED:
                continue
            atoms = heavy_atoms(lines, h, end)
            if not atoms:
                raise SystemExit(f"no heavy atoms parsed for {name}")
            itp = BUILD / f"posre_{name}.itp"
            body = [f"; heavy-atom position restraints for {name}",
                    "; generated by md_05_setup_equilibration.py from the topology's",
                    "; own masses; force constant is set per stage via -DPOSRES_FC",
                    "[ position_restraints ]",
                    "; atom  type       fx       fy       fz"]
            body += [f"{a:6d}     1  POSRES_FC POSRES_FC POSRES_FC" for a in atoms]
            itp.write_text("\n".join(body) + "\n", encoding="utf-8")
            lines[end:end] = ["", "#ifdef POSRES",
                              f'#include "posre_{name}.itp"', "#endif", ""]
            print(f"  {name:9s} {len(atoms):5d} heavy atoms restrained -> {itp.name}")
        top.write_text("\n".join(lines), encoding="utf-8")
        print(f"POSRES includes written into {top.name}")

    # ---- index file: three coupling groups ---------------------------------
    ndx_cmds = (
        "1 | 13 | 14\n"        # Protein | ZMA | NA          -> Solute
        "15 | 16 | 17\n"       # PA | PC | OL                -> Membrane
        "21 | 18 | 19\n"       # Water | Na+ | Cl-           -> Solvent
        "name 24 Solute\nname 25 Membrane\nname 26 Solvent\nq\n"
    )
    (MDP / "make_ndx.in").write_text(ndx_cmds, encoding="utf-8")

    # ---- mdp stages ---------------------------------------------------------
    common = """cutoff-scheme           = Verlet
nstlist                 = 20
coulombtype             = PME
rcoulomb                = 1.0
rvdw                    = 1.0
DispCorr                = EnerPres
constraints             = h-bonds
constraint-algorithm    = LINCS
lincs-iter              = 1
lincs-order             = 4
pbc                     = xyz
; Three baths, not one. A bilayer, its solvent and the solute relax on
; different timescales; a single thermostat over all three lets them trade
; energy through the coupling instead of through physics.
tcoupl                  = V-rescale
tc-grps                 = Solute Membrane Solvent
tau-t                   = 1.0 1.0 1.0
ref-t                   = 310 310 310
nstlog                  = 5000
nstenergy               = 5000
nstxout-compressed      = 25000
"""
    npt = """; Semi-isotropic: xy (membrane plane) and z scale independently. Isotropic
; coupling would force area-per-lipid to follow box height.
pcoupl                  = C-rescale
pcoupltype              = semiisotropic
tau-p                   = 5.0
ref-p                   = 1.0 1.0
compressibility         = 4.5e-5 4.5e-5
; With restraints under a barostat the reference coordinates must scale with
; the box, or the restraints fight the barostat.
refcoord-scaling        = com
"""
    stages = [
        ("nvt",   "define = -DPOSRES -DPOSRES_FC=1000\n", 125000, True,  False,
         "250 ps NVT, heavy atoms restrained at 1000 kJ/mol/nm^2"),
        ("npt1",  "define = -DPOSRES -DPOSRES_FC=1000\n", 500000, False, True,
         "1 ns NPT, restraints still 1000; the box finds its area per lipid"),
        ("npt2",  "define = -DPOSRES -DPOSRES_FC=200\n",  1000000, False, True,
         "2 ns NPT, restraints relaxed to 200"),
        ("npt3",  "define = -DPOSRES -DPOSRES_FC=50\n",   1000000, False, True,
         "2 ns NPT, restraints relaxed to 50"),
        ("prod",  "",                                     50000000, False, True,
         "100 ns production, unrestrained"),
    ]
    for name, define, nsteps, genvel, pressure, note in stages:
        txt = [f"; {note}", "integrator              = md", "dt                      = 0.002",
               f"nsteps                  = {nsteps}"]
        if define:
            txt.append(define.rstrip())
        txt.append(common.rstrip())
        if pressure:
            txt.append(npt.rstrip())
        else:
            txt.append("pcoupl                  = no")
        if genvel:
            txt += ["gen-vel                 = yes", "gen-temp                = 310",
                    "gen-seed                = 20260822"]
        else:
            txt += ["gen-vel                 = no", "continuation            = yes"]
        if name == "prod":
            # Override the shared value rather than appending a second one --
            # grompp rejects a doubly-defined parameter outright, which is how
            # this was found.
            txt = [line.replace("nstxout-compressed      = 25000",
                                "nstxout-compressed      = 50000")
                   for line in txt]
        (MDP / f"{name}.mdp").write_text("\n".join(txt) + "\n", encoding="utf-8")
        print(f"  {name}.mdp  {nsteps * 0.002 / 1000:g} ns  {note}")

    print(f"\nmdp files and make_ndx input written to {MDP}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

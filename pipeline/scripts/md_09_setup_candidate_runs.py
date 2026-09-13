"""
md_09_setup_candidate_runs.py -- position restraints, coupling groups and mdp
files for one candidate, with three production replicates.

Same three prerequisites as the control (md_05_setup_equilibration.py), and
the same reasoning behind each: per-moleculetype position restraints written
from the topology's own masses, three temperature-coupling baths rather than
one, and semi-isotropic pressure coupling for a bilayer.

Two differences from the control, both deliberate:

**Three replicates of 100 ns each, matching the control.** A single
trajectory cannot distinguish a pose that is stable from one that happened not
to move; replicates are the only thing that licenses a statement about
consistency. 100 ns is still shorter than the 300 ns triplicate precedent for
this receptor (Ng et al. 2014), and that limit is reported rather than
glossed. It was raised from an initial 50 ns after checking, not assuming,
what duration this kind of check actually needs: two independent studies
(Bhattarai et al., PLOS Comput Biol 2021, doi:10.1371/journal.pcbi.1008936;
Liu & Watanabe, J Comput Aided Mol Des 2017, doi:10.1007/s10822-016-0005-2)
report that an incorrect pose can remain stable for hundreds of ns to
microseconds -- the earlier "50 ns is enough because a wrong pose leaves
within 20-30 ns" reasoning had no citation behind it and does not survive
contact with the literature. Matching the control's full 100 ns removes an
asymmetry a reviewer would otherwise have to take on faith, without
pretending the underlying limitation -- that pose persistence under MD is
weak evidence of correctness at any duration tested here -- goes away.

**Replicates differ only in the velocity seed**, generated fresh at the start
of each production run from the same equilibrated structure. That isolates the
one source of divergence being sampled: thermal noise, not a different
starting geometry.

Usage:
    python md_09_setup_candidate_runs.py NPLIB01185
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path
from kdss_paths import project_root, resolve, under

RESTRAINED = ("system1", "system2", "LIG")
N_REPLICATES = 3
PROD_NS = 100
SEEDS = (20260822, 20260823, 20260824)


def _root() -> Path:
    # Delegates to kdss_paths, which finds whichever layout this file is inside
    # -- the development tree or a deposit -- by walking up from its own
    # location. It previously tried two absolute paths and exited if neither
    # existed, which on the author's machine meant a copy of the deposit read
    # and wrote the development tree instead of itself.
    return project_root()


ROOT = _root()
DEST = resolve("20_new_runs/md_candidates_2026-08-22")

MEMBRANE = {"PA", "PC", "OL"}
SOLVENT = {"WAT", "HOH", "SOL", "Na+", "Cl-", "NA+", "CL-"}


def split_moleculetypes(lines: list[str]) -> list[tuple[str, int, int]]:
    heads = [i for i, l in enumerate(lines)
             if re.match(r"^\s*\[\s*moleculetype\s*\]", l)]
    out = []
    for k, h in enumerate(heads):
        name = next((l.split(";")[0].strip().split()[0]
                     for l in lines[h + 1:] if l.split(";")[0].strip()), None)
        end = heads[k + 1] if k + 1 < len(heads) else next(
            (i for i, l in enumerate(lines[h:], h)
             if re.match(r"^\s*\[\s*system\s*\]", l)), len(lines))
        out.append((name, h, end))
    return out


def heavy_atoms(lines: list[str], start: int, end: int) -> list[int]:
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
    ap = argparse.ArgumentParser()
    ap.add_argument("compound_id")
    args = ap.parse_args()
    cid = args.compound_id
    run = DEST / cid
    build, mdp = run / "03_build", run / "04_run"
    mdp.mkdir(parents=True, exist_ok=True)

    top = build / "system_gmx.top"
    if not top.is_file():
        raise SystemExit(f"{cid}: no system_gmx.top -- has the build finished?")
    lines = top.read_text(encoding="utf-8").split("\n")

    # ---- position restraints ------------------------------------------------
    if "#ifdef POSRES" in "\n".join(lines):
        print(f"  {cid}: topology already carries POSRES includes")
    else:
        blocks = split_moleculetypes(lines)
        found = {n for n, _, _ in blocks}
        missing = [m for m in RESTRAINED if m not in found]
        if missing:
            raise SystemExit(f"{cid}: moleculetypes not in topology: {missing} "
                             f"(present: {sorted(found)})")
        total = 0
        for name, h, end in sorted(blocks, key=lambda b: -b[1]):
            if name not in RESTRAINED:
                continue
            atoms = heavy_atoms(lines, h, end)
            if not atoms:
                raise SystemExit(f"{cid}: no heavy atoms parsed for {name}")
            total += len(atoms)
            itp = build / f"posre_{name}.itp"
            body = [f"; heavy-atom position restraints for {name}",
                    "; force constant set per stage via -DPOSRES_FC",
                    "[ position_restraints ]",
                    "; atom  type       fx       fy       fz"]
            body += [f"{a:6d}     1  POSRES_FC POSRES_FC POSRES_FC" for a in atoms]
            itp.write_text("\n".join(body) + "\n", encoding="utf-8")
            lines[end:end] = ["", "#ifdef POSRES",
                              f'#include "posre_{name}.itp"', "#endif", ""]
        top.write_text("\n".join(lines), encoding="utf-8")
        print(f"  {cid}: {total} heavy atoms restrained across {len(RESTRAINED)} molecules")

    # ---- coupling groups, built from residue names -------------------------
    # Not via gmx make_ndx: its group renumbering during merges produced a
    # "Solute" group containing the solvent when the control was set up.
    gro = (build / "system_gmx.gro").read_text(encoding="utf-8").split("\n")
    n = int(gro[1])
    groups: dict[str, list[int]] = {"Solute": [], "Membrane": [], "Solvent": []}
    for i in range(n):
        rn = gro[2 + i][5:10].strip()
        key = ("Membrane" if rn in MEMBRANE
               else "Solvent" if rn in SOLVENT else "Solute")
        groups[key].append(i + 1)
    out = []
    for name, ids in groups.items():
        out.append(f"[ {name} ]")
        for j in range(0, len(ids), 15):
            out.append(" ".join(f"{a:6d}" for a in ids[j:j + 15]))
    (mdp / "coupling.ndx").write_text("\n".join(out) + "\n", encoding="utf-8")
    tot = sum(len(v) for v in groups.values())
    if tot != n:
        raise SystemExit(f"{cid}: coupling groups cover {tot} of {n} atoms")
    print(f"  {cid}: Solute {len(groups['Solute'])} / Membrane "
          f"{len(groups['Membrane'])} / Solvent {len(groups['Solvent'])} = {tot}")

    # ---- mdp files ----------------------------------------------------------
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
tcoupl                  = V-rescale
tc-grps                 = Solute Membrane Solvent
tau-t                   = 1.0 1.0 1.0
ref-t                   = 310 310 310
nstlog                  = 5000
nstenergy               = 5000
nstxout-compressed      = 25000
"""
    npt = """pcoupl                  = C-rescale
pcoupltype              = semiisotropic
tau-p                   = 5.0
ref-p                   = 1.0 1.0
compressibility         = 4.5e-5 4.5e-5
refcoord-scaling        = com
"""
    em = """; Steepest descent to the same tolerance as the control build.
integrator      = steep
nsteps          = 5000
emtol           = 1000.0
emstep          = 0.01
cutoff-scheme   = Verlet
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0
DispCorr        = EnerPres
pbc             = xyz
"""
    (mdp / "em.mdp").write_text(em, encoding="utf-8")

    stages = [
        ("nvt",  "define = -DPOSRES -DPOSRES_FC=1000", 125000, True, False),
        ("npt1", "define = -DPOSRES -DPOSRES_FC=1000", 500000, False, True),
        ("npt2", "define = -DPOSRES -DPOSRES_FC=200", 1000000, False, True),
        ("npt3", "define = -DPOSRES -DPOSRES_FC=50", 1000000, False, True),
    ]
    for name, define, nsteps, genvel, pressure in stages:
        t = ["integrator              = md", "dt                      = 0.002",
             f"nsteps                  = {nsteps}", define, common.rstrip()]
        t.append(npt.rstrip() if pressure else "pcoupl                  = no")
        t += (["gen-vel                 = yes", "gen-temp                = 310",
               "gen-seed                = 20260822"] if genvel
              else ["gen-vel                 = no", "continuation            = yes"])
        (mdp / f"{name}.mdp").write_text("\n".join(t) + "\n", encoding="utf-8")

    # Replicates: fresh velocities from the same equilibrated structure, so the
    # only thing that differs between them is thermal noise.
    for r, seed in zip(range(1, N_REPLICATES + 1), SEEDS):
        t = [f"; replicate {r} of {N_REPLICATES}: {PROD_NS} ns, unrestrained",
             "integrator              = md", "dt                      = 0.002",
             f"nsteps                  = {int(PROD_NS * 1000 / 0.002)}",
             common.replace("nstxout-compressed      = 25000",
                            "nstxout-compressed      = 25000    ; 50 ps/frame").rstrip(),
             npt.replace("refcoord-scaling        = com\n", "").rstrip(),
             "gen-vel                 = yes", "gen-temp                = 310",
             f"gen-seed                = {seed}"]
        (mdp / f"prod_r{r}.mdp").write_text("\n".join(t) + "\n", encoding="utf-8")

    print(f"  {cid}: em + 4 equilibration stages + {N_REPLICATES} x {PROD_NS} ns "
          f"replicates (seeds {', '.join(str(s) for s in SEEDS)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

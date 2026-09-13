"""crossdock_zma.py -- does redocking accuracy survive a change of receptor?

Replaces `crossdock_validation.py`, whose output directory was voided. That
script was run three times with three different receptor preparations into one
directory, with `rm -rf` between passes, and the result was a directory whose
files described different experiments and could not be told apart from inside
it. See `20_new_runs/crossdock_validation_2026-09-03/VOID.md`.

**Three rules are enforced here in code, not by intention.**

1. `--out` must not exist, or must be empty. A run directory is written once.
   `20_new_runs/` is immutable under ADR-AIIM-092 and the previous script broke
   that rule three times in eight minutes.
2. Every external input is copied in and hashed into `manifest.json` before
   anything is computed. Nothing is read from another run directory implicitly.
   The voided run pulled a receptor out of `snp_screen_2026-08-04` by hand and
   never declared it.
3. **The within-study control gates the result.** ZM241385 is redocked into
   4EIY, the structure it was solved in, under this run's own preparation. If
   that does not land within `CONTROL_MAX_RMSD_A`, the script writes the control
   and refuses to report any cross-dock number. The preparation, not the
   receptor conformation, is then what is being measured.

**Preparation is the production G2 recipe, transcribed and reused, not
approximated.** `_write_receptor_with_waters()` from `run_snp_screen.py` selects
the protein chain plus the ordered waters within 5 A of the co-crystal ligand;
pdbfixer completes missing atoms and adds hydrogens at pH 7.4; obabel converts
to pdbqt with `-xr` and no `-p`. The waters are not a detail: that function's
own docstring records that dropping them moves rank-1 redocking from 0.28 A to
8.67 A, and the voided run reproduced 8.478 A by dropping them.

**5OLH is excluded.** It was solved with 9XT, so it has no ZM241385 of its own:
its reference pose would have to be constructed rather than observed, and the
water shell would be selected around a different ligand. **3EML is included but
is not a new result** -- ADR-AIIM-036 excluded it on 2026-08-04 for failing
redocking at every water shell tested, and any outcome here confirms that
decision rather than discovering it.

Usage:
    python crossdock_zma.py --out 20_new_runs/crossdock_zma_2026-09-03
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import shlex
import shutil
import subprocess
from pathlib import Path
from kdss_paths import project_root, resolve, under, wsl_root

import numpy as np

# --- pre-registered, ADR-AIIM-162 -------------------------------------------
RMSD_PASS_A = 2.0             # this project's frozen G3 criterion
CONTROL_MAX_RMSD_A = 1.0      # the control must clearly succeed, not merely pass
WATER_CUTOFF_A = 5.0          # frozen config, receptor.structural_waters
PROTONATION_PH = 7.4          # frozen config, receptor.protonation_ph
N_CONFORMERS, EXHAUSTIVENESS, N_POSES, SEED = 8, 32, 10, 20260729
CONTROL = "4EIY"
RECEPTORS = ["3EML", "4EIY", "5OLG", "6MH8"]
CONFIRMS_PRIOR_DECISION = {"3EML": "ADR-AIIM-036 excluded it on 2026-08-04"}
LIGAND_RESN = "ZMA"
RECEPTOR_RESI_MAX = 999       # fusion partners are numbered from 1000: T4
                              # lysozyme in 3EML, BRIL in 4EIY. Matching CA by
                              # residue number across that boundary superposed
                              # 3EML at 8.94 A in the voided run.


def _root() -> Path:
    # Delegates to kdss_paths, which finds whichever layout this file is inside
    # -- the development tree or a deposit -- by walking up from its own
    # location. It previously tried two absolute paths and exited if neither
    # existed, which on the author's machine meant a copy of the deposit read
    # and wrote the development tree instead of itself.
    return project_root()


ROOT = _root()
SCRIPTS = resolve("10_new_system/scripts")
STRUCT = resolve("02_frozen_inputs/a2ar_structures")
BOX_SRC = resolve("20_new_runs/snp_screen_2026-08-04/G2_preparation/docking_box.json")
TEMPLATE_SRC = resolve("20_new_runs/snp_screen_2026-08-04/G3_redocking/reference_template.json")
PDBFIXER_ENV = "mdenv"   # where run_snp_screen.py says PDBFixer lives
WSL_DISTRO, CONDA_ENV = "Ubuntu-22.04", "diffdock"


def _wsl_path(p) -> str:
    """Windows path -> the /mnt/<drive> form the WSL side needs."""
    s = str(p).replace(chr(92), "/")
    return f"/mnt/{s[0].lower()}{s[2:]}" if s[1:3] == ":/" else s


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


MR = _load("multi_reference_ifp_variant", under("10_new_system/scripts", "multi_reference_ifp_variant.py"))
RSS = _load("run_snp_screen", under("10_new_system/scripts", "run_snp_screen.py"))


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def wsl(cmd, timeout=7200):
    full = (f"source ~/miniconda3/etc/profile.d/conda.sh && conda activate {CONDA_ENV} "
            f"&& cd {wsl_root()!r} && {cmd}")
    return subprocess.run(["wsl", "-d", WSL_DISTRO, "--", "bash", "-lc", full],
                          capture_output=True, text=True, timeout=timeout)


def superpose_pdb(src: Path, dst: Path, R, t):
    """Whole file, every record, mapped by (R, t). Waters move with the protein."""
    out = []
    for line in src.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            v = R @ np.array([float(line[30:38]), float(line[38:46]),
                              float(line[46:54])]) + t
            out.append(f"{line[:30]}{v[0]:8.3f}{v[1]:8.3f}{v[2]:8.3f}{line[54:]}")
        elif line.startswith(("TER", "END")):
            out.append(line)
    dst.write_text("\n".join(out) + "\n", encoding="utf-8", newline="\n")


def pdbfix(src: Path, dst: Path) -> str:
    """The production G2 sequence, in the environment production used.

    `run_snp_screen.py` says "PDBFixer lives in mdenv" and that WSL environment
    still exists. An earlier attempt here used a Windows conda environment that
    also has pdbfixer 1.12 and openmm 8.5.2; it loads force fields and reports
    all four platforms, but `addMissingAtoms()` kills the interpreter with exit
    127 the moment it has a side chain to rebuild. 3EML has no missing atoms and
    passed through it; 4EIY has seven and did not. Running where production ran
    reproduces production exactly -- 6,195 atoms and 15 retained waters for
    4EIY, matching `snp_screen_2026-08-04/G2_preparation/4EIY_fixed.pdb`.
    """
    code = (
        "from pdbfixer import PDBFixer\n"
        "from openmm.app import PDBFile\n"
        f"f = PDBFixer(filename={_wsl_path(src)!r})\n"
        "f.findMissingResidues(); f.missingResidues = {}\n"
        "f.findNonstandardResidues(); f.replaceNonstandardResidues()\n"
        "f.removeHeterogens(keepWater=True)\n"
        "f.findMissingAtoms(); f.addMissingAtoms()\n"
        f"f.addMissingHydrogens({PROTONATION_PH})\n"
        f"PDBFile.writeFile(f.topology, f.positions, open({_wsl_path(dst)!r}, 'w'))\n"
        "print('atoms:', sum(1 for _ in f.topology.atoms()), "
        "'waters:', sum(1 for r in f.topology.residues() if r.name in ('HOH','WAT')))\n"
    )
    r = subprocess.run(
        ["wsl", "-d", WSL_DISTRO, "--", "bash", "-lc",
         "source ~/miniconda3/etc/profile.d/conda.sh && conda activate "
         f"{PDBFIXER_ENV} && python -u -c {shlex.quote(code)}"],
        capture_output=True, text=True, timeout=1800)
    if not dst.is_file():
        raise RuntimeError(
            f"pdbfixer failed for {src.name}\nstdout: {r.stdout[-400:]}\n"
            f"stderr: {r.stderr[-400:]}")
    return r.stdout.strip().splitlines()[-1]


def crystal_ligand(src: Path, dst: Path, template_smiles: str):
    """The co-crystal ligand alone, primary altloc, bond orders from the frozen
    template. A PDB HETATM block has no bond orders and inferring them is
    guessing at chemistry."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    keep = [l for l in src.read_text(encoding="utf-8", errors="ignore").splitlines()
            if l.startswith("HETATM") and l[17:20].strip() == LIGAND_RESN
            and l[16] in (" ", "A")]
    dst.write_text("\n".join(keep) + "\nEND\n", encoding="utf-8", newline="\n")
    raw = Chem.MolFromPDBFile(str(dst), removeHs=True, sanitize=False)
    return AllChem.AssignBondOrdersFromTemplate(
        Chem.MolFromSmiles(template_smiles), raw)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    # Absolute from here on: every path below is expressed relative to ROOT for
    # the WSL side, and a relative --out makes that conversion raise.
    out = Path(args.out).resolve()

    # Rule 1, enforced.
    if out.exists() and any(out.iterdir()):
        raise SystemExit(
            f"REFUSED: {out} already exists and is not empty.\n"
            "A run directory is written once (ADR-AIIM-092). Choose a new one.")
    (out / "inputs").mkdir(parents=True, exist_ok=True)
    (out / "receptors").mkdir(exist_ok=True)

    # Rule 2, enforced: copy and hash every external input before computing.
    manifest = {"declared_inputs": {}, "recipe": {
        "preparation": "production G2, transcribed: _write_receptor_with_waters "
                       f"(chain A + waters within {WATER_CUTOFF_A} A of the "
                       "co-crystal ligand) -> pdbfixer (missing atoms, "
                       f"hydrogens at pH {PROTONATION_PH}, keepWater=True) -> "
                       "obabel -xr",
        "docking": f"{N_CONFORMERS} conformers, exhaustiveness {EXHAUSTIVENESS}, "
                   f"{N_POSES} poses, seed {SEED}",
        "excluded": {"5OLH": "solved with 9XT, so it has no ZM241385 of its own: "
                             "the reference pose would be constructed and the "
                             "water shell selected around a different ligand"},
        "confirms_prior_decision": CONFIRMS_PRIOR_DECISION,
    }}
    for src in [BOX_SRC, TEMPLATE_SRC] + [STRUCT / f"{p}.pdb" for p in RECEPTORS]:
        dst = out / "inputs" / src.name
        shutil.copy2(src, dst)
        manifest["declared_inputs"][src.name] = {
            "from": src.relative_to(ROOT).as_posix(), "sha256": sha256(src)}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2),
                                       encoding="utf-8")
    smiles = json.loads(TEMPLATE_SRC.read_text(encoding="utf-8"))["smiles"]

    # --- prepare ------------------------------------------------------------
    ref_ca = MR.ca_by_resi(MR.read_pdb(STRUCT / f"{CONTROL}.pdb"))
    prep = {}
    print("preparing receptors (production G2 recipe)")
    for pdb in RECEPTORS:
        src = STRUCT / f"{pdb}.pdb"
        ca = MR.ca_by_resi(MR.read_pdb(src))
        shared = sorted((r for r in set(ca) & set(ref_ca)
                         if int(r) <= RECEPTOR_RESI_MAX), key=int)
        P = np.array([ca[r] for r in shared]); Q = np.array([ref_ca[r] for r in shared])
        R, t = (MR.kabsch(P, Q) if pdb != CONTROL else (np.eye(3), np.zeros(3)))
        resid = float(np.sqrt((((P @ R.T + t) - Q) ** 2).sum(1).mean()))

        sup = out / "receptors" / f"{pdb}_super.pdb"
        superpose_pdb(src, sup, R, t)
        filt = out / "receptors" / f"{pdb}_protein_plus_waters.pdb"
        wat = RSS._write_receptor_with_waters(sup, filt, chain="A",
                                              cutoff=WATER_CUTOFF_A)
        fixed = out / "receptors" / f"{pdb}_fixed.pdb"
        fixinfo = pdbfix(filt, fixed)
        pdbqt = out / "receptors" / f"{pdb}_prepared.pdbqt"
        wsl(f'obabel -ipdb "{fixed.relative_to(ROOT).as_posix()}" -opdbqt '
            f'-O "{pdbqt.relative_to(ROOT).as_posix()}" -xr')
        lig = out / "receptors" / f"{pdb}_{LIGAND_RESN}_crystal.pdb"
        ref_mol = crystal_ligand(sup, lig, smiles)
        prep[pdb] = {"n_shared_CA": len(shared),
                     "superposition_rmsd_A": round(resid, 3),
                     "waters_retained": wat.get("n_waters_retained"),
                     "pdbfixer": fixinfo,
                     "reference_atoms": ref_mol.GetNumAtoms() if ref_mol else None,
                     "pdbqt_ok": pdbqt.is_file() and pdbqt.stat().st_size > 0}
        p = prep[pdb]
        print(f"  {pdb}: superposition {p['superposition_rmsd_A']:.3f} A, "
              f"waters {p['waters_retained']}, {p['pdbfixer']}, "
              f"reference {p['reference_atoms']} atoms")
    (out / "prep.json").write_text(json.dumps(prep, indent=2), encoding="utf-8")

    # --- dock ---------------------------------------------------------------
    ledger = out / "ligand_ledger.csv"
    ledger.write_text(f"compound_id,SMILES\nZM241385,{smiles}\n",
                      encoding="utf-8", newline="\n")
    box = out / "box.json"
    box.write_text(BOX_SRC.read_text(encoding="utf-8"), encoding="utf-8", newline="\n")
    print("\ndocking ZM241385 into each")
    for pdb in RECEPTORS:
        if not prep[pdb]["pdbqt_ok"]:
            print(f"  {pdb}: no receptor, skipped"); continue
        rec = (out / "receptors" / f"{pdb}_prepared.pdbqt").relative_to(ROOT).as_posix()
        print(f"  {pdb} ...", flush=True)
        wsl(f'python 10_new_system/scripts/g4_dock_library.py '
            f'--ledger "{ledger.relative_to(ROOT).as_posix()}" --receptor "{rec}" '
            f'--box "{box.relative_to(ROOT).as_posix()}" '
            f'--out "{(out / f"scores_{pdb}.csv").relative_to(ROOT).as_posix()}" '
            f'--poses-dir "{(out / f"poses_{pdb}").relative_to(ROOT).as_posix()}" '
            f'--id-col compound_id --smiles-col SMILES '
            f'--n-conformers {N_CONFORMERS} --exhaustiveness {EXHAUSTIVENESS} '
            f'--n-poses {N_POSES} --seed {SEED} --cpu 6')

    # --- score, with rule 3 enforced ----------------------------------------
    import warnings; warnings.filterwarnings("ignore")
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdMolAlign
    RDLogger.DisableLog("rdApp.*")
    res = {}
    for pdb in RECEPTORS:
        sdf = out / f"poses_{pdb}" / "ZM241385.sdf"
        lig = out / "receptors" / f"{pdb}_{LIGAND_RESN}_crystal.pdb"
        if not sdf.is_file():
            res[pdb] = {"available": False, "reason": "no poses written"}; continue
        ref = crystal_ligand(out / "receptors" / f"{pdb}_super.pdb", lig, smiles)
        ms = [m for m in Chem.SDMolSupplier(str(sdf), removeHs=True) if m]
        r = [round(rdMolAlign.CalcRMS(m, ref), 3) for m in ms]
        res[pdb] = {"available": True, "n_poses": len(r), "rank1_rmsd_A": r[0],
                    "best_rmsd_A": min(r), "all_rmsd_A": r,
                    "within_2A": r[0] <= RMSD_PASS_A,
                    "confirms": CONFIRMS_PRIOR_DECISION.get(pdb)}
        print(f"  {pdb}: rank-1 {r[0]} A  (best of {len(r)}: {min(r)})")

    ctrl = res.get(CONTROL, {})
    ctrl_ok = ctrl.get("available") and ctrl["rank1_rmsd_A"] <= CONTROL_MAX_RMSD_A
    report = {
        "question": "Does redocking accuracy survive a change of receptor conformation?",
        "authority": "ADR-AIIM-162",
        "recipe": manifest["recipe"],
        "control": {"receptor": CONTROL, "max_rmsd_A": CONTROL_MAX_RMSD_A,
                    "observed_rmsd_A": ctrl.get("rank1_rmsd_A"), "passed": bool(ctrl_ok)},
        "by_receptor": res,
    }
    if not ctrl_ok:
        report["outcome"] = "CONTROL_FAILED"
        report["cross_dock_results_withheld"] = (
            f"{CONTROL} redocked into its own structure at "
            f"{ctrl.get('rank1_rmsd_A')} A against a ceiling of {CONTROL_MAX_RMSD_A} A. "
            "The preparation, not the receptor conformation, is what this run "
            "measured. No cross-dock number here supports any claim.")
        print(f"\nCONTROL FAILED ({ctrl.get('rank1_rmsd_A')} A) -- results withheld")
    else:
        others = [p for p in RECEPTORS if p != CONTROL and res.get(p, {}).get("available")]
        n = sum(1 for p in others if res[p]["within_2A"])
        report["outcome"] = ("HOLDS" if n == len(others)
                             else "FAILS" if n == 0 else "DEGRADES")
        report["n_passing_of_cross"] = f"{n} of {len(others)}"
        print(f"\ncontrol {ctrl['rank1_rmsd_A']} A -- passed. "
              f"cross-dock {n} of {len(others)} within {RMSD_PASS_A} A "
              f"-> {report['outcome']}")
    (out / "crossdock_zma.json").write_text(json.dumps(report, indent=2),
                                            encoding="utf-8")
    print(f"\nwritten -> {out / 'crossdock_zma.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

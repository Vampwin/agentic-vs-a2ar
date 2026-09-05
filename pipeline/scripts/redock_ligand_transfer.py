"""redock_ligand_transfer.py -- does redocking accuracy survive a change of LIGAND?

`crossdock_zma.py` changed the receptor four times (4EIY, 5OLG, 3EML, 6MH8) and
never changed the ligand: all four structures were solved with ZM241385. That
study's own numbers -- 0.397, 0.318, 0.709 A -- are therefore statements about
one molecule redocked into crystals shaped around it.

This project already holds the reason to doubt they generalise. Under the same
production protocol, none of seven known A2A antagonists is placed in the Asn253
hydrogen bond that ZM241385 reproduces at 2.71 A (ADR-AIIM-163).

So: change the ligand. 5OLH holds **9XT**, a chemically distinct A2A antagonist,
in its own crystallographic frame. `crossdock_zma.py` excluded 5OLH because it
has no ZM241385 and a ZM241385 reference there would have to be constructed --
correct for that question, exactly backwards for this one. Redocking 9XT into
5OLH constructs nothing.

    arm       ligand      receptor   tests
    control   ZM241385    4EIY       the preparation, re-established here
    T1        9XT         5OLH       native redock, LIGAND changed
    T2        9XT         4EIY       ligand AND receptor changed

Criteria are fixed in ADR-AIIM-168, written before this script was run,
including the reading of each outcome.

**The three rules of `crossdock_zma.py` are carried over and enforced in code:**
write-once output directory (ADR-AIIM-092); every input copied in and hashed
into `manifest.json` before anything is computed; and the control gates the
report -- if ZM241385 does not redock into 4EIY within CONTROL_MAX_RMSD_A under
this run's own preparation, no other number is reported.

Usage:
    python redock_ligand_transfer.py --out 20_new_runs/redock_9xt_2026-09-05
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import shlex
import shutil
import subprocess
from pathlib import Path

import numpy as np

# --- pre-registered, ADR-AIIM-168 -------------------------------------------
RMSD_PASS_A = 2.0             # this project's frozen G3 criterion
CONTROL_MAX_RMSD_A = 1.0      # the control must clearly succeed, not merely pass
WATER_CUTOFF_A = 5.0          # frozen config, receptor.structural_waters
PROTONATION_PH = 7.4          # frozen config, receptor.protonation_ph
N_CONFORMERS, EXHAUSTIVENESS, N_POSES, SEED = 8, 32, 10, 20260729
FRAME = "4EIY"                # the frozen docking box is defined in this frame
RECEPTOR_RESI_MAX = 999       # fusion partners are numbered from 1000 (BRIL in
                              # 4EIY, T4 lysozyme in 3EML); matching CA across
                              # that boundary superposed 3EML at 8.94 A once.

# 9XT's template. Already in this repository as
# `np_library_prefilter.get_9xt_smiles()`, used for the §3.4 second-reference
# variant and asserted in `test_control_representativeness.py`. Verified against
# the 5OLH crystal composition: 24 heavy atoms, C16N7O, exactly.
SMILES_9XT = "Nc1nc(-c2ccco2)c2nnn(Cc3ccc(N)c(C)c3)c2n1"

# (arm, compound, crystal resname, source structure of the reference pose,
#  receptor docked into)
ARMS = [
    ("control", "ZM241385", "ZMA", "4EIY", "4EIY"),
    ("T1",      "9XT",      "9XT", "5OLH", "5OLH"),
    ("T2",      "9XT",      "9XT", "5OLH", "4EIY"),
]
# receptor -> the residue name of ITS OWN co-crystal ligand, which is what the
# structural water shell is selected around.
OWN_LIGAND = {"4EIY": "ZMA", "5OLH": "9XT"}


def _root() -> Path:
    for c in (Path(r"D:\Agentic VS Workflow"), Path("/mnt/d/Agentic VS Workflow")):
        if (c / "00_project_control/08_DECISION_LOG.md").is_file():
            return c
    raise SystemExit("project root not found")


ROOT = _root()
SCRIPTS = ROOT / "10_new_system/scripts"
STRUCT = ROOT / "02_frozen_inputs/a2ar_structures"
BOX_SRC = ROOT / "20_new_runs/snp_screen_2026-08-04/G2_preparation/docking_box.json"
TEMPLATE_SRC = ROOT / "20_new_runs/snp_screen_2026-08-04/G3_redocking/reference_template.json"
PDBFIXER_ENV = "mdenv"
WSL_DISTRO, CONDA_ENV = "Ubuntu-22.04", "diffdock"


def _wsl_path(p) -> str:
    s = str(p).replace(chr(92), "/")
    return f"/mnt/{s[0].lower()}{s[2:]}" if s[1:3] == ":/" else s


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


MR = _load("multi_reference_ifp_variant", SCRIPTS / "multi_reference_ifp_variant.py")
RSS = _load("run_snp_screen", SCRIPTS / "run_snp_screen.py")


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def wsl(cmd, timeout=7200):
    full = (f"source ~/miniconda3/etc/profile.d/conda.sh && conda activate {CONDA_ENV} "
            f"&& cd '/mnt/d/Agentic VS Workflow' && {cmd}")
    return subprocess.run(["wsl", "-d", WSL_DISTRO, "--", "bash", "-lc", full],
                          capture_output=True, text=True, timeout=timeout)


def write_receptor_with_waters(pdb_path: Path, out_pdb: Path, chain: str,
                               cutoff: float, ligand_resn: str) -> dict:
    """`run_snp_screen._write_receptor_with_waters`, transcribed, with the
    co-crystal ligand's residue name made a parameter.

    The production function hardcodes "ZMA" and raises when it finds none, so it
    cannot prepare 5OLH at all. Nothing else here differs from it: same primary
    altloc rule, same per-residue minimum distance, same chain filter, same
    output order. `--verify-transcription` checks that claim by running both on
    4EIY and comparing the bytes.
    """
    lines = pdb_path.read_text(encoding="utf-8", errors="replace").splitlines()

    def xyz(l): return (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    def usable(l): return l[16] in (" ", "A")
    def is_hydrogen(l):
        return (l[76:78].strip() or l[12:16].strip()[:1]).upper().startswith("H")

    ligand = [xyz(l) for l in lines
              if l.startswith(("ATOM", "HETATM")) and l[17:20].strip() == ligand_resn
              and usable(l) and not is_hydrogen(l)]
    if not ligand:
        raise RuntimeError(f"No {ligand_resn} heavy atoms in {pdb_path}")

    groups: dict = {}
    for l in lines:
        if not l.startswith("HETATM") or l[17:20].strip() != "HOH" or not usable(l):
            continue
        d = min(math.dist(xyz(l), z) for z in ligand)
        entry = groups.setdefault((l[21], l[22:27]), [d, []])
        entry[0] = min(entry[0], d)
        entry[1].append(l)

    kept_keys = [k for k, (d, _) in groups.items() if d <= cutoff]
    kept = [l for k in kept_keys for l in groups[k][1]]
    protein = [l for l in lines
               if l.startswith("ATOM") and l[21] == chain and usable(l)]
    out_pdb.write_text("\n".join(protein + kept) + "\nEND\n", encoding="utf-8")
    return {"n_waters_in_file": len(groups), "n_waters_retained": len(kept_keys),
            "n_ligand_heavy_atoms": len(ligand), "ligand_resn": ligand_resn}


def superpose_pdb(src: Path, dst: Path, R, t):
    """Whole file, every record, mapped by (R, t). Waters and ligand move with
    the protein."""
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
    """Production G2 sequence, in the environment production used (`mdenv`).
    A Windows conda environment with the same package versions kills the
    interpreter at `addMissingAtoms()`; see `crossdock_zma.pdbfix`."""
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
        raise RuntimeError(f"pdbfixer failed for {src.name}\n"
                           f"stdout: {r.stdout[-400:]}\nstderr: {r.stderr[-400:]}")
    return r.stdout.strip().splitlines()[-1]


def crystal_ligand(src: Path, dst: Path, resn: str, template_smiles: str):
    """The co-crystal ligand alone, primary altloc, bond orders from a template.
    A PDB HETATM block carries no bond orders and inferring them is guessing at
    chemistry."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    keep = [l for l in src.read_text(encoding="utf-8", errors="ignore").splitlines()
            if l.startswith("HETATM") and l[17:20].strip() == resn
            and l[16] in (" ", "A")]
    if not keep:
        raise RuntimeError(f"no {resn} HETATM records in {src}")
    dst.write_text("\n".join(keep) + "\nEND\n", encoding="utf-8", newline="\n")
    raw = Chem.MolFromPDBFile(str(dst), removeHs=True, sanitize=False)
    return AllChem.AssignBondOrdersFromTemplate(
        Chem.MolFromSmiles(template_smiles), raw)


def verify_transcription() -> int:
    """Run the production function and this one on 4EIY and compare the bytes."""
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        a, b = Path(d) / "prod.pdb", Path(d) / "mine.pdb"
        ra = RSS._write_receptor_with_waters(STRUCT / "4EIY.pdb", a,
                                             chain="A", cutoff=WATER_CUTOFF_A)
        rb = write_receptor_with_waters(STRUCT / "4EIY.pdb", b, chain="A",
                                        cutoff=WATER_CUTOFF_A, ligand_resn="ZMA")
        same = a.read_bytes() == b.read_bytes()
        print(f"production : {ra['n_waters_retained']} waters retained")
        print(f"transcribed: {rb['n_waters_retained']} waters retained")
        print("identical output:", same)
        return 0 if same else 1


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out")
    ap.add_argument("--verify-transcription", action="store_true")
    args = ap.parse_args()
    if args.verify_transcription:
        return verify_transcription()
    if not args.out:
        raise SystemExit("--out is required")
    out = Path(args.out).resolve()

    # Rule 1, enforced.
    if out.exists() and any(out.iterdir()):
        raise SystemExit(
            f"REFUSED: {out} already exists and is not empty.\n"
            "A run directory is written once (ADR-AIIM-092). Choose a new one.")
    (out / "inputs").mkdir(parents=True, exist_ok=True)
    (out / "receptors").mkdir(exist_ok=True)

    receptors = sorted({a[4] for a in ARMS} | {a[3] for a in ARMS})

    # Rule 2, enforced.
    manifest = {"declared_inputs": {}, "recipe": {
        "preparation": "production G2, transcribed: chain A + waters within "
                       f"{WATER_CUTOFF_A} A of that structure's OWN co-crystal "
                       "ligand -> pdbfixer (missing atoms, hydrogens at pH "
                       f"{PROTONATION_PH}, keepWater=True) -> obabel -xr",
        "docking": f"{N_CONFORMERS} conformers, exhaustiveness {EXHAUSTIVENESS}, "
                   f"{N_POSES} poses, seed {SEED}",
        "water_shell_ligand": OWN_LIGAND,
        "arms": [{"arm": a, "compound": c, "resn": r, "reference_from": s,
                  "docked_into": d} for a, c, r, s, d in ARMS],
    }}
    for src in ([BOX_SRC, TEMPLATE_SRC]
                + [STRUCT / f"{p}.pdb" for p in receptors]):
        shutil.copy2(src, out / "inputs" / src.name)
        manifest["declared_inputs"][src.name] = {
            "from": src.relative_to(ROOT).as_posix(), "sha256": sha256(src)}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2),
                                       encoding="utf-8")
    smiles_zma = json.loads(TEMPLATE_SRC.read_text(encoding="utf-8"))["smiles"]
    smiles = {"ZM241385": smiles_zma, "9XT": SMILES_9XT}

    # --- prepare, everything expressed in FRAME's coordinates ---------------
    ref_ca = MR.ca_by_resi(MR.read_pdb(STRUCT / f"{FRAME}.pdb"))
    prep = {}
    print("preparing receptors (production G2 recipe, own-ligand water shell)")
    for pdb in receptors:
        src = STRUCT / f"{pdb}.pdb"
        ca = MR.ca_by_resi(MR.read_pdb(src))
        shared = sorted((r for r in set(ca) & set(ref_ca)
                         if int(r) <= RECEPTOR_RESI_MAX), key=int)
        P = np.array([ca[r] for r in shared]); Q = np.array([ref_ca[r] for r in shared])
        R, t = (MR.kabsch(P, Q) if pdb != FRAME else (np.eye(3), np.zeros(3)))
        resid = float(np.sqrt((((P @ R.T + t) - Q) ** 2).sum(1).mean()))

        sup = out / "receptors" / f"{pdb}_super.pdb"
        superpose_pdb(src, sup, R, t)
        filt = out / "receptors" / f"{pdb}_protein_plus_waters.pdb"
        wat = write_receptor_with_waters(sup, filt, chain="A",
                                         cutoff=WATER_CUTOFF_A,
                                         ligand_resn=OWN_LIGAND[pdb])
        fixed = out / "receptors" / f"{pdb}_fixed.pdb"
        fixinfo = pdbfix(filt, fixed)
        pdbqt = out / "receptors" / f"{pdb}_prepared.pdbqt"
        wsl(f'obabel -ipdb "{fixed.relative_to(ROOT).as_posix()}" -opdbqt '
            f'-O "{pdbqt.relative_to(ROOT).as_posix()}" -xr')
        prep[pdb] = {"n_shared_CA": len(shared),
                     "superposition_rmsd_A": round(resid, 3),
                     "water_shell_around": OWN_LIGAND[pdb],
                     "waters_retained": wat["n_waters_retained"],
                     "pdbfixer": fixinfo,
                     "pdbqt_ok": pdbqt.is_file() and pdbqt.stat().st_size > 0}
        p = prep[pdb]
        print(f"  {pdb}: superposition {p['superposition_rmsd_A']:.3f} A, "
              f"waters {p['waters_retained']} around {p['water_shell_around']}, "
              f"{p['pdbfixer']}")
    (out / "prep.json").write_text(json.dumps(prep, indent=2), encoding="utf-8")

    # --- reference poses, all carried into FRAME's coordinates --------------
    refs = {}
    for _, comp, resn, sfrom, _ in ARMS:
        if comp in refs:
            continue
        lig = out / "receptors" / f"{comp}_crystal_from_{sfrom}.pdb"
        refs[comp] = crystal_ligand(out / "receptors" / f"{sfrom}_super.pdb",
                                    lig, resn, smiles[comp])
        print(f"  reference {comp}: {refs[comp].GetNumAtoms()} atoms "
              f"(from {sfrom}, in {FRAME}'s frame)")

    # --- dock ---------------------------------------------------------------
    box = out / "box.json"
    box.write_text(BOX_SRC.read_text(encoding="utf-8"), encoding="utf-8", newline="\n")
    print("\ndocking")
    for arm, comp, _, _, into in ARMS:
        ledger = out / f"ligand_{arm}.csv"
        ledger.write_text(f"compound_id,SMILES\n{comp},{smiles[comp]}\n",
                          encoding="utf-8", newline="\n")
        rec = (out / "receptors" / f"{into}_prepared.pdbqt").relative_to(ROOT).as_posix()
        print(f"  {arm}: {comp} into {into} ...", flush=True)
        wsl(f'python 10_new_system/scripts/g4_dock_library.py '
            f'--ledger "{ledger.relative_to(ROOT).as_posix()}" --receptor "{rec}" '
            f'--box "{box.relative_to(ROOT).as_posix()}" '
            f'--out "{(out / f"scores_{arm}.csv").relative_to(ROOT).as_posix()}" '
            f'--poses-dir "{(out / f"poses_{arm}").relative_to(ROOT).as_posix()}" '
            f'--id-col compound_id --smiles-col SMILES '
            f'--n-conformers {N_CONFORMERS} --exhaustiveness {EXHAUSTIVENESS} '
            f'--n-poses {N_POSES} --seed {SEED} --cpu 6')

    # --- score, with rule 3 enforced ----------------------------------------
    import warnings; warnings.filterwarnings("ignore")
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdMolAlign
    RDLogger.DisableLog("rdApp.*")
    res = {}
    for arm, comp, _, sfrom, into in ARMS:
        sdf = out / f"poses_{arm}" / f"{comp}.sdf"
        if not sdf.is_file():
            res[arm] = {"available": False, "reason": "no poses written"}; continue
        ms = [m for m in Chem.SDMolSupplier(str(sdf), removeHs=True) if m]
        r = [round(rdMolAlign.CalcRMS(m, refs[comp]), 3) for m in ms]
        res[arm] = {"available": True, "compound": comp, "docked_into": into,
                    "reference_from": sfrom, "n_poses": len(r),
                    "rank1_rmsd_A": r[0], "best_rmsd_A": min(r), "all_rmsd_A": r,
                    "within_2A": r[0] <= RMSD_PASS_A}
        print(f"  {arm} ({comp} -> {into}): rank-1 {r[0]} A "
              f"(best of {len(r)}: {min(r)})")

    ctrl = res.get("control", {})
    ctrl_ok = bool(ctrl.get("available")
                   and ctrl["rank1_rmsd_A"] <= CONTROL_MAX_RMSD_A)
    report = {
        "question": "Does redocking accuracy survive a change of LIGAND? "
                    "Every structure in ADR-AIIM-162 was solved with ZM241385.",
        "authority": "ADR-AIIM-168, criteria fixed before the run",
        "recipe": manifest["recipe"],
        "control": {"arm": "control", "compound": "ZM241385", "receptor": FRAME,
                    "max_rmsd_A": CONTROL_MAX_RMSD_A,
                    "observed_rmsd_A": ctrl.get("rank1_rmsd_A"), "passed": ctrl_ok},
        "by_arm": res,
        "declared_reading": {
            "T1_within_2A": "redocking accuracy transfers to a chemically "
                            "distinct ligand; the geometric licence extends "
                            "beyond ZM241385",
            "T1_beyond_2A": "it does not transfer; 0.318-0.709 A describes "
                            "ZM241385, not the protocol",
            "T2": "reported either way, gates nothing",
        },
    }
    if not ctrl_ok:
        report["outcome"] = "CONTROL_FAILED"
        report["results_withheld"] = (
            f"ZM241385 redocked into {FRAME} at {ctrl.get('rank1_rmsd_A')} A "
            f"against a ceiling of {CONTROL_MAX_RMSD_A} A. The preparation, not "
            "the ligand, is what this run measured. No number here supports any "
            "claim.")
        print(f"\nCONTROL FAILED ({ctrl.get('rank1_rmsd_A')} A) -- withheld")
    else:
        t1 = res.get("T1", {})
        report["outcome"] = ("TRANSFERS" if t1.get("within_2A") else
                             "DOES_NOT_TRANSFER" if t1.get("available") else
                             "T1_UNAVAILABLE")
        print(f"\ncontrol {ctrl['rank1_rmsd_A']} A -- passed.  "
              f"T1 (9XT into its own crystal) {t1.get('rank1_rmsd_A')} A "
              f"-> {report['outcome']}")
    (out / "redock_ligand_transfer.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")
    print(f"\nwritten -> {out / 'redock_ligand_transfer.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

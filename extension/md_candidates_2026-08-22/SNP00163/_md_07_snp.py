# ADR-AIIM-316: md_07_build_candidate.py with ONLY the two NP-library
# lookups redirected to SNP-Lib -- the pose directory and the identity
# ledger. Every other line, including the complex.pdb assembly whose
# protein/ligand/Na+/water ordering the original warns is silently
# breakable, is unchanged. 10_new_system is not modified.
"""
md_07_build_candidate.py -- build a membrane system for one Tier-A candidate,
starting from its docked pose.

This is the same protocol as the ZM241385 control (Section 2.9), differing in
exactly one input: the ligand and its coordinates. Keeping everything else
identical is the point -- the control exists so the candidate runs can be
compared against a known binding mode under the same treatment, and any
difference in receptor preparation, bilayer composition, force field or
equilibration schedule would spoil that comparison.

Two things are easier here than they were for the control, and one is harder.

Easier: the docked pose is an SDF, so bond orders are already present -- the
control had to recover them from a SMILES template because PDB does not store
them. And all three candidates are small and rigid (22-25 heavy atoms, two
rotatable bonds), so conformer and parameterisation trouble is unlikely.

Harder: the pose is *predicted*, not crystallographic. It is the rank-1 Vina
pose, and the whole purpose of the simulation is to ask whether it survives.
So the pose is taken verbatim, with no minimisation or manual adjustment that
would quietly improve the thing being tested.

The identity of the molecule is checked against the frozen prioritisation
ledger before anything is built: simulating a different molecule than the one
that was ranked would invalidate the result without failing anywhere.

Usage:
    python md_07_build_candidate.py NPLIB01185
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from pathlib import Path
from kdss_paths import project_root, resolve, under

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

LIGNAME = "LIG"          # residue name inside the built system
WSL_DISTRO = "A2A-MD"
AMBER = "/home/a2amd/A2A_MD/ambertools26_env"


def _root() -> Path:
    # Delegates to kdss_paths, which finds whichever layout this file is inside
    # -- the development tree or a deposit -- by walking up from its own
    # location. It previously tried two absolute paths and exited if neither
    # existed, which on the author's machine meant a copy of the deposit read
    # and wrote the development tree instead of itself.
    return project_root()


ROOT = _root()
FUNNEL = resolve("20_new_runs/snplib_dock_2026-10-06")   # ADR-AIIM-316
CONTROL = resolve("20_new_runs/md_production_2026-08-21")
DEST = resolve("20_new_runs/md_candidates_2026-08-22")


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def wsl(cmd: str) -> tuple[int, str]:
    """Run a bash command inside the AmberTools distro.

    Written to a file rather than passed inline: nested quoting through
    wsl.exe has silently mangled `$VAR` expansions in this project before.
    """
    tmp = under("20_new_runs/md_candidates_2026-08-22/SNP00163", "_cmd_snp00163.sh")
    tmp.parent.mkdir(parents=True, exist_ok=True)
    body = (f"#!/bin/bash\nexport AMBERHOME={AMBER}\n"
            f"export PATH=$AMBERHOME/bin:$PATH\n{cmd}\n")
    tmp.write_bytes(body.encode())            # binary write keeps LF
    wpath = "/mnt/d/" + str(tmp.relative_to(Path("D:/"))).replace("\\", "/")
    p = subprocess.run(["wsl.exe", "-d", WSL_DISTRO, "--", "bash", wpath],
                       capture_output=True, text=True,
                       env={**__import__("os").environ, "MSYS_NO_PATHCONV": "1"})
    return p.returncode, (p.stdout or "") + (p.stderr or "")


def find_pose(cid: str) -> Path:
    for d in sorted(FUNNEL.glob("poses*")):   # ADR-AIIM-316
        f = d / f"{cid}.sdf"
        if f.is_file():
            return f
    raise SystemExit(f"no docked pose found for {cid}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("compound_id")
    args = ap.parse_args()
    cid = args.compound_id
    run = DEST / cid
    build = run / "03_build"
    build.mkdir(parents=True, exist_ok=True)

    # ---- 1. the pose, and a check that it is the right molecule -----------
    sdf = find_pose(cid)
    mols = [m for m in Chem.SDMolSupplier(str(sdf), removeHs=False) if m is not None]
    if not mols:
        raise SystemExit(f"{sdf.name} parsed to zero molecules")
    pose = mols[0]                       # rank-1; the pose the ranking reports

    # ADR-AIIM-316: the SNP ledger in place of the NP prioritisation file.
    ledger = pd.read_csv(under("20_new_runs/snp_screen_2026-08-04", "G1_standardization/snp_source_ledger.csv")).set_index("SNP_ID")
    want = Chem.MolToInchiKey(Chem.MolFromSmiles(ledger.loc[cid, "smiles_canonical"]))
    got = Chem.MolToInchiKey(Chem.RemoveHs(pose))
    if want != got:
        raise SystemExit(
            f"{cid}: the saved pose is not the molecule the ledger ranked\n"
            f"  ledger InChIKey {want}\n  pose   InChIKey {got}")

    charge = Chem.GetFormalCharge(pose)
    ligh = Chem.AddHs(pose, addCoords=True)   # heavy atoms stay exactly where docked
    for a in ligh.GetAtoms():
        ri = a.GetPDBResidueInfo() or Chem.AtomPDBResidueInfo()
        ri.SetResidueName(LIGNAME)
        ri.SetResidueNumber(1)
        ri.SetIsHeteroAtom(True)
        a.SetMonomerInfo(ri)
    lig_pdb = build / "ligand.pdb"
    Chem.MolToPDBFile(ligh, str(lig_pdb))
    lig_sdf = build / "ligand_posed.sdf"
    Chem.SDWriter(str(lig_sdf)).write(ligh)

    print(f"{cid}")
    print(f"  pose        : {sdf.relative_to(ROOT)} (rank 1 of {len(mols)})")
    print(f"  identity    : InChIKey matches the ledger")
    print(f"  atoms       : {pose.GetNumAtoms()} heavy, {ligh.GetNumAtoms()} with H")
    print(f"  formal chg  : {charge:+d}")

    # ---- 2. GAFF2 + AM1-BCC ------------------------------------------------
    # SDF in, not PDB: antechamber segfaulted on PDB input during the control
    # build because it had to perceive bond orders itself. The SDF carries
    # them, so that failure mode is removed rather than worked around.
    rc, out = wsl(
        f'cd "{("/mnt/d/" + str(build.relative_to(Path("D:/"))).replace(chr(92), "/"))}" || exit 1\n'
        f'antechamber -i ligand_posed.sdf -fi sdf -o {LIGNAME}.mol2 -fo mol2 '
        f'-c bcc -nc {charge} -at gaff2 -rn {LIGNAME} -s 0 > antechamber.log 2>&1\n'
        f'echo "antechamber exit=$?"\n'
        f'parmchk2 -i {LIGNAME}.mol2 -f mol2 -o {LIGNAME}.frcmod -s gaff2 '
        f'> parmchk2.log 2>&1\n'
        f'echo "parmchk2 exit=$?"\n'
        f'grep -c ATTN {LIGNAME}.frcmod || true')
    print("  " + "\n  ".join(l for l in out.strip().splitlines() if l.strip())[:400])
    if not (build / f"{LIGNAME}.mol2").is_file():
        raise SystemExit(f"{cid}: antechamber produced no mol2 -- see antechamber.log")

    # parmchk2 flags every parameter it had to guess. GAFF2 has no penalty
    # score, so this count is the only coverage warning available.
    frc = (build / f"{LIGNAME}.frcmod").read_text(encoding="utf-8", errors="replace")
    attn = [l.strip() for l in frc.splitlines() if "ATTN" in l]
    print(f"  parmchk2    : {len(attn)} parameters flagged ATTN"
          + (f" -- first: {attn[0][:60]}" if attn else ""))

    # ---- 3. complex.pdb: shared receptor + this ligand ---------------------
    # The receptor, its caps, its 156 retained waters and the allosteric Na+
    # are byte-identical to the control's. Only the ligand differs.
    receptor_src = under("20_new_runs/md_production_2026-08-21", "03_build/complex.pdb")
    rec_lines = receptor_src.read_text(encoding="utf-8", errors="replace").splitlines()

    # ORDER MATTERS, and getting it wrong is silent. The control's layout is
    # protein / TER / ligand(chain B) / TER / Na+ / waters, and the ligand
    # there is not decoration -- it separates the protein from the
    # crystallographic Na+. Building the file with the ligand appended at the
    # end instead left the Na+ immediately after the last NME, and
    # packmol-memgen absorbed the ion INTO that residue: tleap then died on
    # "Atom .R<NME 304>.A<NA 7> does not have a type". Worse, it would have
    # been easy to "fix" by deleting the ion, which is the allosteric sodium
    # that 4EIY is known for. So the layout below reproduces the control's
    # exactly, and the ion's presence is asserted before the file is written.
    protein, hetero = [], []
    for l in rec_lines:
        if l.startswith("ATOM") or l.startswith("TER"):
            protein.append(l)
        elif l.startswith("HETATM") and l[17:20].strip() in ("NA", "HOH"):
            hetero.append(l)
    while protein and protein[-1].startswith("TER"):
        protein.pop()

    n_na = sum(1 for l in hetero if l[17:20].strip() == "NA")
    n_wat = sum(1 for l in hetero if l[17:20].strip() == "HOH")
    if n_na != 1:
        raise SystemExit(f"{cid}: expected exactly one crystallographic Na+, found {n_na}")

    lig_lines = [l for l in lig_pdb.read_text(encoding="utf-8").splitlines()
                 if l.startswith(("ATOM", "HETATM"))]
    # Chain B, residue 1 -- the chain the control's ligand occupies, so
    # packmol-memgen treats it as a separate molecule rather than folding it
    # into the protein.
    # PDB columns, 0-indexed: 0-5 record, 6-10 serial, 11 space, 12-15 atom
    # name, 16 altLoc, 17-19 resName, 20 space, 21 chainID, 22-25 resSeq.
    # Built by slice rather than by counting spaces in a format string -- the
    # first attempt was one character short and silently produced residue "IG".
    lig_lines = ["HETATM" + l[6:16] + " " + f"{LIGNAME:>3}" + " " + "B"
                 + f"{1:>4}" + l[26:] for l in lig_lines]

    out_lines, serial = [], 1
    def emit(lines):
        nonlocal serial
        for l in lines:
            if l.startswith("TER"):
                out_lines.append("TER")
                continue
            out_lines.append(l[:6] + f"{serial:>5}" + l[11:])
            serial += 1

    emit(protein)
    out_lines.append("TER")
    emit(lig_lines)
    out_lines.append("TER")
    emit(hetero)
    out_lines.append("END")
    complex_pdb = build / "complex.pdb"
    complex_pdb.write_text("\n".join(out_lines) + "\n", encoding="utf-8")

    n_prot = sum(1 for l in protein if not l.startswith("TER"))
    print(f"  complex.pdb : {n_prot} protein + {len(lig_lines)} ligand + "
          f"{n_na} Na+ + {n_wat} waters = {serial - 1} atoms")
    print(f"                order: protein / ligand(chain B) / Na+ / waters "
          f"-- matches the control")

    (run / "build_metadata.json").write_text(json.dumps({
        "compound_id": cid,
        "built": "2026-10-08",
        "pose_source": str(sdf.relative_to(ROOT)),
        "pose_sha256": sha256(sdf),
        "pose_rank": 1,
        "n_poses_available": len(mols),
        "pose_taken_verbatim": "no minimisation or manual adjustment; the "
                               "predicted pose is what the simulation tests",
        "inchikey": got,
        "identity_checked_against": "snp_source_ledger.csv (G1_standardization)",
        "formal_charge": charge,
        "parmchk2_attn_count": len(attn),
        "parmchk2_attn_lines": attn,
        "receptor_source": str(receptor_src.relative_to(ROOT)),
        "receptor_sha256": sha256(receptor_src),
        "protocol": "identical to the ZM241385 control (Section 2.9); only the "
                    "ligand differs",
    }, indent=2), encoding="utf-8")
    print(f"  metadata    -> {(run / 'build_metadata.json').relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""
md_02_prepare_ligand.py — protonate a crystallographic ligand and generate
GAFF2 parameters, entirely from the command line.

The PDB ligand carries heavy-atom coordinates only, and bond orders are not
stored in PDB format. Both are recovered before parameterisation:

  1. Bond orders are assigned from a SMILES template taken from this
     project's own frozen funnel metadata (`funnel_metadata.json`,
     `reference_smiles[0]`) rather than from an external lookup, so the
     template is hash-anchored to an artifact already in the run record.
  2. Hydrogens are added with coordinates, keeping heavy atoms fixed at their
     crystallographic positions.
  3. antechamber assigns AM1-BCC charges and GAFF2 atom types; parmchk2
     reports any parameter it had to estimate.

GAFF2 has no equivalent of CGenFF's penalty score, so the parmchk2 output is
the only warning channel for poor parameter coverage and is therefore parsed
and reported rather than written and ignored.

Run inside the AmberTools environment (needs rdkit + antechamber + parmchk2).

Usage:
    python md_02_prepare_ligand.py --pdb <ligand.pdb> --smiles <SMILES> \
        --name ZMA --charge 0 --out <dir>
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(cmd: list[str], cwd: Path) -> tuple[int, str]:
    p = subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True)
    return p.returncode, (p.stdout or "") + (p.stderr or "")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--smiles", required=True)
    ap.add_argument("--name", default="LIG")
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--out", required=True)
    ap.add_argument("--amberhome", default="/home/a2amd/A2A_MD/ambertools26_env")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    amber = Path(args.amberhome)

    raw = Chem.MolFromPDBFile(args.pdb, removeHs=False, sanitize=False)
    if raw is None:
        print(f"ERROR: could not read {args.pdb}", file=sys.stderr)
        return 2

    template = Chem.MolFromSmiles(args.smiles)
    if template is None:
        print("ERROR: template SMILES did not parse", file=sys.stderr)
        return 2

    try:
        mol = AllChem.AssignBondOrdersFromTemplate(template, raw)
    except Exception as exc:
        print(f"ERROR: bond-order assignment failed: {exc}", file=sys.stderr)
        print("  heavy atoms in PDB :", raw.GetNumAtoms())
        print("  heavy atoms in tmpl:", template.GetNumAtoms())
        return 2

    Chem.SanitizeMol(mol)
    molh = Chem.AddHs(mol, addCoords=True)
    n_h = molh.GetNumAtoms() - mol.GetNumAtoms()

    formal = Chem.GetFormalCharge(molh)
    if formal != args.charge:
        print(f"WARNING: RDKit formal charge {formal} != --charge {args.charge}; "
              f"using {args.charge} for antechamber as instructed")

    lig_h = out / f"{args.name}_H.pdb"
    Chem.MolToPDBFile(molh, str(lig_h), flavor=4)

    # Feed antechamber an SDF, not a PDB. PDB carries no bond orders, so
    # antechamber must perceive them itself via `bondtype`, which segfaults on
    # this ligand. RDKit already holds the correct bond orders (assigned from
    # the frozen SMILES template above), so writing SDF hands them over
    # directly and removes the perception step entirely.
    lig_sdf = out / f"{args.name}_H.sdf"
    w = Chem.SDWriter(str(lig_sdf))
    molh.SetProp("_Name", args.name)
    w.write(molh)
    w.close()

    ac_cmd = [
        str(amber / "bin/antechamber"),
        "-i", lig_sdf.name, "-fi", "sdf",
        "-o", f"{args.name}.mol2", "-fo", "mol2",
        "-c", "bcc", "-nc", str(args.charge),
        "-at", "gaff2", "-rn", args.name, "-s", "2",
    ]
    print("running antechamber (AM1-BCC; this is the slow step) ...")
    rc, log = run(ac_cmd, out)
    (out / "antechamber.log").write_text(log, encoding="utf-8")
    if rc != 0 or not (out / f"{args.name}.mol2").exists():
        print(f"ERROR: antechamber failed (rc={rc}); see antechamber.log", file=sys.stderr)
        print(log[-2000:])
        return 3

    rc, log = run([str(amber / "bin/parmchk2"),
                   "-i", f"{args.name}.mol2", "-f", "mol2",
                   "-o", f"{args.name}.frcmod", "-s", "gaff2"], out)
    (out / "parmchk2.log").write_text(log, encoding="utf-8")
    frcmod = out / f"{args.name}.frcmod"
    if rc != 0 or not frcmod.exists():
        print(f"ERROR: parmchk2 failed (rc={rc})", file=sys.stderr)
        return 4

    # GAFF2 gives no penalty score, so estimated parameters are the only
    # available signal that the ligand is poorly covered. Surface them.
    text = frcmod.read_text(encoding="utf-8")
    flagged = [ln.strip() for ln in text.splitlines()
               if "ATTN" in ln.upper() or "penalty" in ln.lower()]

    meta = {
        "ligand_name": args.name,
        "source_pdb": args.pdb,
        "source_pdb_sha256": sha256(Path(args.pdb)),
        "template_smiles": args.smiles,
        "template_source": "funnel_metadata.json reference_smiles[0] (frozen artifact)",
        "heavy_atoms": mol.GetNumAtoms(),
        "hydrogens_added": n_h,
        "rdkit_formal_charge": formal,
        "antechamber_charge_used": args.charge,
        "atom_type_set": "gaff2",
        "charge_method": "AM1-BCC",
        "parmchk2_flagged_lines": flagged,
        "parmchk2_needs_attention": bool(flagged),
        "antechamber_input_format": "sdf (bond orders supplied by RDKit; "
                                    "PDB input made antechamber's bondtype segfault)",
        "outputs": {
            "ligand_with_H": {"path": str(lig_h), "sha256": sha256(lig_h)},
            "ligand_sdf": {"path": str(lig_sdf), "sha256": sha256(lig_sdf)},
            "mol2": {"path": str(out / f'{args.name}.mol2'),
                     "sha256": sha256(out / f"{args.name}.mol2")},
            "frcmod": {"path": str(frcmod), "sha256": sha256(frcmod)},
        },
    }
    (out / f"{args.name}_ligand_prep_metadata.json").write_text(
        json.dumps(meta, indent=2), encoding="utf-8")

    print(f"heavy atoms      : {mol.GetNumAtoms()}")
    print(f"hydrogens added  : {n_h}")
    print(f"formal charge    : {formal} (antechamber used {args.charge})")
    print(f"parmchk2 ATTN    : {len(flagged)} line(s)")
    for ln in flagged[:10]:
        print(f"   {ln}")
    print(f"\nwritten -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

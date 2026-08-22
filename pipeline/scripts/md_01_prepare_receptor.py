"""
md_01_prepare_receptor.py — scripted preparation of 4EIY for membrane MD.

Implements the decisions recorded in
`00_project_control/25_MD_PROTOCOL_DESIGN.md` §1.1, all of which were read
off the frozen PDB rather than assumed:

  * BRIL (apocytochrome b562RIL) occupies residues 1001-1106 and is deleted;
    the receptor is a fusion construct and simulating the fusion would be
    simulating a crystallisation scaffold.
  * The three engineered mutations in the file all map to BRIL (UniProt
    P0ABE7), not to A2AR (P29274), so the receptor sequence is wild-type and
    nothing is reverted.
  * The allosteric Na+ at Asp2.50 is retained (Liu et al., Science 2012).
  * LCP host lipids (OLA/OLC/OLB) and cryoprotectant (PEG) are removed.
  * Crystallographic cholesterol (CLR) is removed here because cholesterol is
    placed as a bilayer component by packmol-memgen; this is a deliberate
    choice, recorded rather than implicit.
  * ZM241385 (ZMA) is written to its own file for ligand parameterisation.
  * Crystallographic waters are written separately so the retain/discard
    decision stays explicit downstream.

No coordinates are invented. The ICL3 gap left by BRIL removal (residues
209-218) is NOT modelled here; capping or loop modelling is a separate,
recorded decision.

Usage:
    python md_01_prepare_receptor.py --out <dir>
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

ROOT = Path(r"D:\Agentic VS Workflow")
PDB_IN = ROOT / "02_frozen_inputs/a2ar_structures/4EIY.pdb"

BRIL_RANGE = (1001, 1106)          # read from SEQADV/ATOM blocks, see §1.1
DROP_HET = {"OLA", "OLC", "OLB", "PEG", "CLR"}
LIGAND_RES = "ZMA"
KEEP_ION = "NA"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "20_new_runs/md_benchmark_2026-08-21/01_receptor"))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    protein: list[str] = []
    ligand: list[str] = []
    waters: list[str] = []
    ion: list[str] = []
    dropped: Counter[str] = Counter()
    altloc_kept = 0
    altloc_dropped = 0

    for line in PDB_IN.read_text(encoding="utf-8", errors="replace").splitlines():
        rec = line[:6].strip()
        if rec not in ("ATOM", "HETATM"):
            continue

        resname = line[17:20].strip()
        resseq = int(line[22:26])
        altloc = line[16]

        # Keep only the primary alternate conformation; a system cannot hold
        # two conformers of the same atom.
        if altloc not in (" ", "A"):
            altloc_dropped += 1
            continue
        if altloc == "A":
            altloc_kept += 1
            line = line[:16] + " " + line[17:]

        if rec == "ATOM":
            if BRIL_RANGE[0] <= resseq <= BRIL_RANGE[1]:
                dropped["BRIL_protein_atoms"] += 1
            else:
                protein.append(line)
            continue

        # HETATM
        if resname in DROP_HET:
            dropped[resname] += 1
        elif resname == LIGAND_RES:
            ligand.append(line)
        elif resname == "HOH":
            waters.append(line)
        elif resname == KEEP_ION:
            ion.append(line)
        else:
            dropped[f"unexpected:{resname}"] += 1

    def write(name: str, lines: list[str]) -> Path:
        p = out / name
        p.write_text("\n".join(lines) + "\nEND\n", encoding="utf-8")
        return p

    p_prot = write("4eiy_A2AR_noBRIL.pdb", protein + ion)
    p_lig = write("ZMA_ligand.pdb", ligand)
    p_wat = write("crystal_waters.pdb", waters)

    prot_res = sorted({int(l[22:26]) for l in protein})
    blocks: list[tuple[int, int]] = []
    s = p = prot_res[0]
    for n in prot_res[1:]:
        if n != p + 1:
            blocks.append((s, p))
            s = n
        p = n
    blocks.append((s, p))

    meta = {
        "source_pdb": str(PDB_IN),
        "source_sha256": sha256(PDB_IN),
        "bril_range_deleted": list(BRIL_RANGE),
        "het_dropped": dict(dropped),
        "altloc_A_promoted": altloc_kept,
        "altloc_B_or_later_dropped": altloc_dropped,
        "protein_atoms_kept": len(protein),
        "protein_residue_blocks": blocks,
        "icl3_gap": "209-218 (left unmodelled; cap-vs-model is a recorded open decision)",
        "na_ion_atoms": len(ion),
        "ligand_atoms": len(ligand),
        "crystal_water_molecules": len(waters),
        "outputs": {
            "protein_plus_ion": {"path": str(p_prot), "sha256": sha256(p_prot)},
            "ligand": {"path": str(p_lig), "sha256": sha256(p_lig)},
            "crystal_waters": {"path": str(p_wat), "sha256": sha256(p_wat)},
        },
        "note": "Cholesterol (CLR) removed here by design: cholesterol enters "
                "as a bilayer component via packmol-memgen, not as a "
                "crystallographic copy.",
    }
    (out / "receptor_prep_metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

    print(f"protein atoms kept : {len(protein)}")
    print(f"residue blocks     : {blocks}")
    print(f"Na+ atoms          : {len(ion)}")
    print(f"ligand (ZMA) atoms : {len(ligand)}")
    print(f"crystal waters     : {len(waters)}")
    print(f"altLoc A promoted  : {altloc_kept}   (B+ dropped: {altloc_dropped})")
    print(f"dropped            : {dict(dropped)}")
    print(f"\nwritten -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

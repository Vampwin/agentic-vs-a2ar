"""Inputs of ADR-AIIM-353, built from geometry only (no docking result is read).

  inputs/vina_blind_box.json       centre and size of the whole-receptor box (arm A)
  inputs/vina_blind_ledger.csv     the 15 survivors + ZM241385, SMILES as in the existing ledgers
  inputs/4EIY_pocket12A.pdb        4EIY_protein_only.pdb cropped to residues within 12.0 A of the ZM241385 crystal pose (arm B)
  inputs/diffdock_pocket_input.csv DiffDock input for arm B (protein path is made absolute by the run script)

Anchors asserted here are the pre-registered values: 101 residues, box centre/size.
"""
import csv
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem

ROOT = Path(r"D:\Agentic VS Workflow")
RUNS = ROOT / "20_new_runs"
OUT = RUNS / "docking_matched_search_space_2026-10-09"
INP = OUT / "inputs"
PDB = RUNS / "diffdock_snplib_2026-10-07" / "4EIY_protein_only.pdb"
ZMA = RUNS / "snp_screen_2026-08-04" / "G2_preparation" / "ZMA_reference_4EIY.sdf"
LEDGER = RUNS / "snplib_dock_2026-10-06" / "dock_input_ledger.csv"
DDIN = RUNS / "diffdock_snplib_2026-10-07" / "diffdock_input.csv"
ZMA_SMILES = "Nc1nc(NCCc2ccc(O)cc2)nc2nc(-c3ccco3)nn12"
BUFFER, CROP_R = 4.0, 12.0


def main():
    INP.mkdir(exist_ok=True)
    lines = [l for l in PDB.read_text().splitlines() if l.startswith("ATOM")]
    heavy = [(l, np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])) for l in lines if l[76:78].strip() != "H"]
    X = np.array([p for _, p in heavy])
    lo, hi = X.min(0), X.max(0)
    centre = ((hi + lo) / 2).round(3)
    size = [math.ceil((e + 2 * BUFFER) * 2) / 2 for e in (hi - lo)]
    assert np.allclose(centre, [8.647, -22.095, 15.858], atol=0.01) or np.allclose(centre, [8.647, -22.095, 15.86], atol=0.01), centre
    assert size == [69.5, 110.5, 51.0], size
    (INP / "vina_blind_box.json").write_text(json.dumps({
        "center": centre.tolist(), "size": size,
        "centers_by_receptor": {"4EIY": centre.tolist()},
        "derived_from": "4EIY_protein_only.pdb heavy-atom bounding box + 4 A on each side, rounded up to 0.5 A (ADR-AIIM-353)"}, indent=1))

    z = Chem.SDMolSupplier(str(ZMA), removeHs=True)[0].GetConformer().GetPositions()
    keep = set()
    for l, p in heavy:
        if np.min(np.linalg.norm(z - p, axis=1)) <= CROP_R:
            keep.add((l[21], l[22:27]))
    assert len(keep) == 101, len(keep)
    cropped = [l for l in lines if (l[21], l[22:27]) in keep]
    (INP / "4EIY_pocket12A.pdb").write_text("\n".join(cropped) + "\nEND\n")

    dd = pd.read_csv(DDIN)
    led = pd.read_csv(LEDGER)
    m = dd.merge(led, left_on="complex_name", right_on="SNP_ID")
    assert len(m) == 15 and (m.ligand_description == m.smiles_canonical).all()
    rows = [("ZM241385", ZMA_SMILES)] + list(zip(m.SNP_ID, m.smiles_canonical))
    with open(INP / "vina_blind_ledger.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["SNP_ID", "smiles_canonical"]); w.writerows(rows)
    with open(INP / "diffdock_pocket_input.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["complex_name", "protein_path", "ligand_description", "protein_sequence"])
        for n, s in rows:
            w.writerow([n, "4EIY_pocket12A.pdb", s, ""])
    print("box", centre.tolist(), size, "| residues kept", len(keep), "| atoms kept", len(cropped), "| compounds", len(rows))


if __name__ == "__main__":
    main()

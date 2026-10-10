"""Inputs of ADR-AIIM-355, from geometry only (no docking result is read).

  inputs/diffdock_l_input.csv        DiffDock-L blind, 16 compounds, protein path made absolute by the run script
  inputs/config_100samples.yaml      copy of DiffDock's default_inference_args.yaml, samples_per_complex 100 (only change)
  inputs/pocket_C.csv                DiffDock-Pocket: centre = ZM241385 heavy-atom centroid, flexible = 6 residues <= 3.5 A
  inputs/pocket_Cca.csv              same, centre = mean C-alpha of residues with C-alpha <= 5 A of the crystal ligand
  (C-rigid uses pocket_C.csv with --rigid)
Anchors asserted: centre, flexible set, C-alpha centre.
"""
import csv
import re
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = Path(r"D:\Agentic VS Workflow")
RUNS = ROOT / "20_new_runs"
OUT = RUNS / "docking_matched_filter_pocket_2026-10-09"
INP = OUT / "inputs"
PDB = RUNS / "diffdock_snplib_2026-10-07" / "4EIY_protein_only.pdb"
ZMA = RUNS / "snp_screen_2026-08-04" / "G2_preparation" / "ZMA_reference_4EIY.sdf"
LEDGER = RUNS / "docking_matched_search_space_2026-10-09" / "inputs" / "vina_blind_ledger.csv"
YAML_SRC = Path(r"\\wsl.localhost\Ubuntu-22.04\home\win_user\diffdock\default_inference_args.yaml")
FLEX_EXPECT = ["A:168", "A:169", "A:177", "A:250", "A:253", "A:270"]


def main():
    INP.mkdir(exist_ok=True)
    lines = [l for l in PDB.read_text().splitlines() if l.startswith("ATOM")]
    z = Chem.SDMolSupplier(str(ZMA), removeHs=True)[0].GetConformer().GetPositions()
    centroid = z.mean(0)
    res = {}
    for l in lines:
        k = (l[21], int(l[22:26]))
        res.setdefault(k, []).append((l[12:16].strip(), np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]), l[76:78].strip()))
    ca_near, flex = [], []
    for k, a in res.items():
        ca = [p for n, p, e in a if n == "CA"]
        if ca and np.min(np.linalg.norm(z - ca[0], axis=1)) <= 5.0:
            ca_near.append(ca[0])
        hv = [p for n, p, e in a if e != "H"]
        if hv and min(np.min(np.linalg.norm(z - p, axis=1)) for p in hv) <= 3.5:
            flex.append(k)
    flex_s = [f"{c}:{r}" for c, r in sorted(flex, key=lambda x: x[1])]
    ca_c = np.mean(ca_near, axis=0)
    assert np.allclose(centroid, [-0.424, 8.527, 17.127], atol=0.002), centroid
    assert flex_s == FLEX_EXPECT, flex_s
    assert np.allclose(ca_c, [-4.233, 10.800, 20.280], atol=0.002), ca_c

    led = pd.read_csv(LEDGER)                      # ZM241385 first, then the 15 survivors
    assert len(led) == 16 and led.SNP_ID.iloc[0] == "ZM241385"
    with open(INP / "diffdock_l_input.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["complex_name", "protein_path", "ligand_description", "protein_sequence"])
        for n, s in zip(led.SNP_ID, led.smiles_canonical):
            w.writerow([n, "4EIY_protein_only.pdb", s, ""])
    # DiffDock-Pocket's CSV loader accepts ligand FILES only (a SMILES in the ligand column is skipped), so each SMILES
    # is written as an SDF with an RDKit ETKDG conformer (seed 42). The model discards these coordinates
    # (RemoveAllConformers, then its own conformer generation) because --keep_local_structures is not used.
    (INP / "ligands").mkdir(exist_ok=True)
    lig_paths = {}
    for n, s in zip(led.SNP_ID, led.smiles_canonical):
        m = Chem.AddHs(Chem.MolFromSmiles(s)); assert AllChem.EmbedMolecule(m, randomSeed=42) == 0
        m.SetProp("_Name", n)
        with Chem.SDWriter(str(INP / "ligands" / f"{n}.sdf")) as w:
            w.write(m)
        lig_paths[n] = f"/mnt/d/Agentic VS Workflow/20_new_runs/docking_matched_filter_pocket_2026-10-09/inputs/ligands/{n}.sdf"
    fx = "-".join(FLEX_EXPECT)
    for name, c in (("pocket_C.csv", centroid), ("pocket_Cca.csv", ca_c)):
        with open(INP / name, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["complex_name", "experimental_protein", "ligand", "pocket_center_x", "pocket_center_y", "pocket_center_z", "flexible_sidechains"])
            for n, s in zip(led.SNP_ID, led.smiles_canonical):
                w.writerow([n, "4EIY_protein_only.pdb", lig_paths[n], round(float(c[0]), 3), round(float(c[1]), 3), round(float(c[2]), 3), fx])
    y = YAML_SRC.read_text()
    assert re.search(r"^samples_per_complex: 10$", y, re.M)
    (INP / "config_100samples.yaml").write_text(re.sub(r"^samples_per_complex: 10$", "samples_per_complex: 100", y, flags=re.M))
    print("centroid", centroid.round(3), "| flexible", fx, "| C-alpha centre", ca_c.round(3), "| compounds", len(led))


if __name__ == "__main__":
    main()

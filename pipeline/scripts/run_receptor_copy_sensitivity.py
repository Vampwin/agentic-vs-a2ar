"""Does the choice of receptor copy change Vina scores?

All seven copies of 4EIY_prepared.pdbqt in this project differ, by 679
polar-hydrogen coordinates, 72 of which lie inside the 22 A docking box with
the nearest 3.50 A from its centre. The manuscript's 18.05 A statement covers
only the three histidine atom-TYPE differences, not these.

Vina uses a polar-hydrogen model, so hydrogen placement in the site can in
principle change a score. Measured here on four ligands across three receptor
copies, at the production protocol.
"""
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation
from vina import Vina

CEN = [-0.424, 8.527, 17.127]
SIZE = [22, 22, 22]
RECS = {
    "snp_screen_2026-08-04": "20_new_runs/snp_screen_2026-08-04/G2_preparation/4EIY_prepared.pdbqt",
    "repro_run1": "20_new_runs/snp_screen_repro_run1_2026-08-15/G2_preparation/4EIY_prepared.pdbqt",
    "repro_run2": "20_new_runs/snp_screen_repro_run2_2026-08-15/G2_preparation/4EIY_prepared.pdbqt",
}

d = pd.read_csv("02_frozen_inputs/public_benchmarks/dude_a2ar/dude_a2ar_docking_set.csv")
tests = [("ZM241385 co-crystal", "Nc1nc(nc2c1ncn2CCc1ccc(O)cc1)-c1ccco1")]
act = d[d.label == 1]
for i in range(3):
    r = act.iloc[i]
    tests.append(("DUD-E active " + str(r.compound_id), r.SMILES))

hdr = "ligand".ljust(26)
for k in RECS:
    hdr += k.rjust(22)
hdr += "spread".rjust(11)
print(hdr)
print("-" * len(hdr))

worst = 0.0
for name, smi in tests:
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    AllChem.EmbedMolecule(m, randomSeed=20260729)
    AllChem.MMFFOptimizeMolecule(m)
    prep = MoleculePreparation()
    prep.prepare(m)
    lig = prep.write_pdbqt_string()
    row = []
    for rec in RECS.values():
        v = Vina(sf_name="vina", cpu=20, seed=20260729, verbosity=0)
        v.set_receptor(rec)
        v.set_ligand_from_string(lig)
        v.compute_vina_maps(center=CEN, box_size=SIZE)
        v.dock(exhaustiveness=32, n_poses=1)
        row.append(v.energies(n_poses=1)[0][0])
    spread = max(row) - min(row)
    worst = max(worst, spread)
    line = name[:26].ljust(26)
    for e in row:
        line += ("%.4f" % e).rjust(22)
    line += ("%.4f" % spread).rjust(11)
    print(line)

print()
print("largest spread across receptor copies: %.4f kcal/mol" % worst)
print("for scale, the reported score range spans several kcal/mol and the")
print("manuscript reports scores to three decimals")

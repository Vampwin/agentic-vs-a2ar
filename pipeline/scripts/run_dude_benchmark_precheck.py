"""run_dude_benchmark_precheck.py -- is DUD-E A2AR fit for this purpose, before we spend 20 hours?

The benchmarking literature documents DUD, DUD-E and MUV as overestimating
virtual-screening accuracy through obvious and hidden chemical biases, which is
the stated motivation for LIT-PCBA. Taking that at face value would either stop
the run or discount its result in advance.

Two of those biases are measurable from the SMILES alone, with nothing docked,
and the point of this script is to measure them rather than inherit the
reputation:

  1. **Trivial separability.** If a descriptor requiring no structure -- heavy
     atoms, molecular weight, cLogP, TPSA, rotatable bonds -- already separates
     actives from decoys, then enrichment measured on this set is partly or
     wholly that descriptor, and the benchmark cannot answer whether docking
     adds anything. This is the failure mode already documented for the
     project's own ChEMBL benchmark, where heavy-atom count alone reaches 0.584.

  2. **Analogue bias.** If the actives are analogue series from a small number
     of scaffolds, a method that recognises the series wins without modelling
     the binding mode. Measured as Murcko scaffold diversity.

A third, scaffold overlap between actives and decoys, is reported because its
absence is the hazard that discredited DUD-E for *trained* models: a classifier
can learn the decoy distribution rather than protein-ligand interaction. It
does not transfer to a physics-based scoring function, which has no training
set, and the distinction is recorded here rather than assumed.

Output -> 20_new_runs/dude_precheck_2026-09-22/
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kdss_paths import under  # noqa: E402

SET = under("02_frozen_inputs/public_benchmarks/dude_a2ar",
            "dude_a2ar_docking_set.csv")
# What the project's own ChEMBL benchmark scores on the same test, for contrast.
CHEMBL_HEAVY_ATOM_AUC = 0.584


def auc(y, s) -> float:
    y = np.asarray(y)
    s = np.asarray(s, float)
    order = np.argsort(s, kind="mergesort")
    r = np.empty(len(s), float)
    r[order] = np.arange(1, len(s) + 1)
    r = pd.DataFrame({"s": s, "r": r}).groupby("s")["r"].transform("mean").to_numpy()
    npos = int(y.sum())
    nneg = len(y) - npos
    return float((r[y == 1].sum() - npos * (npos + 1) / 2) / (npos * nneg))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(
        under("20_new_runs", "dude_precheck_2026-09-22")))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    raw = pd.read_csv(SET)
    rows = []
    for _, r in raw.iterrows():
        m = Chem.MolFromSmiles(r.SMILES)
        if m is None:
            continue
        try:
            scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=m)
        except Exception:                                   # noqa: BLE001
            scaf = ""
        rows.append({
            "label": int(r.label), "scaffold": scaf,
            "heavy_atoms": m.GetNumHeavyAtoms(),
            "molecular_weight": Descriptors.MolWt(m),
            "clogp": Crippen.MolLogP(m),
            "tpsa": rdMolDescriptors.CalcTPSA(m),
            "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(m),
        })
    d = pd.DataFrame(rows)
    d.to_csv(out / "dude_precheck_descriptors.csv", index=False)

    trivial = {c: round(auc(d.label, d[c]), 4) for c in
               ("heavy_atoms", "molecular_weight", "clogp", "tpsa",
                "rotatable_bonds")}
    worst = max(trivial, key=lambda k: abs(trivial[k] - 0.5))

    act = d[d.label == 1]
    dec = d[d.label == 0]
    series = act.groupby("scaffold").size().sort_values(ascending=False)

    summary = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "purpose": "decide, before docking, whether DUD-E A2AR can answer "
                   "whether docking adds information over trivial descriptors",
        "nothing_docked": True,
        "set": "02_frozen_inputs/public_benchmarks/dude_a2ar/dude_a2ar_docking_set.csv",
        "n_parsed": int(len(d)), "n_actives": int(len(act)),
        "n_decoys": int(len(dec)),
        "trivial_descriptor_auc": trivial,
        "worst_trivial_descriptor": worst,
        "worst_trivial_auc": trivial[worst],
        "chembl_benchmark_heavy_atom_auc_for_contrast": CHEMBL_HEAVY_ATOM_AUC,
        "analogue_bias": {
            "unique_scaffolds_actives": int(act.scaffold.nunique()),
            "n_actives": int(len(act)),
            "median_compounds_per_scaffold": float(series.median()),
            "largest_series": int(series.iloc[0]),
            "largest_series_fraction": round(float(series.iloc[0] / len(act)), 4),
            "top5_scaffold_fraction": round(float(series.head(5).sum() / len(act)), 4),
        },
        "scaffold_overlap_actives_decoys": int(
            len(set(act.scaffold) & set(dec.scaffold))),
        "unique_scaffolds_decoys": int(dec.scaffold.nunique()),
        "reading": {
            "trivial_separability": "none of practical concern"
                if abs(trivial[worst] - 0.5) < 0.12 else "present, and material",
            "analogue_bias": "none"
                if series.iloc[0] / len(act) < 0.05 else "present",
            "zero_scaffold_overlap": "by construction; disqualifying for a "
                                     "TRAINED classifier, which can learn the "
                                     "decoy distribution, but not for a "
                                     "physics-based scoring function with no "
                                     "training set",
            "presumed_inactives": "decoys are presumed, not measured, inactive. "
                                  "A truly active decoy scored well counts as a "
                                  "false positive, so this biases measured "
                                  "enrichment DOWNWARD and is conservative",
        },
    }
    with open(out / "dude_precheck.json", "w", encoding="utf-8",
              newline="\n") as f:
        json.dump(summary, f, indent=1, ensure_ascii=False)
        f.write("\n")

    print("=" * 76)
    print("IS DUD-E A2AR FIT FOR THIS PURPOSE?  (nothing docked)")
    print("=" * 76)
    print(f"  {len(act)} actives / {len(dec)} decoys, {len(d)} parsed\n")
    print("  trivial descriptors, no structure required:")
    for k, v in trivial.items():
        flag = "  <-- worst" if k == worst else ""
        print(f"     {k:18} AUC {v:.4f}{flag}")
    print(f"     {'(ChEMBL benchmark)':18} AUC {CHEMBL_HEAVY_ATOM_AUC:.4f}"
          f"  heavy atoms, for contrast")
    print()
    ab = summary["analogue_bias"]
    print("  analogue bias:")
    print(f"     unique scaffolds among actives : {ab['unique_scaffolds_actives']}"
          f" of {ab['n_actives']}")
    print(f"     largest series                 : {ab['largest_series']} "
          f"({ab['largest_series_fraction']:.1%} of actives)")
    print(f"     top 5 scaffolds cover          : {ab['top5_scaffold_fraction']:.1%}")
    print()
    print(f"  scaffold overlap actives/decoys  : "
          f"{summary['scaffold_overlap_actives_decoys']}")
    print()
    print(f"  VERDICT  trivial separability: {summary['reading']['trivial_separability']}")
    print(f"           analogue bias       : {summary['reading']['analogue_bias']}")
    print()
    print(f"written -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

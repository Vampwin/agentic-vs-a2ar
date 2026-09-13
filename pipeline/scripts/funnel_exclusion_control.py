"""
funnel_exclusion_control.py — build the docking ledger for the funnel-exclusion
control.

The manuscript argues that similarity-based applicability-domain gates wrongly
exclude novel chemotypes. The expansion funnel then uses 2D pharmacophore/ErG
similarity and 3D shape overlap to choose which 5,000 of 25,604 compounds get
docked. A reviewer is entitled to ask whether the funnel reintroduced exactly
the bias the paper criticises.

This builds a random sample of the compounds the funnel *rejected*, so they can
be docked under the identical fast-pass protocol and their score distribution
compared against the 5,000 that were kept. Two arms, because the funnel rejects
in two different ways:

  arm_2d : rejected by the 2D pharmacophore/ErG screen (21,882 -> 8,000).
           The larger cut and the more similarity-like filter.
  arm_3d : rejected by the 3D shape/pharmacophore refinement (8,000 -> 5,000).
           Already standardised, so usable directly.

Sampling is seeded and the seed recorded, so the control is reproducible and
cannot be re-drawn until it gives a convenient answer.

Usage:
    python funnel_exclusion_control.py --n-per-arm 250 --seed 20260821
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from kdss_paths import project_root, resolve, under

import numpy as np
import pandas as pd

def _root() -> Path:
    # Delegates to kdss_paths, which finds whichever layout this file is inside
    # -- the development tree or a deposit -- by walking up from its own
    # location. It previously tried two absolute paths and exited if neither
    # existed, which on the author's machine meant a copy of the deposit read
    # and wrote the development tree instead of itself.
    return project_root()


ROOT = _root()
FUNNEL = resolve("20_new_runs/np_library_funnel_2026-08-15")
POOL = resolve("02_frozen_inputs/lotus_np_library/lotus_pool_biosynthetic_25604.csv")
OUT = resolve("20_new_runs/funnel_exclusion_control_2026-08-21")

MW_MAX = 700.0          # identical to the funnel's stage-1 rule
ROTB_MAX = 10


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def standardize(smiles_list: list[str]) -> list[dict]:
    """Same operations the funnel's stage 1 applied: largest fragment,
    neutralise, canonical SMILES, InChIKey, MW and rotatable-bond count."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors
    from rdkit.Chem.MolStandardize import rdMolStandardize

    RDLogger.DisableLog("rdApp.*")
    lfc = rdMolStandardize.LargestFragmentChooser()
    un = rdMolStandardize.Uncharger()
    rows = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(str(smi))
        if mol is None:
            continue
        try:
            mol = un.uncharge(lfc.choose(mol))
            Chem.SanitizeMol(mol)
        except Exception:
            continue
        rows.append({
            "smiles_standardized": Chem.MolToSmiles(mol),
            "inchikey_standardized": Chem.MolToInchiKey(mol),
            "mw": float(Descriptors.MolWt(mol)),
            "rotatable_bonds": int(Descriptors.NumRotatableBonds(mol)),
        })
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-per-arm", type=int, default=250)
    ap.add_argument("--seed", type=int, default=20260821)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    stage2 = pd.read_csv(under("20_new_runs/np_library_funnel_2026-08-15", "np_library_stage2_2d_ranked.csv"))
    shortlist = pd.read_csv(under("20_new_runs/np_library_funnel_2026-08-15", "np_library_docking_shortlist.csv"))
    kept = set(shortlist["inchikey_standardized"])

    # --- arm_3d: in the 8,000 but not the 5,000 -----------------------------
    excl3 = stage2[~stage2["inchikey_standardized"].isin(kept)].copy()
    take = min(args.n_per_arm, len(excl3))
    arm3 = excl3.iloc[np.sort(rng.choice(len(excl3), size=take, replace=False))].copy()
    arm3["arm"] = "arm_3d"

    # --- arm_2d: in the pool but never reached the 8,000 ---------------------
    pool = pd.read_csv(POOL)
    smi_col = next(c for c in pool.columns if "smiles" in c.lower())
    ik_col = next((c for c in pool.columns if "inchikey" in c.lower()), None)
    reached = set(stage2["structure_inchikey"]) if ik_col else set()
    cand = pool[~pool[ik_col].isin(reached)] if ik_col else pool
    # oversample, because stage-1 property filtering will drop some
    over = min(len(cand), args.n_per_arm * 4)
    draw = cand.iloc[np.sort(rng.choice(len(cand), size=over, replace=False))]

    std = pd.DataFrame(standardize(draw[smi_col].tolist()))
    std = std[(std.mw <= MW_MAX) & (std.rotatable_bonds <= ROTB_MAX)]
    std = std.drop_duplicates("inchikey_standardized")
    std = std[~std.inchikey_standardized.isin(kept)]
    arm2 = std.head(args.n_per_arm).copy()
    arm2["arm"] = "arm_2d"

    cols = ["smiles_standardized", "inchikey_standardized", "mw", "rotatable_bonds", "arm"]
    both = pd.concat([arm2[cols], arm3[cols]], ignore_index=True)
    both.insert(0, "compound_id", [f"EXCL{i:05d}" for i in range(1, len(both) + 1)])
    both = both.rename(columns={"smiles_standardized": "SMILES"})
    ledger = under("20_new_runs/funnel_exclusion_control_2026-08-21", "exclusion_control_ledger.csv")
    both.to_csv(ledger, index=False)

    meta = {
        "purpose": "Control for the funnel's own similarity-based exclusions: dock a "
                   "seeded random sample of rejected compounds under the identical "
                   "fast-pass protocol and compare score distributions against the "
                   "5,000 that were kept.",
        "seed": args.seed,
        "n_per_arm_requested": args.n_per_arm,
        "arm_2d": {"n": int((both.arm == "arm_2d").sum()),
                   "population": "biosynthetic pool entries that never reached the 8,000 "
                                 "2D-screen survivors; standardised and property-filtered "
                                 "here with the same MW<=700 / rotB<=10 rule",
                   "oversampled_before_filter": int(over)},
        "arm_3d": {"n": int((both.arm == "arm_3d").sum()),
                   "population": "the 3,000 stage-2 survivors dropped by the 3D O3A refinement",
                   "population_size": int(len(excl3))},
        "inputs": {"pool": {"path": str(POOL), "sha256": sha256(POOL)},
                   "stage2": str(under("20_new_runs/np_library_funnel_2026-08-15", "np_library_stage2_2d_ranked.csv")),
                   "shortlist": str(under("20_new_runs/np_library_funnel_2026-08-15", "np_library_docking_shortlist.csv"))},
        "ledger": {"path": str(ledger), "sha256": sha256(ledger)},
        "docking_protocol": "fast pass, identical to the expansion screen: "
                            "2 ETKDGv3 conformers, exhaustiveness 8, same receptor and box",
    }
    (under("20_new_runs/funnel_exclusion_control_2026-08-21", "control_metadata.json")).write_text(json.dumps(meta, indent=2), encoding="utf-8")

    print(f"arm_2d : {int((both.arm=='arm_2d').sum())} compounds "
          f"(from {len(cand)} never-reached, oversampled {over})")
    print(f"arm_3d : {int((both.arm=='arm_3d').sum())} compounds (from {len(excl3)} dropped at 3D)")
    print(f"total  : {len(both)}")
    print(f"\nwritten -> {ledger}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

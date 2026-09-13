"""
build_metamorphic_ledgers.py -- construct the input variants for the
metamorphic-invariance tests.

A metamorphic test states a relation that must hold between the outputs of two
runs whose inputs differ in a way that cannot change the answer. It needs no
oracle, which is exactly why it suits agent-written code: it can be checked
without knowing the correct docking score.

Two transformations are produced, each isolating one dependency:

  M1  row-order permutation.  The ledger rows are shuffled at a fixed seed.
      Chemistry is untouched. If scores move, the pipeline carries state
      across compounds -- a shared RNG advanced per row, a reused conformer
      buffer, or a cache keyed by position.

  M2  SMILES re-writing.  Each SMILES is replaced by a different, valid,
      NON-canonical SMILES for the same molecule (RDKit doRandom). The
      molecular graph is identical; only the string is not. If scores move,
      something downstream depends on the input string rather than on the
      molecule -- atom ordering leaking into conformer generation, or a
      string-keyed cache.

Both are run through the identical docking command as production, and both are
compared against the SAME 25-compound subset already docked twice for
seed_repeatability, so the comparison is against an established baseline
rather than against a fresh reference.

An InChIKey identity check is enforced on M2 before anything is docked: if a
rewritten SMILES does not parse back to the same standard InChIKey, the
transformation itself is faulty and the test would measure the wrong thing.

Usage:
    python build_metamorphic_ledgers.py
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from kdss_paths import project_root, resolve, under

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


def _root() -> Path:
    # Delegates to kdss_paths, which finds whichever layout this file is inside
    # -- the development tree or a deposit -- by walking up from its own
    # location. It previously tried two absolute paths and exited if neither
    # existed, which on the author's machine meant a copy of the deposit read
    # and wrote the development tree instead of itself.
    return project_root()


ROOT = _root()
SRC = resolve("20_new_runs/g4_verification_2026-08-21/g4_subset_ledger.csv")
RUN = resolve("20_new_runs/metamorphic_2026-08-21")
SEED = 20260821


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> int:
    RUN.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(SRC)

    # ---- M1 : row-order permutation ---------------------------------------
    rng = np.random.default_rng(SEED)
    order = rng.permutation(len(df))
    m1 = df.iloc[order].reset_index(drop=True)
    p_m1 = under("20_new_runs/metamorphic_2026-08-21", "ledger_M1_rowshuffle.csv")
    m1.to_csv(p_m1, index=False)
    moved = int((np.asarray(order) != np.arange(len(df))).sum())

    # ---- M2 : non-canonical SMILES rewriting ------------------------------
    m2 = df.copy()
    rewritten, unchanged, failures = [], 0, []
    rr = np.random.default_rng(SEED)
    for i, smi in enumerate(df["SMILES"]):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            failures.append({"row": i, "reason": "input SMILES did not parse"})
            rewritten.append(smi)
            continue
        ref_key = Chem.MolToInchiKey(mol)
        new = smi
        # Re-root the traversal at a different atom rather than using doRandom:
        # this RDKit build's MolToSmiles takes no randomSeed, and rooting is
        # deterministic, so the rewritten ledger is reproducible from SEED
        # alone. Try several roots; a molecule with no degrees of freedom in
        # its writing order would otherwise silently be a no-op.
        roots = rr.permutation(mol.GetNumAtoms())
        for a in roots[:25]:
            cand = Chem.MolToSmiles(mol, canonical=False, rootedAtAtom=int(a))
            if cand != smi:
                new = cand
                break
        if new == smi:
            unchanged += 1
        back = Chem.MolFromSmiles(new)
        if back is None or Chem.MolToInchiKey(back) != ref_key:
            failures.append({"row": i, "compound_id": df.compound_id.iloc[i],
                             "reason": "rewritten SMILES is not the same molecule"})
            new = smi
        rewritten.append(new)
    m2["SMILES"] = rewritten
    p_m2 = under("20_new_runs/metamorphic_2026-08-21", "ledger_M2_smiles_rewrite.csv")
    m2.to_csv(p_m2, index=False)

    n_diff = int((m2.SMILES.values != df.SMILES.values).sum())

    meta = {
        "source_ledger": {"path": str(SRC), "sha256": sha256(SRC), "n_rows": int(len(df))},
        "seed": SEED,
        "M1": {"transformation": "row-order permutation at fixed seed",
               "path": str(p_m1), "sha256": sha256(p_m1),
               "rows_moved": moved,
               "invariant": "docking score must be identical per compound"},
        "M2": {"transformation": "SMILES replaced by a different valid "
                                 "non-canonical SMILES for the same molecule",
               "path": str(p_m2), "sha256": sha256(p_m2),
               "rows_rewritten": n_diff,
               "rows_unchanged_no_alternative_found": unchanged,
               "inchikey_identity_enforced": True,
               "transformation_failures": failures,
               "invariant": "docking score must be identical per compound"},
        "baseline_for_comparison":
            "20_new_runs/g4_verification_2026-08-21/rerun_passA_cpu6.csv "
            "(same 25 compounds, same command, already shown identical to "
            "production and to the 20-thread rerun)",
    }
    (under("20_new_runs/metamorphic_2026-08-21", "metamorphic_inputs_metadata.json")).write_text(
        json.dumps(meta, indent=2), encoding="utf-8")

    print(f"M1 rows moved            : {moved}/{len(df)}")
    print(f"M2 SMILES rewritten      : {n_diff}/{len(df)}")
    print(f"M2 no alternative found  : {unchanged}")
    print(f"M2 transformation errors : {len(failures)}")
    if failures:
        print(json.dumps(failures, indent=2))
    print(f"\nwritten -> {p_m1}\n           {p_m2}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

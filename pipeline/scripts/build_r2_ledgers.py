"""build_r2_ledgers.py -- R2 input: the 500 redocked compounds, rewritten.

§3.1 reports that the docking score depends on how the input SMILES is written,
measured on a 25-compound subset. R2 asks the question that matters for the
released result: does the *selection* survive the same transformation? The top
500 and the three Tier-A candidates were produced by the same machinery, and
their stability under re-representation is unmeasured.

The transformation is the one already used for M2 and is reproduced here rather
than reinvented: re-root the SMILES traversal at a different atom, keeping the
molecular graph identical, and refuse any rewrite whose standard InChIKey does
not match the original. A rewrite that changes the molecule would make the test
measure the wrong thing.

Each production chunk is rewritten row-for-row and in the same order, so R2
chunk k pairs exactly with production chunk k. That pairing is what lets the
comparison be per-compound rather than distributional.

Usage:
    python build_r2_ledgers.py [--out DIR]
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


def _root() -> Path:
    marker = "00_project_control/08_DECISION_LOG.md"
    for cand in (Path(r"D:\Agentic VS Workflow"), Path("/mnt/d/Agentic VS Workflow")):
        if (cand / marker).is_file():
            return cand
    raise SystemExit("project root not found")


ROOT = _root()
SRC_DIR = ROOT / "20_new_runs/np_library_funnel_2026-08-15"
SEED = 20260826          # distinct from the M2 seed; R2 is a different draw


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def rewrite(smi: str, rng: np.random.Generator) -> tuple[str, str]:
    """Return (rewritten_smiles, status). Identity-checked by InChIKey."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return smi, "input_did_not_parse"
    ref = Chem.MolToInchiKey(mol)
    for a in rng.permutation(mol.GetNumAtoms())[:25]:
        cand = Chem.MolToSmiles(mol, canonical=False, rootedAtAtom=int(a))
        if cand == smi:
            continue
        back = Chem.MolFromSmiles(cand)
        if back is not None and Chem.MolToInchiKey(back) == ref:
            return cand, "rewritten"
        # a candidate that does not round-trip is discarded, not accepted
    return smi, "no_valid_alternative_found"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "20_new_runs/smiles_redock_500_2026-08-26"))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(SEED)
    chunks, counts = {}, {"rewritten": 0, "no_valid_alternative_found": 0,
                          "input_did_not_parse": 0}
    total = 0
    for k in (1, 2, 3):
        src = SRC_DIR / f"full_ledger_chunk{k}.csv"
        df = pd.read_csv(src)
        new, status = [], []
        for smi in df["SMILES"]:
            s, st = rewrite(str(smi), rng)
            new.append(s)
            status.append(st)
            counts[st] += 1
        out_df = df.copy()
        out_df["SMILES_original"] = df["SMILES"]
        out_df["SMILES"] = new
        out_df["rewrite_status"] = status
        dst = out / f"ledger_R2_chunk{k}.csv"
        out_df.to_csv(dst, index=False)
        chunks[f"chunk{k}"] = {
            "source": str(src.relative_to(ROOT).as_posix()),
            "source_sha256": sha256(src),
            "path": str(dst.relative_to(ROOT).as_posix()),
            "sha256": sha256(dst),
            "n_rows": int(len(df)),
            "n_rewritten": int(sum(s == "rewritten" for s in status)),
        }
        total += len(df)
        print(f"chunk{k}: {len(df)} rows, "
              f"{sum(s == 'rewritten' for s in status)} rewritten -> {dst.name}")

    meta = {
        "purpose": "R2 -- does the top-500 selection and the Tier-A set survive "
                   "re-representation of the input SMILES?",
        "adr": "ADR-AIIM-133",
        "seed": SEED,
        "n_total": total,
        "transformation": "MolToSmiles(canonical=False, rootedAtAtom=a) for a "
                          "random atom a; the molecular graph is unchanged and "
                          "only the string differs",
        "identity_check": "standard InChIKey of the rewritten SMILES must equal "
                          "that of the original, else the rewrite is discarded "
                          "and the original is kept",
        "status_counts": counts,
        "chunks": chunks,
        "pairing": "R2 chunk k is production chunk k rewritten row-for-row in "
                   "the same order, so comparison is per-compound.",
        "baseline_for_comparison":
            "20_new_runs/g4_verification_2026-08-21/full_scores_MERGED.csv",
    }
    (out / "r2_inputs_metadata.json").write_text(json.dumps(meta, indent=2),
                                                 encoding="utf-8")
    print(f"\nstatus counts: {counts}")
    print(f"written -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

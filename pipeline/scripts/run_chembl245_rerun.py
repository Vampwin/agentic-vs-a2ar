"""run_chembl245_rerun.py -- the ChEMBL 245 benchmark under the unified protocol.

Pays the debt ADR-AIIM-215 took on. That amendment changed the docking protocol
from eight pooled conformers to one minimised conformer for the DUD-E run, and
promised that the ChEMBL 245 benchmark would be re-run the same way so the two
benchmarks could be compared without the protocol confounding the comparison.
The author has since fixed the protocol for all docking from here on.

**Protocol identity is structural, not asserted.** Every setting and the docking
function itself are *imported* from `run_dude_enrichment.py`. There is one
implementation, so the two benchmarks cannot drift apart through an edit to one
of them, and a reviewer can check that in a single import line rather than by
diffing two parameter blocks.

Input is the ledger the original run docked from,
`20_new_runs/mmgbsa_pilot_2026-08-10/mmgbsa_full245_ledger.csv`: 245 compounds,
50 active and 195 inactive, whose labels agree with the published benchmark at
245 of 245. The docked poses of the 2026-08-15 run agree with these SMILES for
221 of 245; the remaining 24 differ by protonation or tautomer introduced during
preparation, which is expected and is why the ledger rather than the pose file
is the input.

Nothing here reads the published AUC. This script docks and records; the
comparison is a separate step.

    python run_chembl245_rerun.py --cpu 20

Output -> 20_new_runs/chembl245_rerun_2026-09-22/
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
import time
from pathlib import Path

import pandas as pd
from vina import Vina

sys.path.insert(0, str(Path(__file__).resolve().parent))

# One protocol, one implementation. Importing rather than restating is the
# whole point of this module.
from run_dude_enrichment import (  # noqa: E402
    CENTER,
    EXHAUSTIVENESS,
    N_CONFORMERS_DOCKED,
    N_CONFORMERS_GENERATED,
    N_POSES,
    RECEPTOR,
    SEED,
    SIZE,
    dock_one,
)

ROOT = Path(__file__).resolve().parents[2]
LEDGER = ROOT / "20_new_runs/mmgbsa_pilot_2026-08-10/mmgbsa_full245_ledger.csv"
OUTDIR = ROOT / "20_new_runs/chembl245_rerun_2026-09-22"

FIELDS = ["compound_id", "label", "pchembl", "vina_score",
          "n_conformers_docked", "elapsed_s", "error"]


def done_ids(path: Path) -> set[str]:
    if not path.exists():
        return set()
    try:
        return set(pd.read_csv(path).compound_id.astype(str))
    except Exception:                                        # noqa: BLE001
        return set()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpu", type=int, default=20)
    ap.add_argument("--limit", type=int, default=0, help="smoke test only")
    args = ap.parse_args()

    OUTDIR.mkdir(parents=True, exist_ok=True)
    out = OUTDIR / "scores.csv"

    d = pd.read_csv(LEDGER)
    if args.limit:
        d = d.head(args.limit)
    already = done_ids(out)
    todo = d[~d.compound_id.astype(str).isin(already)]

    print(f"ChEMBL 245 re-run under the unified protocol")
    print(f"  conformers generated {N_CONFORMERS_GENERATED}, "
          f"docked {N_CONFORMERS_DOCKED} (lowest MMFF94 energy)")
    print(f"  exhaustiveness {EXHAUSTIVENESS}, {N_POSES} poses, seed {SEED}")
    print(f"  receptor {RECEPTOR.name}, box centre {CENTER}, size {SIZE}")
    print(f"  {len(d)} compounds  done {len(already)}  to do {len(todo)}  "
          f"cpu={args.cpu}", flush=True)
    if todo.empty:
        print("nothing to do", flush=True)
        return 0

    v = Vina(sf_name="vina", cpu=args.cpu, seed=SEED, verbosity=0)
    v.set_receptor(str(RECEPTOR))
    t = time.time()
    v.compute_vina_maps(center=CENTER, box_size=SIZE)
    print(f"maps computed in {time.time()-t:.1f}s", flush=True)

    new = not out.exists()
    fh = open(out, "a", newline="", encoding="utf-8")
    w = csv.DictWriter(fh, fieldnames=FIELDS)
    if new:
        w.writeheader()
        fh.flush()

    t0 = time.time()
    for i, (_, r) in enumerate(todo.iterrows(), 1):
        c0 = time.time()
        score, n, err = dock_one(v, r.SMILES)
        w.writerow({"compound_id": r.compound_id, "label": int(r.label),
                    "pchembl": r.get("pChEMBL", ""),
                    "vina_score": "" if score is None else round(score, 4),
                    "n_conformers_docked": n,
                    "elapsed_s": round(time.time() - c0, 2),
                    "error": err or ""})
        fh.flush()
        os.fsync(fh.fileno())
        if i % 25 == 0 or i == len(todo):
            el = time.time() - t0
            rate = el / i
            print(f"  {i}/{len(todo)}  {rate:.1f}s/cpd  "
                  f"eta {(len(todo)-i)*rate/60:.0f} min", flush=True)
    fh.close()
    print(f"complete in {(time.time()-t0)/60:.1f} min", flush=True)

    # The plan's standing gate: a mixed-protocol comparison must abort, not warn.
    got = pd.read_csv(out)
    scored = got[got.vina_score.notna()]
    assert (scored.n_conformers_docked == 1).all(), (
        "a row was docked with more than one conformer; this run cannot be "
        "compared with the DUD-E benchmark")
    print(f"protocol gate: all {len(scored)} scored rows docked 1 conformer")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""run_dude_enrichment.py -- dock the DUD-E A2AR benchmark under the frozen protocol.

Executes M5.1, pre-registered in ADR-AIIM-214 before any compound was docked.

Protocol, amended before any result was analysed (ADR-AIIM-215):

    8 ETKDGv3 conformers generated and MMFF94 minimised,
    the LOWEST-ENERGY one docked
    exhaustiveness 32, 10 retained poses
    seed 20260729
    receptor snp_screen_2026-08-04/G2_preparation/4EIY_prepared.pdbqt
    box centre [-0.424, 8.527, 17.127], 22 A cube

The production protocol docked all eight conformers and kept the best score.
Measured here at 22.7 s per dock, that is 121 hours for this benchmark and not
affordable. Docking only the minimum-energy conformer costs 15.1 hours.

The scientific case for the change is not only cost. Vina performs
**flexible-ligand** docking: it samples rotatable-bond torsions internally
during its Monte Carlo search, so supplying eight input conformers is largely
redundant for exactly the degrees of freedom Vina already explores. What input
conformers do fix is ring geometry, which Vina holds rigid -- and for
drug-like molecules with aromatic or small saturated rings a single minimised
conformer is standard practice. Eight pooled conformers is a conservative
choice, above the field norm, not a requirement.

The cost is stated rather than hidden: this validates a standard-practice
protocol that differs in one respect from the eight-conformer protocol that
produced the study's own screening results. `bridge_8conformer_reference.csv`
holds seven compounds scored under the production protocol so the two can be
compared directly.

Runs under WSL `diffdock` (vina 1.2.7). Shardable and resumable: each worker
appends to its own CSV and skips compounds already present, so an interrupted
run continues where it stopped.

    python run_dude_enrichment.py --shard 0 --nshards 4 --cpu 5

Output -> 20_new_runs/dude_enrichment_2026-09-22/
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
import time
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from meeko import MoleculePreparation
from vina import Vina

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[2]
SET = ROOT / "02_frozen_inputs/public_benchmarks/dude_a2ar/dude_a2ar_docking_set.csv"
RECEPTOR = ROOT / "20_new_runs/snp_screen_2026-08-04/G2_preparation/4EIY_prepared.pdbqt"
OUTDIR = ROOT / "20_new_runs/dude_enrichment_2026-09-22"

CENTER = [-0.424, 8.527, 17.127]
SIZE = [22, 22, 22]
N_CONFORMERS_GENERATED = 8      # generated and minimised
N_CONFORMERS_DOCKED = 1         # only the lowest-energy one is docked
EXHAUSTIVENESS = 32
N_POSES = 10
SEED = 20260729

FIELDS = ["compound_id", "label", "vina_score", "n_conformers_docked",
          "elapsed_s", "error"]


def done_ids(path: Path) -> set[str]:
    if not path.exists():
        return set()
    try:
        d = pd.read_csv(path)
        return set(d.compound_id.astype(str))
    except Exception:                                        # noqa: BLE001
        return set()


def dock_one(v: Vina, smiles: str,
             exhaustiveness: int = EXHAUSTIVENESS
             ) -> tuple[float | None, int, str | None]:
    """Dock the minimum-MMFF-energy conformer and return its rank-1 energy.

    `exhaustiveness` defaults to the frozen protocol's 32 and exists so that the
    M2.2 arm (ADR-AIIM-231) can vary the one setting it is about without
    forking this function. Callers that do not pass it -- including
    `run_chembl245_rerun.py`, which imports this function so that the two
    benchmarks cannot drift apart (ADR-AIIM-227) -- are unaffected.
    """
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None, 0, "unparseable SMILES"
    m = Chem.AddHs(m)
    cids = AllChem.EmbedMultipleConfs(m, numConfs=N_CONFORMERS_GENERATED,
                                      randomSeed=SEED)
    if len(cids) == 0:
        return None, 0, "embedding failed"
    # MMFFOptimizeMoleculeConfs returns [(not_converged, energy), ...] in the
    # order of the conformer ids, so the argmin picks the lowest-energy one.
    try:
        res = AllChem.MMFFOptimizeMoleculeConfs(m)
        energies = [e for _conv, e in res]
        best_cid = int(cids[int(min(range(len(energies)),
                                    key=lambda i: energies[i]))])
    except Exception:                                        # noqa: BLE001
        best_cid = int(cids[0])
    try:
        sub = Chem.Mol(m, confId=best_cid)
        prep = MoleculePreparation()
        prep.prepare(sub)
        v.set_ligand_from_string(prep.write_pdbqt_string())
        v.dock(exhaustiveness=exhaustiveness, n_poses=N_POSES)
        return float(v.energies(n_poses=1)[0][0]), N_CONFORMERS_DOCKED, None
    except Exception as exc:                                 # noqa: BLE001
        return None, 0, f"dock failed: {type(exc).__name__}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--cpu", type=int, default=5)
    ap.add_argument("--limit", type=int, default=0, help="smoke test only")
    ap.add_argument("--exhaustiveness", type=int, default=EXHAUSTIVENESS,
                    help="M2.2 arm (ADR-AIIM-231). The frozen protocol is 32; "
                         "Vina's own default is 8. Nothing else changes.")
    ap.add_argument("--out", default=str(OUTDIR),
                    help="output directory; defaults to the frozen M5.1 run, "
                         "so a resume of that run needs no argument")
    args = ap.parse_args()

    outdir = Path(args.out)

    # A row docked at a different exhaustiveness must never land in the M5.1
    # artifact: the whole point of the M2.2 arm is that the two sets of scores
    # are comparable because exactly one thing differs between them, and a
    # mixed directory destroys that. Refuse rather than warn.
    if args.exhaustiveness != EXHAUSTIVENESS and outdir == OUTDIR:
        print(f"refusing: --exhaustiveness {args.exhaustiveness} into the frozen "
              f"M5.1 directory. Give --out a new dated directory.", flush=True)
        return 2

    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / f"scores_shard{args.shard}.csv"

    d = pd.read_csv(SET)
    mine = d.iloc[args.shard::args.nshards].reset_index(drop=True)
    if args.limit:
        mine = mine.head(args.limit)
    already = done_ids(out)
    todo = mine[~mine.compound_id.astype(str).isin(already)]

    print(f"shard {args.shard}/{args.nshards}  assigned {len(mine)}  "
          f"done {len(already)}  to do {len(todo)}  cpu={args.cpu}  "
          f"exhaustiveness={args.exhaustiveness}  out={outdir.name}", flush=True)
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
        score, n, err = dock_one(v, r.SMILES, args.exhaustiveness)
        w.writerow({"compound_id": r.compound_id, "label": int(r.label),
                    "vina_score": "" if score is None else round(score, 4),
                    "n_conformers_docked": n,
                    "elapsed_s": round(time.time() - c0, 2),
                    "error": err or ""})
        fh.flush()
        os.fsync(fh.fileno())
        if i % 25 == 0 or i == len(todo):
            el = time.time() - t0
            rate = el / i
            print(f"  [{args.shard}] {i}/{len(todo)}  {rate:.1f}s/cpd  "
                  f"eta {(len(todo)-i)*rate/3600:.1f}h", flush=True)
    fh.close()
    print(f"shard {args.shard} complete in {(time.time()-t0)/3600:.2f}h", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""
g4_dock_library.py — multi-conformer Vina docking of a compound library.

Runs inside the WSL `diffdock` env (Vina + meeko + RDKit). Invoked by the G4
adapter in run_snp_screen.py, but usable standalone for calibration.

Protocol is the one G3 validated, and it is not optional: the redocking control
only reproduced the ZM241385 crystal pose when structural waters were retained
in the receptor *and* several starting conformers were docked with their poses
pooled by energy. Screening under a looser protocol than the one that passed
validation would make the validation meaningless.

Every compound terminates in exactly one state — `ok` or a recorded failure
reason — so the accountability ledger G4 checks can never silently lose one.

Results are appended per compound and the run resumes from whatever is already
on disk, because a full library takes hours and a crash should not cost all of it.

Usage
-----
    python g4_dock_library.py --ledger L.csv --receptor R.pdbqt --box B.json \
        --out scores.csv [--limit N] [--ids SNP00001,SNP00002]
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolDescriptors
from meeko import MoleculePreparation, PDBQTMolecule
from vina import Vina

RDLogger.DisableLog("rdApp.*")

FIELDS = [
    "compound_id", "status", "failure_reason",
    "vina_score", "vina_mean_top3", "n_conformers_generated", "n_poses_pooled",
    "best_conformer_index", "centroid_dist_A", "on_pocket",
    "mol_formula", "n_heavy_atoms", "n_rotatable_bonds", "elapsed_s",
]


def centroid(mol: Chem.Mol) -> tuple[float, float, float]:
    conf = mol.GetConformer()
    n = mol.GetNumAtoms()
    return tuple(sum(conf.GetAtomPosition(i)[k] for i in range(n)) / n for k in range(3))


def dock_one(smiles: str, receptor: str, box: dict, cfg: dict,
             poses_dir: Path | None = None, compound_id: str = "") -> dict:
    """Dock a single compound across N conformers; return the pooled best.

    When `poses_dir` is given the pooled poses are also written to
    `<poses_dir>/<compound_id>.sdf`, energy-ranked, so that downstream
    pose-dependent analyses (interaction fingerprints, contact profiling)
    can run without re-docking. Purely additive: with `poses_dir=None` the
    returned record is bit-identical to the pre-existing behaviour.
    """
    t0 = time.time()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"status": "failed", "failure_reason": "invalid_smiles"}

    molh = Chem.AddHs(mol)
    ps = AllChem.ETKDGv3()
    ps.randomSeed = int(cfg["seed"])
    ps.pruneRmsThresh = float(cfg["prune_rms"])
    cids = AllChem.EmbedMultipleConfs(molh, numConfs=int(cfg["n_conformers"]), params=ps)
    if not len(cids):
        # Retry once with random coordinates before giving up: ETKDG can fail on
        # strained macrocycles that are otherwise perfectly dockable.
        ps.useRandomCoords = True
        cids = AllChem.EmbedMultipleConfs(molh, numConfs=int(cfg["n_conformers"]), params=ps)
    if not len(cids):
        return {"status": "failed", "failure_reason": "embedding_failure"}

    try:
        AllChem.MMFFOptimizeMoleculeConfs(molh)
    except Exception:
        pass  # unminimised conformers still dock; geometry quality is recorded downstream

    cx, cy, cz = box["center"]
    sx, sy, sz = box["size"]
    pool: list[tuple[float, int, Chem.Mol]] = []

    for ci in cids:
        single = Chem.Mol(molh)
        single.RemoveAllConformers()
        single.AddConformer(molh.GetConformer(ci), assignId=True)
        try:
            prep = MoleculePreparation()
            prep.prepare(single)
            pdbqt = prep.write_pdbqt_string()
        except Exception:
            continue  # this conformer could not be typed; others may still work

        # Vina aborts the whole process — not a catchable Python exception — if
        # handed an empty or atomless ligand file, so one untypeable conformer
        # would otherwise kill the entire library run. Validate before calling it.
        if not pdbqt or not pdbqt.strip() or "ATOM" not in pdbqt:
            continue

        lig_path = "/tmp/_g4_lig.pdbqt"
        with open(lig_path, "w") as fh:
            fh.write(pdbqt)

        try:
            v = Vina(sf_name="vina", cpu=int(cfg.get("cpu", 0)), seed=int(cfg["seed"]), verbosity=0)
            v.set_receptor(receptor)
            v.set_ligand_from_file(lig_path)
            v.compute_vina_maps(center=[cx, cy, cz], box_size=[sx, sy, sz])
            v.dock(exhaustiveness=int(cfg["exhaustiveness"]), n_poses=int(cfg["n_poses"]))
            out_path = "/tmp/_g4_out.pdbqt"
            v.write_poses(out_path, n_poses=int(cfg["n_poses"]), overwrite=True)
            energies = [float(e[0]) for e in v.energies()]
        except Exception:
            continue

        try:
            for k, pose in enumerate(PDBQTMolecule.from_file(out_path, skip_typing=True)):
                mols = pose.export_rdkit_mol()
                m = mols[0] if isinstance(mols, (list, tuple)) else mols
                m = Chem.RemoveHs(m)
                pool.append((energies[k] if k < len(energies) else 0.0, int(ci), m))
        except Exception:
            continue

    if not pool:
        return {"status": "failed", "failure_reason": "docking_failure"}

    # Energy-ranked across every conformer: the pooled best is the compound's score.
    pool.sort(key=lambda t: t[0])
    best_e, best_ci, best_mol = pool[0]
    top3 = [e for e, _, _ in pool[:3]]
    bx, by, bz = centroid(best_mol)
    dist = math.dist((bx, by, bz), (cx, cy, cz))

    if poses_dir is not None:
        # Retaining poses is a provenance decision, not a cosmetic one: the
        # score alone cannot support any later interaction-level analysis.
        n_save = int(cfg.get("n_save_poses", 10))
        try:
            poses_dir.mkdir(parents=True, exist_ok=True)
            with Chem.SDWriter(str(poses_dir / f"{compound_id}.sdf")) as w:
                for rank, (e, ci, m) in enumerate(pool[:n_save], 1):
                    m.SetProp("_Name", f"{compound_id}_pose{rank}")
                    m.SetDoubleProp("vina_energy", float(e))
                    m.SetIntProp("conformer_index", int(ci))
                    m.SetIntProp("pose_rank", rank)
                    w.write(m)
        except Exception as exc:  # never let pose archiving kill a library run
            print(f"  warn: pose write failed for {compound_id}: {exc}", flush=True)

    heavy = Chem.RemoveHs(mol)
    return {
        "status": "ok",
        "failure_reason": "",
        "vina_score": round(best_e, 3),
        "vina_mean_top3": round(sum(top3) / len(top3), 3),
        "n_conformers_generated": len(cids),
        "n_poses_pooled": len(pool),
        "best_conformer_index": best_ci,
        "centroid_dist_A": round(dist, 2),
        "on_pocket": dist <= float(cfg["on_pocket_cutoff"]),
        "mol_formula": rdMolDescriptors.CalcMolFormula(heavy),
        "n_heavy_atoms": heavy.GetNumAtoms(),
        "n_rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(heavy),
        "elapsed_s": round(time.time() - t0, 1),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ledger", required=True)
    ap.add_argument("--receptor", required=True)
    ap.add_argument("--box", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--id-col", default="SNP_ID")
    ap.add_argument("--smiles-col", default="smiles_canonical")
    ap.add_argument("--n-conformers", type=int, default=8)
    ap.add_argument("--exhaustiveness", type=int, default=32)
    ap.add_argument("--n-poses", type=int, default=10)
    ap.add_argument("--seed", type=int, default=20260729)
    ap.add_argument("--prune-rms", type=float, default=0.5)
    ap.add_argument("--on-pocket-cutoff", type=float, default=8.0)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--ids", default="")
    ap.add_argument("--poses-dir", default="",
                    help="if set, archive each compound's energy-ranked poses as "
                         "<dir>/<id>.sdf for later interaction-level analysis")
    ap.add_argument("--n-save-poses", type=int, default=10)
    ap.add_argument("--cpu", type=int, default=0,
                    help="threads passed to Vina; 0 (default) means Vina's own "
                         "'use all available cores' behavior, unchanged from "
                         "prior runs. Set explicitly when running multiple "
                         "chunks in parallel on the same machine, so each "
                         "chunk doesn't independently try to claim every "
                         "core and leave no headroom for the host.")
    args = ap.parse_args()

    import pandas as pd

    df = pd.read_csv(args.ledger)
    if args.ids:
        wanted = {s.strip() for s in args.ids.split(",") if s.strip()}
        df = df[df[args.id_col].astype(str).isin(wanted)]
    if args.limit:
        df = df.head(args.limit)

    box = json.loads(Path(args.box).read_text())
    if "centers_by_receptor" in box:
        box = {"center": box["centers_by_receptor"]["4EIY"], "size": box["size"]}

    cfg = {
        "n_conformers": args.n_conformers, "exhaustiveness": args.exhaustiveness,
        "n_poses": args.n_poses, "seed": args.seed, "prune_rms": args.prune_rms,
        "on_pocket_cutoff": args.on_pocket_cutoff,
        "n_save_poses": args.n_save_poses, "cpu": args.cpu,
    }
    poses_dir = Path(args.poses_dir) if args.poses_dir else None

    out_path = Path(args.out)
    done: set[str] = set()
    if out_path.exists():
        with open(out_path, newline="", encoding="utf-8") as fh:
            done = {r["compound_id"] for r in csv.DictReader(fh)}
        print(f"resuming: {len(done)} compounds already scored", flush=True)
    else:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w", newline="", encoding="utf-8") as fh:
            csv.DictWriter(fh, fieldnames=FIELDS).writeheader()

    total = len(df)
    t_start = time.time()
    for i, row in enumerate(df.itertuples(index=False), 1):
        cid = str(getattr(row, args.id_col))
        if cid in done:
            continue
        smi = getattr(row, args.smiles_col, None)
        if not isinstance(smi, str) or not smi.strip():
            rec = {"status": "failed", "failure_reason": "missing_smiles"}
        else:
            rec = dock_one(smi.strip(), args.receptor, box, cfg,
                           poses_dir=poses_dir, compound_id=cid)
        rec["compound_id"] = cid
        with open(out_path, "a", newline="", encoding="utf-8") as fh:
            csv.DictWriter(fh, fieldnames=FIELDS).writerow(
                {k: rec.get(k, "") for k in FIELDS}
            )
        rate = (time.time() - t_start) / max(1, i)
        print(f"[{i}/{total}] {cid} {rec['status']} "
              f"{rec.get('vina_score', '')} eta={(total - i) * rate / 60:.0f}min",
              flush=True)

    print("done", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

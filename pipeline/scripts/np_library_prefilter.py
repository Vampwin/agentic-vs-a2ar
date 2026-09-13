"""
np_library_prefilter.py — reduce the LOTUS biosynthetic-class pool
(25,604 compounds; ADR-AIIM-079) to a docking-affordable shortlist.

Funnel, per ADR-AIIM-073's original compute budget (~5,000-10,000 docked):

  1. Standardize + property filter (RDKit only)
     - parse SMILES, keep the largest fragment (strip salts/counterions),
       recompute canonical SMILES and InChIKey
     - drop MW > 700 Da (ZM241385 = 337 Da; the A2AR orthosteric pocket is
       compact -- SNP-Lib's own on-pocket exclusion problem was partly a
       size/shape issue, not a reason to carry very large glycosides forward)
     - drop > 10 rotatable bonds (docking pose reliability degrades sharply
       with flexibility, and this benchmark's own baseline AUC is already
       weak -- carrying maximally flexible ligands forward would only add
       docking noise on top of a signal that is already thin)
     - dedupe by InChIKey (post-standardization; some entries in the raw
       LOTUS pool differ only by stereochemistry annotation or salt form)

  2. 2D pharmacophore / reduced-graph similarity screen (Gobbi 2D + ErG,
     the same descriptor spaces validated in E2's multi-space AD sweep)
     against a small reference set of known A2AR ligands: ZM241385, 9XT,
     and the ChEMBL benchmark's actual actives (n=50, pChEMBL >= 6.0) --
     not just the two co-crystal ligands, so the reference captures more of
     the chemotype diversity the benchmark itself already contains. Ranks
     the full standardized pool by max similarity to any reference; this is
     a coarse triage step, not a claim that 2D pharmacophore similarity
     predicts activity (E2 already showed such similarity spaces reject a
     structurally novel library at the far tail -- here it is used only to
     order a candidate pool, not to hard-cut it against a domain threshold).

  3. 3D shape+pharmacophore refinement (RDKit O3A) on the top slice from
     step 2, against the same reference ligands' 3D conformers, to bring in
     a signal step 2's 2D fingerprints cannot see (matching earlier project
     findings that phenanthrene- and xanthine-class scaffolds can overlap
     substantially in 3D while sharing almost no 2D fingerprint bits).

The final shortlist is written with full provenance (which stage kept it,
its rank and score at each stage) so the funnel is auditable end to end --
matches this project's registry discipline even though this is a triage
step rather than a G0-registered pipeline run.

Usage:
    python np_library_prefilter.py --pool 02_frozen_inputs/lotus_np_library/lotus_pool_biosynthetic_25604.csv \
        --n-2d 8000 --n-final 5000 --out 20_new_runs/np_library_funnel_2026-08-15
"""
from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
from kdss_paths import project_root, resolve, under

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger, DataStructs
from rdkit.Chem import (AllChem, Descriptors, rdFingerprintGenerator,
                        rdMolDescriptors, rdReducedGraphs)
from rdkit.Chem.Pharm2D import Generate, Gobbi_Pharm2D
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")

ROOT = project_root()
CHEMBL_TRAINING = resolve("02_frozen_inputs/public_benchmarks/chembl_a2ar/chembl_a2a_training_v3.csv")
ZM241385_SMILES = "Nc1nc(NCCc2ccc(O)cc2)nc2nc(-c3ccco3)nn12"
STRUCT_DIR = resolve("02_frozen_inputs/a2ar_structures")

MW_MAX = 700.0
ROTB_MAX = 10
SEED = 20260815


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def standardize(smi: str) -> tuple[str | None, Chem.Mol | None]:
    mol = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
    if mol is None:
        return None, None
    try:
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
        mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
        uncharger = rdMolStandardize.Uncharger()
        mol = uncharger.uncharge(mol)
        Chem.SanitizeMol(mol)
    except Exception:
        return None, None
    return Chem.MolToSmiles(mol), mol


def get_9xt_smiles() -> str:
    """9XT crystal ligand, from HETNAM (3-[(4-azanyl-3-methyl-phenyl)methyl]-
    7-(furan-2-yl)-[1,2,3]triazolo[4,5-d]pyrimidin-5-amine), matched to the
    5OLH structure already used elsewhere in this project."""
    return "Nc1nc(-c2ccco2)c2nnn(Cc3ccc(N)c(C)c3)c2n1"


def load_references(n_chembl_actives: int = 50) -> list[str]:
    refs = [ZM241385_SMILES, get_9xt_smiles()]
    tr = pd.read_csv(CHEMBL_TRAINING)
    actives = tr[tr["label"] == 1].sort_values("pchembl_value", ascending=False)
    refs += actives["canonical_smiles"].head(n_chembl_actives).tolist()
    return refs


# ---- Stage 2 (2D pharmacophore) worker pool --------------------------------
_WORKER_REF_GOBBI: list = []
_WORKER_REF_ERG: list = []


def _init_2d_worker(ref_smiles: list[str]) -> None:
    global _WORKER_REF_GOBBI, _WORKER_REF_ERG
    from rdkit import Chem as _Chem
    from rdkit.Chem import rdReducedGraphs as _rdReducedGraphs
    from rdkit.Chem.Pharm2D import Generate as _Generate, Gobbi_Pharm2D as _Gobbi_Pharm2D
    mols = [m for m in (_Chem.MolFromSmiles(s) for s in ref_smiles) if m is not None]
    _WORKER_REF_GOBBI = [_Generate.Gen2DFingerprint(m, _Gobbi_Pharm2D.factory) for m in mols]
    _WORKER_REF_ERG = [np.array(_rdReducedGraphs.GetErGFingerprint(m)) for m in mols]


def _score_one_2d(smi: str) -> float:
    from rdkit import Chem as _Chem
    from rdkit import DataStructs as _DataStructs
    from rdkit.Chem import rdReducedGraphs as _rdReducedGraphs
    from rdkit.Chem.Pharm2D import Generate as _Generate, Gobbi_Pharm2D as _Gobbi_Pharm2D
    m = _Chem.MolFromSmiles(smi)
    if m is None:
        return -1.0
    try:
        fp = _Generate.Gen2DFingerprint(m, _Gobbi_Pharm2D.factory)
        gobbi_sims = _DataStructs.BulkTanimotoSimilarity(fp, _WORKER_REF_GOBBI)
        erg = np.array(_rdReducedGraphs.GetErGFingerprint(m))
        erg_sims = [
            float(np.dot(erg, r) / (np.linalg.norm(erg) * np.linalg.norm(r) + 1e-9))
            for r in _WORKER_REF_ERG
        ]
        return max(0.5 * g + 0.5 * e for g, e in zip(gobbi_sims, erg_sims))
    except Exception:
        return -1.0


# ---- Stage 3 (O3A) worker pool ---------------------------------------------
# Stage 2 already triaged the full pool against all 52 references; stage 3's
# job is to refine that ranking with 3D shape+pharmacophore information, so a
# smaller, still-diverse reference subset is sufficient here and keeps the
# per-worker embedding cost down. Reference conformers are embedded ONCE per
# worker process (via the pool initializer) rather than once per compound.
_WORKER_REF_CONFS: list = []


def _init_o3a_worker(ref_smiles: list[str], seed: int) -> None:
    global _WORKER_REF_CONFS
    from rdkit import Chem as _Chem
    from rdkit.Chem import AllChem as _AllChem
    _WORKER_REF_CONFS = []
    for s in ref_smiles:
        m = _Chem.MolFromSmiles(s)
        if m is None:
            continue
        mh = _Chem.AddHs(m)
        if _AllChem.EmbedMolecule(mh, randomSeed=seed) == 0:
            _AllChem.MMFFOptimizeMolecule(mh)
            _WORKER_REF_CONFS.append(mh)


def _score_one_o3a(args: tuple[str, int]) -> float:
    smi, seed = args
    from rdkit import Chem as _Chem
    from rdkit.Chem import AllChem as _AllChem
    m = _Chem.MolFromSmiles(smi)
    if m is None:
        return -1.0
    mh = _Chem.AddHs(m)
    if _AllChem.EmbedMolecule(mh, randomSeed=seed) != 0:
        return -1.0
    try:
        _AllChem.MMFFOptimizeMolecule(mh)
    except Exception:
        pass
    best = -1.0
    for rc in _WORKER_REF_CONFS:
        try:
            o3a = _AllChem.GetO3A(mh, rc)
            best = max(best, o3a.Score())
        except Exception:
            continue
    return best


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool", required=True)
    ap.add_argument("--n-2d", type=int, default=8000)
    ap.add_argument("--n-final", type=int, default=5000)
    ap.add_argument("--n-refs-3d", type=int, default=12,
                    help="reference subset size for the (expensive) 3D O3A stage; "
                         "stage 2 already used the full reference set for 2D triage")
    ap.add_argument("--n-workers", type=int, default=max(1, (os.cpu_count() or 4) - 2))
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    pool = pd.read_csv(args.pool)
    print(f"input pool: {len(pool)} compounds "
          f"(sha256={sha256(Path(args.pool))[:16]}...)")

    # ---- Stage 1: standardize + property filter ----------------------------
    std_smiles, std_mols, keep_idx = [], [], []
    for i, smi in enumerate(pool["structure_smiles"]):
        s, m = standardize(smi)
        if s is None:
            continue
        mw = Descriptors.MolWt(m)
        rotb = rdMolDescriptors.CalcNumRotatableBonds(m)
        if mw > MW_MAX or rotb > ROTB_MAX:
            continue
        std_smiles.append(s)
        std_mols.append(m)
        keep_idx.append(i)

    s1 = pool.iloc[keep_idx].copy()
    s1["smiles_standardized"] = std_smiles
    s1["mw"] = [Descriptors.MolWt(m) for m in std_mols]
    s1["rotatable_bonds"] = [rdMolDescriptors.CalcNumRotatableBonds(m) for m in std_mols]
    s1["inchikey_standardized"] = [Chem.MolToInchiKey(m) for m in std_mols]
    before = len(pool)
    s1 = s1.drop_duplicates("inchikey_standardized").reset_index(drop=True)
    print(f"stage 1 (standardize + MW<={MW_MAX} + rotb<={ROTB_MAX} + dedupe): "
          f"{before} -> {len(s1)}")

    # ---- Stage 2: 2D pharmacophore / ErG similarity screen (parallel) -------
    refs = load_references()
    ref_mols = [m for m in (Chem.MolFromSmiles(s) for s in refs) if m is not None]
    print(f"reference set: {len(ref_mols)} known A2AR ligands "
          f"(ZM241385, 9XT, {len(ref_mols)-2} ChEMBL actives)")
    print(f"running 2D pharmacophore/ErG screen on {len(s1)} compounds using "
          f"{args.n_workers} worker processes...")

    smiles1 = s1["smiles_standardized"].tolist()
    with mp.Pool(processes=args.n_workers, initializer=_init_2d_worker,
                initargs=(refs,)) as mpool:
        max_sim = mpool.map(_score_one_2d, smiles1,
                            chunksize=max(1, len(smiles1) // (args.n_workers * 4)))

    s1["pharmacophore_2d_score"] = np.round(max_sim, 4)
    s2 = s1.sort_values("pharmacophore_2d_score", ascending=False).head(args.n_2d).reset_index(drop=True)
    print(f"stage 2 (2D pharmacophore/ErG top-{args.n_2d}): "
          f"{len(s1)} -> {len(s2)}  score range [{s2['pharmacophore_2d_score'].min():.3f}, "
          f"{s2['pharmacophore_2d_score'].max():.3f}]")

    # ---- Stage 3: 3D shape+pharmacophore refinement (O3A, parallel) ---------
    # ZM241385 and 9XT are always kept (the two structurally distinct A2AR
    # binding modes documented in this project's own work); the remaining
    # slots go to the most potent ChEMBL actives, already sorted that way in
    # load_references().
    refs_3d = [ZM241385_SMILES, get_9xt_smiles()] + refs[2:args.n_refs_3d]
    print(f"stage 3 reference subset: {len(refs_3d)} ligands "
          f"(full 2D triage already used all {len(ref_mols)})")
    print(f"running O3A on {len(s2)} compounds using {args.n_workers} worker processes...")

    tasks = [(smi, SEED) for smi in s2["smiles_standardized"]]
    with mp.Pool(processes=args.n_workers, initializer=_init_o3a_worker,
                initargs=(refs_3d, SEED)) as mpool:
        o3a_scores = mpool.map(_score_one_o3a, tasks,
                               chunksize=max(1, len(tasks) // (args.n_workers * 4)))

    s2["shape_pharm_3d_score"] = np.round(o3a_scores, 4)
    valid = s2[s2["shape_pharm_3d_score"] > 0]
    print(f"stage 3: {len(valid)}/{len(s2)} embedded successfully; "
          f"O3A score range [{valid['shape_pharm_3d_score'].min():.2f}, "
          f"{valid['shape_pharm_3d_score'].max():.2f}]")

    final = valid.sort_values("shape_pharm_3d_score", ascending=False).head(args.n_final).reset_index(drop=True)
    final["funnel_rank"] = range(1, len(final) + 1)
    print(f"\nfinal shortlist: {len(final)} compounds "
          f"({int(final['is_orchidaceae'].sum())} Orchidaceae-linked)")

    out_cols = ["structure_inchikey", "inchikey_standardized", "smiles_standardized",
               "mw", "rotatable_bonds", "structure_taxonomy_npclassifier_02superclass",
               "structure_taxonomy_npclassifier_03class", "organism_name", "is_orchidaceae",
               "pharmacophore_2d_score", "shape_pharm_3d_score", "funnel_rank"]
    final[out_cols].to_csv(outdir / "np_library_docking_shortlist.csv", index=False)
    s2.to_csv(outdir / "np_library_stage2_2d_ranked.csv", index=False)

    meta = dict(
        input_pool=args.pool, input_pool_sha256=sha256(Path(args.pool)),
        n_input=len(pool), n_stage1=len(s1), n_stage2=len(s2), n_final=len(final),
        mw_max=MW_MAX, rotatable_bonds_max=ROTB_MAX,
        reference_set_size=len(ref_mols), reference_smiles=refs,
        combination_rule="0.5*Gobbi2D_Tanimoto + 0.5*ErG_cosine, fixed a priori",
        seed=SEED,
        note="Triage/prefilter only. Ranking by 2D/3D similarity to known A2AR "
             "ligands orders which compounds get scarce docking compute first; "
             "it is not itself a validated activity predictor (see ADR-AIIM-075/076: "
             "IFP/shape-similarity rescoring did not beat plain Vina on the benchmark). "
             "Final activity evidence is the Vina docking step that follows.",
    )
    (outdir / "funnel_metadata.json").write_text(json.dumps(meta, indent=2, default=str))
    print(f"\nwritten -> {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

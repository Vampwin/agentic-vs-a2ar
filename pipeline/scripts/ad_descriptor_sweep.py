"""E2 — Is the applicability domain a property of the library, or of the descriptor?

The frozen G6 rule declares a compound in-domain when its mean Tanimoto to the
5 nearest training-set neighbours reaches the 10th percentile of the training
set's own 5-NN self-similarity. That rule is applied here unchanged, but in
several descriptor spaces rather than ECFP4 alone, so that "0/236 out of domain"
can be read as either a fact about the library or an artefact of one descriptor
choice.

Every space uses the identical rule, identical k, identical percentile; only the
representation changes. No thresholds are tuned and no space is preferred.

Usage:
    python ad_descriptor_sweep.py --out <dir>
"""
from __future__ import annotations

import argparse
import json
import hashlib
from pathlib import Path
from kdss_paths import project_root, resolve, under

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger, DataStructs
from rdkit.Chem import MACCSkeys, rdFingerprintGenerator, rdReducedGraphs, Descriptors
from rdkit.Chem.Pharm2D import Generate, Gobbi_Pharm2D

RDLogger.DisableLog("rdApp.*")

def _root() -> Path:
    # Delegates to kdss_paths, which finds whichever layout this file is inside
    # -- the development tree or a deposit -- by walking up from its own
    # location. It previously tried two absolute paths and exited if neither
    # existed, which on the author's machine meant a copy of the deposit read
    # and wrote the development tree instead of itself.
    return project_root()


ROOT = _root()
TRAIN = resolve("02_frozen_inputs/public_benchmarks/chembl_a2ar/chembl_a2a_training_v3.csv")
SNPLIB = resolve("02_frozen_inputs/snp_library/01_library_validated.csv")

K = 5            # frozen: 5 nearest neighbours
PCTL = 10        # frozen: 10th percentile of training self-similarity

PHYSCHEM = ["MolWt", "MolLogP", "TPSA", "NumHDonors", "NumHAcceptors",
            "NumRotatableBonds", "NumAromaticRings", "FractionCSP3",
            "RingCount", "HeavyAtomCount"]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# ---------------- representations ----------------
def bitvec_matrix(mols, kind):
    """Return an (n, nbits) float32 matrix of binary fingerprints."""
    if kind == "ECFP4":
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        fps = [gen.GetFingerprint(m) for m in mols]
    elif kind == "FCFP4":
        inv = rdFingerprintGenerator.GetMorganFeatureAtomInvGen()
        gen = rdFingerprintGenerator.GetMorganGenerator(
            radius=2, fpSize=2048, atomInvariantsGenerator=inv)
        fps = [gen.GetFingerprint(m) for m in mols]
    elif kind == "MACCS":
        fps = [MACCSkeys.GenMACCSKeys(m) for m in mols]
    elif kind == "AtomPair":
        gen = rdFingerprintGenerator.GetAtomPairGenerator(fpSize=2048)
        fps = [gen.GetFingerprint(m) for m in mols]
    elif kind == "TopTorsion":
        gen = rdFingerprintGenerator.GetTopologicalTorsionGenerator(fpSize=2048)
        fps = [gen.GetFingerprint(m) for m in mols]
    elif kind == "Gobbi_Pharm2D":
        fps = [Generate.Gen2DFingerprint(m, Gobbi_Pharm2D.factory) for m in mols]
    else:
        raise ValueError(kind)
    n = len(fps)
    nb = fps[0].GetNumBits()
    X = np.zeros((n, nb), dtype=np.float32)
    for i, fp in enumerate(fps):
        for b in fp.GetOnBits():
            X[i, b] = 1.0
    return X


def realvec_matrix(mols, kind):
    """Return an (n, d) float32 matrix of real-valued descriptors."""
    if kind == "ErG":
        V = np.array([rdReducedGraphs.GetErGFingerprint(m) for m in mols],
                     dtype=np.float32)
    elif kind == "PhysChem":
        rows = []
        for m in mols:
            rows.append([getattr(Descriptors, d)(m) for d in PHYSCHEM])
        V = np.array(rows, dtype=np.float32)
    else:
        raise ValueError(kind)
    return V


def tanimoto_block(A, B, chunk=512):
    """Binary Tanimoto between rows of A (n,d) and B (m,d). Returns (n,m)."""
    a = A.sum(axis=1)[:, None]
    b = B.sum(axis=1)[None, :]
    out = np.empty((A.shape[0], B.shape[0]), dtype=np.float32)
    for s in range(0, A.shape[0], chunk):
        e = min(s + chunk, A.shape[0])
        inter = A[s:e] @ B.T
        union = a[s:e] + b - inter
        out[s:e] = np.where(union > 0, inter / np.maximum(union, 1e-9), 0.0)
    return out


def cosine_block(A, B, chunk=512):
    """Cosine similarity for real-valued vectors, mapped to [0,1]."""
    An = A / np.maximum(np.linalg.norm(A, axis=1, keepdims=True), 1e-9)
    Bn = B / np.maximum(np.linalg.norm(B, axis=1, keepdims=True), 1e-9)
    out = np.empty((A.shape[0], B.shape[0]), dtype=np.float32)
    for s in range(0, A.shape[0], chunk):
        e = min(s + chunk, A.shape[0])
        out[s:e] = An[s:e] @ Bn.T
    return np.clip(out, 0.0, 1.0)


def knn_mean(S, k, exclude_self):
    """Mean similarity to the k nearest neighbours, row-wise."""
    if exclude_self:
        np.fill_diagonal(S, -np.inf)
    part = np.partition(S, -k, axis=1)[:, -k:]
    return part.mean(axis=1)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "20_new_runs" /
                                         "ad_descriptor_sweep_2026-08-15"))
    # Additive: default reproduces the original SNP-Lib run byte-for-byte.
    # Supplied so the identical rule can be applied to the expanded library,
    # which the original run never covered (declared Limitation #2).
    ap.add_argument("--library", default=str(SNPLIB))
    ap.add_argument("--smiles-col", default="SMILES_std")
    ap.add_argument("--id-col", default="Compound")
    ap.add_argument("--group-col", default="Groups")
    args = ap.parse_args()
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    tr = pd.read_csv(TRAIN)
    sn = pd.read_csv(args.library)

    tr_mols, tr_keep = [], []
    for i, s in enumerate(tr["canonical_smiles"].astype(str)):
        m = Chem.MolFromSmiles(s)
        if m is not None:
            tr_mols.append(m)
            tr_keep.append(i)
    sn_mols, sn_ids, sn_grp = [], [], []
    for _, r in sn.iterrows():
        m = Chem.MolFromSmiles(str(r[args.smiles_col]))
        if m is not None:
            sn_mols.append(m)
            sn_ids.append(r[args.id_col])
            sn_grp.append(r[args.group_col] if args.group_col in sn.columns else "")

    print(f"training molecules parsed : {len(tr_mols)} / {len(tr)}")
    print(f"SNP-Lib molecules parsed  : {len(sn_mols)} / {len(sn)}")
    print(f"rule: in-domain if mean Tanimoto to {K}-NN >= {PCTL}th percentile "
          f"of training self-similarity (frozen G6 rule, unchanged)\n")

    binary = ["ECFP4", "FCFP4", "MACCS", "AtomPair", "TopTorsion", "Gobbi_Pharm2D"]
    real = ["ErG", "PhysChem"]

    rows, per_compound, self_sims = [], {}, {}
    for kind in binary + real:
        if kind in binary:
            A = bitvec_matrix(tr_mols, kind)
            B = bitvec_matrix(sn_mols, kind)
            S_tr = tanimoto_block(A, A)
            S_sn = tanimoto_block(B, A)
        else:
            A = realvec_matrix(tr_mols, kind)
            B = realvec_matrix(sn_mols, kind)
            if kind == "PhysChem":            # standardise before distance
                mu, sd = A.mean(0), A.std(0) + 1e-9
                A = (A - mu) / sd
                B = (B - mu) / sd
            S_tr = cosine_block(A, A)
            S_sn = cosine_block(B, A)

        self_sim = knn_mean(S_tr.copy(), K, exclude_self=True)
        thr = float(np.percentile(self_sim, PCTL))
        # Additive, ADR-AIIM-126. This vector is the only thing that makes
        # coverage comparable ACROSS spaces, and until now it was computed on
        # the line above and thrown away. Collapsed to a single percentile, each
        # space's threshold sits at an unknown position in its own similarity
        # distribution -- so "0 % in ECFP4, 86 % in PhysChem" cannot be
        # separated from those two distributions simply having different
        # shapes. Retained in full below; nothing above this comment changes,
        # and the three legacy outputs are byte-identical as a result.
        self_sims[kind] = self_sim
        lib_sim = knn_mean(S_sn.copy(), K, exclude_self=False)
        inside = lib_sim >= thr

        rows.append(dict(
            descriptor_space=kind,
            threshold=round(thr, 4),
            snplib_mean_knn_sim=round(float(lib_sim.mean()), 4),
            snplib_min=round(float(lib_sim.min()), 4),
            snplib_max=round(float(lib_sim.max()), 4),
            n_inside=int(inside.sum()),
            n_total=len(sn_mols),
            coverage_fraction=round(float(inside.mean()), 4),
        ))
        per_compound[kind] = dict(sim=lib_sim.tolist(), inside=inside.tolist(),
                                  threshold=thr)
        print(f"{kind:>14}  thr={thr:.4f}  lib mean={lib_sim.mean():.4f} "
              f"[{lib_sim.min():.3f}-{lib_sim.max():.3f}]  "
              f"in-domain={int(inside.sum()):>3}/{len(sn_mols)} "
              f"({inside.mean()*100:.1f}%)")

    res = pd.DataFrame(rows)
    res.to_csv(outdir / "ad_coverage_by_descriptor.csv", index=False)

    pc = pd.DataFrame({"compound_id": sn_ids, "group": sn_grp})
    for kind in per_compound:
        pc[f"{kind}_sim"] = np.round(per_compound[kind]["sim"], 4)
        pc[f"{kind}_inside"] = per_compound[kind]["inside"]
    pc["n_spaces_inside"] = pc[[f"{k}_inside" for k in per_compound]].sum(axis=1)
    pc.to_csv(outdir / "ad_per_compound.csv", index=False)

    meta = dict(
        rule=dict(k=K, percentile=PCTL,
                  definition="mean Tanimoto/cosine to k nearest training neighbours "
                             ">= percentile of training self-similarity"),
        inputs={str(TRAIN.name): sha256(TRAIN), Path(args.library).name: sha256(Path(args.library))},
        n_training=len(tr_mols), n_library=len(sn_mols),
        spaces=binary + real,
        note="Identical rule in every space; only the representation differs. "
             "No threshold was tuned and no space was preferred a priori.",
    )
    (outdir / "run_metadata.json").write_text(json.dumps(meta, indent=2))

    # ---- additive outputs (ADR-AIIM-126) ----------------------------------
    # Written AFTER run_metadata.json so the three legacy artifacts are
    # produced by exactly the code path that produced them before.
    np.savez_compressed(outdir / "self_sim_by_space.npz",
                        **{k: v.astype(np.float64) for k, v in self_sims.items()})

    # Where does each library compound sit inside the training set's own
    # self-similarity distribution, in that same space? This is a rank, so it
    # is comparable between spaces in a way the raw similarities are not.
    pct = pd.DataFrame({"compound_id": sn_ids, "group": sn_grp})
    for kind, ss in self_sims.items():
        ss_sorted = np.sort(ss)
        lib = np.asarray(per_compound[kind]["sim"], dtype=np.float64)
        pct[f"{kind}_pctile"] = np.round(
            np.searchsorted(ss_sorted, lib, side="left") / len(ss_sorted) * 100.0, 4)
    pct.to_csv(outdir / "ad_percentile_position.csv", index=False)

    # Provenance the original run did not record. run_metadata.json is left
    # exactly as it was so it stays byte-comparable; this is a separate file.
    import platform, sys as _sys, os as _os, rdkit
    (outdir / "run_provenance.json").write_text(json.dumps(dict(
        argv=_sys.argv,
        python=_sys.version,
        platform=platform.platform(),
        rdkit=rdkit.__version__,
        numpy=np.__version__,
        pandas=pd.__version__,
        thread_env={k: _os.environ.get(k) for k in
                    ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                     "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")},
        note="Recorded because the 2026-08-21 run recorded neither tooling "
             "versions nor the command line, and two RDKit installations "
             "differ across this project.",
    ), indent=2))

    print(f"\nby group, spaces-in-domain (max {len(per_compound)}):")
    print(pc.groupby("group")["n_spaces_inside"].agg(["mean", "max", "count"])
            .sort_values("mean", ascending=False).head(12).to_string())
    print(f"\nwritten -> {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

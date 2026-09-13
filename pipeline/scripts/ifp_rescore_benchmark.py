"""
ifp_rescore_benchmark.py — does structure-based rescoring beat Vina alone on
the corrected-protocol ChEMBL A2AR benchmark?

Pre-registered per ADR-AIIM-073 before this script was run against the full
benchmark: three enrichment techniques with documented A2AR/GPCR precedent
are evaluated retrospectively against the frozen n=245 benchmark, with the
existing corrected-protocol Vina AUC of 0.601 as the comparator.

  1. Asn253 hydrogen-bond filter — A2AR-specific precedent (Congreve et al.,
     J Med Chem 2012; enrichment of known ligands tracks Asn253 contact).
  2. Interaction-fingerprint (IFP) similarity to the ZM241385 reference
     binding mode — general GPCR-VS precedent (Marcou & Rognan, JCIM 2007;
     Sanders et al., JCIM 2012).
  3. An equal-weighted z-score ensemble of (-Vina energy) and IFP similarity
     — the literature's "complementary hypothesis" for combining an
     energy-based score with an independent structural signal. The
     combination weights are fixed a priori (1:1, no tuning against the
     benchmark labels) precisely so this comparison cannot become circular.

Scoring methods compared, all on the identical n (whatever the benchmark
docking run actually produced), all against the same frozen labels:

  vina_only                  existing baseline (best pooled Vina energy)
  vina_at_best_ifp_pose       Vina energy of the pose most similar to the
                             reference binding mode, instead of the lowest-
                             energy pose (does binding-mode fidelity beat
                             pure energy for POSE SELECTION alone?)
  ifp_only                    max IFP-Tanimoto similarity across a compound's
                             saved poses (a purely structural signal)
  asn253_filtered_vina         compounds forming the Asn253 H-bond in any
                             saved pose are ranked above those that do not;
                             within each group, ranked by Vina energy
  combined_zscore              z(-vina_energy) + z(ifp_similarity), 1:1,
                             fixed before this script saw the labels

The interaction fingerprint itself is purely geometric (per-residue contact,
H-bond, and aromaticity bits from atomic coordinates only) — it needs no
bond-order template for the candidate poses, only their 3D coordinates,
which g4_dock_library.py's --poses-dir archives per compound.

Usage:
    python ifp_rescore_benchmark.py \
        --scores 20_new_runs/chembl245_ifp_rescoring_2026-08-15/scores.csv \
        --poses-dir 20_new_runs/chembl245_ifp_rescoring_2026-08-15/poses \
        --ledger 20_new_runs/chembl248_redock_corrected_2026-08-10/chembl248_ledger.csv \
        --out 20_new_runs/chembl245_ifp_rescoring_2026-08-15
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from kdss_paths import project_root, resolve, under

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.ML.Scoring import Scoring

RDLogger.DisableLog("rdApp.*")

ROOT = project_root()
CRYSTAL_PDB = resolve("02_frozen_inputs/a2ar_structures/4EIY.pdb")
RECEPTOR_PDB = (resolve("20_new_runs/routeA_water_shell_2026-08-09")
                / "4EIY_cut5_protein_plus_waters.pdb")
TEMPLATE_SDF = (resolve("20_new_runs/routeA_water_shell_2026-08-09")
               / "ZMA_redocked_4EIY_cut5.sdf")

CONTACT, HB = 4.5, 3.5      # frozen geometric cutoffs, see _domain_gate_bias_probe siblings
BOOTSTRAP_N, BOOTSTRAP_SEED = 1000, 42     # matches the project's existing benchmark convention
EF_FRACTIONS = [0.01, 0.05, 0.10]
BEDROC_ALPHA = 20.0


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# ---------------- geometry ----------------
def read_pdb(path: Path) -> list[dict]:
    atoms = []
    for line in open(path):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        atoms.append(dict(
            name=line[12:16].strip(), resn=line[17:20].strip(),
            chain=line[21].strip(), resi=line[22:26].strip(),
            x=float(line[30:38]), y=float(line[38:46]), z=float(line[46:54]),
            elem=(line[76:78].strip() or line[12:16].strip()[0]),
        ))
    return atoms


def xyz(atoms: list[dict]) -> np.ndarray:
    return np.array([[a["x"], a["y"], a["z"]] for a in atoms], dtype=float)


def ifp_bits(lig_xyz: np.ndarray, lig_elem: list[str], prot_heavy: list[dict],
             prot_xyz: np.ndarray) -> set[str]:
    """Per-residue contact/polar-Hbond/aromatic bits; waters collapsed to one bucket."""
    D = np.linalg.norm(prot_xyz[None, :, :] - lig_xyz[:, None, :], axis=2)
    mind = D.min(axis=0)
    polar = [i for i, e in enumerate(lig_elem) if e in ("N", "O")]
    bits: dict[str, list[int]] = {}
    for k, a in enumerate(prot_heavy):
        if mind[k] > CONTACT:
            continue
        key = "HOH" if a["resn"] in ("HOH", "WAT") else f'{a["resn"]}{a["resi"]}'
        b = bits.setdefault(key, [0, 0, 0])
        b[0] = 1
        if a["elem"] in ("N", "O") and polar and D[polar, k].min() <= HB:
            b[1] = 1
        if a["resn"] in ("PHE", "TYR", "TRP", "HIS"):
            b[2] = 1
    return {f"{r}:{i}" for r, b in bits.items() for i, v in enumerate(b) if v}


def tanimoto(a: set, b: set) -> float:
    return len(a & b) / len(a | b) if (a | b) else 0.0


def asn253_hbond(lig_xyz: np.ndarray, lig_elem: list[str], prot_heavy: list[dict],
                 prot_xyz: np.ndarray) -> bool:
    polar = [i for i, e in enumerate(lig_elem) if e in ("N", "O")]
    if not polar:
        return False
    sc_idx = [k for k, a in enumerate(prot_heavy)
             if a["resn"] == "ASN" and a["resi"] == "253" and a["name"] in ("ND2", "OD1")]
    if not sc_idx:
        return False
    d = np.linalg.norm(prot_xyz[sc_idx][None, :, :] - lig_xyz[polar][:, None, :], axis=2)
    return bool(d.min() <= HB)


# ---------------- scoring metrics ----------------
def _sorted_label_rows(score_label: list[tuple[float, int]]) -> list[list[int]]:
    """RDKit's Scoring module wants rows pre-sorted, index 0 = best-ranked."""
    ordered = sorted(score_label, key=lambda t: -t[0])
    return [[int(lbl)] for _, lbl in ordered]


def compute_metrics(score_label: list[tuple[float, int]]) -> dict:
    rows = _sorted_label_rows(score_label)
    n = len(rows)
    n_pos = sum(r[0] for r in rows)
    out = {"n": n, "n_active": n_pos, "n_inactive": n - n_pos}
    if n_pos == 0 or n_pos == n:
        out["error"] = "degenerate label distribution"
        return out
    out["auc_roc"] = round(float(Scoring.CalcAUC(rows, 0)), 4)
    out["bedroc20"] = round(float(Scoring.CalcBEDROC(rows, 0, BEDROC_ALPHA)), 4)
    ef = Scoring.CalcEnrichment(rows, 0, EF_FRACTIONS)
    for frac, val in zip(EF_FRACTIONS, ef):
        out[f"ef_{int(frac*100)}pct"] = round(float(val), 4)
    return out


def bootstrap_ci(score_label: list[tuple[float, int]], metric: str,
                 n_boot: int = BOOTSTRAP_N, seed: int = BOOTSTRAP_SEED) -> list[float] | None:
    rng = np.random.default_rng(seed)
    n = len(score_label)
    vals = []
    arr = np.array(score_label, dtype=object)
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        sample = [tuple(arr[i]) for i in idx]
        labels = {l for _, l in sample}
        if len(labels) < 2:
            continue
        m = compute_metrics(sample)
        if metric in m:
            vals.append(m[metric])
    if not vals:
        return None
    return [round(float(np.percentile(vals, 2.5)), 4), round(float(np.percentile(vals, 97.5)), 4)]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", required=True)
    ap.add_argument("--poses-dir", required=True)
    ap.add_argument("--ledger", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--no-bootstrap", action="store_true",
                    help="skip CIs for a fast partial-data check")
    args = ap.parse_args()

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    poses_dir = Path(args.poses_dir)

    # ---- reference IFP from the crystal ZM241385 pose (not a redocked one) ----
    zma_lines = [l for l in open(CRYSTAL_PDB)
                if l.startswith("HETATM") and l[17:20].strip() == "ZMA" and l[21] == "A"]
    cryst = Chem.MolFromPDBBlock("".join(zma_lines) + "END\n", removeHs=False, sanitize=False)
    cryst_xyz = np.array(cryst.GetConformer().GetPositions())
    cryst_elem = [a.GetSymbol() for a in cryst.GetAtoms()]

    prot = read_pdb(RECEPTOR_PDB)
    prot_heavy = [a for a in prot if a["elem"] != "H"]
    prot_xyz = xyz(prot_heavy)

    ref_ifp = ifp_bits(cryst_xyz, cryst_elem, prot_heavy, prot_xyz)
    print(f"reference IFP: {len(ref_ifp)} bits from crystal ZM241385 pose "
          f"({RECEPTOR_PDB.name})")
    print(f"reference Asn253 H-bond: {asn253_hbond(cryst_xyz, cryst_elem, prot_heavy, prot_xyz)}\n")

    # ---- labels ----
    ledger = pd.read_csv(args.ledger)[["compound_id", "label"]].drop_duplicates("compound_id")
    labels = dict(zip(ledger["compound_id"].astype(str), ledger["label"].astype(int)))

    # ---- per-compound scores ----
    scores = pd.read_csv(args.scores)
    scores = scores[scores["status"] == "ok"].copy()
    print(f"docked compounds available: {len(scores)}")

    rows = []
    missing_poses = 0
    for _, r in scores.iterrows():
        cid = str(r["compound_id"])
        if cid not in labels:
            continue
        pose_path = poses_dir / f"{cid}.sdf"
        if not pose_path.is_file():
            missing_poses += 1
            continue
        poses = [m for m in Chem.SDMolSupplier(str(pose_path), removeHs=False) if m is not None]
        if not poses:
            missing_poses += 1
            continue

        best_ifp, best_ifp_rank, best_ifp_vina = -1.0, None, None
        any_hbond, hbond_at_rank1, n_hbond_poses = False, False, 0
        for p in poses:
            pos = np.array(p.GetConformer().GetPositions())
            elem = [a.GetSymbol() for a in p.GetAtoms()]
            fp = ifp_bits(pos, elem, prot_heavy, prot_xyz)
            sim = tanimoto(fp, ref_ifp)
            rank = p.GetPropsAsDict().get("pose_rank")
            e = p.GetPropsAsDict().get("vina_energy")
            if sim > best_ifp:
                best_ifp, best_ifp_rank, best_ifp_vina = sim, rank, e
            hb = asn253_hbond(pos, elem, prot_heavy, prot_xyz)
            any_hbond = any_hbond or hb
            n_hbond_poses += int(hb)
            if rank == 1:
                hbond_at_rank1 = hb

        rows.append(dict(
            compound_id=cid, label=labels[cid],
            vina_score=float(r["vina_score"]),
            n_poses=len(poses),
            max_ifp_tanimoto=round(best_ifp, 4),
            best_ifp_pose_rank=best_ifp_rank,
            vina_at_best_ifp_pose=best_ifp_vina,
            any_asn253_hbond=any_hbond,
            asn253_hbond_at_rank1=hbond_at_rank1,
            n_hbond_poses=n_hbond_poses,
            frac_hbond_poses=round(n_hbond_poses / len(poses), 3),
        ))

    df = pd.DataFrame(rows)
    print(f"compounds with archived poses and a label: {len(df)}  "
          f"(missing poses: {missing_poses})")
    if df.empty:
        print("nothing to score yet -- run again once docking has produced poses")
        return 0
    print(f"labels: {int(df['label'].sum())} active / {int((df['label']==0).sum())} inactive\n")

    # ---- z-score ensemble (fixed 1:1 weights, computed before seeing labels) ----
    def z(x: pd.Series) -> pd.Series:
        return (x - x.mean()) / (x.std(ddof=0) + 1e-9)

    df["z_vina"] = z(-df["vina_score"])
    df["z_ifp"] = z(df["max_ifp_tanimoto"])
    df["combined_zscore"] = df["z_vina"] + df["z_ifp"]

    # asn253-filtered vina: hbond-forming compounds ranked above non-forming ones,
    # then by vina energy within each group. Encoded as a single sortable score.
    df["asn253_filtered_score"] = (
        df["any_asn253_hbond"].astype(int) * 1000.0 + (-df["vina_score"])
    )

    methods = {
        "vina_only": -df["vina_score"],
        "vina_at_best_ifp_pose": -df["vina_at_best_ifp_pose"].astype(float),
        "ifp_only": df["max_ifp_tanimoto"],
        "asn253_filtered_vina": df["asn253_filtered_score"],
        "combined_zscore": df["combined_zscore"],
    }

    summary = {}
    print(f"{'method':>22} {'n':>5} {'AUC':>7} {'95% CI':>16} {'BEDROC20':>9} "
          f"{'EF@1%':>7} {'EF@5%':>7} {'EF@10%':>7}")
    print("-" * 100)
    for name, score_col in methods.items():
        sl = list(zip(score_col.tolist(), df["label"].tolist()))
        m = compute_metrics(sl)
        if "error" in m:
            print(f"{name:>22}  -- degenerate labels, skipped --")
            summary[name] = m
            continue
        ci = None if args.no_bootstrap else bootstrap_ci(sl, "auc_roc")
        m["auc_ci95"] = ci
        summary[name] = m
        ci_str = f"[{ci[0]:.3f}, {ci[1]:.3f}]" if ci else "n/a"
        print(f"{name:>22} {m['n']:>5} {m['auc_roc']:>7.3f} {ci_str:>16} "
              f"{m['bedroc20']:>9.3f} {m.get('ef_1pct', float('nan')):>7.2f} "
              f"{m.get('ef_5pct', float('nan')):>7.2f} {m.get('ef_10pct', float('nan')):>7.2f}")

    baseline_auc = summary.get("vina_only", {}).get("auc_roc")
    print(f"\nbaseline (corrected-protocol Vina, prior run): AUC 0.601 "
          f"[0.510-0.683], n=245")
    if baseline_auc is not None:
        print(f"baseline (this run, vina_only, n={summary['vina_only']['n']}): "
              f"AUC {baseline_auc:.3f}")
        for name, m in summary.items():
            if name == "vina_only" or "error" in m:
                continue
            delta = m["auc_roc"] - baseline_auc
            verdict = "BEATS" if delta > 0 else ("TIES" if delta == 0 else "does not beat")
            print(f"  {name:>22}: delta AUC = {delta:+.3f} vs vina_only  ({verdict} baseline)")

    df.to_csv(outdir / "ifp_rescoring_per_compound.csv", index=False)
    meta = dict(
        n_scored=len(df), n_missing_poses=missing_poses,
        reference_receptor=str(RECEPTOR_PDB), reference_receptor_sha256=sha256(RECEPTOR_PDB),
        reference_ifp_n_bits=len(ref_ifp),
        bootstrap_n=BOOTSTRAP_N, bootstrap_seed=BOOTSTRAP_SEED,
        bedroc_alpha=BEDROC_ALPHA, ef_fractions=EF_FRACTIONS,
        combination_rule="combined_zscore = z(-vina_energy) + z(max_ifp_tanimoto), "
                         "1:1 fixed weights, not fit against labels",
        note="Pre-registered per ADR-AIIM-073 before this script was run against the "
             "full n=245 benchmark. Baseline for comparison: corrected-protocol Vina "
             "AUC 0.601 [0.510-0.683] (Section 3.7 of the prior manuscript draft).",
    )
    (outdir / "ifp_rescoring_summary.json").write_text(
        json.dumps({"metadata": meta, "methods": summary}, indent=2, default=str))
    print(f"\nwritten -> {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

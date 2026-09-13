"""
Quick variant of the pre-registered IFP rescoring (ADR-AIIM-075): does a
two-reference IFP beat the single-ZM241385 reference?

ADR-AIIM-075 found max-IFP-Tanimoto to a single ZM241385 crystal reference
carries essentially zero discriminative power (mean identical to 3 decimals
between actives and inactives). One plausible cause: a single reference
binding mode cannot represent the diversity among 50 chemically varied
actives. This project already has a second, independent, structurally
distinct A2AR co-crystal ligand frozen as an input: 9XT in 5OLH, which the
main analysis found binds without any water bridge (unlike ZM241385, which
requires one) -- Section 3.2 of the prior manuscript draft. Using both
references, rather than one, is a quick, well-motivated variant that needs
no new docking: it only needs 9XT's crystal coordinates re-expressed in
4EIY's frame so they are directly comparable to the already-archived
candidate poses.

Method: Kabsch superposition of 5OLH onto 4EIY using matched C-alpha atoms
(chain A, matched by residue number), applied to 9XT's crystal HETATM block.
This is a rigid-body coordinate transform, not a re-dock -- it carries no
new pose uncertainty beyond the crystal structures themselves.

Usage:
    python multi_reference_ifp_variant.py --out <dir>
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from kdss_paths import project_root, resolve, under

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.ML.Scoring import Scoring

RDLogger.DisableLog("rdApp.*")

ROOT = project_root()
STRUCT_DIR = resolve("02_frozen_inputs/a2ar_structures")
CRYSTAL_4EIY = under("02_frozen_inputs/a2ar_structures", "4EIY.pdb")
CRYSTAL_5OLH = under("02_frozen_inputs/a2ar_structures", "5OLH.pdb")
RECEPTOR_PDB = (resolve("20_new_runs/routeA_water_shell_2026-08-09")
                / "4EIY_cut5_protein_plus_waters.pdb")
RUN_DIR = resolve("20_new_runs/chembl245_ifp_rescoring_2026-08-15")

CONTACT, HB = 4.5, 3.5
BOOTSTRAP_N, BOOTSTRAP_SEED = 1000, 42
EF_FRACTIONS = [0.01, 0.05, 0.10]
BEDROC_ALPHA = 20.0


def read_pdb(path: Path, hetres: str | None = None) -> list[dict]:
    out = []
    for line in open(path):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        resn = line[17:20].strip()
        if hetres is not None and resn != hetres:
            continue
        if hetres is None and line.startswith("HETATM"):
            continue
        out.append(dict(name=line[12:16].strip(), resn=resn, chain=line[21].strip(),
                        resi=line[22:26].strip(),
                        x=float(line[30:38]), y=float(line[38:46]), z=float(line[46:54]),
                        elem=(line[76:78].strip() or line[12:16].strip()[0])))
    return out


def ca_by_resi(atoms: list[dict], chain: str = "A") -> dict[str, np.ndarray]:
    return {a["resi"]: np.array([a["x"], a["y"], a["z"]])
           for a in atoms if a["name"] == "CA" and a["chain"] == chain}


def kabsch(P: np.ndarray, Q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rotation+translation mapping P (mobile) onto Q (fixed), least squares."""
    Pc, Qc = P.mean(0), Q.mean(0)
    P0, Q0 = P - Pc, Q - Qc
    H = P0.T @ Q0
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    t = Qc - R @ Pc
    return R, t


def xyz(atoms: list[dict]) -> np.ndarray:
    return np.array([[a["x"], a["y"], a["z"]] for a in atoms], dtype=float)


def ifp_bits(lig_xyz, lig_elem, prot_heavy, prot_xyz) -> set[str]:
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


def _sorted_label_rows(score_label):
    ordered = sorted(score_label, key=lambda t: -t[0])
    return [[int(lbl)] for _, lbl in ordered]


def compute_metrics(score_label):
    rows = _sorted_label_rows(score_label)
    n_pos = sum(r[0] for r in rows)
    out = {"n": len(rows), "n_active": n_pos}
    if n_pos in (0, len(rows)):
        out["error"] = "degenerate labels"
        return out
    out["auc_roc"] = round(float(Scoring.CalcAUC(rows, 0)), 4)
    out["bedroc20"] = round(float(Scoring.CalcBEDROC(rows, 0, BEDROC_ALPHA)), 4)
    for frac, val in zip(EF_FRACTIONS, Scoring.CalcEnrichment(rows, 0, EF_FRACTIONS)):
        out[f"ef_{int(frac*100)}pct"] = round(float(val), 4)
    return out


def bootstrap_ci(score_label, metric, n_boot=BOOTSTRAP_N, seed=BOOTSTRAP_SEED):
    rng = np.random.default_rng(seed)
    arr = np.array(score_label, dtype=object)
    n = len(arr)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        sample = [tuple(arr[i]) for i in idx]
        if len({l for _, l in sample}) < 2:
            continue
        m = compute_metrics(sample)
        if metric in m:
            vals.append(m[metric])
    if not vals:
        return None
    return [round(float(np.percentile(vals, 2.5)), 4), round(float(np.percentile(vals, 97.5)), 4)]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(RUN_DIR))
    args = ap.parse_args()
    outdir = Path(args.out)

    # ---- superpose 5OLH onto 4EIY using matched chain-A CA atoms ----
    ca_4eiy = ca_by_resi(read_pdb(CRYSTAL_4EIY))
    ca_5olh = ca_by_resi(read_pdb(CRYSTAL_5OLH))
    shared = sorted(set(ca_4eiy) & set(ca_5olh))
    P = np.array([ca_5olh[r] for r in shared])   # mobile: 5OLH
    Q = np.array([ca_4eiy[r] for r in shared])    # fixed: 4EIY
    R, t = kabsch(P, Q)
    rmsd = float(np.sqrt(((( R @ P.T).T + t) - Q) ** 2).sum(axis=1).mean() ** 0.5)
    print(f"superposition: {len(shared)} matched C-alpha atoms, RMSD {rmsd:.3f} A "
          f"(prior manuscript reference: 388 atoms, 0.415 A)")

    # ---- 9XT crystal ligand, transformed into 4EIY's frame ----
    xt_atoms = read_pdb(CRYSTAL_5OLH, hetres="9XT")
    xt_xyz_native = xyz(xt_atoms)
    xt_xyz_4eiy_frame = (R @ xt_xyz_native.T).T + t
    xt_elem = [a["elem"] for a in xt_atoms]
    print(f"9XT: {len(xt_atoms)} heavy atoms transformed into 4EIY's frame")

    # ---- ZM241385 crystal ligand (already in 4EIY's own frame) ----
    zma_lines = [l for l in open(CRYSTAL_4EIY)
                if l.startswith("HETATM") and l[17:20].strip() == "ZMA" and l[21] == "A"]
    cryst = Chem.MolFromPDBBlock("".join(zma_lines) + "END\n", removeHs=False, sanitize=False)
    zm_xyz = np.array(cryst.GetConformer().GetPositions())
    zm_elem = [a.GetSymbol() for a in cryst.GetAtoms()]

    prot = read_pdb(RECEPTOR_PDB)
    prot_heavy = [a for a in prot if a["elem"] != "H"]
    prot_xyz = xyz(prot_heavy)

    ref_zm = ifp_bits(zm_xyz, zm_elem, prot_heavy, prot_xyz)
    ref_xt = ifp_bits(xt_xyz_4eiy_frame, xt_elem, prot_heavy, prot_xyz)
    ref_union = ref_zm | ref_xt
    print(f"\nreference IFP bits: ZM241385={len(ref_zm)}  9XT={len(ref_xt)}  "
          f"union={len(ref_union)}  shared={len(ref_zm & ref_xt)}")
    print(f"ZM241385-only bits: {sorted(ref_zm - ref_xt)}")
    print(f"9XT-only bits:      {sorted(ref_xt - ref_zm)}\n")

    # ---- rescore the archived benchmark poses against both references ----
    ledger = pd.read_csv(resolve("20_new_runs/chembl248_redock_corrected_2026-08-10")
                         / "chembl248_ledger.csv")[["compound_id", "label"]].drop_duplicates("compound_id")
    labels = dict(zip(ledger["compound_id"].astype(str), ledger["label"].astype(int)))
    scores = pd.read_csv(under("20_new_runs/chembl245_ifp_rescoring_2026-08-15", "scores.csv"))
    scores = scores[scores["status"] == "ok"]
    poses_dir = under("20_new_runs/chembl245_ifp_rescoring_2026-08-15", "poses")

    rows = []
    for _, r in scores.iterrows():
        cid = str(r["compound_id"])
        if cid not in labels:
            continue
        p = poses_dir / f"{cid}.sdf"
        if not p.is_file():
            continue
        poses = [m for m in Chem.SDMolSupplier(str(p), removeHs=False) if m is not None]
        if not poses:
            continue
        best_zm, best_xt, best_union = 0.0, 0.0, 0.0
        for m in poses:
            pos = np.array(m.GetConformer().GetPositions())
            elem = [a.GetSymbol() for a in m.GetAtoms()]
            fp = ifp_bits(pos, elem, prot_heavy, prot_xyz)
            best_zm = max(best_zm, tanimoto(fp, ref_zm))
            best_xt = max(best_xt, tanimoto(fp, ref_xt))
            best_union = max(best_union, tanimoto(fp, ref_union))
        rows.append(dict(compound_id=cid, label=labels[cid], vina_score=float(r["vina_score"]),
                         ifp_zm_only=round(best_zm, 4), ifp_9xt_only=round(best_xt, 4),
                         ifp_union_ref=round(best_union, 4),
                         ifp_max_of_two=round(max(best_zm, best_xt), 4)))

    df = pd.DataFrame(rows)
    print(f"rescored: {len(df)} compounds\n")

    methods = {
        "ifp_zm_only (baseline, ADR-075)": df["ifp_zm_only"],
        "ifp_9xt_only": df["ifp_9xt_only"],
        "ifp_union_reference": df["ifp_union_ref"],
        "ifp_max_of_two_references": df["ifp_max_of_two"],
    }

    print(f"{'method':>34} {'n':>5} {'AUC':>7} {'95% CI':>16} {'BEDROC20':>9}")
    print("-" * 78)
    summary = {}
    for name, col in methods.items():
        sl = list(zip(col.tolist(), df["label"].tolist()))
        m = compute_metrics(sl)
        if "error" in m:
            print(f"{name:>34}  -- degenerate --")
            summary[name] = m
            continue
        m["auc_ci95"] = bootstrap_ci(sl, "auc_roc")
        summary[name] = m
        ci = m["auc_ci95"]
        ci_str = f"[{ci[0]:.3f}, {ci[1]:.3f}]" if ci else "n/a"
        print(f"{name:>34} {m['n']:>5} {m['auc_roc']:>7.3f} {ci_str:>16} {m['bedroc20']:>9.3f}")

    vina_auc = 0.601
    print(f"\nVina baseline (unchanged): AUC {vina_auc}")
    for name, m in summary.items():
        if "error" in m:
            continue
        print(f"  {name}: delta vs Vina = {m['auc_roc'] - vina_auc:+.3f}")

    df.to_csv(outdir / "multi_reference_ifp_per_compound.csv", index=False)
    (outdir / "multi_reference_ifp_summary.json").write_text(json.dumps(
        {"superposition": {"n_matched_ca": len(shared), "rmsd_A": round(rmsd, 3)},
         "reference_bits": {"zm241385": len(ref_zm), "9xt": len(ref_xt),
                            "union": len(ref_union), "shared": len(ref_zm & ref_xt)},
         "methods": summary}, indent=2, default=str))
    print(f"\nwritten -> {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

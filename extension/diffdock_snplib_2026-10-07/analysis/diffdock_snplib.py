"""ADR-AIIM-317, executed. Predictions frozen before this ran.

Pre-registration SHA-256
  4386e354d4ca7bd1e8a08a36fbfb50ab795f24741853e3e25b66f1957469a265
Decision log at freeze (text/LF)
  3afdea45b0b38d1f4d498025d2a76a01546e86ad2297d34f3dfcaca61b3a3ff6

Pose agreement between DiffDock-L and Vina for the fifteen SNP-Lib survivors. This is a
pose-level check, NOT a ranking channel: DiffDock confidence orders nothing, and the
confidence-versus-score correlation of E4 is a statement about independence only.
Per the diffdock skill, confidence is the model's certainty about a structure and is
not an affinity, and DiffDock is not specifically trained on membrane proteins, which
matters here because the input is the protein alone, without the bilayer.

RMSD is heavy-atom, symmetry-aware and computed IN PLACE (rdMolAlign.CalcRMS, no
alignment), which is correct only because both engines work in the 4EIY frame. That is
checked, not assumed (anchor 1). Nothing is re-docked with Vina and nothing is re-ranked.

Anchors, all required before E1 to E4 may be reported:
  1. the DiffDock input protein's heavy atoms coincide with the Vina receptor's
  2. the Vina centroid distances recomputed from the pose files equal the stored
     centroid_dist_A of ADR-AIIM-313
  3. each Vina pose has RMSD 0 against itself
"""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import rdkit
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolAlign
from scipy.stats import spearmanr

RDLogger.DisableLog("rdApp.*")

ROOT = Path(r"D:\Agentic VS Workflow")
RUNS = ROOT / "20_new_runs"
DD = RUNS / "diffdock_snplib_2026-10-07"
VINA = RUNS / "snplib_dock_2026-10-06"
SNP = RUNS / "snp_screen_2026-08-04" / "G2_preparation"
OUT = DD
PREREG = "4386e354d4ca7bd1e8a08a36fbfb50ab795f24741853e3e25b66f1957469a265"
LOG_AT_FREEZE = "3afdea45b0b38d1f4d498025d2a76a01546e86ad2297d34f3dfcaca61b3a3ff6"

RMSD_AGREE_A = 2.0
SITE_A = 8.0
SHORTLIST = ["SNP00157", "SNP00163", "SNP00154", "SNP00097", "SNP00052",
             "SNP00033", "SNP00142", "SNP00006"]


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def first_mol(path: Path):
    for m in Chem.SDMolSupplier(str(path), removeHs=False):
        if m is not None:
            return m
    return None


def heavy_xyz(m):
    h = Chem.RemoveHs(m)
    return np.array(h.GetConformer().GetPositions()), h


def coords_pdb(path: Path):
    out = []
    for ln in path.read_text(errors="ignore").splitlines():
        if ln.startswith(("ATOM", "HETATM")):
            el = (ln[76:78].strip() or ln[12:16].strip()[:1]).upper()
            if el != "H":
                out.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
    return np.array(out)


def coords_pdbqt(path: Path):
    out = []
    for ln in path.read_text(errors="ignore").splitlines():
        if ln.startswith(("ATOM", "HETATM")):
            at = ln.split()[-1] if len(ln.split()) > 1 else ""
            if at.upper() in ("HD", "H", "HS"):
                continue
            out.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
    return np.array(out)


def main() -> int:
    inputs = {}

    def need(p):
        inputs[str(p.relative_to(ROOT)).replace("\\", "/")] = sha(p)
        return p

    box = json.loads(need(SNP / "docking_box.json").read_text("utf-8"))
    centre = np.array(box["center"], float)                     # ZMA heavy-atom centroid
    scores = pd.read_csv(need(VINA / "scores_with_predictions.csv"))
    scores["SNP_ID"] = scores["SNP_ID"].astype(str)
    surv = pd.read_csv(need(RUNS / "em_survivors_2026-10-06" / "survivors.csv"))
    ids = surv["SNP_ID"].astype(str).tolist()
    assert len(ids) == 15

    # ---- anchor 1: same coordinate frame -----------------------------------
    pdb = coords_pdb(need(DD / "4EIY_protein_only.pdb"))
    qt = coords_pdbqt(need(SNP / "4EIY_prepared.pdbqt"))
    d = np.array([np.min(np.linalg.norm(qt - a, axis=1)) for a in pdb])
    a1 = float(np.percentile(d, 99.5)) < 0.05
    anchor1 = {"n_protein_atoms_pdb": int(len(pdb)), "n_receptor_atoms_pdbqt": int(len(qt)),
               "max_nearest_distance_A": round(float(d.max()), 4),
               "p99.5_nearest_distance_A": round(float(np.percentile(d, 99.5)), 4),
               "verdict": "MET" if a1 else "MISSED"}
    print(f"anchor 1  frame: {len(pdb)} PDB atoms vs {len(qt)} pdbqt atoms, "
          f"p99.5 nearest {anchor1['p99.5_nearest_distance_A']} A -> {anchor1['verdict']}", flush=True)

    rows, problems = [], []
    for sid in ids:
        v = first_mol(need(VINA / "poses" / f"{sid}.sdf"))
        rdir = DD / "results" / sid
        r1 = sorted(rdir.glob("rank1_confidence*.sdf"))
        if v is None or not r1:
            problems.append(sid)
            rows.append({"SNP_ID": sid, "status": "no pose"})
            continue
        f1 = need(r1[0])
        conf = float(re.search(r"confidence(-?[0-9.]+)\.sdf$", f1.name).group(1))
        dm = first_mol(f1)
        vx, vh = heavy_xyz(v)
        dx, dh = heavy_xyz(dm)
        try:
            rmsd = float(rdMolAlign.CalcRMS(dh, vh))
            self0 = float(rdMolAlign.CalcRMS(vh, vh))
            err = ""
        except Exception as e:                                   # noqa: BLE001
            rmsd, self0, err = float("nan"), float("nan"), f"{type(e).__name__}: {e}"[:80]
        # all ten saved poses, as a post hoc descriptive extra
        best = float("nan")
        allr = []
        for f in sorted(rdir.glob("rank*_confidence*.sdf")):
            try:
                m = first_mol(f)
                allr.append(float(rdMolAlign.CalcRMS(Chem.RemoveHs(m), vh)))
            except Exception:                                    # noqa: BLE001
                pass
        if allr:
            best = min(allr)
        rows.append({"SNP_ID": sid, "status": "ok" if not err else "rmsd failed: " + err,
                     "n_heavy_vina": int(vh.GetNumAtoms()), "n_heavy_dd": int(dh.GetNumAtoms()),
                     "dd_confidence": conf,
                     "dd_confidence_band": ("high" if conf > 0 else "moderate" if conf >= -1.5 else "low"),
                     "rmsd_dd_vs_vina_A": rmsd, "vina_self_rmsd_A": self0,
                     "best_of_saved_poses_rmsd_A": best, "n_saved_poses": len(allr),
                     "dd_centroid_dist_A": float(np.linalg.norm(dx.mean(axis=0) - centre)),
                     "vina_centroid_dist_recomputed_A": float(np.linalg.norm(vx.mean(axis=0) - centre))})
    df = pd.DataFrame(rows).merge(
        scores[["SNP_ID", "vina_score", "centroid_dist_A", "on_pocket"]], on="SNP_ID", how="left")
    df["in_shortlist"] = df["SNP_ID"].isin(SHORTLIST)
    df["vina_on_pocket"] = df["on_pocket"].astype(bool)
    df["dd_on_pocket"] = df["dd_centroid_dist_A"] <= SITE_A
    df["agree_within_2A"] = df["rmsd_dd_vs_vina_A"] <= RMSD_AGREE_A

    # ---- anchors 2 and 3 ----------------------------------------------------
    ok = df[df["status"] == "ok"]
    a2 = bool((ok["vina_centroid_dist_recomputed_A"] - ok["centroid_dist_A"]).abs().max() < 0.01)
    a3 = bool((ok["vina_self_rmsd_A"].abs() < 1e-6).all())
    anchors = {"1_frame": anchor1,
               "2_vina_centroid_recomputed_equals_stored": {
                   "max_abs_difference_A": round(float((ok["vina_centroid_dist_recomputed_A"]
                                                        - ok["centroid_dist_A"]).abs().max()), 4),
                   "verdict": "MET" if a2 else "MISSED"},
               "3_vina_self_rmsd_zero": {"verdict": "MET" if a3 else "MISSED"},
               "all_poses_processed": {"n_ok": int(len(ok)), "n_expected": 15, "problems": problems}}
    anchors_ok = a1 and a2 and a3 and len(ok) == 15
    print(f"anchor 2  Vina centroid recomputed vs stored: max diff "
          f"{anchors['2_vina_centroid_recomputed_equals_stored']['max_abs_difference_A']} A -> "
          f"{anchors['2_vina_centroid_recomputed_equals_stored']['verdict']}", flush=True)
    print(f"anchor 3  Vina self-RMSD zero: {'MET' if a3 else 'MISSED'};  processed {len(ok)}/15", flush=True)

    # ---- predictions --------------------------------------------------------
    n_agree = int(ok["agree_within_2A"].sum())
    e1 = n_agree < 8
    off = ok[~ok["vina_on_pocket"]]
    n_off_also = int((off["dd_centroid_dist_A"] > SITE_A).sum())
    e2 = n_off_also >= 5
    r157 = ok[ok["SNP_ID"] == "SNP00157"].iloc[0]
    e3 = bool(r157["rmsd_dd_vs_vina_A"] <= RMSD_AGREE_A)
    rho = float(spearmanr(ok["dd_confidence"], ok["vina_score"]).statistic)
    e4 = abs(rho) < 0.4
    grading = {
        "E1_fewer_than_8_of_15_within_2A": {"n_within_2A": n_agree, "verdict": "MET" if e1 else "MISSED"},
        "E2_at_least_5_of_8_vina_off_pocket_also_off_by_diffdock": {
            "n_vina_off_pocket": int(len(off)), "n_also_off_pocket_by_diffdock": n_off_also,
            "verdict": "MET" if e2 else "MISSED"},
        "E3_SNP00157_within_2A": {"rmsd_A": round(float(r157["rmsd_dd_vs_vina_A"]), 3),
                                  "dd_confidence": float(r157["dd_confidence"]),
                                  "verdict": "MET" if e3 else "MISSED"},
        "E4_abs_spearman_confidence_vs_vina_below_0.4": {"spearman": round(rho, 4), "n": int(len(ok)),
                                                         "verdict": "MET" if e4 else "MISSED"}}
    if not anchors_ok:
        grading = {"NOT_REPORTED": "an anchor failed"}

    df.to_csv(OUT / "pose_comparison.csv", index=False, encoding="utf-8")
    out = {"adr": "ADR-AIIM-317", "prereg_sha256": PREREG, "log_sha256_at_freeze": LOG_AT_FREEZE,
           "rdkit_version": rdkit.__version__,
           "method": "heavy-atom symmetry-aware RMSD computed in place (no alignment) between "
                     "the DiffDock-L rank-1 pose and the Vina rank-1 pose, in the 4EIY frame",
           "thresholds": {"rmsd_agreement_A": RMSD_AGREE_A, "site_A_from_ZMA_centroid": SITE_A},
           "anchors": anchors, "grading": grading,
           "descriptive_not_graded": {
               "dd_confidence_bands": ok["dd_confidence_band"].value_counts().to_dict(),
               "dd_on_pocket_n": int(ok["dd_on_pocket"].sum()),
               "vina_on_pocket_n": int(ok["vina_on_pocket"].sum()),
               "on_off_pocket_classification_agrees_n": int((ok["dd_on_pocket"] == ok["vina_on_pocket"]).sum()),
               "median_rmsd_A": round(float(ok["rmsd_dd_vs_vina_A"].median()), 3),
               "min_max_rmsd_A": [round(float(ok["rmsd_dd_vs_vina_A"].min()), 3),
                                  round(float(ok["rmsd_dd_vs_vina_A"].max()), 3)],
               "best_of_saved_poses_within_2A_n": int((ok["best_of_saved_poses_rmsd_A"] <= RMSD_AGREE_A).sum()),
               "note": "best-of-saved-poses is post hoc and descriptive only"},
           "inputs_sha256": inputs}
    (OUT / "diffdock_snplib.json").write_text(json.dumps(out, indent=1), encoding="utf-8")

    print("\n" + "=" * 90)
    for k, g in grading.items():
        print(f"{k:62s} {g['verdict'] if isinstance(g, dict) else g}")
    print("=" * 90)
    cols = ["SNP_ID", "in_shortlist", "dd_confidence", "dd_confidence_band", "rmsd_dd_vs_vina_A",
            "best_of_saved_poses_rmsd_A", "vina_centroid_dist_recomputed_A", "dd_centroid_dist_A",
            "vina_on_pocket", "dd_on_pocket", "vina_score"]
    print(ok[cols].round(2).to_string(index=False))
    print("\ndescriptive:", json.dumps(out["descriptive_not_graded"]))
    print(f"\nwritten -> {OUT / 'diffdock_snplib.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

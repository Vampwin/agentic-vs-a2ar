"""Grading of ADR-AIIM-355 (pre-registered). Thresholds: agreement RMSD <= 2.0 A; on pocket = centroid within 8.0 A of the
ZM241385 centroid; box = cube of side 22 A on that centroid; G3 >= 14 of 15 survivors with an in-box pose among 100;
G4 >= 13 of 15 on pocket for DiffDock-Pocket; G5, G6 fewer than 8 of 15 within 2.0 A.

  Vina box (existing, ADR-313) : 20_new_runs/snplib_dock_2026-10-06/poses/<ID>.sdf (model 1); control poses_zma/ZM241385.sdf
  arm D  DiffDock-L blind, 100 samples, filtered to the box : diffdock_l_100/results/<ID>/rank*_confidence*.sdf
  arm C  DiffDock-Pocket (flexible, centroid centre) : pocket_C/index*___4EIY_protein_only.pdb___<ID>.sdf/rank1_confidence*.sdf
  variants (descriptive) : pocket_Crigid, pocket_Cca
"""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign

ROOT = Path(r"D:\Agentic VS Workflow")
RUNS = ROOT / "20_new_runs"
OUT = RUNS / "docking_matched_filter_pocket_2026-10-09"
SNP = RUNS / "snp_screen_2026-08-04" / "G2_preparation"
VP = RUNS / "snplib_dock_2026-10-06"
AGREE_A, SITE_A, HALF = 2.0, 8.0, 11.0
SHORTLIST = ["SNP00157", "SNP00163", "SNP00154", "SNP00097", "SNP00052", "SNP00033", "SNP00142", "SNP00006"]
ZMA_SMILES = "Nc1nc(NCCc2ccc(O)cc2)nc2nc(-c3ccco3)nn12"


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def first_heavy(path: Path):
    for m in Chem.SDMolSupplier(str(path), removeHs=False):
        if m is not None:
            return Chem.RemoveHs(m)
    return None


def xyz(m):
    return np.array(m.GetConformer().GetPositions())


def rmsd(a, b):
    try:
        return float(rdMolAlign.CalcRMS(a, b))
    except Exception:                                    # noqa: BLE001
        return float("nan")


def conf_of(p: Path) -> float:
    return float(re.search(r"confidence(-?[0-9.]+)\.sdf$", p.name).group(1))


def pocket_rank1(variant: str, sid: str):
    d = [x for x in (OUT / variant).glob(f"index*___4EIY_protein_only.pdb___{sid}.sdf") if x.is_dir()]
    if not d:
        return None, float("nan")
    f = sorted(d[0].glob("rank1_confidence*.sdf"))
    return (f[0], conf_of(f[0])) if f else (None, float("nan"))


def main():
    have_c = any((OUT / "pocket_C").glob("index*___4EIY_protein_only.pdb___*.sdf"))   # arm C not run or incomplete -> graded as NOT RUN
    hashes = {}

    def need(p):
        hashes[str(Path(p).relative_to(ROOT)).replace("\\", "/")] = sha(Path(p))
        return Path(p)

    box = json.loads(need(SNP / "docking_box.json").read_text("utf-8"))
    centre = np.array(box["center"], float)
    crystal = AllChem.AssignBondOrdersFromTemplate(Chem.MolFromSmiles(ZMA_SMILES), first_heavy(need(SNP / "ZMA_reference_4EIY.sdf")))
    ids = pd.read_csv(need(RUNS / "em_survivors_2026-10-06" / "survivors.csv"))["SNP_ID"].astype(str).tolist()
    assert len(ids) == 15
    order = ["ZM241385"] + ids

    rows, problems = [], []
    for sid in order:
        ctrl = sid == "ZM241385"
        vp = first_heavy(need(VP / ("poses_zma" if ctrl else "poses") / f"{sid}.sdf"))
        # arm D: filter
        poses = sorted((OUT / "diffdock_l_100" / "results" / sid).glob("rank*_confidence*.sdf"),
                       key=lambda p: int(re.match(r"rank(\d+)_", p.name).group(1)))
        n_all, d_mol, d_conf, n_in = len(poses), None, float("nan"), 0
        for p in poses:
            m = first_heavy(p)
            if m is not None and (np.abs(xyz(m) - centre) <= HALF).all():
                n_in += 1
                if d_mol is None:
                    d_mol, d_conf = m, conf_of(p); hashes[str(need(p).relative_to(ROOT)).replace("\\", "/")] = sha(p)
        r = {"compound": sid, "in_shortlist": sid in SHORTLIST, "dd_l_n_samples": n_all, "dd_l_n_in_box": n_in, "dd_l_filtered_confidence": d_conf}
        mols = {"vina_box": vp, "dd_l_filtered": d_mol}
        for v in ("pocket_C", "pocket_Crigid", "pocket_Cca"):
            f, c = pocket_rank1(v, sid)
            mols[v] = first_heavy(need(f)) if f else None
            r[f"{v}_confidence"] = c
        if mols["vina_box"] is None or n_all == 0 or (have_c and mols["pocket_C"] is None):
            problems.append(sid)
        r["status"] = "ok" if sid not in problems else "missing"
        for k, m in mols.items():
            if m is not None:
                r[f"{k}_centroid_dist_A"] = float(np.linalg.norm(xyz(m).mean(0) - centre))
                if ctrl:
                    r[f"{k}_rmsd_to_crystal_A"] = rmsd(m, crystal)
                if k != "vina_box":
                    r[f"{k}_vs_vina_box_rmsd_A"] = rmsd(m, vp)
        rows.append(r)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "filter_pocket_poses.csv", index=False)
    s = df[df["compound"] != "ZM241385"]
    ctrl = df[df["compound"] == "ZM241385"].iloc[0]

    def n_le(col, thr):
        return int((s[col] <= thr).sum()) if col in s else None

    def verdict(cond, val):
        return "UNGRADED (NaN)" if pd.isna(val) else ("MET" if cond else "MISSED")

    anchors = {"arm_C_present": have_c, "all_16_complexes_present_in_arms_C_and_D": {"problems": problems, "verdict": "MET" if not problems else "MISSED"}}
    g = {}
    c1 = ctrl.get("dd_l_filtered_rmsd_to_crystal_A", np.nan)
    g["G1_arm_D_control_within_2A"] = {"rmsd_to_crystal_A": None if pd.isna(c1) else round(float(c1), 3), "verdict": verdict(c1 <= AGREE_A, c1)}
    c2 = ctrl.get("pocket_C_rmsd_to_crystal_A", np.nan)
    g["G2_arm_C_control_within_2A"] = ({"verdict": "NOT RUN"} if not have_c else {"rmsd_to_crystal_A": None if pd.isna(c2) else round(float(c2), 3), "verdict": verdict(c2 <= AGREE_A, c2)})
    n3 = int((s["dd_l_n_in_box"] > 0).sum())
    g["G3_at_least_14_of_15_with_an_in_box_pose"] = {"n": n3, "verdict": "MET" if n3 >= 14 else "MISSED"}
    n4 = n_le("pocket_C_centroid_dist_A", SITE_A)
    g["G4_arm_C_on_pocket_at_least_13"] = {"verdict": "NOT RUN"} if not have_c else {"n_on_pocket": n4, "verdict": "MET" if n4 >= 13 else "MISSED"}
    n5 = n_le("dd_l_filtered_vs_vina_box_rmsd_A", AGREE_A)
    g["G5_box_matched_fewer_than_8_within_2A"] = {"n_within_2A": n5, "n_with_a_pose": int(s["dd_l_filtered_vs_vina_box_rmsd_A"].notna().sum()), "verdict": "MET" if n5 < 8 else "MISSED"}
    n6 = n_le("pocket_C_vs_vina_box_rmsd_A", AGREE_A)
    g["G6_pocket_matched_fewer_than_8_within_2A"] = {"verdict": "NOT RUN"} if not have_c else {"n_within_2A": n6, "verdict": "MET" if n6 < 8 else "MISSED"}

    def med(col):
        return round(float(s[col].median()), 3) if col in s and s[col].notna().any() else None

    desc = {"median_rmsd_to_vina_box_A": {"dd_l_filtered": med("dd_l_filtered_vs_vina_box_rmsd_A"), "pocket_C": med("pocket_C_vs_vina_box_rmsd_A"),
                                          "pocket_Crigid": med("pocket_Crigid_vs_vina_box_rmsd_A"), "pocket_Cca": med("pocket_Cca_vs_vina_box_rmsd_A")},
            "within_2A_of_15": {k: n_le(f"{k}_vs_vina_box_rmsd_A", AGREE_A) for k in ("dd_l_filtered", "pocket_C", "pocket_Crigid", "pocket_Cca")},
            "on_pocket_of_15": {k: n_le(f"{k}_centroid_dist_A", SITE_A) for k in ("vina_box", "dd_l_filtered", "pocket_C", "pocket_Crigid", "pocket_Cca")},
            "dd_l_in_box_poses_of_100": {r["compound"]: r["dd_l_n_in_box"] for r in rows},
            "control_rmsd_to_crystal_A": {k: (None if pd.isna(ctrl.get(f"{k}_rmsd_to_crystal_A", np.nan)) else round(float(ctrl[f"{k}_rmsd_to_crystal_A"]), 3))
                                          for k in ("vina_box", "dd_l_filtered", "pocket_C", "pocket_Crigid", "pocket_Cca")}}
    res = {"adr": "ADR-AIIM-355", "prereg_sha256": (OUT / "PREREG_SHA256.txt").read_text().strip(), "anchors": anchors,
           "grading": g, "descriptive_not_graded": desc, "input_sha256": hashes}
    (OUT / "filter_pocket_result.json").write_text(json.dumps(res, indent=1), "utf-8")
    print(json.dumps({"anchors": anchors, "grading": g, "descriptive": desc}, indent=1))


if __name__ == "__main__":
    main()

"""Grading of ADR-AIIM-353 (pre-registered). Reads four sets of rank-1 poses and writes matched_result.json and
matched_poses.csv in the run folder. Nothing is fitted; thresholds are those of the pre-registration:
agreement RMSD <= 2.0 A, on pocket = centroid within 8.0 A of the ZM241385 centroid, H3 >= 13 of 15, H4/H5 < 8 of 15,
H6 median < 15.3 A.

  Vina pocket   : 20_new_runs/snplib_dock_2026-10-06/poses/<ID>.sdf (model 1); control poses_zma/ZM241385.sdf   [existing, ADR-313]
  DiffDock blind: 20_new_runs/diffdock_snplib_2026-10-07/results/<ID>/rank1*.sdf; control diffdock_orthogonal.../ZM241385_control [existing, ADR-317/327]
  Vina blind    : vina_blind/poses/<ID>.sdf (model 1)                                                          [new, arm A]
  DiffDock pocket: diffdock_pocket/results/<ID>/rank1_confidence*.sdf                                          [new, arm B]
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
OUT = RUNS / "docking_matched_search_space_2026-10-09"
SNP = RUNS / "snp_screen_2026-08-04" / "G2_preparation"
VP = RUNS / "snplib_dock_2026-10-06"
DB = RUNS / "diffdock_snplib_2026-10-07"
DBC = RUNS / "diffdock_orthogonal_2026-08-23" / "results" / "ZM241385_control"
AGREE_A, SITE_A, H3_MIN, H45_MAX, H6_MEDIAN = 2.0, 8.0, 13, 8, 15.3
SHORTLIST = ["SNP00157", "SNP00163", "SNP00154", "SNP00097", "SNP00052", "SNP00033", "SNP00142", "SNP00006"]


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


def dd_rank1(rdir: Path):
    f = sorted(rdir.glob("rank1_confidence*.sdf"))
    if not f:
        return None, float("nan")
    return f[0], float(re.search(r"confidence(-?[0-9.]+)\.sdf$", f[0].name).group(1))


def main():
    hashes = {}

    def need(p):
        hashes[str(p.relative_to(ROOT)).replace("\\", "/")] = sha(p)
        return p

    box = json.loads(need(SNP / "docking_box.json").read_text("utf-8"))
    centre = np.array(box["center"], float)
    # The crystal SDF carries no usable bond orders (its SMILES is saturated), so CalcRMS finds no match against a
    # docked pose. Bond orders are assigned from the ligand's own SMILES; coordinates are untouched.
    crystal = first_heavy(need(SNP / "ZMA_reference_4EIY.sdf"))
    crystal = AllChem.AssignBondOrdersFromTemplate(Chem.MolFromSmiles("Nc1nc(NCCc2ccc(O)cc2)nc2nc(-c3ccco3)nn12"), crystal)
    ids = pd.read_csv(need(RUNS / "em_survivors_2026-10-06" / "survivors.csv"))["SNP_ID"].astype(str).tolist()
    assert len(ids) == 15
    order = ["ZM241385"] + ids

    sets = {}
    for sid in order:
        ctrl = sid == "ZM241385"
        vp = VP / ("poses_zma" if ctrl else "poses") / f"{sid}.sdf"
        db_f, db_c = dd_rank1(DB / "results" / sid) if not ctrl else (DBC / "rank1.sdf", float("nan"))
        if ctrl:
            c1 = sorted(DBC.glob("rank1_confidence*.sdf"))
            db_c = float(re.search(r"confidence(-?[0-9.]+)\.sdf$", c1[0].name).group(1)) if c1 else float("nan")
        vb = OUT / "vina_blind" / "poses" / f"{sid}.sdf"
        dp_f, dp_c = dd_rank1(OUT / "diffdock_pocket" / "results" / sid)
        paths = {"vina_pocket": vp, "dd_blind": db_f, "vina_blind": vb, "dd_pocket": dp_f}
        mols = {k: (first_heavy(need(p)) if p is not None and Path(p).exists() else None) for k, p in paths.items()}
        sets[sid] = (mols, {"dd_blind": db_c, "dd_pocket": dp_c})

    rows, problems = [], []
    for sid in order:
        mols, conf = sets[sid]
        r = {"compound": sid, "in_shortlist": sid in SHORTLIST, "dd_blind_confidence": conf["dd_blind"],
             "dd_pocket_confidence": conf["dd_pocket"]}
        if any(v is None for v in mols.values()):
            problems.append(sid)
            r["status"] = "missing: " + ",".join(k for k, v in mols.items() if v is None)
        else:
            r["status"] = "ok"
        for k, m in mols.items():
            if m is not None:
                r[f"{k}_centroid_dist_A"] = float(np.linalg.norm(xyz(m).mean(0) - centre))
                r[f"{k}_rmsd_to_crystal_A"] = rmsd(m, crystal) if sid == "ZM241385" else float("nan")
        if r["status"] == "ok":
            r["pocket_matched_rmsd_A"] = rmsd(mols["dd_pocket"], mols["vina_pocket"])
            r["blind_matched_rmsd_A"] = rmsd(mols["dd_blind"], mols["vina_blind"])
            r["unmatched_rmsd_A"] = rmsd(mols["dd_blind"], mols["vina_pocket"])
            r["vina_pocket_vs_vina_blind_rmsd_A"] = rmsd(mols["vina_pocket"], mols["vina_blind"])
            r["dd_blind_vs_dd_pocket_rmsd_A"] = rmsd(mols["dd_blind"], mols["dd_pocket"])
        rows.append(r)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "matched_poses.csv", index=False)
    s = df[df["compound"] != "ZM241385"]
    ctrl = df[df["compound"] == "ZM241385"].iloc[0]

    def count_le(col, thr):
        return int((s[col] <= thr).sum())

    # anchors
    qt_ok = int((df["status"] == "ok").sum())
    anchors = {"all_16_poses_present": {"n_ok": qt_ok, "n_expected": 16, "problems": problems,
                                        "verdict": "MET" if qt_ok == 16 else "MISSED"}}
    # unmatched median must reproduce ADR-327's 15.3 (rounded) -- the existing poses are read unchanged
    um = float(s["unmatched_rmsd_A"].median()) if qt_ok == 16 else float("nan")
    anchors["unmatched_median_reproduces_ADR327"] = {"median_A": round(um, 3), "expected": 15.3,
                                                     "verdict": "MET" if abs(um - 15.3) < 0.05 else "MISSED"}

    c_dd = float(ctrl["dd_blind_rmsd_to_crystal_A"])
    anchors["dd_blind_control_reproduces_ADR111_1.917"] = {"rmsd_A": round(c_dd, 3), "verdict": "MET" if abs(c_dd - 1.917) < 0.005 else "MISSED"}
    grading = {}
    grading["H1_vina_blind_control_within_2A"] = {"rmsd_to_crystal_A": round(float(ctrl["vina_blind_rmsd_to_crystal_A"]), 3),
                                                   "verdict": ("UNGRADED (NaN)" if pd.isna(ctrl["vina_blind_rmsd_to_crystal_A"]) else "MET" if ctrl["vina_blind_rmsd_to_crystal_A"] <= AGREE_A else "MISSED")}
    grading["H2_dd_pocket_control_within_2A"] = {"rmsd_to_crystal_A": round(float(ctrl["dd_pocket_rmsd_to_crystal_A"]), 3),
                                                  "verdict": ("UNGRADED (NaN)" if pd.isna(ctrl["dd_pocket_rmsd_to_crystal_A"]) else "MET" if ctrl["dd_pocket_rmsd_to_crystal_A"] <= AGREE_A else "MISSED")}
    n3 = count_le("dd_pocket_centroid_dist_A", SITE_A)
    grading["H3_dd_pocket_on_pocket_at_least_13"] = {"n_on_pocket": n3, "verdict": "MET" if n3 >= H3_MIN else "MISSED"}
    n4 = count_le("pocket_matched_rmsd_A", AGREE_A)
    grading["H4_pocket_matched_agree_fewer_than_8"] = {"n_within_2A": n4, "verdict": "MET" if n4 < H45_MAX else "MISSED"}
    n5 = count_le("blind_matched_rmsd_A", AGREE_A)
    grading["H5_blind_matched_agree_fewer_than_8"] = {"n_within_2A": n5, "verdict": "MET" if n5 < H45_MAX else "MISSED"}
    m6 = float(s["pocket_matched_rmsd_A"].median())
    grading["H6_pocket_matched_median_below_15.3"] = {"median_A": round(m6, 3), "verdict": "MET" if m6 < H6_MEDIAN else "MISSED"}

    desc = {"on_pocket_counts_2x2_of_15": {
                "vina_pocket": count_le("vina_pocket_centroid_dist_A", SITE_A), "vina_blind": count_le("vina_blind_centroid_dist_A", SITE_A),
                "dd_blind": count_le("dd_blind_centroid_dist_A", SITE_A), "dd_pocket": count_le("dd_pocket_centroid_dist_A", SITE_A)},
            "median_rmsd_A": {"pocket_matched": round(m6, 3), "blind_matched": round(float(s["blind_matched_rmsd_A"].median()), 3),
                              "unmatched_adr327": round(um, 3)},
            "within_2A_of_15": {"pocket_matched": n4, "blind_matched": n5, "unmatched": count_le("unmatched_rmsd_A", AGREE_A)},
            "vina_pocket_vs_vina_blind": {"median_rmsd_A": round(float(s["vina_pocket_vs_vina_blind_rmsd_A"].median()), 3),
                                          "within_2A": count_le("vina_pocket_vs_vina_blind_rmsd_A", AGREE_A)},
            "dd_blind_vs_dd_pocket": {"median_rmsd_A": round(float(s["dd_blind_vs_dd_pocket_rmsd_A"].median()), 3),
                                      "within_2A": count_le("dd_blind_vs_dd_pocket_rmsd_A", AGREE_A)},
            "control": {k: (None if pd.isna(ctrl[k]) else round(float(ctrl[k]), 3)) for k in
                        ("vina_pocket_rmsd_to_crystal_A", "vina_blind_rmsd_to_crystal_A", "dd_blind_rmsd_to_crystal_A",
                         "dd_pocket_rmsd_to_crystal_A")}}
    res = {"adr": "ADR-AIIM-353", "prereg_sha256": (OUT / "PREREG_SHA256.txt").read_text().strip(),
           "thresholds": {"agree_A": AGREE_A, "site_A": SITE_A, "H3_min": H3_MIN, "H4_H5_max_exclusive": H45_MAX, "H6_median_A": H6_MEDIAN},
           "anchors": anchors, "grading": grading, "descriptive_not_graded": desc, "input_sha256": hashes}
    (OUT / "matched_result.json").write_text(json.dumps(res, indent=1), "utf-8")
    print(json.dumps({"anchors": anchors, "grading": grading, "descriptive": desc}, indent=1))


if __name__ == "__main__":
    main()

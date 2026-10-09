"""ADR-AIIM-313, executed. Predictions frozen before this ran.

Pre-registration SHA-256
  d161e191a4e193a08d741559cbcbcdbff719819d765c65cff01f7afcf31cca95
Decision log at freeze (text/LF)
  62f8f5d77755418c7297b6e9464aa6bf15589e5cecb8f75467fd5af89946ca95

The SNP-Lib arm is ligand-based only and molecular dynamics consumes a posed
ligand, so the 69 in-domain compounds are docked here. The docking itself is
`10_new_system/scripts/g4_dock_library.py`, the frozen script, INVOKED AT ITS
VALIDATED DEFAULTS AND NOT MODIFIED. This file builds its input, runs the
ZM241385 anchor, computes the predictions the grading needs, and grades.

No number from `exploratory_sbdd_ranking_2026-08-09` is used: that run is
DO-NOT-REUSE (ADR-AIIM-313 constraint 4).

Vina is capped at 8 of 20 threads so ADR-AIIM-305's DiffDock run keeps host
threads; the two do not contend for the GPU.
"""
from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "4")

import numpy as np
import pandas as pd
import rdkit
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor
from sklearn.feature_selection import VarianceThreshold
from sklearn.pipeline import Pipeline

RDLogger.DisableLog("rdApp.*")

ROOT = Path(r"D:\Agentic VS Workflow")
SCRIPTS = ROOT / "10_new_system" / "scripts"
CUR_RUN = ROOT / "20_new_runs" / "chembl251_full_curation_2026-10-05"
REG_RUN = ROOT / "20_new_runs" / "pki_regression_2026-10-05"
EM_RUN = ROOT / "20_new_runs" / "pki_censored_em_2026-10-05"
PC_RUN = ROOT / "20_new_runs" / "pool_comparison_2026-10-06"
SNP_SCREEN = ROOT / "20_new_runs" / "snp_screen_2026-08-04"
SNPLIB = SNP_SCREEN / "G1_standardization" / "snp_source_ledger.csv"
RECEPTOR = SNP_SCREEN / "G2_preparation" / "4EIY_prepared.pdbqt"
BOX = SNP_SCREEN / "G2_preparation" / "docking_box.json"
ZMA_REF = SNP_SCREEN / "G3_redocking" / "ZMA_redocked_4EIY.sdf"
OUTDIR = ROOT / "20_new_runs" / "snplib_dock_2026-10-06"
sys.path.insert(0, str(SCRIPTS))

PREREG = "d161e191a4e193a08d741559cbcbcdbff719819d765c65cff01f7afcf31cca95"
LOG_AT_FREEZE = "62f8f5d77755418c7297b6e9464aa6bf15589e5cecb8f75467fd5af89946ca95"
RS, VINA_CPU, WSL_DISTRO = 42, 8, "Ubuntu-22.04"

ZMA_SMILES = "Nc1nc(NCCc2ccc(O)cc2)nc2nc(-c3ccco3)nn12"
SHORTLIST = ["SNP00157", "SNP00163", "SNP00154", "SNP00097",
             "SNP00052", "SNP00033", "SNP00142", "SNP00006"]
NINE = SHORTLIST + ["SNP00223"]
Z1_RMSD_MAX, Z2_MIN_ON_POCKET, Z3_MIN_SHORTLIST, Z4_MAX_RHO = 2.0, 6, 2, 0.30


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def wslpath(p: Path) -> str:
    s = str(p).replace("\\", "/")
    return "/mnt/" + s[0].lower() + s[2:]


def dock(ledger: Path, out_csv: Path, poses: Path, tag: str):
    """Invoke the FROZEN script in the WSL diffdock env. Nothing is installed."""
    poses.mkdir(parents=True, exist_ok=True)
    cmd = (f"source ~/miniconda3/etc/profile.d/conda.sh && conda activate diffdock && "
           f"python '{wslpath(SCRIPTS / 'g4_dock_library.py')}' "
           f"--ledger '{wslpath(ledger)}' "
           f"--receptor '{wslpath(RECEPTOR)}' "
           f"--box '{wslpath(BOX)}' "
           f"--out '{wslpath(out_csv)}' "
           f"--poses-dir '{wslpath(poses)}' "
           f"--cpu {VINA_CPU}")
    log = OUTDIR / f"dock_{tag}.log"
    with log.open("w", encoding="utf-8") as fh:
        return subprocess.Popen(["wsl", "-d", WSL_DISTRO, "--", "bash", "-lc", cmd],
                                stdout=fh, stderr=subprocess.STDOUT)


def pipe1():
    return Pipeline([("var", VarianceThreshold(threshold=0.01)),
                     ("rf", RandomForestRegressor(random_state=RS, n_jobs=1))])


def best_pose(sdf: Path):
    for m in Chem.SDMolSupplier(str(sdf), removeHs=True):
        if m is not None:
            return m
    return None


def main() -> int:
    (OUTDIR / "analysis").mkdir(parents=True, exist_ok=True)
    qsar = load("_qsar", SCRIPTS / "run_qsar_lbdd_2026-08-07.py")
    cur = load("_cur", CUR_RUN / "analysis" / "chembl251_curate.py")
    rg = load("_reg", REG_RUN / "analysis" / "pki_regression.py")
    em = load("_em", EM_RUN / "analysis" / "pki_censored_em.py")
    pc = load("_pc", PC_RUN / "analysis" / "pool_comparison.py")
    em.pipe = pipe1
    t0 = time.time()

    # ---- the modelling set and the 69-compound in-domain pool ------------
    tbl, _ = rg.rebuild_table(cur)
    b = tbl[tbl["endpoint_class"] == "binding"]
    b = b[b["basis"].isin(["point_median", "point_median_gap",
                           "censored_upper", "censored_lower"])]
    b = b[np.isfinite(b["px"])].reset_index(drop=True)
    kind = np.where(b["basis"] == "censored_upper", 1,
                    np.where(b["basis"] == "censored_lower", -1, 0)).astype(int)
    msm = b["std_smiles"].tolist()
    y0 = b["px"].to_numpy(float)
    iks = [qsar._smiles_to_inchikey(s) or f"__N{i}__" for i, s in enumerate(msm)]
    o = sorted(range(len(msm)), key=lambda j: iks[j])
    msm = [msm[i] for i in o]
    y0, kind = y0[o], kind[o]
    scaf = [qsar._murcko(s) for s in msm]
    X = qsar.compute_raw_features(msm)

    sl_raw = pd.read_csv(SNPLIB)
    sl_raw = sl_raw[sl_raw["smiles_canonical"].notna()].reset_index(drop=True)
    ref = [pc._morgan.GetFingerprint(Chem.MolFromSmiles(s)) for s in msm
           if Chem.MolFromSmiles(s)]
    A, wf = pc.gate(sl_raw, "SNP_ID", "smiles_canonical", ref, qsar, "SNP-Lib")
    A = A[["uid", "_can", "ad_score", "murcko"]].reset_index(drop=True)
    print(f"in-domain pool {len(A)}  [{time.time()-t0:.0f}s]", flush=True)

    ledger = OUTDIR / "dock_input_ledger.csv"
    pd.DataFrame({"SNP_ID": A["uid"], "smiles_canonical": A["_can"]}).to_csv(
        ledger, index=False, encoding="utf-8")
    zma_ledger = OUTDIR / "dock_input_zma.csv"
    pd.DataFrame({"SNP_ID": ["ZM241385"],
                  "smiles_canonical": [ZMA_SMILES]}).to_csv(
        zma_ledger, index=False, encoding="utf-8")

    # ---- docking starts now, in the background, while the model refits ---
    p_lib = dock(ledger, OUTDIR / "scores.csv", OUTDIR / "poses", "library")
    p_zma = dock(zma_ledger, OUTDIR / "zma_scores.csv", OUTDIR / "poses_zma", "zma")
    print("docking launched (frozen script, WSL diffdock env, Vina 1.2.7)",
          flush=True)

    # ---- predictions for all 69, saved this time ------------------------
    cells = [em.cell_name(c) for c in em.CELLS]
    X_A = qsar.compute_raw_features(A["_can"].tolist())
    preds = {}
    for cfg in em.CELLS:
        nm = em.cell_name(cfg)
        tr, _te = em.split(scaf, len(y0), **cfg)
        mdl, n_it, _c, _ = em.censored_em(X, y0, kind, tr)
        preds[nm] = mdl.predict(X_A)
        print(f"  {nm:40s} EM {n_it}it  [{time.time()-t0:.0f}s]", flush=True)
    mean_pki = np.vstack([preds[c] for c in cells]).mean(axis=0)
    pred_df = pd.DataFrame({"SNP_ID": A["uid"], "ad_score": A["ad_score"],
                            "murcko": A["murcko"],
                            "mean_pred_pki": np.round(mean_pki, 4)})
    for c in cells:
        pred_df[f"pred_{c}"] = np.round(preds[c], 4)
    pred_df.to_csv(OUTDIR / "pool69_predictions.csv", index=False,
                   encoding="utf-8")
    print(f"predictions for {len(pred_df)} saved  [{time.time()-t0:.0f}s]",
          flush=True)

    for p, nm in ((p_zma, "zma"), (p_lib, "library")):
        rc = p.wait()
        print(f"docking [{nm}] exit={rc}  [{time.time()-t0:.0f}s]", flush=True)

    # ---- Z1 anchor -------------------------------------------------------
    anchor = {"rmsd": None, "note": ""}
    zp = OUTDIR / "poses_zma" / "ZM241385.sdf"
    if zp.exists() and ZMA_REF.exists():
        a, r = best_pose(zp), best_pose(ZMA_REF)
        if a is not None and r is not None:
            try:
                ref_noh = Chem.RemoveHs(r)
                a_noh = Chem.RemoveHs(a)
                anchor["rmsd"] = round(float(AllChem.GetBestRMS(a_noh, ref_noh)), 4)
            except Exception as e:                       # noqa: BLE001
                anchor["note"] = f"RMSD not computable: {e}"
        else:
            anchor["note"] = "a pose could not be read"
    else:
        anchor["note"] = f"missing {zp if not zp.exists() else ZMA_REF}"
    (OUTDIR / "zma_anchor.json").write_text(json.dumps(anchor, indent=1),
                                            encoding="utf-8")
    print(f"Z1 anchor RMSD {anchor['rmsd']}  {anchor['note']}", flush=True)

    # ---- grading ---------------------------------------------------------
    sc = pd.read_csv(OUTDIR / "scores.csv")
    sc["SNP_ID"] = sc[[c for c in sc.columns
                       if c.lower() in ("compound_id", "snp_id", "id")][0]].astype(str)
    m = sc.merge(pred_df, on="SNP_ID", how="left")
    m.to_csv(OUTDIR / "scores_with_predictions.csv", index=False, encoding="utf-8")

    ok = m[m["status"] == "ok"] if "status" in m.columns else m
    onp = ok[ok["on_pocket"] == True] if "on_pocket" in ok.columns else ok.iloc[0:0]
    nine_on = sorted(set(onp["SNP_ID"]) & set(NINE))
    sl_on = sorted(set(onp["SNP_ID"]) & set(SHORTLIST))
    sub = ok.dropna(subset=["vina_score", "mean_pred_pki"])
    rho = (float(spearmanr(sub["vina_score"], sub["mean_pred_pki"]).statistic)
           if len(sub) > 2 else float("nan"))

    out = {"adr": "ADR-AIIM-313", "prereg_sha256": PREREG,
           "log_sha256_at_freeze": LOG_AT_FREEZE,
           "rdkit_version": rdkit.__version__, "vina": "1.2.7 (WSL diffdock env)",
           "frozen_script": "10_new_system/scripts/g4_dock_library.py, unmodified",
           "receptor": str(RECEPTOR.relative_to(ROOT)).replace("\\", "/"),
           "box": json.loads(BOX.read_text(encoding="utf-8"))["center"],
           "pool_size": int(len(A)), "waterfall": wf,
           "n_docked_ok": int(len(ok)), "n_on_pocket": int(len(onp)),
           "on_pocket_of_nine": nine_on,
           "on_pocket_of_shortlist": sl_on,
           "vina_score_min": round(float(ok["vina_score"].min()), 3) if len(ok) else None,
           "vina_score_median": round(float(ok["vina_score"].median()), 3) if len(ok) else None,
           "anchor": anchor,
           "grading": {
               "Z1_zma_redock_rmsd_at_most_2A": {
                   "rmsd": anchor["rmsd"], "threshold": Z1_RMSD_MAX,
                   "verdict": ("MET" if (anchor["rmsd"] is not None
                                         and anchor["rmsd"] <= Z1_RMSD_MAX)
                               else "MISSED")},
               "Z2_at_least_6_of_nine_on_pocket": {
                   "n": len(nine_on), "threshold": Z2_MIN_ON_POCKET,
                   "which": nine_on,
                   "verdict": "MET" if len(nine_on) >= Z2_MIN_ON_POCKET else "MISSED"},
               "Z3_at_least_2_shortlisted_on_pocket": {
                   "n": len(sl_on), "threshold": Z3_MIN_SHORTLIST, "which": sl_on,
                   "verdict": "MET" if len(sl_on) >= Z3_MIN_SHORTLIST else "MISSED"},
               "Z4_channels_not_redundant": {
                   "spearman_vina_vs_pki": round(rho, 4) if rho == rho else None,
                   "n": int(len(sub)), "threshold_abs": Z4_MAX_RHO,
                   "verdict": ("MET" if (rho == rho and abs(rho) < Z4_MAX_RHO)
                               else "MISSED")}},
           "seconds": round(time.time() - t0, 1)}
    (OUTDIR / "snplib_dock.json").write_text(json.dumps(out, indent=1),
                                             encoding="utf-8")

    print("\n" + "=" * 76)
    for k, g in out["grading"].items():
        print(f"{k:48s} {g['verdict']}")
    print("=" * 76)
    print(f"docked ok {len(ok)}/{len(A)}   on-pocket {len(onp)}")
    print(f"on-pocket among the nine ({len(nine_on)}): {nine_on}")
    print(f"on-pocket among the shortlist ({len(sl_on)}): {sl_on}")
    print(f"Spearman vina vs predicted pKi: {rho:.4f}  (n={len(sub)})")
    if len(onp):
        cols = [c for c in ("SNP_ID", "vina_score", "centroid_dist_A",
                            "mean_pred_pki", "ad_score") if c in onp.columns]
        print("\ntop on-pocket by Vina:")
        print(onp.sort_values("vina_score")[cols].head(12).to_string(index=False))
    print(f"\nwritten -> {OUTDIR / 'snplib_dock.json'}  [{out['seconds']}s]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

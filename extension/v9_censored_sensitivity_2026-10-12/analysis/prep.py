"""ADR-AIIM-412: build the common inputs once (modelling set, features, scaffolds, pool) exactly as em_survivors.py does, and check
the pool against the stored 69 predictions. Writes cache.npz and pool.csv beside this folder. Nothing in 10_new_system is changed."""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

ROOT = Path(r"D:\Agentic VS Workflow")
RUNS = ROOT / "20_new_runs"
SCRIPTS = ROOT / "10_new_system" / "scripts"
OUT = RUNS / "v9_censored_sensitivity_2026-10-12"
sys.path.insert(0, str(SCRIPTS))


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def main():
    t0 = time.time()
    qsar = load("_qsar", SCRIPTS / "run_qsar_lbdd_2026-08-07.py")
    cur = load("_cur", RUNS / "chembl251_full_curation_2026-10-05" / "analysis" / "chembl251_curate.py")
    rg = load("_reg", RUNS / "pki_regression_2026-10-05" / "analysis" / "pki_regression.py")
    pc = load("_pc", RUNS / "pool_comparison_2026-10-06" / "analysis" / "pool_comparison.py")

    tbl, _ = rg.rebuild_table(cur)
    b = tbl[tbl["endpoint_class"] == "binding"]
    b = b[b["basis"].isin(["point_median", "point_median_gap", "censored_upper", "censored_lower"])]
    b = b[np.isfinite(b["px"])].reset_index(drop=True)
    kind = np.where(b["basis"] == "censored_upper", 1, np.where(b["basis"] == "censored_lower", -1, 0)).astype(int)
    msm = b["std_smiles"].tolist()
    y0 = b["px"].to_numpy(float)
    iks = [qsar._smiles_to_inchikey(s) or f"__N{i}__" for i, s in enumerate(msm)]
    o = sorted(range(len(msm)), key=lambda j: iks[j])
    msm = [msm[i] for i in o]
    y0, kind = y0[o], kind[o]
    scaf = [qsar._murcko(s) for s in msm]
    X = qsar.compute_raw_features(msm)
    print(f"modelling set {len(y0)}  X {X.shape} {X.dtype}  point {(kind==0).sum()} up {(kind==1).sum()} lo {(kind==-1).sum()} [{time.time()-t0:.0f}s]", flush=True)

    sl_raw = pd.read_csv(RUNS / "snp_screen_2026-08-04" / "G1_standardization" / "snp_source_ledger.csv")
    sl_raw = sl_raw[sl_raw["smiles_canonical"].notna()].reset_index(drop=True)
    ref = [pc._morgan.GetFingerprint(Chem.MolFromSmiles(s)) for s in msm if Chem.MolFromSmiles(s)]
    A, _wf = pc.gate(sl_raw, "SNP_ID", "smiles_canonical", ref, qsar, "SNP-Lib")
    A = A[["uid", "_can", "ad_score", "murcko", "source"]].reset_index(drop=True)
    X_A = qsar.compute_raw_features(A["_can"].tolist())
    print(f"pool (all SNP-Lib, gated table) {len(A)}  X_A {X_A.shape}  [{time.time()-t0:.0f}s]", flush=True)

    stored = pd.read_csv(RUNS / "snplib_dock_2026-10-06" / "pool69_predictions.csv")
    in_pool = A["uid"].isin(stored["SNP_ID"])
    print("pool uids present in the stored 69:", int(in_pool.sum()), "of", len(stored), flush=True)
    A.to_csv(OUT / "pool_all_gated.csv", index=False)
    np.savez_compressed(OUT / "cache.npz", X=X, y0=y0, kind=kind, X_A=X_A,
                        scaf=np.array(scaf, dtype=object), uid=A["uid"].to_numpy(dtype=object),
                        murcko=A["murcko"].to_numpy(dtype=object))
    (OUT / "prep_summary.json").write_text(json.dumps({
        "n": int(len(y0)), "d_raw": int(X.shape[1]), "point": int((kind == 0).sum()),
        "upper": int((kind == 1).sum()), "lower": int((kind == -1).sum()),
        "pool_rows": int(len(A)), "pool_in_stored69": int(in_pool.sum()), "stored69": int(len(stored)),
        "seconds": round(time.time() - t0, 1)}, indent=1), encoding="utf-8")
    print("done", flush=True)


if __name__ == "__main__":
    main()

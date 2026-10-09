"""ADR-AIIM-305 stage 1, executed. Predictions and the stop rule frozen before
this ran.

Pre-registration SHA-256
  da3da56f9a217940bb010c6e5664245694a2bec5cc67253642010995ebaf7ec2
Decision log at freeze (text/LF, per the 2026-10-05 addendum)
  b2d269a04cd0d4a52194f06e902abc8c6aae8eb8ade96e4ff937a3f0ffc646d1

Does DiffDock confidence rank actives above decoys at all?

The diffdock SKILL.md says three times that it does not predict affinity
("Not for binding affinity prediction", "Confidence != Affinity", and an
explicit "DOES NOT Predict: Binding affinity"). It also lists virtual
screening as an in-scope use. Those are compatible: screening by pose
plausibility is a different and weaker signal than affinity, and whether it
still enriches is the empirical question W2 asks and the documentation does
not answer.

STOP RULE, frozen: if DiffDock's AUC does not exceed Vina's AUC ON THE SAME
593-COMPOUND SUBSAMPLE, stage 2 is not run and no DiffDock channel is added to
anything. Vina at 0.7005 already dilutes the ligand channel (ADR-AIIM-301), so
a structure-based channel weaker than Vina cannot help.

Parameters are ADR-AIIM-111's verbatim -- samples_per_complex 20,
default_inference_args.yaml -- so the two runs are comparable. Halving the cost
by dropping to 10 samples would change the method and is deliberately not done
(ADR-AIIM-305 constraint 3).

A complex yielding no pose is recorded as MISSING, never imputed as a low
score: imputing a bad confidence for a failed generation manufactures a signal
from a crash (constraint 8).

`10_new_system` is not modified. DiffDock is not installed, updated or
reconfigured; the existing env and commit are used as they stand.

Usage:  python diffdock_dude.py --prep      build the input CSV
        python diffdock_dude.py --analyse   parse results and grade
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import re
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "3")

import numpy as np
import pandas as pd
from rdkit import RDLogger
from sklearn.metrics import roc_auc_score

RDLogger.DisableLog("rdApp.*")

ROOT = Path(r"D:\Agentic VS Workflow")
SCRIPTS = ROOT / "10_new_system" / "scripts"
CUR_RUN = ROOT / "20_new_runs" / "chembl251_full_curation_2026-10-05"
REG_RUN = ROOT / "20_new_runs" / "pki_regression_2026-10-05"
DUDE = ROOT / "20_new_runs" / "dude_enrichment_2026-09-22" / "analysis" / "dude_enrichment_per_compound.csv"
PRIOR_DD = ROOT / "20_new_runs" / "diffdock_orthogonal_2026-08-23"
OUTDIR = ROOT / "20_new_runs" / "diffdock_dude_2026-10-05"
sys.path.insert(0, str(SCRIPTS))

PREREG = "da3da56f9a217940bb010c6e5664245694a2bec5cc67253642010995ebaf7ec2"
LOG_AT_FREEZE = "b2d269a04cd0d4a52194f06e902abc8c6aae8eb8ade96e4ff937a3f0ffc646d1"
RS = 42
N_DECOYS = 400              # ADR-AIIM-305: 193 actives + 400 decoys
SAMPLES_PER_COMPLEX = 20    # ADR-AIIM-111 verbatim
RECEPTOR = "4EIY_protein_only.pdb"
VINA_FULL_AUC = 0.7005      # ADR-AIIM-301, full 2,112 set, for context only
W3_FAIL_FLOOR = 0.01


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def build_subsample():
    """DUD-E with the ADR-AIIM-301 leakage control, then 193 actives + 400
    seeded decoys. Returns the frame and the leakage bookkeeping."""
    qsar = load("_qsar", SCRIPTS / "run_qsar_lbdd_2026-08-07.py")
    cur = load("_cur", CUR_RUN / "analysis" / "chembl251_curate.py")
    rg = load("_reg", REG_RUN / "analysis" / "pki_regression.py")

    tbl, _ = rg.rebuild_table(cur)
    b = tbl[tbl["endpoint_class"] == "binding"].reset_index(drop=True)
    keep = b["basis"].isin(["point_median", "point_median_gap",
                            "censored_upper", "censored_lower"])
    b = b[keep & np.isfinite(b["px"])].reset_index(drop=True)
    model_ik = {qsar._smiles_to_inchikey(s) for s in b["std_smiles"]}
    model_ik.discard(None)

    du = pd.read_csv(DUDE)
    du = du[du["SMILES"].notna() & du["vina_score"].notna()].reset_index(drop=True)
    rows_in = len(du)
    du["std_smiles"] = [cur.standardize(s) for s in du["SMILES"]]
    du = du[du["std_smiles"].notna()].reset_index(drop=True)
    du["ik"] = [qsar._smiles_to_inchikey(s) for s in du["std_smiles"]]
    leak = du["ik"].isin(model_ik)
    removed = int(leak.sum())
    du = du[~leak].reset_index(drop=True)

    act = du[du["label"] == 1].reset_index(drop=True)
    dec = du[du["label"] == 0].sample(n=min(N_DECOYS, (du["label"] == 0).sum()),
                                      random_state=RS).reset_index(drop=True)
    sub = pd.concat([act, dec], ignore_index=True)
    meta = {"rows_in": rows_in, "leakage_removed": removed,
            "surviving_full": int(len(du)),
            "surviving_actives": int((du["label"] == 1).sum()),
            "surviving_decoys": int((du["label"] == 0).sum()),
            "subsample_n": int(len(sub)),
            "subsample_actives": int((sub["label"] == 1).sum()),
            "subsample_decoys": int((sub["label"] == 0).sum()),
            "decoy_sample_seed": RS}
    return sub, meta


def prep():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    sub, meta = build_subsample()
    # DiffDock resolves complex_name into a directory name; keep it filesystem
    # safe and unique, and keep the mapping so results can be joined back.
    sub["complex_name"] = [f"DUDE_{i:05d}" for i in range(len(sub))]
    inp = pd.DataFrame({"complex_name": sub["complex_name"],
                        "protein_path": RECEPTOR,
                        "ligand_description": sub["std_smiles"],
                        "protein_sequence": ""})
    inp.to_csv(OUTDIR / "diffdock_input.csv", index=False)
    sub[["complex_name", "compound_id", "label", "vina_score",
         "std_smiles"]].to_csv(OUTDIR / "subsample_map.csv", index=False)
    import shutil
    shutil.copy(PRIOR_DD / RECEPTOR, OUTDIR / RECEPTOR)
    (OUTDIR / "prep_meta.json").write_text(json.dumps(meta, indent=1),
                                           encoding="utf-8")
    print(json.dumps(meta, indent=1))
    print(f"\n{len(inp)} complexes -> {OUTDIR / 'diffdock_input.csv'}")
    print(f"estimated wall time at 65.0 s/complex: "
          f"{len(inp)*65/3600:.1f} h")
    return 0


def parse_results(results_dir: Path):
    """rank1 confidence per complex from the filename DiffDock writes.
    A complex with no rank1_confidence*.sdf is MISSING, not a low score."""
    conf, missing = {}, []
    for d in sorted(results_dir.iterdir()):
        if not d.is_dir():
            continue
        hit = None
        for f in d.iterdir():
            m = re.match(r"rank1_confidence(-?\d+\.?\d*)\.sdf$", f.name)
            if m:
                hit = float(m.group(1))
                break
        if hit is None:
            missing.append(d.name)
        else:
            conf[d.name] = hit
    return conf, missing


def analyse():
    t0 = time.time()
    sub = pd.read_csv(OUTDIR / "subsample_map.csv")
    meta = json.loads((OUTDIR / "prep_meta.json").read_text(encoding="utf-8"))
    results = OUTDIR / "results"
    if not results.exists():
        print(f"no results at {results}; run DiffDock first")
        return 1
    conf, missing = parse_results(results)
    sub["diffdock_conf"] = sub["complex_name"].map(conf)
    n_total = len(sub)
    n_missing = int(sub["diffdock_conf"].isna().sum())
    scored = sub[sub["diffdock_conf"].notna()].reset_index(drop=True)
    y = scored["label"].to_numpy(dtype=int)

    auc_dd = float(roc_auc_score(y, scored["diffdock_conf"].to_numpy(float)))
    # Vina on EXACTLY the same scored compounds, so the comparison is like for like
    auc_vina = float(roc_auc_score(y, -scored["vina_score"].to_numpy(float)))

    stop_passes = auc_dd > auc_vina
    out = {"adr": "ADR-AIIM-305", "stage": 1, "prereg_sha256": PREREG,
           "log_sha256_at_freeze": LOG_AT_FREEZE,
           "diffdock": {"samples_per_complex": SAMPLES_PER_COMPLEX,
                        "config": "default_inference_args.yaml",
                        "receptor": RECEPTOR,
                        "env": "WSL Ubuntu-22.04, conda env diffdock",
                        "repo_commit": os.environ.get("DIFFDOCK_COMMIT", "85c49b6")},
           "subsample": meta,
           "completion": {"complexes_submitted": n_total,
                          "scored": int(len(scored)),
                          "missing": n_missing,
                          "missing_frac": round(n_missing / n_total, 4),
                          "missing_names": missing[:50]},
           "scored_composition": {"actives": int(y.sum()),
                                  "decoys": int((y == 0).sum()),
                                  "active_rate": round(float(y.mean()), 4)},
           "auc": {"diffdock_confidence": round(auc_dd, 4),
                   "vina_same_subsample": round(auc_vina, 4),
                   "vina_full_set_for_context": VINA_FULL_AUC},
           "confidence_distribution": {
               "actives": {"min": round(float(scored.loc[y == 1, "diffdock_conf"].min()), 3),
                           "median": round(float(scored.loc[y == 1, "diffdock_conf"].median()), 3),
                           "max": round(float(scored.loc[y == 1, "diffdock_conf"].max()), 3)},
               "decoys": {"min": round(float(scored.loc[y == 0, "diffdock_conf"].min()), 3),
                          "median": round(float(scored.loc[y == 0, "diffdock_conf"].median()), 3),
                          "max": round(float(scored.loc[y == 0, "diffdock_conf"].max()), 3)}},
           "stop_rule": {
               "condition": "DiffDock AUC must exceed Vina AUC on the same subsample",
               "passes": bool(stop_passes),
               "stage_2_authorised": bool(stop_passes)}}

    out["grading"] = {
        "W1_diffdock_below_vina": {
            "auc_diffdock": round(auc_dd, 4), "auc_vina": round(auc_vina, 4),
            "verdict": "MET" if auc_dd < auc_vina else "MISSED"},
        "W2_diffdock_above_chance": {
            "auc_diffdock": round(auc_dd, 4),
            "verdict": "MET" if auc_dd > 0.50 else "MISSED"},
        "W3_failures_above_1pct": {
            "missing_frac": round(n_missing / n_total, 4),
            "verdict": "MET" if (n_missing / n_total) > W3_FAIL_FLOOR else "MISSED"},
    }
    out["caveats"] = [
        "diffdock SKILL.md states three times that confidence is not an "
        "affinity estimate; this measures whether it enriches regardless, "
        "which the documentation does not address.",
        "DiffDock is 3D-generative and ADR-AIIM-250 measured this project's 3D "
        "stage reproducing 1 of 40 across interpreters; this run does not test "
        "its own reproducibility.",
        "A validated channel would still not license calling anything active "
        "(ADR-AIIM-299: conformal count 0 at pKi >= 6.0 under every framing).",
        "Missing complexes are excluded from AUC, never imputed."]
    out["seconds"] = round(time.time() - t0, 1)
    (OUTDIR / "stage1_auc.json").write_text(json.dumps(out, indent=1),
                                            encoding="utf-8")

    print("=" * 76)
    for t, g in out["grading"].items():
        print(f"{t:34s} {g['verdict']}")
    print("=" * 76)
    print(f"submitted {n_total}  scored {len(scored)}  missing {n_missing} "
          f"({n_missing/n_total*100:.2f} %)")
    print(f"scored composition: {int(y.sum())} active / {int((y==0).sum())} decoy "
          f"({y.mean()*100:.1f} %)")
    print(f"\nAUC  DiffDock confidence {auc_dd:.4f}")
    print(f"     Vina, same subsample {auc_vina:.4f}")
    print(f"     (Vina full set, context only: {VINA_FULL_AUC})")
    cd = out["confidence_distribution"]
    print(f"\nconfidence  actives  min {cd['actives']['min']}  "
          f"median {cd['actives']['median']}  max {cd['actives']['max']}")
    print(f"            decoys   min {cd['decoys']['min']}  "
          f"median {cd['decoys']['median']}  max {cd['decoys']['max']}")
    print(f"\nSTOP RULE: {'PASSES -- stage 2 authorised' if stop_passes else 'FIRES -- stage 2 NOT run, no DiffDock channel added'}")
    print(f"\nwritten -> {OUTDIR / 'stage1_auc.json'}")
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--prep", action="store_true")
    ap.add_argument("--analyse", action="store_true")
    a = ap.parse_args()
    if a.prep:
        raise SystemExit(prep())
    if a.analyse:
        raise SystemExit(analyse())
    ap.print_help()
    raise SystemExit(2)

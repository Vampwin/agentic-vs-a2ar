"""ADR-AIIM-305 stage 1, analysed on the 335 complexes that finished.

Pre-registration SHA-256 of ADR-AIIM-305
  da3da56f9a217940bb010c6e5664245694a2bec5cc67253642010995ebaf7ec2
Decision log at freeze (text/LF)
  b2d269a04cd0d4a52194f06e902abc8c6aae8eb8ade96e4ff937a3f0ffc646d1

WHAT HAPPENED. The 593-complex run stopped at 335 when the host GPU passthrough
faulted (ADR-AIIM-310 / 322). ADR-AIIM-305 does not say what to do on an
interrupted run, and the frozen script `diffdock_dude.py` would count the 258
complexes that never ran as FAILED, which would fabricate a failure rate. This
script is therefore a separate, declared deviation. Its handling was fixed before
any confidence value was read:

  1. Denominator. Attempted complexes are DUDE_00000 to DUDE_00334 (335). The 258
     later complexes were never run and are "not run", not "missing". W3 (failure
     rate) is graded over the 335 attempted only.
  2. W1, W2 graded exactly as frozen (DiffDock AUC vs Vina AUC on the SAME scored
     compounds; DiffDock AUC above 0.50), on the 335. The stop rule is applied as
     written, and reported as resting on a partial subsample.
  3. Composition. The subsample map lists the 193 actives first and the 400 seeded
     random decoys after them, so the 335 are all 193 actives plus the first 142
     decoys. Those 142 are a random subset of the 400 (DataFrame.sample returns
     shuffled rows). The active rate is therefore 57.6 % against a design of
     32.5 %. ROC-AUC does not depend on the base rate in expectation, but the
     decoy side is smaller, so confidence intervals are added (descriptive,
     stratified bootstrap, 2000 resamples, seed 20261008).
  4. Not imputed: a complex with no pose is missing and excluded.

Anchors, all required before grading:
  A1  subsample_map.csv and diffdock_input.csv agree row for row
  A2  the scored complexes are exactly the contiguous block 00000 to 00334, and
      every one has rank1.sdf and ten saved poses
  A3  the heavy-atom count of each DiffDock rank-1 pose equals that of the input
      SMILES (so results are joined to the right compound)

Descriptive, post hoc, not graded: AUC of ligand heavy-atom count alone and its
Spearman correlation with confidence (DUD-E decoys are property-matched, so size
should carry little signal), and Vina's AUC on all 593 for context.
The diffdock skill (v1.1) says confidence is not an affinity; this measures
whether it enriches regardless. `10_new_system` is not modified; nothing re-docked.
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
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

RDLogger.DisableLog("rdApp.*")

OUT = Path(r"D:\Agentic VS Workflow\20_new_runs\diffdock_dude_2026-10-05")
PREREG = "da3da56f9a217940bb010c6e5664245694a2bec5cc67253642010995ebaf7ec2"
LOG_AT_FREEZE = "b2d269a04cd0d4a52194f06e902abc8c6aae8eb8ade96e4ff937a3f0ffc646d1"
VINA_FULL_AUC = 0.7005
W3_FLOOR = 0.01
N_BOOT, SEED = 2000, 20261008


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def rank1_conf(d: Path):
    for f in d.iterdir():
        m = re.match(r"rank1_confidence(-?\d+\.?\d*)\.sdf$", f.name)
        if m:
            return float(m.group(1))
    return None


def heavy_sdf(path: Path):
    for m in Chem.SDMolSupplier(str(path), removeHs=False):
        if m is not None:
            return Chem.RemoveHs(m).GetNumAtoms()
    return None


def boot_ci(y, a, b, rng):
    """95 % percentile CI for AUC(a), AUC(b) and AUC(a) - AUC(b), stratified."""
    pos, neg = np.where(y == 1)[0], np.where(y == 0)[0]
    A, B = [], []
    for _ in range(N_BOOT):
        i = np.concatenate([rng.choice(pos, len(pos)), rng.choice(neg, len(neg))])
        A.append(roc_auc_score(y[i], a[i]))
        B.append(roc_auc_score(y[i], b[i]))
    A, B = np.array(A), np.array(B)
    q = lambda v: [round(float(np.percentile(v, 2.5)), 4), round(float(np.percentile(v, 97.5)), 4)]
    return q(A), q(B), q(A - B)


def main() -> int:
    inputs = {}

    def need(p):
        inputs[p.name] = sha(p)
        return p

    sub = pd.read_csv(need(OUT / "subsample_map.csv"))
    inp = pd.read_csv(need(OUT / "diffdock_input.csv"))
    meta = json.loads(need(OUT / "prep_meta.json").read_text("utf-8"))
    n_total = len(sub)
    assert n_total == 593 == len(inp)

    # ---- A1 ---------------------------------------------------------------
    a1 = bool((sub["complex_name"].to_numpy() == inp["complex_name"].to_numpy()).all()
              and (sub["std_smiles"].to_numpy() == inp["ligand_description"].to_numpy()).all())
    print(f"A1  subsample_map == diffdock_input row for row: {'MET' if a1 else 'MISSED'}", flush=True)

    # ---- attempted vs not run ---------------------------------------------
    res = OUT / "results"
    dirs = {d.name: d for d in res.iterdir() if d.is_dir()}
    conf = {n: rank1_conf(d) for n, d in dirs.items()}
    has = sorted(n for n, c in conf.items() if c is not None)
    last = int(has[-1].split("_")[1])
    attempted = [f"DUDE_{i:05d}" for i in range(last + 1)]
    not_run = [n for n in sub["complex_name"] if n not in set(attempted)]
    assert all(rank1_conf(dirs[n]) is None and not any(dirs[n].iterdir()) for n in not_run), \
        "a later complex has content"
    failed = [n for n in attempted if conf.get(n) is None]

    # ---- A2 ---------------------------------------------------------------
    n_poses = {n: len(list(dirs[n].glob("rank*_confidence*.sdf"))) for n in has}
    a2 = (has == attempted[: len(has)] and len(failed) == 0 and
          all(v == 10 and (dirs[n] / "rank1.sdf").exists() for n, v in n_poses.items()))
    print(f"A2  contiguous block 00000-{last:05d}, rank1.sdf and 10 poses each: "
          f"{'MET' if a2 else 'MISSED'}  (attempted {len(attempted)}, scored {len(has)}, "
          f"failed {len(failed)}, not run {len(not_run)})", flush=True)

    sc = sub[sub["complex_name"].isin(has)].copy().reset_index(drop=True)
    sc["conf"] = sc["complex_name"].map(conf)

    # ---- A3 ---------------------------------------------------------------
    bad = []
    nh_in = {}
    for r in sc.itertuples():
        m = Chem.MolFromSmiles(r.std_smiles)
        nh_in[r.complex_name] = m.GetNumHeavyAtoms() if m is not None else None
        nh_out = heavy_sdf(dirs[r.complex_name] / "rank1.sdf")
        if nh_in[r.complex_name] is None or nh_out != nh_in[r.complex_name]:
            bad.append((r.complex_name, nh_in[r.complex_name], nh_out))
    a3 = not bad
    print(f"A3  rank-1 pose heavy-atom count == input SMILES for all {len(sc)}: "
          f"{'MET' if a3 else 'MISSED ' + str(bad[:5])}", flush=True)
    sc["n_heavy"] = sc["complex_name"].map(nh_in)

    anchors_ok = a1 and a2 and a3
    y = sc["label"].to_numpy(int)
    dd = sc["conf"].to_numpy(float)
    vina = -sc["vina_score"].to_numpy(float)
    auc_dd, auc_v = float(roc_auc_score(y, dd)), float(roc_auc_score(y, vina))
    rng = np.random.default_rng(SEED)
    ci_dd, ci_v, ci_diff = boot_ci(y, dd, vina, rng)

    # descriptive, post hoc
    auc_size = float(roc_auc_score(y, sc["n_heavy"].to_numpy(float)))
    rho_size = float(spearmanr(dd, sc["n_heavy"]).statistic)
    full_v = float(roc_auc_score(sub["label"], -sub["vina_score"]))
    rho_dv = float(spearmanr(dd, sc["vina_score"]).statistic)
    sc_act, sc_dec = sc[sc.label == 1], sc[sc.label == 0]
    q = lambda s: {"min": round(float(s.min()), 3), "median": round(float(s.median()), 3),
                   "max": round(float(s.max()), 3)}
    fail_frac = len(failed) / len(attempted)

    grading = {
        "W1_diffdock_below_vina": {"auc_diffdock": round(auc_dd, 4), "auc_vina": round(auc_v, 4),
                                   "verdict": "MET" if auc_dd < auc_v else "MISSED"},
        "W2_diffdock_above_chance": {"auc_diffdock": round(auc_dd, 4),
                                     "verdict": "MET" if auc_dd > 0.5 else "MISSED"},
        "W3_failures_above_1pct": {"failed": len(failed), "attempted": len(attempted),
                                   "failed_frac": round(fail_frac, 4),
                                   "verdict": "MET" if fail_frac > W3_FLOOR else "MISSED"}}
    if not anchors_ok:
        grading = {"NOT_REPORTED": "an anchor failed"}
    stop_passes = auc_dd > auc_v
    out = {"adr": "ADR-AIIM-305", "stage": 1, "analysis": "partial, 335 of 593 complexes",
           "prereg_sha256": PREREG, "log_sha256_at_freeze": LOG_AT_FREEZE,
           "rdkit_version": rdkit.__version__,
           "completion": {"submitted": n_total, "attempted": len(attempted),
                          "scored": len(sc), "failed_no_pose": len(failed),
                          "not_run": len(not_run), "last_complex": has[-1]},
           "subsample_design": meta,
           "scored_composition": {"actives": int(y.sum()), "decoys": int((y == 0).sum()),
                                  "active_rate": round(float(y.mean()), 4),
                                  "design_active_rate": round(193 / 593, 4),
                                  "decoys_of_400_run": int((y == 0).sum())},
           "anchors": {"A1_map_equals_input": a1, "A2_contiguous_complete": a2,
                       "A3_pose_matches_smiles": a3},
           "auc": {"diffdock_confidence": round(auc_dd, 4), "diffdock_ci95_bootstrap": ci_dd,
                   "vina_same_compounds": round(auc_v, 4), "vina_ci95_bootstrap": ci_v,
                   "difference_dd_minus_vina": round(auc_dd - auc_v, 4),
                   "difference_ci95_bootstrap": ci_diff, "n_boot": N_BOOT, "seed": SEED},
           "grading": grading,
           "stop_rule": {"condition": "DiffDock AUC must exceed Vina AUC on the same subsample",
                         "passes": bool(stop_passes), "stage_2_authorised": bool(stop_passes),
                         "basis": "335 of 593 complexes; a partial subsample"},
           "confidence_distribution": {"actives": q(sc_act["conf"]), "decoys": q(sc_dec["conf"])},
           "descriptive_post_hoc": {
               "auc_heavy_atom_count_alone": round(auc_size, 4),
               "spearman_confidence_vs_heavy_atoms": round(rho_size, 4),
               "spearman_confidence_vs_vina_score": round(rho_dv, 4),
               "vina_auc_all_593_context": round(full_v, 4),
               "vina_auc_full_2112_set_context": VINA_FULL_AUC},
           "inputs_sha256": inputs}
    (OUT / "stage1_auc_partial335.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    sc[["complex_name", "compound_id", "label", "vina_score", "n_heavy", "conf"]].to_csv(
        OUT / "stage1_scored_335.csv", index=False)

    print("\n" + "=" * 78)
    for k, v in grading.items():
        print(f"{k:34s} {v['verdict'] if isinstance(v, dict) else v}")
    print("=" * 78)
    print(f"scored {len(sc)} ({int(y.sum())} active / {int((y == 0).sum())} decoy, "
          f"{y.mean() * 100:.1f}% active; design 32.5%)   failed {len(failed)}   not run {len(not_run)}")
    print(f"AUC  DiffDock {auc_dd:.4f}  95% CI {ci_dd}")
    print(f"     Vina     {auc_v:.4f}  95% CI {ci_v}   (same {len(sc)} compounds)")
    print(f"     DD - Vina {auc_dd - auc_v:+.4f}  95% CI {ci_diff}")
    print(f"context: Vina on all 593 {full_v:.4f}; full-set 0.7005")
    print(f"post hoc: AUC of heavy-atom count alone {auc_size:.4f}; rho(conf, size) {rho_size:+.3f}; "
          f"rho(conf, Vina score) {rho_dv:+.3f}")
    print(f"confidence actives {q(sc_act['conf'])}  decoys {q(sc_dec['conf'])}")
    print(f"STOP RULE: {'PASSES' if stop_passes else 'FIRES: stage 2 not run, no DiffDock channel added'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

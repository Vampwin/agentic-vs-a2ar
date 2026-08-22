"""
run_confidence_tier_positive_control_2026-08-08.py

Confidence-tier framework validation, positive-control arm (ADR-AIIM-051).

The retroactive SNP-Lib check (ADR-AIIM-045) reproduced VERY LOW from first
principles as a *negative* control: a real natural-product library correctly
rated low confidence because 0/236 compounds fall inside the QSAR model's
applicability domain. This script is the paired *positive* control: does the
same framework, same frozen threshold, same fixed model, correctly rate a
small set of natural products HIGHER confidence when there is an a priori,
mechanistic reason to expect them in-domain?

The n=7 set is deliberately small and deliberately NOT presented as a second
screening-scale library (see 08_DECISION_LOG.md ADR-AIIM-051 for the review
of why a "second library" framing does not survive scrutiny at this n).
It is 7 real natural purine/methylxanthine alkaloids and nucleosides --
xanthine, 1-methylxanthine, 3-methylxanthine, theobromine, theophylline,
7-methylpurin-6-amine, deoxyadenosine -- with real literature-sourced
Ki/IC50/EC50 activity against human adenosine A2A receptor from NPASS
(bidd.group/NPASS), confirmed NOT present in the exact ChEMBL rows used to
train the QSAR model (chembl_a2a_raw_v2.csv, ADR-AIIM-047/049).

Reuses the frozen AD method unchanged from run_qsar_lbdd_2026-08-07.py:
ECFP4 (radius 2, 2048 bit), 5-NN mean Tanimoto similarity, 10th-percentile
threshold derived from that run's own training set (0.6811, frozen). This
script recomputes the training-set 5-NN self-similarity threshold from the
already-registered chembl_a2a_training_set.csv and asserts it reproduces the
frozen value -- if it does not, something about the training set or the
scaffold split has drifted, and the run halts rather than silently using a
different threshold.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

_SYSTEM_ROOT = Path(__file__).resolve().parent.parent
if str(_SYSTEM_ROOT) not in sys.path:
    sys.path.insert(0, str(_SYSTEM_ROOT))

import numpy as np
import pandas as pd

from kdss_audit.artifact_registry import ArtifactEntry, ArtifactRegistry, FileAccessEntry, sha256_file
from evaluators._confidence_tier import assess_lbdd
from evaluators._assay_blindness import blocked_columns_in

# Reuse the frozen AD/scaffold functions unchanged rather than reimplementing them.
sys.path.insert(0, str(Path(__file__).resolve().parent))
import importlib.util as _ilu
_spec = _ilu.spec_from_file_location(
    "run_qsar_lbdd_2026_08_07", Path(__file__).resolve().parent / "run_qsar_lbdd_2026-08-07.py"
)
_qsar = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(_qsar)

ROOT = _SYSTEM_ROOT.parent
TRAINING_SET_CSV = ROOT / "20_new_runs" / "qsar_lbdd_2026-08-07" / "chembl_a2a_training_set.csv"
PRIOR_SUMMARY_JSON = ROOT / "20_new_runs" / "qsar_lbdd_2026-08-07" / "qsar_lbdd_summary.json"
OUT_DIR = ROOT / "20_new_runs" / "confidence_tier_positive_control_2026-08-08"
APPROVED_BY = "adisorn.siric@gmail.com"

# ── The n=7 positive-control set ────────────────────────────────────────────
# Source: NPASS (bidd.group/NPASS), target NPT217 = human adenosine A2A
# receptor (UniProt P29274). All 7 confirmed absent from chembl_a2a_raw_v2.csv
# (the exact 8320-row pull used to train the QSAR model, ADR-AIIM-047/049).
# SMILES from NPASS3.0_naturalproducts_structure.txt.
POSITIVE_CONTROL_SET = [
    {"np_id": "NPC287876", "name": "Xanthine", "chembl_id": "CHEMBL1424",
     "smiles": "Oc1nc(O)c2c(n1)[nH]cn2"},
    {"np_id": "NPC278549", "name": "1-Methylxanthine", "chembl_id": "CHEMBL1250",
     "smiles": "Cn1c(O)nc2c(c1=O)nc[nH]2"},
    {"np_id": "NPC41958", "name": "3-Methylxanthine", "chembl_id": "CHEMBL619",
     "smiles": "Cn1c2c(c(nc1=O)O)nc[nH]2"},
    {"np_id": "NPC119133", "name": "Theobromine", "chembl_id": "CHEMBL1114",
     "smiles": "Cn1cnc2c1c(O)nc(=O)n2C"},
    {"np_id": "NPC252603", "name": "Theophylline", "chembl_id": "CHEMBL190",
     "smiles": "Cn1c(=O)n(C)c2c(c1=O)nc[nH]2"},
    {"np_id": "NPC174114", "name": "7-Methylpurin-6-Amine", "chembl_id": "CHEMBL449346",
     "smiles": "Cn1cnc2c1c(N)ncn2"},
    {"np_id": "NPC209525", "name": "Deoxyadenosine", "chembl_id": "CHEMBL449329",
     "smiles": "OC[C@H]1O[C@H](C[C@@H]1O)n1cnc2c1ncnc2N"},
]


def log(msg: str) -> None:
    print(f"[posctrl] {msg}")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    registry = ArtifactRegistry(OUT_DIR / "artifact_registry.jsonl")

    def record_read(path: Path, purpose: str) -> Path:
        registry.record_file_access(FileAccessEntry(reader="run_confidence_tier_positive_control",
                                                      file_path=str(path), purpose=purpose))
        return path

    def register(artifact_id: str, path: Path, description: str, tool: str, decision: str = "") -> None:
        registry.register_artifact(ArtifactEntry(
            artifact_id=artifact_id, path=str(path), sha256=sha256_file(path), stage="G8_reporting",
            description=description, tool=tool, tool_version="", conda_env="", command="",
            random_seed=None, agent_decision_summary=decision,
        ))

    if not TRAINING_SET_CSV.is_file():
        raise FileNotFoundError(f"QSAR training set not found at {TRAINING_SET_CSV} -- run "
                                 "run_qsar_lbdd_2026-08-07.py first.")
    record_read(TRAINING_SET_CSV, "real ChEMBL A2A training set (ADR-AIIM-047), to recompute the "
                                   "frozen AD threshold and 5-NN reference pool")
    combined = pd.read_csv(TRAINING_SET_CSV)
    smiles_all = combined["canonical_smiles"].tolist()
    y_all = combined["label"].to_numpy()
    log(f"loaded {len(combined)} training-set rows (sha256={sha256_file(TRAINING_SET_CSV)[:16]}...)")

    prior = json.loads(PRIOR_SUMMARY_JSON.read_text(encoding="utf-8"))
    frozen_ad_threshold = prior["applicability_domain"]["threshold"]
    frozen_cv_auc = prior["cv_results"]["roc_auc_mean"]
    frozen_external_auc = prior["held_out_test"]["roc_auc"]
    frozen_n_positive = prior["training_set"]["n_active"]

    # ---- Reproduce the same scaffold train/test split, recompute AD threshold
    train_idx, _test_idx = _qsar.scaffold_train_test_split(smiles_all, y_all)
    smiles_train = [smiles_all[i] for i in train_idx]
    recomputed_threshold = _qsar.compute_ad_threshold(smiles_train)
    threshold_matches = abs(recomputed_threshold - frozen_ad_threshold) < 1e-6
    log(f"recomputed AD threshold={recomputed_threshold:.6f} vs frozen={frozen_ad_threshold:.6f} "
        f"-> {'MATCH' if threshold_matches else 'MISMATCH -- HALT'}")
    if not threshold_matches:
        raise RuntimeError(
            "Recomputed AD threshold does not match the frozen value from "
            "qsar_lbdd_summary.json. The training set or scaffold split has drifted since "
            "the Phase A rerun -- refusing to score the positive-control set against a "
            "threshold that is not actually the frozen one."
        )

    # ---- Assay-blindness: this is a control-set CSV we are about to write; make
    # sure it carries no forbidden activity-annotation columns before registering it.
    control_df = pd.DataFrame(POSITIVE_CONTROL_SET)
    blocked = blocked_columns_in(control_df.columns)
    if blocked:
        raise RuntimeError(f"positive-control CSV would carry blocked columns: {blocked}")

    query_smiles = [c["smiles"] for c in POSITIVE_CONTROL_SET]
    ad_scores = _qsar.compute_ad_scores(query_smiles, smiles_train)
    ad_flag_outside = ad_scores < recomputed_threshold
    ad_coverage_fraction = float((~ad_flag_outside).mean())

    control_df["ad_mean_tanimoto"] = np.round(ad_scores, 4)
    control_df["ad_threshold"] = round(recomputed_threshold, 4)
    control_df["outside_domain"] = ad_flag_outside
    control_path = OUT_DIR / "positive_control_ad_scores.csv"
    control_df.to_csv(control_path, index=False)
    register("positive_control_ad_scores", control_path,
              f"AD scores for the n=7 methylxanthine/purine positive-control set. "
              f"Coverage={ad_coverage_fraction:.1%}",
              "rdkit BulkTanimotoSimilarity (frozen method from run_qsar_lbdd_2026-08-07.py)")

    for row in POSITIVE_CONTROL_SET:
        idx = [c["np_id"] for c in POSITIVE_CONTROL_SET].index(row["np_id"])
        log(f"  {row['np_id']} {row['name']:<24} ad_score={ad_scores[idx]:.3f} "
            f"{'OUTSIDE' if ad_flag_outside[idx] else 'inside'} domain")
    log(f"applicability domain: {int((~ad_flag_outside).sum())}/{len(query_smiles)} "
        f"({ad_coverage_fraction:.1%}) of the positive-control set within domain "
        f"[threshold={recomputed_threshold:.4f}, same frozen method as SNP-Lib]")

    # ---- Feed into the same assess_lbdd() used for SNP-Lib. Internal validation,
    # external validation, and imprecision are properties of the fixed QSAR model
    # and are therefore identical to the SNP-Lib run -- only AD coverage can differ.
    lbdd_report = assess_lbdd(
        ad_coverage_fraction=ad_coverage_fraction,
        cv_auc=frozen_cv_auc,
        n_positive=int(frozen_n_positive),
        external_auc=frozen_external_auc,
        external_benchmark_representative=True,
    )
    log(f"confidence tier (LBDD, n=7 positive control): {lbdd_report['final_tier']}")
    for reason in lbdd_report["downgrade_reasons"]:
        log(f"  downgrade: {reason}")

    summary = {
        "run_date": "2026-08-08",
        "purpose": "confidence-tier framework positive control (ADR-AIIM-051), paired with the "
                   "SNP-Lib negative control (ADR-AIIM-045/047)",
        "n_compounds": len(POSITIVE_CONTROL_SET),
        "source": "NPASS (bidd.group/NPASS), target NPT217 = human adenosine A2A receptor "
                  "(UniProt P29274); real literature Ki/IC50/EC50 activity",
        "independence": "all 7 confirmed absent from chembl_a2a_raw_v2.csv, the exact ChEMBL "
                        "pull used to train this QSAR model (ADR-AIIM-047/049)",
        "applicability_domain": {
            "threshold": recomputed_threshold, "threshold_frozen_reference": frozen_ad_threshold,
            "threshold_matches_frozen": threshold_matches,
            "n": len(query_smiles), "coverage_fraction": ad_coverage_fraction,
        },
        "model_properties_reused_unchanged_from_snplib_run": {
            "cv_auc": frozen_cv_auc, "external_auc": frozen_external_auc,
            "n_positive": frozen_n_positive,
        },
        "confidence_tier_lbdd": lbdd_report,
        "per_compound": [
            {**row, "ad_score": round(float(ad_scores[i]), 4), "outside_domain": bool(ad_flag_outside[i])}
            for i, row in enumerate(POSITIVE_CONTROL_SET)
        ],
    }
    (OUT_DIR / "positive_control_summary.json").write_text(
        json.dumps(summary, indent=2, default=str), encoding="utf-8")
    register("positive_control_summary", OUT_DIR / "positive_control_summary.json",
              f"Confidence-tier positive-control result: {lbdd_report['final_tier']}",
              "evaluators._confidence_tier.assess_lbdd")
    log(f"wrote {OUT_DIR / 'positive_control_summary.json'}")
    log("DONE")


if __name__ == "__main__":
    main()

"""
run_qsar_lbdd_2026-08-07.py — QSAR (ligand-based) path, rerun under audit.

Phase A of the QSAR/MD/MM-GBSA rerun decided after the confidence-tier
framework was frozen (ADR-AIIM-046): redo the QSAR/applicability-domain
analysis on real, already-registered data so the LBDD arm of the
confidence-tier assessment has fresh numbers instead of only the retroactive
legacy figures. Phase B (MD/MM-GBSA) is explicitly NOT a rerun -- the
existing 15-system results are audited, not recomputed (human decision,
matches the original project plan's "do not rerun MD" stop rule).

Modelling logic is adapted from the legacy, already leakage-fixed script
(`D:\\kdss_SNP v2\\...\\scripts\\09_chembl_qsar.py`): ECFP4+FCFP4 features,
VarianceThreshold inside the CV pipeline (not before it), scaffold-grouped
StratifiedGroupKFold so no scaffold crosses the train/validation boundary,
isotonic-calibrated RandomForest, and applicability domain as the 10th
percentile of 5-NN Tanimoto similarity *within this run's own training set*
-- not a fixed literature constant. ADR-AIIM-046's G6 contract field
`ad_reference_threshold: 0.681` documented that number as "matching what was
used ad hoc for the 172/172 finding," which was the *value from the legacy
run*, not a portable constant; this run derives its own value from the same
frozen ChEMBL data and reports whatever it actually gets.

Real, already-hashed inputs used (no live ChEMBL API calls):
  - 02_frozen_inputs/public_benchmarks/chembl_a2ar/chembl_a2ar_validation_set.csv
    (3372 rows, G0-registered)
  - 20_new_runs/snp_screen_2026-08-04/G1_standardization/snp_source_ledger.csv
    (236 real SNP-Lib compounds, G1-registered)

Label cutoffs (ACTIVE_CUTOFF=6.0, INACTIVE_CUTOFF=5.5 pChEMBL) are carried
over unchanged from the legacy run rather than re-tuned now, precisely
because changing a threshold after seeing what a run looks like is the one
thing this project's discipline forbids.
"""

from __future__ import annotations

import random
import sys
import warnings
from collections import Counter
from pathlib import Path
from typing import Any

warnings.filterwarnings("ignore", category=UserWarning)

_SYSTEM_ROOT = Path(__file__).resolve().parent.parent
if str(_SYSTEM_ROOT) not in sys.path:
    sys.path.insert(0, str(_SYSTEM_ROOT))

import joblib
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, inchi as rdInchi
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.calibration import CalibratedClassifierCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import VarianceThreshold
from sklearn.metrics import (
    average_precision_score, balanced_accuracy_score, matthews_corrcoef, roc_auc_score,
)
from sklearn.model_selection import RandomizedSearchCV, StratifiedGroupKFold, StratifiedKFold
from sklearn.pipeline import Pipeline

from kdss_audit.artifact_registry import ArtifactEntry, ArtifactRegistry, FileAccessEntry, sha256_file
from evaluators._confidence_tier import assess_lbdd

ROOT = _SYSTEM_ROOT.parent
FROZEN = ROOT / "02_frozen_inputs"
# chembl_a2ar_validation_set.csv (chembl_benchmark_csv) turned out to contain
# only potent/very_potent actives (pChEMBL 7.0-11.0, 0 rows <=5.5) -- a
# curated actives-only reference set, not usable for classifier training.
# The raw ChEMBL pull with a genuine active/inactive spread (pChEMBL
# 4.0-11.0, 843 rows <=5.5) is now a formal, hash-anchored G0 input in its
# own right (chembl_a2ar_raw_csv / hash_chembl_a2ar_raw_csv, ADR-AIIM-049) --
# no longer read from the legacy project location this script originally
# used.
CHEMBL_RAW_SOURCE = FROZEN / "public_benchmarks" / "chembl_a2ar" / "chembl_a2a_raw_v2.csv"
SNPLIB_LEDGER = ROOT / "20_new_runs" / "snp_screen_2026-08-04" / "G1_standardization" / "snp_source_ledger.csv"
OUT_DIR = ROOT / "20_new_runs" / "qsar_lbdd_2026-08-07"
APPROVED_BY = "adisorn.siric@gmail.com"

# Frozen parameters, carried over unchanged from the legacy run.
ACTIVE_CUTOFF = 6.0
INACTIVE_CUTOFF = 5.5
N_OUTER_FOLDS = 5
N_INNER_FOLDS = 3
N_HP_ITER = 20
VAR_THRESHOLD = 0.01
AD_PERCENTILE = 10
AD_K = 5
TEST_FRACTION = 0.20
RANDOM_STATE = 42

HP_SEARCH_SPACE = {
    "rf__n_estimators": [100, 200, 500],
    "rf__max_depth": [None, 15, 25],
    "rf__min_samples_leaf": [1, 5, 10],
    "rf__max_features": ["sqrt", 0.1, 0.2],
}


def log(msg: str) -> None:
    print(f"[qsar] {msg}")


# ── Featurization ───────────────────────────────────────────────────────────

def _ecfp4(smi: str) -> np.ndarray:
    mol = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
    if mol is None:
        return np.zeros(2048, dtype=float)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048, useFeatures=False)
    return np.frombuffer(fp.ToBitString().encode(), dtype="u1") - ord("0")


def _fcfp4(smi: str) -> np.ndarray:
    mol = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
    if mol is None:
        return np.zeros(2048, dtype=float)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048, useFeatures=True)
    return np.frombuffer(fp.ToBitString().encode(), dtype="u1") - ord("0")


def compute_raw_features(smiles_list: list[str]) -> np.ndarray:
    ecfp4 = np.array([_ecfp4(s) for s in smiles_list], dtype=float)
    fcfp4 = np.array([_fcfp4(s) for s in smiles_list], dtype=float)
    return np.concatenate([ecfp4, fcfp4], axis=1)


def _smiles_to_inchikey(smi: str) -> str | None:
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None
        ichi = rdInchi.MolToInchi(mol)
        return rdInchi.InchiToInchiKey(ichi) if ichi else None
    except Exception:
        return None


def _murcko(smi: str) -> str:
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return "__NO_SCAFFOLD__"
        s = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
        return s if s else "__GENERIC__"
    except Exception:
        return "__NO_SCAFFOLD__"


def get_scaffold_groups(smiles: list[str]) -> np.ndarray:
    scaffolds = [_murcko(s) for s in smiles]
    unique = sorted(set(scaffolds))
    sc_map = {s: i for i, s in enumerate(unique)}
    return np.array([sc_map[s] for s in scaffolds])


def scaffold_train_test_split(smiles: list[str], labels: np.ndarray,
                               test_frac: float = TEST_FRACTION,
                               random_state: int = RANDOM_STATE):
    rng = random.Random(random_state)
    scaffolds = [_murcko(s) for s in smiles]
    scaf_to_idx: dict[str, list[int]] = {}
    for i, sc in enumerate(scaffolds):
        scaf_to_idx.setdefault(sc, []).append(i)
    sorted_scafs = sorted(scaf_to_idx.items(), key=lambda x: len(x[1]), reverse=True)
    rng.shuffle(sorted_scafs)
    n_test_target = int(len(labels) * test_frac)
    test_set, train_set = set(), set()
    for sc, idxs in sorted_scafs:
        (test_set if len(test_set) < n_test_target else train_set).update(idxs)
    return np.array(sorted(train_set)), np.array(sorted(test_set))


# ── Applicability domain ────────────────────────────────────────────────────

def _ecfp4_fp(smi: str):
    mol = Chem.MolFromSmiles(smi)
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048) if mol else None


def compute_ad_threshold(train_smiles: list[str], k: int = AD_K) -> float:
    fps = [f for f in (_ecfp4_fp(s) for s in train_smiles) if f is not None]
    sample = random.Random(RANDOM_STATE).sample(fps, min(500, len(fps)))
    nn_scores = []
    for i, fp in enumerate(sample):
        others = [f for j, f in enumerate(fps) if j != i]
        if others:
            sims = DataStructs.BulkTanimotoSimilarity(fp, others)
            nn_scores.append(np.mean(sorted(sims, reverse=True)[:k]))
    arr = np.array(nn_scores)
    threshold = float(np.percentile(arr, AD_PERCENTILE))
    log(f"AD: training self-similarity mean={arr.mean():.3f} SD={arr.std():.3f} "
        f"{AD_PERCENTILE}th-pct={threshold:.3f}")
    return threshold


def compute_ad_scores(query_smiles: list[str], train_smiles: list[str], k: int = AD_K) -> np.ndarray:
    train_fps = [f for f in (_ecfp4_fp(s) for s in train_smiles) if f is not None]
    scores = []
    for smi in query_smiles:
        fp = _ecfp4_fp(smi)
        if fp is None or not train_fps:
            scores.append(0.0)
        else:
            sims = DataStructs.BulkTanimotoSimilarity(fp, train_fps)
            scores.append(float(np.mean(sorted(sims, reverse=True)[:k])))
    return np.array(scores)


# ── Model ────────────────────────────────────────────────────────────────────

def _make_pipe(random_state: int = RANDOM_STATE, **rf_kwargs) -> Pipeline:
    return Pipeline([
        ("var", VarianceThreshold(threshold=VAR_THRESHOLD)),
        ("rf", RandomForestClassifier(class_weight="balanced", random_state=random_state,
                                       n_jobs=-1, **rf_kwargs)),
    ])


def nested_cv_evaluate(X_raw: np.ndarray, y: np.ndarray, scaffold_groups: np.ndarray) -> dict[str, Any]:
    outer_cv = StratifiedGroupKFold(n_splits=N_OUTER_FOLDS, shuffle=True, random_state=RANDOM_STATE)
    inner_cv = StratifiedKFold(n_splits=N_INNER_FOLDS, shuffle=True, random_state=RANDOM_STATE)
    fold_metrics: dict[str, list] = {k: [] for k in ("roc_auc", "mcc", "bal_acc", "best_params")}

    for fold, (tr_idx, val_idx) in enumerate(outer_cv.split(X_raw, y, groups=scaffold_groups)):
        X_tr, X_val = X_raw[tr_idx], X_raw[val_idx]
        y_tr, y_val = y[tr_idx], y[val_idx]
        search = RandomizedSearchCV(_make_pipe(), HP_SEARCH_SPACE, n_iter=N_HP_ITER,
                                     cv=inner_cv, scoring="roc_auc", random_state=RANDOM_STATE,
                                     n_jobs=1, refit=True)
        search.fit(X_tr, y_tr)
        best_pipe = search.best_estimator_
        fold_metrics["best_params"].append(search.best_params_)
        y_prob = best_pipe.predict_proba(X_val)[:, 1]
        y_pred = best_pipe.predict(X_val)
        fold_metrics["roc_auc"].append(roc_auc_score(y_val, y_prob))
        fold_metrics["mcc"].append(matthews_corrcoef(y_val, y_pred))
        fold_metrics["bal_acc"].append(balanced_accuracy_score(y_val, y_pred))
        log(f"fold {fold + 1}/{N_OUTER_FOLDS}: AUC={fold_metrics['roc_auc'][-1]:.3f} "
            f"MCC={fold_metrics['mcc'][-1]:.3f}")

    results: dict[str, Any] = {}
    for key in ("roc_auc", "mcc", "bal_acc"):
        arr = np.array(fold_metrics[key])
        results[f"{key}_mean"] = float(arr.mean())
        results[f"{key}_sd"] = float(arr.std())
    results["best_params"] = fold_metrics["best_params"]
    return results


def get_modal_params(best_params_list: list[dict]) -> dict:
    rf_keys = [k.replace("rf__", "") for k in HP_SEARCH_SPACE]
    modal = {}
    for key in rf_keys:
        pipe_key = f"rf__{key}"
        vals = [str(p.get(pipe_key)) for p in best_params_list]
        modal_str = Counter(vals).most_common(1)[0][0]
        orig = next((p[pipe_key] for p in best_params_list if pipe_key in p), None)
        if modal_str == "None":
            modal[key] = None
        elif isinstance(orig, int):
            modal[key] = int(modal_str)
        elif isinstance(orig, float):
            modal[key] = float(modal_str)
        else:
            modal[key] = modal_str
    return modal


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    registry = ArtifactRegistry(OUT_DIR / "artifact_registry.jsonl")

    def record_read(path: Path, purpose: str) -> Path:
        registry.record_file_access(FileAccessEntry(reader="run_qsar_lbdd", file_path=str(path), purpose=purpose))
        return path

    def register(artifact_id: str, path: Path, description: str, tool: str, decision: str = "") -> None:
        registry.register_artifact(ArtifactEntry(
            artifact_id=artifact_id, path=str(path), sha256=sha256_file(path), stage="G5_LBDD",
            description=description, tool=tool, tool_version="", conda_env="", command="",
            random_seed=None, agent_decision_summary=decision,
        ))

    # ---- 1. Load the real raw ChEMBL data (formal G0 input, ADR-AIIM-049) --
    if not CHEMBL_RAW_SOURCE.is_file():
        raise FileNotFoundError(f"chembl_a2ar_raw_csv G0 input not found at {CHEMBL_RAW_SOURCE}")
    chembl_local = OUT_DIR / "chembl_a2a_raw_actives_and_inactives.csv"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    chembl_local.write_bytes(CHEMBL_RAW_SOURCE.read_bytes())
    record_read(CHEMBL_RAW_SOURCE, "real raw ChEMBL A2A activities (actives+inactives), "
                                    "G0 input chembl_a2ar_raw_csv (ADR-AIIM-049)")
    register("chembl_raw_actives_inactives", chembl_local,
              "Real raw ChEMBL A2A IC50/Ki pull (8320 rows, pChEMBL 4.0-11.0). "
              "Formal G0 input as of ADR-AIIM-049 (hash_chembl_a2ar_raw_csv); distinct from "
              "chembl_benchmark_csv, which is actives-only (0 rows <=5.5 pChEMBL) and cannot "
              "train a classifier.",
              "legacy ChEMBL REST API pull (unchanged, not a live call)",
              "Uses the real active/inactive spread the legacy QSAR script used, now formally "
              "hash-anchored under G0 rather than read ad hoc from the legacy project location.")

    raw = pd.read_csv(chembl_local)
    log(f"loaded {len(raw)} real ChEMBL A2A activity records "
        f"(sha256={sha256_file(chembl_local)[:16]}...)")

    record_read(SNPLIB_LEDGER, "SNP-Lib library for AD/prediction (real, G1-registered)")
    snplib = pd.read_csv(SNPLIB_LEDGER)
    log(f"loaded {len(snplib)} real SNP-Lib compounds "
        f"(sha256={sha256_file(SNPLIB_LEDGER)[:16]}...)")

    # ---- 2. Build training set (frozen cutoffs) ------------------------------
    df = raw.copy()
    act_rel = df["standard_relation"].isin(["=", "<=", "<"])
    inact_rel = df["standard_relation"].isin(["=", ">=", ">"])
    actives = df[act_rel & (df["pchembl_value"] >= ACTIVE_CUTOFF)].copy()
    inactives = df[inact_rel & (df["pchembl_value"] <= INACTIVE_CUTOFF)].copy()
    actives["label"] = 1
    inactives["label"] = 0
    combined = pd.concat([actives, inactives], ignore_index=True)
    combined["inchikey"] = combined["canonical_smiles"].apply(_smiles_to_inchikey)
    combined = combined.dropna(subset=["inchikey", "canonical_smiles"])
    combined = combined.sort_values("pchembl_value", ascending=False).drop_duplicates(subset=["inchikey"]).reset_index(drop=True)
    n_act, n_ina = (combined["label"] == 1).sum(), (combined["label"] == 0).sum()
    log(f"training set: {n_act} actives, {n_ina} inactives (ratio {n_act / max(n_ina, 1):.1f}:1) "
        f"[ACTIVE_CUTOFF={ACTIVE_CUTOFF}, INACTIVE_CUTOFF={INACTIVE_CUTOFF}, frozen from legacy run]")

    train_path = OUT_DIR / "chembl_a2a_training_set.csv"
    combined.to_csv(train_path, index=False)
    register("qsar_training_set", train_path,
              f"ChEMBL A2A training set: {n_act} actives / {n_ina} inactives, deduplicated by InChIKey",
              "pandas + rdkit", "Frozen cutoffs (6.0/5.5 pChEMBL) carried over unchanged from the legacy run")

    smiles_all = combined["canonical_smiles"].tolist()
    y_all = combined["label"].to_numpy()

    # ---- 3. Scaffold-grouped train/test split, then nested CV on train ------
    train_idx, test_idx = scaffold_train_test_split(smiles_all, y_all)
    smiles_train = [smiles_all[i] for i in train_idx]
    smiles_test = [smiles_all[i] for i in test_idx]
    y_train, y_test = y_all[train_idx], y_all[test_idx]
    log(f"scaffold split: {len(train_idx)} train / {len(test_idx)} test "
        f"(active rate train={y_train.mean():.3f} test={y_test.mean():.3f})")

    X_raw_train = compute_raw_features(smiles_train)
    scaffold_groups_train = get_scaffold_groups(smiles_train)

    log(f"running nested CV: {N_OUTER_FOLDS} outer x ({N_HP_ITER} iter x {N_INNER_FOLDS} inner)...")
    cv_results = nested_cv_evaluate(X_raw_train, y_train, scaffold_groups_train)
    log(f"CV ROC-AUC = {cv_results['roc_auc_mean']:.3f} +/- {cv_results['roc_auc_sd']:.3f}  "
        f"MCC = {cv_results['mcc_mean']:.3f}")

    # ---- 4. Final model on full training set, evaluate on held-out test -----
    modal_params = get_modal_params(cv_results["best_params"])
    base_rf = RandomForestClassifier(**modal_params, class_weight="balanced",
                                      random_state=RANDOM_STATE, n_jobs=-1)
    cal = CalibratedClassifierCV(base_rf, method="isotonic", cv=5)
    final_pipe = Pipeline([("var", VarianceThreshold(threshold=VAR_THRESHOLD)), ("cal", cal)])
    final_pipe.fit(X_raw_train, y_train)

    X_raw_test = compute_raw_features(smiles_test)
    y_prob_test = final_pipe.predict_proba(X_raw_test)[:, 1]
    test_auc = float(roc_auc_score(y_test, y_prob_test))
    test_pr_inactive = float(average_precision_score(1 - y_test, 1 - y_prob_test))
    log(f"held-out scaffold-split test: AUC={test_auc:.3f} PR-AUC(inactive)={test_pr_inactive:.3f} "
        f"n={len(y_test)} ({int((y_test == 0).sum())} inactives)")

    model_path = OUT_DIR / "qsar_pipeline_a2a.pkl"
    joblib.dump(final_pipe, model_path)
    register("qsar_calibrated_pipeline", model_path,
              f"Calibrated RF QSAR pipeline. CV AUC={cv_results['roc_auc_mean']:.3f}, "
              f"held-out test AUC={test_auc:.3f}",
              "scikit-learn RandomForestClassifier + CalibratedClassifierCV(isotonic)")

    # ---- 5. Applicability domain: derive from training set, apply to real SNP-Lib
    ad_threshold = compute_ad_threshold(smiles_train)
    smi_col = "smiles_canonical"
    snp_smiles = [s for s in snplib[smi_col].tolist() if isinstance(s, str) and s]
    ad_scores = compute_ad_scores(snp_smiles, smiles_train)
    ad_flag_outside = ad_scores < ad_threshold
    ad_coverage_fraction = float((~ad_flag_outside).mean())
    log(f"applicability domain: {int((~ad_flag_outside).sum())}/{len(snp_smiles)} "
        f"({ad_coverage_fraction:.1%}) of the real SNP-Lib library within domain "
        f"[threshold={ad_threshold:.3f}, derived from this run's own training set]")

    p_active = final_pipe.predict_proba(compute_raw_features(snp_smiles))[:, 1]
    pred_df = pd.DataFrame({
        "SNP_ID": snplib.loc[snplib[smi_col].isin(snp_smiles), "SNP_ID"].values[: len(snp_smiles)],
        "smiles": snp_smiles,
        "p_active_calibrated": np.round(p_active, 4),
        "ad_mean_tanimoto": np.round(ad_scores, 4),
        "ad_threshold": round(ad_threshold, 4),
        "outside_domain": ad_flag_outside,
    })
    pred_path = OUT_DIR / "snplib_qsar_predictions.csv"
    pred_df.to_csv(pred_path, index=False)
    register("qsar_snplib_predictions", pred_path,
              f"QSAR P(active) and AD status for {len(snp_smiles)} real SNP-Lib compounds. "
              f"AD coverage={ad_coverage_fraction:.1%}",
              "scikit-learn + rdkit BulkTanimotoSimilarity")

    # ---- 6. Feed real numbers into the confidence-tier LBDD assessment ------
    lbdd_report = assess_lbdd(
        ad_coverage_fraction=ad_coverage_fraction,
        cv_auc=cv_results["roc_auc_mean"],
        n_positive=int(n_act),
        external_auc=test_auc,
        external_benchmark_representative=True,  # scaffold-held-out split of the same ChEMBL set
    )
    log(f"confidence tier (LBDD, real 2026-08-07 rerun): {lbdd_report['final_tier']}")
    for reason in lbdd_report["downgrade_reasons"]:
        log(f"  downgrade: {reason}")

    summary = {
        "run_date": "2026-08-07",
        "training_set": {"n_active": int(n_act), "n_inactive": int(n_ina),
                          "active_cutoff": ACTIVE_CUTOFF, "inactive_cutoff": INACTIVE_CUTOFF},
        "cv_results": {"roc_auc_mean": cv_results["roc_auc_mean"], "roc_auc_sd": cv_results["roc_auc_sd"],
                       "mcc_mean": cv_results["mcc_mean"], "bal_acc_mean": cv_results["bal_acc_mean"]},
        "held_out_test": {"n": len(y_test), "roc_auc": test_auc, "pr_auc_inactive": test_pr_inactive},
        "applicability_domain": {"threshold": ad_threshold, "percentile": AD_PERCENTILE, "k": AD_K,
                                  "snplib_n": len(snp_smiles), "snplib_coverage_fraction": ad_coverage_fraction},
        "confidence_tier_lbdd": lbdd_report,
    }
    import json
    (OUT_DIR / "qsar_lbdd_summary.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    log(f"wrote {OUT_DIR / 'qsar_lbdd_summary.json'}")
    log("DONE")


if __name__ == "__main__":
    main()

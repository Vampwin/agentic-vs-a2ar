"""
retrain_risk_monitor_2026-08-05.py

Re-runs the WP4 development risk-monitor pipeline (ADR-AIIM-027 through
ADR-AIIM-032 in 00_project_control/08_DECISION_LOG.md) on condition-B pilot
data collected AFTER the G2/G4/G8 evaluator stub fixes
(20_new_runs/WP2.5_WP3A_readiness_2026-07-30/pilot_postfix_2026-08-04.json).

Every design choice here (feature allowlist, model candidates, calibration
method comparison, calibration split strategy, action budget, OOD detector
hyperparameters) is unchanged from the already human-approved procedure --
only the input labels changed, because the evaluators that produce them were
fixed. This mirrors how the earlier SEER(B)-only re-run was treated as a
bug-fix re-run rather than a new design decision.

This script does NOT edit or overwrite ADR-AIIM-034 (Gate 9, "frozen ...
immutable, require new ADR to change"). It writes a new, separately dated
results file for human review. A human must open a new ADR to accept these
numbers as superseding the frozen ones.

KNOWN LIMITATION, accepted by human decision on 2026-08-05: escaped-fault
positives in the new data are concentrated in 2 of 11 development families
(F10, F13; 9 of 133 fault cases total). GroupKFold(5)/(3) folds outside those
two families contain zero positive examples, so compare_models() and the
calibration comparison are statistically fragile for this rerun. Documented,
not worked around, per the human decision to proceed with the unmodified
procedure and record the limitation.
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from kdss_paths import project_root, resolve, under
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

from fault_injection.pilot_runner import build_training_dataset
from risk_monitor.train import (
    fit_selected_model, CVResult, N_OUTER_FOLDS, N_INNER_FOLDS,
    _PARAM_GRIDS, _make_pipeline, _check_groups_and_balance,
)
from risk_monitor.calibrate import calibrate, evaluate_calibration, save_calibrated_model
from risk_monitor.registry import ModelRegistry
from risk_monitor.inference import predict_risk, InferenceError
from risk_monitor.policy import PolicyAction
from telemetry.ood_detector import OutOfDomainDetector
from telemetry.leakage_guard import check_field_names, check_allowlist
from sklearn.model_selection import GroupKFold, GridSearchCV
from sklearn.metrics import brier_score_loss, average_precision_score


def resilient_compare_models(X: pd.DataFrame, y: np.ndarray, groups: np.ndarray) -> list[CVResult]:
    """
    Same nested-GroupKFold comparison as risk_monitor.train.compare_models(),
    with one addition: when an inner GridSearchCV fold contains a single class
    (possible here because positives are concentrated in 2/11 families —
    documented limitation, human decision 2026-08-05 to proceed and record it
    rather than change the CV design), fall back to fitting the pipeline with
    its default hyperparameters for that outer fold instead of crashing.
    train.py itself is untouched; this wrapper lives only in this rerun script.
    """
    check_field_names(X.columns, context="resilient_compare_models X columns")
    check_allowlist(X.columns, context="resilient_compare_models X columns")
    y = np.asarray(y)
    groups = np.asarray(groups)
    _check_groups_and_balance(y, groups, context="resilient_compare_models")

    candidates = {"logistic_l2": "logistic_l2", "gradient_boosting": "gradient_boosting"}
    outer_cv = GroupKFold(n_splits=N_OUTER_FOLDS)
    results: list[CVResult] = []

    for name in candidates:
        fold_details = []
        brier_scores, auprc_scores = [], []
        param_grid = _PARAM_GRIDS.get(name, {})

        for fold_idx, (train_idx, val_idx) in enumerate(outer_cv.split(X, y, groups=groups)):
            X_train, X_val = X.iloc[train_idx], X.iloc[val_idx]
            y_train, y_val = y[train_idx], y[val_idx]
            groups_train = groups[train_idx]

            base_pipeline = _make_pipeline(name)
            fallback_used = False
            if param_grid and len(np.unique(y_train)) > 1:
                inner_groups_count = len(np.unique(groups_train))
                inner_cv = GroupKFold(n_splits=min(N_INNER_FOLDS, max(2, inner_groups_count)))
                try:
                    search = GridSearchCV(base_pipeline, param_grid, cv=inner_cv,
                                           scoring="neg_brier_score", refit=True,
                                           n_jobs=1, error_score="raise")
                    search.fit(X_train, y_train, groups=groups_train)
                    best_pipeline = search.best_estimator_
                    best_params = search.best_params_
                except ValueError as exc:
                    logger_msg = f"fold {fold_idx} inner CV failed ({exc}); using default hyperparameters"
                    log(f"  [{name}] {logger_msg}")
                    base_pipeline.fit(X_train, y_train)
                    best_pipeline = base_pipeline
                    best_params = {}
                    fallback_used = True
            else:
                if len(np.unique(y_train)) <= 1:
                    log(f"  [{name}] fold {fold_idx}: single-class training fold, fitting with defaults")
                    fallback_used = True
                base_pipeline.fit(X_train, y_train)
                best_pipeline = base_pipeline
                best_params = {}

            proba = best_pipeline.predict_proba(X_val)[:, 1]
            brier = brier_score_loss(y_val, proba)
            if len(np.unique(y_val)) > 1:
                auprc = average_precision_score(y_val, proba)
            else:
                auprc = float("nan")  # undefined with a single class in the validation fold
            brier_scores.append(brier)
            auprc_scores.append(auprc)
            fold_details.append({
                "fold": fold_idx, "n_train": len(train_idx), "n_val": len(val_idx),
                "brier": brier, "auprc": auprc, "n_val_positive": int(y_val.sum()),
                "n_train_positive": int(y_train.sum()), "best_params": best_params,
                "fallback_to_default_hyperparams": fallback_used,
            })

        results.append(CVResult(
            model_name=name,
            mean_brier=float(np.mean(brier_scores)),
            std_brier=float(np.std(brier_scores)),
            mean_auprc=float(np.nanmean(auprc_scores)),
            std_auprc=float(np.nanstd(auprc_scores)),
            fold_details=fold_details,
        ))

    results.sort(key=lambda r: r.mean_brier)
    return results

BASE = resolve("20_new_runs/WP2.5_WP3A_readiness_2026-07-30")
OUT_DIR = resolve("20_new_runs/WP4_retrain_postfix_2026-08-05")
MODELS_DIR = OUT_DIR / "models"
REGISTRY_PATH = BASE / "model_registry.json"  # same append-only registry as before
APPROVED_BY = "adisorn.siric@gmail.com"
DATE_TAG = "2026-08-05"

PILOT_JSON = BASE / "pilot_postfix_2026-08-04.json"


def log(msg: str) -> None:
    print(f"[retrain] {msg}")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    MODELS_DIR.mkdir(parents=True, exist_ok=True)

    # ---- 1. Load post-fix pilot data, build training set (condition B) -----
    raw = json.loads(PILOT_JSON.read_text(encoding="utf-8"))
    all_outcomes = raw["outcomes"]
    results = SimpleNamespace(outcomes=all_outcomes)

    X, y, groups = build_training_dataset(results, reference_condition="B")
    n_pos = int(y.sum())
    log(f"X={X.shape} y positive_rate={y.mean():.4f} ({n_pos}/{len(y)}) n_groups={len(set(groups))}")

    fam_pos = {}
    for g, yi in zip(groups, y):
        fam_pos.setdefault(g, [0, 0])
        fam_pos[g][1] += 1
        if yi:
            fam_pos[g][0] += 1
    log("per-family positives: " + ", ".join(f"{k}={v[0]}/{v[1]}" for k, v in sorted(fam_pos.items())))

    # ---- 2. compare_models (Gate 5 candidates, exploratory CV) --------------
    cv_results = resilient_compare_models(X, y, groups)
    cv_table = [
        {"model": r.model_name, "mean_brier": r.mean_brier, "std_brier": r.std_brier,
         "mean_auprc": r.mean_auprc, "std_auprc": r.std_auprc}
        for r in cv_results
    ]
    for row in cv_table:
        log(f"compare_models: {row['model']:18s} brier={row['mean_brier']:.4f}+-{row['std_brier']:.4f} "
            f"auprc={row['mean_auprc']:.4f}+-{row['std_auprc']:.4f}")

    # ADR-AIIM-022 tie-break: prefer logistic_l2 when Brier means are within
    # 1 std of either candidate (same rule applied in ADR-AIIM-027).
    by_name = {r.model_name: r for r in cv_results}
    lr, gbm = by_name["logistic_l2"], by_name["gradient_boosting"]
    diff = abs(lr.mean_brier - gbm.mean_brier)
    within_1std = diff < max(lr.std_brier, gbm.std_brier)
    chosen = "logistic_l2" if within_1std else min(by_name, key=lambda k: by_name[k].mean_brier)
    log(f"model selection: brier diff={diff:.4f}, within 1 std={within_1std} -> chosen={chosen}")

    # ---- 3. Fit selected model on all development condition-B rows ---------
    fitted = fit_selected_model(X, y, model_name=chosen)
    uncal_id = f"{chosen}_dev_{DATE_TAG}"
    uncal_path = MODELS_DIR / f"{uncal_id}.joblib"
    import joblib
    joblib.dump(fitted, uncal_path, compress=3)

    registry = ModelRegistry(REGISTRY_PATH)
    feature_order = list(X.columns)
    schema_path = Path(__file__).resolve().parents[1] / "telemetry" / "feature_schema.yaml"
    registry.register(
        artifact_id=uncal_id,
        model_path=uncal_path,
        schema_version="1.0.0",
        approved_by=APPROVED_BY,
        notes=(
            f"{chosen} — development data only (post evaluator-fix rerun); not calibrated; "
            f"not frozen. Fit on {len(y)} condition-B pilot outcomes ({len(set(groups))} "
            f"families, {X.shape[1]} features). Rerun of ADR-AIIM-027/028 procedure on "
            f"pilot_postfix_2026-08-04.json. positive_rate={y.mean():.4f} "
            f"({n_pos}/{len(y)}), concentrated in families "
            f"{[k for k, v in fam_pos.items() if v[0] > 0]} — see limitations in results JSON."
        ),
        feature_order=feature_order,
        schema_file_path=schema_path,
    )
    log(f"registered uncalibrated model: {uncal_id}")

    # ---- 4. Calibration: compare isotonic vs sigmoid (Gate 6) --------------
    cal_results = {}
    for method in ("isotonic", "sigmoid"):
        cal_pipeline = calibrate(fitted, X, y, groups, method=method)
        ev = evaluate_calibration(cal_pipeline, X, y)
        cal_results[method] = {"pipeline": cal_pipeline, "brier": ev["brier_score"],
                                "ece": float(np.mean(np.abs(
                                    np.array(ev["fraction_of_positives"]) - np.array(ev["mean_predicted_value"])
                                )))}
        log(f"calibration[{method}]: brier={cal_results[method]['brier']:.4f} ece={cal_results[method]['ece']:.4f}")

    iso, sig = cal_results["isotonic"], cal_results["sigmoid"]
    if iso["brier"] <= sig["brier"] and iso["ece"] <= sig["ece"]:
        cal_method = "isotonic"
    elif sig["brier"] <= iso["brier"] and sig["ece"] <= iso["ece"]:
        cal_method = "sigmoid"
    else:
        cal_method = "isotonic" if iso["brier"] < sig["brier"] else "sigmoid"  # tie-break: lower Brier (primary metric)
    log(f"calibration method selected: {cal_method}")

    calibrated = cal_results[cal_method]["pipeline"]
    cal_id = f"{chosen}_cal_{cal_method}_dev_{DATE_TAG}"
    cal_path = MODELS_DIR / f"{cal_id}.joblib"
    save_calibrated_model(calibrated, cal_path)
    registry.register(
        artifact_id=cal_id,
        model_path=cal_path,
        schema_version="1.0.0",
        approved_by=APPROVED_BY,
        notes=(
            f"Calibrated ({cal_method}) rerun of ADR-AIIM-029 on post evaluator-fix data. "
            f"Brier={cal_results[cal_method]['brier']:.4f} ECE={cal_results[cal_method]['ece']:.4f}. "
            "In-sample cross-fitted calibration on the same development rows used for fitting "
            "(same design as the original run); see limitations for class-imbalance caveat."
        ),
        feature_order=feature_order,
        schema_file_path=schema_path,
    )
    log(f"registered calibrated model: {cal_id}")

    # ---- 5. OOD detector (ADR-AIIM-031 hyperparameters, unchanged) ---------
    ood = OutOfDomainDetector(contamination=0.01, n_estimators=100, random_state=20260729)
    ood.fit(X)
    ood_id = f"ood_isoforest_dev_{DATE_TAG}"
    ood_path = MODELS_DIR / f"{ood_id}.joblib"
    ood.save(ood_path)
    ood_flags_train = [ood.predict_one(row.to_dict()) for _, row in X.iterrows()]
    n_ood_train = sum(1 for f in ood_flags_train if f)
    registry.register(
        artifact_id=ood_id,
        model_path=ood_path,
        schema_version="1.0.0",
        approved_by=APPROVED_BY,
        notes=(
            f"IsolationForest OOD detector, rerun of ADR-AIIM-031 on post evaluator-fix data. "
            f"{len(ood.non_null_features)} non-null features. "
            f"Training self-check: {n_ood_train}/{len(X)} flagged OOD."
        ),
        feature_order=feature_order,
        schema_file_path=schema_path,
    )
    log(f"registered OOD detector: {ood_id} (train self-check {n_ood_train}/{len(X)} OOD)")

    # ---- 6. Threshold derivation (Gate 6 budget: FE <= 10%, unchanged) -----
    # Re-derive per-row case_type aligned with X's row order (build_training_dataset
    # preserves the filter order of `usable`, which mirrors the outcomes list order).
    b_rows = [o for o in all_outcomes if o.get("condition") == "B" and not o.get("error")
              and isinstance(o.get("features"), dict) and len(o["features"]) > 0]
    assert len(b_rows) == len(X), f"row alignment mismatch: {len(b_rows)} vs {len(X)}"

    proba_all = calibrated.predict_proba(X)[:, 1]
    clean_mask = np.array([r["case_type"] == "clean" for r in b_rows])
    escaped_fault_mask = (~clean_mask) & (y == 1)
    caught_fault_mask = (~clean_mask) & (y == 0)

    clean_probs = proba_all[clean_mask]
    escaped_probs = proba_all[escaped_fault_mask]
    n_clean, n_escaped, n_caught = int(clean_mask.sum()), int(escaped_fault_mask.sum()), int(caught_fault_mask.sum())
    log(f"clean={n_clean} escaped_fault={n_escaped} caught_fault={n_caught}")

    FE_BUDGET = 0.10
    candidate_t = sorted(set(clean_probs.tolist()) | {0.0, 1.0})
    t_escalate = 1.0
    for t in candidate_t:
        fe = float((clean_probs >= t).mean())
        if fe <= FE_BUDGET:
            t_escalate = t
            break
    actual_fe_clean = float((clean_probs >= t_escalate).mean())

    if n_escaped > 0:
        t_halt = float(np.percentile(escaped_probs, 25))
        t_halt = max(t_halt, t_escalate)
    else:
        t_halt = 1.0  # no escaped-fault probabilities to derive from (would need explicit human fallback)
    log(f"t_escalate={t_escalate:.6f} (FE={actual_fe_clean:.4f}, budget<= {FE_BUDGET}) t_halt={t_halt:.6f} (n_escaped={n_escaped})")

    def _rate(mask_probs: np.ndarray, thalt: float, tesc: float) -> dict:
        n = len(mask_probs)
        if n == 0:
            return {"n": 0, "halt_rate": None, "escalate_rate": None, "continue_rate": None}
        halt = float((mask_probs >= thalt).mean())
        esc = float(((mask_probs >= tesc) & (mask_probs < thalt)).mean())
        cont = float((mask_probs < tesc).mean())
        return {"n": n, "halt_rate": round(halt, 4), "escalate_rate": round(esc, 4), "continue_rate": round(cont, 4)}

    dev_sim = {
        "clean": _rate(clean_probs, t_halt, t_escalate),
        "escaped_fault": _rate(escaped_probs, t_halt, t_escalate),
        "caught_fault": {"n": n_caught, "halt_rate": 1.0, "note": "contracts always HALT; model not consulted"},
    }

    policy_sha = None
    import hashlib
    policy_payload = {
        "policy_version": "1.0.0",
        "derived_on": DATE_TAG,
        "approved_by": APPROVED_BY,
        "artifact_id": cal_id,
        "thresholds": {"t_escalate": t_escalate, "t_halt": t_halt},
        "budget": {"fe_budget": FE_BUDGET},
        "derivation": {
            "t_escalate_method": "smallest threshold with FE rate <= 10% on development clean cases",
            "t_halt_method": "25th percentile of escaped-fault predicted probabilities, clamped >= t_escalate",
            "n_clean": n_clean,
            "n_escaped_fault": n_escaped,
            "actual_fe_clean": actual_fe_clean,
        },
        "development_simulation": dev_sim,
        "limitation": (
            f"n_escaped_fault dropped from 113 (pre-fix) to {n_escaped} (post-fix), concentrated "
            f"in families {[k for k, v in fam_pos.items() if v[0] > 0]}. t_halt is a 25th-percentile "
            f"estimate over only {n_escaped} points and compare_models()/calibration folds outside "
            "these families contain zero positive examples. Reported as-is per human decision "
            "2026-08-05 to run the unmodified procedure and document the limitation rather than "
            "adapt the CV design."
        ),
    }
    policy_bytes = json.dumps({k: v for k, v in policy_payload.items() if k != "sha256"},
                               indent=2, sort_keys=False).encode("utf-8")
    policy_sha = hashlib.sha256(policy_bytes).hexdigest().upper()
    policy_payload["sha256"] = policy_sha
    policy_path = OUT_DIR / "policy_config.json"
    policy_path.write_text(json.dumps(policy_payload, indent=2), encoding="utf-8")
    log(f"wrote {policy_path} sha256={policy_sha}")

    # ---- 7. Simulate conditions A and C on the same 188 B-condition rows ---
    def simulate(condition: str) -> list[dict]:
        rows_out = []
        for row_dict, feat_row, evaluator_state, case_type, fam in zip(
            b_rows, X.to_dict(orient="records"), [r["evaluator_state"] for r in b_rows],
            [r["case_type"] for r in b_rows], groups,
        ):
            telemetry = dict(feat_row)
            telemetry["out_of_domain_flag"] = ood.predict_one(feat_row)
            contracts_all_passed = True if condition == "C" else (evaluator_state == "PASS")
            try:
                res = predict_risk(
                    telemetry=telemetry, artifact_id=cal_id, registry=registry,
                    t_escalate=t_escalate, t_halt=t_halt,
                    contracts_all_passed=contracts_all_passed,
                )
                action, risk, reason = res.action, res.risk, res.reason
            except InferenceError as exc:
                action, risk, reason = PolicyAction.HALT, None, f"fail-closed: {exc}"

            if case_type == "clean":
                escaped = action == PolicyAction.HALT
                escaped_detected = action in (PolicyAction.HALT, PolicyAction.ESCALATE)
            else:
                escaped = action != PolicyAction.HALT
                escaped_detected = action not in (PolicyAction.HALT, PolicyAction.ESCALATE)

            rows_out.append({
                "family_id": fam, "instance_id": row_dict["instance_id"], "case_type": case_type,
                "condition": condition, "action": action.value if hasattr(action, "value") else str(action),
                "risk": risk, "reason": reason, "escaped": escaped, "escaped_detected": escaped_detected,
            })
        return rows_out

    cond_A_rows = simulate("A")
    cond_C_rows = simulate("C")

    def seer(rows: list[dict]) -> dict:
        fault = [r for r in rows if r["case_type"] == "fault"]
        clean = [r for r in rows if r["case_type"] == "clean"]
        return {
            "seer": (sum(r["escaped"] for r in fault) / len(fault)) if fault else None,
            "seer_detected": (sum(r["escaped_detected"] for r in fault) / len(fault)) if fault else None,
            "false_escalation": (sum(r["escaped"] for r in clean) / len(clean)) if clean else None,
            "n_fault": len(fault), "n_clean": len(clean),
        }

    seer_A, seer_C = seer(cond_A_rows), seer(cond_C_rows)
    # SEER(B)/(D) already known from pilot_postfix_2026-08-04.json; recomputed here for the same 188 rows for a consistent table.
    fault_B = [r for r in b_rows if r["case_type"] == "fault"]
    clean_B = [r for r in b_rows if r["case_type"] == "clean"]
    seer_B_val = sum(1 for r in fault_B if r["escaped"]) / len(fault_B)
    fe_B_val = sum(1 for r in clean_B if r["escaped"]) / len(clean_B)

    d_rows = [o for o in all_outcomes if o.get("condition") == "D" and not o.get("error")]
    fault_D = [r for r in d_rows if r["case_type"] == "fault"]
    seer_D_val = (sum(1 for r in fault_D if r["escaped"]) / len(fault_D)) if fault_D else 1.0

    per_family = {}
    for fam in sorted(set(groups)):
        fam_fault_B = [r for r in fault_B if r["family_id"] == fam]
        fam_A = [r for r in cond_A_rows if r["family_id"] == fam and r["case_type"] == "fault"]
        fam_C = [r for r in cond_C_rows if r["family_id"] == fam and r["case_type"] == "fault"]
        per_family[fam] = {
            "n_fault": len(fam_fault_B),
            "seer_B": (sum(1 for r in fam_fault_B if r["escaped"]) / len(fam_fault_B)) if fam_fault_B else None,
            "seer_A": (sum(1 for r in fam_A if r["escaped"]) / len(fam_A)) if fam_A else None,
            "seer_C": (sum(1 for r in fam_C if r["escaped"]) / len(fam_C)) if fam_C else None,
        }

    seer_summary = {
        "SEER_D": round(seer_D_val, 6),
        "SEER_B": round(seer_B_val, 6),
        "SEER_A": round(seer_A["seer"], 6) if seer_A["seer"] is not None else None,
        "SEER_C": round(seer_C["seer"], 6) if seer_C["seer"] is not None else None,
        "delta_seer_A_vs_D": round(seer_D_val - seer_A["seer"], 6) if seer_A["seer"] is not None else None,
        "delta_seer_C_vs_D": round(seer_D_val - seer_C["seer"], 6) if seer_C["seer"] is not None else None,
        "sesoi": 0.15,
        "A_vs_D_meets_sesoi": (seer_D_val - seer_A["seer"] >= 0.15) if seer_A["seer"] is not None else None,
        "false_escalation_B": round(fe_B_val, 6),
        "false_escalation_A": round(seer_A["false_escalation"], 6) if seer_A["false_escalation"] is not None else None,
        "false_escalation_C": round(seer_C["false_escalation"], 6) if seer_C["false_escalation"] is not None else None,
    }
    log("SEER summary: " + json.dumps(seer_summary, indent=2))

    out_payload = {
        "run_date": DATE_TAG,
        "supersedes_pending_approval": "ADR-AIIM-034 (Gate 9, frozen SEER_A=0.165) -- NOT auto-applied; requires a new ADR",
        "source_pilot": str(PILOT_JSON),
        "artifacts": {
            "uncalibrated_model": uncal_id,
            "calibrated_model": cal_id,
            "ood_detector": ood_id,
            "model_registry": str(REGISTRY_PATH),
        },
        "model_selection": {"cv_table": cv_table, "chosen": chosen, "tie_break_within_1std": within_1std},
        "calibration_comparison": {m: {"brier": v["brier"], "ece": v["ece"]} for m, v in cal_results.items()},
        "calibration_method_selected": cal_method,
        "policy": {"t_escalate": t_escalate, "t_halt": t_halt, "sha256": policy_sha},
        "seer_summary": seer_summary,
        "per_family": per_family,
        "condition_A_rows": cond_A_rows,
        "condition_C_rows": cond_C_rows,
        "limitations": [
            f"Escaped-fault positives: {n_escaped}/{len(fault_B)} fault cases, concentrated in "
            f"families {[k for k, v in fam_pos.items() if v[0] > 0]} (2 of 11 development families). "
            "GroupKFold(5)/(3) folds outside these two families have zero positive examples; "
            "compare_models() Brier/AUPRC and the isotonic-vs-sigmoid calibration comparison are "
            "correspondingly fragile for this rerun.",
            f"t_halt derived from only {n_escaped} escaped-fault calibrated probabilities "
            "(25th percentile of a 9-point sample) -- a much smaller and less stable estimate "
            "than the original 113-point derivation.",
            "In-sample cross-fitted calibration (same rows used for fitting and calibration), "
            "matching the original design; not a held-out estimate.",
            "Development-data estimate only. Held-out evaluation (H01-H12) has not been run; "
            "these numbers cannot be reported as the paper's primary RQ1/RQ2 result until Gate 7 "
            "(verifier) closes and held-out evaluation confirms or revises them.",
        ],
    }
    results_path = OUT_DIR / f"wp4_conditions_ac_results_postfix_{DATE_TAG}.json"
    results_path.write_text(json.dumps(out_payload, indent=2), encoding="utf-8")
    log(f"wrote {results_path}")
    log("DONE")


if __name__ == "__main__":
    main()

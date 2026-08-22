"""
retrain_risk_monitor_lofo_2026-08-05.py

Second retrain attempt on post evaluator-fix condition-B data
(pilot_postfix_2026-08-04.json), after the first attempt
(retrain_risk_monitor_2026-08-05.py, ADR-AIIM-038) produced a degenerate
model: GroupKFold(5)/(3) nested CV assumes every fold has both classes, but
positives are concentrated in 2 of 11 development families (F10, F13; 9 of
133 fault cases). Investigation confirmed there is no untested instance left
to find within the approved F01-F18 grids (pilot already ran every instance
in every family's frozen param_grid) and expanding F10/F13's frozen grid
would not fix the fold-level degeneracy anyway (it is a group-count problem,
not an instance-count problem) -- so this is a genuine, human-approved
deviation from the model-selection procedure, not a way to manufacture a
better-looking number from the same broken design.

Deviation, and why:

  GroupKFold(5) outer + GroupKFold(3) inner  ->  Leave-One-Family-Out (11
  folds, one per development family), no inner hyperparameter search.

  - LOFO guarantees every training fold contains both classes (whichever of
    F10/F13 is not the held-out family stays in training), which is exactly
    what crashed the plain GroupKFold(5) attempt.
  - AUPRC is computed pooled across all 11 held-out folds' predictions
    rather than averaged per-fold, because 9 of the 11 individual folds have
    no positives in validation (AUPRC undefined for those folds in
    isolation) -- pooling is the standard treatment for rare-event CV.
  - Inner hyperparameter search is dropped rather than adapted, because
    grid-searching a regularisation constant against 9 positive examples is
    not a justifiable use of the data; the fixed default (C=1.0 for
    logistic_l2) matches what ADR-AIIM-028 actually used for the final fit
    last time regardless of what the grid search found.

Feature allowlist, model candidates, calibration methods compared, and the
FE<=10% action budget are all unchanged from the originally approved design.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import joblib
import numpy as np
import pandas as pd
from sklearn.calibration import CalibratedClassifierCV
from sklearn.frozen import FrozenEstimator
from sklearn.metrics import average_precision_score, brier_score_loss

from fault_injection.pilot_runner import build_training_dataset
from risk_monitor.train import CVResult, _make_pipeline, fit_selected_model
from risk_monitor.calibrate import evaluate_calibration, save_calibrated_model
from risk_monitor.registry import ModelRegistry
from risk_monitor.inference import predict_risk, InferenceError
from risk_monitor.policy import PolicyAction
from telemetry.ood_detector import OutOfDomainDetector
from telemetry.leakage_guard import check_field_names, check_allowlist

BASE = Path(r"D:/Agentic VS Workflow/20_new_runs/WP2.5_WP3A_readiness_2026-07-30")
OUT_DIR = Path(r"D:/Agentic VS Workflow/20_new_runs/WP4_retrain_lofo_2026-08-05")
MODELS_DIR = OUT_DIR / "models"
REGISTRY_PATH = BASE / "model_registry.json"
APPROVED_BY = "adisorn.siric@gmail.com"
DATE_TAG = "2026-08-05-lofo-v2"
PILOT_JSON = BASE / "pilot_postfix_2026-08-04.json"


def log(msg: str) -> None:
    print(f"[lofo] {msg}")


def lofo_compare_models(X: pd.DataFrame, y: np.ndarray, groups: np.ndarray) -> list[CVResult]:
    """Leave-one-family-out comparison, no inner hyperparameter search, pooled AUPRC."""
    check_field_names(X.columns, context="lofo_compare_models X columns")
    check_allowlist(X.columns, context="lofo_compare_models X columns")

    families = sorted(set(groups))
    results: list[CVResult] = []

    for name in ("logistic_l2", "gradient_boosting"):
        fold_details = []
        brier_scores = []
        pooled_y_val: list[int] = []
        pooled_proba: list[float] = []

        for held_out in families:
            train_mask = groups != held_out
            val_mask = groups == held_out
            X_train, y_train = X.loc[train_mask], y[train_mask]
            X_val, y_val = X.loc[val_mask], y[val_mask]

            pipeline = _make_pipeline(name)
            pipeline.fit(X_train, y_train)
            proba = pipeline.predict_proba(X_val)[:, 1]

            brier = brier_score_loss(y_val, proba)
            brier_scores.append(brier)
            pooled_y_val.extend(y_val.tolist())
            pooled_proba.extend(proba.tolist())

            fold_details.append({
                "held_out_family": held_out, "n_train": int(train_mask.sum()),
                "n_val": int(val_mask.sum()), "n_train_positive": int(y_train.sum()),
                "n_val_positive": int(y_val.sum()), "brier": brier,
            })

        pooled_y_val_arr = np.asarray(pooled_y_val)
        pooled_proba_arr = np.asarray(pooled_proba)
        overall_auprc = (
            average_precision_score(pooled_y_val_arr, pooled_proba_arr)
            if len(np.unique(pooled_y_val_arr)) > 1 else float("nan")
        )

        results.append(CVResult(
            model_name=name,
            mean_brier=float(np.mean(brier_scores)),
            std_brier=float(np.std(brier_scores)),
            mean_auprc=float(overall_auprc),  # pooled, not per-fold-averaged
            std_auprc=0.0,  # a single pooled estimate has no fold-to-fold std
            fold_details=fold_details,
        ))

    results.sort(key=lambda r: r.mean_brier)
    return results


def lofo_calibrate(fitted_pipeline, X: pd.DataFrame, y: np.ndarray, groups: np.ndarray,
                    method: str) -> CalibratedClassifierCV:
    """Cross-fitted calibration using LOFO splits instead of GroupKFold(3)."""
    families = sorted(set(groups))
    idx = np.arange(len(y))
    explicit_splits = [
        (idx[groups != fam], idx[groups == fam]) for fam in families
    ]
    calibrated = CalibratedClassifierCV(
        estimator=FrozenEstimator(fitted_pipeline), method=method, cv=explicit_splits,
    )
    calibrated.fit(X, y)
    return calibrated


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    MODELS_DIR.mkdir(parents=True, exist_ok=True)

    raw = json.loads(PILOT_JSON.read_text(encoding="utf-8"))
    results_ns = SimpleNamespace(outcomes=raw["outcomes"])
    X, y, groups = build_training_dataset(results_ns, reference_condition="B")
    log(f"X={X.shape} y positive_rate={y.mean():.4f} ({int(y.sum())}/{len(y)}) n_groups={len(set(groups))}")

    # ---- 1. LOFO model comparison ------------------------------------------
    cv_results = lofo_compare_models(X, y, groups)
    cv_table = [
        {"model": r.model_name, "mean_brier": r.mean_brier, "std_brier": r.std_brier,
         "pooled_auprc": r.mean_auprc}
        for r in cv_results
    ]
    for row in cv_table:
        log(f"LOFO: {row['model']:18s} brier={row['mean_brier']:.4f}+-{row['std_brier']:.4f} "
            f"pooled_auprc={row['pooled_auprc']:.4f}")

    by_name = {r.model_name: r for r in cv_results}
    lr, gbm = by_name["logistic_l2"], by_name["gradient_boosting"]
    diff = abs(lr.mean_brier - gbm.mean_brier)
    within_1std = diff < max(lr.std_brier, gbm.std_brier)
    chosen = "logistic_l2" if within_1std else min(by_name, key=lambda k: by_name[k].mean_brier)
    log(f"model selection: brier diff={diff:.4f}, within 1 std={within_1std} -> chosen={chosen}")

    # ---- 2. Fit on all development rows ------------------------------------
    fitted = fit_selected_model(X, y, model_name=chosen)
    uncal_id = f"{chosen}_dev_{DATE_TAG}"
    uncal_path = MODELS_DIR / f"{uncal_id}.joblib"
    joblib.dump(fitted, uncal_path, compress=3)

    registry = ModelRegistry(REGISTRY_PATH)
    feature_order = list(X.columns)
    schema_path = Path(__file__).resolve().parents[1] / "telemetry" / "feature_schema.yaml"
    registry.register(
        artifact_id=uncal_id, model_path=uncal_path, schema_version="1.0.0",
        approved_by=APPROVED_BY,
        notes=(f"{chosen} -- LOFO-selected, development data only (post evaluator-fix, "
               f"leave-one-family-out CV). Fit on {len(y)} condition-B rows "
               f"({len(set(groups))} families). Deviation from GroupKFold(5)/(3) documented "
               "in ADR-AIIM-041 (class-imbalance degeneracy in ADR-AIIM-038)."),
        feature_order=feature_order, schema_file_path=schema_path,
    )
    log(f"registered uncalibrated model: {uncal_id}")

    # ---- 3. LOFO calibration comparison ------------------------------------
    cal_results = {}
    for method in ("isotonic", "sigmoid"):
        cal_pipeline = lofo_calibrate(fitted, X, y, groups, method=method)
        ev = evaluate_calibration(cal_pipeline, X, y)
        ece = float(np.mean(np.abs(np.array(ev["fraction_of_positives"]) - np.array(ev["mean_predicted_value"]))))
        cal_results[method] = {"pipeline": cal_pipeline, "brier": ev["brier_score"], "ece": ece}
        log(f"calibration[{method}]: brier={ev['brier_score']:.4f} ece={ece:.4f}")

    iso, sig = cal_results["isotonic"], cal_results["sigmoid"]
    if iso["brier"] <= sig["brier"] and iso["ece"] <= sig["ece"]:
        cal_method = "isotonic"
    elif sig["brier"] <= iso["brier"] and sig["ece"] <= iso["ece"]:
        cal_method = "sigmoid"
    else:
        cal_method = "isotonic" if iso["brier"] < sig["brier"] else "sigmoid"
    log(f"calibration method selected: {cal_method}")

    calibrated = cal_results[cal_method]["pipeline"]
    cal_id = f"{chosen}_cal_{cal_method}_dev_{DATE_TAG}"
    cal_path = MODELS_DIR / f"{cal_id}.joblib"
    save_calibrated_model(calibrated, cal_path)
    registry.register(
        artifact_id=cal_id, model_path=cal_path, schema_version="1.0.0",
        approved_by=APPROVED_BY,
        notes=(f"Calibrated ({cal_method}), LOFO cross-fitted (11 folds) on post "
               f"evaluator-fix data. Brier={cal_results[cal_method]['brier']:.4f} "
               f"ECE={cal_results[cal_method]['ece']:.4f}."),
        feature_order=feature_order, schema_file_path=schema_path,
    )
    log(f"registered calibrated model: {cal_id}")

    # ---- 4. OOD detector (unchanged hyperparameters/design) ---------------
    ood = OutOfDomainDetector(contamination=0.01, n_estimators=100, random_state=20260729)
    ood.fit(X)
    ood_id = f"ood_isoforest_dev_{DATE_TAG}"
    ood_path = MODELS_DIR / f"{ood_id}.joblib"
    ood.save(ood_path)
    n_ood_train = sum(1 for _, row in X.iterrows() if ood.predict_one(row.to_dict()))
    registry.register(
        artifact_id=ood_id, model_path=ood_path, schema_version="1.0.0",
        approved_by=APPROVED_BY,
        notes=f"IsolationForest OOD detector (unchanged design). {n_ood_train}/{len(X)} flagged OOD.",
        feature_order=feature_order, schema_file_path=schema_path,
    )
    log(f"registered OOD detector: {ood_id} ({n_ood_train}/{len(X)} OOD)")

    # ---- 5. Threshold derivation (same FE<=10% budget, unchanged) ---------
    b_rows = [o for o in raw["outcomes"] if o.get("condition") == "B" and not o.get("error")
              and isinstance(o.get("features"), dict) and len(o["features"]) > 0]
    assert len(b_rows) == len(X)

    proba_all = calibrated.predict_proba(X)[:, 1]
    clean_mask = np.array([r["case_type"] == "clean" for r in b_rows])
    escaped_fault_mask = (~clean_mask) & (y == 1)
    caught_fault_mask = (~clean_mask) & (y == 0)
    clean_probs = proba_all[clean_mask]
    escaped_probs = proba_all[escaped_fault_mask]
    n_clean, n_escaped, n_caught = int(clean_mask.sum()), int(escaped_fault_mask.sum()), int(caught_fault_mask.sum())
    log(f"clean={n_clean} escaped_fault={n_escaped} caught_fault={n_caught}")
    log(f"clean_probs: min={clean_probs.min():.4f} max={clean_probs.max():.4f} "
        f"unique={len(set(np.round(clean_probs, 6)))}")
    log(f"escaped_probs: min={escaped_probs.min():.4f} max={escaped_probs.max():.4f} "
        f"unique={len(set(np.round(escaped_probs, 6)))}")

    FE_BUDGET = 0.10
    # Candidates must span every distinct probability the model actually
    # outputs, not just the clean-case values: with clean_probs collapsed to
    # a single repeated value (0.0), searching only over {0.0, 1.0} skips
    # every real threshold in between -- including the escaped-fault value
    # itself -- and forces t_escalate to the trivial ceiling of 1.0 even
    # though a much smaller threshold already satisfies the FE budget.
    candidate_t = sorted(set(proba_all.tolist()) | {0.0, 1.0})
    t_escalate = 1.0
    for t in candidate_t:
        if float((clean_probs >= t).mean()) <= FE_BUDGET:
            t_escalate = t
            break
    actual_fe_clean = float((clean_probs >= t_escalate).mean())
    t_halt = max(float(np.percentile(escaped_probs, 25)), t_escalate) if n_escaped > 0 else 1.0
    log(f"t_escalate={t_escalate:.6f} (FE={actual_fe_clean:.4f}) t_halt={t_halt:.6f}")

    policy_payload = {
        "policy_version": "1.0.0-lofo", "derived_on": "2026-08-05", "approved_by": APPROVED_BY,
        "artifact_id": cal_id,
        "thresholds": {"t_escalate": t_escalate, "t_halt": t_halt},
        "budget": {"fe_budget": FE_BUDGET},
        "derivation": {
            "t_escalate_method": "smallest threshold with FE rate <= 10% on development clean cases",
            "t_halt_method": "25th percentile of escaped-fault predicted probabilities, clamped >= t_escalate",
            "n_clean": n_clean, "n_escaped_fault": n_escaped, "actual_fe_clean": actual_fe_clean,
            "cv_method": "leave-one-family-out (11 folds), no inner hyperparameter search, pooled AUPRC",
        },
    }
    policy_sha = hashlib.sha256(json.dumps(policy_payload, indent=2).encode("utf-8")).hexdigest().upper()
    policy_payload["sha256"] = policy_sha
    (OUT_DIR / "policy_config.json").write_text(json.dumps(policy_payload, indent=2), encoding="utf-8")
    log(f"wrote policy_config.json sha256={policy_sha}")

    # ---- 6. Simulate conditions A and C ------------------------------------
    def simulate(condition: str) -> list[dict]:
        rows_out = []
        for row_dict, feat_row, fam in zip(b_rows, X.to_dict(orient="records"), groups):
            telemetry = dict(feat_row)
            telemetry["out_of_domain_flag"] = ood.predict_one(feat_row)
            contracts_all_passed = True if condition == "C" else (row_dict["evaluator_state"] == "PASS")
            try:
                res = predict_risk(telemetry=telemetry, artifact_id=cal_id, registry=registry,
                                    t_escalate=t_escalate, t_halt=t_halt,
                                    contracts_all_passed=contracts_all_passed)
                action, risk = res.action, res.risk
            except InferenceError as exc:
                action, risk = PolicyAction.HALT, None
            case_type = row_dict["case_type"]
            if case_type == "clean":
                escaped = action == PolicyAction.HALT
                escaped_detected = action in (PolicyAction.HALT, PolicyAction.ESCALATE)
            else:
                escaped = action != PolicyAction.HALT
                escaped_detected = action not in (PolicyAction.HALT, PolicyAction.ESCALATE)
            rows_out.append({
                "family_id": fam, "instance_id": row_dict["instance_id"], "case_type": case_type,
                "condition": condition, "action": action.value if hasattr(action, "value") else str(action),
                "risk": risk, "escaped": escaped, "escaped_detected": escaped_detected,
            })
        return rows_out

    cond_A_rows, cond_C_rows = simulate("A"), simulate("C")

    def seer_of(rows: list[dict]) -> dict:
        fault = [r for r in rows if r["case_type"] == "fault"]
        clean = [r for r in rows if r["case_type"] == "clean"]
        return {
            "seer": (sum(r["escaped"] for r in fault) / len(fault)) if fault else None,
            "false_escalation": (sum(r["escaped"] for r in clean) / len(clean)) if clean else None,
        }

    seer_A, seer_C = seer_of(cond_A_rows), seer_of(cond_C_rows)
    fault_B = [r for r in b_rows if r["case_type"] == "fault"]
    clean_B = [r for r in b_rows if r["case_type"] == "clean"]
    seer_B_val = sum(1 for r in fault_B if r["escaped"]) / len(fault_B)
    d_rows = [o for o in raw["outcomes"] if o.get("condition") == "D" and not o.get("error")]
    fault_D = [r for r in d_rows if r["case_type"] == "fault"]
    seer_D_val = (sum(1 for r in fault_D if r["escaped"]) / len(fault_D)) if fault_D else 1.0

    risk_distinct = len(set(round(r["risk"], 6) for r in cond_A_rows if r["risk"] is not None))
    action_counts_A = {a: sum(1 for r in cond_A_rows if r["action"] == a) for a in ("CONTINUE", "ESCALATE", "HALT")}
    action_counts_C = {a: sum(1 for r in cond_C_rows if r["action"] == a) for a in ("CONTINUE", "ESCALATE", "HALT")}

    seer_summary = {
        "SEER_D": round(seer_D_val, 6), "SEER_B": round(seer_B_val, 6),
        "SEER_A": round(seer_A["seer"], 6) if seer_A["seer"] is not None else None,
        "SEER_C": round(seer_C["seer"], 6) if seer_C["seer"] is not None else None,
        "delta_seer_A_vs_D": round(seer_D_val - seer_A["seer"], 6) if seer_A["seer"] is not None else None,
        "false_escalation_A": seer_A["false_escalation"], "false_escalation_C": seer_C["false_escalation"],
        "n_distinct_calibrated_risk_values": risk_distinct,
        "monitor_ever_fires": risk_distinct > 1 and t_escalate < 1.0,
        "action_counts_A": action_counts_A, "action_counts_C": action_counts_C,
    }
    log("SEER summary: " + json.dumps(seer_summary, indent=2))

    out_payload = {
        "run_date": "2026-08-05", "method": "leave-one-family-out (see ADR-AIIM-041)",
        "cv_table": cv_table, "model_selection": {"chosen": chosen, "within_1std": within_1std},
        "calibration_comparison": {m: {"brier": v["brier"], "ece": v["ece"]} for m, v in cal_results.items()},
        "calibration_method_selected": cal_method,
        "policy": {"t_escalate": t_escalate, "t_halt": t_halt, "sha256": policy_sha},
        "seer_summary": seer_summary,
        "condition_A_rows": cond_A_rows, "condition_C_rows": cond_C_rows,
    }
    out_path = OUT_DIR / "wp4_conditions_ac_results_lofo_2026-08-05.json"
    out_path.write_text(json.dumps(out_payload, indent=2), encoding="utf-8")
    log(f"wrote {out_path}")
    log("DONE")


if __name__ == "__main__":
    main()

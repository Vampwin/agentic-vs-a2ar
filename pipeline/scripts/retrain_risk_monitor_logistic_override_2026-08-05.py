"""
retrain_risk_monitor_logistic_override_2026-08-05.py

Third pass on the development risk monitor, addressing the overfitting
caveat raised in ADR-AIIM-041: the LOFO comparison mechanically selected
gradient_boosting by Brier score, but (a) pooled AUPRC actually favoured
logistic_l2 instead, (b) GBM's calibrated output collapsed to two values
(0.0 / 0.1111) across all 188 rows -- a binary flag, not a risk gradient --
and (c) the final fit sees all 9 positive examples with no held-out
protection, which is exactly the regime where a 100-tree depth-3 ensemble
can memorise rather than generalise.

This script does not change the LOFO methodology (still the fix from
ADR-AIIM-041: leave-one-family-out, no inner hyperparameter search, pooled
AUPRC, and the corrected full-landscape threshold search). It changes one
thing: the model selected after comparison. Human decision (2026-08-05):
override the mechanical Brier-only pick and select logistic_l2 instead,
because with 9 positive examples an ensemble's Brier advantage is more
plausibly explained by memorisation than by genuine discrimination, the
prespecified secondary metric (AUPRC) already preferred logistic_l2, and
logistic regression's linear decision boundary and L2 regularisation are
structurally far less able to memorise 9 points than a 100-tree ensemble.

Reuses lofo_compare_models() and lofo_calibrate() from
retrain_risk_monitor_lofo_2026-08-05.py unchanged -- only the selection step
after comparison differs.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path
from kdss_paths import project_root, resolve, under
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import joblib
import numpy as np
import pandas as pd

from fault_injection.pilot_runner import build_training_dataset
from risk_monitor.train import fit_selected_model
from risk_monitor.calibrate import evaluate_calibration, save_calibrated_model
from risk_monitor.registry import ModelRegistry
from risk_monitor.inference import predict_risk, InferenceError
from risk_monitor.policy import PolicyAction
from telemetry.ood_detector import OutOfDomainDetector

_lofo_path = Path(__file__).resolve().parent / "retrain_risk_monitor_lofo_2026-08-05.py"
_spec = importlib.util.spec_from_file_location("_lofo_module", _lofo_path)
_lofo = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_lofo)  # module-level code only defines functions; no main() runs

lofo_compare_models = _lofo.lofo_compare_models
lofo_calibrate = _lofo.lofo_calibrate

BASE = resolve("20_new_runs/WP2.5_WP3A_readiness_2026-07-30")
OUT_DIR = resolve("20_new_runs/WP4_retrain_logistic_override_2026-08-05")
MODELS_DIR = OUT_DIR / "models"
REGISTRY_PATH = BASE / "model_registry.json"
APPROVED_BY = "adisorn.siric@gmail.com"
DATE_TAG = "2026-08-05-logistic-sigmoid-override"
PILOT_JSON = BASE / "pilot_postfix_2026-08-04.json"
FORCED_MODEL = "logistic_l2"


def log(msg: str) -> None:
    print(f"[override] {msg}")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    MODELS_DIR.mkdir(parents=True, exist_ok=True)

    raw = json.loads(PILOT_JSON.read_text(encoding="utf-8"))
    results_ns = SimpleNamespace(outcomes=raw["outcomes"])
    X, y, groups = build_training_dataset(results_ns, reference_condition="B")
    log(f"X={X.shape} y positive_rate={y.mean():.4f} ({int(y.sum())}/{len(y)}) n_groups={len(set(groups))}")

    cv_results = lofo_compare_models(X, y, groups)
    for r in cv_results:
        log(f"LOFO (for the record): {r.model_name:18s} brier={r.mean_brier:.4f}+-{r.std_brier:.4f} "
            f"pooled_auprc={r.mean_auprc:.4f}")
    log(f"OVERRIDE: selecting {FORCED_MODEL} regardless of Brier ranking -- see module docstring "
        "and ADR-AIIM-042 for the human-approved rationale.")
    chosen = FORCED_MODEL

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
        notes=(f"{chosen} -- human-overridden selection (ADR-AIIM-042), not the LOFO Brier "
               f"winner (gradient_boosting). Fit on {len(y)} condition-B rows "
               f"({len(set(groups))} families)."),
        feature_order=feature_order, schema_file_path=schema_path,
    )
    log(f"registered uncalibrated model: {uncal_id}")

    cal_results = {}
    for method in ("isotonic", "sigmoid"):
        cal_pipeline = lofo_calibrate(fitted, X, y, groups, method=method)
        ev = evaluate_calibration(cal_pipeline, X, y)
        ece = float(np.mean(np.abs(np.array(ev["fraction_of_positives"]) - np.array(ev["mean_predicted_value"]))))
        cal_results[method] = {"pipeline": cal_pipeline, "brier": ev["brier_score"], "ece": ece}
        log(f"calibration[{method}]: brier={ev['brier_score']:.4f} ece={ece:.4f}")

    # OVERRIDE: isotonic mechanically wins on Brier/ECE, but ADR-AIIM-029
    # already flagged in the original run that ECE=0.0000 on a small
    # calibration set is a signature of isotonic overfitting the calibration
    # curve rather than genuine perfect calibration -- and the calibration
    # set here (9 positives) is smaller than ADR-029's (113). Isotonic's step
    # function collapses the output to 2 values across all 188 rows; sigmoid
    # (a smooth 2-parameter fit) keeps a continuous ~15-value gradient at an
    # essentially indistinguishable Brier cost. Same reasoning as the model
    # override above: prefer the option less able to overfit 9 points.
    cal_method = "sigmoid"
    log(f"calibration method OVERRIDDEN to {cal_method} (isotonic mechanically wins Brier/ECE "
        f"but collapses to few output values on this small a calibration set -- see ADR-AIIM-042)")

    calibrated = cal_results[cal_method]["pipeline"]
    cal_id = f"{chosen}_cal_{cal_method}_dev_{DATE_TAG}"
    cal_path = MODELS_DIR / f"{cal_id}.joblib"
    save_calibrated_model(calibrated, cal_path)
    registry.register(
        artifact_id=cal_id, model_path=cal_path, schema_version="1.0.0",
        approved_by=APPROVED_BY,
        notes=(f"Calibrated ({cal_method}), human-overridden model selection (ADR-AIIM-042). "
               f"Brier={cal_results[cal_method]['brier']:.4f} ECE={cal_results[cal_method]['ece']:.4f}."),
        feature_order=feature_order, schema_file_path=schema_path,
    )
    log(f"registered calibrated model: {cal_id}")

    ood = OutOfDomainDetector(contamination=0.01, n_estimators=100, random_state=20260729)
    ood.fit(X)
    ood_id = f"ood_isoforest_dev_{DATE_TAG}"
    ood_path = MODELS_DIR / f"{ood_id}.joblib"
    ood.save(ood_path)
    n_ood_train = sum(1 for _, row in X.iterrows() if ood.predict_one(row.to_dict()))
    registry.register(
        artifact_id=ood_id, model_path=ood_path, schema_version="1.0.0",
        approved_by=APPROVED_BY, notes=f"IsolationForest OOD detector. {n_ood_train}/{len(X)} flagged OOD.",
        feature_order=feature_order, schema_file_path=schema_path,
    )
    log(f"registered OOD detector: {ood_id} ({n_ood_train}/{len(X)} OOD)")

    b_rows = [o for o in raw["outcomes"] if o.get("condition") == "B" and not o.get("error")
              and isinstance(o.get("features"), dict) and len(o["features"]) > 0]
    assert len(b_rows) == len(X)

    proba_all = calibrated.predict_proba(X)[:, 1]
    log(f"calibrated risk distribution: min={proba_all.min():.4f} max={proba_all.max():.4f} "
        f"n_distinct={len(set(np.round(proba_all, 6)))}")

    clean_mask = np.array([r["case_type"] == "clean" for r in b_rows])
    escaped_fault_mask = (~clean_mask) & (y == 1)
    clean_probs = proba_all[clean_mask]
    escaped_probs = proba_all[escaped_fault_mask]
    n_clean, n_escaped = int(clean_mask.sum()), int(escaped_fault_mask.sum())

    FE_BUDGET = 0.10
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
        "policy_version": "1.0.0-logistic-override", "derived_on": "2026-08-05",
        "approved_by": APPROVED_BY, "artifact_id": cal_id,
        "thresholds": {"t_escalate": t_escalate, "t_halt": t_halt},
        "budget": {"fe_budget": FE_BUDGET},
        "derivation": {
            "t_escalate_method": "smallest threshold with FE rate <= 10% on development clean cases "
                                  "(candidate set = full observed probability landscape, ADR-AIIM-041 fix)",
            "t_halt_method": "25th percentile of escaped-fault predicted probabilities, clamped >= t_escalate",
            "n_clean": n_clean, "n_escaped_fault": n_escaped,
            "cv_method": "leave-one-family-out (11 folds), model forced to logistic_l2 (ADR-AIIM-042)",
        },
    }
    policy_sha = hashlib.sha256(json.dumps(policy_payload, indent=2).encode("utf-8")).hexdigest().upper()
    policy_payload["sha256"] = policy_sha
    (OUT_DIR / "policy_config.json").write_text(json.dumps(policy_payload, indent=2), encoding="utf-8")
    log(f"wrote policy_config.json sha256={policy_sha}")

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
            except InferenceError:
                action, risk = PolicyAction.HALT, None
            case_type = row_dict["case_type"]
            if case_type == "clean":
                escaped = action == PolicyAction.HALT
            else:
                escaped = action != PolicyAction.HALT
            rows_out.append({
                "family_id": fam, "instance_id": row_dict["instance_id"], "case_type": case_type,
                "condition": condition, "action": action.value if hasattr(action, "value") else str(action),
                "risk": risk, "escaped": escaped,
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
    seer_B_val = sum(1 for r in fault_B if r["escaped"]) / len(fault_B)
    d_rows = [o for o in raw["outcomes"] if o.get("condition") == "D" and not o.get("error")]
    fault_D = [r for r in d_rows if r["case_type"] == "fault"]
    seer_D_val = (sum(1 for r in fault_D if r["escaped"]) / len(fault_D)) if fault_D else 1.0

    seer_summary = {
        "SEER_D": round(seer_D_val, 6), "SEER_B": round(seer_B_val, 6),
        "SEER_A": round(seer_A["seer"], 6) if seer_A["seer"] is not None else None,
        "SEER_C": round(seer_C["seer"], 6) if seer_C["seer"] is not None else None,
        "delta_seer_A_vs_D": round(seer_D_val - seer_A["seer"], 6) if seer_A["seer"] is not None else None,
        "false_escalation_A": seer_A["false_escalation"], "false_escalation_C": seer_C["false_escalation"],
        "n_distinct_calibrated_risk_values": len(set(np.round(proba_all, 6))),
    }
    log("SEER summary: " + json.dumps(seer_summary, indent=2))

    out_payload = {
        "run_date": "2026-08-05", "method": "leave-one-family-out CV, model forced to logistic_l2 (ADR-AIIM-042)",
        "model_selection": {"chosen": chosen, "override": True, "lofo_brier_winner": "gradient_boosting"},
        "calibration_comparison": {m: {"brier": v["brier"], "ece": v["ece"]} for m, v in cal_results.items()},
        "calibration_method_selected": cal_method,
        "policy": {"t_escalate": t_escalate, "t_halt": t_halt, "sha256": policy_sha},
        "seer_summary": seer_summary,
        "condition_A_rows": cond_A_rows, "condition_C_rows": cond_C_rows,
    }
    out_path = OUT_DIR / "wp4_conditions_ac_results_logistic_override_2026-08-05.json"
    out_path.write_text(json.dumps(out_payload, indent=2), encoding="utf-8")
    log(f"wrote {out_path}")
    log("DONE")


if __name__ == "__main__":
    main()

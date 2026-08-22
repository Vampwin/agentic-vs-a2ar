"""
inference.py — Fail-closed risk-monitor inference (ADR-AIIM-018, ADR-AIIM-020).

Loads a trusted registered scikit-learn pipeline and returns a calibrated
unsafe-release probability for the current pre-release telemetry record.

Fail-closed conditions — any of these raises InferenceError:
  - model artifact not in registry
  - hash mismatch on the registered artifact
  - telemetry contains a forbidden or unknown field (LeakageError)
  - required non-nullable feature is missing (None when not nullable)
  - out_of_domain_flag is None (unresolved — domain check not run)
  - out_of_domain_flag is True (record is out-of-domain; prediction not trusted)
  - telemetry schema version does not match the model's training schema
  - feature order in registry does not match pipeline.feature_names_in_
  - prediction error or NaN output from the pipeline

The inference engine cannot:
  - alter contract outcomes
  - write to the artifact registry
  - read held-out labels, fault identity, or experimental condition
  - override a report_lock block
"""

from __future__ import annotations

import importlib
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .policy import PolicyAction, PolicyResult, evaluate_policy
from .registry import ModelRegistry, RegistryError
from telemetry.leakage_guard import (
    LeakageError,
    validate_telemetry,
    get_allowed_names,
    get_schema_version,
    get_canonical_feature_order,
)

logger = logging.getLogger(__name__)

# Features that must be non-null (None not allowed) for any inference to proceed.
# Stage-dependent nullable features (scientific_summaries, shadow_rmsd*, shadow_score_rank_rho)
# and out_of_domain_flag (handled separately) are not in this set.
_REQUIRED_NON_NULL = frozenset({
    "n_missing_input_hashes",
    "n_unregistered_artifacts",
    "has_unexpected_lineage",
    "n_tool_warnings",
    "n_retries",
    "elapsed_seconds",
    "n_nonzero_exit_codes",
    "n_missing_outputs",
    "n_schema_drift_fields",
    "n_count_discrepancies",
    "n_duplicate_ids",
    "n_lost_records",
    "n_unexpected_records",
    "timestamp_staleness_seconds",
    "shadow_id_match_fraction",
    "shadow_hash_match_all",
    "n_null_features",
    # out_of_domain_flag is NOT here — it is nullable (None = unresolved)
    # and handled with its own fail-closed logic below.
})


class InferenceError(RuntimeError):
    """Raised when inference cannot proceed safely; the caller must ESCALATE or HALT."""


@dataclass(frozen=True)
class InferenceResult:
    risk: float
    action: PolicyAction
    reason: str
    artifact_id: str
    schema_version_match: bool


def _import_joblib():
    try:
        return importlib.import_module("joblib")
    except ImportError as exc:
        raise InferenceError(
            "joblib is required for model loading but is not installed"
        ) from exc


def _import_numpy():
    try:
        return importlib.import_module("numpy")
    except ImportError as exc:
        raise InferenceError("numpy is required for inference but is not installed") from exc


def predict_risk(
    telemetry: dict[str, Any],
    artifact_id: str,
    registry: ModelRegistry,
    t_escalate: float,
    t_halt: float,
    contracts_all_passed: bool,
) -> InferenceResult:
    """
    Return a calibrated risk estimate and the resulting policy action.

    Parameters
    ----------
    telemetry:
        Output of telemetry.feature_extractor.extract(), with out_of_domain_flag
        resolved (False = in-domain, True = out-of-domain, None = unresolved).
        Inference fails closed when out_of_domain_flag is None or True.
    artifact_id:
        Registry key for the frozen trusted model artifact.
    registry:
        ModelRegistry instance; used to verify and locate the model file.
    t_escalate, t_halt:
        Frozen policy thresholds (must come from the human-approved policy config).
    contracts_all_passed:
        True only if every required hard contract for this stage passed.
    """
    # 1. Leakage / allowlist check — fail closed
    try:
        validate_telemetry(telemetry, context="inference.predict_risk")
    except LeakageError as exc:
        raise InferenceError(f"Leakage guard rejected telemetry: {exc}") from exc

    # 2. Required non-null features
    missing = [f for f in _REQUIRED_NON_NULL if telemetry.get(f) is None]
    if missing:
        raise InferenceError(
            f"Required non-nullable features missing: {missing}.  "
            "Cannot proceed with inference — failing closed."
        )

    # 3. out_of_domain_flag: fail closed on unresolved (None) or out-of-domain (True)
    ood = telemetry.get("out_of_domain_flag")
    if ood is None:
        raise InferenceError(
            "out_of_domain_flag is None (unresolved): the applicability-domain check "
            "has not been run. Inference fails closed until domain status is resolved."
        )
    if ood is True:
        raise InferenceError(
            "out_of_domain_flag is True: this telemetry record is outside the "
            "monitor's training distribution. Inference fails closed for out-of-domain input."
        )

    # 4. Schema version match
    current_schema_version = get_schema_version()
    schema_version_match = True
    try:
        entry = registry.get_entry(artifact_id)
        if entry.schema_version != current_schema_version:
            logger.warning(
                "Schema version mismatch: model trained on %s, current schema is %s. "
                "Failing closed.",
                entry.schema_version,
                current_schema_version,
            )
            schema_version_match = False
            raise InferenceError(
                f"Telemetry schema version mismatch: model={entry.schema_version!r}, "
                f"current={current_schema_version!r}"
            )
    except RegistryError as exc:
        raise InferenceError(f"Registry lookup failed: {exc}") from exc

    # 5. Load and verify the model artifact (hash check)
    joblib = _import_joblib()
    np = _import_numpy()
    try:
        model_path = registry.verify_and_get_path(artifact_id)
    except RegistryError as exc:
        raise InferenceError(f"Model artifact verification failed: {exc}") from exc

    try:
        pipeline = joblib.load(model_path)
    except Exception as exc:
        raise InferenceError(f"Failed to load model from {model_path}: {exc}") from exc

    # 6. Determine canonical feature order
    # Prefer the order stored in the registry (from training), then canonical schema order.
    feature_order = list(entry.feature_order) if entry.feature_order else get_canonical_feature_order()

    # 7. Verify pipeline feature names against stored order (if available)
    if hasattr(pipeline, "feature_names_in_"):
        pipeline_order = list(pipeline.feature_names_in_)
        if pipeline_order != feature_order:
            raise InferenceError(
                f"Feature order mismatch: pipeline.feature_names_in_ has "
                f"{pipeline_order[:5]}… but registry feature_order has {feature_order[:5]}…. "
                "Model was likely trained with a different feature order."
            )

    # 8. Build feature vector in the training order (as DataFrame to match training names)
    try:
        import pandas as pd
        X = pd.DataFrame(
            [[_safe_float(telemetry.get(f)) for f in feature_order]],
            columns=feature_order,
        )
    except Exception as exc:
        raise InferenceError(f"Feature vector construction failed: {exc}") from exc

    # 9. Predict calibrated probability
    try:
        proba = pipeline.predict_proba(X)
        if proba.shape != (1, 2):
            raise InferenceError(
                f"Unexpected predict_proba shape {proba.shape}; expected (1, 2)"
            )
        risk = float(proba[0, 1])
    except InferenceError:
        raise
    except Exception as exc:
        raise InferenceError(f"Pipeline prediction failed: {exc}") from exc

    if not (0.0 <= risk <= 1.0) or np.isnan(float(risk)):
        raise InferenceError(
            f"Pipeline returned invalid risk probability: {risk!r}"
        )

    # 10. Apply frozen policy (hard-contract precedence enforced inside evaluate_policy)
    policy_result: PolicyResult = evaluate_policy(
        risk=risk,
        t_escalate=t_escalate,
        t_halt=t_halt,
        contracts_all_passed=contracts_all_passed,
    )

    return InferenceResult(
        risk=risk,
        action=policy_result.action,
        reason=policy_result.reason,
        artifact_id=artifact_id,
        schema_version_match=schema_version_match,
    )


def _safe_float(value: Any) -> float:
    """Convert a telemetry value to float, using NaN for None/missing."""
    import math
    if value is None:
        return float("nan")
    if isinstance(value, bool):
        return 1.0 if value else 0.0
    return float(value)

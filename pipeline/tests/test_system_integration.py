"""
test_system_integration.py — System integration tests (ADR-AIIM-018, ADR-AIIM-021).

Verifies end-to-end behavior using real ArtifactRegistry and ModelRegistry objects:

  1. ML artifact records are written and recoverable from ArtifactRegistry.
  2. Non-empty approved_by is required for ML model registration.
  3. Failed contract overrides low monitor risk → HALT action.
  4. OOD and missing artifact conditions fail closed through the full pipeline.
  5. ReportLock cannot be released based on low monitor risk alone.
"""

import numpy as np
import pandas as pd
import pytest
import joblib
from pathlib import Path

from telemetry.leakage_guard import get_canonical_feature_order, get_schema_version
from risk_monitor.inference import InferenceError, predict_risk
from risk_monitor.registry import ModelRegistry, RegistryError
from risk_monitor.policy import PolicyAction


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_and_register_pipeline(tmp_path: Path, registry: ModelRegistry) -> object:
    """Train a dummy pipeline and register it in the ModelRegistry."""
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    feature_order = get_canonical_feature_order()
    rng = np.random.default_rng(42)
    X = pd.DataFrame(rng.random((20, len(feature_order))), columns=feature_order)
    y = rng.integers(0, 2, size=20)

    pipeline = Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(max_iter=100, random_state=0)),
    ])
    pipeline.fit(X, y)

    model_file = tmp_path / "model.joblib"
    joblib.dump(pipeline, model_file)

    entry = registry.register(
        artifact_id="integration_model_v1",
        model_path=model_file,
        schema_version=get_schema_version(),
        approved_by="integration_test_fixture",
        feature_order=feature_order,
    )
    return entry, registry


def _clean_telemetry(out_of_domain: bool = False) -> dict:
    return {
        "n_missing_input_hashes": 0,
        "n_unregistered_artifacts": 0,
        "has_unexpected_lineage": False,
        "n_tool_warnings": 0,
        "n_retries": 0,
        "elapsed_seconds": 5.0,
        "n_nonzero_exit_codes": 0,
        "n_missing_outputs": 0,
        "n_schema_drift_fields": 0,
        "n_count_discrepancies": 0,
        "n_duplicate_ids": 0,
        "n_lost_records": 0,
        "n_unexpected_records": 0,
        "timestamp_staleness_seconds": 0.0,
        "pose_failure_fraction": None,
        "score_distribution_shift": None,
        "n_steric_clashes": None,
        "convergence_flag": None,
        "shadow_rmsd_best_pose": None,
        "shadow_score_rank_rho": None,
        "shadow_id_match_fraction": 1.0,
        "shadow_hash_match_all": True,
        "n_null_features": 6,
        "out_of_domain_flag": out_of_domain,
    }


# ---------------------------------------------------------------------------
# 1. ArtifactRegistry ML records
# ---------------------------------------------------------------------------

def test_artifact_registry_records_ml_schema(tmp_path):
    """register_ml_schema_artifact() writes a recoverable record."""
    from kdss_audit.artifact_registry import ArtifactRegistry
    reg = ArtifactRegistry(tmp_path / "registry.jsonl")

    reg.register_ml_schema_artifact(
        schema_path="telemetry/feature_schema.yaml",
        sha256="AABBCC",
        schema_version="PENDING_HUMAN_APPROVAL",
        approved_by="test_user",
        feature_count=24,
        notes="integration test",
    )

    ml_records = reg.get_ml_records(record_type="ml_schema_artifact")
    assert len(ml_records) == 1
    assert ml_records[0]["feature_count"] == 24
    assert ml_records[0]["approved_by"] == "test_user"


def test_artifact_registry_records_ml_model(tmp_path):
    """register_ml_model_artifact() writes a recoverable record."""
    from kdss_audit.artifact_registry import ArtifactRegistry
    reg = ArtifactRegistry(tmp_path / "registry.jsonl")

    reg.register_ml_model_artifact(
        artifact_id="model_v1",
        model_path="/path/to/model.joblib",
        sha256="DEADBEEF",
        schema_version="PENDING_HUMAN_APPROVAL",
        approved_by="test_reviewer",
        model_name="logistic_l2",
        mean_brier=0.25,
    )

    ml_records = reg.get_ml_records(record_type="ml_model_artifact")
    assert len(ml_records) == 1
    assert ml_records[0]["model_name"] == "logistic_l2"
    assert ml_records[0]["approved_by"] == "test_reviewer"


def test_artifact_registry_rejects_empty_approved_by_for_model(tmp_path):
    """register_ml_model_artifact() raises ValueError when approved_by is empty."""
    from kdss_audit.artifact_registry import ArtifactRegistry
    reg = ArtifactRegistry(tmp_path / "registry.jsonl")

    with pytest.raises(ValueError, match="approved_by"):
        reg.register_ml_model_artifact(
            artifact_id="model_v1",
            model_path="/path/to/model.joblib",
            sha256="DEADBEEF",
            schema_version="PENDING_HUMAN_APPROVAL",
            approved_by="",   # empty — must be rejected
            model_name="logistic_l2",
            mean_brier=0.25,
        )


def test_artifact_registry_records_prediction_event(tmp_path):
    """record_ml_prediction_event() writes a recoverable audit record."""
    from kdss_audit.artifact_registry import ArtifactRegistry
    reg = ArtifactRegistry(tmp_path / "registry.jsonl")

    reg.record_ml_prediction_event(
        artifact_id="model_v1",
        risk_score=0.45,
        action="ESCALATE",
        reason="risk above t_escalate",
        schema_version_match=True,
        out_of_domain_flag=False,
        workflow_stage="G3_pose_generation",
    )

    pred_records = reg.get_ml_records(record_type="ml_prediction_event")
    assert len(pred_records) == 1
    assert pred_records[0]["risk_score"] == pytest.approx(0.45)
    assert pred_records[0]["action"] == "ESCALATE"


# ---------------------------------------------------------------------------
# 2. Contract precedence overrides monitor
# ---------------------------------------------------------------------------

def test_failed_contract_halts_regardless_of_low_risk(tmp_path):
    """
    When contracts_all_passed=False, predict_risk returns HALT even when
    the monitor's calibrated risk is below t_escalate.
    This enforces ADR-AIIM-018: hard contracts take precedence over the learned monitor.
    """
    registry = ModelRegistry(tmp_path / "registry.json")
    entry, registry = _make_and_register_pipeline(tmp_path, registry)

    result = predict_risk(
        telemetry=_clean_telemetry(out_of_domain=False),
        artifact_id=entry.artifact_id,
        registry=registry,
        t_escalate=0.5,
        t_halt=0.8,
        contracts_all_passed=False,   # contract failed
    )
    assert result.action == PolicyAction.HALT, (
        "A failed contract must produce HALT regardless of monitor risk score. "
        f"Got action={result.action} with risk={result.risk:.4f}."
    )


# ---------------------------------------------------------------------------
# 3. OOD and missing artifact fail closed through full pipeline
# ---------------------------------------------------------------------------

def test_ood_true_fails_closed_integration(tmp_path):
    """Full pipeline: OOD=True raises InferenceError (fail closed)."""
    registry = ModelRegistry(tmp_path / "registry.json")
    entry, registry = _make_and_register_pipeline(tmp_path, registry)

    telemetry = _clean_telemetry(out_of_domain=False)
    telemetry["out_of_domain_flag"] = True

    with pytest.raises(InferenceError):
        predict_risk(
            telemetry=telemetry,
            artifact_id=entry.artifact_id,
            registry=registry,
            t_escalate=0.3,
            t_halt=0.7,
            contracts_all_passed=True,
        )


def test_missing_artifact_fails_closed_integration(tmp_path):
    """Full pipeline: deleted model file raises InferenceError (fail closed)."""
    registry = ModelRegistry(tmp_path / "registry.json")
    entry, registry = _make_and_register_pipeline(tmp_path, registry)

    # Delete the model file after registration
    model_path = Path(entry.model_path)
    model_path.unlink()

    with pytest.raises(InferenceError):
        predict_risk(
            telemetry=_clean_telemetry(out_of_domain=False),
            artifact_id=entry.artifact_id,
            registry=registry,
            t_escalate=0.3,
            t_halt=0.7,
            contracts_all_passed=True,
        )


# ---------------------------------------------------------------------------
# 4. ReportLock cannot be released by monitor risk alone
# ---------------------------------------------------------------------------

def test_report_lock_not_releasable_by_monitor(tmp_path):
    """
    ReportLock.can_report() returns False when the artifact lacks a contract PASS.
    A monitor risk score has no mechanism to release the lock (no accept_risk method).
    This enforces ADR-AIIM-018: the learned monitor cannot override contract outcomes.
    """
    from kdss_audit.artifact_registry import ArtifactRegistry
    from kdss_audit.report_lock import ReportLock

    # Create registry with NO contract pass records
    registry = ArtifactRegistry(tmp_path / "audit_registry.jsonl")
    lock_file = tmp_path / "report_lock.json"
    lock = ReportLock(registry=registry, lock_file=lock_file)

    # With no contract PASS, can_report must return False
    assert not lock.can_report("some_artifact"), (
        "can_report() must return False when the artifact has no passed contract record."
    )

    # assert_can_report must raise PermissionError (not silently pass)
    with pytest.raises(PermissionError):
        lock.assert_can_report(["some_artifact"])

    # ReportLock must have no method that accepts a risk score as release criterion
    assert not hasattr(lock, "accept_risk"), (
        "ReportLock must not have an accept_risk() method; "
        "monitor risk cannot release the report lock (ADR-AIIM-018)."
    )
    assert not hasattr(lock, "release_if_low_risk"), (
        "ReportLock must not auto-release based on risk score."
    )

    # Confirm the lock's condition is still not met after simulating monitor output
    simulated_low_risk = 0.05  # noqa: F841 — never passed to ReportLock
    assert not lock.can_report("some_artifact"), (
        "can_report() must remain False regardless of monitor risk score; "
        "only a contract PASS or human adjudication record can release it."
    )

"""
test_monitor_fail_closed.py — Fail-closed behavior tests (ADR-AIIM-018).

Verifies that inference.predict_risk raises InferenceError (never returns a
result that could be mistaken for a valid low-risk assessment) for every
defined failure mode:

  1. Model artifact not in registry.
  2. Hash mismatch between registered and on-disk artifact.
  3. Required non-nullable feature missing from telemetry.
  4. Telemetry contains a forbidden field (LeakageError → InferenceError).
  5. Telemetry contains an unknown non-forbidden field (LeakageError → InferenceError).
  6. Telemetry schema version does not match the model's training version.
  7. out_of_domain_flag is True (OOD — prediction not trusted).
  8. out_of_domain_flag is None (unresolved — domain check not run).
  9. Feature order in pipeline does not match registry's stored feature_order.
"""

import json
import tempfile
from pathlib import Path

import joblib
import numpy as np
import pytest
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer

from risk_monitor.inference import InferenceError, InferenceResult, predict_risk
from risk_monitor.registry import ModelRegistry, RegistryError
from risk_monitor.policy import PolicyAction

from telemetry.leakage_guard import (
    get_allowed_names,
    get_canonical_feature_order,
    get_schema_version,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _minimal_clean_telemetry(out_of_domain: bool = False) -> dict:
    """
    A fully valid telemetry record that passes all guards.
    out_of_domain_flag defaults to False (in-domain, domain check resolved).
    """
    return {
        "n_missing_input_hashes": 0,
        "n_unregistered_artifacts": 0,
        "has_unexpected_lineage": False,
        "n_tool_warnings": 0,
        "n_retries": 0,
        "elapsed_seconds": 10.0,
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


def _make_dummy_pipeline(feature_order: list[str] | None = None) -> Pipeline:
    """
    Create a minimal sklearn Pipeline for testing.
    Trained on canonical feature order so pipeline.feature_names_in_ matches
    the order that inference.py derives from the schema.
    """
    import pandas as pd
    if feature_order is None:
        feature_order = get_canonical_feature_order()
    n_features = len(feature_order)

    pipeline = Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(max_iter=100, random_state=0)),
    ])
    rng = np.random.default_rng(0)
    X_dummy = pd.DataFrame(
        rng.random((20, n_features)), columns=feature_order
    )
    y_dummy = rng.integers(0, 2, size=20)
    pipeline.fit(X_dummy, y_dummy)
    return pipeline


def _current_schema_version() -> str:
    return get_schema_version()


@pytest.fixture
def tmp_registry(tmp_path):
    """A ModelRegistry backed by a temporary directory."""
    return ModelRegistry(tmp_path / "test_registry.json")


@pytest.fixture
def registered_model(tmp_registry, tmp_path):
    """A minimal trained Pipeline registered with the current schema version."""
    model_file = tmp_path / "model.joblib"
    pipeline = _make_dummy_pipeline()
    joblib.dump(pipeline, model_file)
    schema_version = _current_schema_version()
    entry = tmp_registry.register(
        artifact_id="test_model_v1",
        model_path=model_file,
        schema_version=schema_version,
        approved_by="test_fixture",
        notes="pytest fixture",
        feature_order=get_canonical_feature_order(),
    )
    return entry, tmp_registry


# ---------------------------------------------------------------------------
# Failure mode 1: artifact not in registry
# ---------------------------------------------------------------------------

def test_fail_closed_unknown_artifact(tmp_registry):
    """predict_risk raises InferenceError when artifact_id is not in the registry."""
    with pytest.raises(InferenceError):
        predict_risk(
            telemetry=_minimal_clean_telemetry(),
            artifact_id="does_not_exist",
            registry=tmp_registry,
            t_escalate=0.3,
            t_halt=0.7,
            contracts_all_passed=True,
        )


# ---------------------------------------------------------------------------
# Failure mode 2: hash mismatch
# ---------------------------------------------------------------------------

def test_fail_closed_hash_mismatch(tmp_path):
    """predict_risk raises InferenceError when the on-disk model hash has changed."""
    registry = ModelRegistry(tmp_path / "registry.json")
    model_file = tmp_path / "model.joblib"
    pipeline = _make_dummy_pipeline()
    joblib.dump(pipeline, model_file)
    registry.register(
        artifact_id="v1",
        model_path=model_file,
        schema_version=_current_schema_version(),
        approved_by="test",
        feature_order=get_canonical_feature_order(),
    )

    # Tamper with the file after registration
    model_file.write_bytes(model_file.read_bytes() + b"\x00TAMPERED")

    with pytest.raises(InferenceError):
        predict_risk(
            telemetry=_minimal_clean_telemetry(),
            artifact_id="v1",
            registry=registry,
            t_escalate=0.3,
            t_halt=0.7,
            contracts_all_passed=True,
        )


# ---------------------------------------------------------------------------
# Failure mode 3: required non-nullable feature missing
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("missing_key", [
    "n_missing_input_hashes",
    "elapsed_seconds",
    "n_tool_warnings",
    "shadow_hash_match_all",
    "n_null_features",
])
def test_fail_closed_missing_required_feature(missing_key, registered_model):
    """predict_risk raises InferenceError when a required non-nullable feature is None."""
    entry, registry = registered_model
    telemetry = _minimal_clean_telemetry()
    telemetry[missing_key] = None

    with pytest.raises(InferenceError):
        predict_risk(
            telemetry=telemetry,
            artifact_id=entry.artifact_id,
            registry=registry,
            t_escalate=0.3,
            t_halt=0.7,
            contracts_all_passed=True,
        )


# ---------------------------------------------------------------------------
# Failure mode 4: telemetry contains a forbidden field
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("bad_key", [
    "camp_value",
    "fault_id",
    "condition_label",
    "held_out_label",
    "contract_decision",
    "assay_result",
])
def test_fail_closed_forbidden_field(bad_key, registered_model):
    """predict_risk raises InferenceError when telemetry contains a forbidden field."""
    entry, registry = registered_model
    telemetry = _minimal_clean_telemetry()
    telemetry[bad_key] = 0

    with pytest.raises(InferenceError):
        predict_risk(
            telemetry=telemetry,
            artifact_id=entry.artifact_id,
            registry=registry,
            t_escalate=0.3,
            t_halt=0.7,
            contracts_all_passed=True,
        )


# ---------------------------------------------------------------------------
# Failure mode 5: unknown non-forbidden field
# ---------------------------------------------------------------------------

def test_fail_closed_unknown_field(registered_model):
    """predict_risk raises InferenceError when telemetry contains an unknown field."""
    entry, registry = registered_model
    telemetry = _minimal_clean_telemetry()
    telemetry["my_undeclared_feature"] = 99

    with pytest.raises(InferenceError):
        predict_risk(
            telemetry=telemetry,
            artifact_id=entry.artifact_id,
            registry=registry,
            t_escalate=0.3,
            t_halt=0.7,
            contracts_all_passed=True,
        )


# ---------------------------------------------------------------------------
# Failure mode 6: schema version mismatch
# ---------------------------------------------------------------------------

def test_fail_closed_schema_version_mismatch(tmp_path):
    """predict_risk raises InferenceError when the model was trained on a different schema."""
    registry = ModelRegistry(tmp_path / "registry.json")
    model_file = tmp_path / "model.joblib"
    pipeline = _make_dummy_pipeline()
    joblib.dump(pipeline, model_file)
    registry.register(
        artifact_id="old_schema_model",
        model_path=model_file,
        schema_version="0.9.0",   # intentionally wrong
        approved_by="test",
        feature_order=get_canonical_feature_order(),
    )

    with pytest.raises(InferenceError):
        predict_risk(
            telemetry=_minimal_clean_telemetry(),
            artifact_id="old_schema_model",
            registry=registry,
            t_escalate=0.3,
            t_halt=0.7,
            contracts_all_passed=True,
        )


# ---------------------------------------------------------------------------
# Failure mode 7: out_of_domain_flag is True (OOD)
# ---------------------------------------------------------------------------

def test_fail_closed_ood_flag_true(registered_model):
    """
    predict_risk raises InferenceError when out_of_domain_flag is True.
    An out-of-domain record must never produce a trusted prediction.
    """
    entry, registry = registered_model
    telemetry = _minimal_clean_telemetry()
    telemetry["out_of_domain_flag"] = True

    with pytest.raises(InferenceError, match="out.of.domain"):
        predict_risk(
            telemetry=telemetry,
            artifact_id=entry.artifact_id,
            registry=registry,
            t_escalate=0.3,
            t_halt=0.7,
            contracts_all_passed=True,
        )


# ---------------------------------------------------------------------------
# Failure mode 8: out_of_domain_flag is None (unresolved)
# ---------------------------------------------------------------------------

def test_fail_closed_ood_flag_none(registered_model):
    """
    predict_risk raises InferenceError when out_of_domain_flag is None.
    The domain check must be resolved before inference can proceed.
    """
    entry, registry = registered_model
    telemetry = _minimal_clean_telemetry()
    telemetry["out_of_domain_flag"] = None

    with pytest.raises(InferenceError, match="unresolved|out.of.domain"):
        predict_risk(
            telemetry=telemetry,
            artifact_id=entry.artifact_id,
            registry=registry,
            t_escalate=0.3,
            t_halt=0.7,
            contracts_all_passed=True,
        )


# ---------------------------------------------------------------------------
# Failure mode 9: feature order mismatch
# ---------------------------------------------------------------------------

def test_fail_closed_feature_order_mismatch(tmp_path):
    """
    predict_risk raises InferenceError when pipeline.feature_names_in_ does not
    match the canonical feature order stored in the registry.
    A permuted training order would cause silently wrong predictions.
    """
    registry = ModelRegistry(tmp_path / "registry.json")
    model_file = tmp_path / "model_wrong_order.joblib"

    # Train on REVERSED canonical order
    canonical = get_canonical_feature_order()
    reversed_order = list(reversed(canonical))
    pipeline = _make_dummy_pipeline(feature_order=reversed_order)
    joblib.dump(pipeline, model_file)

    # Register with the CORRECT canonical order (mismatch with training)
    registry.register(
        artifact_id="mismatched_order",
        model_path=model_file,
        schema_version=_current_schema_version(),
        approved_by="test",
        feature_order=canonical,   # stored order differs from pipeline's fitted order
    )

    with pytest.raises(InferenceError, match="[Ff]eature.order"):
        predict_risk(
            telemetry=_minimal_clean_telemetry(),
            artifact_id="mismatched_order",
            registry=registry,
            t_escalate=0.3,
            t_halt=0.7,
            contracts_all_passed=True,
        )


# ---------------------------------------------------------------------------
# Positive control: valid inference succeeds
# ---------------------------------------------------------------------------

def test_valid_inference_returns_result(registered_model):
    """
    predict_risk returns an InferenceResult (not InferenceError) when all
    conditions are met: correct schema version, OOD=False, canonical feature order,
    matching hash, and all required fields present.
    """
    entry, registry = registered_model
    telemetry = _minimal_clean_telemetry(out_of_domain=False)

    result = predict_risk(
        telemetry=telemetry,
        artifact_id=entry.artifact_id,
        registry=registry,
        t_escalate=0.3,
        t_halt=0.7,
        contracts_all_passed=True,
    )
    assert isinstance(result, InferenceResult)
    assert 0.0 <= result.risk <= 1.0
    assert result.artifact_id == entry.artifact_id
    assert result.schema_version_match is True


# ---------------------------------------------------------------------------
# Registry-level fail-closed tests
# ---------------------------------------------------------------------------

def test_registry_fails_closed_on_missing_file(tmp_path):
    """verify_and_get_path raises RegistryError if the registered file was deleted."""
    registry = ModelRegistry(tmp_path / "registry.json")
    model_file = tmp_path / "model.joblib"
    pipeline = _make_dummy_pipeline()
    joblib.dump(pipeline, model_file)
    registry.register(
        artifact_id="v1",
        model_path=model_file,
        schema_version=_current_schema_version(),
        approved_by="test",
    )
    model_file.unlink()

    with pytest.raises(RegistryError):
        registry.verify_and_get_path("v1")


def test_registry_fails_closed_on_unregistered_id(tmp_registry):
    """verify_and_get_path raises RegistryError for an unknown artifact_id."""
    with pytest.raises(RegistryError):
        tmp_registry.verify_and_get_path("nonexistent")


def test_registry_append_only_rejects_duplicate(tmp_path):
    """Registering the same artifact_id twice raises RegistryError (append-only)."""
    registry = ModelRegistry(tmp_path / "registry.json")
    model_file = tmp_path / "model.joblib"
    joblib.dump(_make_dummy_pipeline(), model_file)
    registry.register(
        artifact_id="v1",
        model_path=model_file,
        schema_version=_current_schema_version(),
        approved_by="test",
    )
    with pytest.raises(RegistryError, match="already registered"):
        registry.register(
            artifact_id="v1",
            model_path=model_file,
            schema_version=_current_schema_version(),
            approved_by="test",
        )


def test_registry_rejects_empty_approved_by(tmp_path):
    """Registering with empty approved_by raises RegistryError."""
    registry = ModelRegistry(tmp_path / "registry.json")
    model_file = tmp_path / "model.joblib"
    joblib.dump(_make_dummy_pipeline(), model_file)
    with pytest.raises(RegistryError, match="approved_by"):
        registry.register(
            artifact_id="v1",
            model_path=model_file,
            schema_version=_current_schema_version(),
            approved_by="",
        )

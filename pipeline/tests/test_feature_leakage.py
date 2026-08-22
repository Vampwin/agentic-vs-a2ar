"""
test_feature_leakage.py — Forbidden-field and allowlist enforcement tests.

Verifies that leakage_guard raises LeakageError for all forbidden-feature categories
defined in ADR-AIIM-020 and feature_schema.yaml, and that known permitted features
pass without error.  No sklearn or shap required.
"""

import pytest
from telemetry.leakage_guard import (
    LeakageError,
    check_field_names,
    check_allowlist,
    validate_telemetry,
    validate_record,
    get_canonical_feature_order,
    get_allowed_names,
)


# ---------------------------------------------------------------------------
# Forbidden pattern tests
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("bad_field", [
    # cAMP / assay
    "camp_value",
    "cAMP_label",
    "assay_result",
    "assay_activity",
    # benchmark labels
    "benchmark_label",
    "known_active_flag",
    "known_inactive_flag",
    # fault identity
    "fault_id",
    "fault_family",
    "fault_severity",
    "expected_gate",
    # experimental condition
    "condition_id",
    "condition_label",
    "experimental_condition",
    # held-out / post-release
    "held_out_label",
    "post_release_outcome",
    "evaluator_outcome",
    # categorical contract decision
    "contract_decision",
    "contract_outcome",
    "contract_pass",
    "contract_fail",
    "gate_pass",
    "gate_fail",
    # assay-informed candidate selection
    "old_candidate",
    "assay_informed_threshold",
    "similarity_to_active",
])
def test_check_field_names_rejects_forbidden(bad_field):
    """check_field_names raises LeakageError for every forbidden pattern."""
    with pytest.raises(LeakageError):
        check_field_names([bad_field], context="test")


def test_check_field_names_case_insensitive():
    """Pattern matching is case-insensitive."""
    for variant in ["CAMP_VALUE", "Camp_Value", "cAmP_VaLuE"]:
        with pytest.raises(LeakageError):
            check_field_names([variant], context="test")


def test_check_field_names_accepts_known_permitted():
    """Known permitted field names pass without error."""
    permitted = [
        "n_missing_input_hashes",
        "n_tool_warnings",
        "elapsed_seconds",
        "shadow_rmsd_best_pose",
        "shadow_hash_match_all",
        "n_null_features",
        "out_of_domain_flag",
    ]
    check_field_names(permitted, context="test")  # must not raise


# ---------------------------------------------------------------------------
# Allowlist tests
# ---------------------------------------------------------------------------

def test_check_allowlist_rejects_unknown_field():
    """Fields not in feature_schema.yaml are rejected regardless of content."""
    with pytest.raises(LeakageError):
        check_allowlist(["totally_unknown_field"], context="test")


def test_check_allowlist_accepts_schema_names():
    """Every name in the schema allowlist passes the allowlist check."""
    from telemetry.leakage_guard import get_allowed_names
    names = list(get_allowed_names())
    check_allowlist(names, context="test")  # must not raise


# ---------------------------------------------------------------------------
# validate_telemetry integration tests
# ---------------------------------------------------------------------------

def _minimal_clean_record():
    """Minimal telemetry record that passes both checks."""
    return {
        "n_missing_input_hashes": 0,
        "n_unregistered_artifacts": 0,
        "has_unexpected_lineage": False,
        "n_tool_warnings": 0,
        "n_retries": 0,
        "elapsed_seconds": 1.0,
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
        "out_of_domain_flag": False,
    }


def test_validate_telemetry_clean_record_passes():
    """A clean record with all permitted fields passes validation."""
    validate_telemetry(_minimal_clean_record(), context="test")


@pytest.mark.parametrize("forbidden_key", [
    "camp_value",
    "fault_id",
    "condition_id",
    "held_out_label",
    "contract_decision",
    "assay_result",
])
def test_validate_telemetry_rejects_forbidden_key(forbidden_key):
    """validate_telemetry raises LeakageError when a forbidden key is present."""
    record = _minimal_clean_record()
    record[forbidden_key] = 0
    with pytest.raises(LeakageError):
        validate_telemetry(record, context="test")


def test_validate_telemetry_rejects_unknown_key():
    """validate_telemetry raises LeakageError for an unrecognised non-forbidden key."""
    record = _minimal_clean_record()
    record["my_custom_feature"] = 42
    with pytest.raises(LeakageError):
        validate_telemetry(record, context="test")


# ---------------------------------------------------------------------------
# Feature extractor integration tests
# ---------------------------------------------------------------------------

def test_extractor_output_passes_leakage_guard():
    """The extractor's output passes validate_telemetry without raising."""
    from telemetry.feature_extractor import extract
    record = extract(
        registered_inputs=[],
        registered_outputs=[],
        tool_logs=[],
        n_retries=0,
        elapsed_seconds=5.0,
    )
    validate_telemetry(record, context="extractor integration test")


def test_extractor_contains_only_schema_fields():
    """The extractor's output contains exactly the allowlisted schema field names."""
    from telemetry.feature_extractor import extract
    from telemetry.leakage_guard import get_allowed_names
    record = extract(
        registered_inputs=[],
        registered_outputs=[],
        tool_logs=[],
        n_retries=0,
        elapsed_seconds=5.0,
    )
    allowed = get_allowed_names()
    extra = set(record.keys()) - allowed
    assert not extra, f"Extractor emitted fields not in schema: {extra}"


def test_extractor_identical_across_conditions():
    """
    The extractor produces identical field names for all A-D conditions.
    (Conditions differ in enforcement, not in telemetry schema.)
    """
    from telemetry.feature_extractor import extract
    shared_kwargs = dict(
        registered_inputs=[],
        registered_outputs=[],
        tool_logs=[],
        n_retries=0,
        elapsed_seconds=5.0,
    )
    r_a = extract(**shared_kwargs)  # condition A
    r_b = extract(**shared_kwargs)  # condition B
    r_c = extract(**shared_kwargs)  # condition C
    r_d = extract(**shared_kwargs)  # condition D
    assert set(r_a.keys()) == set(r_b.keys()) == set(r_c.keys()) == set(r_d.keys())


def test_extractor_out_of_domain_flag_is_none():
    """
    The extractor emits out_of_domain_flag=None (unresolved).
    A separate applicability-domain step must set it before inference.
    """
    from telemetry.feature_extractor import extract
    record = extract(
        registered_inputs=[],
        registered_outputs=[],
        tool_logs=[],
        n_retries=0,
        elapsed_seconds=5.0,
    )
    assert record["out_of_domain_flag"] is None, (
        "out_of_domain_flag must be None (unresolved) from the extractor; "
        "it is set by a separate domain-check step, not by the extractor."
    )


# ---------------------------------------------------------------------------
# validate_record tests
# ---------------------------------------------------------------------------

def _full_clean_record():
    """A complete telemetry record with all 24 fields, all values valid."""
    return {
        "n_missing_input_hashes": 0,
        "n_unregistered_artifacts": 0,
        "has_unexpected_lineage": False,
        "n_tool_warnings": 0,
        "n_retries": 0,
        "elapsed_seconds": 1.0,
        "n_nonzero_exit_codes": 0,
        "n_missing_outputs": 0,
        "n_schema_drift_fields": 0,
        "n_count_discrepancies": 0,
        "n_duplicate_ids": 0,
        "n_lost_records": 0,
        "n_unexpected_records": 0,
        "timestamp_staleness_seconds": 0.0,
        "pose_failure_fraction": None,      # nullable — None OK
        "score_distribution_shift": None,   # nullable
        "n_steric_clashes": None,           # nullable
        "convergence_flag": None,           # nullable
        "shadow_rmsd_best_pose": None,      # nullable
        "shadow_score_rank_rho": None,      # nullable
        "shadow_id_match_fraction": 1.0,
        "shadow_hash_match_all": True,
        "n_null_features": 6,
        "out_of_domain_flag": None,         # nullable (None = unresolved)
    }


def test_validate_record_clean_passes():
    """A complete, valid record passes validate_record."""
    validate_record(_full_clean_record(), context="test")


def test_validate_record_requires_all_fields():
    """validate_record raises LeakageError when a required field is absent."""
    rec = _full_clean_record()
    del rec["elapsed_seconds"]
    with pytest.raises(LeakageError, match="elapsed_seconds"):
        validate_record(rec, context="test_missing_key")


def test_validate_record_rejects_non_nullable_none():
    """validate_record raises LeakageError when a non-nullable field is None."""
    rec = _full_clean_record()
    rec["elapsed_seconds"] = None   # elapsed_seconds is non-nullable
    with pytest.raises(LeakageError):
        validate_record(rec, context="test_null_non_nullable")


def test_validate_record_allows_nullable_none():
    """validate_record allows None for explicitly nullable fields."""
    rec = _full_clean_record()
    rec["pose_failure_fraction"] = None     # nullable=true → allowed
    rec["out_of_domain_flag"] = None        # nullable=true → allowed
    validate_record(rec, context="test_nullable_none")


def test_validate_record_rejects_wrong_int_dtype():
    """validate_record raises LeakageError when an int field receives a float."""
    rec = _full_clean_record()
    rec["n_tool_warnings"] = 1.5  # must be int
    with pytest.raises(LeakageError, match="int"):
        validate_record(rec, context="test_dtype_int")


def test_validate_record_rejects_bool_as_int():
    """bool is a subclass of int in Python; validate_record rejects it for int fields."""
    rec = _full_clean_record()
    rec["n_tool_warnings"] = True  # bool not accepted for int fields
    with pytest.raises(LeakageError, match="int"):
        validate_record(rec, context="test_bool_as_int")


def test_validate_record_rejects_wrong_float_dtype():
    """validate_record raises LeakageError when a float field receives a string."""
    rec = _full_clean_record()
    rec["elapsed_seconds"] = "10.0"  # must be numeric
    with pytest.raises(LeakageError):
        validate_record(rec, context="test_dtype_float")


def test_validate_record_rejects_nonfinite_float():
    """Non-finite float values (inf, nan) are rejected for non-nullable float fields."""
    import math
    rec = _full_clean_record()
    rec["elapsed_seconds"] = float("inf")
    with pytest.raises(LeakageError, match="non.finite"):
        validate_record(rec, context="test_nonfinite_inf")

    rec2 = _full_clean_record()
    rec2["elapsed_seconds"] = float("nan")
    with pytest.raises(LeakageError, match="non.finite"):
        validate_record(rec2, context="test_nonfinite_nan")


def test_validate_record_rejects_fraction_exceeding_one():
    """pose_failure_fraction above 1.0 violates range_max=1.0."""
    rec = _full_clean_record()
    rec["pose_failure_fraction"] = 1.5   # range is [0.0, 1.0]
    with pytest.raises(LeakageError, match="range_max"):
        validate_record(rec, context="test_range_fraction")


def test_validate_record_rejects_negative_count():
    """Negative count fields violate range_min=0."""
    rec = _full_clean_record()
    rec["n_tool_warnings"] = -1   # range_min=0
    with pytest.raises(LeakageError, match="range_min"):
        validate_record(rec, context="test_range_negative_count")


def test_validate_record_rejects_correlation_out_of_range():
    """shadow_score_rank_rho outside [-1.0, 1.0] violates range bounds."""
    rec = _full_clean_record()
    rec["shadow_score_rank_rho"] = 1.5   # range_max=1.0
    with pytest.raises(LeakageError, match="range_max"):
        validate_record(rec, context="test_range_correlation_high")

    rec2 = _full_clean_record()
    rec2["shadow_score_rank_rho"] = -1.5  # range_min=-1.0
    with pytest.raises(LeakageError, match="range_min"):
        validate_record(rec2, context="test_range_correlation_low")


def test_validate_record_rejects_negative_elapsed():
    """elapsed_seconds must be non-negative (range_min=0.0)."""
    rec = _full_clean_record()
    rec["elapsed_seconds"] = -1.0
    with pytest.raises(LeakageError, match="range_min"):
        validate_record(rec, context="test_range_elapsed")


# ---------------------------------------------------------------------------
# get_canonical_feature_order tests
# ---------------------------------------------------------------------------

def test_canonical_order_length():
    """The canonical feature order has exactly 24 entries."""
    order = get_canonical_feature_order()
    assert len(order) == 24, f"Expected 24 features, got {len(order)}"


def test_canonical_order_matches_allowlist():
    """The canonical order contains exactly the same names as the allowlist."""
    order = get_canonical_feature_order()
    allowed = get_allowed_names()
    assert set(order) == allowed, "Canonical order and allowlist must contain identical names."


def test_canonical_order_is_stable():
    """Two calls to get_canonical_feature_order return the same list."""
    order1 = get_canonical_feature_order()
    order2 = get_canonical_feature_order()
    assert order1 == order2


def test_canonical_order_starts_with_provenance_family():
    """The first family declared in the schema is provenance_completeness."""
    order = get_canonical_feature_order()
    assert order[0] == "n_missing_input_hashes", (
        f"First feature should be 'n_missing_input_hashes', got '{order[0]}'"
    )


def test_canonical_order_ends_with_uncertainty_family():
    """The last two features belong to the uncertainty_context family."""
    order = get_canonical_feature_order()
    assert order[-2] == "n_null_features", (
        f"Second-to-last should be 'n_null_features', got '{order[-2]}'"
    )
    assert order[-1] == "out_of_domain_flag", (
        f"Last feature should be 'out_of_domain_flag', got '{order[-1]}'"
    )

"""
feature_extractor.py — Deterministic common telemetry extractor (ADR-AIIM-020).

Produces the same pre-release feature record in all four factorial conditions A–D.
This extractor is read-only: it inspects workflow state but never writes artifacts,
triggers tool calls, or reads fault identity, condition identity, held-out labels,
post-release outcomes, or assay data.

Call extract() with the current workflow state.  The return value is a validated
dict whose keys are exactly those listed in feature_schema.yaml.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any

from .leakage_guard import validate_telemetry


def _sha256_file(path: Path) -> str | None:
    if not path.is_file():
        return None
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest().upper()


def extract(
    *,
    # Provenance: list of {path: str, expected_hash: str|None, registered: bool, lineage_complete: bool}
    registered_inputs: list[dict[str, Any]],
    # Artifact presence: list of {path: str, registered: bool}
    registered_outputs: list[dict[str, Any]],
    # Tool invocation logs: list of {exit_code: int, warnings: int}
    tool_logs: list[dict[str, Any]],
    n_retries: int,
    elapsed_seconds: float,
    # Identifier sets for cross-stage consistency
    reference_id_set: set[str] | None = None,
    output_id_set: set[str] | None = None,
    # Output file timestamps (epoch seconds)
    output_timestamp: float | None = None,
    latest_input_timestamp: float | None = None,
    # Schema drift: expected vs. observed field names
    frozen_schema_fields: set[str] | None = None,
    output_schema_fields: set[str] | None = None,
    # List of (expected_count, actual_count) pairs, one per count-check in this stage
    count_checks: list[tuple[int, int]] | None = None,
    # Stage-dependent scientific metrics; absent keys produce null features
    stage_metrics: dict[str, Any] | None = None,
    # Optional fitted OutOfDomainDetector; when provided, resolves out_of_domain_flag.
    # When None (default), out_of_domain_flag stays None and inference.py fails closed.
    ood_detector=None,
) -> dict[str, Any]:
    """
    Extract and validate a telemetry record for the current pre-release state.

    Returns a dict whose keys are exactly those in feature_schema.yaml.
    Raises LeakageError (from validate_telemetry) if the record contains any
    forbidden or unknown field name.
    """
    sm = stage_metrics or {}
    ref_ids = reference_id_set or set()
    out_ids = output_id_set or set()

    # --- Provenance completeness ---
    n_missing_input_hashes = sum(
        1
        for inp in registered_inputs
        if not inp.get("expected_hash")
        or _sha256_file(Path(inp["path"])) != inp.get("expected_hash")
    )
    n_unregistered_artifacts = sum(
        1 for out in registered_outputs if not out.get("registered", False)
    )
    has_unexpected_lineage = any(
        not inp.get("lineage_complete", True) for inp in registered_inputs
    )

    # --- Execution telemetry ---
    n_tool_warnings = sum(int(log.get("warnings", 0)) for log in tool_logs)
    n_nonzero_exit_codes = sum(
        1 for log in tool_logs if int(log.get("exit_code", 0)) != 0
    )

    # --- Artifact integrity ---
    n_missing_outputs = sum(
        1 for out in registered_outputs if not Path(out["path"]).is_file()
    )
    if frozen_schema_fields is not None and output_schema_fields is not None:
        n_schema_drift_fields = len(output_schema_fields - frozen_schema_fields)
    else:
        n_schema_drift_fields = 0

    n_count_discrepancies = sum(
        1 for (expected, actual) in (count_checks or []) if expected != actual
    )
    n_duplicate_ids = int(sm.get("n_duplicate_ids", 0))

    # --- Cross-stage consistency ---
    n_lost_records = len(ref_ids - out_ids)
    n_unexpected_records = len(out_ids - ref_ids)
    timestamp_staleness_seconds = float(
        (output_timestamp or 0.0) - (latest_input_timestamp or 0.0)
    )

    # --- Scientific summaries (stage-dependent, nullable) ---
    pose_failure_fraction = sm.get("pose_failure_fraction")
    score_distribution_shift = sm.get("score_distribution_shift")
    n_steric_clashes_raw = sm.get("n_steric_clashes")
    n_steric_clashes = int(n_steric_clashes_raw) if n_steric_clashes_raw is not None else None
    convergence_flag = sm.get("convergence_flag")

    # --- Shadow check context (continuous raw values, NOT contract decisions) ---
    shadow_rmsd_best_pose = sm.get("shadow_rmsd_best_pose")
    shadow_score_rank_rho = sm.get("shadow_score_rank_rho")
    shadow_id_match_fraction = (
        len(ref_ids & out_ids) / max(len(ref_ids), 1) if ref_ids else 1.0
    )
    shadow_hash_match_all = bool(n_missing_input_hashes == 0)

    record: dict[str, Any] = {
        # provenance_completeness
        "n_missing_input_hashes": n_missing_input_hashes,
        "n_unregistered_artifacts": n_unregistered_artifacts,
        "has_unexpected_lineage": has_unexpected_lineage,
        # execution_telemetry
        "n_tool_warnings": n_tool_warnings,
        "n_retries": n_retries,
        "elapsed_seconds": float(elapsed_seconds),
        "n_nonzero_exit_codes": n_nonzero_exit_codes,
        # artifact_integrity
        "n_missing_outputs": n_missing_outputs,
        "n_schema_drift_fields": n_schema_drift_fields,
        "n_count_discrepancies": n_count_discrepancies,
        "n_duplicate_ids": n_duplicate_ids,
        # cross_stage_consistency
        "n_lost_records": n_lost_records,
        "n_unexpected_records": n_unexpected_records,
        "timestamp_staleness_seconds": timestamp_staleness_seconds,
        # scientific_summaries (nullable)
        "pose_failure_fraction": pose_failure_fraction,
        "score_distribution_shift": score_distribution_shift,
        "n_steric_clashes": n_steric_clashes,
        "convergence_flag": convergence_flag,
        # shadow_check_context
        "shadow_rmsd_best_pose": shadow_rmsd_best_pose,
        "shadow_score_rank_rho": shadow_score_rank_rho,
        "shadow_id_match_fraction": float(shadow_id_match_fraction),
        "shadow_hash_match_all": shadow_hash_match_all,
        # uncertainty_context (n_null_features must be computed from the record itself)
        # None = unresolved (domain check not yet run); set to True/False by
        # a separate applicability-domain step before calling predict_risk().
        # Inference fails closed when this remains None.
        "out_of_domain_flag": None,
    }

    # Count nulls across all features already computed (out_of_domain_flag still None here).
    # n_null_features is intentionally computed BEFORE resolving out_of_domain_flag so that
    # its value matches the training distribution, where out_of_domain_flag was always None.
    record["n_null_features"] = sum(1 for v in record.values() if v is None)

    # Resolve out_of_domain_flag AFTER n_null_features is set.
    # When ood_detector is None the flag stays None; inference.py will fail closed.
    if ood_detector is not None:
        record["out_of_domain_flag"] = ood_detector.predict_one(record)

    validate_telemetry(record, context="feature_extractor.extract()")
    return record

"""
Append-only provenance log for all artifacts, contracts, and human interventions.

Written as JSON Lines (one record per line) to a file in the run directory.
The registry is never modified after a record is written — only appended to.
"""

from __future__ import annotations

import hashlib
import json
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest().upper()


@dataclass
class ArtifactRecord:
    record_type: str          # "artifact", "contract_outcome", "human_intervention", "access_log"
    record_id: str            # Unique within this registry
    timestamp: float = field(default_factory=time.time)
    data: dict[str, Any] = field(default_factory=dict)


@dataclass
class ArtifactEntry:
    artifact_id: str
    path: str
    sha256: str
    stage: str
    description: str
    tool: str
    tool_version: str
    conda_env: str
    command: str
    random_seed: int | None
    timestamp: float = field(default_factory=time.time)
    resource_use: dict[str, Any] = field(default_factory=dict)
    agent_decision_summary: str = ""
    extra: dict[str, Any] = field(default_factory=dict)


@dataclass
class ContractOutcomeEntry:
    contract_id: str
    gate: str
    artifact_ids: list[str]
    outcome: str                  # AuditState value
    checks_passed: list[str]
    checks_failed: list[str]
    checks_warned: list[str]
    state_history: list[dict[str, Any]]
    timestamp: float = field(default_factory=time.time)
    human_adjudicated: bool = False
    adjudication_note: str = ""


@dataclass
class HumanInterventionEntry:
    intervention_id: str
    gate: str
    information_visible_to_human: str
    decision: str
    prespecified: bool
    effect_on_result: str
    # Artifacts this adjudication covers. ReportLock releases an artifact only
    # when it has a passed contract OR appears here, so an intervention that
    # names no artifact cannot unlock anything.
    artifact_ids: list[str] = field(default_factory=list)
    timestamp: float = field(default_factory=time.time)


@dataclass
class FileAccessEntry:
    """Records which files the orchestrator or evaluators READ — enables label-isolation audit."""
    reader: str           # "orchestrator", "contract_engine", or evaluator module name
    file_path: str
    purpose: str
    timestamp: float = field(default_factory=time.time)


class ArtifactRegistry:
    def __init__(self, registry_path: Path) -> None:
        self._path = registry_path
        self._path.parent.mkdir(parents=True, exist_ok=True)
        # Touch the file if new; never truncate an existing one.
        if not self._path.exists():
            self._path.touch()

    def _append(self, record_type: str, record_id: str, data: dict[str, Any]) -> None:
        record = {"record_type": record_type, "record_id": record_id,
                  "timestamp": time.time(), **data}
        with open(self._path, "a", encoding="utf-8") as f:
            f.write(json.dumps(record, ensure_ascii=False) + "\n")

    def register_artifact(self, entry: ArtifactEntry) -> None:
        self._append("artifact", entry.artifact_id, asdict(entry))

    def record_contract_outcome(self, entry: ContractOutcomeEntry) -> None:
        self._append("contract_outcome", entry.contract_id, asdict(entry))

    def record_human_intervention(self, entry: HumanInterventionEntry) -> None:
        self._append("human_intervention", entry.intervention_id, asdict(entry))

    def record_file_access(self, entry: FileAccessEntry) -> None:
        self._append("access_log", f"access_{int(time.time()*1e6)}", asdict(entry))

    def record_benchmark_role_declaration(
        self,
        role: str,          # "Option_A_development" or "Option_B_heldout"
        declared_by: str,   # human identifier
        rationale: str,
    ) -> None:
        self._append(
            "benchmark_role_declaration",
            f"benchmark_role_{int(time.time())}",
            {"role": role, "declared_by": declared_by, "rationale": rationale},
        )

    def get_all_records(self) -> list[dict[str, Any]]:
        records = []
        with open(self._path, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line:
                    records.append(json.loads(line))
        return records

    def get_benchmark_role_declaration(self) -> dict[str, Any] | None:
        for rec in self.get_all_records():
            if rec["record_type"] == "benchmark_role_declaration":
                return rec
        return None

    def get_file_access_log(self) -> list[dict[str, Any]]:
        return [r for r in self.get_all_records() if r["record_type"] == "access_log"]

    def count_human_interventions(self) -> int:
        return sum(1 for r in self.get_all_records() if r["record_type"] == "human_intervention")

    def artifact_has_passed_contract(self, artifact_id: str) -> bool:
        """True if the artifact has at least one passed contract record."""
        for rec in self.get_all_records():
            if (rec["record_type"] == "contract_outcome"
                    and artifact_id in rec.get("artifact_ids", [])
                    and rec.get("outcome") == "PASS"):
                return True
        return False

    # ── ML risk-monitor artifact and event records ────────────────────────────

    def register_ml_schema_artifact(
        self,
        schema_path: str,
        sha256: str,
        schema_version: str,
        approved_by: str,
        feature_count: int,
        notes: str = "",
    ) -> None:
        """Record the frozen feature schema artifact (Gate 1 approval)."""
        self._append(
            "ml_schema_artifact",
            f"ml_schema_{int(time.time()*1e6)}",
            {
                "schema_path": schema_path,
                "sha256": sha256,
                "schema_version": schema_version,
                "approved_by": approved_by,
                "feature_count": feature_count,
                "notes": notes,
            },
        )

    def register_ml_feature_order(
        self,
        feature_order: list[str],
        feature_order_sha256: str,
        schema_version: str,
        approved_by: str,
    ) -> None:
        """Record the canonical feature order artifact (must be immutable after approval)."""
        self._append(
            "ml_feature_order",
            f"ml_feature_order_{int(time.time()*1e6)}",
            {
                "feature_order": feature_order,
                "feature_order_sha256": feature_order_sha256,
                "schema_version": schema_version,
                "approved_by": approved_by,
                "n_features": len(feature_order),
            },
        )

    def register_ml_model_artifact(
        self,
        artifact_id: str,
        model_path: str,
        sha256: str,
        schema_version: str,
        approved_by: str,
        model_name: str,
        mean_brier: float,
        feature_order_sha256: str = "",
        notes: str = "",
    ) -> None:
        """
        Record a trained model artifact.  approved_by must be non-empty (human approval).
        This mirrors the ModelRegistry entry but in the provenance audit trail.
        """
        if not approved_by or not approved_by.strip():
            raise ValueError(
                "approved_by must be non-empty for ML model registration "
                "(human approval required, ADR-AIIM-021)."
            )
        self._append(
            "ml_model_artifact",
            f"ml_model_{artifact_id}_{int(time.time()*1e6)}",
            {
                "artifact_id": artifact_id,
                "model_path": model_path,
                "sha256": sha256,
                "schema_version": schema_version,
                "approved_by": approved_by,
                "model_name": model_name,
                "mean_brier": mean_brier,
                "feature_order_sha256": feature_order_sha256,
                "notes": notes,
            },
        )

    def register_ml_calibration_artifact(
        self,
        artifact_id: str,
        calibrated_model_path: str,
        sha256: str,
        base_model_artifact_id: str,
        brier_score: float,
        method: str,
        n_cal_folds: int,
        approved_by: str,
        notes: str = "",
    ) -> None:
        """Record a calibrated model artifact."""
        self._append(
            "ml_calibration_artifact",
            f"ml_cal_{artifact_id}_{int(time.time()*1e6)}",
            {
                "artifact_id": artifact_id,
                "calibrated_model_path": calibrated_model_path,
                "sha256": sha256,
                "base_model_artifact_id": base_model_artifact_id,
                "brier_score": brier_score,
                "method": method,
                "n_cal_folds": n_cal_folds,
                "approved_by": approved_by,
                "notes": notes,
            },
        )

    def record_ml_prediction_event(
        self,
        artifact_id: str,
        risk_score: float,
        action: str,
        reason: str,
        schema_version_match: bool,
        out_of_domain_flag: bool | None,
        workflow_stage: str,
    ) -> None:
        """Record a single risk-monitor inference event for audit purposes."""
        self._append(
            "ml_prediction_event",
            f"ml_pred_{int(time.time()*1e6)}",
            {
                "artifact_id": artifact_id,
                "risk_score": risk_score,
                "action": action,
                "reason": reason,
                "schema_version_match": schema_version_match,
                "out_of_domain_flag": out_of_domain_flag,
                "workflow_stage": workflow_stage,
            },
        )

    def get_ml_records(self, record_type: str | None = None) -> list[dict[str, Any]]:
        """Return all ML-specific audit records, optionally filtered by record_type."""
        ml_types = {
            "ml_schema_artifact", "ml_feature_order", "ml_model_artifact",
            "ml_calibration_artifact", "ml_prediction_event",
        }
        return [
            r for r in self.get_all_records()
            if r["record_type"] in ml_types
            and (record_type is None or r["record_type"] == record_type)
        ]

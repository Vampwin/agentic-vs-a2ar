"""
Workflow orchestrator — drives the virtual screening pipeline stage by stage.

The orchestrator maintains workflow state, chooses only allowed tools and
recovery actions, and records its reasoning. It cannot bypass a failed contract,
alter a frozen threshold, or read held-out fault labels.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .artifact_registry import ArtifactEntry, ArtifactRegistry, FileAccessEntry, sha256_file
from .contract_engine import ContractEngine, ContractResult
from .report_lock import ReportLock
from .state_machine import AuditState

logger = logging.getLogger(__name__)


@dataclass
class StageResult:
    stage: str
    contract_result: ContractResult
    artifacts_produced: list[str] = field(default_factory=list)


class WorkflowOrchestrator:
    """
    Coordinates the G0→G8 pipeline. Each stage:
    1. Invokes the deterministic domain tool.
    2. Registers resulting artifacts.
    3. Runs the corresponding contract gate.
    4. Proceeds, repairs, escalates, or halts based on the contract outcome.

    The orchestrator records all file reads via the registry access log,
    enabling the G5a label-isolation check.
    """

    def __init__(
        self,
        run_dir: Path,
        contracts_dir: Path,
        config_path: Path,
        registry: ArtifactRegistry,
    ) -> None:
        self.run_dir = run_dir
        self.engine = ContractEngine(contracts_dir=contracts_dir, registry=registry)
        self.registry = registry
        self.config_path = config_path
        self.stage_results: list[StageResult] = []
        self._halted = False
        self._escalation_pending = False

    # ── Public interface ──────────────────────────────────────────────────────

    def run_stage(
        self,
        stage_name: str,
        contract_file: str,
        tool_fn: Any,           # callable: (run_dir, config, registry) -> dict
        tool_outputs: dict[str, Any] | None = None,
    ) -> StageResult:
        """
        Run one pipeline stage. tool_fn is a domain tool (docking, standardization,
        etc.) that the orchestrator calls deterministically. tool_outputs, if
        provided, skips calling tool_fn (for testing or re-entry).
        """
        if self._halted:
            raise RuntimeError(f"Workflow halted before {stage_name}; cannot proceed.")

        logger.info("Starting stage: %s", stage_name)

        if tool_outputs is None:
            tool_outputs = tool_fn(
                run_dir=self.run_dir,
                config_path=self.config_path,
                registry=self.registry,
            )

        artifact_ids = self._register_outputs(stage_name, tool_outputs)
        contract_result = self.engine.run(
            contract_file=contract_file,
            inputs=tool_outputs,
            artifact_ids=artifact_ids,
        )

        stage_result = StageResult(
            stage=stage_name,
            contract_result=contract_result,
            artifacts_produced=artifact_ids,
        )
        self.stage_results.append(stage_result)

        self._handle_outcome(stage_name, contract_result)
        return stage_result

    def record_read(self, reader: str, file_path: str, purpose: str) -> None:
        """
        Call this whenever the orchestrator or a domain tool reads a file.
        Used by the G5a evaluator to verify label isolation.
        """
        self.registry.record_file_access(
            FileAccessEntry(reader=reader, file_path=file_path, purpose=purpose)
        )

    def request_human_adjudication(
        self,
        gate: str,
        information_visible: str,
        decision_needed: str,
    ) -> str:
        """
        Pause the workflow and request human input. Returns the human's decision
        string (provided externally). The caller must record the intervention via
        registry.record_human_intervention before calling run_stage again.

        This method does NOT implement an interactive prompt — it raises an
        exception that the runner script catches to surface to the human.
        """
        raise EscalationRequired(
            gate=gate,
            information_visible=information_visible,
            decision_needed=decision_needed,
        )

    def build_report_lock(self, lock_file: Path) -> ReportLock:
        return ReportLock(registry=self.registry, lock_file=lock_file)

    # ── Internal ──────────────────────────────────────────────────────────────

    def _register_outputs(self, stage: str, tool_outputs: dict[str, Any]) -> list[str]:
        artifact_ids = []
        for key, value in tool_outputs.items():
            if not isinstance(value, (str, Path)):
                continue
            p = Path(value)
            if not p.exists():
                continue
            artifact_id = f"{stage}::{key}"
            entry = ArtifactEntry(
                artifact_id=artifact_id,
                path=str(p),
                sha256=sha256_file(p),
                stage=stage,
                description=key,
                tool="",          # Populated by domain tool registration
                tool_version="",
                conda_env="",
                command="",
                random_seed=None,
            )
            self.registry.register_artifact(entry)
            artifact_ids.append(artifact_id)
        return artifact_ids

    def _handle_outcome(self, stage: str, result: ContractResult) -> None:
        state = result.state
        if state == AuditState.PASS:
            logger.info("%s: PASS", stage)
        elif state == AuditState.SKIPPED_BY_PROTOCOL:
            logger.info("%s: SKIPPED_BY_PROTOCOL", stage)
        elif state == AuditState.AUTO_REPAIR:
            logger.warning("%s: AUTO_REPAIR triggered", stage)
            # Repair was attempted inside ContractEngine; re-check state
            if result.state != AuditState.PASS:
                self._escalation_pending = True
        elif state == AuditState.ESCALATE:
            logger.warning("%s: ESCALATE — human adjudication required", stage)
            self._escalation_pending = True
        elif state == AuditState.HALT:
            logger.error("%s: HALT — workflow stopped", stage)
            self._halted = True


class EscalationRequired(Exception):
    def __init__(self, gate: str, information_visible: str, decision_needed: str) -> None:
        self.gate = gate
        self.information_visible = information_visible
        self.decision_needed = decision_needed
        super().__init__(
            f"Human adjudication required at gate {gate}: {decision_needed}"
        )

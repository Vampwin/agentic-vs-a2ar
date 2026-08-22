"""
Loads versioned YAML contract definitions, invokes the specified evaluator,
and records outcomes in the artifact registry.

The contract engine runs outside the acting agent's control: it is invoked by
the orchestrator but cannot be instructed by the acting agent to skip checks
or relax thresholds.
"""

from __future__ import annotations

import importlib
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from .artifact_registry import ArtifactRegistry, ContractOutcomeEntry
from .state_machine import AuditState, StateMachine


@dataclass
class ContractResult:
    gate: str
    contract_id: str
    state: AuditState
    checks_passed: list[str]
    checks_failed: list[str]
    checks_warned: list[str]
    evidence: dict[str, Any]
    state_history: list[dict[str, Any]]
    elapsed_seconds: float
    error: str | None = None


class ContractEngine:
    def __init__(
        self,
        contracts_dir: Path,
        registry: ArtifactRegistry,
    ) -> None:
        self._contracts_dir = contracts_dir
        self._registry = registry

    def load_contract(self, contract_file: str) -> dict[str, Any]:
        path = self._contracts_dir / contract_file
        with open(path, encoding="utf-8") as f:
            return yaml.safe_load(f)

    def run(
        self,
        contract_file: str,
        inputs: dict[str, Any],
        artifact_ids: list[str],
    ) -> ContractResult:
        t0 = time.time()
        contract = self.load_contract(contract_file)
        gate = contract["contract_id"]
        evaluator_module_path = contract["evaluator"]  # e.g. "evaluators.G0_evaluator"
        sm = StateMachine(gate_id=gate)

        try:
            mod = importlib.import_module(evaluator_module_path)
            evaluator_fn = getattr(mod, "evaluate")
            eval_result = evaluator_fn(contract=contract, inputs=inputs, registry=self._registry)
        except Exception as exc:
            sm.set_initial(
                AuditState.HALT,
                reason=f"Evaluator raised exception: {exc}",
                evidence=["evaluator_exception"],
            )
            result = ContractResult(
                gate=gate,
                contract_id=gate,
                state=AuditState.HALT,
                checks_passed=[],
                checks_failed=[],
                checks_warned=[],
                evidence={},
                state_history=sm.full_history(),
                elapsed_seconds=time.time() - t0,
                error=str(exc),
            )
            self._record(result, contract, artifact_ids)
            return result

        passed = eval_result.get("passed", [])
        failed = eval_result.get("failed", [])
        warned = eval_result.get("warned", [])
        evidence = eval_result.get("evidence", {})

        # Determine initial state from evaluator output
        initial_state_str = eval_result.get("state", None)
        if initial_state_str:
            initial_state = AuditState(initial_state_str)
        elif failed:
            # Default: any critical failure → HALT; any non-critical → ESCALATE
            critical_failures = eval_result.get("critical_failures", failed)
            initial_state = AuditState.HALT if critical_failures else AuditState.ESCALATE
        else:
            initial_state = AuditState.PASS

        sm.set_initial(
            initial_state,
            reason=eval_result.get("summary", "Evaluator completed"),
            evidence=list(evidence.keys()),
        )

        # Auto-repair if the contract specifies it and state is AUTO_REPAIR
        if initial_state == AuditState.AUTO_REPAIR:
            repair_result = self._attempt_auto_repair(contract, inputs, sm)
            passed = repair_result.get("passed", passed)
            failed = repair_result.get("failed", failed)

        result = ContractResult(
            gate=gate,
            contract_id=gate,
            state=sm.current,
            checks_passed=passed,
            checks_failed=failed,
            checks_warned=warned,
            evidence=evidence,
            state_history=sm.full_history(),
            elapsed_seconds=time.time() - t0,
        )
        self._record(result, contract, artifact_ids)
        return result

    def _attempt_auto_repair(
        self,
        contract: dict[str, Any],
        inputs: dict[str, Any],
        sm: StateMachine,
    ) -> dict[str, Any]:
        """
        Only executes prespecified, reversible repair actions from the contract.
        Cannot alter thresholds, change evaluation sets, or reinterpret failures.
        """
        repair_actions = contract.get("on_fail", [])
        # Only run if "rerun_with_same_frozen_parameters" is listed
        if "rerun_with_same_frozen_parameters" not in repair_actions:
            sm.transition(
                AuditState.ESCALATE,
                reason="AUTO_REPAIR requested but no reversible action available",
                evidence=["no_reversible_repair"],
            )
            return {}
        # Repair logic is evaluator-specific; signal the evaluator to rerun
        evaluator_module_path = contract["evaluator"]
        mod = importlib.import_module(evaluator_module_path)
        repair_fn = getattr(mod, "repair", None)
        if repair_fn is None:
            sm.transition(
                AuditState.ESCALATE,
                reason="Evaluator has no repair() function",
                evidence=["no_repair_function"],
            )
            return {}
        repair_result = repair_fn(contract=contract, inputs=inputs)
        if repair_result.get("success"):
            sm.transition(
                AuditState.PASS,
                reason="Auto-repair succeeded",
                evidence=repair_result.get("evidence_keys", []),
            )
        else:
            sm.transition(
                AuditState.ESCALATE,
                reason=f"Auto-repair failed: {repair_result.get('reason', 'unknown')}",
                evidence=repair_result.get("evidence_keys", []),
            )
        return repair_result

    def _record(
        self,
        result: ContractResult,
        contract: dict[str, Any],
        artifact_ids: list[str],
    ) -> None:
        entry = ContractOutcomeEntry(
            contract_id=result.contract_id,
            gate=result.gate,
            artifact_ids=artifact_ids,
            outcome=result.state.value,
            checks_passed=result.checks_passed,
            checks_failed=result.checks_failed,
            checks_warned=result.checks_warned,
            state_history=result.state_history,
        )
        self._registry.record_contract_outcome(entry)

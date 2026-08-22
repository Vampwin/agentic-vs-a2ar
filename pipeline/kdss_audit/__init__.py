"""
kdss_audit — contract engine and audit infrastructure for the AIIM agentic workflow.

Components:
    state_machine      — PASS/AUTO_REPAIR/ESCALATE/HALT transitions
    artifact_registry  — append-only provenance log (SHA-256, seeds, versions)
    contract_engine    — loads YAML contracts, invokes evaluators, records outcomes
    report_lock        — releases report artifacts only from passed-contract runs
    orchestrator       — drives the full workflow stage by stage
"""

from .state_machine import AuditState, StateTransition, StateMachine
from .artifact_registry import ArtifactRegistry
from .contract_engine import ContractEngine, ContractResult
from .report_lock import ReportLock
from .orchestrator import WorkflowOrchestrator

__all__ = [
    "AuditState",
    "StateTransition",
    "StateMachine",
    "ArtifactRegistry",
    "ContractEngine",
    "ContractResult",
    "ReportLock",
    "WorkflowOrchestrator",
]

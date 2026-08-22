"""
Four-state machine for audit gate outcomes.

Allowed states: PASS, AUTO_REPAIR, ESCALATE, HALT.
Every transition carries a machine-readable reason and links to evidence.
The acting agent cannot trigger a PASS — only the contract engine can emit it.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from enum import Enum
from typing import Any


class AuditState(str, Enum):
    PASS = "PASS"
    AUTO_REPAIR = "AUTO_REPAIR"
    ESCALATE = "ESCALATE"
    HALT = "HALT"
    SKIPPED_BY_PROTOCOL = "SKIPPED_BY_PROTOCOL"  # G7 when MD not run


# Which transitions are allowed per current state.
# HALT and SKIPPED_BY_PROTOCOL are terminal — no transitions out.
_ALLOWED_TRANSITIONS: dict[AuditState, set[AuditState]] = {
    AuditState.PASS: set(),
    AuditState.AUTO_REPAIR: {AuditState.PASS, AuditState.ESCALATE, AuditState.HALT},
    AuditState.ESCALATE: {AuditState.PASS, AuditState.HALT},
    AuditState.HALT: set(),
    AuditState.SKIPPED_BY_PROTOCOL: set(),
}


@dataclass
class StateTransition:
    from_state: AuditState
    to_state: AuditState
    reason: str
    evidence: list[str] = field(default_factory=list)
    timestamp: float = field(default_factory=time.time)
    human_initiated: bool = False

    def as_dict(self) -> dict[str, Any]:
        return {
            "from_state": self.from_state.value,
            "to_state": self.to_state.value,
            "reason": self.reason,
            "evidence": self.evidence,
            "timestamp": self.timestamp,
            "human_initiated": self.human_initiated,
        }


class StateMachine:
    def __init__(self, gate_id: str) -> None:
        self.gate_id = gate_id
        self.history: list[StateTransition] = []
        self._current: AuditState | None = None

    @property
    def current(self) -> AuditState | None:
        return self._current

    def set_initial(self, state: AuditState, reason: str, evidence: list[str] | None = None) -> None:
        if self._current is not None:
            raise RuntimeError(f"{self.gate_id}: initial state already set")
        self._current = state
        self.history.append(
            StateTransition(
                from_state=state,  # no prior state; from == to marks the initial entry
                to_state=state,
                reason=reason,
                evidence=evidence or [],
            )
        )

    def transition(
        self,
        to_state: AuditState,
        reason: str,
        evidence: list[str] | None = None,
        human_initiated: bool = False,
    ) -> AuditState:
        if self._current is None:
            raise RuntimeError(f"{self.gate_id}: call set_initial before transition")
        allowed = _ALLOWED_TRANSITIONS.get(self._current, set())
        if to_state not in allowed:
            raise ValueError(
                f"{self.gate_id}: transition {self._current} → {to_state} is not allowed. "
                f"Allowed: {allowed}"
            )
        t = StateTransition(
            from_state=self._current,
            to_state=to_state,
            reason=reason,
            evidence=evidence or [],
            human_initiated=human_initiated,
        )
        self.history.append(t)
        self._current = to_state
        return self._current

    def is_terminal(self) -> bool:
        return self._current in (AuditState.HALT, AuditState.SKIPPED_BY_PROTOCOL)

    def full_history(self) -> list[dict[str, Any]]:
        return [t.as_dict() for t in self.history]

"""
policy.py — Frozen selective risk policy (ADR-AIIM-018, 14_ML_RISK_CONTROL_PROTOCOL.md).

Two probability thresholds divide risk into three action regions:

    risk < t_escalate          → CONTINUE (only if all required contracts passed)
    t_escalate <= risk < t_halt → ESCALATE (human adjudication)
    risk >= t_halt              → HALT     (safe halt or prespecified reversible repair)

Thresholds are selected using development/calibration data to satisfy prespecified
clean-case escalation and high-severity miss-risk budgets.  They require human
approval before the held-out run and are frozen in a separate config that can be
registered and hashed.  This file provides the action logic only; it never stores
the thresholds itself so that threshold values cannot be silently changed by editing
this file.

Hard contracts have ABSOLUTE precedence (ADR-AIIM-018):
- A HALT from a failed contract cannot be downgraded to CONTINUE by any risk score.
- A low risk score cannot release an artifact blocked by report_lock.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class PolicyAction(str, Enum):
    CONTINUE = "CONTINUE"
    ESCALATE = "ESCALATE"
    HALT = "HALT"


@dataclass(frozen=True)
class PolicyResult:
    action: PolicyAction
    risk: float
    t_escalate: float
    t_halt: float
    reason: str


def evaluate_policy(
    risk: float,
    t_escalate: float,
    t_halt: float,
    contracts_all_passed: bool,
) -> PolicyResult:
    """
    Map a calibrated risk probability to a policy action.

    Parameters
    ----------
    risk:
        Calibrated probability that the current pre-release state is unsafe to release.
    t_escalate, t_halt:
        Frozen probability thresholds.  Must satisfy 0 < t_escalate <= t_halt <= 1.
    contracts_all_passed:
        True only if every required hard contract for this stage passed.
        If False, the result is HALT regardless of risk score (ADR-AIIM-018).
    """
    if not (0.0 < t_escalate <= t_halt <= 1.0):
        raise ValueError(
            f"Invalid thresholds: t_escalate={t_escalate}, t_halt={t_halt}. "
            "Must satisfy 0 < t_escalate <= t_halt <= 1."
        )
    if not (0.0 <= risk <= 1.0):
        raise ValueError(f"risk={risk!r} is not a probability in [0, 1]")

    # Hard-contract precedence — ADR-AIIM-018
    if not contracts_all_passed:
        return PolicyResult(
            action=PolicyAction.HALT,
            risk=risk,
            t_escalate=t_escalate,
            t_halt=t_halt,
            reason="Hard contract failed — contract precedence overrides model (ADR-AIIM-018)",
        )

    if risk >= t_halt:
        reason = f"risk={risk:.4f} >= t_halt={t_halt:.4f}"
        action = PolicyAction.HALT
    elif risk >= t_escalate:
        reason = f"t_escalate={t_escalate:.4f} <= risk={risk:.4f} < t_halt={t_halt:.4f}"
        action = PolicyAction.ESCALATE
    else:
        reason = f"risk={risk:.4f} < t_escalate={t_escalate:.4f}"
        action = PolicyAction.CONTINUE

    return PolicyResult(
        action=action,
        risk=risk,
        t_escalate=t_escalate,
        t_halt=t_halt,
        reason=reason,
    )

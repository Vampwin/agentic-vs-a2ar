"""
test_contract_precedence.py — Hard-contract precedence over model predictions (ADR-AIIM-018).

Verifies that:
1. A failed hard contract produces HALT regardless of model risk score.
2. A low risk score does not downgrade a contract HALT to CONTINUE.
3. The model cannot change a contract FAIL to PASS.
4. evaluate_policy enforces contract precedence via contracts_all_passed=False.
5. Inference from a trusted model is still overridden by a contract failure.
"""

import pytest
from risk_monitor.policy import PolicyAction, evaluate_policy


# ---------------------------------------------------------------------------
# Policy-level contract precedence tests
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("risk", [0.0, 0.001, 0.05, 0.4, 0.99, 1.0])
def test_failed_contract_always_halts(risk):
    """Any risk score produces HALT when contracts_all_passed=False (ADR-AIIM-018)."""
    result = evaluate_policy(
        risk=risk,
        t_escalate=0.3,
        t_halt=0.7,
        contracts_all_passed=False,
    )
    assert result.action == PolicyAction.HALT, (
        f"risk={risk}: expected HALT due to failed contract, got {result.action}"
    )
    assert "contract" in result.reason.lower(), (
        "Policy reason should mention the contract precedence rule"
    )


def test_low_risk_cannot_release_blocked_artifact():
    """
    A very low risk score must not produce CONTINUE when contracts failed.
    This directly tests ADR-AIIM-018: 'A low-risk prediction cannot release
    an artifact blocked by report_lock.'
    """
    result = evaluate_policy(
        risk=0.001,      # near-zero risk
        t_escalate=0.3,
        t_halt=0.7,
        contracts_all_passed=False,  # contract blocked the artifact
    )
    assert result.action == PolicyAction.HALT
    assert result.action != PolicyAction.CONTINUE


def test_high_risk_with_passing_contracts_halts():
    """When contracts pass but risk is high, action is HALT from risk policy."""
    result = evaluate_policy(
        risk=0.95,
        t_escalate=0.3,
        t_halt=0.7,
        contracts_all_passed=True,
    )
    assert result.action == PolicyAction.HALT


def test_mid_risk_with_passing_contracts_escalates():
    """When contracts pass and risk is mid-range, action is ESCALATE."""
    result = evaluate_policy(
        risk=0.5,
        t_escalate=0.3,
        t_halt=0.7,
        contracts_all_passed=True,
    )
    assert result.action == PolicyAction.ESCALATE


def test_low_risk_with_passing_contracts_continues():
    """When contracts pass and risk is low, action is CONTINUE."""
    result = evaluate_policy(
        risk=0.1,
        t_escalate=0.3,
        t_halt=0.7,
        contracts_all_passed=True,
    )
    assert result.action == PolicyAction.CONTINUE


def test_policy_action_at_threshold_boundary_escalate():
    """Risk exactly at t_escalate produces ESCALATE, not CONTINUE."""
    result = evaluate_policy(
        risk=0.3,
        t_escalate=0.3,
        t_halt=0.7,
        contracts_all_passed=True,
    )
    assert result.action == PolicyAction.ESCALATE


def test_policy_action_at_threshold_boundary_halt():
    """Risk exactly at t_halt produces HALT, not ESCALATE."""
    result = evaluate_policy(
        risk=0.7,
        t_escalate=0.3,
        t_halt=0.7,
        contracts_all_passed=True,
    )
    assert result.action == PolicyAction.HALT


def test_invalid_thresholds_raise():
    """Nonsensical threshold values raise ValueError before any action is taken."""
    with pytest.raises(ValueError):
        evaluate_policy(risk=0.5, t_escalate=0.8, t_halt=0.3, contracts_all_passed=True)

    with pytest.raises(ValueError):
        evaluate_policy(risk=0.5, t_escalate=0.0, t_halt=0.7, contracts_all_passed=True)


def test_invalid_risk_raises():
    """Out-of-range risk values raise ValueError."""
    with pytest.raises(ValueError):
        evaluate_policy(risk=-0.1, t_escalate=0.3, t_halt=0.7, contracts_all_passed=True)

    with pytest.raises(ValueError):
        evaluate_policy(risk=1.1, t_escalate=0.3, t_halt=0.7, contracts_all_passed=True)


# ---------------------------------------------------------------------------
# PolicyResult immutability — frozen dataclass
# ---------------------------------------------------------------------------

def test_policy_result_is_frozen():
    """PolicyResult is frozen: model output cannot be mutated after the fact."""
    result = evaluate_policy(
        risk=0.5,
        t_escalate=0.3,
        t_halt=0.7,
        contracts_all_passed=True,
    )
    with pytest.raises(Exception):
        result.action = PolicyAction.CONTINUE  # type: ignore[misc]


# ---------------------------------------------------------------------------
# Contract-precedence narrative consistency
# ---------------------------------------------------------------------------

def test_risk_stored_even_when_contract_overrides():
    """
    The risk probability is still recorded in the PolicyResult when contracts failed.
    This supports audit logging without changing the action.
    """
    result = evaluate_policy(
        risk=0.05,
        t_escalate=0.3,
        t_halt=0.7,
        contracts_all_passed=False,
    )
    assert result.action == PolicyAction.HALT
    assert result.risk == pytest.approx(0.05)  # risk value preserved for logging

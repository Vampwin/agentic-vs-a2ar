"""
Unit tests for StateMachine — verifies all allowed and forbidden transitions.
"""

import pytest

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from kdss_audit.state_machine import AuditState, StateMachine


def make_sm(initial: AuditState = AuditState.PASS) -> StateMachine:
    sm = StateMachine("test_gate")
    sm.set_initial(initial, reason="test")
    return sm


def test_initial_state_recorded():
    sm = make_sm(AuditState.PASS)
    assert sm.current == AuditState.PASS
    assert len(sm.history) == 1


def test_pass_is_terminal():
    sm = make_sm(AuditState.PASS)
    with pytest.raises(ValueError):
        sm.transition(AuditState.HALT, reason="should not happen")


def test_halt_is_terminal():
    sm = make_sm(AuditState.HALT)
    with pytest.raises(ValueError):
        sm.transition(AuditState.ESCALATE, reason="should not happen")


def test_auto_repair_to_pass():
    sm = make_sm(AuditState.AUTO_REPAIR)
    sm.transition(AuditState.PASS, reason="repair succeeded")
    assert sm.current == AuditState.PASS


def test_auto_repair_to_escalate():
    sm = make_sm(AuditState.AUTO_REPAIR)
    sm.transition(AuditState.ESCALATE, reason="repair failed")
    assert sm.current == AuditState.ESCALATE


def test_auto_repair_to_halt():
    sm = make_sm(AuditState.AUTO_REPAIR)
    sm.transition(AuditState.HALT, reason="unrecoverable")
    assert sm.current == AuditState.HALT
    assert sm.is_terminal()


def test_escalate_to_pass():
    sm = make_sm(AuditState.ESCALATE)
    sm.transition(AuditState.PASS, reason="human approved", human_initiated=True)
    assert sm.current == AuditState.PASS


def test_escalate_to_halt():
    sm = make_sm(AuditState.ESCALATE)
    sm.transition(AuditState.HALT, reason="human stopped workflow", human_initiated=True)
    assert sm.is_terminal()


def test_escalate_cannot_go_to_auto_repair():
    sm = make_sm(AuditState.ESCALATE)
    with pytest.raises(ValueError):
        sm.transition(AuditState.AUTO_REPAIR, reason="invalid")


def test_history_length():
    sm = make_sm(AuditState.AUTO_REPAIR)
    sm.transition(AuditState.ESCALATE, reason="r1")
    sm.transition(AuditState.PASS, reason="r2", human_initiated=True)
    assert len(sm.history) == 3  # initial + 2 transitions


def test_full_history_serializable():
    sm = make_sm(AuditState.AUTO_REPAIR)
    sm.transition(AuditState.PASS, reason="ok")
    h = sm.full_history()
    import json
    serialized = json.dumps(h)
    assert "PASS" in serialized


def test_double_initial_raises():
    sm = StateMachine("g")
    sm.set_initial(AuditState.PASS, reason="first")
    with pytest.raises(RuntimeError):
        sm.set_initial(AuditState.PASS, reason="second")


def test_transition_without_initial_raises():
    sm = StateMachine("g")
    with pytest.raises(RuntimeError):
        sm.transition(AuditState.PASS, reason="no initial set")


def test_skipped_by_protocol_is_terminal():
    sm = StateMachine("G7")
    sm.set_initial(AuditState.SKIPPED_BY_PROTOCOL, reason="md_decision = no_MD")
    assert sm.is_terminal()
    with pytest.raises(ValueError):
        sm.transition(AuditState.PASS, reason="invalid")

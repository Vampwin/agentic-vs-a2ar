"""test_orchestrator.py — branch coverage of the workflow orchestrator (ADR-AIIM-199).

Companion to `test_contract_engine.py`. The external audit recorded 0 of 22
branches executed here, and no test in the suite imported the module
(ADR-AIIM-198), although `scripts/run_snp_screen.py` drives the production
screen through it.

What is under test is the orchestrator's own control flow, so the contract
engine is stubbed: these tests decide what state comes back and assert what the
orchestrator does with it. The engine's own behaviour is covered separately.

Two branches are deliberately not covered, and the reason is recorded in
`test_unhandled_state_is_silently_ignored` and
`test_auto_repair_branch_is_unreachable_from_the_engine` below. They are
unreachable rather than untested, which is a different finding and is reported
as such rather than papered over with a contrived assertion.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from kdss_audit.artifact_registry import ArtifactRegistry
from kdss_audit.contract_engine import ContractResult
from kdss_audit.orchestrator import (EscalationRequired, StageResult,
                                     WorkflowOrchestrator)
from kdss_audit.report_lock import ReportLock
from kdss_audit.state_machine import AuditState


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _result(state, gate="g0", **kw):
    return ContractResult(
        gate=gate,
        contract_id=gate,
        state=state,
        checks_passed=kw.get("passed", []),
        checks_failed=kw.get("failed", []),
        checks_warned=kw.get("warned", []),
        evidence=kw.get("evidence", {}),
        state_history=kw.get("history", []),
        elapsed_seconds=0.0,
        error=kw.get("error"),
    )


@pytest.fixture
def orch(tmp_path):
    registry = ArtifactRegistry(tmp_path / "audit_registry.jsonl")
    o = WorkflowOrchestrator(
        run_dir=tmp_path / "run",
        contracts_dir=tmp_path / "contracts",
        config_path=tmp_path / "config.yaml",
        registry=registry,
    )
    return o


def _stub_engine(orchestrator, state):
    """Make the stubbed engine return `state`, and record what it was called with."""
    calls = []

    def run(contract_file, inputs, artifact_ids):
        calls.append({"contract_file": contract_file, "inputs": inputs,
                      "artifact_ids": artifact_ids})
        return _result(state)

    orchestrator.engine.run = run
    return calls


# ---------------------------------------------------------------------------
# The halt guard — the line that makes a HALT actually stop the workflow
# ---------------------------------------------------------------------------

def test_halt_sets_the_flag_and_blocks_the_next_stage(orch):
    _stub_engine(orch, AuditState.HALT)

    orch.run_stage("G0", "g0.yaml", tool_fn=None, tool_outputs={})

    assert orch._halted is True
    with pytest.raises(RuntimeError, match="halted before G1"):
        orch.run_stage("G1", "g1.yaml", tool_fn=None, tool_outputs={})


def test_a_halted_workflow_records_no_further_stages(orch):
    _stub_engine(orch, AuditState.HALT)
    orch.run_stage("G0", "g0.yaml", tool_fn=None, tool_outputs={})

    with pytest.raises(RuntimeError):
        orch.run_stage("G1", "g1.yaml", tool_fn=None, tool_outputs={})

    assert [s.stage for s in orch.stage_results] == ["G0"]


def test_pass_leaves_the_workflow_runnable(orch):
    _stub_engine(orch, AuditState.PASS)

    first = orch.run_stage("G0", "g0.yaml", tool_fn=None, tool_outputs={})
    second = orch.run_stage("G1", "g1.yaml", tool_fn=None, tool_outputs={})

    assert orch._halted is False
    assert orch._escalation_pending is False
    assert isinstance(first, StageResult) and second.stage == "G1"


# ---------------------------------------------------------------------------
# Outcome dispatch
# ---------------------------------------------------------------------------

def test_escalate_sets_escalation_pending_without_halting(orch):
    _stub_engine(orch, AuditState.ESCALATE)

    orch.run_stage("G0", "g0.yaml", tool_fn=None, tool_outputs={})

    assert orch._escalation_pending is True
    assert orch._halted is False


def test_skipped_by_protocol_sets_neither_flag(orch):
    _stub_engine(orch, AuditState.SKIPPED_BY_PROTOCOL)

    orch.run_stage("G7", "g7.yaml", tool_fn=None, tool_outputs={})

    assert orch._halted is False
    assert orch._escalation_pending is False


def test_auto_repair_branch_is_unreachable_from_the_engine(orch):
    """Reaching `_handle_outcome` with AUTO_REPAIR requires bypassing the engine.

    `ContractEngine._attempt_auto_repair` transitions to PASS, ESCALATE or HALT
    on every path, so a real `ContractResult` never carries AUTO_REPAIR out of
    the engine (see `test_auto_repair_always_resolves`). The branch is therefore
    dead in production, and the guard inside it — `if result.state !=
    AuditState.PASS` — is a tautology: the branch was only entered because the
    state equals AUTO_REPAIR, so the comparison cannot be false.

    The test drives it directly to document the behaviour, and the tautology's
    false arc stays uncovered because no input can reach it (ADR-AIIM-199).
    """
    orch.engine.run = lambda contract_file, inputs, artifact_ids: _result(
        AuditState.AUTO_REPAIR)

    orch.run_stage("G0", "g0.yaml", tool_fn=None, tool_outputs={})

    assert orch._escalation_pending is True
    assert orch._halted is False


def test_unhandled_state_is_silently_ignored(orch, caplog):
    """An unrecognised state falls through every branch and sets no flag.

    `_handle_outcome` is a chain of `if/elif` with no trailing `else`. A state
    outside the enum — or `None` — therefore produces no log line, no halt and
    no escalation: the workflow continues as though the gate had passed. That
    is a fail-*open* default in a component whose stated design is fail-closed.

    No production path currently produces such a state, so this is a latent
    gap rather than an active defect, and it is reported rather than patched:
    changing the module would change the system that produced this study's
    results (ADR-AIIM-199).
    """
    orch.engine.run = lambda contract_file, inputs, artifact_ids: _result(None)

    with caplog.at_level(logging.DEBUG, logger="kdss_audit.orchestrator"):
        orch.run_stage("G0", "g0.yaml", tool_fn=None, tool_outputs={})

    assert orch._halted is False
    assert orch._escalation_pending is False
    assert not [r for r in caplog.records if r.levelno >= logging.WARNING]


# ---------------------------------------------------------------------------
# Tool invocation
# ---------------------------------------------------------------------------

def test_tool_fn_is_called_when_no_outputs_are_supplied(orch, tmp_path):
    calls = _stub_engine(orch, AuditState.PASS)
    seen = {}

    def tool_fn(run_dir, config_path, registry):
        seen["run_dir"] = run_dir
        seen["config_path"] = config_path
        return {"note": "produced"}

    orch.run_stage("G0", "g0.yaml", tool_fn=tool_fn)

    assert seen["run_dir"] == tmp_path / "run"
    assert seen["config_path"] == tmp_path / "config.yaml"
    assert calls[0]["inputs"] == {"note": "produced"}


def test_supplied_outputs_bypass_the_tool(orch):
    calls = _stub_engine(orch, AuditState.PASS)

    def tool_fn(**kw):  # pragma: no cover - must not be called
        raise AssertionError("tool_fn must not be called when outputs are given")

    orch.run_stage("G0", "g0.yaml", tool_fn=tool_fn, tool_outputs={"k": "v"})

    assert calls[0]["inputs"] == {"k": "v"}


# ---------------------------------------------------------------------------
# Artifact registration — what gets into the audit trail, and what does not
# ---------------------------------------------------------------------------

def test_existing_files_are_registered_with_a_hash(orch, tmp_path):
    calls = _stub_engine(orch, AuditState.PASS)
    f = tmp_path / "out.csv"
    f.write_text("a,b\n1,2\n", encoding="utf-8")

    stage = orch.run_stage("G1", "g1.yaml", tool_fn=None,
                           tool_outputs={"table": str(f)})

    assert stage.artifacts_produced == ["G1::table"]
    assert calls[0]["artifact_ids"] == ["G1::table"]
    import json
    recs = [json.loads(l) for l in
            (tmp_path / "audit_registry.jsonl").read_text(
                encoding="utf-8").strip().splitlines()]
    art = [r for r in recs if r["record_type"] == "artifact"][0]
    assert art["artifact_id"] == "G1::table"
    assert len(art["sha256"]) == 64


def test_non_path_values_are_not_registered(orch):
    _stub_engine(orch, AuditState.PASS)

    stage = orch.run_stage("G1", "g1.yaml", tool_fn=None,
                           tool_outputs={"n_rows": 42, "ok": True,
                                         "scores": [1, 2, 3]})

    assert stage.artifacts_produced == []


def test_paths_that_do_not_exist_are_silently_skipped(orch, tmp_path):
    """A declared output that was never written leaves no registry entry.

    Worth pinning: the orchestrator does not fail when a stage claims an
    artifact it did not produce. Detecting that is left to the gate contracts.
    """
    _stub_engine(orch, AuditState.PASS)

    stage = orch.run_stage("G1", "g1.yaml", tool_fn=None,
                           tool_outputs={"missing": str(tmp_path / "nope.csv")})

    assert stage.artifacts_produced == []


def test_empty_outputs_register_nothing(orch):
    _stub_engine(orch, AuditState.PASS)

    stage = orch.run_stage("G1", "g1.yaml", tool_fn=None, tool_outputs={})

    assert stage.artifacts_produced == []


def test_path_objects_are_accepted_alongside_strings(orch, tmp_path):
    _stub_engine(orch, AuditState.PASS)
    f = tmp_path / "out.txt"
    f.write_text("x", encoding="utf-8")

    stage = orch.run_stage("G1", "g1.yaml", tool_fn=None,
                           tool_outputs={"as_path": Path(f)})

    assert stage.artifacts_produced == ["G1::as_path"]


# ---------------------------------------------------------------------------
# Access log and escalation surface
# ---------------------------------------------------------------------------

def test_record_read_writes_an_access_log_entry(orch, tmp_path):
    orch.record_read(reader="G5a", file_path="/data/labels.csv",
                     purpose="label isolation check")

    import json
    recs = [json.loads(l) for l in
            (tmp_path / "audit_registry.jsonl").read_text(
                encoding="utf-8").strip().splitlines()]
    access = [r for r in recs if r["record_type"] == "access_log"]
    assert len(access) == 1
    assert access[0]["reader"] == "G5a"
    assert access[0]["purpose"] == "label isolation check"


def test_request_human_adjudication_raises_rather_than_prompting(orch):
    with pytest.raises(EscalationRequired) as exc:
        orch.request_human_adjudication(
            gate="G4",
            information_visible="score distribution only",
            decision_needed="accept or reject the shift")

    assert exc.value.gate == "G4"
    assert "G4" in str(exc.value)
    assert exc.value.decision_needed == "accept or reject the shift"


def test_build_report_lock_returns_a_lock_bound_to_the_registry(orch, tmp_path):
    lock = orch.build_report_lock(tmp_path / "report_lock.json")

    assert isinstance(lock, ReportLock)

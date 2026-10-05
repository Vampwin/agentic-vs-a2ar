"""test_contract_engine.py — branch coverage of the contract engine (ADR-AIIM-199).

Written because the external audit (ADR-AIIM-189) recorded 0 of 12 branches
executed in `kdss_audit/contract_engine.py`, and a re-reading of that artifact
(ADR-AIIM-198) established that no test in the suite imported the module at all.
The module is not dead code: `scripts/run_snp_screen.py` reaches it through
`WorkflowOrchestrator`.

These tests exercise the state-selection rule and the whole auto-repair path,
which is where the engine decides between PASS, ESCALATE and HALT. They assert
the documented behaviour; they do not assert that the behaviour is correct by
construction, and where the code and its comment disagree the test records what
the code does.

Evaluators are injected into `sys.modules` rather than written to disk, because
the engine resolves them with `importlib.import_module`, which consults
`sys.modules` first. That keeps the test independent of the evaluator package.
"""

from __future__ import annotations

import sys
import types

import pytest
import yaml

from kdss_audit.artifact_registry import ArtifactRegistry
from kdss_audit.contract_engine import ContractEngine, ContractResult
from kdss_audit.state_machine import AuditState


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def engine(tmp_path):
    contracts = tmp_path / "contracts"
    contracts.mkdir()
    registry = ArtifactRegistry(tmp_path / "audit_registry.jsonl")
    return ContractEngine(contracts_dir=contracts, registry=registry), contracts


@pytest.fixture
def fake_evaluator(request):
    """Register a throwaway evaluator module and remove it afterwards."""
    created: list[str] = []

    def _make(name, evaluate=None, repair=None):
        mod = types.ModuleType(name)
        if evaluate is not None:
            mod.evaluate = evaluate
        if repair is not None:
            mod.repair = repair
        sys.modules[name] = mod
        created.append(name)
        return name

    yield _make

    for name in created:
        sys.modules.pop(name, None)


def _write_contract(contracts_dir, name, evaluator, **extra):
    body = {"contract_id": name.replace(".yaml", ""), "evaluator": evaluator}
    body.update(extra)
    path = contracts_dir / name
    path.write_text(yaml.safe_dump(body), encoding="utf-8")
    return name


def _ev(**payload):
    """An evaluate() that returns a fixed payload."""
    def evaluate(contract, inputs, registry):
        return payload
    return evaluate


# ---------------------------------------------------------------------------
# State selection: the rule the whole contract layer rests on
# ---------------------------------------------------------------------------

def test_evaluator_exception_halts_and_is_recorded(engine, fake_evaluator):
    """An evaluator that raises must fail closed to HALT, not propagate."""
    eng, contracts = engine

    def boom(contract, inputs, registry):
        raise ValueError("evaluator exploded")

    name = fake_evaluator("fake_eval_boom", evaluate=boom)
    cf = _write_contract(contracts, "g_boom.yaml", name)

    result = eng.run(contract_file=cf, inputs={}, artifact_ids=[])

    assert result.state is AuditState.HALT
    assert "evaluator exploded" in result.error
    assert result.checks_passed == [] and result.checks_failed == []
    # The failure is durable, not just returned.
    assert "evaluator_exception" in result.state_history[0]["evidence"]


def test_explicit_state_from_evaluator_is_honoured(engine, fake_evaluator):
    """An evaluator naming its own state overrides the default rule."""
    eng, contracts = engine
    name = fake_evaluator(
        "fake_eval_explicit",
        evaluate=_ev(state="ESCALATE", failed=["c1"], summary="explicit"))
    cf = _write_contract(contracts, "g_explicit.yaml", name)

    result = eng.run(contract_file=cf, inputs={}, artifact_ids=[])

    assert result.state is AuditState.ESCALATE


def test_failure_without_declared_criticality_halts(engine, fake_evaluator):
    """`critical_failures` defaults to `failed`, so an undeclared failure HALTs.

    This is the fail-closed direction: an evaluator that reports a failure but
    forgets to classify it is treated as critical rather than as a warning.
    """
    eng, contracts = engine
    name = fake_evaluator("fake_eval_undeclared",
                          evaluate=_ev(failed=["c1"], passed=["c0"]))
    cf = _write_contract(contracts, "g_undeclared.yaml", name)

    result = eng.run(contract_file=cf, inputs={}, artifact_ids=[])

    assert result.state is AuditState.HALT
    assert result.checks_failed == ["c1"]


def test_failure_declared_non_critical_escalates(engine, fake_evaluator):
    """An explicitly empty `critical_failures` downgrades HALT to ESCALATE."""
    eng, contracts = engine
    name = fake_evaluator("fake_eval_noncritical",
                          evaluate=_ev(failed=["c1"], critical_failures=[]))
    cf = _write_contract(contracts, "g_noncritical.yaml", name)

    result = eng.run(contract_file=cf, inputs={}, artifact_ids=[])

    assert result.state is AuditState.ESCALATE


def test_no_failures_passes(engine, fake_evaluator):
    eng, contracts = engine
    name = fake_evaluator("fake_eval_clean",
                          evaluate=_ev(passed=["c0", "c1"], warned=["w0"],
                                       evidence={"k": 1}))
    cf = _write_contract(contracts, "g_clean.yaml", name)

    result = eng.run(contract_file=cf, inputs={}, artifact_ids=[])

    assert result.state is AuditState.PASS
    assert result.checks_passed == ["c0", "c1"]
    assert result.checks_warned == ["w0"]


# ---------------------------------------------------------------------------
# Auto-repair: every exit from this path is a fail-closed decision
# ---------------------------------------------------------------------------

def test_auto_repair_without_reversible_action_escalates(engine, fake_evaluator):
    """AUTO_REPAIR with no prespecified reversible action must not self-clear."""
    eng, contracts = engine
    name = fake_evaluator("fake_eval_ar_noaction",
                          evaluate=_ev(state="AUTO_REPAIR", failed=["c1"]))
    cf = _write_contract(contracts, "g_ar_noaction.yaml", name,
                         on_fail=["notify_human"])

    result = eng.run(contract_file=cf, inputs={}, artifact_ids=[])

    assert result.state is AuditState.ESCALATE
    reasons = [h["reason"] for h in result.state_history]
    assert any("no reversible action available" in r for r in reasons)


def test_auto_repair_without_repair_function_escalates(engine, fake_evaluator):
    """The action is allowed but the evaluator cannot perform it."""
    eng, contracts = engine
    name = fake_evaluator("fake_eval_ar_nofn",
                          evaluate=_ev(state="AUTO_REPAIR", failed=["c1"]))
    cf = _write_contract(contracts, "g_ar_nofn.yaml", name,
                         on_fail=["rerun_with_same_frozen_parameters"])

    result = eng.run(contract_file=cf, inputs={}, artifact_ids=[])

    assert result.state is AuditState.ESCALATE
    reasons = [h["reason"] for h in result.state_history]
    assert any("no repair() function" in r for r in reasons)


def test_auto_repair_success_passes_and_reports_repaired_checks(engine,
                                                                fake_evaluator):
    eng, contracts = engine

    def repair(contract, inputs):
        return {"success": True, "passed": ["c1"], "failed": [],
                "evidence_keys": ["rerun_ok"]}

    name = fake_evaluator("fake_eval_ar_ok",
                          evaluate=_ev(state="AUTO_REPAIR", failed=["c1"]),
                          repair=repair)
    cf = _write_contract(contracts, "g_ar_ok.yaml", name,
                         on_fail=["rerun_with_same_frozen_parameters"])

    result = eng.run(contract_file=cf, inputs={}, artifact_ids=[])

    assert result.state is AuditState.PASS
    # The repaired outcome, not the pre-repair one, is what gets recorded.
    assert result.checks_passed == ["c1"]
    assert result.checks_failed == []


def test_auto_repair_failure_escalates(engine, fake_evaluator):
    """A repair that reports failure must escalate, never silently pass."""
    eng, contracts = engine

    def repair(contract, inputs):
        return {"success": False, "reason": "rerun diverged",
                "evidence_keys": ["rerun_diverged"]}

    name = fake_evaluator("fake_eval_ar_fail",
                          evaluate=_ev(state="AUTO_REPAIR", failed=["c1"]),
                          repair=repair)
    cf = _write_contract(contracts, "g_ar_fail.yaml", name,
                         on_fail=["rerun_with_same_frozen_parameters"])

    result = eng.run(contract_file=cf, inputs={}, artifact_ids=[])

    assert result.state is AuditState.ESCALATE
    reasons = [h["reason"] for h in result.state_history]
    assert any("rerun diverged" in r for r in reasons)


# ---------------------------------------------------------------------------
# The engine never leaves a gate sitting in AUTO_REPAIR
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("on_fail,repair,expected", [
    (["notify_human"], None, AuditState.ESCALATE),
    (["rerun_with_same_frozen_parameters"], None, AuditState.ESCALATE),
    (["rerun_with_same_frozen_parameters"],
     lambda contract, inputs: {"success": True}, AuditState.PASS),
    (["rerun_with_same_frozen_parameters"],
     lambda contract, inputs: {"success": False}, AuditState.ESCALATE),
])
def test_auto_repair_always_resolves(engine, fake_evaluator, on_fail, repair,
                                     expected):
    """Every auto-repair exit resolves to a terminal-for-the-gate state.

    This is what makes `WorkflowOrchestrator._handle_outcome`'s AUTO_REPAIR
    branch unreachable in production (ADR-AIIM-199): the engine cannot return a
    result still in AUTO_REPAIR.
    """
    eng, contracts = engine
    name = fake_evaluator(f"fake_eval_resolve_{id(repair)}",
                          evaluate=_ev(state="AUTO_REPAIR", failed=["c1"]),
                          repair=repair)
    cf = _write_contract(contracts, f"g_resolve_{id(repair)}.yaml", name,
                         on_fail=on_fail)

    result = eng.run(contract_file=cf, inputs={}, artifact_ids=[])

    assert result.state is expected
    assert result.state is not AuditState.AUTO_REPAIR


def test_outcome_is_written_to_the_registry(engine, fake_evaluator, tmp_path):
    """The engine records outside its own return value, or nothing is auditable."""
    eng, contracts = engine
    name = fake_evaluator("fake_eval_recorded", evaluate=_ev(passed=["c0"]))
    cf = _write_contract(contracts, "g_recorded.yaml", name)

    eng.run(contract_file=cf, inputs={}, artifact_ids=["a::1"])

    lines = (tmp_path / "audit_registry.jsonl").read_text(
        encoding="utf-8").strip().splitlines()
    assert len(lines) == 1
    import json
    rec = json.loads(lines[0])
    assert rec["record_type"] == "contract_outcome"
    assert rec["outcome"] == "PASS"
    assert rec["artifact_ids"] == ["a::1"]


def test_load_contract_reads_from_the_configured_directory(engine):
    eng, contracts = engine
    _write_contract(contracts, "g_load.yaml", "evaluators.G0_evaluator",
                    threshold=0.5)

    loaded = eng.load_contract("g_load.yaml")

    assert loaded["contract_id"] == "g_load"
    assert loaded["threshold"] == 0.5


def test_result_is_a_contract_result_with_timing(engine, fake_evaluator):
    eng, contracts = engine
    name = fake_evaluator("fake_eval_timing", evaluate=_ev(passed=[]))
    cf = _write_contract(contracts, "g_timing.yaml", name)

    result = eng.run(contract_file=cf, inputs={}, artifact_ids=[])

    assert isinstance(result, ContractResult)
    assert result.elapsed_seconds >= 0.0
    assert result.error is None

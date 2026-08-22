"""
Tests for retrospective_enrichment and its reconciliation with
control_representativeness (ADR-AIIM-073).

control_representativeness is a ligand-similarity proxy that cannot in
principle be satisfied by a structurally novel library. retrospective_
enrichment is the direct, structure-based alternative: dock a labeled
benchmark under the identical protocol and measure whether it actually
separates actives from inactives. These tests pin three things: the check
fires correctly on its own evidence, a passing result supersedes a failed
control_representativeness rather than leaving an unresolved ESCALATE, and
contracts that do not declare retrospective_enrichment behave exactly as
before this amendment.
"""

from __future__ import annotations

import pytest

from evaluators.G3_evaluator import evaluate

ZM241385 = "Nc1nc(NCCc2ccc(O)cc2)nc2nc(-c3ccco3)nn12"
ORCHID_NPS = [
    "COc1cc(O)cc2c1-c1c(O)cc(O)cc1CC2",
    "COc1c(O)cc2ccc3cc(O)ccc3c2c1OC",
    "O=C1C[C@@H](c2ccc(O)cc2)Oc2cc(O)cc(O)c21",
    "O=c1cc(-c2ccc(O)cc2)oc2cc(O)cc(O)c12",
]

pytest.importorskip("rdkit")

_BOTH_DECLARED = {
    "checks": [
        {"id": "control_representativeness", "similarity_cutoff": 0.35,
         "min_similarity_coverage": 0.50, "min_scaffold_coverage": 0.10,
         "severity": "high"},
        {"id": "retrospective_enrichment",
         "threshold_auc_min": 0.60, "threshold_bedroc20_min": 0.40,
         "severity": "critical"},
    ]
}

_ONLY_CONTROL_REP = {
    "checks": [
        {"id": "control_representativeness", "similarity_cutoff": 0.35,
         "min_similarity_coverage": 0.50, "min_scaffold_coverage": 0.10,
         "severity": "high"},
    ]
}


def _run(contract, extra_inputs=None):
    inputs = {"control_ligand_smiles": [ZM241385], "screening_library_smiles": ORCHID_NPS}
    inputs.update(extra_inputs or {})
    return evaluate(contract=contract, inputs=inputs, registry=None)


def test_passing_enrichment_alone():
    result = evaluate(
        contract={"checks": [{"id": "retrospective_enrichment",
                              "threshold_auc_min": 0.60,
                              "threshold_bedroc20_min": 0.40,
                              "severity": "critical"}]},
        inputs={"retrospective_enrichment": {"auc_roc": 0.72, "bedroc20": 0.45,
                                             "protocol_matches_screen": True}},
        registry=None,
    )
    assert "retrospective_enrichment" in result["passed"]
    assert result["evidence"]["retrospective_enrichment"]["auc_pass"] is True


def test_failing_enrichment_alone_is_critical():
    result = evaluate(
        contract={"checks": [{"id": "retrospective_enrichment",
                              "threshold_auc_min": 0.60,
                              "threshold_bedroc20_min": 0.40,
                              "severity": "critical"}]},
        inputs={"retrospective_enrichment": {"auc_roc": 0.51,
                                             "protocol_matches_screen": True}},
        registry=None,
    )
    assert "retrospective_enrichment" in result["failed"]
    assert "retrospective_enrichment" in result["critical_failures"]
    assert result["state"] == "HALT"


def test_protocol_mismatch_is_not_valid_evidence_even_if_auc_looks_good():
    """
    Benchmark evidence from a different protocol (different water shell,
    different box) says nothing about the protocol actually being used.
    """
    result = evaluate(
        contract={"checks": [{"id": "retrospective_enrichment",
                              "threshold_auc_min": 0.60,
                              "severity": "critical"}]},
        inputs={"retrospective_enrichment": {"auc_roc": 0.95,
                                             "protocol_matches_screen": False}},
        registry=None,
    )
    assert "retrospective_enrichment" in result["failed"]
    assert "retrospective_enrichment" in result["critical_failures"]
    assert "not valid evidence" in result["evidence"]["retrospective_enrichment"]["error"]


def test_missing_enrichment_warns_when_not_required_by_a_failed_proxy():
    """Simply absent, with nothing else forcing the issue: a warning, not a failure."""
    result = evaluate(
        contract={"checks": [{"id": "retrospective_enrichment", "severity": "critical"}]},
        inputs={},
        registry=None,
    )
    assert "retrospective_enrichment" in result["warned"]


def test_passing_enrichment_supersedes_a_failed_control_representativeness():
    """
    The real case: ZM241385 does not resemble the orchid library
    (control_representativeness fails), but direct structure-based evidence
    that the protocol separates actives from inactives closes that gap.
    """
    result = _run(_BOTH_DECLARED,
                  {"retrospective_enrichment": {"auc_roc": 0.65, "bedroc20": 0.42,
                                                "protocol_matches_screen": True}})
    assert "control_representativeness" in result["failed"]      # still true, reported
    assert "retrospective_enrichment" in result["passed"]
    assert result["evidence"]["retrospective_enrichment"]["supersedes_control_representativeness"]
    assert result["evidence"]["control_representativeness"]["superseded_by"] == "retrospective_enrichment"
    assert "retrospective_enrichment" not in result["critical_failures"]


def test_missing_enrichment_becomes_critical_once_control_representativeness_fails():
    """
    Both declared, but no retrospective_enrichment evidence was supplied: the
    similarity proxy failed and nothing structure-based corroborates
    transferability. That combination must not be allowed to sit at ESCALATE
    indefinitely -- it is a critical gap in the evidence.
    """
    result = _run(_BOTH_DECLARED)
    assert "control_representativeness" in result["failed"]
    assert "retrospective_enrichment" in result["critical_failures"]
    assert result["state"] == "HALT"


def test_contracts_without_retrospective_enrichment_are_unaffected():
    """
    Backward compatibility: a contract that never declares
    retrospective_enrichment must behave exactly as it did before this
    amendment -- ESCALATE, not HALT, on a failed similarity proxy alone.
    """
    result = _run(_ONLY_CONTROL_REP)
    assert "control_representativeness" in result["failed"]
    assert "retrospective_enrichment" not in result["evidence"]
    assert result["state"] == "ESCALATE"
    assert result["critical_failures"] == []

"""
Tests for the validation-representativeness check.

The check exists because a redocking control can pass while saying nothing about
the compounds a pipeline goes on to rank. These tests pin both directions: it
must fire on a control that is chemically unrelated to the library, and it must
stay quiet when the control sits inside the library's chemical space. A check
that only ever fires is as useless as one that never does.
"""

from __future__ import annotations

import pytest

from evaluators._control_representativeness import assess
from evaluators.G3_evaluator import evaluate

# ZM241385 — the A2AR co-crystal ligand used as the redocking control.
ZM241385 = "Nc1nc(NCCc2ccc(O)cc2)nc2nc(-c3ccco3)nn12"

# Representative orchid natural products: phenanthrenes, bibenzyls, flavonoids.
ORCHID_NPS = [
    "COc1cc(O)cc2c1-c1c(O)cc(O)cc1CC2",          # dihydrophenanthrene
    "COc1c(O)cc2ccc3cc(O)ccc3c2c1OC",            # nudol-like phenanthrene
    "O=C1C[C@@H](c2ccc(O)cc2)Oc2cc(O)cc(O)c21",  # naringenin
    "O=c1cc(-c2ccc(O)cc2)oc2cc(O)cc(O)c12",      # apigenin
    "COc1cc2c(c(OC)c1O)-c1ccc(O)cc1CC2",         # flavanthridin-like
    "COC1=CC(=O)c2c(ccc3cc(O)ccc23)C1=O",        # densiflorol-like quinone
    "Oc1ccc(CCc2cc(O)cc(O)c2)cc1",               # bibenzyl
    "COc1c(OC)c2c3c(c1O)c(=O)oc1cc(O)cc(c13)CC2",  # coeloginin-like
]

# Close analogues of the control — a library the control genuinely speaks for.
TRIAZOLOPYRIMIDINES = [
    "Nc1nc(NCCc2ccc(O)cc2)nc2nc(-c3ccco3)nn12",   # the control itself
    "Nc1nc(-c2ccco2)c2nnn(Cc3ccc(N)c(C)c3)c2n1",  # 9XT
    "Nc1nc(NCCc2ccccc2)nc2nc(-c3ccco3)nn12",
    "Nc1nc(NCCc2ccc(F)cc2)nc2nc(-c3ccco3)nn12",
    "Nc1nc(NCCc2ccc(O)cc2)nc2nc(-c3cccs3)nn12",
]

_CONTRACT = {
    "checks": [
        {"id": "control_representativeness",
         "similarity_cutoff": 0.35,
         "min_similarity_coverage": 0.50,
         "min_scaffold_coverage": 0.10,
         "severity": "high"},
    ]
}

pytest.importorskip("rdkit")


def _run(controls, library):
    return evaluate(
        contract=_CONTRACT,
        inputs={"control_ligand_smiles": controls,
                "screening_library_smiles": library},
        registry=None,
    )


def test_unrelated_control_is_flagged():
    """
    ZM241385 against orchid natural products: the real case. A synthetic
    aminotriazolopyrimidine shares no framework with phenanthrenes, bibenzyls or
    flavonoids, so a passing redocking result cannot license ranking them.
    """
    result = _run([ZM241385], ORCHID_NPS)
    assert "control_representativeness" in result["failed"]
    ev = result["evidence"]["control_representativeness"]
    assert ev["similarity_coverage"] < 0.5
    assert ev["scaffold_coverage"] < 0.1
    assert "does not demonstrate the protocol transfers" in ev["interpretation"]


def test_related_control_passes():
    """The check must not fire when the control sits inside the library's space."""
    result = _run([ZM241385], TRIAZOLOPYRIMIDINES)
    assert "control_representativeness" in result["passed"]
    ev = result["evidence"]["control_representativeness"]
    assert ev["similarity_coverage"] >= 0.5 or ev["scaffold_coverage"] >= 0.1


def test_scaffold_route_alone_can_pass():
    """
    Either route suffices. A library sharing the control's framework passes on
    scaffold coverage even where fingerprint similarity is diluted by
    substituents.
    """
    report = assess([ZM241385], TRIAZOLOPYRIMIDINES)
    assert report["available"]
    assert report["scaffold_coverage"] > 0.0


def test_missing_inputs_warn_rather_than_pass():
    """
    Representativeness that could not be measured must never be reported as
    established — the same rule applied to unverifiable hash anchors in G0.
    """
    result = _run([], ORCHID_NPS)
    assert "control_representativeness" in result["warned"]
    assert "control_representativeness" not in result["passed"]


def test_report_exposes_the_numbers():
    """
    The thresholds are conventions, so the underlying measurements must always
    be on the record for a human to judge, not just a pass/fail verdict.
    """
    report = assess([ZM241385], ORCHID_NPS)
    for key in ("similarity_coverage", "scaffold_coverage",
                "nearest_similarity_max", "nearest_similarity_median",
                "n_controls", "n_library", "descriptors"):
        assert key in report, key
    desc = report["descriptors"][ZM241385]
    assert desc["n_descriptors"] == 8
    assert set(desc["inside_library_band"]) | set(desc["outside_library_band"])


def test_state_is_escalate_not_halt():
    """
    An unrepresentative control does not invalidate the pose result; it means a
    human must decide whether ranking is justified. High severity, so the stage
    escalates rather than halting.
    """
    result = _run([ZM241385], ORCHID_NPS)
    assert result["state"] == "ESCALATE"
    assert "control_representativeness" not in result["critical_failures"]

"""
Tests for the domain-gate bias probe (ADR-AIIM-073, E2).

E2 found SNP-Lib 0% in-domain across six independent structural descriptor
spaces while 64.4% in-domain physicochemically. These tests pin the module
that turns that comparison into a checkable signal: it must fire when
structural coverage is near zero and physicochemical coverage is
substantial (the "structurally novel, not chemical garbage" signature), and
it must stay quiet when there is nothing to distinguish -- a library that is
out of domain everywhere is a different, equally real finding this probe
must not paper over.
"""

from __future__ import annotations

from evaluators import _domain_gate_bias_probe as probe
from evaluators.G6_evaluator import evaluate

# The actual E2 result (2026-08-15), reproduced here as fixed known-answer data.
SNPLIB_E2_RESULT = {
    "ECFP4": {"coverage_fraction": 0.0},
    "FCFP4": {"coverage_fraction": 0.0},
    "MACCS": {"coverage_fraction": 0.025},
    "AtomPair": {"coverage_fraction": 0.008},
    "TopTorsion": {"coverage_fraction": 0.0},
    "Gobbi_Pharm2D": {"coverage_fraction": 0.0},
    "ErG": {"coverage_fraction": 0.004},
    "PhysChem": {"coverage_fraction": 0.644},
}

_CONTRACT = {"checks": [{"id": "applicability_domain_multi_space",
                         "min_spaces": 3, "severity": "high"}]}


def _run(coverage_by_space):
    return evaluate(_CONTRACT,
                    {"applicability_domain_by_space": coverage_by_space},
                    registry=None)


def test_divergence_detected_on_the_real_e2_result():
    report = probe.assess(SNPLIB_E2_RESULT)
    assert report["available"]
    assert report["n_spaces"] == 8
    assert report["structural_physchem_divergence"] is True
    assert report["structural_max_coverage"] <= 0.05
    assert report["physicochemical_max_coverage"] >= 0.30
    assert "structurally novel" in report["interpretation"]


def test_uniformly_out_of_domain_library_is_not_flagged_as_divergent():
    """
    A library that is out of domain everywhere -- structural AND
    physicochemical -- is a real, different finding. The probe must not
    claim divergence it did not find.
    """
    uniform = {"ECFP4": {"coverage_fraction": 0.0},
              "MACCS": {"coverage_fraction": 0.01},
              "PhysChem": {"coverage_fraction": 0.02}}
    report = probe.assess(uniform)
    assert report["structural_physchem_divergence"] is False


def test_ordinary_in_domain_library_is_not_flagged():
    """A library that resembles training data in every space: nothing to see."""
    ordinary = {"ECFP4": {"coverage_fraction": 0.85},
               "MACCS": {"coverage_fraction": 0.90},
               "PhysChem": {"coverage_fraction": 0.95}}
    report = probe.assess(ordinary)
    assert report["structural_physchem_divergence"] is False


def test_too_few_spaces_warns_rather_than_asserts():
    result = _run({"ECFP4": {"coverage_fraction": 0.0}})
    assert "applicability_domain_multi_space" in result["warned"]
    assert result["evidence"]["applicability_domain_multi_space"]["enough_spaces"] is False


def test_missing_data_warns_rather_than_passes():
    result = _run({})
    assert "applicability_domain_multi_space" in result["warned"]
    assert result["evidence"]["applicability_domain_multi_space"]["available"] is False


def test_enough_spaces_passes_and_carries_the_divergence_finding_forward():
    """
    Reaching PASS must never mean the divergence finding is discarded --
    only that there was enough evidence to make a judgement at all. The
    finding itself travels in evidence for the human decision downstream.
    """
    result = _run(SNPLIB_E2_RESULT)
    assert "applicability_domain_multi_space" in result["passed"]
    ev = result["evidence"]["applicability_domain_multi_space"]
    assert ev["structural_physchem_divergence"] is True
    assert result["state"] in ("PASS", "ESCALATE")  # never silently HALT on this alone
    assert "applicability_domain_multi_space" not in result.get("critical_failures", [])

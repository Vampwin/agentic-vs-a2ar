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


# --- false-refusal probe (ADR-AIIM-147) ------------------------------------

def test_false_refusal_requires_a_separation_statement():
    """A rate measured on the compounds that define the domain is zero by
    construction. The evaluator must refuse rather than report it."""
    r = probe.assess_false_refusal(
        {"ECFP4": [True, False, True]}, "", "some panel")
    assert r["available"] is False
    assert "define the domain" in r["reason"]


def test_false_refusal_rejects_mismatched_panel_sizes():
    r = probe.assess_false_refusal(
        {"ECFP4": [True, False], "MACCS": [True]},
        "held out by scaffold", "panel")
    assert r["available"] is False
    assert "disagree on panel size" in r["reason"]


def test_false_refusal_rejects_an_empty_panel():
    r = probe.assess_false_refusal(
        {"ECFP4": []}, "held out by scaffold", "panel")
    assert r["available"] is False


def test_false_refusal_counts_per_space_and_never_averages():
    r = probe.assess_false_refusal(
        {"ECFP4": [True, True, False, False],      # 50 %
         "PhysChem": [True, True, True, False]},   # 25 %
        "held out by scaffold", "four known actives")
    assert r["available"] is True
    assert r["by_space"]["ECFP4"]["false_refusal_rate"] == 0.5
    assert r["by_space"]["PhysChem"]["false_refusal_rate"] == 0.25
    assert r["worst_space"] == "ECFP4"
    assert r["best_space"] == "PhysChem"
    assert r["spread"] == 0.25
    # The spaces are not commensurable; no pooled rate is offered.
    assert not any("overall" in k or "mean" in k for k in r)


def test_false_refusal_flags_only_spaces_at_or_above_the_threshold():
    r = probe.assess_false_refusal(
        {"ECFP4": [False, False, True, True],   # 50 %, at the threshold
         "MACCS": [False, True, True, True]},   # 25 %
        "held out by scaffold", "panel")
    assert r["spaces_above_threshold"] == ["ECFP4"]


def test_false_refusal_separates_structural_from_physicochemical():
    r = probe.assess_false_refusal(
        {"ECFP4": [False, False, False, True],
         "PhysChem": [True, True, True, True]},
        "held out by scaffold", "panel")
    sp = r["structural_vs_physicochemical"]
    assert sp["structural_worst"] == 0.75
    assert sp["physicochemical_worst"] == 0.0
    assert sp["difference"] == 0.75


# --- the guard against a contract check with no evaluator branch -----------
# ADR-AIIM-159. Adding applicability_domain_false_refusal to G6_selection.yaml
# v1.3.0 before the evaluator had a branch for it made the check vanish: the
# if/elif chain fell through and the id appeared in neither passed, failed nor
# warned. The gate would have gone on reporting PASS while silently not
# evaluating a check its own contract declares.

def test_a_contract_check_with_no_evaluator_branch_fails_critically():
    from evaluators import G6_evaluator
    contract = {
        "contract_id": "G6_selection_v1",
        "checks": [{"id": "a_check_no_evaluator_implements", "severity": "high"}],
    }

    class _Registry:
        def record_contract_outcome(self, *a, **k):
            pass

    r = G6_evaluator.evaluate(contract, {}, _Registry())
    assert "a_check_no_evaluator_implements" in r["failed"]
    assert r["state"] == "HALT", "an unimplemented check must not let the gate pass"
    ev = r["evidence"]["a_check_no_evaluator_implements"]
    assert ev["available"] is False
    assert "not implemented" in ev["reason"]

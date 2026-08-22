"""
Integration test: the confidence_tier_assessment check as wired into the real
G8_reporting.yaml contract and G8_evaluator.py, not just the isolated
_confidence_tier module (see test_confidence_tier.py for the unit-level
tests, including the retroactive SNP-Lib check).
"""

from __future__ import annotations

import yaml

from evaluators import G8_evaluator

_CONTRACT = yaml.safe_load(open("contracts/G8_reporting.yaml", encoding="utf-8"))


def _minimal_report_inputs(confidence_tier_inputs: dict | None = None) -> dict:
    return {
        "manuscript_result_snippets": [],
        "figure_manifest_json": [],
        "table_manifest_json": [],
        "artifact_registry_snapshot": [],
        "manuscript_metadata": {
            "human_intervention_count": 0,
            "ai_use_declaration": "Claude Sonnet 5 used under human direction; see Methods.",
        },
        "confidence_tier_inputs": confidence_tier_inputs,
    }


def test_missing_confidence_tier_inputs_warns():
    result = G8_evaluator.evaluate(_CONTRACT, _minimal_report_inputs(None), registry=None)
    assert "confidence_tier_assessment" in result["warned"]
    assert "confidence_tier_assessment" not in result["passed"]


def test_confidence_tier_check_reproduces_snplib_very_low_end_to_end():
    ct_inputs = {
        "paths_attempted": ["SBDD", "LBDD"],
        "sbdd": {
            "rmsd_angstrom": 0.379,
            "control_representativeness_report": {
                "available": True, "similarity_coverage": 0.000, "scaffold_coverage": 0.000,
            },
            "ligand_dependent_assumption": True,
            "structural_indirectness_detail": "water network shown ligand-specific by cross-docking",
        },
        "lbdd": {
            "ad_coverage_fraction": 0.0,
            "cv_auc": 0.907,
            "n_positive": None,
        },
    }
    result = G8_evaluator.evaluate(_CONTRACT, _minimal_report_inputs(ct_inputs), registry=None)
    assert "confidence_tier_assessment" in result["passed"]
    ev = result["evidence"]["confidence_tier_assessment"]
    assert ev["overall_tier"] == "VERY LOW"
    assert ev["sbdd"]["final_tier"] == "VERY LOW"
    assert ev["lbdd"]["final_tier"] == "VERY LOW"
    # High-severity ESCALATE, not a critical HALT -- a VERY LOW tier is a
    # valid, reportable outcome, not a contract violation.
    assert "confidence_tier_assessment" not in result["critical_failures"]


def test_threshold_immutability_warns_with_no_registered_anchor():
    inputs = _minimal_report_inputs(None)
    result = G8_evaluator.evaluate(_CONTRACT, inputs, registry=None)
    assert "confidence_tier_threshold_immutability" in result["warned"]


def test_threshold_immutability_passes_when_anchor_matches():
    check = next(c for c in _CONTRACT["checks"] if c["id"] == "confidence_tier_assessment")
    current = {k: v for k, v in check.items() if k not in ("id", "description", "severity")}
    inputs = _minimal_report_inputs(None)
    inputs["registered_confidence_tier_thresholds"] = current
    result = G8_evaluator.evaluate(_CONTRACT, inputs, registry=None)
    assert "confidence_tier_threshold_immutability" in result["passed"]


def test_threshold_immutability_halts_on_mismatch():
    inputs = _minimal_report_inputs(None)
    inputs["registered_confidence_tier_thresholds"] = {"sbdd_pose_pass_threshold_angstrom": 3.0}
    result = G8_evaluator.evaluate(_CONTRACT, inputs, registry=None)
    assert "confidence_tier_threshold_immutability" in result["failed"]
    assert "confidence_tier_threshold_immutability" in result["critical_failures"]
    assert result["state"] == "HALT"


def test_confidence_tier_check_reports_high_when_evidence_is_strong():
    ct_inputs = {
        "paths_attempted": ["SBDD"],
        "sbdd": {
            "rmsd_angstrom": 0.3,
            "control_representativeness_report": {
                "available": True, "similarity_coverage": 0.75, "scaffold_coverage": 0.4,
            },
            "auc": 0.8, "ef1pct": 6.0, "bedroc20": 0.55,
        },
    }
    result = G8_evaluator.evaluate(_CONTRACT, _minimal_report_inputs(ct_inputs), registry=None)
    assert "confidence_tier_assessment" in result["passed"]
    assert result["evidence"]["confidence_tier_assessment"]["overall_tier"] == "HIGH"

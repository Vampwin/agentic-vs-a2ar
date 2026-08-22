"""
Tests for the GRADE-style confidence-tier framework
(00_project_control/17_CONFIDENCE_TIER_FRAMEWORK.md).

The most important test here is the retroactive check: applying the
framework to the real, already-known SNP-Lib/A2AR numbers must reproduce
ADR-AIIM-037's actual conclusion (no ranking should have been released).
A framework that disagreed with a conclusion already reached through careful
manual investigation would be a framework not to trust on a library where
the answer isn't already known.
"""

from __future__ import annotations

import yaml

from evaluators import _confidence_tier as ct
from evaluators._confidence_tier import (
    assess_lbdd,
    assess_sbdd,
    combine,
    rate_applicability_domain,
    rate_control_representativeness,
    rate_external_scoring_validity,
    rate_imprecision,
    rate_pose_validation_margin,
)


def test_confidence_tier_thresholds_frozen():
    """
    The contract (contracts/G8_reporting.yaml) is the authoritative, frozen
    copy of every cut point (ADR-AIIM-046); evaluators/_confidence_tier.py's
    DEFAULT_* constants are the fallback used when a caller doesn't pass the
    contract's values through. If these ever drift apart, a caller that
    forgets to wire the contract values would silently get different
    behaviour depending on the code path -- fail closed here instead.
    """
    contract = yaml.safe_load(open("contracts/G8_reporting.yaml", encoding="utf-8"))
    check = next(c for c in contract["checks"] if c["id"] == "confidence_tier_assessment")

    expected = {
        "sbdd_pose_pass_threshold_angstrom": ct.DEFAULT_POSE_PASS_THRESHOLD,
        "sbdd_pose_margin_threshold_angstrom": ct.DEFAULT_POSE_MARGIN_THRESHOLD,
        "sbdd_min_similarity_coverage": ct.DEFAULT_MIN_SIMILARITY_COVERAGE,
        "sbdd_min_scaffold_coverage": ct.DEFAULT_MIN_SCAFFOLD_COVERAGE,
        "sbdd_thin_margin_multiplier": ct.DEFAULT_THIN_MARGIN_MULTIPLIER,
        "sbdd_external_auc_min": ct.DEFAULT_EXTERNAL_AUC_MIN,
        "sbdd_external_ef1pct_min": ct.DEFAULT_EXTERNAL_EF1PCT_MIN,
        "sbdd_external_bedroc20_min": ct.DEFAULT_EXTERNAL_BEDROC20_MIN,
        "lbdd_ad_coverage_high": ct.DEFAULT_AD_COVERAGE_HIGH,
        "lbdd_ad_coverage_low": ct.DEFAULT_AD_COVERAGE_LOW,
        "lbdd_cv_auc_high": ct.DEFAULT_CV_AUC_HIGH,
        "lbdd_cv_auc_low": ct.DEFAULT_CV_AUC_LOW,
        "lbdd_external_auc_min": ct.DEFAULT_LBDD_EXTERNAL_AUC_MIN,
        "lbdd_imprecision_min_recommended": ct.DEFAULT_IMPRECISION_MIN_RECOMMENDED,
        "lbdd_imprecision_min_floor": ct.DEFAULT_IMPRECISION_MIN_FLOOR,
    }
    for key, module_value in expected.items():
        assert check[key] == module_value, (
            f"{key}: contract has {check[key]!r}, module default has {module_value!r}"
        )


# ── Individual domain raters ────────────────────────────────────────────────

def test_pose_validation_margin_not_serious_with_wide_margin():
    r = rate_pose_validation_margin(0.379)
    assert r["rating"] == "not_serious"


def test_pose_validation_margin_serious_when_thin():
    r = rate_pose_validation_margin(1.8)
    assert r["rating"] == "serious"


def test_pose_validation_margin_very_serious_when_failing():
    r = rate_pose_validation_margin(4.32)
    assert r["rating"] == "very_serious"


def test_pose_validation_margin_missing_is_very_serious():
    assert rate_pose_validation_margin(None)["rating"] == "very_serious"


def test_control_representativeness_very_serious_on_real_snplib_numbers():
    # The real numbers from the production run against the 236-compound library.
    report = {"available": True, "similarity_coverage": 0.000, "scaffold_coverage": 0.000}
    r = rate_control_representativeness(report)
    assert r["rating"] == "very_serious"


def test_control_representativeness_not_serious_when_well_covered():
    report = {"available": True, "similarity_coverage": 0.72, "scaffold_coverage": 0.35}
    r = rate_control_representativeness(report)
    assert r["rating"] == "not_serious"


def test_applicability_domain_very_serious_at_zero_coverage():
    # The real 172/172-outside-domain finding.
    r = rate_applicability_domain(0.0)
    assert r["rating"] == "very_serious"


def test_applicability_domain_not_serious_at_high_coverage():
    assert rate_applicability_domain(0.9)["rating"] == "not_serious"


def test_imprecision_very_serious_for_tiny_positive_class():
    # The actual n_positive=9 the development monitor was fit on.
    r = rate_imprecision(9)
    assert r["rating"] == "very_serious"


def test_imprecision_not_serious_for_adequate_sample():
    assert rate_imprecision(50)["rating"] == "not_serious"


def test_external_scoring_skipped_reports_not_serious():
    r = rate_external_scoring_validity(None, None, None, skip=True, skip_reason="test")
    assert r["rating"] == "not_serious"
    assert "not evaluated" in r["detail"]


# ── Weakest-link aggregation ────────────────────────────────────────────────

def test_single_very_serious_domain_floors_at_low_not_very_low():
    """One very-serious domain downgrades two tiers: HIGH -> LOW, not VERY LOW."""
    result = assess_sbdd(
        rmsd_angstrom=0.379,
        control_representativeness_report={"available": True, "similarity_coverage": 0.0,
                                             "scaffold_coverage": 0.0},
        ligand_dependent_assumption=False,
    )
    assert result["final_tier"] == "LOW"


def test_good_result_elsewhere_cannot_buy_back_a_very_serious_domain():
    """A strong external AUC must not rescue a very-serious control-representativeness finding."""
    result = assess_sbdd(
        rmsd_angstrom=0.2,
        control_representativeness_report={"available": True, "similarity_coverage": 0.0,
                                             "scaffold_coverage": 0.0},
        ligand_dependent_assumption=False,
        auc=0.95, ef1pct=10.0, bedroc20=0.9,
    )
    # External scoring is skipped (not scored as "not_serious" via good numbers) because
    # control representativeness already failed -- the AUC is not informative here.
    external = next(d for d in result["domains"] if d["domain"] == "external_scoring_validity")
    assert "not evaluated" in external["detail"]
    assert result["final_tier"] != "HIGH"


def test_all_domains_clean_gives_high_tier():
    result = assess_sbdd(
        rmsd_angstrom=0.3,
        control_representativeness_report={"available": True, "similarity_coverage": 0.7,
                                             "scaffold_coverage": 0.4},
        ligand_dependent_assumption=False,
        auc=0.75, ef1pct=8.0, bedroc20=0.6,
    )
    assert result["final_tier"] == "HIGH"


# ── Retroactive check against the real SNP-Lib / A2AR result ──────────────

def test_snplib_a2ar_retroactive_check_reproduces_very_low():
    """
    Reproduces ADR-AIIM-037's actual conclusion (no ranking released) from the
    real numbers already on record: control_representativeness = 0.000/0.000
    against the real 236-compound library, and the water network shown
    ligand-dependent by the cross-docking control (ADR-AIIM-036).
    """
    sbdd = assess_sbdd(
        rmsd_angstrom=0.379,  # ZM241385 -> 4EIY, own receptor, real value
        control_representativeness_report={
            "available": True, "similarity_coverage": 0.000, "scaffold_coverage": 0.000,
        },
        ligand_dependent_assumption=True,
        structural_indirectness_detail=(
            "pose validity depends on the full 15-water bridging network, shown "
            "ligand-specific by the cross-docking control (9.294 A in 5OLH)"
        ),
    )
    lbdd = assess_lbdd(
        ad_coverage_fraction=0.0,  # 172/172 outside domain
        cv_auc=0.907, n_positive=None,
    )
    overall = combine(sbdd, lbdd)

    assert sbdd["final_tier"] == "VERY LOW"
    assert lbdd["final_tier"] == "VERY LOW"
    assert overall["overall_tier"] == "VERY LOW"

    # The reasons must be traceable to the real evidence, not a black-box score.
    reasons_text = " ".join(sbdd["downgrade_reasons"] + lbdd["downgrade_reasons"])
    assert "control_representativeness" in reasons_text
    assert "applicability_domain_coverage" in reasons_text


def test_combine_uses_the_worse_of_attempted_paths():
    strong_sbdd = assess_sbdd(
        rmsd_angstrom=0.3,
        control_representativeness_report={"available": True, "similarity_coverage": 0.8,
                                             "scaffold_coverage": 0.5},
        auc=0.8, ef1pct=6.0, bedroc20=0.5,
    )
    weak_lbdd = assess_lbdd(ad_coverage_fraction=0.0, cv_auc=0.5, n_positive=5)
    overall = combine(strong_sbdd, weak_lbdd)
    assert strong_sbdd["final_tier"] == "HIGH"
    assert weak_lbdd["final_tier"] == "VERY LOW"
    assert overall["overall_tier"] == "VERY LOW"


def test_combine_with_no_paths_is_very_low():
    assert combine()["overall_tier"] == "VERY LOW"

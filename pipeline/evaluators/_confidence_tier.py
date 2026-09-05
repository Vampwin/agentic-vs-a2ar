"""
Confidence tier assessment — GRADE-style evidence grading for a virtual-
screening result, adapted from clinical medicine's Grading of Recommendations
Assessment, Development and Evaluation framework.

Design record: 00_project_control/17_CONFIDENCE_TIER_FRAMEWORK.md

Why a tier and not a percentage: a calibrated probability needs outcome data
to check it against (the way a weather forecast is checkable against how
often it actually rained). This project has no ground-truth binding data to
calibrate a confidence percentage against, and building one anyway repeats a
mistake already made once this session — isotonic calibration on 9 positive
examples reported ECE = 0.0000, which looked precise and was overfit
(ADR-AIIM-041/042). A GRADE-style tier, downgraded from a starting point by
named, itemized concerns computed from checks the audit system already runs,
does not carry that false-precision claim.

Aggregation is weakest-link, not a weighted average: a "very serious" finding
in one domain cannot be offset by a good result in another. This mirrors
hard-contract precedence used everywhere else in this system (a HALT cannot
be bought back by a good score elsewhere) and avoids the averaging-masks-a-
critical-weakness failure mode the audit layer otherwise exists to catch.
"""

from __future__ import annotations

from typing import Any, Literal

Rating = Literal["not_serious", "serious", "very_serious"]
Tier = Literal["HIGH", "MODERATE", "LOW", "VERY LOW"]

# ── Frozen cut points (ADR-AIIM-046, frozen 2026-08-07) ─────────────────────
# These are conventions, not derived constants -- exactly like
# _control_representativeness.py's thresholds. The values here are the
# fallback defaults; contracts/G8_reporting.yaml's confidence_tier_assessment
# check carries the authoritative copy, and a threshold_immutability-style
# check (G8_evaluator.py) fails closed if the two ever disagree. Do not edit
# either copy without a new ADR -- that is fault family F07
# (undeclared threshold movement) in this project's own taxonomy.

# SBDD path
DEFAULT_POSE_PASS_THRESHOLD = 2.0        # Angstrom, matches G3 rmsd_rank1 threshold
DEFAULT_POSE_MARGIN_THRESHOLD = 1.0      # Angstrom, half the pass threshold
DEFAULT_MIN_SIMILARITY_COVERAGE = 0.50   # matches G3 control_representativeness
DEFAULT_MIN_SCAFFOLD_COVERAGE = 0.10     # matches G3 control_representativeness
DEFAULT_THIN_MARGIN_MULTIPLIER = 1.5     # "thin pass" = within 1.5x the min coverage
DEFAULT_EXTERNAL_AUC_MIN = 0.60          # matches G5b auc_roc threshold
DEFAULT_EXTERNAL_EF1PCT_MIN = 3.0        # matches G5b ef_1pct threshold
DEFAULT_EXTERNAL_BEDROC20_MIN = 0.40     # matches G5b bedroc20 threshold

# LBDD path
DEFAULT_AD_COVERAGE_HIGH = 0.80          # >= this: not_serious
DEFAULT_AD_COVERAGE_LOW = 0.40           # >= this: serious; below: very_serious
DEFAULT_CV_AUC_HIGH = 0.75               # >= this: not_serious
DEFAULT_CV_AUC_LOW = 0.60                # >= this: serious; below: very_serious
DEFAULT_LBDD_EXTERNAL_AUC_MIN = 0.60
DEFAULT_IMPRECISION_MIN_RECOMMENDED = 30  # positive-class size
DEFAULT_IMPRECISION_MIN_FLOOR = 10

_TIER_ORDER: list[Tier] = ["HIGH", "MODERATE", "LOW", "VERY LOW"]
_DOWNGRADE_STEPS = {"not_serious": 0, "serious": 1, "very_serious": 2}


def _domain(name: str, rating: Rating, detail: str, value: Any = None) -> dict[str, Any]:
    return {"domain": name, "rating": rating, "detail": detail, "value": value}


def _apply_downgrades(domains: list[dict[str, Any]]) -> tuple[Tier, list[str]]:
    idx = 0
    reasons = []
    for d in domains:
        step = _DOWNGRADE_STEPS[d["rating"]]
        if step:
            idx = min(idx + step, len(_TIER_ORDER) - 1)
            reasons.append(f"{d['domain']}: {d['rating']} ({d['detail']})")
    return _TIER_ORDER[idx], reasons


# ── SBDD domain raters ──────────────────────────────────────────────────────

def rate_pose_validation_margin(rmsd_angstrom: float | None,
                                 pass_threshold: float = DEFAULT_POSE_PASS_THRESHOLD,
                                 margin_threshold: float = DEFAULT_POSE_MARGIN_THRESHOLD) -> dict[str, Any]:
    if rmsd_angstrom is None:
        return _domain("pose_validation_margin", "very_serious",
                        "no redocking RMSD available")
    if rmsd_angstrom <= margin_threshold:
        return _domain("pose_validation_margin", "not_serious",
                        f"RMSD {rmsd_angstrom:.3f} A, well inside the {pass_threshold} A threshold",
                        rmsd_angstrom)
    if rmsd_angstrom <= pass_threshold:
        return _domain("pose_validation_margin", "serious",
                        f"RMSD {rmsd_angstrom:.3f} A passes but with a thin margin "
                        f"to the {pass_threshold} A threshold", rmsd_angstrom)
    return _domain("pose_validation_margin", "very_serious",
                    f"RMSD {rmsd_angstrom:.3f} A exceeds the {pass_threshold} A threshold",
                    rmsd_angstrom)


def rate_control_representativeness(
    report: dict[str, Any] | None,
    min_similarity_coverage: float = DEFAULT_MIN_SIMILARITY_COVERAGE,
    min_scaffold_coverage: float = DEFAULT_MIN_SCAFFOLD_COVERAGE,
) -> dict[str, Any]:
    if not report or not report.get("available"):
        return _domain("control_representativeness", "very_serious",
                        "could not be measured (missing control or library SMILES)")
    sim = report.get("similarity_coverage", 0.0)
    scaf = report.get("scaffold_coverage", 0.0)
    if sim >= min_similarity_coverage or scaf >= min_scaffold_coverage:
        thin = (sim < min_similarity_coverage * DEFAULT_THIN_MARGIN_MULTIPLIER
                and scaf < min_scaffold_coverage * DEFAULT_THIN_MARGIN_MULTIPLIER)
        rating = "serious" if thin else "not_serious"
        return _domain("control_representativeness", rating,
                        f"similarity_coverage={sim:.3f}, scaffold_coverage={scaf:.3f}",
                        {"similarity_coverage": sim, "scaffold_coverage": scaf})
    return _domain("control_representativeness", "very_serious",
                    f"similarity_coverage={sim:.3f}, scaffold_coverage={scaf:.3f} — "
                    "the redocking control does not resemble the screening library",
                    {"similarity_coverage": sim, "scaffold_coverage": scaf})


def rate_structural_indirectness(ligand_dependent_assumption: bool,
                                  detail: str = "") -> dict[str, Any]:
    """
    `ligand_dependent_assumption` is a human/evaluator-recorded judgment (e.g.
    "pose validity depends on a structural water network shown to be
    ligand-specific," ADR-AIIM-036) — not something derivable purely from a
    single number, so this rater takes an explicit flag rather than
    recomputing anything.
    """
    if ligand_dependent_assumption:
        return _domain("structural_model_indirectness", "very_serious",
                        detail or "pose validation depends on an assumption known to be "
                                   "ligand-specific")
    return _domain("structural_model_indirectness", "not_serious",
                    detail or "no undocumented ligand-specific assumption identified")


def rate_external_scoring_validity(
    auc: float | None, ef1pct: float | None, bedroc20: float | None,
    auc_min: float = DEFAULT_EXTERNAL_AUC_MIN,
    ef1pct_min: float = DEFAULT_EXTERNAL_EF1PCT_MIN,
    bedroc20_min: float = DEFAULT_EXTERNAL_BEDROC20_MIN,
    skip: bool = False, skip_reason: str = "",
) -> dict[str, Any]:
    if skip:
        return _domain("external_scoring_validity", "not_serious",
                        f"not evaluated — {skip_reason}", None)
    if auc is None:
        return _domain("external_scoring_validity", "very_serious", "no external benchmark available")
    passed = [m for m, v, t in (("auc", auc, auc_min), ("ef1pct", ef1pct, ef1pct_min),
                                 ("bedroc20", bedroc20, bedroc20_min)) if v is not None and v >= t]
    n_metrics = sum(1 for v in (auc, ef1pct, bedroc20) if v is not None)
    if n_metrics and len(passed) == n_metrics:
        return _domain("external_scoring_validity", "not_serious",
                        f"AUC={auc}, EF1%={ef1pct}, BEDROC20={bedroc20} all clear thresholds",
                        {"auc": auc, "ef1pct": ef1pct, "bedroc20": bedroc20})
    if len(passed) > 0:
        return _domain("external_scoring_validity", "serious",
                        f"AUC={auc}, EF1%={ef1pct}, BEDROC20={bedroc20} — only some metrics clear thresholds",
                        {"auc": auc, "ef1pct": ef1pct, "bedroc20": bedroc20})
    return _domain("external_scoring_validity", "very_serious",
                    f"AUC={auc}, EF1%={ef1pct}, BEDROC20={bedroc20} — thresholds not met",
                    {"auc": auc, "ef1pct": ef1pct, "bedroc20": bedroc20})


def assess_sbdd(
    rmsd_angstrom: float | None,
    control_representativeness_report: dict[str, Any] | None,
    ligand_dependent_assumption: bool = False,
    structural_indirectness_detail: str = "",
    auc: float | None = None, ef1pct: float | None = None, bedroc20: float | None = None,
    pose_pass_threshold: float = DEFAULT_POSE_PASS_THRESHOLD,
    pose_margin_threshold: float = DEFAULT_POSE_MARGIN_THRESHOLD,
    min_similarity_coverage: float = DEFAULT_MIN_SIMILARITY_COVERAGE,
    min_scaffold_coverage: float = DEFAULT_MIN_SCAFFOLD_COVERAGE,
    external_auc_min: float = DEFAULT_EXTERNAL_AUC_MIN,
    external_ef1pct_min: float = DEFAULT_EXTERNAL_EF1PCT_MIN,
    external_bedroc20_min: float = DEFAULT_EXTERNAL_BEDROC20_MIN,
) -> dict[str, Any]:
    control_domain = rate_control_representativeness(
        control_representativeness_report, min_similarity_coverage, min_scaffold_coverage,
    )
    # An external-AUC computed on a benchmark chemically unlike the screening
    # library is not informative about the library, so skip rather than
    # score it once control representativeness has already floored the tier.
    skip_external = control_domain["rating"] == "very_serious"
    domains = [
        rate_pose_validation_margin(rmsd_angstrom, pose_pass_threshold, pose_margin_threshold),
        control_domain,
        rate_structural_indirectness(ligand_dependent_assumption, structural_indirectness_detail),
        rate_external_scoring_validity(
            auc, ef1pct, bedroc20, external_auc_min, external_ef1pct_min, external_bedroc20_min,
            skip=skip_external,
            skip_reason="control does not represent the screening library "
                        "(external AUC would not be informative)",
        ),
    ]
    tier, reasons = _apply_downgrades(domains)
    return {"path": "SBDD", "starting_tier": "HIGH", "final_tier": tier,
            "domains": domains, "downgrade_reasons": reasons}


# ── LBDD domain raters ──────────────────────────────────────────────────────

def rate_applicability_domain(coverage_fraction: float | None,
                               high: float = DEFAULT_AD_COVERAGE_HIGH,
                               low: float = DEFAULT_AD_COVERAGE_LOW) -> dict[str, Any]:
    if coverage_fraction is None:
        return _domain("applicability_domain_coverage", "very_serious",
                        "applicability domain not assessed")
    if coverage_fraction >= high:
        return _domain("applicability_domain_coverage", "not_serious",
                        f"{coverage_fraction:.1%} of the library within the training "
                        "chemical space", coverage_fraction)
    if coverage_fraction >= low:
        return _domain("applicability_domain_coverage", "serious",
                        f"{coverage_fraction:.1%} of the library within domain",
                        coverage_fraction)
    return _domain("applicability_domain_coverage", "very_serious",
                    f"only {coverage_fraction:.1%} of the library within the training "
                    "chemical space", coverage_fraction)


def rate_internal_validation(cv_auc: float | None, n_positive: int | None,
                              fold_degenerate: bool = False,
                              high: float = DEFAULT_CV_AUC_HIGH,
                              low: float = DEFAULT_CV_AUC_LOW) -> dict[str, Any]:
    if cv_auc is None:
        return _domain("internal_validation_quality", "very_serious", "no CV performance available")
    if fold_degenerate:
        return _domain("internal_validation_quality", "very_serious",
                        f"CV AUC={cv_auc:.3f} but cross-validation folds were degenerate "
                        f"(n_positive={n_positive})", cv_auc)
    if cv_auc >= high:
        return _domain("internal_validation_quality", "not_serious", f"CV AUC={cv_auc:.3f}", cv_auc)
    if cv_auc >= low:
        return _domain("internal_validation_quality", "serious", f"CV AUC={cv_auc:.3f}", cv_auc)
    return _domain("internal_validation_quality", "very_serious", f"CV AUC={cv_auc:.3f}", cv_auc)


def rate_external_validation_lbdd(auc: float | None, benchmark_representative: bool = True,
                                   auc_min: float = DEFAULT_LBDD_EXTERNAL_AUC_MIN) -> dict[str, Any]:
    if auc is None:
        return _domain("external_validation", "very_serious", "no independent benchmark available")
    if not benchmark_representative:
        return _domain("external_validation", "very_serious",
                        f"AUC={auc:.3f} but the benchmark is not chemically representative "
                        "of the screening library", auc)
    if auc >= auc_min:
        return _domain("external_validation", "not_serious", f"AUC={auc:.3f}", auc)
    return _domain("external_validation", "very_serious", f"AUC={auc:.3f} below {auc_min}", auc)


def rate_imprecision(n_positive: int | None,
                      min_recommended: int = DEFAULT_IMPRECISION_MIN_RECOMMENDED,
                      min_floor: int = DEFAULT_IMPRECISION_MIN_FLOOR) -> dict[str, Any]:
    if n_positive is None:
        return _domain("imprecision", "very_serious", "positive-class size unknown")
    if n_positive >= min_recommended:
        return _domain("imprecision", "not_serious", f"n_positive={n_positive}", n_positive)
    if n_positive >= min_floor:
        return _domain("imprecision", "serious",
                        f"n_positive={n_positive}, below the {min_recommended} rule-of-thumb minimum",
                        n_positive)
    return _domain("imprecision", "very_serious",
                    f"n_positive={n_positive} is too small for the metric to be trustworthy",
                    n_positive)


def assess_lbdd(
    ad_coverage_fraction: float | None,
    cv_auc: float | None, n_positive: int | None, fold_degenerate: bool = False,
    external_auc: float | None = None, external_benchmark_representative: bool = True,
    ad_coverage_high: float = DEFAULT_AD_COVERAGE_HIGH,
    ad_coverage_low: float = DEFAULT_AD_COVERAGE_LOW,
    cv_auc_high: float = DEFAULT_CV_AUC_HIGH,
    cv_auc_low: float = DEFAULT_CV_AUC_LOW,
    external_auc_min: float = DEFAULT_LBDD_EXTERNAL_AUC_MIN,
    imprecision_min_recommended: int = DEFAULT_IMPRECISION_MIN_RECOMMENDED,
    imprecision_min_floor: int = DEFAULT_IMPRECISION_MIN_FLOOR,
) -> dict[str, Any]:
    domains = [
        rate_applicability_domain(ad_coverage_fraction, ad_coverage_high, ad_coverage_low),
        rate_internal_validation(cv_auc, n_positive, fold_degenerate, cv_auc_high, cv_auc_low),
        rate_external_validation_lbdd(external_auc, external_benchmark_representative, external_auc_min),
        rate_imprecision(n_positive, imprecision_min_recommended, imprecision_min_floor),
    ]
    tier, reasons = _apply_downgrades(domains)
    return {"path": "LBDD", "starting_tier": "HIGH", "final_tier": tier,
            "domains": domains, "downgrade_reasons": reasons}


def combine(*assessments: dict[str, Any]) -> dict[str, Any]:
    """
    Overall tier across every path actually attempted. Still weakest-link:
    if only one path was attempted, its tier is the overall tier; if more
    than one, the worse of the attempted paths' tiers governs, because a
    strong LBDD result does not make an unvalidated SBDD ranking trustworthy
    or vice versa — they answer different questions and both are being
    reported.
    """
    present = [a for a in assessments if a]
    if not present:
        return {"overall_tier": "VERY LOW", "paths": [], "reason": "no path was assessed"}
    worst_idx = max(_TIER_ORDER.index(a["final_tier"]) for a in present)
    return {"overall_tier": _TIER_ORDER[worst_idx], "paths": present}

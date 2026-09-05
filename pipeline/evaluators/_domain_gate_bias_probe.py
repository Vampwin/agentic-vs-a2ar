"""
Domain-gate bias probe — is a single-descriptor-space "0% coverage" verdict
being over-interpreted?

E2 (2026-08-15, ADR-AIIM-073) found that SNP-Lib's 0% applicability-domain
coverage against the frozen ChEMBL A2AR training set reproduces across six
independent structural descriptor spaces (ECFP4, FCFP4, topological torsion,
Gobbi 2D pharmacophore, ErG, AtomPair) — including two designed specifically
for scaffold hopping — while the same library sits inside the training set's
*physicochemical* space (64.4% coverage). That pattern is the expected
signature of a structurally novel library, not evidence of a defective
fingerprint choice or an unscreenable library: a genuine analogue of a known
actives class would show high coverage in at least one structural space, not
zero in six.

A verdict computed in a single descriptor space cannot tell these two cases
apart. This module does not decide whether a ranking should be released — it
flags when a domain-gate conclusion rests on one space, or when a structural/
physicochemical divergence exists and has not been surfaced alongside the
coverage number, so the human decision at G6/G8 is made with the full
picture rather than with whatever fingerprint was computed first.
"""

from __future__ import annotations

from typing import Any

# A verdict computed in only one space cannot rule out a descriptor artefact.
DEFAULT_MIN_SPACES = 3
# Cheap divergence heuristic: near-zero everywhere structural, but well inside
# domain physicochemically. Neither cutoff is a hard scientific constant; both
# are reported alongside the raw per-space numbers so a human can re-judge.
_STRUCTURAL_NEAR_ZERO = 0.05
_PHYSCHEM_SUBSTANTIAL = 0.30


def _is_physchem(space_name: str) -> bool:
    return "phys" in space_name.lower()


def assess(
    coverage_by_space: dict[str, dict[str, Any]],
    min_spaces: int = DEFAULT_MIN_SPACES,
) -> dict[str, Any]:
    """
    coverage_by_space: {space_name: {"coverage_fraction": float, ...}, ...},
    one entry per descriptor space (the shape produced by
    10_new_system/scripts/ad_descriptor_sweep.py's per-space rows).

    Returns a report; the caller (a G6-style evaluator) applies whatever
    pass/warn policy its contract declares. This function never itself
    decides pass/fail — only what is true about the evidence supplied.
    """
    if not coverage_by_space:
        return {"available": False,
                "reason": "no per-space applicability-domain data supplied"}

    fractions = {
        name: float(v.get("coverage_fraction", 0.0))
        for name, v in coverage_by_space.items()
    }
    n_spaces = len(fractions)

    structural = {k: v for k, v in fractions.items() if not _is_physchem(k)}
    physchem = {k: v for k, v in fractions.items() if _is_physchem(k)}

    structural_max = max(structural.values()) if structural else None
    physchem_max = max(physchem.values()) if physchem else None

    # The signature this module exists to detect: structurally novel but not
    # physicochemically unusual. Divergence — not either number alone — is
    # the signal. A library that is genuinely out-of-domain everywhere,
    # structural and physicochemical alike, is a different and equally real
    # finding, and this check must not obscure that case by always firing.
    divergent = (
        structural_max is not None and physchem_max is not None
        and structural_max <= _STRUCTURAL_NEAR_ZERO
        and physchem_max >= _PHYSCHEM_SUBSTANTIAL
    )

    report: dict[str, Any] = {
        "available": True,
        "n_spaces": n_spaces,
        "spaces": sorted(fractions),
        "coverage_by_space": fractions,
        "has_structural_space": bool(structural),
        "has_physicochemical_space": bool(physchem),
        "structural_max_coverage": structural_max,
        "physicochemical_max_coverage": physchem_max,
        "enough_spaces": n_spaces >= min_spaces,
        "min_spaces_required": min_spaces,
        "structural_physchem_divergence": divergent,
    }
    if divergent:
        report["interpretation"] = (
            "Coverage is near zero in every structural descriptor space checked "
            "while the library sits inside the training set's physicochemical "
            "range. This is the expected signature of a structurally novel "
            "library, not evidence that the library is chemically unscreenable. "
            "A HALT or exclusion decision resting on this applicability-domain "
            "result alone requires corroborating structure-based evidence "
            "(e.g. G3 retrospective_enrichment) before being treated as final."
        )
    elif structural_max is not None and structural_max <= _STRUCTURAL_NEAR_ZERO:
        report["interpretation"] = (
            "Coverage is near zero in the structural space(s) checked, and no "
            "physicochemical space was supplied (or it is also low), so no "
            "structural/physicochemical divergence can be confirmed. This "
            "result alone does not distinguish a genuinely out-of-domain "
            "library from a structurally novel one -- report both possibilities."
        )
    return report


# --------------------------------------------------------------------------
# False refusal on a known-active panel.
#
# The G6 contract has always named this mechanism -- run a panel of compounds
# known to bind the target through the applicability-domain gate and measure how
# many it refuses -- and until 2026-08-29 no code implemented it. `assess()`
# above audits how many descriptor spaces a verdict rests on; it says nothing
# about whether the gate is right. An internal review recorded that gap as an
# instance of this paper's own defect class, a contract asserting more than its
# evaluator does, and this closes it (M3, ADR-AIIM-147).
#
# The measurement only means anything if the panel is disjoint from the set that
# defines the domain. Every compound in the reference set is inside the domain by
# construction, so scoring the reference set against itself returns a refusal
# rate near zero and reads as a pass. The caller is responsible for that
# disjointness; this function refuses to report a rate without being told, in
# writing, how it was achieved.
FALSE_REFUSAL_CONCERNING = 0.50   # a gate refusing most known actives is broken


def assess_false_refusal(
    in_domain_by_space: dict[str, list[bool]],
    panel_disjoint_from_reference: str,
    panel_description: str,
    concerning_at: float = FALSE_REFUSAL_CONCERNING,
) -> dict[str, Any]:
    """
    in_domain_by_space: {space_name: [bool, ...]}, one entry per descriptor
        space, each list one flag per panel compound, produced by applying the
        gate's own frozen rule.
    panel_disjoint_from_reference: prose stating how the panel was kept out of
        the set that defines the domain. Required and non-empty; a rate measured
        without that separation is not a false-refusal rate.
    panel_description: what the compounds are and where they came from.

    A false refusal is a compound known to be active that the gate places
    out of domain. The rate is reported per space and never averaged across
    spaces: the spaces are not commensurable, which is the finding this whole
    probe exists to support.
    """
    if not panel_disjoint_from_reference.strip():
        return {"available": False,
                "reason": "no statement of how the panel was separated from the "
                          "reference set; a refusal rate measured on compounds "
                          "that define the domain is zero by construction"}
    if not in_domain_by_space:
        return {"available": False, "reason": "no per-space flags supplied"}

    sizes = {len(v) for v in in_domain_by_space.values()}
    if len(sizes) != 1:
        return {"available": False,
                "reason": f"spaces disagree on panel size: {sorted(sizes)}"}
    n = sizes.pop()
    if n == 0:
        return {"available": False, "reason": "panel is empty"}

    by_space = {}
    for name, flags in in_domain_by_space.items():
        refused = sum(1 for f in flags if not f)
        by_space[name] = {
            "n_panel": n,
            "n_refused": refused,
            "false_refusal_rate": round(refused / n, 4),
            "is_physicochemical": _is_physchem(name),
        }

    structural = {k: v for k, v in by_space.items() if not v["is_physicochemical"]}
    physchem = {k: v for k, v in by_space.items() if v["is_physicochemical"]}
    worst = max(by_space, key=lambda k: by_space[k]["false_refusal_rate"])
    best = min(by_space, key=lambda k: by_space[k]["false_refusal_rate"])

    report: dict[str, Any] = {
        "available": True,
        "panel": panel_description,
        "panel_size": n,
        "separation": panel_disjoint_from_reference,
        "by_space": by_space,
        "worst_space": worst,
        "best_space": best,
        "worst_rate": by_space[worst]["false_refusal_rate"],
        "best_rate": by_space[best]["false_refusal_rate"],
        "spread": round(by_space[worst]["false_refusal_rate"]
                        - by_space[best]["false_refusal_rate"], 4),
        "concerning_threshold": concerning_at,
        "spaces_above_threshold": sorted(
            k for k, v in by_space.items()
            if v["false_refusal_rate"] >= concerning_at),
    }

    n_bad = len(report["spaces_above_threshold"])
    if n_bad == len(by_space):
        report["interpretation"] = (
            "Every descriptor space refuses at least half of a panel of known "
            "actives that it has never seen. The gate is not distinguishing "
            "unfamiliar from inactive, and an exclusion resting on it carries "
            "no evidence about binding.")
    elif n_bad:
        report["interpretation"] = (
            f"{n_bad} of {len(by_space)} spaces refuse at least half of a panel "
            f"of known actives; {report['best_space']} refuses "
            f"{report['best_rate']:.1%}. A verdict computed in one space is a "
            f"verdict about that space, and which space was chosen matters more "
            f"than the number it returns.")
    else:
        report["interpretation"] = (
            "No space refuses half or more of the known-active panel. The gate "
            "admits unfamiliar actives at a rate that does not by itself "
            "explain an exclusion.")
    if structural and physchem:
        s_worst = max(v["false_refusal_rate"] for v in structural.values())
        p_worst = max(v["false_refusal_rate"] for v in physchem.values())
        report["structural_vs_physicochemical"] = {
            "structural_worst": s_worst,
            "physicochemical_worst": p_worst,
            "difference": round(s_worst - p_worst, 4),
        }
    return report

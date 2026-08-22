"""
G5b evaluator — external validation performance criteria.
Conditional on G5a passing.
"""

from __future__ import annotations
from pathlib import Path
from typing import Any

import numpy as np


def evaluate(contract: dict[str, Any], inputs: dict[str, Any], registry: Any) -> dict[str, Any]:
    passed, failed, warned, evidence = [], [], [], {}
    critical_failures = []

    thresholds = contract.get("thresholds", {})
    auc_min = thresholds.get("pass_auc_min", 0.60)
    ef1pct_min = thresholds.get("pass_ef1pct_min", 3.0)
    bedroc_min = thresholds.get("pass_bedroc20_min", 0.40)

    # Verify threshold immutability
    cid = "threshold_immutability"
    registered_auc = inputs.get("registered_auc_min_threshold")
    if registered_auc is not None and abs(registered_auc - auc_min) > 1e-9:
        failed.append(cid)
        critical_failures.append(cid)
        evidence[cid] = {
            "registered": registered_auc, "current": auc_min,
            "error": "Threshold changed since registration",
        }
    else:
        passed.append(cid)
        evidence[cid] = {"auc_min": auc_min, "ef1pct_min": ef1pct_min, "bedroc_min": bedroc_min}

    # Performance metrics
    auc = inputs.get("auc_roc")
    ef1 = inputs.get("ef_1pct")
    bedroc = inputs.get("bedroc20")

    for cid, value, threshold, name in [
        ("auc_roc", auc, auc_min, "AUC"),
        ("ef_1pct", ef1, ef1pct_min, "EF@1%"),
        ("bedroc20", bedroc, bedroc_min, "BEDROC20"),
    ]:
        if value is None:
            warned.append(cid)
            evidence[cid] = {"warning": f"{name} not in inputs"}
        else:
            ok = value >= threshold
            evidence[cid] = {
                "value": round(float(value), 4),
                "threshold": threshold,
                "pass": ok,
                "ci": inputs.get(f"{cid}_ci"),
            }
            (passed if ok else failed).append(cid)

    # Reporting language.
    # Declared critical in the contract: if benchmark_role cannot be resolved we
    # cannot know whether "held-out external validation" language is permitted,
    # so this fails closed. A critical check that could not be evaluated must
    # never be reported as passed.
    cid = "reporting_language_check"
    role = inputs.get("benchmark_role", "")
    if "Option_A" in str(role) or "development" in str(role).lower():
        evidence[cid] = {
            "required_caveat": (
                "This AUC/EF was computed on a development benchmark "
                "(ChEMBL A2AR set, Options A). Results are exploratory and not a "
                "held-out external validation."
            )
        }
        passed.append(cid)
    elif "Option_B" in str(role) or "heldout" in str(role).lower():
        evidence[cid] = {"reporting": "held-out external validation"}
        passed.append(cid)
    else:
        failed.append(cid)
        critical_failures.append(cid)
        evidence[cid] = {
            "error": "benchmark_role not set; cannot determine permitted reporting language"
        }

    # Critical breaches (threshold tampering, unresolvable reporting language)
    # halt outright; metric shortfalls are high/medium severity and escalate to
    # a human per the contract's on_fail policy.
    if critical_failures:
        state = "HALT"
    elif failed:
        state = "ESCALATE"
    else:
        state = "PASS"
    return {
        "state": state,
        "passed": passed,
        "failed": failed,
        "warned": warned,
        "critical_failures": critical_failures,
        "evidence": evidence,
        "summary": (
            f"{len(passed)} checks passed, {len(failed)} failed, {len(warned)} warned. "
            f"State: {state}"
        ),
    }

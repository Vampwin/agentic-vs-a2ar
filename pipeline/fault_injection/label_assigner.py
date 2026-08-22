"""
label_assigner.py — Development-set label assignment for risk-monitor training.

Converts feasibility pilot outcomes to binary `unsafe_release` labels suitable
for training the risk monitor under the approved Gate 4 label criteria
(ADR-AIIM-019, Gate 4 approved in session log §10):

    unsafe_release = 1  if fault case AND escaped the reference condition
    unsafe_release = 0  if clean case, OR fault case caught in reference condition

The reference condition defaults to "B" (contracts only), so the monitor learns
to predict which cases slip through contracts.  This supports the primary
A-vs-D contrast: under condition A (hybrid), the monitor flags cases that
condition B (contracts alone) would miss.

FORBIDDEN inputs — this module must never receive:
  - cAMP values, assay outcomes, or experimental-condition identifiers
  - Held-out manifest labels (H01–H12; sealed until definitive model frozen)
  - Post-release outcomes or any measurement taken after the audit decision

Only development-family pilot outcomes (from run_pilot()) and ground-truth
injection metadata (escaped / not escaped, family_id, case_type) may be used.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any

import pandas as pd

logger = logging.getLogger(__name__)

# Gate 4 approved reference condition for deriving unsafe_release labels.
# Changing this requires a new ADR and human approval.
_DEFAULT_REFERENCE_CONDITION: str = "B"

_FORBIDDEN_TERMS: frozenset[str] = frozenset({
    "camp", "assay", "ic50", "ec50", "ki", "active", "inactive",
    "held_out", "heldout", "h01", "h02", "h03", "h04", "h05",
    "h06", "h07", "h08", "h09", "h10", "h11", "h12",
})


@dataclass(frozen=True)
class LabelRecord:
    instance_id: str
    family_id: str
    condition: str
    case_type: str        # "fault" | "clean"
    unsafe_release: int   # 0 or 1
    group: str            # family_id used as GroupKFold key


def _check_no_assay_leakage(outcomes: list[dict[str, Any]]) -> None:
    for o in outcomes:
        for key, val in o.items():
            key_lower = str(key).lower()
            val_lower = str(val).lower() if val is not None else ""
            for term in _FORBIDDEN_TERMS:
                if term in key_lower or term in val_lower:
                    raise ValueError(
                        f"label_assigner: possible assay leakage — "
                        f"term {term!r} found in outcome field {key!r}={val!r}. "
                        "Held-out labels and assay data are forbidden here."
                    )


def assign_labels(
    outcomes: list[dict[str, Any]],
    reference_condition: str = _DEFAULT_REFERENCE_CONDITION,
) -> pd.DataFrame:
    """
    Assign Gate 4 `unsafe_release` labels from pilot outcomes.

    Parameters
    ----------
    outcomes:
        List of outcome dicts from PilotResults.outcomes (as produced by
        run_pilot() and serialised in the pilot JSON).  Each dict must contain:
          - ``case_type``: "fault" or "clean"
          - ``condition``: e.g. "B", "D"
          - ``family_id``: development fault family string
          - ``instance_id``: unique case identifier
          - ``escaped``: bool — True if the case was not contained
          - ``error``: bool — True if the evaluator raised an exception
    reference_condition:
        Condition used to derive labels (default "B": contracts only).
        Must appear in the provided outcomes.

    Returns
    -------
    pd.DataFrame with columns:
        instance_id, family_id, condition, case_type, unsafe_release (int 0/1), group.
    Rows with ``error=True`` are excluded (evaluator could not determine outcome;
    omitting them is the fail-closed choice).

    Raises
    ------
    ValueError
        If no outcomes match the reference condition, or if forbidden terms
        (assay data, held-out label identifiers) are detected in the outcomes.
    """
    _check_no_assay_leakage(outcomes)

    ref_outcomes = [
        o for o in outcomes
        if o.get("condition") == reference_condition and not o.get("error", False)
    ]
    if not ref_outcomes:
        raise ValueError(
            f"label_assigner: no non-errored outcomes found for "
            f"reference condition {reference_condition!r}. "
            "Ensure the pilot was run with this condition and produced results."
        )

    records: list[LabelRecord] = []
    for o in ref_outcomes:
        case_type = o.get("case_type", "clean")
        escaped = bool(o.get("escaped", False))

        # Gate 4 label criteria (human-approved):
        # unsafe_release = 1  iff  fault case AND escaped reference condition
        # unsafe_release = 0  iff  clean case, OR fault case caught by reference condition
        if case_type == "fault" and escaped:
            label = 1
        else:
            label = 0

        family_id = str(o.get("family_id", ""))
        records.append(LabelRecord(
            instance_id=str(o.get("instance_id", "")),
            family_id=family_id,
            condition=reference_condition,
            case_type=case_type,
            unsafe_release=label,
            group=family_id,
        ))

    df = pd.DataFrame([r.__dict__ for r in records])

    n_positive = int((df["unsafe_release"] == 1).sum())
    n_total = len(df)
    logger.info(
        "label_assigner: assigned %d labels (%d positive / %d negative) "
        "from condition %r outcomes.",
        n_total, n_positive, n_total - n_positive, reference_condition,
    )
    if n_positive == 0:
        logger.warning(
            "label_assigner: no positive (unsafe_release=1) labels assigned. "
            "The monitor cannot be trained with all-zero labels."
        )
    if n_positive == n_total:
        logger.warning(
            "label_assigner: all labels are positive (unsafe_release=1). "
            "Check that clean cases are included in the pilot outcomes."
        )

    return df


def summarise_labels(df: pd.DataFrame) -> dict[str, Any]:
    """
    Return a summary dict for logging or reporting.

    Does not expose family-level escape rates (those could encode condition
    identifiers indirectly); only aggregate statistics are returned.
    """
    n_total = len(df)
    n_positive = int((df["unsafe_release"] == 1).sum())
    n_fault = int((df["case_type"] == "fault").sum())
    n_clean = int((df["case_type"] == "clean").sum())
    n_groups = int(df["group"].nunique())

    return {
        "n_total": n_total,
        "n_positive": n_positive,
        "n_negative": n_total - n_positive,
        "positive_rate": round(n_positive / max(n_total, 1), 4),
        "n_fault_cases": n_fault,
        "n_clean_cases": n_clean,
        "n_groups": n_groups,
    }

"""
SEER (Silent Error Escape Rate) calculator.

SEER = fault_cases_that_escaped_detection / all_evaluable_fault_cases

Pseudoreplication rule (from experimental-design skill):
    The experimental unit is the FAULT CASE, not the individual run.
    If a fault case is run N times (for stochasticity estimates), it still
    contributes 1 to the denominator, not N. Multiple runs of the same case
    characterize within-case variability only.

The acting agent does not compute SEER for the held-out test set — that is
the custodian's responsibility via the heldout_interface.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd


@dataclass
class FaultCaseResult:
    fault_id: str
    fault_class: str
    gate_triggered: str          # Which gate (G0–G8) the fault should trigger
    gate_detected: str | None    # Which gate actually flagged it; None = escaped
    detected: bool
    n_runs: int = 1              # How many times this fault case was run
    run_outcomes: list[str] = field(default_factory=list)  # Per-run state values


@dataclass
class SEERResult:
    total_fault_cases: int
    detected: int
    escaped: int
    seer: float
    seer_ci_lower: float
    seer_ci_upper: float
    bootstrap_n: int
    per_gate_detection: dict[str, dict[str, int]]
    per_class_detection: dict[str, dict[str, int]]
    note: str = ""


def compute_seer(
    fault_results: list[FaultCaseResult],
    bootstrap_n: int = 1000,
    rng_seed: int = 20260729,
) -> SEERResult:
    """
    Compute SEER from a list of fault case outcomes.

    One fault case = one experimental unit regardless of n_runs.
    Bootstrap CI uses case-level resampling (not run-level).
    """
    import numpy as np

    n = len(fault_results)
    if n == 0:
        raise ValueError("No fault results provided")

    detected_flags = [int(r.detected) for r in fault_results]
    n_detected = sum(detected_flags)
    n_escaped = n - n_detected
    seer = n_escaped / n

    # Bootstrap CI for SEER (case-level resampling)
    rng = np.random.default_rng(rng_seed)
    boot_seers = []
    for _ in range(bootstrap_n):
        sample = rng.choice(detected_flags, size=n, replace=True)
        boot_seers.append(1.0 - float(sample.mean()))
    ci_lower = float(np.percentile(boot_seers, 2.5))
    ci_upper = float(np.percentile(boot_seers, 97.5))

    # Per-gate breakdown
    per_gate: dict[str, dict[str, int]] = {}
    for r in fault_results:
        gate = r.gate_triggered
        per_gate.setdefault(gate, {"detected": 0, "escaped": 0})
        if r.detected:
            per_gate[gate]["detected"] += 1
        else:
            per_gate[gate]["escaped"] += 1

    # Per-class breakdown
    per_class: dict[str, dict[str, int]] = {}
    for r in fault_results:
        cls = r.fault_class
        per_class.setdefault(cls, {"detected": 0, "escaped": 0})
        if r.detected:
            per_class[cls]["detected"] += 1
        else:
            per_class[cls]["escaped"] += 1

    return SEERResult(
        total_fault_cases=n,
        detected=n_detected,
        escaped=n_escaped,
        seer=round(seer, 4),
        seer_ci_lower=round(ci_lower, 4),
        seer_ci_upper=round(ci_upper, 4),
        bootstrap_n=bootstrap_n,
        per_gate_detection=per_gate,
        per_class_detection=per_class,
    )


def load_fault_results_csv(path: Path) -> list[FaultCaseResult]:
    """
    Load fault results from a CSV produced by the heldout custodian.
    Required columns: fault_id, fault_class, gate_triggered, gate_detected, detected.
    """
    df = pd.read_csv(path)
    required = {"fault_id", "fault_class", "gate_triggered", "gate_detected", "detected"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns in fault results CSV: {missing}")

    results = []
    for _, row in df.iterrows():
        results.append(FaultCaseResult(
            fault_id=str(row["fault_id"]),
            fault_class=str(row["fault_class"]),
            gate_triggered=str(row["gate_triggered"]),
            gate_detected=None if pd.isna(row["gate_detected"]) else str(row["gate_detected"]),
            detected=bool(row["detected"]),
            n_runs=int(row.get("n_runs", 1)),
        ))
    return results


def save_seer_result(result: SEERResult, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    data = {
        "total_fault_cases": result.total_fault_cases,
        "detected": result.detected,
        "escaped": result.escaped,
        "seer": result.seer,
        "seer_95ci": [result.seer_ci_lower, result.seer_ci_upper],
        "bootstrap_n": result.bootstrap_n,
        "per_gate_detection": result.per_gate_detection,
        "per_class_detection": result.per_class_detection,
        "note": result.note,
        "pseudoreplication_note": (
            "SEER denominator = unique fault cases. "
            "Repeated runs of the same case are NOT separate cases."
        ),
    }
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)

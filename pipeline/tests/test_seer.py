"""
Known-answer tests for the SEER calculator.
Verifies pseudoreplication rule and bootstrap CI correctness.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))
from benchmark.scoring.seer_calculator import FaultCaseResult, compute_seer

numpy = pytest.importorskip("numpy", reason="numpy not installed")


def make_results(n_total: int, n_escaped: int) -> list[FaultCaseResult]:
    results = []
    for i in range(n_total):
        detected = i >= n_escaped  # first n_escaped are escaped
        results.append(FaultCaseResult(
            fault_id=f"F{i:03d}",
            fault_class="test_class",
            gate_triggered="G0",
            gate_detected=None if not detected else "G0",
            detected=detected,
        ))
    return results


def test_all_detected_seer_zero():
    results = make_results(10, 0)
    seer = compute_seer(results, bootstrap_n=200)
    assert seer.seer == 0.0
    assert seer.escaped == 0
    assert seer.detected == 10


def test_all_escaped_seer_one():
    results = make_results(10, 10)
    seer = compute_seer(results, bootstrap_n=200)
    assert seer.seer == 1.0
    assert seer.escaped == 10
    assert seer.detected == 0


def test_partial_escape():
    # 8 detected, 2 escaped out of 10
    results = make_results(10, 2)
    seer = compute_seer(results, bootstrap_n=500)
    assert abs(seer.seer - 0.2) < 1e-9
    assert seer.seer_ci_lower < seer.seer < seer.seer_ci_upper


def test_pseudoreplication_unit_is_case():
    """n_runs does not inflate the denominator."""
    results = [
        FaultCaseResult("F001", "hash_mismatch", "G0", "G0", detected=True, n_runs=10),
        FaultCaseResult("F002", "hash_mismatch", "G0", None, detected=False, n_runs=10),
    ]
    seer = compute_seer(results, bootstrap_n=200)
    assert seer.total_fault_cases == 2  # NOT 20
    assert seer.escaped == 1
    assert abs(seer.seer - 0.5) < 1e-9


def test_bootstrap_ci_width_shrinks_with_more_cases():
    small = make_results(5, 1)
    large = make_results(50, 10)
    s_small = compute_seer(small, bootstrap_n=500)
    s_large = compute_seer(large, bootstrap_n=500)
    ci_small = s_small.seer_ci_upper - s_small.seer_ci_lower
    ci_large = s_large.seer_ci_upper - s_large.seer_ci_lower
    assert ci_small > ci_large  # larger n → narrower CI


def test_per_gate_breakdown():
    results = [
        FaultCaseResult("F001", "hash_mismatch", "G0", "G0", detected=True),
        FaultCaseResult("F002", "hash_mismatch", "G0", None, detected=False),
        FaultCaseResult("F003", "label_leak", "G5a", "G5a", detected=True),
    ]
    seer = compute_seer(results, bootstrap_n=200)
    assert seer.per_gate_detection["G0"]["detected"] == 1
    assert seer.per_gate_detection["G0"]["escaped"] == 1
    assert seer.per_gate_detection["G5a"]["detected"] == 1


def test_empty_raises():
    with pytest.raises(ValueError):
        compute_seer([], bootstrap_n=200)

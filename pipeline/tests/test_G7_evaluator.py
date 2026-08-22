"""
Tests for G7_evaluator.py -- the MD/energy refinement audit. Written
alongside the first real implementation (previously a stub that only
handled the no_MD skip case), exercising both the skip path and a real
audit against a small, real-shaped MD summary table.
"""

from __future__ import annotations

import pandas as pd
import pytest

from evaluators import G7_evaluator

_CONTRACT = {
    "checks": [
        {"id": "trajectory_completion", "severity": "critical"},
        {"id": "equilibration_verified", "severity": "high"},
        {"id": "replicate_consistency", "severity": "medium"},
        {"id": "claim_boundary", "severity": "critical"},
        {"id": "md_not_used_for_ranking", "severity": "critical"},
        {"id": "seed_recorded", "severity": "high"},
    ]
}


def test_no_md_is_skipped():
    result = G7_evaluator.evaluate(_CONTRACT, {"md_decision_recorded": "no_MD"}, registry=None)
    assert result["state"] == "SKIPPED_BY_PROTOCOL"


def test_trajectory_completion_passes_when_all_complete(tmp_path):
    df = pd.DataFrame({"compound_id": ["A", "B"], "status": ["completed", "completed"],
                        "prod_ns": [50.0, 50.0]})
    csv = tmp_path / "traj.csv"
    df.to_csv(csv, index=False)
    result = G7_evaluator.evaluate(_CONTRACT, {
        "md_decision_recorded": "MD_run", "md_completion_log": str(csv),
    }, registry=None)
    assert "trajectory_completion" in result["passed"]


def test_trajectory_completion_fails_on_incomplete_run(tmp_path):
    df = pd.DataFrame({"compound_id": ["A", "B"], "status": ["completed", "failed"],
                        "prod_ns": [50.0, 12.0]})
    csv = tmp_path / "traj.csv"
    df.to_csv(csv, index=False)
    result = G7_evaluator.evaluate(_CONTRACT, {
        "md_decision_recorded": "MD_run", "md_completion_log": str(csv),
    }, registry=None)
    assert "trajectory_completion" in result["failed"]
    assert "trajectory_completion" in result["critical_failures"]
    assert result["state"] == "HALT"


def test_equilibration_fails_closed_without_timeseries_data():
    """
    A trajectory-level summary statistic (mean RMSD) is not the same
    evidence as a per-frame plateau check -- this must fail closed, not
    pass just because *some* MD data exists.
    """
    result = G7_evaluator.evaluate(_CONTRACT, {"md_decision_recorded": "MD_run"}, registry=None)
    assert "equilibration_verified" in result["failed"]


def test_replicate_consistency_warns_without_replicates(tmp_path):
    df = pd.DataFrame({"compound_id": ["A", "B"], "status": ["completed", "completed"],
                        "prod_ns": [50.0, 50.0]})
    csv = tmp_path / "traj.csv"
    df.to_csv(csv, index=False)
    result = G7_evaluator.evaluate(_CONTRACT, {
        "md_decision_recorded": "MD_run", "md_completion_log": str(csv),
    }, registry=None)
    assert "replicate_consistency" in result["warned"]
    assert "replicate_consistency" not in result["failed"]


def test_claim_boundary_fails_without_declaration():
    result = G7_evaluator.evaluate(_CONTRACT, {"md_decision_recorded": "MD_run"}, registry=None)
    assert "claim_boundary" in result["failed"]
    assert "claim_boundary" in result["critical_failures"]


def test_claim_boundary_passes_with_explicit_declaration():
    result = G7_evaluator.evaluate(_CONTRACT, {
        "md_decision_recorded": "MD_run",
        "required_claim_boundary": "MD results are binding-stability annotations only "
                                    "and do not constitute experimental validation.",
    }, registry=None)
    assert "claim_boundary" in result["passed"]


def test_md_not_used_for_ranking_vacuously_passes_with_no_selection():
    result = G7_evaluator.evaluate(_CONTRACT, {"md_decision_recorded": "MD_run"}, registry=None)
    assert "md_not_used_for_ranking" in result["passed"]


def test_md_not_used_for_ranking_fails_if_reranked():
    result = G7_evaluator.evaluate(_CONTRACT, {
        "md_decision_recorded": "MD_run",
        "selected_candidates_csv": "some/path.csv",
        "selection_reranked_by_md": True,
    }, registry=None)
    assert "md_not_used_for_ranking" in result["failed"]


def test_seed_recorded_fails_without_seeds(tmp_path):
    df = pd.DataFrame({"compound_id": ["A", "B"], "status": ["completed", "completed"],
                        "prod_ns": [50.0, 50.0]})
    csv = tmp_path / "traj.csv"
    df.to_csv(csv, index=False)
    result = G7_evaluator.evaluate(_CONTRACT, {
        "md_decision_recorded": "MD_run", "md_completion_log": str(csv),
    }, registry=None)
    assert "seed_recorded" in result["failed"]


def test_refuses_to_read_a_file_with_assay_columns(tmp_path):
    df = pd.DataFrame({"compound_id": ["A"], "status": ["completed"], "prod_ns": [50.0],
                        "known_active_cAMP": [1]})
    csv = tmp_path / "tainted.csv"
    df.to_csv(csv, index=False)
    result = G7_evaluator.evaluate(_CONTRACT, {
        "md_decision_recorded": "MD_run", "md_completion_log": str(csv),
    }, registry=None)
    assert "trajectory_completion" in result["failed"]
    assert "forbidden columns" in result["evidence"]["trajectory_completion"]["error"]

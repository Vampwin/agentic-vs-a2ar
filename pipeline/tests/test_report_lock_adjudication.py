"""
Regression tests for the human-adjudication release path of ReportLock.

This path was previously unreachable for two independent reasons:

  1. HumanInterventionEntry carried no artifact_ids field, so an adjudication
     record could not name the artifact it covered.
  2. ReportLock.can_report() read rec["data"]["artifact_ids"], but
     ArtifactRegistry._append flattens the entry into the record's top level,
     so the lookup silently matched nothing.

Together these meant every human-adjudicated artifact stayed permanently
blocked. These tests pin the corrected behaviour, and equally importantly pin
that adjudication remains *scoped* — it must not become a blanket unlock.
"""

from __future__ import annotations

import pytest

from kdss_audit.artifact_registry import ArtifactRegistry, HumanInterventionEntry
from kdss_audit.report_lock import ReportLock


def _make_lock(tmp_path):
    registry = ArtifactRegistry(tmp_path / "audit_registry.jsonl")
    lock = ReportLock(registry=registry, lock_file=tmp_path / "report_lock.json")
    return registry, lock


def _intervention(intervention_id: str, artifact_ids: list[str]) -> HumanInterventionEntry:
    return HumanInterventionEntry(
        intervention_id=intervention_id,
        gate="G5b_performance_v1",
        information_visible_to_human="EF@1% below frozen threshold; AUC above threshold",
        decision="proceed",
        prespecified=True,
        effect_on_result="metric reported with shortfall disclosed",
        artifact_ids=artifact_ids,
    )


def test_adjudication_releases_named_artifact(tmp_path):
    """A human intervention naming an artifact must release exactly that artifact."""
    registry, lock = _make_lock(tmp_path)

    assert not lock.can_report("ranked_scores_csv"), (
        "Artifact must start blocked when it has no passed contract."
    )

    registry.record_human_intervention(_intervention("hi_001", ["ranked_scores_csv"]))

    assert lock.can_report("ranked_scores_csv"), (
        "A human adjudication record naming the artifact must release it. "
        "If this fails, the artifact_ids lookup in ReportLock.can_report() is "
        "reading the wrong nesting level again."
    )
    lock.assert_can_report(["ranked_scores_csv"])  # must not raise


def test_adjudication_is_scoped_to_named_artifacts(tmp_path):
    """Adjudicating one artifact must not unlock a different one."""
    registry, lock = _make_lock(tmp_path)
    registry.record_human_intervention(_intervention("hi_002", ["ranked_scores_csv"]))

    assert not lock.can_report("selection_panel_csv"), (
        "Adjudication must be scoped to the artifacts it names; it is not a "
        "blanket unlock."
    )
    with pytest.raises(PermissionError):
        lock.assert_can_report(["ranked_scores_csv", "selection_panel_csv"])


def test_adjudication_without_artifact_ids_releases_nothing(tmp_path):
    """An intervention that names no artifact must not release anything."""
    registry, lock = _make_lock(tmp_path)
    registry.record_human_intervention(_intervention("hi_003", []))

    assert not lock.can_report("ranked_scores_csv"), (
        "An adjudication record with an empty artifact_ids list must not "
        "release any artifact."
    )


def test_artifact_ids_persist_at_record_top_level(tmp_path):
    """
    Pins the storage contract the release path depends on: _append flattens the
    entry, so artifact_ids lives at the record's top level, not under "data".
    """
    registry, _ = _make_lock(tmp_path)
    registry.record_human_intervention(_intervention("hi_004", ["artifact_x"]))

    records = [
        r for r in registry.get_all_records()
        if r["record_type"] == "human_intervention"
    ]
    assert len(records) == 1
    rec = records[0]
    assert rec["artifact_ids"] == ["artifact_x"], (
        "artifact_ids must be readable at the record top level."
    )
    assert "data" not in rec, (
        "Records are flattened by _append; a nested 'data' key would mean the "
        "storage format changed and ReportLock.can_report() must be updated too."
    )

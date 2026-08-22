"""
Tests for fault_injection/label_assigner.py.

Verifies Gate 4 label criteria: unsafe_release = 1 iff fault AND escaped
reference condition; 0 for clean cases or caught faults.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from fault_injection.label_assigner import assign_labels, summarise_labels


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_outcome(
    instance_id: str = "F01-001",
    family_id: str = "F01",
    condition: str = "B",
    case_type: str = "fault",
    escaped: bool = True,
    error: bool = False,
) -> dict:
    return {
        "instance_id": instance_id,
        "family_id": family_id,
        "condition": condition,
        "case_type": case_type,
        "escaped": escaped,
        "error": error,
    }


def _make_clean(condition: str = "B", idx: int = 0) -> dict:
    return _make_outcome(
        instance_id=f"CLEAN-{idx:03d}",
        family_id="",
        condition=condition,
        case_type="clean",
        escaped=False,
    )


# ---------------------------------------------------------------------------
# TestLabelCriteria
# ---------------------------------------------------------------------------

class TestLabelCriteria:
    """Gate 4 label logic: fault+escaped=1, everything else=0."""

    def test_fault_escaped_is_positive(self):
        outcomes = [_make_outcome(escaped=True)]
        df = assign_labels(outcomes)
        assert df.iloc[0]["unsafe_release"] == 1

    def test_fault_caught_is_negative(self):
        outcomes = [_make_outcome(escaped=False)]
        df = assign_labels(outcomes)
        assert df.iloc[0]["unsafe_release"] == 0

    def test_clean_case_is_negative(self):
        outcomes = [_make_clean()]
        df = assign_labels(outcomes)
        assert df.iloc[0]["unsafe_release"] == 0

    def test_clean_escaped_is_still_negative(self):
        # A clean case should never be labelled 1 even if escaped=True (anomaly in data)
        outcomes = [_make_outcome(case_type="clean", escaped=True)]
        df = assign_labels(outcomes)
        assert df.iloc[0]["unsafe_release"] == 0

    def test_mixed_cases(self):
        outcomes = [
            _make_outcome("F01-001", "F01", "B", "fault", escaped=True),
            _make_outcome("F01-002", "F01", "B", "fault", escaped=False),
            _make_clean("B", 0),
        ]
        df = assign_labels(outcomes)
        assert len(df) == 3
        pos = df[df["unsafe_release"] == 1]
        assert len(pos) == 1
        assert pos.iloc[0]["instance_id"] == "F01-001"


# ---------------------------------------------------------------------------
# TestErrorFiltering
# ---------------------------------------------------------------------------

class TestErrorFiltering:
    """Outcomes with error=True must be excluded (fail-closed)."""

    def test_error_outcome_excluded(self):
        outcomes = [
            _make_outcome("F02-001", "F02", "B", "fault", escaped=True, error=True),
            _make_outcome("F01-001", "F01", "B", "fault", escaped=True, error=False),
        ]
        df = assign_labels(outcomes)
        assert len(df) == 1
        assert df.iloc[0]["instance_id"] == "F01-001"

    def test_all_errors_raises(self):
        outcomes = [
            _make_outcome(error=True),
        ]
        with pytest.raises(ValueError, match="no non-errored outcomes"):
            assign_labels(outcomes)


# ---------------------------------------------------------------------------
# TestReferenceCondition
# ---------------------------------------------------------------------------

class TestReferenceCondition:
    """Only outcomes from the reference condition are included."""

    def test_default_condition_is_B(self):
        outcomes = [
            _make_outcome("F01-001", "F01", "B", "fault", True),
            _make_outcome("F01-002", "F01", "D", "fault", True),
        ]
        df = assign_labels(outcomes)
        assert len(df) == 1
        assert all(df["condition"] == "B")

    def test_custom_reference_condition(self):
        outcomes = [
            _make_outcome("F01-001", "F01", "B", "fault", True),
            _make_outcome("F01-002", "F01", "D", "fault", True),
        ]
        df = assign_labels(outcomes, reference_condition="D")
        assert len(df) == 1
        assert all(df["condition"] == "D")

    def test_missing_reference_condition_raises(self):
        outcomes = [_make_outcome(condition="B")]
        with pytest.raises(ValueError, match="no non-errored outcomes"):
            assign_labels(outcomes, reference_condition="A")


# ---------------------------------------------------------------------------
# TestOutputSchema
# ---------------------------------------------------------------------------

class TestOutputSchema:
    """Return type and column names must match the specification."""

    def setup_method(self):
        outcomes = [
            _make_outcome("F01-001", "F01", "B", "fault", True),
            _make_clean("B", 0),
        ]
        self.df = assign_labels(outcomes)

    def test_returns_dataframe(self):
        assert isinstance(self.df, pd.DataFrame)

    def test_required_columns_present(self):
        required = {"instance_id", "family_id", "condition", "case_type",
                    "unsafe_release", "group"}
        assert required.issubset(set(self.df.columns))

    def test_unsafe_release_is_integer(self):
        assert self.df["unsafe_release"].dtype in (int, "int64", "int32")

    def test_group_equals_family_id(self):
        for _, row in self.df.iterrows():
            assert row["group"] == row["family_id"]

    def test_no_nans_in_unsafe_release(self):
        assert not self.df["unsafe_release"].isna().any()


# ---------------------------------------------------------------------------
# TestLeakageGuard
# ---------------------------------------------------------------------------

class TestLeakageGuard:
    """Forbidden terms in outcomes raise ValueError."""

    def test_camp_in_key_raises(self):
        outcomes = [{"camp_value": 5.0, "condition": "B", "case_type": "clean",
                     "escaped": False, "error": False, "instance_id": "X",
                     "family_id": ""}]
        with pytest.raises(ValueError, match="assay leakage"):
            assign_labels(outcomes)

    def test_held_out_label_id_raises(self):
        outcomes = [{"h01_label": 1, "condition": "B", "case_type": "clean",
                     "escaped": False, "error": False, "instance_id": "X",
                     "family_id": ""}]
        with pytest.raises(ValueError, match="assay leakage"):
            assign_labels(outcomes)

    def test_assay_in_value_raises(self):
        outcomes = [{"some_field": "assay_data", "condition": "B", "case_type": "clean",
                     "escaped": False, "error": False, "instance_id": "X",
                     "family_id": ""}]
        with pytest.raises(ValueError, match="assay leakage"):
            assign_labels(outcomes)

    def test_clean_outcome_does_not_raise(self):
        outcomes = [_make_clean("B", 0)]
        df = assign_labels(outcomes)
        assert len(df) == 1


# ---------------------------------------------------------------------------
# TestSummariseLabels
# ---------------------------------------------------------------------------

class TestSummariseLabels:
    """summarise_labels() returns aggregate counts."""

    def test_returns_dict(self):
        df = pd.DataFrame([
            {"instance_id": "A", "family_id": "F01", "condition": "B",
             "case_type": "fault", "unsafe_release": 1, "group": "F01"},
            {"instance_id": "B", "family_id": "F01", "condition": "B",
             "case_type": "clean", "unsafe_release": 0, "group": ""},
        ])
        summary = summarise_labels(df)
        assert isinstance(summary, dict)

    def test_counts_correct(self):
        df = pd.DataFrame([
            {"unsafe_release": 1, "case_type": "fault", "group": "F01"},
            {"unsafe_release": 1, "case_type": "fault", "group": "F01"},
            {"unsafe_release": 0, "case_type": "clean", "group": ""},
        ])
        summary = summarise_labels(df)
        assert summary["n_total"] == 3
        assert summary["n_positive"] == 2
        assert summary["n_negative"] == 1
        assert summary["n_fault_cases"] == 2
        assert summary["n_clean_cases"] == 1

    def test_positive_rate(self):
        df = pd.DataFrame([
            {"unsafe_release": 1, "case_type": "fault", "group": "F01"},
            {"unsafe_release": 0, "case_type": "clean", "group": ""},
            {"unsafe_release": 0, "case_type": "fault", "group": "F02"},
            {"unsafe_release": 0, "case_type": "fault", "group": "F02"},
        ])
        summary = summarise_labels(df)
        assert abs(summary["positive_rate"] - 0.25) < 1e-6

    def test_n_groups(self):
        df = pd.DataFrame([
            {"unsafe_release": 1, "case_type": "fault", "group": "F01"},
            {"unsafe_release": 1, "case_type": "fault", "group": "F02"},
            {"unsafe_release": 0, "case_type": "clean", "group": ""},
        ])
        summary = summarise_labels(df)
        assert summary["n_groups"] == 3  # F01, F02, ""

    def test_empty_dataframe(self):
        df = pd.DataFrame(columns=["unsafe_release", "case_type", "group"])
        summary = summarise_labels(df)
        assert summary["n_total"] == 0
        assert summary["positive_rate"] == 0.0


# ---------------------------------------------------------------------------
# TestIntegrationWithPilotOutput
# ---------------------------------------------------------------------------

class TestIntegrationWithPilotOutput:
    """Round-trip test using pilot-style outcome list."""

    def _pilot_outcomes(self) -> list[dict]:
        outcomes = []
        for fid in ["F01", "F02", "F03"]:
            for i in range(3):
                escaped = fid != "F01"
                outcomes.append({
                    "instance_id": f"{fid}-{i:03d}",
                    "family_id": fid,
                    "condition": "B",
                    "case_type": "fault",
                    "escaped": escaped,
                    "error": False,
                })
        for i in range(2):
            outcomes.append({
                "instance_id": f"CLEAN-{i:03d}",
                "family_id": "",
                "condition": "B",
                "case_type": "clean",
                "escaped": False,
                "error": False,
            })
        return outcomes

    def test_correct_label_count(self):
        df = assign_labels(self._pilot_outcomes())
        assert len(df) == 11  # 9 fault + 2 clean

    def test_f01_all_caught(self):
        df = assign_labels(self._pilot_outcomes())
        f01 = df[df["family_id"] == "F01"]
        assert all(f01["unsafe_release"] == 0)

    def test_f02_f03_all_escaped(self):
        df = assign_labels(self._pilot_outcomes())
        f23 = df[df["family_id"].isin(["F02", "F03"])]
        assert all(f23["unsafe_release"] == 1)

    def test_clean_cases_all_negative(self):
        df = assign_labels(self._pilot_outcomes())
        clean = df[df["case_type"] == "clean"]
        assert all(clean["unsafe_release"] == 0)

    def test_summary_round_trip(self):
        df = assign_labels(self._pilot_outcomes())
        s = summarise_labels(df)
        assert s["n_positive"] == 6   # F02 + F03: 3 each
        assert s["n_negative"] == 5   # F01 (3) + clean (2)

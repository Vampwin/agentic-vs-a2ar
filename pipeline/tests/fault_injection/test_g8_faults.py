"""
Tests for G8 fault injectors: F17 (STALE_FIGURE), F18 (VALUE_MISMATCH).
"""

import csv
import os
import sys
import time
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from fault_injection.g8.f17_stale_figure import F17StaleFigure
from fault_injection.g8.f18_value_mismatch import F18ValueMismatch


# ---------------------------------------------------------------------------
# F17 STALE_FIGURE
# ---------------------------------------------------------------------------

class TestF17Spec:
    def setup_method(self):
        self.inj = F17StaleFigure()

    def test_family_id(self):
        assert self.inj.spec.family_id == "F17"

    def test_gate(self):
        assert self.inj.spec.gate == "G8"

    def test_severity(self):
        assert self.inj.spec.severity == "medium"

    def test_param_grid_length(self):
        assert len(self.inj.spec.param_grid) == 10

    def test_staleness_hours_ascending(self):
        hours = [p["staleness_hours"] for p in self.inj.spec.param_grid]
        assert hours == sorted(hours)

    def test_no_assay_leakage(self):
        self.inj.verify_no_assay_leakage()


class TestF17Inject:
    def setup_method(self):
        self.inj = F17StaleFigure()

    def test_figure_content_changed(self, tmp_path):
        fig = tmp_path / "main_figure.png"
        original = b"\x89PNG\r\n\x1a\nORIGINAL_CONTENT_BYTES"
        fig.write_bytes(original)

        params = {"staleness_hours": 24, "target_figure": "main_figure"}
        rec = self.inj.inject(tmp_path, params, seed=0, instance_index=3)

        new_content = fig.read_bytes()
        assert new_content != original
        assert rec.family_id == "F17"

    def test_mtime_set_in_past(self, tmp_path):
        fig = tmp_path / "main_figure.png"
        fig.write_bytes(b"\x89PNG\r\n\x1a\nSOME_CONTENT")
        before = time.time()

        params = {"staleness_hours": 48, "target_figure": "main_figure"}
        self.inj.inject(tmp_path, params, seed=0, instance_index=4)

        mtime = os.stat(fig).st_mtime
        expected_oldest = before - 48 * 3600 - 5  # 5 s tolerance
        assert mtime < before
        assert mtime > expected_oldest

    def test_creates_figure_if_missing(self, tmp_path):
        # No figure present — injector should create a synthetic one
        params = {"staleness_hours": 1, "target_figure": "main_figure"}
        rec = self.inj.inject(tmp_path, params, seed=0, instance_index=0)
        assert len(rec.affected_items) == 1
        assert (tmp_path / "main_figure.png").exists()

    def test_record_ground_truth_mentions_staleness(self, tmp_path):
        params = {"staleness_hours": 72, "target_figure": "main_figure"}
        rec = self.inj.inject(tmp_path, params, seed=0, instance_index=0)
        assert "72" in rec.ground_truth or "stale" in rec.ground_truth.lower()

    def test_all_staleness_values_covered(self):
        hours = [p["staleness_hours"] for p in self.inj.spec.param_grid]
        assert 1 in hours
        assert 2160 in hours   # 90 days


# ---------------------------------------------------------------------------
# F18 VALUE_MISMATCH
# ---------------------------------------------------------------------------

class TestF18Spec:
    def setup_method(self):
        self.inj = F18ValueMismatch()

    def test_family_id(self):
        assert self.inj.spec.family_id == "F18"

    def test_gate(self):
        assert self.inj.spec.gate == "G8"

    def test_severity(self):
        assert self.inj.spec.severity == "medium"

    def test_param_grid_not_empty(self):
        assert len(self.inj.spec.param_grid) >= 1

    def test_no_assay_leakage(self):
        self.inj.verify_no_assay_leakage()

    def test_sign_values(self):
        signs = {p["sign"] for p in self.inj.spec.param_grid}
        assert all(abs(s) == 1 for s in signs)


def _make_results_table(tmp_path: Path) -> Path:
    """Write F18's expected results_table.csv schema."""
    p = tmp_path / "results_table.csv"
    p.write_text(
        "condition,seer_pct,brier_score,n_compounds\n"
        "A_hybrid,0.15,0.12,172\n"
        "B_contracts_only,0.35,0.28,172\n"
        "C_monitor_only,0.42,0.31,172\n"
        "D_neither,0.58,0.45,172\n",
        encoding="utf-8",
    )
    return p


class TestF18Inject:
    def setup_method(self):
        self.inj = F18ValueMismatch()

    def test_value_modified(self, tmp_path):
        original = _make_results_table(tmp_path).read_text()
        params = self.inj.spec.param_grid[0]
        rec = self.inj.inject(tmp_path, params, seed=0, instance_index=0)
        assert rec.family_id == "F18"
        modified = (tmp_path / "results_table.csv").read_text()
        assert original != modified

    def test_returns_injection_record(self, tmp_path):
        _make_results_table(tmp_path)
        params = self.inj.spec.param_grid[0]
        rec = self.inj.inject(tmp_path, params, seed=0, instance_index=0)
        assert isinstance(rec.affected_items, list)
        assert len(rec.affected_items) >= 1

    def test_creates_file_if_missing(self, tmp_path):
        params = self.inj.spec.param_grid[0]
        rec = self.inj.inject(tmp_path, params, seed=0, instance_index=0)
        assert (tmp_path / "results_table.csv").exists()
        assert rec.family_id == "F18"

    def test_different_deltas_produce_different_values(self, tmp_path):
        import shutil
        ws1 = tmp_path / "ws1"; ws2 = tmp_path / "ws2"
        ws1.mkdir(); ws2.mkdir()
        _make_results_table(ws1); _make_results_table(ws2)

        params_small = {"delta_fraction": 0.01, "sign": 1, "target_column": "seer_pct"}
        params_large = {"delta_fraction": 0.50, "sign": 1, "target_column": "seer_pct"}
        self.inj.inject(ws1, params_small, seed=0, instance_index=0)
        self.inj.inject(ws2, params_large, seed=0, instance_index=0)

        assert (ws1 / "results_table.csv").read_text() != (ws2 / "results_table.csv").read_text()

    def test_no_camp_in_record(self, tmp_path):
        _make_results_table(tmp_path)
        params = self.inj.spec.param_grid[0]
        rec = self.inj.inject(tmp_path, params, seed=0, instance_index=0)
        assert "camp" not in rec.ground_truth.lower()
        assert "assay" not in rec.ground_truth.lower()

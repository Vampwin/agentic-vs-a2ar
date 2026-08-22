"""
Tests for G4 fault injectors: F09 (SCORE_DIRECTION), F10 (MW_SCORE_BIAS).
"""

import csv
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from fault_injection.g4.f09_score_direction import F09ScoreDirection
from fault_injection.g4.f10_mw_score_bias import F10MwScoreBias


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_scores_csv(tmp_path: Path, n: int = 50, col: str = "vina_score") -> Path:
    p = tmp_path / "docking_scores.csv"
    rows = [f"compound_id,SMILES,{col}"]
    for i in range(n):
        score = -5.0 - i * 0.1
        rows.append(f"CMP{i:05d},CCO,{score:.2f}")
    p.write_text("\n".join(rows), encoding="utf-8")
    return p


# ---------------------------------------------------------------------------
# F09 SCORE_DIRECTION
# ---------------------------------------------------------------------------

class TestF09Spec:
    def setup_method(self):
        self.inj = F09ScoreDirection()

    def test_family_id(self):
        assert self.inj.spec.family_id == "F09"

    def test_gate(self):
        assert self.inj.spec.gate == "G4"

    def test_severity(self):
        assert self.inj.spec.severity == "critical"

    def test_param_grid_length(self):
        assert len(self.inj.spec.param_grid) == 15

    def test_no_assay_leakage(self):
        self.inj.verify_no_assay_leakage()

    def test_score_columns_in_params(self):
        cols = {p["score_column"] for p in self.inj.spec.param_grid}
        assert "vina_score" in cols

    def test_instance_id_format(self):
        assert self.inj.instance_id(14) == "F09_015"


class TestF09Inject:
    def setup_method(self):
        self.inj = F09ScoreDirection()

    def test_scores_negated(self, tmp_path):
        _make_scores_csv(tmp_path, n=50, col="vina_score")
        params = {"n_affected": None, "score_column": "vina_score"}
        rec = self.inj.inject(tmp_path, params, seed=0, instance_index=0)

        with open(tmp_path / "docking_scores.csv", newline="") as f:
            rows = list(csv.DictReader(f))
        scores = [float(r["vina_score"]) for r in rows]

        # Original scores were all negative; after negation all should be positive
        assert all(s > 0 for s in scores)
        assert rec.family_id == "F09"

    def test_partial_negation(self, tmp_path):
        _make_scores_csv(tmp_path, n=50, col="vina_score")
        params = {"n_affected": 10, "score_column": "vina_score"}
        rec = self.inj.inject(tmp_path, params, seed=42, instance_index=3)

        with open(tmp_path / "docking_scores.csv", newline="") as f:
            rows = list(csv.DictReader(f))
        scores = [float(r["vina_score"]) for r in rows]

        positive_count = sum(1 for s in scores if s > 0)
        negative_count = sum(1 for s in scores if s < 0)
        assert positive_count == 10
        assert negative_count == 40

    def test_affected_ids_in_record(self, tmp_path):
        _make_scores_csv(tmp_path, n=30, col="vina_score")
        params = {"n_affected": 5, "score_column": "vina_score"}
        rec = self.inj.inject(tmp_path, params, seed=1, instance_index=0)
        assert len(rec.affected_items) == 5

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            self.inj.inject(tmp_path, {"n_affected": None, "score_column": "vina_score"},
                            seed=0, instance_index=0)

    def test_no_matching_col_raises(self, tmp_path):
        _make_scores_csv(tmp_path, n=10, col="binding_energy")
        # Remove the default score column name
        params = {"n_affected": None, "score_column": "nonexistent_xyz_score"}
        # Flexible match falls back to any "score" column — should still work
        # because "binding_energy" does not contain "score" either
        with pytest.raises(ValueError, match="no score column"):
            # Overwrite to a CSV with no "score" keyword at all
            (tmp_path / "docking_scores.csv").write_text("compound_id,energy\nA,-5.0\n")
            self.inj.inject(tmp_path, params, seed=0, instance_index=0)


# ---------------------------------------------------------------------------
# F10 MW_SCORE_BIAS
# ---------------------------------------------------------------------------

class TestF10Spec:
    def setup_method(self):
        self.inj = F10MwScoreBias()

    def test_family_id(self):
        assert self.inj.spec.family_id == "F10"

    def test_gate(self):
        assert self.inj.spec.gate == "G4"

    def test_severity(self):
        assert self.inj.spec.severity == "high"

    def test_param_grid_length(self):
        assert len(self.inj.spec.param_grid) == 12

    def test_no_assay_leakage(self):
        self.inj.verify_no_assay_leakage()

    def test_mw_direction_values(self):
        directions = {p["mw_direction"] for p in self.inj.spec.param_grid}
        assert "heavy_favoured" in directions or "high_mw" in directions or len(directions) >= 1


class TestF10Inject:
    def setup_method(self):
        self.inj = F10MwScoreBias()

    def _make_scores_with_smiles(self, tmp_path: Path, n: int = 30) -> Path:
        smiles_list = [
            "CCO", "CCCO", "CCCCO", "CCCCCO",
            "c1ccccc1", "c1ccccc1CC", "c1ccccc1CCCC",
            "CC(=O)O", "CC(=O)OCC", "CC(N)=O",
        ]
        p = tmp_path / "docking_scores.csv"
        rows = ["compound_id,SMILES,vina_score"]
        for i in range(n):
            smi = smiles_list[i % len(smiles_list)]
            rows.append(f"CMP{i:05d},{smi},{-5.0 - i * 0.1:.2f}")
        p.write_text("\n".join(rows), encoding="utf-8")
        return p

    def test_returns_injection_record(self, tmp_path):
        self._make_scores_with_smiles(tmp_path)
        params = self.inj.spec.param_grid[0]
        rec = self.inj.inject(tmp_path, params, seed=0, instance_index=0)
        assert rec.family_id == "F10"
        assert rec.gate == "G4"

    def test_scores_modified(self, tmp_path):
        self._make_scores_with_smiles(tmp_path)
        original_content = (tmp_path / "docking_scores.csv").read_text()
        params = self.inj.spec.param_grid[0]
        self.inj.inject(tmp_path, params, seed=0, instance_index=0)
        new_content = (tmp_path / "docking_scores.csv").read_text()
        assert original_content != new_content

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            self.inj.inject(tmp_path, self.inj.spec.param_grid[0], seed=0, instance_index=0)

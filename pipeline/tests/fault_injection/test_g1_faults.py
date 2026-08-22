"""
Tests for G1 fault injectors: F03 (STEREO_LOSS), F04 (SALT_RETAINED),
F05 (INVALID_VALENCE).
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from fault_injection.g1.f03_stereo_loss import F03StereoLoss
from fault_injection.g1.f04_salt_retained import F04SaltRetained
from fault_injection.g1.f05_invalid_valence import F05InvalidValence


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_ledger(tmp_path: Path, smiles_rows: list[str]) -> Path:
    p = tmp_path / "standardised_ledger.csv"
    header = "compound_id,SMILES,name"
    rows = [header] + smiles_rows
    p.write_text("\n".join(rows), encoding="utf-8")
    return p


def _stereo_smiles(n: int = 20) -> list[str]:
    return [f"CMP{i:05d},[C@@H](O)(F)Cl,stereo_{i}" for i in range(n)]


def _simple_smiles(n: int = 20) -> list[str]:
    return [f"CMP{i:05d},CCO,cmpd_{i}" for i in range(n)]


# ---------------------------------------------------------------------------
# F03 STEREO_LOSS
# ---------------------------------------------------------------------------

class TestF03Spec:
    def setup_method(self):
        self.inj = F03StereoLoss()

    def test_family_id(self):
        assert self.inj.spec.family_id == "F03"

    def test_gate(self):
        assert self.inj.spec.gate == "G1"

    def test_severity(self):
        assert self.inj.spec.severity == "high"

    def test_param_grid_length(self):
        assert len(self.inj.spec.param_grid) == 12

    def test_stereo_fraction_range(self):
        fracs = [p["stereo_fraction"] for p in self.inj.spec.param_grid]
        assert min(fracs) >= 0.0
        assert max(fracs) <= 1.0

    def test_no_assay_leakage(self):
        self.inj.verify_no_assay_leakage()


class TestF03Inject:
    def setup_method(self):
        self.inj = F03StereoLoss()

    def test_stereocentres_removed(self, tmp_path):
        _write_ledger(tmp_path, _stereo_smiles(20))
        params = {"stereo_fraction": 1.0}
        rec = self.inj.inject(tmp_path, params, seed=0, instance_index=0)
        content = (tmp_path / "standardised_ledger.csv").read_text()
        assert "@" not in content
        assert rec.family_id == "F03"

    def test_partial_stereo_removal(self, tmp_path):
        _write_ledger(tmp_path, _stereo_smiles(20))
        params = {"stereo_fraction": 0.5}
        rec = self.inj.inject(tmp_path, params, seed=7, instance_index=5)
        content = (tmp_path / "standardised_ledger.csv").read_text()
        # Some rows should still have @ (50% kept)
        assert "@" in content
        assert len(rec.affected_items) > 0

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            self.inj.inject(tmp_path, {"stereo_fraction": 0.5}, seed=0, instance_index=0)

    def test_containment_criterion_mentions_g1(self, tmp_path):
        _write_ledger(tmp_path, _stereo_smiles(10))
        rec = self.inj.inject(tmp_path, {"stereo_fraction": 1.0}, seed=0, instance_index=0)
        assert "G1" in rec.containment_criterion or "stereo" in rec.containment_criterion.lower()


# ---------------------------------------------------------------------------
# F04 SALT_RETAINED
# ---------------------------------------------------------------------------

class TestF04Spec:
    def setup_method(self):
        self.inj = F04SaltRetained()

    def test_family_id(self):
        assert self.inj.spec.family_id == "F04"

    def test_gate(self):
        assert self.inj.spec.gate == "G1"

    def test_severity(self):
        assert self.inj.spec.severity == "medium"

    def test_param_grid_length(self):
        assert len(self.inj.spec.param_grid) == 12

    def test_no_assay_leakage(self):
        self.inj.verify_no_assay_leakage()

    def test_salt_types_in_param_grid(self):
        salt_types = {p["salt_type"] for p in self.inj.spec.param_grid}
        assert "Cl" in salt_types or "Na" in salt_types


class TestF04Inject:
    def setup_method(self):
        self.inj = F04SaltRetained()

    def test_salt_appended(self, tmp_path):
        _write_ledger(tmp_path, _simple_smiles(20))
        params = {"salt_type": "Cl", "n_compounds": 5}
        rec = self.inj.inject(tmp_path, params, seed=0, instance_index=0)
        content = (tmp_path / "standardised_ledger.csv").read_text()
        assert "[Cl-]" in content
        assert rec.family_id == "F04"
        assert len(rec.affected_items) == 5

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            self.inj.inject(tmp_path, {"salt_type": "Na", "n_compounds": 3},
                            seed=0, instance_index=0)

    def test_dot_separator_present(self, tmp_path):
        _write_ledger(tmp_path, _simple_smiles(10))
        params = {"salt_type": "Na", "n_compounds": 3}
        self.inj.inject(tmp_path, params, seed=1, instance_index=1)
        content = (tmp_path / "standardised_ledger.csv").read_text()
        # Salt appended with . separator: SMILES.[Na+]
        assert ".[Na+]" in content


# ---------------------------------------------------------------------------
# F05 INVALID_VALENCE
# ---------------------------------------------------------------------------

class TestF05Spec:
    def setup_method(self):
        self.inj = F05InvalidValence()

    def test_family_id(self):
        assert self.inj.spec.family_id == "F05"

    def test_gate(self):
        assert self.inj.spec.gate == "G1"

    def test_severity(self):
        assert self.inj.spec.severity == "high"

    def test_param_grid_length(self):
        assert len(self.inj.spec.param_grid) == 10

    def test_no_assay_leakage(self):
        self.inj.verify_no_assay_leakage()


class TestF05Inject:
    def setup_method(self):
        self.inj = F05InvalidValence()

    def test_invalid_smiles_written(self, tmp_path):
        _write_ledger(tmp_path, _simple_smiles(20))
        params = self.inj.spec.param_grid[0]
        rec = self.inj.inject(tmp_path, params, seed=0, instance_index=0)
        content = (tmp_path / "standardised_ledger.csv").read_text()
        # Should contain some indication of broken structure
        assert rec.family_id == "F05"
        assert len(rec.affected_items) >= 1

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            self.inj.inject(tmp_path, self.inj.spec.param_grid[0], seed=0, instance_index=0)

    def test_record_ground_truth_not_empty(self, tmp_path):
        _write_ledger(tmp_path, _simple_smiles(20))
        rec = self.inj.inject(tmp_path, self.inj.spec.param_grid[0], seed=0, instance_index=0)
        assert len(rec.ground_truth) > 10

"""
Tests for G0 fault injectors: F01 (ID_DUPLICATE), F02 (HASH_MISMATCH).
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from fault_injection.g0.f01_id_duplicate import F01IdDuplicate
from fault_injection.g0.f02_hash_mismatch import F02HashMismatch


# ---------------------------------------------------------------------------
# F01 ID_DUPLICATE
# ---------------------------------------------------------------------------

class TestF01Spec:
    def setup_method(self):
        self.inj = F01IdDuplicate()

    def test_family_id(self):
        assert self.inj.spec.family_id == "F01"

    def test_gate(self):
        assert self.inj.spec.gate == "G0"

    def test_severity(self):
        assert self.inj.spec.severity == "medium"

    def test_param_grid_length(self):
        assert len(self.inj.spec.param_grid) == 10

    def test_param_fractions_ascending(self):
        fracs = [p["duplicate_fraction"] for p in self.inj.spec.param_grid]
        assert fracs == sorted(fracs)

    def test_no_assay_leakage(self):
        self.inj.verify_no_assay_leakage()  # must not raise

    def test_instance_id_format(self):
        assert self.inj.instance_id(0) == "F01_001"
        assert self.inj.instance_id(9) == "F01_010"


def _make_compounds_csv(tmp_path: Path, n: int = 50) -> Path:
    p = tmp_path / "compounds.csv"
    rows = ["compound_id,SMILES,name"] + [
        f"CMP{i:05d},CCO,compound_{i}" for i in range(n)
    ]
    p.write_text("\n".join(rows), encoding="utf-8")
    return p


class TestF01Inject:
    def setup_method(self):
        self.inj = F01IdDuplicate()

    def test_returns_injection_record(self, tmp_path):
        _make_compounds_csv(tmp_path)
        params = {"duplicate_fraction": 0.05}
        rec = self.inj.inject(tmp_path, params, seed=42, instance_index=4)
        assert rec.family_id == "F01"
        assert rec.gate == "G0"
        assert isinstance(rec.affected_items, list)
        assert len(rec.affected_items) >= 1

    def test_file_grows_by_duplicated_rows(self, tmp_path):
        _make_compounds_csv(tmp_path, n=100)
        import csv as csvmod
        with open(tmp_path / "compounds.csv", newline="") as f:
            original_count = sum(1 for _ in csvmod.DictReader(f))

        params = {"duplicate_fraction": 0.10}
        rec = self.inj.inject(tmp_path, params, seed=0, instance_index=9)

        with open(tmp_path / "compounds.csv", newline="") as f:
            new_count = sum(1 for _ in csvmod.DictReader(f))

        assert new_count > original_count
        assert len(rec.affected_items) == new_count - original_count

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            self.inj.inject(tmp_path, {"duplicate_fraction": 0.01}, seed=0, instance_index=0)

    def test_different_seeds_give_different_items(self, tmp_path):
        _make_compounds_csv(tmp_path, n=100)
        params = {"duplicate_fraction": 0.05}

        import shutil
        ws1 = tmp_path / "ws1"
        ws2 = tmp_path / "ws2"
        ws1.mkdir(); ws2.mkdir()
        shutil.copy(tmp_path / "compounds.csv", ws1 / "compounds.csv")
        shutil.copy(tmp_path / "compounds.csv", ws2 / "compounds.csv")

        rec1 = self.inj.inject(ws1, params, seed=1, instance_index=0)
        rec2 = self.inj.inject(ws2, params, seed=2, instance_index=0)
        # Different seeds → different affected items (probabilistically)
        assert rec1.affected_items != rec2.affected_items


# ---------------------------------------------------------------------------
# F02 HASH_MISMATCH
# ---------------------------------------------------------------------------

class TestF02Spec:
    def setup_method(self):
        self.inj = F02HashMismatch()

    def test_family_id(self):
        assert self.inj.spec.family_id == "F02"

    def test_gate(self):
        assert self.inj.spec.gate == "G0"

    def test_severity(self):
        assert self.inj.spec.severity == "high"

    def test_param_grid_length(self):
        assert len(self.inj.spec.param_grid) == 15

    def test_no_assay_leakage(self):
        self.inj.verify_no_assay_leakage()

    def test_instance_id_format(self):
        assert self.inj.instance_id(0) == "F02_001"


class TestF02Inject:
    def setup_method(self):
        self.inj = F02HashMismatch()

    def _make_workspace(self, tmp_path: Path) -> Path:
        (tmp_path / "registered_input.bin").write_bytes(b"A" * 256)
        return tmp_path

    def test_returns_injection_record(self, tmp_path):
        self._make_workspace(tmp_path)
        params = {"offset_fraction": 0.1, "flip_bit": 0}
        rec = self.inj.inject(tmp_path, params, seed=0, instance_index=0)
        assert rec.family_id == "F02"
        assert rec.gate == "G0"

    def test_file_is_modified(self, tmp_path):
        self._make_workspace(tmp_path)
        original = (tmp_path / "registered_input.bin").read_bytes()
        params = {"offset_fraction": 0.5, "flip_bit": 3}
        self.inj.inject(tmp_path, params, seed=0, instance_index=7)
        modified = (tmp_path / "registered_input.bin").read_bytes()
        assert original != modified
        assert len(original) == len(modified)

    def test_exactly_one_byte_differs(self, tmp_path):
        self._make_workspace(tmp_path)
        original = (tmp_path / "registered_input.bin").read_bytes()
        params = {"offset_fraction": 0.2, "flip_bit": 7}
        self.inj.inject(tmp_path, params, seed=0, instance_index=3)
        modified = (tmp_path / "registered_input.bin").read_bytes()
        differences = sum(1 for a, b in zip(original, modified) if a != b)
        assert differences == 1

    def test_missing_workspace_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            params = {"offset_fraction": 0.1, "flip_bit": 0}
            self.inj.inject(tmp_path, params, seed=0, instance_index=0)

    def test_empty_file_raises(self, tmp_path):
        (tmp_path / "registered_input.bin").write_bytes(b"")
        params = {"offset_fraction": 0.5, "flip_bit": 0}
        with pytest.raises(ValueError, match="empty"):
            self.inj.inject(tmp_path, params, seed=0, instance_index=0)

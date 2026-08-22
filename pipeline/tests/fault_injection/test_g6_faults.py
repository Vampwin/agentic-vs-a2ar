"""
Tests for G6 fault injectors:
  F13 (UNDECLARED_EXTERNAL_LABEL_CONSTRAINT),
  F14 (SCAFFOLD_DUPLICATION).
"""

import csv
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from fault_injection.g6.f13_undeclared_external_label_constraint import (
    F13UndeclaredExternalLabelConstraint,
)
from fault_injection.g6.f14_scaffold_duplication import F14ScaffoldDuplication


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_RANKED_SMILES = [
    "CCO", "CCCO", "CCCCO", "CCCCCO",
    "c1ccccc1", "c1ccccc1CC", "c1ccccc1CCCC",
    "CC(=O)O", "CC(=O)OCC", "CC(N)=O",
    "CNCCO", "CN(C)CCO", "CCNCC", "c1ccncc1",
    "c1ccc(O)cc1", "c1ccc(N)cc1", "c1ccc(F)cc1", "c1ccc(Cl)cc1",
    "CC1=CC=CC=C1", "CC1=CC=CC=C1F",
]


def _make_ranked_csv(tmp_path: Path, n: int = 20) -> Path:
    p = tmp_path / "ranked_compounds.csv"
    rows = ["compound_id,SMILES,vina_score,MW,rotatable_bonds,HBA,HBD,ring_count"]
    for i in range(n):
        smi = _RANKED_SMILES[i % len(_RANKED_SMILES)]
        mw = 40 + i * 10
        rows.append(f"CMP{i:05d},{smi},{-5.0 - i * 0.1:.2f},{mw},{i%6},{i%5},{i%4},{i%3}")
    p.write_text("\n".join(rows), encoding="utf-8")
    return p


def _make_selected_panel(tmp_path: Path, n: int = 15) -> Path:
    p = tmp_path / "selected_panel.csv"
    rows = ["compound_id,SMILES,vina_score"]
    for i in range(n):
        smi = _RANKED_SMILES[i % len(_RANKED_SMILES)]
        rows.append(f"CMP{i:05d},{smi},{-5.0 - i * 0.1:.2f}")
    p.write_text("\n".join(rows), encoding="utf-8")
    return p


# ---------------------------------------------------------------------------
# F13 UNDECLARED_EXTERNAL_LABEL_CONSTRAINT
# ---------------------------------------------------------------------------

class TestF13Spec:
    def setup_method(self):
        self.inj = F13UndeclaredExternalLabelConstraint()

    def test_family_id(self):
        assert self.inj.spec.family_id == "F13"

    def test_gate(self):
        assert self.inj.spec.gate == "G6"

    def test_severity(self):
        assert self.inj.spec.severity == "critical"

    def test_param_grid_length(self):
        assert len(self.inj.spec.param_grid) == 20

    def test_no_assay_leakage(self):
        self.inj.verify_no_assay_leakage()

    def test_filter_names_structural_only(self):
        allowed = {"MW", "rotatable_bonds", "HBA", "HBD", "ring_count"}
        forbidden_assay = {"camp", "active", "assay", "ic50", "ec50", "ki", "activity"}
        for p in self.inj.spec.param_grid:
            fname = p["filter_name"].lower()
            for term in forbidden_assay:
                assert term not in fname, f"Assay term '{term}' found in filter_name '{fname}'"
            # The filter must be from the structural-only set
            assert p["filter_name"] in allowed, (
                f"Non-structural filter '{p['filter_name']}' in F13 params"
            )

    def test_instance_id_format(self):
        assert self.inj.instance_id(0) == "F13_001"
        assert self.inj.instance_id(19) == "F13_020"


class TestF13Inject:
    def setup_method(self):
        self.inj = F13UndeclaredExternalLabelConstraint()

    def test_rows_removed_by_filter(self, tmp_path):
        _make_ranked_csv(tmp_path, n=20)
        # Use MW filter with a very tight threshold to remove some rows
        params = {"filter_name": "MW", "threshold": 100.0}
        rec = self.inj.inject(tmp_path, params, seed=0, instance_index=0)

        with open(tmp_path / "ranked_compounds.csv", newline="") as f:
            remaining = list(csv.DictReader(f))

        assert len(remaining) < 20
        assert rec.family_id == "F13"
        assert rec.gate == "G6"
        assert len(rec.affected_items) > 0

    def test_ground_truth_describes_filter(self, tmp_path):
        _make_ranked_csv(tmp_path, n=20)
        params = {"filter_name": "MW", "threshold": 150.0}
        rec = self.inj.inject(tmp_path, params, seed=0, instance_index=0)
        assert "MW" in rec.ground_truth or "filter" in rec.ground_truth.lower()

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            self.inj.inject(tmp_path, {"filter_name": "MW", "threshold": 100.0},
                            seed=0, instance_index=0)

    def test_no_camp_in_ground_truth(self, tmp_path):
        _make_ranked_csv(tmp_path, n=20)
        params = self.inj.spec.param_grid[0]
        rec = self.inj.inject(tmp_path, params, seed=0, instance_index=0)
        # "camp" must not appear as an assay reference (note: "no assay data" is an allowed disclosure phrase)
        assert "camp" not in rec.ground_truth.lower()
        assert "ic50" not in rec.ground_truth.lower()
        assert "ec50" not in rec.ground_truth.lower()


# ---------------------------------------------------------------------------
# F14 SCAFFOLD_DUPLICATION
# ---------------------------------------------------------------------------

class TestF14Spec:
    def setup_method(self):
        self.inj = F14ScaffoldDuplication()

    def test_family_id(self):
        assert self.inj.spec.family_id == "F14"

    def test_gate(self):
        assert self.inj.spec.gate == "G6"

    def test_severity(self):
        assert self.inj.spec.severity == "medium"

    def test_param_grid_length(self):
        assert len(self.inj.spec.param_grid) == 12

    def test_no_assay_leakage(self):
        self.inj.verify_no_assay_leakage()


class TestF14Inject:
    def setup_method(self):
        self.inj = F14ScaffoldDuplication()

    def test_panel_grows(self, tmp_path):
        _make_selected_panel(tmp_path, n=15)
        params = self.inj.spec.param_grid[0]
        n_dup = params.get("n_duplicates", 2)
        rec = self.inj.inject(tmp_path, params, seed=0, instance_index=0)

        with open(tmp_path / "selected_panel.csv", newline="") as f:
            new_count = sum(1 for _ in csv.DictReader(f))

        assert new_count > 15
        assert rec.family_id == "F14"

    def test_returns_injection_record(self, tmp_path):
        _make_selected_panel(tmp_path, n=15)
        params = self.inj.spec.param_grid[0]
        rec = self.inj.inject(tmp_path, params, seed=42, instance_index=0)
        assert isinstance(rec.affected_items, list)

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            self.inj.inject(tmp_path, self.inj.spec.param_grid[0], seed=0, instance_index=0)

"""
Tests for stub injectors (F06–F08, F11–F12, F15–F16).

Each stub must:
  - Have a valid spec (family_id, gate, severity, param_grid)
  - Pass verify_no_assay_leakage()
  - Raise NotImplementedError on inject()
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from fault_injection.g2.f06_wrong_protonation import F06WrongProtonation
from fault_injection.g3.f07_off_pocket_pose import F07OffPocketPose
from fault_injection.g3.f08_rmsd_exceeded import F08RmsdExceeded
from fault_injection.g5.f11_label_leakage import F11LabelLeakage
from fault_injection.g5.f12_eval_set_mismatch import F12EvalSetMismatch
from fault_injection.g7.f15_truncated_traj import F15TruncatedTraj
from fault_injection.g7.f16_residue_index_error import F16ResidueIndexError


_STUBS = [
    ("F06", "G2", F06WrongProtonation),
    ("F07", "G3", F07OffPocketPose),
    ("F08", "G3", F08RmsdExceeded),
    ("F11", "G5", F11LabelLeakage),
    ("F12", "G5", F12EvalSetMismatch),
    ("F15", "G7", F15TruncatedTraj),
    ("F16", "G7", F16ResidueIndexError),
]


@pytest.mark.parametrize("family_id,gate,cls", _STUBS)
class TestStubSpec:
    def test_family_id(self, family_id, gate, cls, tmp_path):
        inj = cls()
        assert inj.spec.family_id == family_id

    def test_gate(self, family_id, gate, cls, tmp_path):
        inj = cls()
        assert inj.spec.gate == gate

    def test_param_grid_not_empty(self, family_id, gate, cls, tmp_path):
        inj = cls()
        assert len(inj.spec.param_grid) > 0

    def test_no_assay_leakage(self, family_id, gate, cls, tmp_path):
        inj = cls()
        inj.verify_no_assay_leakage()  # must not raise

    def test_inject_raises_not_implemented(self, family_id, gate, cls, tmp_path):
        inj = cls()
        params = inj.spec.param_grid[0]
        with pytest.raises(NotImplementedError):
            inj.inject(tmp_path, params, seed=0, instance_index=0)

    def test_instance_id_format(self, family_id, gate, cls, tmp_path):
        inj = cls()
        iid = inj.instance_id(0)
        assert iid == f"{family_id}_001"

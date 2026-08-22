"""
manifest.py — Registry of all 18 approved development fault families.

Provides a single dict mapping family_id → FaultInjector instance.
Implemented families raise on inject(); stub families raise NotImplementedError.
"""

from __future__ import annotations

from fault_injection.base import FaultInjector

from fault_injection.g0.f01_id_duplicate import F01IdDuplicate
from fault_injection.g0.f02_hash_mismatch import F02HashMismatch
from fault_injection.g1.f03_stereo_loss import F03StereoLoss
from fault_injection.g1.f04_salt_retained import F04SaltRetained
from fault_injection.g1.f05_invalid_valence import F05InvalidValence
from fault_injection.g2.f06_wrong_protonation import F06WrongProtonation
from fault_injection.g3.f07_off_pocket_pose import F07OffPocketPose
from fault_injection.g3.f08_rmsd_exceeded import F08RmsdExceeded
from fault_injection.g4.f09_score_direction import F09ScoreDirection
from fault_injection.g4.f10_mw_score_bias import F10MwScoreBias
from fault_injection.g5.f11_label_leakage import F11LabelLeakage
from fault_injection.g5.f12_eval_set_mismatch import F12EvalSetMismatch
from fault_injection.g6.f13_undeclared_external_label_constraint import (
    F13UndeclaredExternalLabelConstraint,
)
from fault_injection.g6.f14_scaffold_duplication import F14ScaffoldDuplication
from fault_injection.g7.f15_truncated_traj import F15TruncatedTraj
from fault_injection.g7.f16_residue_index_error import F16ResidueIndexError
from fault_injection.g8.f17_stale_figure import F17StaleFigure
from fault_injection.g8.f18_value_mismatch import F18ValueMismatch

FAULT_REGISTRY: dict[str, FaultInjector] = {
    "F01": F01IdDuplicate(),
    "F02": F02HashMismatch(),
    "F03": F03StereoLoss(),
    "F04": F04SaltRetained(),
    "F05": F05InvalidValence(),
    "F06": F06WrongProtonation(),    # stub — requires WP5
    "F07": F07OffPocketPose(),       # stub — requires WP5
    "F08": F08RmsdExceeded(),        # stub — requires WP5
    "F09": F09ScoreDirection(),
    "F10": F10MwScoreBias(),
    "F11": F11LabelLeakage(),        # stub — requires WP5
    "F12": F12EvalSetMismatch(),     # stub — requires WP5
    "F13": F13UndeclaredExternalLabelConstraint(),
    "F14": F14ScaffoldDuplication(),
    "F15": F15TruncatedTraj(),       # stub — requires WP7
    "F16": F16ResidueIndexError(),   # stub — requires WP7
    "F17": F17StaleFigure(),
    "F18": F18ValueMismatch(),
}

# Verify all 18 families are registered
assert len(FAULT_REGISTRY) == 18, f"Expected 18 families, got {len(FAULT_REGISTRY)}"

# Verify no assay leakage in any implemented family (stubs pass trivially)
for fid, injector in FAULT_REGISTRY.items():
    injector.verify_no_assay_leakage()


def get_injector(family_id: str) -> FaultInjector:
    """Return the FaultInjector for a given family_id."""
    if family_id not in FAULT_REGISTRY:
        raise KeyError(f"Unknown fault family: {family_id!r}. Valid: {sorted(FAULT_REGISTRY)}")
    return FAULT_REGISTRY[family_id]


def implementable_families() -> list[str]:
    """Return family IDs that have full inject() (not stubs)."""
    stubs = {"F06", "F07", "F08", "F11", "F12", "F15", "F16"}
    return [fid for fid in sorted(FAULT_REGISTRY) if fid not in stubs]


def stub_families() -> list[str]:
    """Return family IDs that are stubs pending WP5/WP7."""
    return sorted({"F06", "F07", "F08", "F11", "F12", "F15", "F16"})

"""
Tests for fault_injection/manifest.py.

Verifies registry completeness, correct class mapping, leakage checks,
and helper functions.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from fault_injection.manifest import (
    FAULT_REGISTRY,
    get_injector,
    implementable_families,
    stub_families,
)


class TestRegistryCompleteness:
    def test_eighteen_families(self):
        assert len(FAULT_REGISTRY) == 18

    def test_all_family_ids_present(self):
        expected = {f"F{i:02d}" for i in range(1, 19)}
        assert set(FAULT_REGISTRY.keys()) == expected

    def test_each_value_is_injector(self):
        from fault_injection.base import FaultInjector
        for fid, inj in FAULT_REGISTRY.items():
            assert isinstance(inj, FaultInjector), f"{fid} is not a FaultInjector"


class TestRegistryFamilyIds:
    """Each injector's spec.family_id must match its registry key."""

    def test_spec_family_id_matches_key(self):
        for key, inj in FAULT_REGISTRY.items():
            assert inj.spec.family_id == key, (
                f"Registry key '{key}' but spec.family_id='{inj.spec.family_id}'"
            )


class TestNoAssayLeakage:
    """manifest-level leakage check runs at import; explicit here for clarity."""

    def test_all_pass_leakage_check(self):
        for fid, inj in FAULT_REGISTRY.items():
            inj.verify_no_assay_leakage()  # must not raise


class TestGetInjector:
    def test_returns_correct_injector(self):
        inj = get_injector("F01")
        assert inj.spec.family_id == "F01"

    def test_unknown_key_raises(self):
        with pytest.raises(KeyError, match="Unknown fault family"):
            get_injector("F99")

    def test_returns_same_instance_each_time(self):
        # Registry is a module-level dict; same object returned
        assert get_injector("F09") is get_injector("F09")


class TestImplementableFamilies:
    def test_returns_eleven_families(self):
        # 18 total − 7 stubs = 11
        assert len(implementable_families()) == 11

    def test_stubs_excluded(self):
        impl = set(implementable_families())
        for sid in stub_families():
            assert sid not in impl

    def test_sorted(self):
        impl = implementable_families()
        assert impl == sorted(impl)


class TestStubFamilies:
    def test_returns_seven_stubs(self):
        assert len(stub_families()) == 7

    def test_correct_stubs(self):
        expected_stubs = {"F06", "F07", "F08", "F11", "F12", "F15", "F16"}
        assert set(stub_families()) == expected_stubs

    def test_sorted(self):
        stubs = stub_families()
        assert stubs == sorted(stubs)


class TestGateDistribution:
    """Sanity-check that gates are distributed across G0–G8."""

    def test_g0_has_two_families(self):
        g0 = [fid for fid, inj in FAULT_REGISTRY.items() if inj.spec.gate == "G0"]
        assert len(g0) == 2

    def test_g4_has_two_families(self):
        g4 = [fid for fid, inj in FAULT_REGISTRY.items() if inj.spec.gate == "G4"]
        assert len(g4) == 2

    def test_g8_has_two_families(self):
        g8 = [fid for fid, inj in FAULT_REGISTRY.items() if inj.spec.gate == "G8"]
        assert len(g8) == 2

    def test_all_gates_are_known(self):
        known_gates = {"G0", "G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8"}
        for fid, inj in FAULT_REGISTRY.items():
            assert inj.spec.gate in known_gates, (
                f"{fid}: unknown gate '{inj.spec.gate}'"
            )

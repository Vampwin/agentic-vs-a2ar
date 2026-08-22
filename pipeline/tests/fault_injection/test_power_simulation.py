"""
Tests for fault_injection/power_simulation.py

Skill: statistical-power v1.0
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

# ── path setup ────────────────────────────────────────────────────────────────
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from fault_injection.power_simulation import (
    J_MAX,
    PowerSimConfig,
    PowerSimResults,
    PowerPoint,
    _family_seer_B,
    _one_replicate,
    _wilson_ci,
    build_sensitivity_table,
    estimate_power,
    find_mde,
    run_power_simulation,
)


# ─── Helpers ──────────────────────────────────────────────────────────────────

def _fast_config(**kwargs) -> PowerSimConfig:
    """Return a PowerSimConfig with small n_sims for fast unit tests."""
    defaults = dict(
        mu_B=0.60,
        icc=0.30,
        k_instances=10,
        n_sims=200,
        seed=7,
        J_range=list(range(2, J_MAX + 1)),
        sensitivity_icc_grid=[0.2, 0.8],
        sensitivity_mu_B_grid=[0.5, 0.8],
    )
    defaults.update(kwargs)
    return PowerSimConfig(**defaults)


def _rng(seed: int = 0) -> np.random.Generator:
    return np.random.default_rng(seed)


# ═════════════════════════════════════════════════════════════════════════════
# TestWilsonCI
# ═════════════════════════════════════════════════════════════════════════════

class TestWilsonCI:
    def test_all_successes(self):
        lo, hi = _wilson_ci(100, 100)
        assert lo > 0.9
        assert hi >= 0.999

    def test_all_failures(self):
        lo, hi = _wilson_ci(0, 100)
        assert lo == 0.0
        assert hi < 0.1

    def test_zero_n(self):
        lo, hi = _wilson_ci(0, 0)
        assert lo == 0.0
        assert hi == 1.0

    def test_half(self):
        lo, hi = _wilson_ci(50, 100)
        assert 0.4 <= lo <= 0.5
        assert 0.5 <= hi <= 0.6

    def test_ordered(self):
        lo, hi = _wilson_ci(30, 100)
        assert lo < hi

    def test_bounds(self):
        lo, hi = _wilson_ci(1, 1000)
        assert 0.0 <= lo <= hi <= 1.0


# ═════════════════════════════════════════════════════════════════════════════
# TestFamilySeerB
# ═════════════════════════════════════════════════════════════════════════════

class TestFamilySeerB:
    def test_shape(self):
        rng = _rng(0)
        result = _family_seer_B(J=7, k=10, mu_B=0.5, icc=0.4, rng=rng)
        assert result.shape == (7,)

    def test_values_in_unit_interval(self):
        rng = _rng(1)
        result = _family_seer_B(J=20, k=10, mu_B=0.6, icc=0.5, rng=rng)
        assert np.all(result >= 0.0)
        assert np.all(result <= 1.0)

    def test_icc_zero_deterministic(self):
        """ICC=0: all families get the same escape rate mu_B."""
        rng = _rng(2)
        result = _family_seer_B(J=50, k=1, mu_B=0.7, icc=0.0, rng=rng)
        # All families have exactly mu_B for their draw probability
        # The mean of iid Bernoulli(0.7) k=1 → each entry is 0 or 1,
        # but the generating p_j = 0.7 for all j (no between-family variation).
        # We can't test exact equality with k=1, but variance should be small
        # relative to ICC=1.0 scenario.
        rng2 = _rng(3)
        result_icc1 = _family_seer_B(J=50, k=1, mu_B=0.7, icc=1.0, rng=rng2)
        # ICC=0 gives consistent p_j → low family-mean variance
        # ICC=1 gives random p_j ~ Bernoulli → high family-mean variance
        assert result.var() <= result_icc1.var() + 0.05  # allow small tolerance

    def test_icc_one_binary_outcomes(self):
        """ICC=1: each family draw from Bernoulli(mu_B), then all instances same."""
        rng = _rng(4)
        result = _family_seer_B(J=100, k=20, mu_B=0.5, icc=1.0, rng=rng)
        # Each family should have seer = 0 or 1 only
        # (since p_j ~ Bernoulli(0.5) → either p_j=0 or p_j=1 → all instances same)
        for val in result:
            assert val == 0.0 or val == 1.0

    def test_icc_zero_mean_approx_mu_B(self):
        rng = _rng(5)
        result = _family_seer_B(J=200, k=20, mu_B=0.4, icc=0.0, rng=rng)
        assert abs(result.mean() - 0.4) < 0.1

    def test_icc_midrange_mean_approx_mu_B(self):
        rng = _rng(6)
        result = _family_seer_B(J=500, k=10, mu_B=0.65, icc=0.3, rng=rng)
        assert abs(result.mean() - 0.65) < 0.05

    def test_icc_one_mean_approx_mu_B(self):
        rng = _rng(7)
        result = _family_seer_B(J=500, k=5, mu_B=0.7, icc=1.0, rng=rng)
        assert abs(result.mean() - 0.7) < 0.1

    def test_numerical_guard_small_icc(self):
        """Very small ICC should not raise an error."""
        rng = _rng(8)
        result = _family_seer_B(J=5, k=10, mu_B=0.5, icc=1e-6, rng=rng)
        assert result.shape == (5,)

    def test_numerical_guard_large_icc_below_one(self):
        """ICC just below 1.0 should not raise an error."""
        rng = _rng(9)
        result = _family_seer_B(J=5, k=10, mu_B=0.5, icc=0.9999, rng=rng)
        assert result.shape == (5,)


# ═════════════════════════════════════════════════════════════════════════════
# TestOneReplicate
# ═════════════════════════════════════════════════════════════════════════════

class TestOneReplicate:
    def test_perfect_containment_always_rejects(self):
        """mu_B=0.0 → all contained → ΔSEER=1.0 >> SESOI → always reject."""
        rng = _rng(10)
        results = [
            _one_replicate(J=5, k=10, mu_B=0.0, mu_D=1.0,
                           icc=0.0, sesoi=0.15, alpha=0.05, rng=rng)
            for _ in range(50)
        ]
        assert all(results)

    def test_no_effect_never_rejects(self):
        """mu_B=1.0 = mu_D → ΔSEER=0 < SESOI → never reject."""
        rng = _rng(11)
        results = [
            _one_replicate(J=5, k=10, mu_B=1.0, mu_D=1.0,
                           icc=0.0, sesoi=0.15, alpha=0.05, rng=rng)
            for _ in range(50)
        ]
        assert not any(results)

    def test_returns_bool(self):
        rng = _rng(12)
        result = _one_replicate(J=4, k=5, mu_B=0.6, mu_D=1.0,
                                icc=0.3, sesoi=0.15, alpha=0.05, rng=rng)
        assert isinstance(result, bool)

    def test_large_effect_rejects_often(self):
        """True ΔSEER=0.6 >> SESOI=0.15: high reject rate."""
        rng = _rng(13)
        hits = sum(
            _one_replicate(J=8, k=15, mu_B=0.4, mu_D=1.0,
                           icc=0.1, sesoi=0.15, alpha=0.05, rng=rng)
            for _ in range(200)
        )
        assert hits / 200 > 0.70

    def test_effect_at_boundary_low_power(self):
        """True ΔSEER ≈ SESOI: power near alpha level."""
        rng = _rng(14)
        hits = sum(
            _one_replicate(J=4, k=10, mu_B=0.85, mu_D=1.0,
                           icc=0.5, sesoi=0.15, alpha=0.05, rng=rng)
            for _ in range(200)
        )
        assert hits / 200 < 0.50

    def test_zero_variance_deterministic_reject(self):
        """ICC=0 + extreme mu_B → std_delta=0 path → deterministic."""
        rng = _rng(15)
        result = _one_replicate(J=4, k=100, mu_B=0.0, mu_D=1.0,
                                icc=0.0, sesoi=0.15, alpha=0.05, rng=rng)
        assert result is True

    def test_zero_variance_deterministic_fail(self):
        rng = _rng(16)
        result = _one_replicate(J=4, k=100, mu_B=1.0, mu_D=1.0,
                                icc=0.0, sesoi=0.15, alpha=0.05, rng=rng)
        assert result is False


# ═════════════════════════════════════════════════════════════════════════════
# TestEstimatePower
# ═════════════════════════════════════════════════════════════════════════════

class TestEstimatePower:
    def test_returns_power_point(self):
        cfg = _fast_config(mu_B=0.50, icc=0.20)
        pt = estimate_power(J=5, config=cfg)
        assert isinstance(pt, PowerPoint)
        assert pt.J == 5

    def test_power_in_unit_interval(self):
        cfg = _fast_config()
        pt = estimate_power(J=6, config=cfg)
        assert 0.0 <= pt.power <= 1.0

    def test_ci_ordered(self):
        cfg = _fast_config()
        pt = estimate_power(J=5, config=cfg)
        assert pt.ci_low <= pt.power <= pt.ci_high

    def test_high_power_large_effect_low_icc(self):
        """Large effect + low ICC + many families → power ≥ 0.80."""
        cfg = _fast_config(mu_B=0.30, icc=0.05, n_sims=400)
        pt = estimate_power(J=J_MAX, config=cfg)
        assert pt.power >= 0.70  # allow MC variance

    def test_low_power_small_effect_high_icc(self):
        """Effect barely above SESOI + max clustering + few families → low power."""
        cfg = _fast_config(mu_B=0.84, icc=1.0, n_sims=300)
        pt = estimate_power(J=3, config=cfg)
        assert pt.power < 0.5

    def test_power_increases_with_J(self):
        """Power should be non-decreasing in J (loosely, over large range)."""
        cfg = _fast_config(mu_B=0.65, icc=0.30, n_sims=300)
        pt_small = estimate_power(J=2, config=cfg)
        pt_large = estimate_power(J=J_MAX, config=cfg)
        # Allow some MC noise; trend should hold
        assert pt_large.power >= pt_small.power - 0.10

    def test_n_sims_recorded(self):
        cfg = _fast_config(n_sims=150)
        pt = estimate_power(J=4, config=cfg)
        assert pt.n_sims == 150


# ═════════════════════════════════════════════════════════════════════════════
# TestFindMDE
# ═════════════════════════════════════════════════════════════════════════════

class TestFindMDE:
    def test_mde_gt_sesoi(self):
        """MDE must exceed SESOI (otherwise any effect > 0 would be sufficient)."""
        cfg = _fast_config(icc=0.30)
        mde = find_mde(J=J_MAX, config=cfg, tol=0.05, n_sims_mde=200)
        assert mde > cfg.sesoi

    def test_mde_lt_one(self):
        cfg = _fast_config()
        mde = find_mde(J=J_MAX, config=cfg, tol=0.05, n_sims_mde=200)
        assert mde < 1.0

    def test_mde_decreases_with_more_families(self):
        """More families → smaller MDE."""
        cfg = _fast_config(icc=0.30, n_sims=100)
        mde_small = find_mde(J=3, config=cfg, tol=0.05, n_sims_mde=200)
        mde_large = find_mde(J=J_MAX, config=cfg, tol=0.05, n_sims_mde=200)
        assert mde_large <= mde_small + 0.05  # allow MC noise

    def test_mde_increases_with_icc(self):
        """Higher ICC → worse precision → larger MDE needed."""
        cfg_lo = _fast_config(icc=0.10)
        cfg_hi = _fast_config(icc=0.90)
        mde_lo = find_mde(J=J_MAX, config=cfg_lo, tol=0.05, n_sims_mde=200)
        mde_hi = find_mde(J=J_MAX, config=cfg_hi, tol=0.05, n_sims_mde=200)
        assert mde_hi >= mde_lo - 0.05

    def test_mde_returns_float(self):
        cfg = _fast_config()
        mde = find_mde(J=J_MAX, config=cfg, tol=0.1, n_sims_mde=100)
        assert isinstance(mde, float)


# ═════════════════════════════════════════════════════════════════════════════
# TestBuildSensitivityTable
# ═════════════════════════════════════════════════════════════════════════════

class TestBuildSensitivityTable:
    @pytest.fixture
    def small_cfg(self):
        return _fast_config(
            sensitivity_icc_grid=[0.20, 0.80],
            sensitivity_mu_B_grid=[0.50, 0.80],
        )

    def test_row_count(self, small_cfg):
        table = build_sensitivity_table(small_cfg, n_sims_grid=100)
        assert len(table) == 4  # 2 ICC × 2 mu_B

    def test_required_keys(self, small_cfg):
        table = build_sensitivity_table(small_cfg, n_sims_grid=100)
        for row in table:
            for key in ("icc", "mu_B", "true_delta_SEER", "J_max", "power",
                        "ci_low", "ci_high"):
                assert key in row

    def test_j_max_in_all_rows(self, small_cfg):
        table = build_sensitivity_table(small_cfg, n_sims_grid=100)
        for row in table:
            assert row["J_max"] == J_MAX

    def test_power_in_unit_interval(self, small_cfg):
        table = build_sensitivity_table(small_cfg, n_sims_grid=100)
        for row in table:
            assert 0.0 <= row["power"] <= 1.0

    def test_true_delta_computed_correctly(self, small_cfg):
        table = build_sensitivity_table(small_cfg, n_sims_grid=100)
        for row in table:
            expected = round(1.0 - row["mu_B"], 3)
            assert abs(row["true_delta_SEER"] - expected) < 1e-6

    def test_low_icc_higher_power_than_high_icc(self, small_cfg):
        """ICC=0.20 row should generally have >= power vs ICC=0.80 for same mu_B."""
        table = build_sensitivity_table(small_cfg, n_sims_grid=300)
        # Group by mu_B
        for mu_B in [0.50, 0.80]:
            rows = [r for r in table if r["mu_B"] == mu_B]
            lo = next(r["power"] for r in rows if r["icc"] == 0.20)
            hi = next(r["power"] for r in rows if r["icc"] == 0.80)
            assert lo >= hi - 0.15  # allow MC variance

    def test_lower_mu_B_higher_power(self, small_cfg):
        """mu_B=0.50 (larger delta) should have >= power vs mu_B=0.80."""
        table = build_sensitivity_table(small_cfg, n_sims_grid=300)
        for icc in [0.20, 0.80]:
            rows = [r for r in table if r["icc"] == icc]
            p_low_mu = next(r["power"] for r in rows if r["mu_B"] == 0.50)
            p_high_mu = next(r["power"] for r in rows if r["mu_B"] == 0.80)
            assert p_low_mu >= p_high_mu - 0.15


# ═════════════════════════════════════════════════════════════════════════════
# TestRunPowerSimulation — integration
# ═════════════════════════════════════════════════════════════════════════════

class TestRunPowerSimulation:
    @pytest.fixture(scope="class")
    def results_small_effect(self):
        cfg = _fast_config(
            mu_B=0.70, icc=0.50, n_sims=200,
            sensitivity_icc_grid=[0.30],
            sensitivity_mu_B_grid=[0.70],
        )
        return run_power_simulation(cfg)

    @pytest.fixture(scope="class")
    def results_large_effect(self):
        cfg = _fast_config(
            mu_B=0.30, icc=0.10, n_sims=300,
            sensitivity_icc_grid=[0.10],
            sensitivity_mu_B_grid=[0.30],
        )
        return run_power_simulation(cfg)

    def test_returns_power_sim_results(self, results_small_effect):
        assert isinstance(results_small_effect, PowerSimResults)

    def test_power_curve_length(self, results_small_effect):
        cfg = _fast_config()
        assert len(results_small_effect.power_curve) == len(cfg.J_range)

    def test_power_curve_has_required_keys(self, results_small_effect):
        for row in results_small_effect.power_curve:
            for key in ("J", "power", "ci_low", "ci_high", "n_sims"):
                assert key in row

    def test_power_at_j_max_in_unit_interval(self, results_small_effect):
        assert 0.0 <= results_small_effect.power_at_J_max <= 1.0

    def test_mde_at_j_max_gt_sesoi(self, results_small_effect):
        cfg = _fast_config()
        assert results_small_effect.mde_at_J_max > cfg.sesoi

    def test_recommendation_contains_gate9(self, results_small_effect):
        assert "Gate 9" in results_small_effect.recommendation

    def test_icc_note_present(self, results_small_effect):
        assert len(results_small_effect.icc_note) > 20

    def test_sensitivity_table_present(self, results_small_effect):
        assert isinstance(results_small_effect.sensitivity_table, list)
        assert len(results_small_effect.sensitivity_table) >= 1

    def test_config_in_output(self, results_small_effect):
        assert "mu_B" in results_small_effect.config
        assert "icc" in results_small_effect.config

    def test_min_j_for_target_is_none_or_int(self, results_small_effect):
        mjt = results_small_effect.min_J_for_target
        assert mjt is None or isinstance(mjt, int)

    def test_min_j_within_j_range_when_found(self, results_small_effect):
        mjt = results_small_effect.min_J_for_target
        cfg = _fast_config()
        if mjt is not None:
            assert cfg.J_range[0] <= mjt <= cfg.J_range[-1]

    def test_large_effect_reaches_target(self, results_large_effect):
        """Large effect + low ICC: should achieve target power at J_MAX."""
        assert results_large_effect.power_at_J_max > 0.65  # allow MC noise

    def test_large_effect_min_j_found(self, results_large_effect):
        """Large effect + low ICC: min_J_for_target should be found."""
        assert results_large_effect.min_J_for_target is not None

    def test_power_curve_j_values_match_range(self, results_small_effect):
        cfg = _fast_config()
        curve_j = [row["J"] for row in results_small_effect.power_curve]
        assert curve_j == list(cfg.J_range)

    def test_small_effect_high_icc_recommendation_mentions_infeasible_or_power(
        self, results_small_effect
    ):
        rec = results_small_effect.recommendation
        # Should mention either not achievable or power value
        assert "power" in rec.lower() or "feasible" in rec.lower()

    def test_power_curve_values_monotone_tendency(self, results_large_effect):
        """Power should not decrease sharply from J=8 to J_MAX (large effect)."""
        curve = results_large_effect.power_curve
        p_mid = next(r["power"] for r in curve if r["J"] == 8)
        p_max = next(r["power"] for r in curve if r["J"] == J_MAX)
        assert p_max >= p_mid - 0.15  # allow MC noise

    def test_all_power_values_in_unit_interval(self, results_small_effect):
        for row in results_small_effect.power_curve:
            assert 0.0 <= row["power"] <= 1.0

    def test_all_ci_ordered(self, results_small_effect):
        for row in results_small_effect.power_curve:
            assert row["ci_low"] <= row["power"] <= row["ci_high"]


# ═════════════════════════════════════════════════════════════════════════════
# TestProtocolConstants
# ═════════════════════════════════════════════════════════════════════════════

class TestProtocolConstants:
    def test_sesoi_frozen(self):
        from fault_injection.power_simulation import SESOI
        assert SESOI == 0.15

    def test_alpha_frozen(self):
        from fault_injection.power_simulation import ALPHA
        assert ALPHA == 0.05

    def test_power_target_frozen(self):
        from fault_injection.power_simulation import POWER_TARGET
        assert POWER_TARGET == 0.80

    def test_j_max_frozen(self):
        from fault_injection.power_simulation import J_MAX
        assert J_MAX == 11

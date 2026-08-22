"""
power_simulation.py — Monte Carlo power simulation for Gate 9.

Primary analysis
----------------
Does the contracts component reduce SEER by ≥ SESOI = 0.15?

Design: 2×2 blocked factorial, blocking variable = fault family.
Test unit: family (not individual instance).
Outcome: binary escaped (0 = caught, 1 = escaped) per instance.

Test (pre-specified, Gate 9)
----------------------------
Paired family-level one-sample t-test on ΔSEER = SEER_D − SEER_B.
  H0: mean(ΔSEER) ≤ SESOI = 0.15
  H1: mean(ΔSEER) > SESOI
Reject H0 if t > t_{α, J−1} (one-tailed, α = 0.05).

Data-generating model (Beta-Binomial)
--------------------------------------
For each simulation replicate with J families and k instances per family:

  1. p_B_j ~ Beta(μ_B·c, (1−μ_B)·c)   where c = (1−ICC)/ICC
       μ_B = true treatment SEER (condition B: contracts only)
       ICC = intraclass correlation (between-family variance in escape rate)
       Special cases: ICC=0 → p_B_j = μ_B (no clustering);
                      ICC≥1 → p_B_j ~ Bernoulli(μ_B) (max clustering)
  2. For each family j: y_{ij} ~ Bernoulli(p_B_j), i = 1..k
     seer_B_j = mean(y_j)
  3. SEER_D_j = 1.0 (all escape under no defense — deterministic)
  4. delta_j  = SEER_D_j − seer_B_j

Test statistic:
  t = (mean(delta) − SESOI) / (std(delta) / sqrt(J))
  Reject if t > t_{0.05, J−1}

Outputs
-------
• Power curve: power vs. J for fixed ICC and μ_B
• Sensitivity table: power at J_max across (ICC × μ_B) grid
• Minimum J for 80% power (or "infeasible" if > J_max)
• MDE at J_max (minimum detectable SEER reduction at 80% power)
• Recommendation text for Gate 9 report

Skill: statistical-power v1.0  (simulation-based power, cluster design)
ADR: ADR-AIIM-019, ADR-AIIM-023, ADR-AIIM-024
"""

from __future__ import annotations

import json
import logging
import math
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
from scipy import stats

logger = logging.getLogger(__name__)

# ─── Protocol constants (frozen, ADR-AIIM-023 / Gate 9) ──────────────────────

SESOI: float = 0.15   # smallest effect of interest (absolute SEER reduction)
ALPHA: float = 0.05   # one-tailed significance level
POWER_TARGET: float = 0.80
J_MAX: int = 11       # implementable development families (F01–F18 minus 7 stubs)
J_MAX_LABEL = "11 (all implementable development families)"


# ─── Dataclasses ─────────────────────────────────────────────────────────────

@dataclass
class PowerSimConfig:
    mu_B: float = 0.688        # pilot SEER under condition B (3-family estimate)
    mu_D: float = 1.0          # baseline SEER under D (deterministic, no defense)
    icc: float = 0.50          # intraclass correlation (conservatively below pilot 1.0)
    k_instances: int = 10      # instances per family (pilot default)
    sesoi: float = SESOI
    alpha: float = ALPHA
    power_target: float = POWER_TARGET
    n_sims: int = 5_000        # Monte Carlo replicates
    seed: int = 42
    J_range: list[int] = field(default_factory=lambda: list(range(2, J_MAX + 1)))
    sensitivity_icc_grid: list[float] = field(
        default_factory=lambda: [0.20, 0.40, 0.60, 0.80, 1.00]
    )
    sensitivity_mu_B_grid: list[float] = field(
        default_factory=lambda: [0.50, 0.60, 0.70, 0.80, 0.90]
    )


@dataclass
class PowerPoint:
    J: int
    power: float
    ci_low: float
    ci_high: float
    n_sims: int

    def __str__(self) -> str:
        return (
            f"J={self.J}: power={self.power:.3f} "
            f"(95% MC CI {self.ci_low:.3f}–{self.ci_high:.3f})"
        )


@dataclass
class PowerSimResults:
    config: dict
    power_curve: list[dict]           # [{J, power, ci_low, ci_high}, ...]
    min_J_for_target: Optional[int]   # smallest J reaching power_target, or None
    power_at_J_max: float
    mde_at_J_max: float               # MDE at J=J_MAX for the primary ICC
    sensitivity_table: list[dict]     # [{icc, mu_B, J_max_power, mu_D_J_max}, ...]
    recommendation: str
    icc_note: str


# ─── Wilson CI for a proportion ───────────────────────────────────────────────

def _wilson_ci(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    denom = 1 + z ** 2 / n
    center = (p + z ** 2 / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / denom
    return (max(0.0, center - half), min(1.0, center + half))


# ─── Beta-Binomial data generator ─────────────────────────────────────────────

def _family_seer_B(
    J: int,
    k: int,
    mu_B: float,
    icc: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Generate family-level SEER_B estimates under the Beta-Binomial model.

    Returns array of shape (J,) with seer_B[j] = mean of k binary outcomes.
    """
    if icc >= 1.0:
        # Maximum clustering: each family is a single Bernoulli draw
        p_j = rng.binomial(1, mu_B, J).astype(float)
    elif icc <= 0.0:
        # No clustering: all families have the same escape rate
        p_j = np.full(J, mu_B)
    else:
        # Beta distribution: c = (1-ICC)/ICC
        c = (1.0 - icc) / icc
        alpha_B = mu_B * c
        beta_B = (1.0 - mu_B) * c
        # Numerical guard: Beta params must be > 0
        alpha_B = max(alpha_B, 1e-6)
        beta_B = max(beta_B, 1e-6)
        p_j = rng.beta(alpha_B, beta_B, J)

    # Instance-level outcomes for each family
    seer = np.empty(J)
    for j in range(J):
        outcomes = rng.binomial(1, float(p_j[j]), k)
        seer[j] = outcomes.mean()
    return seer


# ─── Single simulation replicate ─────────────────────────────────────────────

def _one_replicate(
    J: int,
    k: int,
    mu_B: float,
    mu_D: float,
    icc: float,
    sesoi: float,
    alpha: float,
    rng: np.random.Generator,
) -> bool:
    """
    Simulate one experiment and return True iff H0 is rejected.

    H0: mean(ΔSEER) ≤ sesoi  (contracts don't reduce SEER enough)
    H1: mean(ΔSEER) > sesoi
    Test: one-tailed one-sample t-test on family-level ΔSEER.
    """
    seer_B = _family_seer_B(J, k, mu_B, icc, rng)
    seer_D = np.full(J, mu_D)
    delta = seer_D - seer_B  # ΔSEER: how much contracts reduced escape rate

    mean_d = delta.mean()
    std_d = delta.std(ddof=1)

    if std_d == 0.0:
        # Deterministic outcome (e.g., ICC=0, all families identical)
        return bool(mean_d > sesoi)

    t_stat = (mean_d - sesoi) / (std_d / math.sqrt(J))
    t_crit = stats.t.ppf(1.0 - alpha, df=J - 1)
    return bool(t_stat > t_crit)


# ─── Power at fixed J ────────────────────────────────────────────────────────

def estimate_power(
    J: int,
    config: PowerSimConfig,
) -> PowerPoint:
    """Estimate power at a given number of families J."""
    rng = np.random.default_rng(config.seed + J)  # stable per J
    hits = sum(
        1 for _ in range(config.n_sims)
        if _one_replicate(
            J=J, k=config.k_instances,
            mu_B=config.mu_B, mu_D=config.mu_D,
            icc=config.icc, sesoi=config.sesoi,
            alpha=config.alpha, rng=rng,
        )
    )
    lo, hi = _wilson_ci(hits, config.n_sims)
    return PowerPoint(
        J=J,
        power=hits / config.n_sims,
        ci_low=lo,
        ci_high=hi,
        n_sims=config.n_sims,
    )


# ─── MDE at fixed J ──────────────────────────────────────────────────────────

def _estimate_power_at_effect(
    J: int,
    k: int,
    mu_B: float,       # true treatment SEER (= 1 - true_delta)
    icc: float,
    sesoi: float,
    alpha: float,
    n_sims: int,
    seed: int,
) -> float:
    rng = np.random.default_rng(seed)
    hits = sum(
        1 for _ in range(n_sims)
        if _one_replicate(J=J, k=k, mu_B=mu_B, mu_D=1.0,
                          icc=icc, sesoi=sesoi, alpha=alpha, rng=rng)
    )
    return hits / n_sims


def find_mde(
    J: int,
    config: PowerSimConfig,
    tol: float = 0.01,
    n_sims_mde: int = 2_000,
) -> float:
    """
    Find the minimum detectable SEER reduction (true_delta) that gives
    power ≥ power_target at the given J.

    Binary search over mu_B = 1 - true_delta in [sesoi, 1.0].
    Returns the true_delta (SEER reduction) at the 80% power threshold.
    """
    # Binary search: lower mu_B → larger delta → higher power
    # We want the highest mu_B (smallest delta) that still gives target power
    lo_mu, hi_mu = config.sesoi, 1.0 - config.sesoi

    for _ in range(30):
        mid_mu = (lo_mu + hi_mu) / 2
        pwr = _estimate_power_at_effect(
            J=J, k=config.k_instances,
            mu_B=mid_mu, icc=config.icc,
            sesoi=config.sesoi, alpha=config.alpha,
            n_sims=n_sims_mde, seed=config.seed + 9999,
        )
        if pwr >= config.power_target:
            lo_mu = mid_mu  # can afford a larger mu_B (smaller delta)
        else:
            hi_mu = mid_mu  # need smaller mu_B (larger delta)
        if hi_mu - lo_mu < tol:
            break

    mde = 1.0 - lo_mu  # minimum detectable SEER reduction
    return round(mde, 3)


# ─── Sensitivity table ────────────────────────────────────────────────────────

def build_sensitivity_table(
    config: PowerSimConfig,
    n_sims_grid: int = 2_000,
) -> list[dict]:
    """
    Power at J_MAX for each (ICC, μ_B) combination in the sensitivity grid.
    """
    rows: list[dict] = []
    for icc in config.sensitivity_icc_grid:
        for mu_B in config.sensitivity_mu_B_grid:
            cfg_cell = PowerSimConfig(
                mu_B=mu_B,
                mu_D=config.mu_D,
                icc=icc,
                k_instances=config.k_instances,
                sesoi=config.sesoi,
                alpha=config.alpha,
                power_target=config.power_target,
                n_sims=n_sims_grid,
                seed=config.seed,
            )
            pt = estimate_power(J=J_MAX, config=cfg_cell)
            true_delta = round(config.mu_D - mu_B, 3)
            rows.append({
                "icc": icc,
                "mu_B": mu_B,
                "true_delta_SEER": true_delta,
                "J_max": J_MAX,
                "power": round(pt.power, 3),
                "ci_low": round(pt.ci_low, 3),
                "ci_high": round(pt.ci_high, 3),
            })
            logger.debug(
                "Sensitivity: ICC=%.2f μ_B=%.2f → power=%.3f",
                icc, mu_B, pt.power,
            )
    return rows


# ─── Main simulation ──────────────────────────────────────────────────────────

def run_power_simulation(config: PowerSimConfig) -> PowerSimResults:
    """
    Run the full Gate 9 power simulation.

    1. Power curve: estimate power at each J in config.J_range.
    2. Find minimum J for target power (or report infeasible).
    3. MDE at J_MAX.
    4. Sensitivity table over (ICC × μ_B) grid.
    5. Draft recommendation for Gate 9.
    """
    logger.info(
        "Power simulation: μ_B=%.3f ICC=%.2f k=%d α=%.2f n_sims=%d",
        config.mu_B, config.icc, config.k_instances,
        config.alpha, config.n_sims,
    )

    # ── 1. Power curve ────────────────────────────────────────────────────────
    power_curve: list[PowerPoint] = []
    for J in config.J_range:
        pt = estimate_power(J, config)
        power_curve.append(pt)
        logger.info("  %s", pt)

    # ── 2. Minimum J ──────────────────────────────────────────────────────────
    min_J: Optional[int] = None
    for pt in power_curve:
        if pt.power >= config.power_target:
            min_J = pt.J
            break

    # ── 3. Power at J_MAX ─────────────────────────────────────────────────────
    jmax_curve = [pt for pt in power_curve if pt.J == J_MAX]
    power_at_jmax = jmax_curve[0].power if jmax_curve else float("nan")

    # ── 4. MDE at J_MAX ───────────────────────────────────────────────────────
    mde = find_mde(J_MAX, config, n_sims_mde=2_000)
    logger.info("MDE at J_MAX=%d: %.3f SEER reduction", J_MAX, mde)

    # ── 5. Sensitivity table ──────────────────────────────────────────────────
    logger.info("Building sensitivity table (%d cells)...",
                len(config.sensitivity_icc_grid) * len(config.sensitivity_mu_B_grid))
    sens_table = build_sensitivity_table(config, n_sims_grid=2_000)

    # ── 6. Recommendation text ────────────────────────────────────────────────
    true_delta = config.mu_D - config.mu_B
    if min_J is not None and min_J <= J_MAX:
        reach_str = (
            f"80 % power is achievable at J = {min_J} families "
            f"(ICC = {config.icc:.2f}, μ_B = {config.mu_B:.3f})."
        )
    else:
        reach_str = (
            f"80 % power is NOT achievable within J_max = {J_MAX} families "
            f"at ICC = {config.icc:.2f}.  "
            f"Power at J = {J_MAX} is {power_at_jmax:.3f}.  "
            f"Consider: (a) accepting lower power for the pilot, "
            f"(b) re-estimating ICC with more families, or "
            f"(c) expanding the development family set with Gate 2 approval."
        )

    icc_note = (
        f"ICC = {config.icc:.2f} was the primary simulation value.  "
        f"The feasibility pilot (3 families) gave ICC = 1.00 — likely an overestimate "
        f"due to extreme family-level heterogeneity with few families.  "
        f"A full pilot using all {J_MAX} implementable families will yield a "
        f"more reliable ICC estimate before the Gate 9 sample-size freeze."
    )

    recommendation = "\n".join([
        f"Gate 9 Power Analysis Summary",
        f"─────────────────────────────",
        f"Primary ICC assumption:  {config.icc:.2f}",
        f"True SEER reduction:     {true_delta:.3f} (μ_D={config.mu_D:.2f} − μ_B={config.mu_B:.3f})",
        f"SESOI:                   {config.sesoi:.2f}",
        f"α (one-tailed):          {config.alpha:.2f}",
        f"Power target:            {config.power_target:.0%}",
        f"Max families (J_max):    {J_MAX}",
        f"k instances/family:      {config.k_instances}",
        f"",
        f"Result:  Power at J={J_MAX} = {power_at_jmax:.3f}",
        f"         MDE at J={J_MAX}   = {mde:.3f} absolute SEER reduction",
        f"",
        reach_str,
        f"",
        f"Interpretation:",
        f"  With ICC={config.icc:.2f} and μ_B={config.mu_B:.3f} (pilot estimate),",
        (
            f"  ΔSEER={true_delta:.3f} EXCEEDS SESOI={config.sesoi:.2f}: "
            f"the effect is real but variance is too high at this ICC."
            if true_delta > config.sesoi else
            f"  ΔSEER={true_delta:.3f} IS BELOW SESOI={config.sesoi:.2f}: "
            f"the B-vs-D contrast cannot reject H0 by design — "
            f"contracts alone do not reduce SEER by ≥ {config.sesoi:.2f}. "
            f"Consider the A-vs-D (hybrid vs. neither) contrast as the primary analysis."
        ),
        f"  ICC={config.icc:.2f} → between-family variance dominates; each additional",
        f"  family adds limited independent information.",
        f"",
        f"Next steps (Gate 9 closure):",
        f"  1. Run full pilot with all {J_MAX} implementable families to re-estimate ICC.",
        f"  2. Present this power table to human author for sample-size freeze.",
        f"  3. Gate 8 (compute budget) decision depends on final J and n_sims.",
    ])

    return PowerSimResults(
        config=asdict(config),
        power_curve=[
            {
                "J": pt.J,
                "power": round(pt.power, 4),
                "ci_low": round(pt.ci_low, 4),
                "ci_high": round(pt.ci_high, 4),
                "n_sims": pt.n_sims,
            }
            for pt in power_curve
        ],
        min_J_for_target=min_J,
        power_at_J_max=round(power_at_jmax, 4),
        mde_at_J_max=mde,
        sensitivity_table=sens_table,
        recommendation=recommendation,
        icc_note=icc_note,
    )


# ─── CLI ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import argparse
    import sys

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    parser = argparse.ArgumentParser(
        description="Gate 9 Monte Carlo power simulation (statistical-power skill v1.0)"
    )
    parser.add_argument("--mu-B",    type=float, default=0.688, metavar="F",
                        help="Treatment SEER under condition B (default: 0.688 from pilot)")
    parser.add_argument("--icc",     type=float, default=0.50,  metavar="F",
                        help="Intraclass correlation (default: 0.50; pilot estimate: 1.00)")
    parser.add_argument("--k",       type=int,   default=10,    metavar="N",
                        help="Instances per family (default: 10)")
    parser.add_argument("--n-sims",  type=int,   default=5_000, metavar="N",
                        help="Monte Carlo replicates (default: 5000)")
    parser.add_argument("--seed",    type=int,   default=42,
                        help="Random seed (default: 42)")
    parser.add_argument("--output",  type=Path,  default=None,
                        help="Write JSON results to FILE")
    args = parser.parse_args()

    cfg = PowerSimConfig(
        mu_B=args.mu_B,
        icc=args.icc,
        k_instances=args.k,
        n_sims=args.n_sims,
        seed=args.seed,
    )
    results = run_power_simulation(cfg)

    print(results.recommendation)
    print()
    print("Power curve:")
    for row in results.power_curve:
        bar = "#" * int(row["power"] * 40)
        print(f"  J={row['J']:2d}  {row['power']:.3f}  |{bar}")
    print()
    print(results.icc_note)

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w", encoding="utf-8") as fh:
            json.dump(asdict(results), fh, indent=2, default=str)
        print(f"\nFull results -> {args.output}", file=sys.stderr)

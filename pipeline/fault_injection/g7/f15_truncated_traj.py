"""
F15 TRUNCATED_TRAJ — G7 fault, severity: high

Injection: Truncate the MD trajectory to a fraction of the production length
before analysis, so convergence and energy statistics are computed on insufficient
sampling.  STUB — requires real MD trajectory files (WP5/WP7).

Parameterisation (10 instances):
  fraction_kept: [0.05, 0.10, 0.15, 0.20, 0.25, 0.33, 0.40, 0.50, 0.67, 0.80]

No assay data: only the trajectory frame count is reduced.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from fault_injection.base import FaultInjector, FaultSpec, InjectionRecord

_PARAM_GRID: list[dict[str, Any]] = [
    {"fraction_kept": f}
    for f in [0.05, 0.10, 0.15, 0.20, 0.25, 0.33, 0.40, 0.50, 0.67, 0.80]
]


class F15TruncatedTraj(FaultInjector):
    """Truncate MD trajectory to a fraction of production length (STUB)."""

    @property
    def spec(self) -> FaultSpec:
        return FaultSpec(
            family_id="F15",
            family_name="TRUNCATED_TRAJ",
            gate="G7",
            severity="high",
            description=(
                "Truncate the MD trajectory to a specified fraction of full length "
                "before RMSD/energy analysis.  STUB — requires real trajectory files (WP7)."
            ),
            param_grid=_PARAM_GRID,
        )

    def verify_no_assay_leakage(self) -> None:
        for p in _PARAM_GRID:
            assert "camp" not in str(p).lower()
            assert "assay" not in str(p).lower()

    def inject(
        self,
        workspace: Path,
        params: dict[str, Any],
        seed: int,
        instance_index: int,
    ) -> InjectionRecord:
        raise NotImplementedError(
            "F15 requires real MD trajectory files.  Implement at WP7."
        )


F15TruncatedTraj().verify_no_assay_leakage()

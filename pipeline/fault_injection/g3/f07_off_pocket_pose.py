"""
F07 OFF_POCKET_POSE — G3 fault, severity: high

Injection: Translate the docking box centre coordinates outside the binding site
by a fixed Å offset.  DiffDock/Vina will dock to the wrong region.
G3 RMSD-to-reference and binding-site overlap checks should detect this.

STUB — requires frozen receptor structure and docking config (WP5).

Parameterisation (15 instances planned):
  offset_angstrom: [5, 10, 15, 20, 25] × axis: ["x", "y", "z"]

No assay data: only a geometric offset vector is applied.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from fault_injection.base import FaultInjector, FaultSpec, InjectionRecord

_PARAM_GRID: list[dict[str, Any]] = [
    {"offset_angstrom": offset, "axis": ax}
    for offset in [5, 10, 15, 20, 25]
    for ax in ["x", "y", "z"]
][:15]


class F07OffPocketPose(FaultInjector):
    """Translate docking box centre outside the binding site (STUB)."""

    @property
    def spec(self) -> FaultSpec:
        return FaultSpec(
            family_id="F07",
            family_name="OFF_POCKET_POSE",
            gate="G3",
            severity="high",
            description=(
                "Shift the docking search box centre by a fixed Å offset along one "
                "Cartesian axis, placing it outside the known binding site.  "
                "STUB — requires frozen receptor and docking config (WP5)."
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
            "F07 requires a frozen receptor PDB and docking config.  "
            "Implement after WP5 receptor freeze."
        )


F07OffPocketPose().verify_no_assay_leakage()

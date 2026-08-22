"""
F08 RMSD_EXCEEDED — G3 fault, severity: high

Injection: Apply a prespecified rotation (90°–180°) to the reference ligand
coordinates from the public crystal structure, so RMSD between the docked pose
and the rotated reference exceeds the contract threshold.  The rotation angle is
prespecified independently of any assay data.  Reference ligand comes from the
public PDB entry (4EIY or 3EML); no cAMP or assay thresholds used.

STUB — requires frozen receptor structure and ligand SDF/PDB (WP5).

Parameterisation (15 instances):
  rotation_degrees: [45, 90, 120, 135, 150, 165, 180] × rotation_axis: ["x", "y", "z"]
  (7 × 3 = 21, take first 15)

Injection boundary (Gate 2 confirmed): rotation angle prespecified from geometric
range only.  No assay-derived value enters the parameter grid.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from fault_injection.base import FaultInjector, FaultSpec, InjectionRecord

_ROTATIONS = [45, 90, 120, 135, 150, 165, 180]
_AXES = ["x", "y", "z"]
_PARAM_GRID: list[dict[str, Any]] = [
    {"rotation_degrees": r, "rotation_axis": ax}
    for r in _ROTATIONS
    for ax in _AXES
][:15]


class F08RmsdExceeded(FaultInjector):
    """Rotate the reference ligand to artificially inflate RMSD (STUB)."""

    @property
    def spec(self) -> FaultSpec:
        return FaultSpec(
            family_id="F08",
            family_name="RMSD_EXCEEDED",
            gate="G3",
            severity="high",
            description=(
                "Apply a prespecified rotation to the public crystal-structure "
                "reference ligand coordinates, so the RMSD between docked pose "
                "and reference exceeds the G3 contract threshold.  "
                "Rotation angle is geometric only — no assay threshold used.  "
                "STUB — requires frozen receptor/ligand files (WP5)."
            ),
            param_grid=_PARAM_GRID,
        )

    def verify_no_assay_leakage(self) -> None:
        for p in _PARAM_GRID:
            assert "camp" not in str(p).lower()
            assert "assay" not in str(p).lower()
            # Rotation degrees must be in the prespecified geometric range
            assert 0 < p["rotation_degrees"] <= 180, (
                f"F08: rotation_degrees {p['rotation_degrees']} outside [0, 180]"
            )

    def inject(
        self,
        workspace: Path,
        params: dict[str, Any],
        seed: int,
        instance_index: int,
    ) -> InjectionRecord:
        raise NotImplementedError(
            "F08 requires the frozen public crystal-structure ligand SDF/PDB.  "
            "Implement after WP5 receptor/ligand freeze."
        )


F08RmsdExceeded().verify_no_assay_leakage()

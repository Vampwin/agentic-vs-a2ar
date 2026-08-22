"""
F16 RESIDUE_INDEX_ERROR — G7 fault, severity: high

Injection: Replace the binding-site residue index used for RMSD calculation with
an incorrect index, so the computed binding-site RMSD reflects the wrong region.
STUB — requires real receptor structure and residue definitions (WP7).

No assay data: only a residue integer index is changed.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from fault_injection.base import FaultInjector, FaultSpec, InjectionRecord

_PARAM_GRID: list[dict[str, Any]] = [
    {"residue_offset": offset}
    for offset in [-50, -30, -20, -10, -5, 5, 10, 20, 30, 50]
]


class F16ResidueIndexError(FaultInjector):
    """Use wrong residue index for binding-site RMSD calculation (STUB)."""

    @property
    def spec(self) -> FaultSpec:
        return FaultSpec(
            family_id="F16",
            family_name="RESIDUE_INDEX_ERROR",
            gate="G7",
            severity="high",
            description=(
                "Shift the binding-site residue index by a fixed offset, so RMSD "
                "is computed for the wrong region of the receptor.  "
                "STUB — requires frozen receptor residue definitions (WP7)."
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
            "F16 requires the frozen binding-site residue definition.  Implement at WP7."
        )


F16ResidueIndexError().verify_no_assay_leakage()

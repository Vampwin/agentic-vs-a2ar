"""
F06 WRONG_PROTONATION — G2 fault, severity: medium

Injection: Flip the protonation state of histidine residues in the receptor PDB
file from HID (Nδ-H) to HIE (Nε-H) or vice versa, silently changing the
receptor preparation without updating the preparation log.

STUB — requires a real receptor PDB file (4EIY or 3EML prepared structure).
Implementation deferred to WP5 when receptor structures are frozen.

Parameterisation (10 instances planned):
  n_his_flipped: [1, 2, 3]  × residue_selection: ["first", "random", "all"]

No assay data: only residue names in PDB ATOM records are changed.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from fault_injection.base import FaultInjector, FaultSpec, InjectionRecord

_PARAM_GRID: list[dict[str, Any]] = [
    {"n_his_flipped": n, "residue_selection": sel}
    for n in [1, 2, 3]
    for sel in ["first", "random", "all"]
][:10]


class F06WrongProtonation(FaultInjector):
    """Flip histidine protonation state in the receptor PDB (STUB)."""

    @property
    def spec(self) -> FaultSpec:
        return FaultSpec(
            family_id="F06",
            family_name="WRONG_PROTONATION",
            gate="G2",
            severity="medium",
            description=(
                "Flip HID↔HIE protonation on histidine residues in the prepared "
                "receptor PDB file.  STUB — requires frozen receptor structure (WP5)."
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
            "F06 requires a real receptor PDB file (4EIY/3EML).  "
            "Implement after receptor structures are frozen at WP5."
        )


F06WrongProtonation().verify_no_assay_leakage()

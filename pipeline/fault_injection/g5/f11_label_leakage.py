"""
F11 LABEL_LEAKAGE — G5 fault, severity: critical

Injection: Expose public ChEMBL active IDs in the ranking step input, so the
ranking agent can implicitly learn which structures are known actives.  Uses only
PUBLIC benchmark identifiers — not assay measurements or cAMP values.

STUB — requires the frozen public A2A benchmark file (human approval pending, WP5
Gate freeze of public benchmark labels).

Parameterisation (10 instances planned):
  n_leaked_ids: [5, 10, 15, 20, 25, 30, 40, 50, 75, 100]

No private assay data: only public ChEMBL IDs from a frozen benchmark set.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from fault_injection.base import FaultInjector, FaultSpec, InjectionRecord

_PARAM_GRID: list[dict[str, Any]] = [
    {"n_leaked_ids": n} for n in [5, 10, 15, 20, 25, 30, 40, 50, 75, 100]
]


class F11LabelLeakage(FaultInjector):
    """Expose public benchmark active IDs in ranking input (STUB)."""

    @property
    def spec(self) -> FaultSpec:
        return FaultSpec(
            family_id="F11",
            family_name="LABEL_LEAKAGE",
            gate="G5",
            severity="critical",
            description=(
                "Inject public ChEMBL active IDs into the ranking step input, "
                "leaking benchmark label information.  STUB — requires frozen "
                "public A2A benchmark (WP5)."
            ),
            param_grid=_PARAM_GRID,
        )

    def verify_no_assay_leakage(self) -> None:
        # Note: this fault SIMULATES leakage using PUBLIC benchmark IDs only.
        # No in-house cAMP values may be used. Verified at Gate 2.
        for p in _PARAM_GRID:
            assert "camp" not in str(p).lower()

    def inject(
        self,
        workspace: Path,
        params: dict[str, Any],
        seed: int,
        instance_index: int,
    ) -> InjectionRecord:
        raise NotImplementedError(
            "F11 requires the frozen public A2A benchmark active ID list.  "
            "Implement after benchmark freeze at WP5."
        )


F11LabelLeakage().verify_no_assay_leakage()

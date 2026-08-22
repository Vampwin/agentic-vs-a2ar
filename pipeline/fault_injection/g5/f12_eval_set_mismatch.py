"""
F12 EVAL_SET_MISMATCH — G5 fault, severity: critical

Injection: Swap the actives and decoys lists for benchmark evaluation, so
enrichment is computed against the wrong reference set.  STUB — requires the
frozen public A2A benchmark (WP5).

No private assay data: only public benchmark file roles are swapped.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from fault_injection.base import FaultInjector, FaultSpec, InjectionRecord

_PARAM_GRID: list[dict[str, Any]] = [
    {"swap_mode": mode} for mode in ["full_swap", "partial_50pct", "partial_25pct"]
]


class F12EvalSetMismatch(FaultInjector):
    """Swap actives/decoys roles in the benchmark evaluation (STUB)."""

    @property
    def spec(self) -> FaultSpec:
        return FaultSpec(
            family_id="F12",
            family_name="EVAL_SET_MISMATCH",
            gate="G5",
            severity="critical",
            description=(
                "Swap or partially swap actives and decoys in the benchmark "
                "evaluation set.  STUB — requires frozen public benchmark (WP5)."
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
            "F12 requires the frozen public A2A benchmark actives/decoys files.  "
            "Implement after benchmark freeze at WP5."
        )


F12EvalSetMismatch().verify_no_assay_leakage()

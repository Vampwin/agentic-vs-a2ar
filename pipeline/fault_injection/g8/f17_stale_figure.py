"""
F17 STALE_FIGURE — G8 fault, severity: medium

Injection: Replace the current figure file with a version from a prior run
directory (or a synthetic older version with an earlier timestamp and different
content hash).  G8 reporting contract should detect that the figure's hash does
not match the registered artifact hash for the current run.

Parameterisation (10 instances):
  staleness_hours: [1, 6, 12, 24, 48, 72, 168, 336, 720, 2160]
    — how many hours "old" the stale figure appears to be (via mtime manipulation)
  target_figure: "main_figure"  (single figure name; real path resolved at runtime)

No assay data: only file content and modification time are changed.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any

from fault_injection.base import FaultInjector, FaultSpec, InjectionRecord

_STALENESS_HOURS = [1, 6, 12, 24, 48, 72, 168, 336, 720, 2160]
_PARAM_GRID: list[dict[str, Any]] = [
    {"staleness_hours": h, "target_figure": "main_figure"}
    for h in _STALENESS_HOURS
]


class F17StaleFigure(FaultInjector):
    """Replace the current figure with a stale version from a prior run."""

    @property
    def spec(self) -> FaultSpec:
        return FaultSpec(
            family_id="F17",
            family_name="STALE_FIGURE",
            gate="G8",
            severity="medium",
            description=(
                "Replace the registered figure file with an older version "
                "(different content hash, earlier mtime).  "
                "Expected detection: G8 contract SHA-256 mismatch between "
                "registered artifact hash and current file."
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
        # Look for any figure file in workspace
        figure_candidates = (
            list(workspace.glob("*.png")) +
            list(workspace.glob("*.pdf")) +
            list(workspace.glob("figures/*.png")) +
            list(workspace.glob("figures/*.pdf"))
        )
        if not figure_candidates:
            # Create a synthetic figure placeholder for testing
            fig_path = workspace / "main_figure.png"
            fig_path.write_bytes(b"\x89PNG\r\n\x1a\nCURRENT_VERSION_CONTENT")
        else:
            fig_path = sorted(figure_candidates)[0]

        current_content = fig_path.read_bytes()
        staleness_hours = params["staleness_hours"]

        # Write "stale" content (different bytes → different hash)
        stale_content = b"STALE_FIGURE_CONTENT_FROM_PRIOR_RUN_" + str(staleness_hours).encode()
        fig_path.write_bytes(stale_content)

        # Set mtime to simulate staleness
        stale_mtime = time.time() - staleness_hours * 3600
        import os
        os.utime(fig_path, (stale_mtime, stale_mtime))

        return InjectionRecord(
            instance_id=self.instance_id(instance_index),
            family_id="F17",
            family_name="STALE_FIGURE",
            gate="G8",
            severity="medium",
            params=params,
            affected_items=[str(fig_path.name)],
            ground_truth=(
                f"Replaced {fig_path.name} with stale content "
                f"({staleness_hours}h old mtime).  "
                f"Original size: {len(current_content)} bytes; "
                f"stale size: {len(stale_content)} bytes."
            ),
            containment_criterion=(
                "G8 contract: figure SHA-256 must match registered artifact hash "
                "for the current run"
            ),
        )


F17StaleFigure().verify_no_assay_leakage()

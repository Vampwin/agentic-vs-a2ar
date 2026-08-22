"""
base.py — Common interface for fault injection scripts.

Every fault family implements FaultInjector.  The manifest registers all families
and generates per-instance parameter grids for the feasibility pilot and definitive
benchmark runs.

Design constraints (ADR-AIIM-020, ADR-AIIM-023):
  - No cAMP values, assay results, or known-active IDs may enter any injection script.
  - Each instance is fully reproducible from (family_id, instance_index, seed).
  - Ground truth is recorded per instance for the independent evaluator.
  - unsafe_release label is set by the independent benchmark scorer, not the injector.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass
class FaultSpec:
    """Static description of one fault family (approved at Gate 2)."""
    family_id: str       # "F01" … "F18"
    family_name: str     # "ID_DUPLICATE" etc.
    gate: str            # "G0" … "G8"
    severity: str        # "critical" | "high" | "medium"
    description: str
    param_grid: list[dict[str, Any]]  # one dict per instance (10–30 entries)


@dataclass
class InjectionRecord:
    """
    Per-instance record written after injection.

    The independent evaluator reads this to score whether the fault was detected.
    The acting agent must not read held-out InjectionRecords before evaluation.
    """
    instance_id: str              # "F01_001", "F01_002", …
    family_id: str
    family_name: str
    gate: str
    severity: str
    params: dict[str, Any]        # parameters actually used
    affected_items: list[str]     # compound IDs, file paths, or record keys changed
    ground_truth: str             # human-readable description of what was injected
    containment_criterion: str    # which contract/check should catch this
    # unsafe_release is assigned by the independent benchmark scorer after evaluation
    unsafe_release: int | None = None
    extra: dict[str, Any] = field(default_factory=dict)


class FaultInjector(ABC):
    """Abstract base class for all fault families."""

    @property
    @abstractmethod
    def spec(self) -> FaultSpec:
        """Return the static FaultSpec for this family."""

    @abstractmethod
    def inject(
        self,
        workspace: Path,
        params: dict[str, Any],
        seed: int,
        instance_index: int,
    ) -> InjectionRecord:
        """
        Apply the fault to a copy of the workflow workspace.

        Parameters
        ----------
        workspace:
            Path to the isolated run directory for this instance.
            The injector may read and write files inside workspace.
            Must not read any file outside workspace that contains assay data.
        params:
            Parameter dict from spec.param_grid[instance_index].
        seed:
            Random seed for any non-deterministic step.  Pass numpy RNG as
            np.random.default_rng(seed) inside the implementation.
        instance_index:
            0-based index into spec.param_grid; used to form the instance_id.

        Returns
        -------
        InjectionRecord with ground_truth and affected_items populated.
        unsafe_release is left None; the independent evaluator sets it.
        """

    def instance_id(self, index: int) -> str:
        """Return canonical instance ID, e.g. 'F03_007'."""
        return f"{self.spec.family_id}_{index + 1:03d}"

    @abstractmethod
    def verify_no_assay_leakage(self) -> None:
        """
        Raise AssertionError if any parameter in param_grid references
        assay data (cAMP, known-active IDs, assay-derived thresholds).

        Called at import time for each injector.  Prevents accidental leakage
        during development (ADR-AIIM-020, Gate 2 human confirmation).
        """

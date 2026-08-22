"""
F18 VALUE_MISMATCH — G8 fault, severity: medium

Injection: Change one numeric value in the manuscript results table (e.g. a SEER
percentage, a Brier score, or a compound count) by a small but detectable amount.
G8 reporting contract should detect that the table value does not match the
registered result artifact.

Parameterisation (15 instances):
  delta_fraction: [0.01, 0.05, 0.10, 0.20, 0.50]  — fractional change to value
  sign: [+1, -1, +1]  — direction of change (3 sign values recycled across 5 deltas)
  target_column: ["seer_pct", "brier_score", "n_compounds"]  — 5 combinations
  → 15 instances

No assay data: only computed result values (SEER %, Brier, compound counts) changed.
These are workflow-output metrics, not assay measurements.
"""

from __future__ import annotations

import csv
import io
from pathlib import Path
from typing import Any

import numpy as np

from fault_injection.base import FaultInjector, FaultSpec, InjectionRecord

_DELTAS = [0.01, 0.05, 0.10, 0.20, 0.50]
_SIGNS = [1, -1, 1, -1, 1]  # alternating per delta
_COLS = ["seer_pct", "brier_score", "n_compounds"]
_PARAM_GRID: list[dict[str, Any]] = [
    {"delta_fraction": d, "sign": s, "target_column": col}
    for (d, s), col in zip(zip(_DELTAS, _SIGNS), _COLS * 5)
][:15]


class F18ValueMismatch(FaultInjector):
    """Alter one numeric value in the manuscript results table."""

    @property
    def spec(self) -> FaultSpec:
        return FaultSpec(
            family_id="F18",
            family_name="VALUE_MISMATCH",
            gate="G8",
            severity="medium",
            description=(
                "Change one numeric value in the manuscript results table by a "
                "fractional amount, so the table no longer matches the registered "
                "result artifact.  Expected detection: G8 contract value-match "
                "check between manuscript table and frozen result file."
            ),
            param_grid=_PARAM_GRID,
        )

    def verify_no_assay_leakage(self) -> None:
        for p in _PARAM_GRID:
            assert "camp" not in str(p).lower()
            assert "assay" not in str(p).lower()
            assert "active" not in str(p).lower()

    def inject(
        self,
        workspace: Path,
        params: dict[str, Any],
        seed: int,
        instance_index: int,
    ) -> InjectionRecord:
        table_path = workspace / "results_table.csv"
        if not table_path.exists():
            # Create a synthetic results table for testing
            table_path.write_text(
                "condition,seer_pct,brier_score,n_compounds\n"
                "A_hybrid,0.15,0.12,172\n"
                "B_contracts_only,0.35,0.28,172\n"
                "C_monitor_only,0.42,0.31,172\n"
                "D_neither,0.58,0.45,172\n",
                encoding="utf-8",
            )

        rng = np.random.default_rng(seed)
        delta_frac = params["delta_fraction"]
        sign = params["sign"]
        target_col = params["target_column"]

        with open(table_path, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            fieldnames = list(reader.fieldnames or [])

        # Find actual column (flexible match)
        actual_col = next(
            (c for c in fieldnames if target_col.lower() in c.lower()), None
        ) or next(
            (c for c in fieldnames if any(k in c.lower() for k in ["seer", "brier", "n_comp"])),
            None
        )
        if actual_col is None:
            actual_col = fieldnames[1] if len(fieldnames) > 1 else fieldnames[0]

        row_idx = int(rng.integers(0, len(rows)))
        original_val_str = rows[row_idx].get(actual_col, "0")
        try:
            original_val = float(original_val_str)
        except ValueError:
            original_val = 0.0

        delta = sign * delta_frac * abs(original_val) if original_val != 0 else sign * delta_frac
        new_val = original_val + delta
        rows[row_idx][actual_col] = f"{new_val:.6g}"

        buf = io.StringIO()
        writer = csv.DictWriter(buf, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
        table_path.write_text(buf.getvalue(), encoding="utf-8")

        condition = rows[row_idx].get("condition", str(row_idx))
        return InjectionRecord(
            instance_id=self.instance_id(instance_index),
            family_id="F18",
            family_name="VALUE_MISMATCH",
            gate="G8",
            severity="medium",
            params=params,
            affected_items=[f"{condition}/{actual_col}"],
            ground_truth=(
                f"Changed '{actual_col}' for condition '{condition}': "
                f"{original_val:.6g} → {new_val:.6g} "
                f"(delta_fraction={delta_frac}, sign={sign:+d})"
            ),
            containment_criterion=(
                "G8 contract: every numeric value in the manuscript results table "
                "must match the registered frozen result artifact to within tolerance"
            ),
        )


F18ValueMismatch().verify_no_assay_leakage()

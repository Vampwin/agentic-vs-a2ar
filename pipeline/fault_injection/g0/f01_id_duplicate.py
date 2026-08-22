"""
F01 ID_DUPLICATE — G0 fault, severity: medium

Injection: Duplicate one row in the compound input CSV, giving two rows with the
same compound ID.  G0 contract checks for duplicate IDs and should catch this.

Parameterisation (10 instances): duplicate_fraction in {0.01, 0.02, …, 0.10},
producing a different number of duplicated rows each time.  The specific rows to
duplicate are chosen deterministically from the seed.

No assay data used: only the CSV row index and compound ID column are touched.
"""

from __future__ import annotations

import csv
import io
from pathlib import Path
from typing import Any

import numpy as np

from fault_injection.base import FaultInjector, FaultSpec, InjectionRecord

_PARAM_GRID: list[dict[str, Any]] = [
    {"duplicate_fraction": round(0.01 * (i + 1), 2)} for i in range(10)
]


class F01IdDuplicate(FaultInjector):
    """Duplicate rows in the compound input CSV."""

    @property
    def spec(self) -> FaultSpec:
        return FaultSpec(
            family_id="F01",
            family_name="ID_DUPLICATE",
            gate="G0",
            severity="medium",
            description=(
                "Duplicate a fraction of rows in the compound input CSV so that "
                "multiple rows share the same compound identifier.  "
                "Expected detection: G0 contract duplicate-ID check."
            ),
            param_grid=_PARAM_GRID,
        )

    def verify_no_assay_leakage(self) -> None:
        for p in _PARAM_GRID:
            assert "camp" not in str(p).lower(), "cAMP reference in F01 params"
            assert "assay" not in str(p).lower(), "assay reference in F01 params"

    def inject(
        self,
        workspace: Path,
        params: dict[str, Any],
        seed: int,
        instance_index: int,
    ) -> InjectionRecord:
        csv_path = workspace / "compounds.csv"
        if not csv_path.exists():
            raise FileNotFoundError(f"F01 expects {csv_path}")

        rng = np.random.default_rng(seed)
        fraction = params["duplicate_fraction"]

        with open(csv_path, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            fieldnames = reader.fieldnames or []

        n_to_dup = max(1, int(len(rows) * fraction))
        indices = rng.choice(len(rows), size=n_to_dup, replace=False).tolist()
        duplicated_ids = [rows[i].get("compound_id", str(i)) for i in indices]

        # Insert duplicated rows immediately after their original
        new_rows = list(rows)
        for offset, idx in enumerate(sorted(indices)):
            new_rows.insert(idx + 1 + offset, dict(rows[idx]))

        buf = io.StringIO()
        writer = csv.DictWriter(buf, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(new_rows)
        csv_path.write_text(buf.getvalue(), encoding="utf-8")

        return InjectionRecord(
            instance_id=self.instance_id(instance_index),
            family_id="F01",
            family_name="ID_DUPLICATE",
            gate="G0",
            severity="medium",
            params=params,
            affected_items=duplicated_ids,
            ground_truth=(
                f"Duplicated {n_to_dup} row(s) "
                f"(fraction={fraction}) in compounds.csv: {duplicated_ids}"
            ),
            containment_criterion="G0 contract duplicate-compound-ID check",
        )


# Verify at import time
F01IdDuplicate().verify_no_assay_leakage()

"""
F09 SCORE_DIRECTION — G4 fault, severity: critical

Injection: Negate all Vina docking scores so that more-negative (better) scores
become more-positive (worse) and vice versa.  The ranking is completely inverted.
G4 score-direction contract (comparison against reference scores) should catch this.

Parameterisation (15 instances): n_compounds_affected × score_column_name.
  n_compounds_affected: [10, 20, 30, 50, None (all)]  — 5 values
  score_column_name: ["vina_score", "docking_score", "score"]  — 3 values
  → 15 instances

No assay data used: only numeric score values are negated.
"""

from __future__ import annotations

import csv
import io
from pathlib import Path
from typing import Any

import numpy as np

from fault_injection.base import FaultInjector, FaultSpec, InjectionRecord

_N_AFFECTED = [10, 20, 30, 50, None]
_SCORE_COLS = ["vina_score", "docking_score", "score"]
_PARAM_GRID: list[dict[str, Any]] = [
    {"n_affected": n, "score_column": col}
    for n in _N_AFFECTED
    for col in _SCORE_COLS
]


class F09ScoreDirection(FaultInjector):
    """Negate docking scores to invert the ranking direction."""

    @property
    def spec(self) -> FaultSpec:
        return FaultSpec(
            family_id="F09",
            family_name="SCORE_DIRECTION",
            gate="G4",
            severity="critical",
            description=(
                "Negate all (or a subset of) Vina docking scores so the ranking "
                "order is fully inverted.  Expected detection: G4 contract "
                "shadow score-direction check (correlation with reference scores)."
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
        score_path = workspace / "docking_scores.csv"
        if not score_path.exists():
            raise FileNotFoundError(f"F09 expects {score_path}")

        rng = np.random.default_rng(seed)
        target_col = params["score_column"]
        n_affected = params["n_affected"]

        with open(score_path, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            fieldnames = list(reader.fieldnames or [])

        # Find actual score column (flexible match)
        actual_col = next(
            (c for c in fieldnames if target_col.lower() in c.lower()), None
        ) or next(
            (c for c in fieldnames if "score" in c.lower()), None
        )
        if actual_col is None:
            raise ValueError(f"F09: no score column matching '{target_col}' found")

        id_col = next((c for c in fieldnames if "id" in c.lower()), fieldnames[0])

        numeric_indices = [
            i for i, r in enumerate(rows)
            if r.get(actual_col, "").strip() not in ("", "NA", "None")
        ]

        if n_affected is None:
            selected = numeric_indices
        else:
            selected = rng.choice(
                numeric_indices, size=min(n_affected, len(numeric_indices)), replace=False
            ).tolist()

        affected_ids = []
        for idx in selected:
            try:
                original = float(rows[idx][actual_col])
                rows[idx][actual_col] = str(-original)
                affected_ids.append(rows[idx].get(id_col, str(idx)))
            except ValueError:
                continue

        buf = io.StringIO()
        writer = csv.DictWriter(buf, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
        score_path.write_text(buf.getvalue(), encoding="utf-8")

        return InjectionRecord(
            instance_id=self.instance_id(instance_index),
            family_id="F09",
            family_name="SCORE_DIRECTION",
            gate="G4",
            severity="critical",
            params=params,
            affected_items=affected_ids[:20],
            ground_truth=(
                f"Negated '{actual_col}' for {len(affected_ids)} compounds "
                f"(n_affected={n_affected}); ranking fully inverted for affected rows."
            ),
            containment_criterion=(
                "G4 contract: shadow_score_rank_rho must be positive "
                "(Spearman ρ ≥ threshold with reference ranking)"
            ),
        )


F09ScoreDirection().verify_no_assay_leakage()

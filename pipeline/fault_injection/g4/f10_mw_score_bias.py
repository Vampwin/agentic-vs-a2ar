"""
F10 MW_SCORE_BIAS — G4 fault, severity: high

Injection: Replace docking scores with molecular-weight ranks, so larger molecules
are ranked higher regardless of docking quality.  This simulates a silent scoring
bias that conflates molecular weight with affinity.

MW is computed from the SMILES column using RDKit — no assay data required.

Parameterisation (12 instances):
  mw_direction: ["descending", "ascending"]  (larger MW = rank 1, or rank N)
  n_replaced: [20, 40, 60, 80, 100, None(all)]
  → 12 instances

No assay data: only molecular weight (a structural descriptor) and score ranks used.
"""

from __future__ import annotations

import csv
import io
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors

from fault_injection.base import FaultInjector, FaultSpec, InjectionRecord

_DIRECTIONS = ["descending", "ascending"]
_N_REPLACED = [20, 40, 60, 80, 100, None]
_PARAM_GRID: list[dict[str, Any]] = [
    {"mw_direction": d, "n_replaced": n}
    for n in _N_REPLACED
    for d in _DIRECTIONS
]


class F10MwScoreBias(FaultInjector):
    """Replace docking scores with molecular-weight ranks."""

    @property
    def spec(self) -> FaultSpec:
        return FaultSpec(
            family_id="F10",
            family_name="MW_SCORE_BIAS",
            gate="G4",
            severity="high",
            description=(
                "Replace docking scores with molecular-weight ranks, silently "
                "biasing selection toward larger (or smaller) molecules.  "
                "Expected detection: G4 contract shadow_score_rank_rho drop "
                "vs reference scores."
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
            raise FileNotFoundError(f"F10 expects {score_path}")

        rng = np.random.default_rng(seed)

        with open(score_path, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            fieldnames = list(reader.fieldnames or [])

        score_col = next((c for c in fieldnames if "score" in c.lower()), None)
        smiles_col = next((c for c in fieldnames if "smiles" in c.lower()), None)
        id_col = next((c for c in fieldnames if "id" in c.lower()), fieldnames[0])

        if score_col is None:
            raise ValueError("F10: no score column found")
        if smiles_col is None:
            raise ValueError("F10: no SMILES column found — needed for MW calc")

        n_replaced = params["n_replaced"]
        direction = params["mw_direction"]

        numeric_indices = [
            i for i, r in enumerate(rows)
            if r.get(score_col, "").strip() not in ("", "NA", "None")
            and r.get(smiles_col, "").strip()
        ]

        if n_replaced is None:
            selected = numeric_indices
        else:
            selected = rng.choice(
                numeric_indices, size=min(n_replaced, len(numeric_indices)), replace=False
            ).tolist()

        # Compute MW for selected rows using RDKit
        mws = []
        for idx in selected:
            mol = Chem.MolFromSmiles(rows[idx][smiles_col])
            mw = Descriptors.ExactMolWt(mol) if mol else 0.0
            mws.append((idx, mw))

        # Rank by MW
        reverse = direction == "descending"
        mws_sorted = sorted(mws, key=lambda x: x[1], reverse=reverse)
        mw_ranks = {idx: rank + 1 for rank, (idx, _) in enumerate(mws_sorted)}

        affected_ids = []
        for idx in selected:
            rows[idx][score_col] = str(-float(mw_ranks[idx]))  # negative = "better" for Vina
            affected_ids.append(rows[idx].get(id_col, str(idx)))

        buf = io.StringIO()
        writer = csv.DictWriter(buf, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
        score_path.write_text(buf.getvalue(), encoding="utf-8")

        return InjectionRecord(
            instance_id=self.instance_id(instance_index),
            family_id="F10",
            family_name="MW_SCORE_BIAS",
            gate="G4",
            severity="high",
            params=params,
            affected_items=affected_ids[:20],
            ground_truth=(
                f"Replaced '{score_col}' with MW-rank scores ({direction} MW = rank 1) "
                f"for {len(selected)} compounds."
            ),
            containment_criterion=(
                "G4 contract: shadow_score_rank_rho between MW-replaced scores and "
                "reference docking scores should fall below acceptable threshold"
            ),
        )


F10MwScoreBias().verify_no_assay_leakage()

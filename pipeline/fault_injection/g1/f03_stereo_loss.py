"""
F03 STEREO_LOSS — G1 fault, severity: high

Injection: Remove '@' and '@@' stereocentre tokens from a SMILES string, causing
the standardised structure to lose chirality information.  G1 stereochemistry
contract should flag the mismatch between input and standardised InChI.

Parameterisation (12 instances): stereo_fraction (fraction of chiral SMILES to
strip, from 0.25 to 1.0 in steps of ~0.07).  Molecule selection is deterministic
from seed.  SMILES manipulation uses RDKit (rdkit skill v1.1).

No assay data used: only SMILES tokens are modified.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem

from fault_injection.base import FaultInjector, FaultSpec, InjectionRecord

_PARAM_GRID: list[dict[str, Any]] = [
    {"stereo_fraction": round(0.25 + 0.068 * i, 3)} for i in range(12)
]

_STEREO_RE = re.compile(r"@{1,2}")


def _strip_stereo(smi: str) -> str:
    """Remove all @ and @@ tokens from a SMILES string."""
    return _STEREO_RE.sub("", smi)


class F03StereoLoss(FaultInjector):
    """Strip stereocentre tokens from chiral SMILES entries."""

    @property
    def spec(self) -> FaultSpec:
        return FaultSpec(
            family_id="F03",
            family_name="STEREO_LOSS",
            gate="G1",
            severity="high",
            description=(
                "Remove '@'/'@@' stereocentre tokens from a fraction of chiral "
                "SMILES in the standardised ledger.  Expected detection: G1 "
                "contract stereocentre-count comparison between input and "
                "standardised structure."
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
        import csv, io

        ledger_path = workspace / "standardised_ledger.csv"
        if not ledger_path.exists():
            raise FileNotFoundError(f"F03 expects {ledger_path}")

        rng = np.random.default_rng(seed)
        fraction = params["stereo_fraction"]

        with open(ledger_path, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            fieldnames = list(reader.fieldnames or [])

        smiles_col = next(
            (c for c in fieldnames if "smiles" in c.lower()), None
        )
        id_col = next(
            (c for c in fieldnames if "id" in c.lower()), fieldnames[0]
        )
        if smiles_col is None:
            raise ValueError("F03: no SMILES column found in ledger")

        chiral_indices = [
            i for i, r in enumerate(rows)
            if "@" in r.get(smiles_col, "")
        ]
        n_to_strip = max(1, int(len(chiral_indices) * fraction))
        selected = rng.choice(chiral_indices, size=min(n_to_strip, len(chiral_indices)),
                              replace=False).tolist()

        affected_ids = []
        for idx in selected:
            original = rows[idx][smiles_col]
            rows[idx][smiles_col] = _strip_stereo(original)
            affected_ids.append(rows[idx].get(id_col, str(idx)))

        buf = io.StringIO()
        writer = csv.DictWriter(buf, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
        ledger_path.write_text(buf.getvalue(), encoding="utf-8")

        return InjectionRecord(
            instance_id=self.instance_id(instance_index),
            family_id="F03",
            family_name="STEREO_LOSS",
            gate="G1",
            severity="high",
            params=params,
            affected_items=affected_ids,
            ground_truth=(
                f"Stripped '@'/'@@' from {len(selected)} chiral SMILES "
                f"(stereo_fraction={fraction}): {affected_ids[:5]}{'…' if len(affected_ids) > 5 else ''}"
            ),
            containment_criterion=(
                "G1 contract: stereocentre count in standardised InChI must match "
                "stereocentre count in source SMILES"
            ),
        )


F03StereoLoss().verify_no_assay_leakage()

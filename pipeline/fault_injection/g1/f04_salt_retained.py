"""
F04 SALT_RETAINED — G1 fault, severity: medium

Injection: Keep the counter-ion fragment in the SMILES string by writing the
multi-component form (parent.counter_ion) into the standardised ledger instead of
the de-salted parent.  G1 salt-stripping contract should detect that the
standardised structure still contains a disconnected fragment.

Parameterisation (12 instances): salt_type × n_compounds.
  salt_types: [Cl, Na, K, Ca, Br, Mg] (6 common counter-ions).
  n_compounds per type: [5, 10] → 12 instances.

Counter-ion SMILES are standard chemistry — no assay data.
"""

from __future__ import annotations

import csv
import io
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem

from fault_injection.base import FaultInjector, FaultSpec, InjectionRecord

_SALTS = {
    "Cl":  "[Cl-]",
    "Na":  "[Na+]",
    "K":   "[K+]",
    "Ca":  "[Ca+2]",
    "Br":  "[Br-]",
    "Mg":  "[Mg+2]",
}
_PARAM_GRID: list[dict[str, Any]] = [
    {"salt_type": s, "n_compounds": n}
    for s in _SALTS
    for n in [5, 10]
]


class F04SaltRetained(FaultInjector):
    """Append a counter-ion fragment to de-salted SMILES entries."""

    @property
    def spec(self) -> FaultSpec:
        return FaultSpec(
            family_id="F04",
            family_name="SALT_RETAINED",
            gate="G1",
            severity="medium",
            description=(
                "Append a common counter-ion (Cl⁻, Na⁺, etc.) to de-salted SMILES "
                "entries, simulating a salt-stripping step that was silently skipped.  "
                "Expected detection: G1 contract fragment-count check on standardised SMILES."
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
        ledger_path = workspace / "standardised_ledger.csv"
        if not ledger_path.exists():
            raise FileNotFoundError(f"F04 expects {ledger_path}")

        rng = np.random.default_rng(seed)
        salt_smi = _SALTS[params["salt_type"]]
        n = params["n_compounds"]

        with open(ledger_path, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            fieldnames = list(reader.fieldnames or [])

        smiles_col = next((c for c in fieldnames if "smiles" in c.lower()), None)
        id_col = next((c for c in fieldnames if "id" in c.lower()), fieldnames[0])
        if smiles_col is None:
            raise ValueError("F04: no SMILES column in ledger")

        valid_indices = [
            i for i, r in enumerate(rows)
            if r.get(smiles_col, "") and "." not in r[smiles_col]
        ]
        selected = rng.choice(valid_indices, size=min(n, len(valid_indices)),
                              replace=False).tolist()

        affected_ids = []
        for idx in selected:
            rows[idx][smiles_col] = rows[idx][smiles_col] + "." + salt_smi
            affected_ids.append(rows[idx].get(id_col, str(idx)))

        buf = io.StringIO()
        writer = csv.DictWriter(buf, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
        ledger_path.write_text(buf.getvalue(), encoding="utf-8")

        return InjectionRecord(
            instance_id=self.instance_id(instance_index),
            family_id="F04",
            family_name="SALT_RETAINED",
            gate="G1",
            severity="medium",
            params=params,
            affected_items=affected_ids,
            ground_truth=(
                f"Appended '{salt_smi}' counter-ion to {len(selected)} compounds: "
                f"{affected_ids[:5]}{'…' if len(affected_ids) > 5 else ''}"
            ),
            containment_criterion=(
                "G1 contract: standardised SMILES must contain exactly one "
                "connected fragment (no '.' separators)"
            ),
        )


F04SaltRetained().verify_no_assay_leakage()

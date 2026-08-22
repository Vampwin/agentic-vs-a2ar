"""
F05 INVALID_VALENCE — G1 fault, severity: high

Injection: Introduce a pentavalent carbon atom in a SMILES string (e.g. replace a
ring carbon token with 'C(F)(F)(F)(F)(F)' inline or force valence 5 via bracket
notation '[C@@H](F)(F)(F)(F)' — but the simplest reproducible approach is to
swap a '[CH2]' or 'C' token for '[C]5' which forces a valence RDKit cannot
normalise).  In practice: write '[C]' (no implicit H) attached to more bonds than
carbon's maximum valence allows.  RDKit sanitisation or the G1 valence contract
should reject the structure.

Parameterisation (10 instances): target_position (fraction along molecule SMILES
string, 5 values) × replacement_pattern (2 patterns).

No assay data used.
"""

from __future__ import annotations

import csv
import io
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from fault_injection.base import FaultInjector, FaultSpec, InjectionRecord

# Two patterns that reliably create invalid valence without needing assay info:
#   pattern A: append "(F)(F)(F)(F)(F)" to a methyl group → 5-valent carbon
#   pattern B: prefix a heavy-atom carbon with extra explicit bonds via
#              a bracket atom that forces valence error
_PATTERNS = ["append_5F_to_methyl", "bracket_overvalent"]
_POSITIONS = [0.1, 0.3, 0.5, 0.7, 0.9]
_PARAM_GRID: list[dict[str, Any]] = [
    {"target_position": pos, "pattern": pat}
    for pos in _POSITIONS
    for pat in _PATTERNS
]


def _make_overvalent_smiles(smiles: str, position: float, pattern: str) -> str | None:
    """
    Return a SMILES that will fail RDKit valence check, or None if not possible.

    Strategy:
      - Parse with sanitize=False so we can manipulate atoms.
      - Find a carbon at roughly `position` fraction along the atom list.
      - Depending on pattern, either:
          A: add enough F substituents to push valence above 4, or
          B: write bracket atom [C](bonds) that exceeds normal valence.

    Simpler and more reproducible: build the broken SMILES by text substitution
    of the first occurrence of 'C(' or 'C' with a known-broken fragment.
    """
    if pattern == "append_5F_to_methyl":
        # Find 'C' that is followed by non-ring-closure char and replace
        idx = int(len(smiles) * position)
        for i in range(max(0, idx - 5), min(len(smiles), idx + 5)):
            if smiles[i] == "C" and (i + 1 >= len(smiles) or smiles[i + 1] not in "0123456789"):
                # Insert pentavalent carbon fragment
                return smiles[:i] + "C(F)(F)(F)(F)(F)" + smiles[i + 1:]
        return None
    elif pattern == "bracket_overvalent":
        # Replace first 'C' with [CH6] — 6 explicit hydrogens on carbon → valence error
        idx = smiles.find("C")
        if idx == -1:
            return None
        return smiles[:idx] + "[CH6]" + smiles[idx + 1:]
    return None


class F05InvalidValence(FaultInjector):
    """Introduce a chemically impossible valence into SMILES entries."""

    @property
    def spec(self) -> FaultSpec:
        return FaultSpec(
            family_id="F05",
            family_name="INVALID_VALENCE",
            gate="G1",
            severity="high",
            description=(
                "Inject a pentavalent carbon (or similarly overvalent atom) into a "
                "SMILES string so that RDKit sanitisation fails.  "
                "Expected detection: G1 contract RDKit sanitisation check."
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
        ledger_path = workspace / "standardised_ledger.csv"
        if not ledger_path.exists():
            raise FileNotFoundError(f"F05 expects {ledger_path}")

        rng = np.random.default_rng(seed)

        with open(ledger_path, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            fieldnames = list(reader.fieldnames or [])

        # Use the same column-priority order as the G1 evaluator's smiles_roundtrip
        # check so the injected fault lands in the column the evaluator actually reads.
        _SMILES_PRIORITY = ["smiles_canonical", "parent_smiles", "smiles_input", "SMILES", "smiles"]
        smiles_col = next(
            (c for c in _SMILES_PRIORITY if c in fieldnames),
            next((c for c in fieldnames if "smiles" in c.lower()), None),
        )
        id_col = next((c for c in fieldnames if "id" in c.lower()), fieldnames[0])
        if smiles_col is None:
            raise ValueError("F05: no SMILES column in ledger")

        valid_indices = [i for i, r in enumerate(rows) if r.get(smiles_col, "")]
        if not valid_indices:
            raise ValueError("F05: no SMILES entries in ledger")

        target_idx = int(rng.choice(valid_indices))
        original_smi = rows[target_idx][smiles_col]
        broken = _make_overvalent_smiles(
            original_smi, params["target_position"], params["pattern"]
        )
        if broken is None:
            # Fallback: prepend a bracket atom that forces valence error
            broken = "[CH6]" + original_smi

        rows[target_idx][smiles_col] = broken
        affected_id = rows[target_idx].get(id_col, str(target_idx))

        buf = io.StringIO()
        writer = csv.DictWriter(buf, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
        ledger_path.write_text(buf.getvalue(), encoding="utf-8")

        # Confirm the broken SMILES actually fails RDKit sanitisation
        mol = Chem.MolFromSmiles(broken, sanitize=False)
        try:
            if mol is not None:
                Chem.SanitizeMol(mol)
                # Sanitisation unexpectedly passed — record this for the evaluator
                broken_confirmed = False
            else:
                broken_confirmed = True
        except Exception:
            broken_confirmed = True

        return InjectionRecord(
            instance_id=self.instance_id(instance_index),
            family_id="F05",
            family_name="INVALID_VALENCE",
            gate="G1",
            severity="high",
            params=params,
            affected_items=[affected_id],
            ground_truth=(
                f"Replaced SMILES of {affected_id} with invalid-valence form "
                f"(pattern={params['pattern']}); "
                f"RDKit sanitisation fails: {broken_confirmed}"
            ),
            containment_criterion=(
                "G1 contract: all standardised SMILES must pass RDKit full sanitisation"
            ),
        )


F05InvalidValence().verify_no_assay_leakage()

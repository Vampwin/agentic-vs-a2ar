"""
F13 UNDECLARED_EXTERNAL_LABEL_CONSTRAINT — G6 fault, severity: critical

Injection: Silently apply a hidden structural filter during the compound selection
step, biasing the panel toward compounds that satisfy an undeclared medicinal-
chemistry property constraint.  The filter uses only standard Lipinski-style
structural descriptors — NO cAMP values, NO assay results, NO known-active IDs.

This simulates a reviewer silently preferring drug-like compounds without
documenting the preference, so the final selection deviates unexpectedly from
what the declared scoring alone would produce.

Parameterisation (20 instances — 5 filters × 4 threshold variants):

  Structural filters (Lipinski-style, no assay connection):
    F_MW:  MW < threshold      thresholds: [350, 400, 450, 500] Da
    F_ROT: rotatable_bonds ≤ threshold  thresholds: [5, 7, 9, 12]
    F_HBA: HBA ≤ threshold     thresholds: [5, 7, 8, 10]
    F_HBD: HBD ≤ threshold     thresholds: [2, 3, 4, 5]
    F_RING: ring_count ≤ threshold  thresholds: [2, 3, 4, 5]

Each instance applies one filter with one threshold to the selection input.

Scientific-critical-thinking check: MW, rotatable bonds, HBA, HBD, ring count are
ALL computable from SMILES alone — they are structural properties with NO connection
to biological activity, cAMP, or any assay outcome (ADR-AIIM-020, Gate 2 confirmed).
"""

from __future__ import annotations

import csv
import io
from pathlib import Path
from typing import Any

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from fault_injection.base import FaultInjector, FaultSpec, InjectionRecord

# Structural filter definitions — all Lipinski-style, zero assay connection
_FILTERS: list[dict[str, Any]] = [
    {"filter_name": "MW",             "threshold": 350},
    {"filter_name": "MW",             "threshold": 400},
    {"filter_name": "MW",             "threshold": 450},
    {"filter_name": "MW",             "threshold": 500},
    {"filter_name": "rotatable_bonds","threshold": 5},
    {"filter_name": "rotatable_bonds","threshold": 7},
    {"filter_name": "rotatable_bonds","threshold": 9},
    {"filter_name": "rotatable_bonds","threshold": 12},
    {"filter_name": "HBA",            "threshold": 5},
    {"filter_name": "HBA",            "threshold": 7},
    {"filter_name": "HBA",            "threshold": 8},
    {"filter_name": "HBA",            "threshold": 10},
    {"filter_name": "HBD",            "threshold": 2},
    {"filter_name": "HBD",            "threshold": 3},
    {"filter_name": "HBD",            "threshold": 4},
    {"filter_name": "HBD",            "threshold": 5},
    {"filter_name": "ring_count",     "threshold": 2},
    {"filter_name": "ring_count",     "threshold": 3},
    {"filter_name": "ring_count",     "threshold": 4},
    {"filter_name": "ring_count",     "threshold": 5},
]
_PARAM_GRID = _FILTERS


def _compute_descriptor(mol: Any, filter_name: str) -> float:
    if filter_name == "MW":
        return Descriptors.ExactMolWt(mol)
    elif filter_name == "rotatable_bonds":
        return float(rdMolDescriptors.CalcNumRotatableBonds(mol))
    elif filter_name == "HBA":
        return float(rdMolDescriptors.CalcNumHBA(mol))
    elif filter_name == "HBD":
        return float(rdMolDescriptors.CalcNumHBD(mol))
    elif filter_name == "ring_count":
        return float(rdMolDescriptors.CalcNumRings(mol))
    raise ValueError(f"Unknown filter: {filter_name}")


class F13UndeclaredExternalLabelConstraint(FaultInjector):
    """Silently apply a structural filter to compound selection (no assay data)."""

    @property
    def spec(self) -> FaultSpec:
        return FaultSpec(
            family_id="F13",
            family_name="UNDECLARED_EXTERNAL_LABEL_CONSTRAINT",
            gate="G6",
            severity="critical",
            description=(
                "Silently apply an undeclared structural property filter (MW, "
                "rotatable bonds, HBA, HBD, or ring count) during compound "
                "selection, biasing the panel toward drug-like structures.  "
                "Filter uses ONLY standard medicinal-chemistry descriptors — "
                "no cAMP, no assay activity, no known-active IDs (ADR-AIIM-020).  "
                "Expected detection: G6 contract unexplained distribution shift "
                "in selected-panel properties."
            ),
            param_grid=_PARAM_GRID,
        )

    def verify_no_assay_leakage(self) -> None:
        forbidden = {"camp", "assay", "active", "inactive", "ic50", "ec50", "ki"}
        for p in _PARAM_GRID:
            for v in str(p).lower().split():
                assert v not in forbidden, f"F13 param contains forbidden term: {v}"
        # Confirm filter names are purely structural
        for p in _PARAM_GRID:
            assert p["filter_name"] in {
                "MW", "rotatable_bonds", "HBA", "HBD", "ring_count"
            }, f"F13: non-structural filter name: {p['filter_name']}"

    def inject(
        self,
        workspace: Path,
        params: dict[str, Any],
        seed: int,
        instance_index: int,
    ) -> InjectionRecord:
        selection_path = workspace / "ranked_compounds.csv"
        if not selection_path.exists():
            raise FileNotFoundError(f"F13 expects {selection_path}")

        with open(selection_path, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            fieldnames = list(reader.fieldnames or [])

        smiles_col = next((c for c in fieldnames if "smiles" in c.lower()), None)
        id_col = next((c for c in fieldnames if "id" in c.lower()), fieldnames[0])
        if smiles_col is None:
            raise ValueError("F13: no SMILES column found")

        filter_name = params["filter_name"]
        threshold = params["threshold"]

        # Compute descriptor for each row; filter to those satisfying the constraint
        passing_rows, removed_rows = [], []
        for row in rows:
            smi = row.get(smiles_col, "")
            mol = Chem.MolFromSmiles(smi) if smi else None
            if mol is None:
                passing_rows.append(row)  # keep unparseable rows (conservative)
                continue
            val = _compute_descriptor(mol, filter_name)
            if val <= threshold:
                passing_rows.append(row)
            else:
                removed_rows.append(row)

        removed_ids = [r.get(id_col, "?") for r in removed_rows]

        buf = io.StringIO()
        writer = csv.DictWriter(buf, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(passing_rows)
        selection_path.write_text(buf.getvalue(), encoding="utf-8")

        return InjectionRecord(
            instance_id=self.instance_id(instance_index),
            family_id="F13",
            family_name="UNDECLARED_EXTERNAL_LABEL_CONSTRAINT",
            gate="G6",
            severity="critical",
            params=params,
            affected_items=removed_ids[:20],
            ground_truth=(
                f"Applied undeclared structural filter: {filter_name} ≤ {threshold}.  "
                f"Removed {len(removed_rows)} compounds from selection input.  "
                f"Filter uses ONLY RDKit structural descriptors — no assay data."
            ),
            containment_criterion=(
                "G6 contract: unexplained shift in selected-panel property distribution "
                "(MW, rotatable bonds, etc.) versus unfiltered ranking"
            ),
        )


F13UndeclaredExternalLabelConstraint().verify_no_assay_leakage()

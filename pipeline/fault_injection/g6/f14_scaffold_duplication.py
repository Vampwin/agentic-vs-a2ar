"""
F14 SCAFFOLD_DUPLICATION — G6 fault, severity: medium

Injection: Select duplicate Murcko scaffolds into the final compound panel,
violating the structural diversity constraint.  Implemented by replacing diverse
panel members with additional compounds that share a Murcko scaffold already
represented in the panel.

Parameterisation (12 instances):
  n_duplicates: [2, 3, 4, 5, 6, 8] — number of extra scaffold-redundant compounds added
  strategy: ["top_scored", "random"]  — which scaffold to duplicate
  → 12 instances

No assay data: Murcko scaffold computed from SMILES via RDKit.
"""

from __future__ import annotations

import csv
import io
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

from fault_injection.base import FaultInjector, FaultSpec, InjectionRecord

_N_DUPS = [2, 3, 4, 5, 6, 8]
_STRATEGIES = ["top_scored", "random"]
_PARAM_GRID: list[dict[str, Any]] = [
    {"n_duplicates": n, "strategy": s}
    for n in _N_DUPS
    for s in _STRATEGIES
]


def _get_murcko(smi: str) -> str | None:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    try:
        core = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(core)
    except Exception:
        return None


class F14ScaffoldDuplication(FaultInjector):
    """Add scaffold-redundant compounds to the selection panel."""

    @property
    def spec(self) -> FaultSpec:
        return FaultSpec(
            family_id="F14",
            family_name="SCAFFOLD_DUPLICATION",
            gate="G6",
            severity="medium",
            description=(
                "Replace diverse panel members with extra compounds sharing an "
                "already-represented Murcko scaffold, violating the diversity rule.  "
                "Expected detection: G6 contract Murcko scaffold diversity check."
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
        panel_path = workspace / "selected_panel.csv"
        ranked_path = workspace / "ranked_compounds.csv"
        if not panel_path.exists():
            raise FileNotFoundError(f"F14 expects {panel_path}")

        rng = np.random.default_rng(seed)
        n_dup = params["n_duplicates"]
        strategy = params["strategy"]

        with open(panel_path, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            panel_rows = list(reader)
            fieldnames = list(reader.fieldnames or [])

        smiles_col = next((c for c in fieldnames if "smiles" in c.lower()), None)
        id_col = next((c for c in fieldnames if "id" in c.lower()), fieldnames[0])
        if smiles_col is None:
            raise ValueError("F14: no SMILES column in panel")

        # Compute scaffold for each panel member
        scaffold_map: dict[str, list[int]] = {}
        for i, row in enumerate(panel_rows):
            sc = _get_murcko(row.get(smiles_col, ""))
            if sc:
                scaffold_map.setdefault(sc, []).append(i)

        # Pick a scaffold that appears once to duplicate
        single_scaffolds = {sc: idxs for sc, idxs in scaffold_map.items()
                            if len(idxs) == 1}
        if not single_scaffolds:
            raise ValueError("F14: no unique scaffolds to duplicate in panel")

        if strategy == "top_scored":
            target_scaffold = next(iter(single_scaffolds))
        else:
            target_scaffold = rng.choice(list(single_scaffolds.keys()))

        original_idx = single_scaffolds[target_scaffold][0]
        original_row = panel_rows[original_idx]

        # Duplicate the row n_dup times with modified IDs
        added_ids = []
        for k in range(n_dup):
            new_row = dict(original_row)
            dup_id = f"{original_row.get(id_col, 'X')}_dup{k+1}"
            new_row[id_col] = dup_id
            panel_rows.append(new_row)
            added_ids.append(dup_id)

        buf = io.StringIO()
        writer = csv.DictWriter(buf, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(panel_rows)
        panel_path.write_text(buf.getvalue(), encoding="utf-8")

        return InjectionRecord(
            instance_id=self.instance_id(instance_index),
            family_id="F14",
            family_name="SCAFFOLD_DUPLICATION",
            gate="G6",
            severity="medium",
            params=params,
            affected_items=added_ids,
            ground_truth=(
                f"Added {n_dup} scaffold-redundant copies of panel member "
                f"'{original_row.get(id_col)}' (scaffold: {target_scaffold[:40]}…)"
            ),
            containment_criterion=(
                "G6 contract: Murcko scaffold diversity check — each scaffold "
                "may appear at most once in the final panel"
            ),
        )


F14ScaffoldDuplication().verify_no_assay_leakage()

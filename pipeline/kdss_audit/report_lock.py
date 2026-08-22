"""
Report lock — prevents generation of figures, tables, or manuscript snippets
from artifacts whose required contracts did not pass (or lack adjudication).

The acting agent may not bypass the lock. The lock is released only when
ContractEngine has recorded a PASS for the required artifact.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any

from .artifact_registry import ArtifactRegistry


class ReportLock:
    def __init__(self, registry: ArtifactRegistry, lock_file: Path) -> None:
        self._registry = registry
        self._lock_file = lock_file

    def can_report(self, artifact_id: str) -> bool:
        """True only if the artifact has a PASS contract record or human adjudication."""
        if self._registry.artifact_has_passed_contract(artifact_id):
            return True
        # Check for human adjudication record.
        # ArtifactRegistry._append flattens the entry into the record's top
        # level, so artifact_ids is read from the record directly — reading it
        # from a nested "data" key silently matches nothing and would leave
        # every adjudicated artifact permanently blocked.
        for rec in self._registry.get_all_records():
            if (rec["record_type"] == "human_intervention"
                    and artifact_id in rec.get("artifact_ids", [])):
                return True
        return False

    def assert_can_report(self, artifact_ids: list[str]) -> None:
        blocked = [aid for aid in artifact_ids if not self.can_report(aid)]
        if blocked:
            raise PermissionError(
                f"ReportLock: cannot release report — the following artifacts "
                f"have no passed contract or adjudication record: {blocked}"
            )

    def lock_report(
        self,
        figures: list[dict[str, str]],
        tables: list[dict[str, str]],
        result_snippets: list[dict[str, Any]],
        run_id: str,
    ) -> Path:
        """
        Freeze the report manifest. Raises PermissionError if any source artifact
        is unregistered or lacks a passed contract.

        figures: list of {"figure_id": ..., "artifact_id": ..., "caption": ...}
        tables:  list of {"table_id":  ..., "artifact_id": ..., "caption": ...}
        result_snippets: list of {"claim": ..., "value": ..., "artifact_id": ...}
        """
        all_artifact_ids = (
            [f["artifact_id"] for f in figures]
            + [t["artifact_id"] for t in tables]
            + [s["artifact_id"] for s in result_snippets]
        )
        self.assert_can_report(all_artifact_ids)

        manifest = {
            "run_id": run_id,
            "locked_at": time.time(),
            "human_intervention_count": self._registry.count_human_interventions(),
            "figures": figures,
            "tables": tables,
            "result_snippets": result_snippets,
        }
        self._lock_file.parent.mkdir(parents=True, exist_ok=True)
        with open(self._lock_file, "w", encoding="utf-8") as f:
            json.dump(manifest, f, indent=2, ensure_ascii=False)
        return self._lock_file

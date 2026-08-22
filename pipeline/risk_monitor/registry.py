"""
registry.py — Trusted model artifact registry (ADR-AIIM-021).

Only a model whose SHA-256 hash matches a registered entry may be loaded.
Missing hash, hash mismatch, or unregistered path all fail closed.

The registry is a JSON file at a path chosen by the caller.  Each entry records:
  - artifact_id         : human-readable label
  - model_path          : absolute path to the serialised sklearn Pipeline (joblib)
  - sha256              : expected uppercase hex SHA-256 of the file
  - schema_version      : telemetry schema version at training time
  - feature_order       : canonical feature name list used at training time
  - feature_order_sha256: SHA-256 of the JSON-serialised feature_order list
  - schema_sha256       : SHA-256 of the feature_schema.yaml file at registration time
  - registered_at       : ISO-8601 UTC timestamp
  - approved_by         : human identifier (required; non-empty string)
  - notes               : optional string

The registry is append-only: re-registering an existing artifact_id is an error.
The acting agent cannot write new registry entries; only the model developer or
custodian may do so before unsealing the held-out partition.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


class RegistryError(RuntimeError):
    """Raised when a model artifact cannot be trusted."""


@dataclass
class RegistryEntry:
    artifact_id: str
    model_path: str
    sha256: str
    schema_version: str
    registered_at: str
    approved_by: str
    notes: str = ""
    feature_order: list[str] = field(default_factory=list)
    feature_order_sha256: str = ""
    schema_sha256: str = ""

    def as_dict(self) -> dict[str, Any]:
        return {
            "artifact_id": self.artifact_id,
            "model_path": self.model_path,
            "sha256": self.sha256,
            "schema_version": self.schema_version,
            "registered_at": self.registered_at,
            "approved_by": self.approved_by,
            "notes": self.notes,
            "feature_order": self.feature_order,
            "feature_order_sha256": self.feature_order_sha256,
            "schema_sha256": self.schema_sha256,
        }


def _sha256_file(path: Path) -> str:
    if not path.is_file():
        raise RegistryError(f"Model file not found: {path}")
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest().upper()


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest().upper()


class ModelRegistry:
    """Persistent, append-only registry of trusted model artifacts."""

    def __init__(self, registry_path: Path | str) -> None:
        self._path = Path(registry_path)
        self._entries: dict[str, RegistryEntry] = {}
        if self._path.is_file():
            self._load()

    def _load(self) -> None:
        with open(self._path, encoding="utf-8") as fh:
            data = json.load(fh)
        for e in data.get("entries", []):
            # Backward compatibility: older entries may lack new fields
            e.setdefault("feature_order", [])
            e.setdefault("feature_order_sha256", "")
            e.setdefault("schema_sha256", "")
            self._entries[e["artifact_id"]] = RegistryEntry(**e)

    def _save(self) -> None:
        self._path.parent.mkdir(parents=True, exist_ok=True)
        payload = {"entries": [e.as_dict() for e in self._entries.values()]}
        with open(self._path, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2)

    def register(
        self,
        artifact_id: str,
        model_path: Path | str,
        schema_version: str,
        approved_by: str,
        notes: str = "",
        feature_order: list[str] | None = None,
        schema_file_path: Path | str | None = None,
    ) -> RegistryEntry:
        """
        Compute the file hash and register the model artifact.

        Parameters
        ----------
        artifact_id:
            Unique human-readable key.  Registering the same ID twice is an error
            (append-only; prevents silent overwrites).
        model_path:
            Path to the serialised sklearn Pipeline (.joblib).
        schema_version:
            Telemetry schema_version at training time (from feature_schema.yaml).
        approved_by:
            Non-empty human identifier; required for all registrations.
        feature_order:
            Ordered list of feature names used during training.  When provided,
            the SHA-256 of the JSON-serialised list is also stored.
        schema_file_path:
            Path to feature_schema.yaml.  When provided, its SHA-256 is stored.

        Raises RegistryError if the file is missing, approved_by is empty, or the
        artifact_id already exists.
        """
        if not approved_by or not approved_by.strip():
            raise RegistryError(
                f"approved_by must be a non-empty string for artifact '{artifact_id}'. "
                "Only a named human may approve model registration."
            )
        if artifact_id in self._entries:
            raise RegistryError(
                f"Artifact '{artifact_id}' is already registered. "
                "The registry is append-only; choose a new artifact_id for each model."
            )

        p = Path(model_path)
        sha = _sha256_file(p)

        # Feature order hash
        fo = feature_order or []
        fo_sha = _sha256_bytes(json.dumps(fo, sort_keys=False).encode()) if fo else ""

        # Schema file hash
        schema_sha = ""
        if schema_file_path is not None:
            sp = Path(schema_file_path)
            if sp.is_file():
                schema_sha = _sha256_bytes(sp.read_bytes())

        entry = RegistryEntry(
            artifact_id=artifact_id,
            model_path=str(p.resolve()),
            sha256=sha,
            schema_version=schema_version,
            registered_at=datetime.now(timezone.utc).isoformat(),
            approved_by=approved_by,
            notes=notes,
            feature_order=fo,
            feature_order_sha256=fo_sha,
            schema_sha256=schema_sha,
        )
        self._entries[artifact_id] = entry
        self._save()
        return entry

    def verify_and_get_path(self, artifact_id: str) -> Path:
        """
        Return the trusted model path only if it exists and its hash matches.
        Raises RegistryError on any failure — always fails closed.
        """
        if artifact_id not in self._entries:
            raise RegistryError(
                f"Model artifact '{artifact_id}' is not in the registry. "
                "Only registered, hash-verified artifacts may be loaded."
            )
        entry = self._entries[artifact_id]
        p = Path(entry.model_path)
        actual_sha = _sha256_file(p)
        if actual_sha != entry.sha256.upper():
            raise RegistryError(
                f"Hash mismatch for artifact '{artifact_id}': "
                f"expected {entry.sha256!r}, got {actual_sha!r}. "
                "Model file may have been tampered with or replaced."
            )
        return p

    def get_entry(self, artifact_id: str) -> RegistryEntry:
        if artifact_id not in self._entries:
            raise RegistryError(f"Artifact '{artifact_id}' not found in registry")
        return self._entries[artifact_id]

    def list_artifact_ids(self) -> list[str]:
        return list(self._entries.keys())

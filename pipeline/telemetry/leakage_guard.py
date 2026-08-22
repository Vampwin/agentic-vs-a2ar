"""
leakage_guard.py — Forbidden-feature policy enforcement (ADR-AIIM-020).

Loads forbidden patterns and the allowed-name set from feature_schema.yaml at import
time.  Call validate_telemetry() on field names and validate_record() on full records
before they can reach the model.  Raises LeakageError on any violation so the caller
must deal with it explicitly — the guard never silently drops a field.
"""

from __future__ import annotations

import math
import re
import yaml
from pathlib import Path
from typing import Any, Iterable

SCHEMA_PATH = Path(__file__).parent / "feature_schema.yaml"


class LeakageError(ValueError):
    """Raised when a forbidden or unknown field is detected in telemetry."""


def _compile_forbidden(schema: dict) -> list[re.Pattern]:
    return [
        re.compile(p, re.IGNORECASE)
        for p in schema.get("forbidden_patterns", [])
    ]


def _collect_feature_info(schema: dict) -> tuple[frozenset[str], list[str], dict[str, dict]]:
    """
    Returns (allowed_names frozenset, canonical_order list, feature_metadata dict).

    canonical_order preserves schema YAML declaration order across all families.
    feature_metadata maps name → full feature dict (dtype, nullable, range_min, etc.).
    """
    names: list[str] = []
    seen: set[str] = set()
    metadata: dict[str, dict] = {}
    for family in schema.get("feature_families", {}).values():
        for feat in family.get("features", []):
            name = feat["name"]
            names.append(name)
            seen.add(name)
            metadata[name] = feat
    return frozenset(seen), names, metadata


def _load_schema() -> dict:
    with open(SCHEMA_PATH, encoding="utf-8") as fh:
        return yaml.safe_load(fh)


_SCHEMA: dict = _load_schema()
_FORBIDDEN_PATTERNS: list[re.Pattern] = _compile_forbidden(_SCHEMA)
_ALLOWED_NAMES, _CANONICAL_ORDER, _FEATURE_METADATA = _collect_feature_info(_SCHEMA)


# ── Name-only guards ──────────────────────────────────────────────────────────

def check_field_names(fields: Iterable[str], context: str = "") -> None:
    """Raise LeakageError if any field name matches a forbidden pattern."""
    for field in fields:
        for pat in _FORBIDDEN_PATTERNS:
            if pat.search(field):
                raise LeakageError(
                    f"Forbidden field {field!r} in {context!r}: "
                    f"matches pattern '{pat.pattern}' (ADR-AIIM-020)"
                )


def check_allowlist(fields: Iterable[str], context: str = "") -> None:
    """Raise LeakageError if any field name is not in the schema allowlist."""
    unknown = [f for f in fields if f not in _ALLOWED_NAMES]
    if unknown:
        raise LeakageError(
            f"Unknown fields {unknown!r} in {context!r}: "
            f"not in feature_schema.yaml allowlist (ADR-AIIM-020)"
        )


def validate_telemetry(record: dict, context: str = "telemetry") -> None:
    """Name-only validation: forbidden patterns then allowlist.  Call before any model input."""
    check_field_names(record.keys(), context)
    check_allowlist(record.keys(), context)


# ── Full record validation ────────────────────────────────────────────────────

def _check_dtype(name: str, value: Any, expected: str, context: str) -> None:
    """Raise LeakageError if value does not match expected dtype string."""
    if expected == "bool":
        if not isinstance(value, bool):
            raise LeakageError(
                f"Field {name!r} in {context!r}: expected bool, "
                f"got {type(value).__name__} = {value!r}"
            )
    elif expected == "int":
        # bool is a subclass of int in Python — reject it for int fields
        if isinstance(value, bool) or not isinstance(value, int):
            raise LeakageError(
                f"Field {name!r} in {context!r}: expected int, "
                f"got {type(value).__name__} = {value!r}"
            )
    elif expected == "float":
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise LeakageError(
                f"Field {name!r} in {context!r}: expected numeric (int or float), "
                f"got {type(value).__name__} = {value!r}"
            )
        # Non-finite values are always an error; use None for missing (nullable) fields
        if not math.isfinite(float(value)):
            raise LeakageError(
                f"Field {name!r} in {context!r}: non-finite value {value!r}. "
                "Use None for missing values on nullable fields."
            )


def validate_record(record: dict, context: str = "validate_record") -> None:
    """
    Full record validation: field names, key completeness, dtype, range, finiteness.

    Checks enforced in order:
      1. Forbidden-pattern check on all field names.
      2. Allowlist membership check on all field names.
      3. Key completeness: every schema field must be present (nullable may be None).
      4. Per-field: nullable contract (non-nullable field must not be None).
      5. Per-field: dtype check (int, float, bool).
      6. Per-field: range check (range_min, range_max from schema).

    Raises LeakageError on the first violation found.
    """
    # 1 & 2: name checks
    check_field_names(record.keys(), context)
    check_allowlist(record.keys(), context)

    # 3: key completeness
    missing_keys = _ALLOWED_NAMES - set(record.keys())
    if missing_keys:
        raise LeakageError(
            f"Record in {context!r} is missing required schema fields: "
            f"{sorted(missing_keys)}. "
            "All 24 schema fields must be present (nullable fields may be None)."
        )

    # 4–6: per-field validation
    for name in _CANONICAL_ORDER:
        meta = _FEATURE_METADATA[name]
        value = record[name]
        nullable = meta.get("nullable", False)

        # 4: nullable contract
        if value is None:
            if not nullable:
                raise LeakageError(
                    f"Field {name!r} in {context!r} is None but nullable=false."
                )
            continue  # nullable None passes all further checks

        # 5: dtype
        expected_dtype = meta.get("dtype")
        if expected_dtype:
            _check_dtype(name, value, expected_dtype, context)

        # 6: range
        range_min = meta.get("range_min")
        range_max = meta.get("range_max")
        numeric_value = float(value) if not isinstance(value, bool) else None
        if numeric_value is not None:
            if range_min is not None and numeric_value < float(range_min):
                raise LeakageError(
                    f"Field {name!r} in {context!r}: value {value!r} "
                    f"is below range_min={range_min!r}."
                )
            if range_max is not None and numeric_value > float(range_max):
                raise LeakageError(
                    f"Field {name!r} in {context!r}: value {value!r} "
                    f"exceeds range_max={range_max!r}."
                )


# ── Schema accessors ──────────────────────────────────────────────────────────

def get_allowed_names() -> frozenset[str]:
    """Return the frozen set of permitted feature names from the schema."""
    return _ALLOWED_NAMES


def get_canonical_feature_order() -> list[str]:
    """Return feature names in schema declaration order (families then within each family)."""
    return list(_CANONICAL_ORDER)


def get_schema_version() -> str:
    return _SCHEMA.get("schema_version", "unknown")

"""
telemetry — Common pre-release telemetry extraction and leakage enforcement.

All four factorial conditions A–D use the same extractor; only enforcement
and model invocation differ between conditions.
"""

from .leakage_guard import (
    LeakageError,
    check_field_names,
    check_allowlist,
    validate_telemetry,
    get_allowed_names,
    get_schema_version,
)
from .feature_extractor import extract

__all__ = [
    "LeakageError",
    "check_field_names",
    "check_allowlist",
    "validate_telemetry",
    "get_allowed_names",
    "get_schema_version",
    "extract",
]

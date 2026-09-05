"""
Canonical assay-blindness blocklist (ADR-AIIM-004).

Single source of truth for the column names that must never appear in any
pipeline input, intermediate, or output. The policy is enforced at several
gates — G0 (inputs), G4 (scoring outputs), G6 (selection) — and those gates
must agree exactly: a column that is blocked at one gate but permitted at
another is a leak path, and the paper's claim of mechanical enforcement rests
on there being one definition rather than several copies that can drift.

Do not redefine these lists in an evaluator. Import from here.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable

import pandas as pd

# ── Scope ────────────────────────────────────────────────────────────────────
#
# Strictness increases as data flows toward selection, because that is where a
# leaked label actually changes a scientific conclusion:
#
#   source library   in-house assay markers forbidden. Literature activity
#                    annotations (IC50, "Biological activities") are permitted
#                    here — the SNP-Lib source is a curated literature library
#                    and those columns are part of it — but flagged, because
#                    they must not survive standardization.
#
#   pipeline data    ledger, scores, ranking, selected panel: NO activity
#                    column of any provenance. These files feed ranking, so a
#                    literature IC50 riding along is as damaging as an in-house
#                    one.
#
#   public benchmark exempt from column-absence checks. A benchmark without
#                    labels is not a benchmark; its labels are guarded at read
#                    time by G5a label isolation, not by pretending they are
#                    absent.

# In-house assay markers (ADR-AIIM-004). Forbidden everywhere in the pipeline.
BLOCKED_EXACT: frozenset[str] = frozenset({
    "known_active", "a_sim_actives", "pct_inhibition",
    "cAMP", "assay_label", "is_active",
})
BLOCKED_SUBSTRINGS: tuple[str, ...] = ("camp", "assay", "inhibit")

# Activity annotations of any provenance, including public literature values.
# Forbidden in pipeline data; permitted-but-flagged in the source library.
ACTIVITY_SUBSTRINGS: tuple[str, ...] = (
    "ic50", "ec50", "pchembl", "potency", "bioactiv",
    "biological activit", "activity", "active",
)


def is_blocked_column(column: str, strict: bool = False) -> bool:
    """
    True if this column violates the assay-blindness policy.

    strict=False checks in-house assay markers only (source-library scope).
    strict=True additionally rejects any activity annotation (pipeline scope).
    """
    col_lower = column.lower()
    if col_lower in {b.lower() for b in BLOCKED_EXACT}:
        return True
    if any(sub in col_lower for sub in BLOCKED_SUBSTRINGS):
        return True
    if strict and any(sub in col_lower for sub in ACTIVITY_SUBSTRINGS):
        return True
    return False


def blocked_columns_in(columns: Iterable[str], strict: bool = False) -> list[str]:
    """Return every blocked name in an iterable of column names."""
    return [c for c in columns if is_blocked_column(c, strict=strict)]


def activity_annotations_in(columns: Iterable[str]) -> list[str]:
    """
    Activity annotations that are tolerated in the source library but must not
    propagate downstream. Reported so the source's contents are on the record.
    """
    out = []
    for c in columns:
        if is_blocked_column(c, strict=False):
            continue  # already a hard violation; reported separately
        if any(sub in c.lower() for sub in ACTIVITY_SUBSTRINGS):
            out.append(c)
    return out


def scan_csv_header(csv_path: Path | str) -> list[str]:
    """
    Return blocked column names in a CSV's header.

    A file that cannot be read returns no findings — callers that treat this
    check as critical must verify the file exists separately, because
    "unreadable" is not the same as "clean".
    """
    path = Path(csv_path) if csv_path else None
    if not path or not path.is_file():
        return []
    try:
        header = pd.read_csv(path, nrows=0)
    except Exception:
        return []
    return blocked_columns_in(header.columns)


def scan_dataframe(df: "pd.DataFrame") -> list[str]:
    """Return blocked column names present in an in-memory dataframe."""
    return blocked_columns_in(df.columns)

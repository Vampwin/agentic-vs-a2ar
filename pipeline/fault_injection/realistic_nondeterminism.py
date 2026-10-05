"""realistic_nondeterminism.py -- run-to-run variation that changes no science.

The seeded-defect harness builds every workspace from a fixed seed
(`_REAL_DATA_SEED`), so calling a workspace factory twice yields byte-identical
files. ADR-AIIM-194 recorded the consequence: a golden-file comparator detects
128 of 128 seeded defects, because a byte-exact diff against a re-built
reference catches any injected change by construction. The benchmark could not
separate a detector that models faults from one that compares bytes.

This module removes that property, by making the workspace factories behave the
way the production pipeline behaves.

WHAT IS INJECTED, AND WHY EACH ONE IS NEUTRAL
---------------------------------------------
Every perturbation here must satisfy two conditions: it must occur in the real
pipeline, and it must leave every scientific fact in the file unchanged. A
perturbation that changes a value, an ordering that carries meaning, or a
relationship between files is a *fault*, not noise, and belongs in the fault
registry rather than here.

  row_order   Row order in an unordered table. Real: parallel writes, dict
              iteration order, database result order, multiprocessing pool
              completion order. Neutral: the table is a set of records.
              **Applied only to files whose order carries no meaning.**
              `ranked_compounds.csv`, `selected_panel.csv` and
              `results_table.csv` are ORDERED OUTPUTS -- their row order is a
              scientific claim -- and are never shuffled. Shuffling them would
              be fault F-like behaviour that contracts SHOULD catch.

  float_fmt   Decimal representation of a float column. Real: numpy/pandas
              version differences, locale, writer settings, round-tripping
              through JSON. Neutral: only ever widened, never narrowed, so no
              significant digit is lost. -9.571 may be written -9.571000 but
              never -9.57.

  quoting     CSV quoting style. Real: different writers (pandas, csv, R,
              awk). Neutral: the parsed values are identical.

  mtime       File modification time. Real: always. Neutral: no content.

WHAT IS NOT INJECTED
--------------------
Nothing that alters a value, a row's membership, a column's name, a
cross-file relationship, or the order of an ordered output. Those are the
things the contracts exist to detect, and putting them here would rig the
comparison in the contracts' favour rather than testing it.

Pre-registered in ADR-AIIM-201 before any arm was run.
"""
from __future__ import annotations

import csv
import os
import time
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd

# Files whose row order is a scientific claim. Never shuffled.
ORDERED_OUTPUTS: frozenset[str] = frozenset({
    "ranked_compounds.csv",   # ranked by docking score
    "selected_panel.csv",     # the top-N selection, in order
    "results_table.csv",      # the reported table
})

# Float formats. Every one preserves at least the three decimals the harness
# writes, so no significant digit is ever lost.
_FLOAT_FORMATS: tuple[Optional[str], ...] = (None, "%.4f", "%.6f", "%.8f")

_QUOTING = (csv.QUOTE_MINIMAL, csv.QUOTE_NONNUMERIC)


def _perturb_csv(path: Path, rng: np.random.Generator) -> dict[str, Any]:
    """Rewrite one CSV with neutral run-to-run variation. Returns what was done."""
    applied: dict[str, Any] = {"row_order": False, "float_fmt": None,
                               "quoting": "minimal"}
    try:
        df = pd.read_csv(path)
    except Exception:
        return applied
    if df.empty:
        return applied

    if path.name not in ORDERED_OUTPUTS and len(df) > 1:
        perm = rng.permutation(len(df))
        # A shuffle that happens to be the identity would leave the file
        # byte-identical and silently weaken the experiment, so reject it.
        if np.array_equal(perm, np.arange(len(df))) and len(df) > 1:
            perm = np.roll(perm, 1)
        df = df.iloc[perm].reset_index(drop=True)
        applied["row_order"] = True

    fmt = _FLOAT_FORMATS[int(rng.integers(len(_FLOAT_FORMATS)))]
    quot = _QUOTING[int(rng.integers(len(_QUOTING)))]
    applied["float_fmt"] = fmt or "default"
    applied["quoting"] = "nonnumeric" if quot == csv.QUOTE_NONNUMERIC else "minimal"

    df.to_csv(path, index=False, float_format=fmt, quoting=quot)
    return applied


def _jitter_mtime(path: Path, rng: np.random.Generator) -> None:
    """Move a file's mtime by up to an hour. Content is untouched."""
    delta = float(rng.uniform(-3600, 3600))
    t = max(1.0, time.time() + delta)
    os.utime(path, (t, t))


def perturb_workspace(ws: Path, rng: np.random.Generator) -> dict[str, Any]:
    """Apply neutral run-to-run variation to every file in a built workspace.

    Called once after a workspace factory runs. Two calls with different `rng`
    states produce workspaces that are semantically identical and almost
    certainly not byte-identical -- which is the condition the production
    pipeline is in, and the condition the benchmark previously lacked.
    """
    record: dict[str, Any] = {"files": {}}
    for p in sorted(ws.rglob("*")):
        if not p.is_file():
            continue
        if p.suffix.lower() == ".csv":
            record["files"][p.name] = _perturb_csv(p, rng)
        else:
            record["files"][p.name] = {"row_order": False, "float_fmt": None,
                                       "quoting": None}
        _jitter_mtime(p, rng)
    return record


def make_setup_map(base_setup: dict[str, Any], rng: np.random.Generator
                   ) -> dict[str, Any]:
    """Wrap each workspace factory so its output carries run-to-run variation.

    The returned map is a drop-in replacement for `pilot_runner._WORKSPACE_SETUP`
    and for the `setup_map` argument of
    `conventional_test_baseline.build_reference`, so neither the harness nor
    the comparator needs to know this module exists.
    """
    def wrap(fn):
        def inner(ws: Path, *args, **kwargs):
            out = fn(ws, *args, **kwargs)
            perturb_workspace(Path(ws), rng)
            return out
        return inner

    return {gate: wrap(fn) for gate, fn in base_setup.items()}


def semantic_equal(a: Path, b: Path) -> tuple[bool, str]:
    """Do two workspaces carry the same science, ignoring the neutral noise?

    Used to verify the injector before the benchmark runs: if this ever returns
    False for two clean workspaces, the perturbation is not neutral and the
    experiment is invalid.
    """
    fa = {p.name for p in a.rglob("*") if p.is_file()}
    fb = {p.name for p in b.rglob("*") if p.is_file()}
    if fa != fb:
        return False, f"file sets differ: {fa ^ fb}"
    for name in sorted(fa):
        pa, pb = a / name, b / name
        if not name.lower().endswith(".csv"):
            if pa.read_bytes() != pb.read_bytes():
                return False, f"{name}: non-CSV bytes differ"
            continue
        da, db = pd.read_csv(pa), pd.read_csv(pb)
        if sorted(da.columns) != sorted(db.columns):
            return False, f"{name}: columns differ"
        db = db[da.columns]
        if len(da) != len(db):
            return False, f"{name}: row count {len(da)} vs {len(db)}"
        if name in ORDERED_OUTPUTS:
            ka, kb = da, db
        else:
            by = [c for c in da.columns
                  if da[c].dtype == object or "id" in c.lower()][:2] or list(da.columns[:1])
            ka = da.sort_values(by).reset_index(drop=True)
            kb = db.sort_values(by).reset_index(drop=True)
        for c in da.columns:
            sa, sb = ka[c], kb[c]
            if pd.api.types.is_numeric_dtype(sa) and pd.api.types.is_numeric_dtype(sb):
                if not np.allclose(sa.fillna(-9e9), sb.fillna(-9e9), rtol=0, atol=5e-4):
                    return False, f"{name}: column {c} differs numerically"
            else:
                if not sa.astype(str).equals(sb.astype(str)):
                    return False, f"{name}: column {c} differs"
    return True, "semantically identical"

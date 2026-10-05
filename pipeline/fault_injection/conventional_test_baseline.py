"""
conventional_test_baseline.py — a conventional software-testing comparator.

What a competent engineering team writes WITHOUT a contract framework:
snapshot (golden-file) comparison, schema assertions, and referential
integrity between a gate's own outputs. Scope pre-registered in ADR-AIIM-193
BEFORE this was run against any fault family, and not adjusted afterwards.

Explicitly excluded, by design — each of these IS the contract framework
rather than ordinary practice:
  - any comparison against an independently registered, hash-anchored
    reference file;
  - statistical calibration, correlation or bias detection;
  - timestamp, file-order or provenance forensics.

Returns the same PASS/HALT vocabulary as `sanity_check_baseline.py`, so the
pilot harness's scoring applies unchanged. There is no ESCALATE: a team
without a formal framework has no human-adjudication channel.

Findings are tagged with the practice that produced them (`snapshot:`,
`schema:`, `refint:`) so the result can report how much of any detection is
the snapshot comparison alone. ADR-AIIM-193 records why that matters.
"""
from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any

import pandas as pd

# Files each gate emits. The reference is built by the harness's own workspace
# factory, so this list is the factories' output, not a guess.
GATE_FILES: dict[str, list[str]] = {
    "G0": ["compounds.csv"],
    "G1": ["standardised_ledger.csv"],
    "G4": ["docking_scores.csv"],
    "G6": ["ranked_compounds.csv", "selected_panel.csv"],
    "G8": ["results_table.csv", "main_figure.png"],
}


def _sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _read(p: Path):
    if not p.is_file() or p.stat().st_size == 0:
        return None
    try:
        return pd.read_csv(p)
    except Exception:
        return None


def evaluate(gate: str, ws: Path, golden: Path) -> dict[str, Any]:
    """Run the conventional-test baseline against a workspace.

    `golden` is a directory holding the reference outputs for this gate,
    produced by the same workspace factory from a known-good run.
    """
    findings: list[str] = []
    files = GATE_FILES.get(gate, [])

    # -- practice 1: snapshot / golden-file comparison -----------------------
    for fn in files:
        got, ref = ws / fn, golden / fn
        if not ref.is_file():
            continue                      # nothing to compare against
        if not got.is_file():
            findings.append(f"snapshot: {fn} missing from output")
            continue
        if _sha(got) != _sha(ref):
            findings.append(f"snapshot: {fn} differs from the reference output")

    # -- practice 2: schema assertions ---------------------------------------
    for fn in [f for f in files if f.endswith(".csv")]:
        got_df, ref_df = _read(ws / fn), _read(golden / fn)
        if ref_df is None:
            continue
        if got_df is None:
            findings.append(f"schema: {fn} missing or unparseable")
            continue
        if list(got_df.columns) != list(ref_df.columns):
            findings.append(
                f"schema: {fn} columns {list(got_df.columns)} != "
                f"reference {list(ref_df.columns)}")
            continue
        for col in ref_df.columns:
            gk, rk = got_df[col].dtype.kind, ref_df[col].dtype.kind
            # numeric kinds are interchangeable (int vs float); a numeric
            # column turning into text is not.
            num = set("iuf")
            if (gk in num) != (rk in num):
                findings.append(
                    f"schema: {fn}[{col}] dtype kind {gk!r} != reference {rk!r}")

    # -- practice 3: referential integrity between a gate's own files --------
    if gate == "G6":
        ranked, panel = _read(ws / "ranked_compounds.csv"), _read(ws / "selected_panel.csv")
        if ranked is not None and panel is not None:
            idc = "compound_id"
            if idc in ranked.columns and idc in panel.columns:
                orphans = set(panel[idc]) - set(ranked[idc])
                if orphans:
                    findings.append(
                        f"refint: {len(orphans)} selected_panel id(s) absent from "
                        f"ranked_compounds (e.g. {sorted(orphans)[0]})")
            if idc in panel.columns:
                dup = int(panel[idc].duplicated().sum())
                if dup:
                    findings.append(f"refint: selected_panel has {dup} duplicate {idc}")

    return {"state": "HALT" if findings else "PASS", "findings": findings}


def build_reference(gate: str, dest: Path, setup_map: dict) -> Path:
    """Produce the reference output for a gate from a known-good run."""
    dest.mkdir(parents=True, exist_ok=True)
    setup_map[gate](dest)
    return dest


# ---------------------------------------------------------------------------
# Tolerant comparison (ADR-AIIM-202 addendum)
# ---------------------------------------------------------------------------

def _canonical(path: Path, sort_rows: bool, decimals: int = 3):
    """Parse a CSV into a comparable canonical form.

    Normalises the two things a naive byte diff trips over under realistic
    non-determinism: row order and float formatting.
    """
    import pandas as _pd
    df = _pd.read_csv(path)
    for c in df.columns:
        if _pd.api.types.is_float_dtype(df[c]):
            df[c] = df[c].round(decimals)
    df = df[sorted(df.columns)]
    if sort_rows and len(df) > 1:
        df = df.sort_values(list(df.columns), kind="mergesort").reset_index(drop=True)
    return df


def evaluate_tolerant(gate: str, ws: Path, golden: Path,
                      sort_all: bool = True) -> dict[str, Any]:
    """The strong-man version of golden-file comparison.

    A competent team facing a non-deterministic pipeline would not diff bytes;
    it would normalise first. This does that: rows sorted, floats rounded to
    the precision the pipeline reports, column order ignored.

    `sort_all` is the knob that matters. Sorting every table removes
    false escalations caused by benign row-order variation, and simultaneously
    blinds the comparator to any fault whose signature *is* an ordering
    change. Knowing which files may be sorted and which may not is exactly the
    per-file semantic knowledge a contract encodes, which is the point
    ADR-AIIM-202 tests.
    """
    import pandas as _pd
    findings: list[str] = []
    for fname in GATE_FILES.get(gate, []):
        a, b = ws / fname, golden / fname
        if not b.exists():
            continue
        if not a.exists():
            findings.append(f"tolerant:{fname} missing")
            continue
        if not fname.lower().endswith(".csv"):
            if _sha(a) != _sha(b):
                findings.append(f"tolerant:{fname} bytes differ")
            continue
        try:
            da = _canonical(a, sort_all)
            db = _canonical(b, sort_all)
        except Exception as exc:                      # noqa: BLE001
            findings.append(f"tolerant:{fname} unreadable ({exc})")
            continue
        if list(da.columns) != list(db.columns):
            findings.append(f"tolerant:{fname} columns differ")
            continue
        if len(da) != len(db):
            findings.append(f"tolerant:{fname} row count {len(da)} vs {len(db)}")
            continue
        if not da.equals(db):
            diff = [c for c in da.columns if not da[c].equals(db[c])]
            findings.append(f"tolerant:{fname} values differ in {','.join(diff[:4])}")
    return {"state": "HALT" if findings else "PASS", "findings": findings}

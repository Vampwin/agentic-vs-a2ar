"""
sanity_check_baseline.py — deterministic "D-prime" naive-analyst baseline.

Represents a competent-but-unaided analyst without a frozen contract
framework: existence/parse checks, loose row-count plausibility, basic
type/range sanity (including RDKit parseability, which is ordinary
cheminformatics hygiene, not fault-specific), and duplicate-ID checks.

Explicitly excluded, by design (this is the line between "naive" and
"audited"):
  - No comparison against any independently-registered/hash-anchored
    reference file.
  - No timestamp or file-creation-order forensics.
  - No statistical calibration (e.g. correlation/bias detection).
  - No scaffold-diversity or provenance-column scanning (these duplicate
    what the real G6 checks do; including them would not be "naive").

These checks and their scope were fixed BEFORE being run against any fault
family, and are not adjusted after seeing which faults they do or do not
catch (see the design-review discussion preceding this file).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd


def _try_read_csv(path: Path) -> "pd.DataFrame | None":
    if not path.exists() or path.stat().st_size == 0:
        return None
    try:
        return pd.read_csv(path)
    except Exception:
        return None


def _smiles_cols(df: pd.DataFrame) -> list[str]:
    """All SMILES-like columns, not just the first: a naive but thorough
    analyst checking 'do the structures parse' would not arbitrarily stop
    at whichever SMILES-named column happens to come first."""
    return [c for c in df.columns if "smiles" in c.lower()]


def _id_col(df: pd.DataFrame) -> str | None:
    for c in df.columns:
        cl = c.lower()
        if cl in ("compound_id", "snp_id", "id"):
            return c
    for c in df.columns:
        if "id" in c.lower():
            return c
    return None


def _check_smiles_parse(df: pd.DataFrame, col: str) -> list[str]:
    """Return list of findings for SMILES that fail RDKit parse or full sanitisation."""
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    findings = []
    for _, row in df.iterrows():
        smi = row.get(col)
        if not isinstance(smi, str) or not smi:
            continue
        mol = Chem.MolFromSmiles(smi, sanitize=False)
        if mol is None:
            findings.append(f"unparseable SMILES: {smi[:40]}")
            continue
        try:
            Chem.SanitizeMol(mol)
        except Exception as exc:
            findings.append(f"sanitisation failed ({exc.__class__.__name__}): {smi[:40]}")
    return findings


# Expected clean row count per gate's MAIN file, from the workspace-setup
# functions' own known design (G0/G1 sample _N_COMPOUNDS=50; G4 samples up
# to 50 real scored compounds; G6's ranked_compounds.csv is the full
# candidate list, also 50; G8's results_table.csv has exactly one row per
# condition, always 4, by construction -- not a fault, a fixed table shape).
# A naive analyst calibrates their own expectation to what a clean run of
# that specific workspace produces, not one blind global number.
_EXPECTED_N: dict[str, int] = {"G0": 50, "G1": 50, "G4": 50, "G6": 50, "G8": 4}


def evaluate(gate: str, ws: Path, expected_n: "int | None" = None) -> dict[str, Any]:
    """
    Run the naive sanity-check baseline against a workspace.

    Returns {"state": "PASS"|"HALT", "findings": [...]}. Mirrors the real
    evaluators' PASS/HALT vocabulary so the same escaped/detected scoring
    logic in the pilot harness applies unchanged; this baseline never
    returns ESCALATE (a naive analyst without a formal framework has no
    human-adjudication channel to escalate to).
    """
    if expected_n is None:
        expected_n = _EXPECTED_N.get(gate, 50)
    findings: list[str] = []

    gate_files = {
        "G0": [("compounds.csv", True)],
        "G1": [("standardised_ledger.csv", True)],
        "G4": [("docking_scores.csv", True)],
        "G6": [("ranked_compounds.csv", True), ("selected_panel.csv", False)],
        "G8": [("results_table.csv", True)],
    }.get(gate, [])

    dfs: dict[str, pd.DataFrame] = {}
    for filename, required in gate_files:
        df = _try_read_csv(ws / filename)
        if df is None:
            if required:
                findings.append(f"{filename}: missing or unparseable")
            continue
        dfs[filename] = df

    if not dfs:
        return {"state": "HALT" if findings else "PASS", "findings": findings}

    main_filename = gate_files[0][0]
    main_df = dfs.get(main_filename)
    if main_df is None:
        return {"state": "HALT", "findings": findings}

    # Check 2: row-count plausibility (+/- 5% of a loosely-remembered expectation)
    n = len(main_df)
    if expected_n > 0 and abs(n - expected_n) / expected_n > 0.05:
        findings.append(f"{main_filename}: row count {n} far from expected ~{expected_n}")

    # Check 3a: duplicate-ID hygiene (ordinary habit, not a formal provenance check)
    idc = _id_col(main_df)
    if idc is not None:
        n_dup = int(main_df[idc].duplicated().sum())
        if n_dup > 0:
            findings.append(f"{main_filename}: {n_dup} duplicate values in ID column {idc!r}")

    # Check 3b: SMILES parseability (ordinary cheminformatics hygiene) --
    # every SMILES-like column, not just the first one found.
    for smc in _smiles_cols(main_df):
        bad = _check_smiles_parse(main_df, smc)
        if bad:
            findings.append(f"{main_filename}[{smc}]: {len(bad)} SMILES failed RDKit "
                             f"parse/sanitisation (e.g. {bad[0]})")

    # Check 3c: numeric range/type sanity, gate-specific but generic domain knowledge
    if gate == "G4" and "vina_score" in main_df.columns:
        scores = pd.to_numeric(main_df["vina_score"], errors="coerce")
        if scores.isna().any():
            findings.append("docking_scores.csv: non-numeric or NaN vina_score values")
        if (scores.dropna() > 0).any():
            findings.append("docking_scores.csv: positive Vina score(s) present "
                             "(Vina scores are conventionally negative)")
        if (scores.dropna() < -15).any() or (scores.dropna() > 0).any():
            findings.append("docking_scores.csv: vina_score outside plausible range [-15, 0]")

    if gate == "G8":
        for col in ("seer_pct", "brier_score"):
            if col in main_df.columns:
                vals = pd.to_numeric(main_df[col], errors="coerce")
                if vals.isna().any() or (vals.dropna() < 0).any() or (vals.dropna() > 1).any():
                    findings.append(f"results_table.csv: {col} outside [0, 1]")

    if gate == "G6":
        for fname, df in dfs.items():
            for col in ("MW", "HBA", "HBD", "ring_count", "rotatable_bonds"):
                if col in df.columns:
                    vals = pd.to_numeric(df[col], errors="coerce")
                    if vals.isna().any() or (vals.dropna() < 0).any():
                        findings.append(f"{fname}: {col} has negative or non-numeric values")

    state = "HALT" if findings else "PASS"
    return {"state": state, "findings": findings}

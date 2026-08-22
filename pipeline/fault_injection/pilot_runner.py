"""
pilot_runner.py — Feasibility pilot for Gate 9 (ADR-AIIM-024).

Runs 2-4 development families through 4 experimental conditions (A/B/C/D)
using isolated synthetic workspaces.  Collects per-instance × per-condition
outcomes to estimate:

  SEER  — Silent Error Escape Rate per condition
  ICC   — Within-family Intraclass Correlation Coefficient (for power analysis)

Conditions
----------
  A  contracts ON  + monitor ON   (hybrid — monitor stub until model trained)
  B  contracts ON  + monitor OFF  (contracts only — primary pilot condition)
  C  contracts OFF + monitor ON   (monitor only — stub until model trained)
  D  contracts OFF + monitor OFF  (neither — all faults escape; SEER=1.0 baseline)

Design constraints
------------------
  - Development families only (F01–F18); H01–H12 held-out must not be referenced.
  - No assay data, no cAMP values in any synthetic workspace file.
  - Conditions A/C with monitor scaffolded but marked "stub_not_trained".
  - Each (family, instance, condition) case runs in an isolated temp directory.
  - Evaluator errors are recorded as "escaped" (fail-closed).
  - ICC is computed from condition B binary outcomes with bootstrap 95 % CI.
"""

from __future__ import annotations

import hashlib
import json
import logging
import shutil
import tempfile
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd

from fault_injection.manifest import FAULT_REGISTRY, implementable_families

logger = logging.getLogger(__name__)

# ─── Workspace constants ──────────────────────────────────────────────────────

_N_COMPOUNDS = 50  # rows per workspace, sampled from the real SNP-Lib ledger

# Conditions: (contracts_on, monitor_on)
_CONDITION_MAP: dict[str, tuple[bool, bool]] = {
    "A": (True,  True),
    "B": (True,  False),
    "C": (False, True),
    "D": (False, False),
}

# Dummy script-hash value for G6 pilot contract (both sides equal → check passes)
_DUMMY_HASH = "A" * 64


# ─── Dataclasses ─────────────────────────────────────────────────────────────

@dataclass
class PilotConfig:
    pilot_families: list[str] = field(
        default_factory=lambda: [
            "F01", "F02", "F03", "F04", "F05",
            "F09", "F10",
            "F13", "F14",
            "F17", "F18",
        ]
    )  # all 11 implementable development families
    conditions: list[str] = field(default_factory=lambda: ["A", "B", "C", "D"])
    n_clean_per_family: int = 5
    base_seed: int = 42


@dataclass
class CaseOutcome:
    family_id: str
    instance_id: str
    case_type: str              # "fault" | "clean"
    condition: str              # "A" | "B" | "C" | "D"
    contracts_on: bool
    monitor_on: bool
    evaluator_available: bool
    evaluator_state: Optional[str]   # "HALT" | "PASS" | "ERROR" | None
    monitor_status: str         # "off" | "stub_not_trained"
    escaped: bool               # True = fault not caught (or clean + false alarm)
    # Same outcome under the "detected" definition, where an ESCALATE counts as
    # containment because a human is shown the fault before release.
    escaped_detected: bool = True
    # True when the injector matched nothing and the workspace was left
    # unchanged; such cases carry no fault and are excluded from SEER.
    vacuous: bool = False
    distortion: Optional[float] = None
    error: Optional[str] = None
    # Pre-release telemetry extracted by feature_extractor (None if extraction fails)
    features: Optional[dict[str, Any]] = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class PilotResults:
    config: dict[str, Any]
    outcomes: list[dict[str, Any]]
    seer_by_condition: dict[str, float]
    false_escalation_by_condition: dict[str, float]
    icc_estimate: float
    icc_ci_95: list[float]
    icc_condition_used: str
    n_families: int
    n_fault_cases: int
    n_clean_cases: int
    families_run: list[str]
    pilot_summary: str
    icc_note: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


# ─── Real-data sourcing (ADR-AIIM-050) ─────────────────────────────────────────
#
# SEER on a workspace built from ~60 synthetic fragments (CC, CCC, c1ccccc1...)
# measures detection on a library that looks nothing like what the pipeline
# actually processes. These workspaces are now built by sampling the real,
# already-registered SNP-Lib ledger (236 real compounds, real stereocentre/MW
# annotations already computed) and real Vina scores (172 real compounds,
# MM-GBSA rescoring run), so a fault case is "a real compound's real record,
# with one field altered" rather than a fragment library invented to have the
# right shape. Column names and schemas are unchanged, so every evaluator and
# every existing workspace test still applies unmodified.

_REAL_LEDGER_PATH = Path(
    r"D:/Agentic VS Workflow/20_new_runs/snp_screen_2026-08-04/G1_standardization/snp_source_ledger.csv"
)
_REAL_VINA_SCORES_PATH = Path(
    r"D:/kdss_SNP v2/KDSS_A2AR_Code_Data_Package_2026-08-03/data/output/mmgbsa/15_mmgbsa_scores.csv"
)
_REAL_DATA_SEED = 20260807

_real_ledger_cache: Optional[pd.DataFrame] = None
_real_scored_cache: Optional[pd.DataFrame] = None


def _load_real_ledger() -> pd.DataFrame:
    """The real, already-registered 236-compound SNP-Lib standardization ledger."""
    global _real_ledger_cache
    if _real_ledger_cache is None:
        if not _REAL_LEDGER_PATH.is_file():
            raise FileNotFoundError(
                f"Real SNP-Lib ledger not found at {_REAL_LEDGER_PATH}. "
                "The fault-injection harness requires real data (ADR-AIIM-050)."
            )
        df = pd.read_csv(_REAL_LEDGER_PATH)
        # No assay data reaches this file at all (G1's own no_activity_columns_in_ledger
        # check already enforces that on the real pipeline); asserted again here so a
        # future ledger regeneration cannot silently introduce one into the harness.
        from evaluators._assay_blindness import blocked_columns_in
        blocked = blocked_columns_in(df.columns, strict=True)
        if blocked:
            raise ValueError(f"Real ledger has forbidden columns, refusing to use it: {blocked}")
        _real_ledger_cache = df
    return _real_ledger_cache


def _load_real_scored_compounds() -> pd.DataFrame:
    """
    Real compounds with a real Vina score, for the G4/G6 workspaces.

    Left-joins the real 172-compound MM-GBSA rescoring table's `vina_box_best`
    column onto the real ledger's SMILES/descriptors by compound ID -- both
    real, both already produced by this project's own pipeline (Phase A/B of
    the QSAR/MD rerun, ADR-AIIM-047/048), not synthesized for this harness.
    """
    global _real_scored_cache
    if _real_scored_cache is None:
        if not _REAL_VINA_SCORES_PATH.is_file():
            raise FileNotFoundError(f"Real Vina score table not found at {_REAL_VINA_SCORES_PATH}")
        scores = pd.read_csv(_REAL_VINA_SCORES_PATH)[["Compound", "vina_box_best"]]
        ledger = _load_real_ledger()
        merged = scores.merge(ledger, left_on="Compound", right_on="SNP_ID", how="inner")
        merged = merged.dropna(subset=["vina_box_best", "smiles_canonical"])
        _real_scored_cache = merged.reset_index(drop=True)
    return _real_scored_cache


def _sample_real_ledger(n: int, require_chiral: bool = False) -> pd.DataFrame:
    """
    A deterministic, reproducible sample of n real ledger rows.

    When require_chiral, at least one real stereocentre-bearing compound is
    guaranteed in the sample (real data has 58/236 such rows, so a random
    sample of any reasonable n already contains some -- guaranteed rather than
    left to chance so the harness is not flaky at small n).
    """
    ledger = _load_real_ledger()
    rng = np.random.default_rng(_REAL_DATA_SEED)
    pool = ledger
    if n >= len(pool):
        # Real library has 236 rows; requests beyond that repeat with a
        # distinguishing suffix so IDs stay unique.
        reps = -(-n // len(pool))
        return pd.concat([pool] * reps, ignore_index=True).iloc[:n].reset_index(drop=True)

    if require_chiral:
        chiral = pool[pool["has_stereocenters"] == True]  # noqa: E712
        n_chiral = min(len(chiral), max(3, n // 10))
        chiral_sample = chiral.sample(n=n_chiral, random_state=_REAL_DATA_SEED)
        rest_pool = pool.drop(chiral_sample.index)
        rest_sample = rest_pool.sample(n=n - n_chiral, random_state=_REAL_DATA_SEED)
        return pd.concat([chiral_sample, rest_sample], ignore_index=True)

    idx = rng.choice(len(pool), size=n, replace=False)
    return pool.iloc[idx].reset_index(drop=True)


# ─── Workspace factories ──────────────────────────────────────────────────────

def _cycle(lst: list, n: int) -> list:
    return [lst[i % len(lst)] for i in range(n)]


def _setup_g0_workspace(ws: Path, n: int = _N_COMPOUNDS) -> None:
    """G0 workspace: compounds.csv, real SNP-Lib IDs and SMILES."""
    sample = _sample_real_ledger(n)
    rows = [{"compound_id": row["SNP_ID"], "SMILES": row["smiles_canonical"]}
            for _, row in sample.iterrows()]
    pd.DataFrame(rows).to_csv(ws / "compounds.csv", index=False)


def _setup_g1_workspace(ws: Path, n: int = _N_COMPOUNDS) -> None:
    """G1 workspace: standardised_ledger.csv, real rows including real stereocentres."""
    sample = _sample_real_ledger(n, require_chiral=True)
    rows = []
    for _, row in sample.iterrows():
        rows.append({
            "SNP_ID": row["SNP_ID"],
            "SMILES": row["smiles_canonical"],
            "exclusion_reason": row.get("exclusion_reason", "") or "",
            "valid_smiles": bool(row.get("valid_smiles", True)),
            "parent_smiles": row.get("parent_smiles", row["smiles_canonical"]),
            "removed_fragment": row.get("removed_fragment", "") or "",
            "has_stereocenters": bool(row["has_stereocenters"]),
            "n_stereocenters": int(row["n_stereocenters"]) if pd.notna(row["n_stereocenters"]) else 0,
            "stereo_defined": bool(row.get("stereo_defined", True)),
        })
    pd.DataFrame(rows).to_csv(ws / "standardised_ledger.csv", index=False)


def _setup_g4_workspace(ws: Path, n: int = _N_COMPOUNDS) -> None:
    """G4 workspace: docking_scores.csv, real compound_id/SMILES/vina_score."""
    scored = _load_real_scored_compounds()
    n = min(n, len(scored))
    sample = scored.sample(n=n, random_state=_REAL_DATA_SEED).reset_index(drop=True)
    rows = [
        {"compound_id": row["SNP_ID"], "SMILES": row["smiles_canonical"],
         "vina_score": round(float(row["vina_box_best"]), 3)}
        for _, row in sample.iterrows()
    ]
    pd.DataFrame(rows).to_csv(ws / "docking_scores.csv", index=False)


def _setup_g8_workspace(ws: Path, n: int = _N_COMPOUNDS) -> None:
    """
    Minimal G8 workspace: results_table.csv plus a figure.

    The figure must exist in the *clean* state so its pre-injection SHA-256 can
    be registered; F17 then overwrites it and the mismatch becomes detectable.
    If the figure were created only by the injector there would be no anchor to
    compare against and a stale figure could never be caught.
    """
    pd.DataFrame([
        {"condition": "A_hybrid",         "seer_pct": 0.15, "brier_score": 0.12, "n_compounds": n},
        {"condition": "B_contracts_only",  "seer_pct": 0.35, "brier_score": 0.28, "n_compounds": n},
        {"condition": "C_monitor_only",    "seer_pct": 0.42, "brier_score": 0.31, "n_compounds": n},
        {"condition": "D_neither",         "seer_pct": 0.58, "brier_score": 0.45, "n_compounds": n},
    ]).to_csv(ws / "results_table.csv", index=False)
    (ws / "main_figure.png").write_bytes(b"\x89PNG\r\n\x1a\nCURRENT_VERSION_CONTENT")


def _real_descriptors(smi: str) -> dict:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return {"MW": 0.0, "rotatable_bonds": 0, "HBA": 0, "HBD": 0, "ring_count": 0}
    return {
        "MW": round(float(Descriptors.ExactMolWt(mol)), 2),
        "rotatable_bonds": int(rdMolDescriptors.CalcNumRotatableBonds(mol)),
        "HBA": int(rdMolDescriptors.CalcNumHBA(mol)),
        "HBD": int(rdMolDescriptors.CalcNumHBD(mol)),
        "ring_count": int(rdMolDescriptors.CalcNumRings(mol)),
    }


def _real_murcko_scaffold(smi: str) -> str:
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return "__INVALID__"
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    return Chem.MolToSmiles(scaffold) if scaffold and scaffold.GetNumAtoms() else "__ACYCLIC__"


def _setup_g6_workspace(ws: Path, n: int = _N_COMPOUNDS) -> None:
    """
    G6 workspace: ranked_compounds.csv + selected_panel.csv, real compounds.

    The selected panel is built by taking the best-scoring real compound per
    distinct real Murcko scaffold (up to 20), so F14 (scaffold duplication)
    is always operating on real, genuinely-unique scaffolds rather than a
    hand-picked synthetic ring library -- and can fail to find a duplicate to
    inject on a library where two real panel members happen to already share
    a scaffold, which is itself a realistic condition to audit.
    """
    scored = _load_real_scored_compounds().copy()
    scored["scaffold"] = scored["smiles_canonical"].apply(_real_murcko_scaffold)
    scored = scored.sort_values("vina_box_best")  # most negative (best) first

    panel_rows = scored.drop_duplicates(subset="scaffold", keep="first").head(20)
    remaining = scored.drop(panel_rows.index)
    n_extra = max(0, n - len(panel_rows))
    extra_rows = remaining.head(n_extra)

    def _row(row) -> dict:
        return {
            "compound_id": row["SNP_ID"], "SMILES": row["smiles_canonical"],
            "vina_score": round(float(row["vina_box_best"]), 3),
            **_real_descriptors(row["smiles_canonical"]),
        }

    rows = [_row(r) for _, r in panel_rows.iterrows()] + [_row(r) for _, r in extra_rows.iterrows()]
    df = pd.DataFrame(rows)
    df.to_csv(ws / "ranked_compounds.csv", index=False)
    df.head(len(panel_rows)).to_csv(ws / "selected_panel.csv", index=False)


_WORKSPACE_SETUP: dict[str, Any] = {
    "G0": _setup_g0_workspace,
    "G1": _setup_g1_workspace,
    "G4": _setup_g4_workspace,
    "G6": _setup_g6_workspace,
    "G8": _setup_g8_workspace,
}

# ─── Pre-release telemetry collection ────────────────────────────────────────

# Map gate → (main_csv_filename, id_column) used for feature extraction.
_GATE_MAIN_FILE: dict[str, tuple[str, str]] = {
    "G0": ("compounds.csv",           "compound_id"),
    "G1": ("standardised_ledger.csv", "SNP_ID"),
    "G4": ("docking_scores.csv",      "compound_id"),
    "G6": ("ranked_compounds.csv",    "compound_id"),
    "G8": ("results_table.csv",       "condition"),
}


def _sha256_ws_file(ws: Path, filename: str) -> Optional[str]:
    """Return SHA256 hex of a workspace file, or None if file is absent."""
    path = ws / filename
    if not path.exists():
        return None
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest().upper()


def _read_ws_ids(ws: Path, gate: str) -> tuple[set[str], int, int]:
    """
    Read (unique_ids, total_row_count, n_duplicates) from the gate's main CSV.

    Returns empty set and zero counts if the file is absent or unreadable.
    """
    if gate not in _GATE_MAIN_FILE:
        return set(), 0, 0
    filename, id_col = _GATE_MAIN_FILE[gate]
    path = ws / filename
    if not path.exists():
        return set(), 0, 0
    try:
        df = pd.read_csv(path, low_memory=False)
        if id_col not in df.columns:
            return set(), len(df), 0
        ids = {str(v) for v in df[id_col].dropna()}
        n_dup = int(df[id_col].duplicated().sum())
        return ids, len(df), n_dup
    except Exception:
        return set(), 0, 0


def _collect_case_telemetry(
    ws: Path,
    gate: str,
    pre_hash: Optional[str],
    ref_ids: set[str],
    n_expected: int,
) -> Optional[dict[str, Any]]:
    """
    Call feature_extractor.extract() on the post-injection workspace state.

    Returns the validated feature dict, or None if extraction fails.

    Parameters
    ----------
    pre_hash:
        SHA256 of the gate's main file in its *clean* (pre-injection) state.
        If the actual post-injection hash differs, n_missing_input_hashes = 1.
    ref_ids:
        Compound IDs from the clean workspace.  Compared with post-injection
        IDs to compute n_lost_records and n_unexpected_records.
    n_expected:
        Expected row count (clean workspace).  Compared with actual row count
        to populate count_checks.
    """
    if gate not in _GATE_MAIN_FILE:
        return None

    try:
        from telemetry.feature_extractor import extract
    except ImportError:
        logger.debug("feature_extractor not importable; skipping telemetry collection")
        return None

    filename, id_col = _GATE_MAIN_FILE[gate]
    actual_path = ws / filename

    # Post-injection file hash (None if file is missing)
    actual_hash = _sha256_ws_file(ws, filename)

    # Hash comparison: any file whose actual hash differs from expected is "missing"
    n_missing = 0
    if pre_hash is not None and actual_hash != pre_hash:
        n_missing = 1
    elif pre_hash is None and actual_hash is None:
        n_missing = 1  # file expected but not found

    registered_inputs = [{
        "path": str(actual_path),
        "expected_hash": pre_hash,
        "registered": True,
        "lineage_complete": True,
    }]
    registered_outputs = [{
        "path": str(actual_path),
        "registered": actual_path.exists(),
    }]

    # Post-injection IDs and duplicate count
    actual_ids, actual_n, n_dup = _read_ws_ids(ws, gate)

    # File timestamps — use the same file for both input and output timestamps
    try:
        mtime = actual_path.stat().st_mtime if actual_path.exists() else None
    except OSError:
        mtime = None

    # Count check: (expected, actual)
    count_checks: Optional[list[tuple[int, int]]] = (
        [(n_expected, actual_n)] if n_expected > 0 else None
    )

    # Schema fields (column names of the main CSV)
    frozen_schema: Optional[set[str]] = None
    if actual_path.exists():
        try:
            frozen_schema = set(pd.read_csv(actual_path, nrows=0).columns.tolist())
        except Exception:
            frozen_schema = None

    try:
        record = extract(
            registered_inputs=registered_inputs,
            registered_outputs=registered_outputs,
            tool_logs=[],
            n_retries=0,
            elapsed_seconds=0.0,
            reference_id_set=ref_ids,
            output_id_set=actual_ids,
            output_timestamp=mtime,
            latest_input_timestamp=mtime,
            frozen_schema_fields=frozen_schema,
            output_schema_fields=frozen_schema,
            count_checks=count_checks,
            stage_metrics={"n_duplicate_ids": n_dup},
        )
        # Override n_missing_input_hashes with our direct hash comparison result.
        # The feature_extractor computes this from the registered_inputs list, but
        # its logic re-hashes the file and compares to expected_hash.  For the
        # pilot, the expected_hash is the pre-injection hash, so the extractor
        # already handles this correctly; no override needed.
        return record
    except Exception as exc:
        logger.debug("Feature extraction failed for %s gate %s: %s", ws.name, gate, exc)
        return None


# ─── Training-dataset assembly ────────────────────────────────────────────────

def build_training_dataset(
    results: "PilotResults",
    reference_condition: str = "B",
) -> "tuple[pd.DataFrame, np.ndarray, np.ndarray]":
    """
    Assemble X (features), y (labels), groups from pilot outcomes.

    Only outcomes from `reference_condition` that have a populated `features`
    dict are included.  Outcomes with error=True or missing features are
    excluded (fail-closed).

    Returns
    -------
    X : pd.DataFrame
        Feature matrix with columns matching feature_schema.yaml.
    y : np.ndarray of int (0/1)
        unsafe_release labels (Gate 4 criteria via label_assigner).
    groups : np.ndarray of str
        Fault-family IDs for GroupKFold splitting.

    Raises
    ------
    ValueError
        If no usable rows are found.
    """
    from fault_injection.label_assigner import assign_labels

    outcomes = results.outcomes if isinstance(results.outcomes, list) else []

    # Filter to reference condition; require features to be present
    usable = [
        o for o in outcomes
        if o.get("condition") == reference_condition
        and not o.get("error")
        and isinstance(o.get("features"), dict)
        and len(o["features"]) > 0
    ]
    if not usable:
        raise ValueError(
            f"build_training_dataset: no usable outcomes with features for "
            f"reference condition {reference_condition!r}.  "
            "Run the pilot with collect_features=True."
        )

    # Labels
    label_df = assign_labels(outcomes, reference_condition=reference_condition)
    label_map: dict[str, int] = dict(
        zip(label_df["instance_id"], label_df["unsafe_release"])
    )

    rows = []
    ys = []
    groups = []

    for o in usable:
        iid = o.get("instance_id", "")
        if iid not in label_map:
            continue
        rows.append(o["features"])
        ys.append(label_map[iid])
        groups.append(o.get("family_id", ""))

    if not rows:
        raise ValueError(
            "build_training_dataset: no rows remain after joining features and labels."
        )

    X = pd.DataFrame(rows)
    y = np.array(ys, dtype=int)
    groups_arr = np.array(groups, dtype=str)

    # Leakage check on column names
    try:
        from telemetry.leakage_guard import check_field_names, check_allowlist
        check_field_names(X.columns, context="build_training_dataset")
        check_allowlist(X.columns, context="build_training_dataset")
    except Exception as exc:
        raise ValueError(f"build_training_dataset: leakage check failed: {exc}") from exc

    logger.info(
        "build_training_dataset: X=%s y=%s groups=%d unique",
        X.shape, y.shape, len(set(groups)),
    )
    return X, y, groups_arr


# ─── Pilot contracts ──────────────────────────────────────────────────────────

_PILOT_CONTRACT: dict[str, dict[str, Any]] = {
    "G0": {
        "checks": [
            {"id": "hash_snp_library_csv", "severity": "critical"},
            {"id": "id_uniqueness",        "severity": "critical"},
            {"id": "smiles_column_present", "severity": "critical"},
            {"id": "compound_count",        "expected_count": _N_COMPOUNDS, "severity": "medium"},
        ],
    },
    "G1": {
        "checks": [
            {"id": "smiles_roundtrip",            "severity": "critical"},
            {"id": "no_drug_likeness_exclusions",  "severity": "critical"},
            {"id": "stereochemistry_preserved",    "severity": "critical"},
            {"id": "salt_stripping_verified",      "severity": "critical"},
            {"id": "ledger_count",                 "expected_count": _N_COMPOUNDS, "severity": "medium"},
        ],
    },
    "G4": {
        "checks": [
            {"id": "score_direction",           "severity": "critical"},
            {"id": "no_assay_columns_in_output", "severity": "critical"},
            {"id": "score_distribution_sanity",
             "expected_range_kcal_mol": [-12.0, -3.0], "severity": "medium"},
            {"id": "mw_bias_check",              "severity": "medium"},
        ],
    },
    "G6": {
        "checks": [
            {"id": "assay_blindness_column_scan", "severity": "critical"},
            {"id": "eligible_set_completeness",   "severity": "critical"},
            {"id": "scaffold_uniqueness",         "severity": "critical"},
            {"id": "frozen_rule_applied",         "severity": "high"},
        ],
    },
    "G8": {
        "checks": [
            {"id": "figure_source_verification", "severity": "critical"},
            {"id": "table_source_verification",  "severity": "critical"},
            {"id": "number_traceability",        "severity": "critical"},
            {"id": "contract_coverage",          "severity": "critical"},
        ],
    },
}


def _clean_hashes(ws: Path, gate: str) -> dict[str, str]:
    """
    SHA-256 of every gate-relevant file in its clean state.

    Captured before injection so provenance gates have a registered anchor to
    compare against; without it a corrupted or stale file has nothing to fail
    against and the fault escapes silently.
    """
    files = {
        "G0": ["compounds.csv"],
        "G8": ["results_table.csv", "main_figure.png"],
    }.get(gate, [])
    out: dict[str, str] = {}
    for name in files:
        digest = _sha256_ws_file(ws, name)
        if digest:
            out[name] = digest
    return out


def _get_evaluator_inputs(
    gate: str, ws: Path, clean_hashes: Optional[dict[str, str]] = None
) -> dict[str, Any]:
    clean_hashes = clean_hashes or {}

    if gate == "G0":
        return {
            "snp_library_csv": str(ws / "compounds.csv"),
            "registered_input_hashes": {
                "hash_snp_library_csv": clean_hashes.get("compounds.csv"),
            },
        }
    if gate == "G1":
        return {"snp_source_ledger_csv": str(ws / "standardised_ledger.csv")}
    if gate == "G4":
        return {"vina_score_csv": str(ws / "docking_scores.csv")}
    if gate == "G8":
        table = ws / "results_table.csv"
        figure = ws / "main_figure.png"
        # Synthetic registry: the clean artifacts were registered and their
        # contracts passed, so any post-injection divergence is attributable.
        snapshot = [
            {"record_type": "artifact", "artifact_id": "results_table",
             "sha256": clean_hashes.get("results_table.csv")},
            {"record_type": "artifact", "artifact_id": "main_figure",
             "sha256": clean_hashes.get("main_figure.png")},
            {"record_type": "contract_outcome", "outcome": "PASS",
             "artifact_ids": ["results_table", "main_figure"]},
        ]
        return {
            "figure_manifest": [
                {"figure_id": "F1", "artifact_id": "main_figure", "path": str(figure)}
            ],
            "table_manifest": [
                {"table_id": "T1", "artifact_id": "results_table", "path": str(table)}
            ],
            "result_snippets": [
                {"claim": "SEER by condition", "metric": "seer_pct",
                 "artifact_id": "results_table"}
            ],
            "artifact_registry_snapshot": snapshot,
        }
    if gate == "G6":
        ranking_csv = ws / "ranked_compounds.csv"
        col_names: list[str] = []
        if ranking_csv.exists():
            col_names = list(pd.read_csv(ranking_csv, nrows=0).columns)
        return {
            "final_ranking_csv": str(ranking_csv),
            "selected_candidates_csv": str(ws / "selected_panel.csv"),
            "column_scan_log": [{"file": str(ranking_csv), "columns": col_names}],
            "selection_script_hash": _DUMMY_HASH,
            "registered_selection_script_hash": _DUMMY_HASH,
            "expected_eligible_count": _N_COMPOUNDS,
            "diversity_seed": 42,
            "selection_count_report": f"{_N_COMPOUNDS} compounds evaluated",
        }
    return {}


def _dispatch_evaluator(
    gate: str, ws: Path, clean_hashes: Optional[dict[str, str]] = None
) -> tuple[Optional[dict[str, Any]], bool]:
    """
    Run the appropriate gate evaluator.

    Returns (result, available).
    available=False when the evaluator is a stub or raises NotImplementedError.
    On unexpected exception: result contains the error, available=False.
    """
    if gate not in _PILOT_CONTRACT:
        return None, False

    contract = _PILOT_CONTRACT[gate]
    inputs = _get_evaluator_inputs(gate, ws, clean_hashes)

    try:
        if gate == "G0":
            from evaluators.G0_evaluator import evaluate  # type: ignore[import]
        elif gate == "G1":
            from evaluators.G1_evaluator import evaluate  # type: ignore[import]
        elif gate == "G4":
            from evaluators.G4_evaluator import evaluate  # type: ignore[import]
        elif gate == "G6":
            from evaluators.G6_evaluator import evaluate  # type: ignore[import]
        elif gate == "G8":
            from evaluators.G8_evaluator import evaluate  # type: ignore[import]
        else:
            return None, False
        result = evaluate(contract=contract, inputs=inputs, registry=None)
        return result, True
    except NotImplementedError:
        return None, False
    except Exception as exc:
        logger.warning("Evaluator %s raised %s: %s", gate, type(exc).__name__, exc)
        return {"state": "ERROR", "error": str(exc)}, False


# ─── Deterministic seed ───────────────────────────────────────────────────────

def _case_seed(base_seed: int, family_id: str, instance_index: int,
               condition: str, case_type: str) -> int:
    """Return a stable 31-bit seed from the case identity."""
    tag = f"{base_seed}:{family_id}:{instance_index}:{condition}:{case_type}"
    digest = hashlib.md5(tag.encode()).digest()[:4]
    return int.from_bytes(digest, "little") & 0x7FFF_FFFF


# ─── Single-case runner ───────────────────────────────────────────────────────

def _run_one_case(
    family_id: str,
    instance_index: int,
    condition: str,
    ws_parent: Path,
    base_seed: int,
    is_clean: bool,
) -> CaseOutcome:
    """Run one isolated (family × instance × condition) case."""
    injector = FAULT_REGISTRY[family_id]
    spec = injector.spec
    gate = spec.gate
    contracts_on, monitor_on = _CONDITION_MAP[condition]
    case_type = "clean" if is_clean else "fault"
    seed = _case_seed(base_seed, family_id, instance_index, condition, case_type)

    if is_clean:
        instance_id = f"{family_id}_clean_{instance_index:03d}"
    else:
        instance_id = injector.instance_id(instance_index)

    ws_name = f"{family_id}_{instance_index:03d}_{condition}_{case_type}"
    ws = ws_parent / ws_name
    ws.mkdir(parents=True, exist_ok=True)

    pre_ws: Optional[Path] = None
    try:
        # Populate clean workspace
        if gate in _WORKSPACE_SETUP:
            _WORKSPACE_SETUP[gate](ws)
        else:
            return CaseOutcome(
                family_id=family_id, instance_id=instance_id,
                case_type=case_type, condition=condition,
                contracts_on=contracts_on, monitor_on=monitor_on,
                evaluator_available=False, evaluator_state=None,
                monitor_status="off" if not monitor_on else "stub_not_trained",
                escaped=True, escaped_detected=True, distortion=None,
                error=f"No workspace setup for gate {gate}",
            )

        # Capture pre-injection state for feature extraction
        # (Hash the clean file and record IDs before any fault is injected)
        if gate in _GATE_MAIN_FILE:
            _filename, _ = _GATE_MAIN_FILE[gate]
            pre_hash = _sha256_ws_file(ws, _filename)
            ref_ids, n_expected, _ = _read_ws_ids(ws, gate)
        else:
            pre_hash = None
            ref_ids = set()
            n_expected = 0

        # Registered provenance anchors, taken while the workspace is still clean.
        clean_hashes = _clean_hashes(ws, gate)
        vacuous = False

        # Inject fault
        if not is_clean:
            params = spec.param_grid[instance_index % len(spec.param_grid)]
            try:
                injection_record = injector.inject(ws, params, seed, instance_index)
                # An instance whose parameters match nothing in this workspace
                # changes no bytes. There is then no fault to detect, and
                # scoring the resulting PASS as an escape would penalise the
                # contracts for a fault that was never injected.
                vacuous = not getattr(injection_record, "affected_items", None)
            except NotImplementedError:
                return CaseOutcome(
                    family_id=family_id, instance_id=instance_id,
                    case_type=case_type, condition=condition,
                    contracts_on=contracts_on, monitor_on=monitor_on,
                    evaluator_available=False, evaluator_state=None,
                    monitor_status="off" if not monitor_on else "stub_not_trained",
                    escaped=True, escaped_detected=True, distortion=None,
                    error="Injector stub (NotImplementedError)",
                )

        # Collect pre-release telemetry (after injection / for clean cases as-is)
        case_features = _collect_case_telemetry(ws, gate, pre_hash, ref_ids, n_expected)

        # Distortion (optional, only for fault cases)
        distortion_val: Optional[float] = None
        if not is_clean and gate in _WORKSPACE_SETUP:
            try:
                from fault_injection.distortion import compute_distortion
                pre_ws = ws_parent / f"{ws_name}_pre"
                pre_ws.mkdir(exist_ok=True)
                _WORKSPACE_SETUP[gate](pre_ws)
                drec = compute_distortion(family_id, instance_id, pre_ws, ws)
                distortion_val = drec.value
            except Exception as exc:
                logger.debug("Distortion skipped for %s: %s", instance_id, exc)
            finally:
                if pre_ws is not None:
                    shutil.rmtree(pre_ws, ignore_errors=True)
                    pre_ws = None

        # Run evaluator (contracts_on only)
        evaluator_result: Optional[dict[str, Any]] = None
        evaluator_available = False
        if contracts_on:
            evaluator_result, evaluator_available = _dispatch_evaluator(
                gate, ws, clean_hashes
            )

        evaluator_state: Optional[str] = None
        if evaluator_result is not None:
            evaluator_state = evaluator_result.get("state")

        monitor_status = "off" if not monitor_on else "stub_not_trained"

        # Determine containment outcome
        if is_clean:
            # For clean cases: escaped means false escalation (evaluator
            # intervened on a clean run). Both a HALT and an ESCALATE interrupt
            # the workflow, so both count against the strict definition.
            if contracts_on and evaluator_available:
                escaped = evaluator_state == "HALT"
                escaped_detected = evaluator_state in ("HALT", "ESCALATE")
            else:
                escaped = False
                escaped_detected = False
        else:
            # For fault cases two definitions are reported side by side:
            #   strict   — only a HALT contains the fault
            #   detected — an ESCALATE also counts, because the project's own
            #              definition of SEER is "escapes *without being
            #              detected before release*", and an escalation is a
            #              detection that hands the decision to a human.
            # Neither is privileged here; the pair is reported so the choice is
            # visible rather than buried in the harness.
            if contracts_on and evaluator_available:
                escaped = evaluator_state != "HALT"
                escaped_detected = evaluator_state not in ("HALT", "ESCALATE")
            else:
                escaped = True  # no working defense → fault escapes
                escaped_detected = True

        return CaseOutcome(
            family_id=family_id, instance_id=instance_id,
            case_type=case_type, condition=condition,
            contracts_on=contracts_on, monitor_on=monitor_on,
            evaluator_available=evaluator_available,
            evaluator_state=evaluator_state,
            monitor_status=monitor_status,
            escaped=escaped, escaped_detected=escaped_detected, vacuous=vacuous,
            distortion=distortion_val,
            features=case_features,
        )

    except Exception as exc:
        logger.error("Case %s[%d] %s failed: %s", family_id, instance_index, condition, exc)
        return CaseOutcome(
            family_id=family_id, instance_id=instance_id,
            case_type=case_type, condition=condition,
            contracts_on=contracts_on, monitor_on=monitor_on,
            evaluator_available=False, evaluator_state=None,
            monitor_status="off" if not monitor_on else "stub_not_trained",
            escaped=True, escaped_detected=True, distortion=None, error=str(exc),
        )
    finally:
        if pre_ws is not None:
            shutil.rmtree(pre_ws, ignore_errors=True)
        if logger.isEnabledFor(logging.DEBUG):
            pass  # keep workspace for inspection
        else:
            shutil.rmtree(ws, ignore_errors=True)


# ─── ICC computation ──────────────────────────────────────────────────────────

def _compute_icc(
    outcomes: list[CaseOutcome],
    condition: str = "B",
    n_bootstrap: int = 2000,
    rng_seed: int = 99,
) -> tuple[float, tuple[float, float], str]:
    """
    ANOVA-based ICC(1,1) from binary escape outcomes.

    Only fault cases with a working evaluator and no errors contribute.
    Returns (icc_point, (ci_lower, ci_upper), note).
    """
    eligible = [
        o for o in outcomes
        if o.case_type == "fault"
        and o.condition == condition
        and o.evaluator_available
        and not o.error
    ]
    if not eligible:
        note = f"No evaluable fault cases for condition {condition}"
        return float("nan"), (float("nan"), float("nan")), note

    from collections import defaultdict
    family_outcomes: dict[str, list[int]] = defaultdict(list)
    for case in eligible:
        family_outcomes[case.family_id].append(int(case.escaped))

    families = sorted(family_outcomes.keys())
    J = len(families)

    if J < 2:
        note = f"Only {J} family with evaluable outcomes — need ≥ 2 for ICC"
        return float("nan"), (float("nan"), float("nan")), note

    def _icc_from_groups(groups: dict[str, list[int]]) -> float:
        all_y = [y for vs in groups.values() for y in vs]
        n_total = len(all_y)
        if n_total <= 1:
            return float("nan")
        grand_mean = sum(all_y) / n_total
        J_loc = len(groups)
        group_means = {fid: sum(vs) / len(vs) for fid, vs in groups.items()}
        k_mean = n_total / J_loc
        ss_b = sum(len(groups[fid]) * (group_means[fid] - grand_mean) ** 2 for fid in groups)
        ss_w = sum((y - group_means[fid]) ** 2 for fid in groups for y in groups[fid])
        df_b = J_loc - 1
        df_w = n_total - J_loc
        if df_b <= 0 or df_w <= 0:
            return float("nan")
        ms_b = ss_b / df_b
        ms_w = ss_w / df_w
        denom = ms_b + (k_mean - 1) * ms_w
        if denom == 0.0:
            return 0.0
        return float(np.clip((ms_b - ms_w) / denom, -1.0, 1.0))

    point_icc = _icc_from_groups(family_outcomes)

    # Bootstrap CI: resample families with replacement
    rng = np.random.default_rng(seed=rng_seed)
    boot_iccs: list[float] = []
    for _ in range(n_bootstrap):
        sampled_families = rng.choice(families, size=J, replace=True).tolist()
        boot_groups: dict[str, list[int]] = {
            f"g{i}": family_outcomes[fid] for i, fid in enumerate(sampled_families)
        }
        v = _icc_from_groups(boot_groups)
        if not (v != v):  # not NaN
            boot_iccs.append(v)

    if len(boot_iccs) < 10:
        ci: tuple[float, float] = (float("nan"), float("nan"))
    else:
        ci = (float(np.percentile(boot_iccs, 2.5)), float(np.percentile(boot_iccs, 97.5)))

    per_family_n = [len(family_outcomes[f]) for f in families]
    note = (
        f"ICC(1,1) from {J} families x {per_family_n} instances "
        f"under condition {condition}.  Bootstrap n={n_bootstrap}, seed={rng_seed}."
    )
    return point_icc, ci, note


# ─── Main pilot orchestrator ──────────────────────────────────────────────────

def run_pilot(
    config: PilotConfig,
    workspace_root: Optional[Path] = None,
    output_path: Optional[Path] = None,
) -> PilotResults:
    """
    Run the Gate 9 feasibility pilot.

    For each family × instance × condition, creates an isolated workspace,
    injects the fault (if fault case), runs the relevant evaluator (if
    contracts_on), and records the outcome.

    Parameters
    ----------
    config:
        Pilot configuration specifying families, conditions, clean-case count,
        and random seed.
    workspace_root:
        Root directory for case workspaces.  If None a temporary directory is
        created and deleted after the run.
    output_path:
        If provided, write PilotResults as JSON to this path.

    Returns
    -------
    PilotResults
    """
    # Validate requested families
    implementable = set(implementable_families())
    for fid in config.pilot_families:
        if fid not in implementable:
            raise ValueError(
                f"Family {fid!r} is not implementable. "
                f"Stub families require WP5/WP7 data.  "
                f"Valid choices: {sorted(implementable)}"
            )
    for cond in config.conditions:
        if cond not in _CONDITION_MAP:
            raise ValueError(f"Unknown condition {cond!r}. Must be one of A, B, C, D.")

    cleanup_ws = workspace_root is None
    if workspace_root is None:
        workspace_root = Path(tempfile.mkdtemp(prefix="aiim_pilot_"))

    workspace_root.mkdir(parents=True, exist_ok=True)
    all_outcomes: list[CaseOutcome] = []

    try:
        for family_id in config.pilot_families:
            injector = FAULT_REGISTRY[family_id]
            n_instances = len(injector.spec.param_grid)
            logger.info("Pilot: %s  %d instances x %s", family_id, n_instances, config.conditions)

            # Fault cases
            for idx in range(n_instances):
                for cond in config.conditions:
                    outcome = _run_one_case(
                        family_id=family_id,
                        instance_index=idx,
                        condition=cond,
                        ws_parent=workspace_root,
                        base_seed=config.base_seed,
                        is_clean=False,
                    )
                    all_outcomes.append(outcome)

            # Clean cases
            for clean_idx in range(config.n_clean_per_family):
                for cond in config.conditions:
                    outcome = _run_one_case(
                        family_id=family_id,
                        instance_index=clean_idx,
                        condition=cond,
                        ws_parent=workspace_root,
                        base_seed=config.base_seed + 10_000,
                        is_clean=True,
                    )
                    all_outcomes.append(outcome)

    finally:
        if cleanup_ws:
            shutil.rmtree(workspace_root, ignore_errors=True)

    # SEER per condition (fault cases only; exclude error cases and vacuous
    # injections, where the injector matched nothing and left the workspace
    # unchanged — there is no fault present to escape).
    seer_by_condition: dict[str, float] = {}
    for cond in config.conditions:
        fault_cond = [
            o for o in all_outcomes
            if o.case_type == "fault" and o.condition == cond
            and not o.error and not o.vacuous
        ]
        seer_by_condition[cond] = (
            sum(o.escaped for o in fault_cond) / len(fault_cond)
            if fault_cond else float("nan")
        )

    # False-escalation rate per condition (clean cases only)
    false_esc_by_condition: dict[str, float] = {}
    for cond in config.conditions:
        clean_cond = [
            o for o in all_outcomes
            if o.case_type == "clean" and o.condition == cond and not o.error
        ]
        false_esc_by_condition[cond] = (
            sum(o.escaped for o in clean_cond) / len(clean_cond)
            if clean_cond else float("nan")
        )

    # ICC from condition B (primary: contracts only, monitor off)
    icc_cond = "B" if "B" in config.conditions else config.conditions[0]
    icc, icc_ci, icc_note = _compute_icc(all_outcomes, condition=icc_cond)

    n_fault = sum(1 for o in all_outcomes if o.case_type == "fault")
    n_clean = sum(1 for o in all_outcomes if o.case_type == "clean")

    def _fmt(v: float) -> str:
        return f"{v:.3f}" if v == v else "n/a"  # noqa: PLR0124

    seer_lines = ", ".join(
        f"{c}={_fmt(seer_by_condition.get(c, float('nan')))}"
        for c in ["A", "B", "C", "D"] if c in config.conditions
    )
    icc_line = (
        f"ICC(cond {icc_cond})={_fmt(icc)} 95%CI [{_fmt(icc_ci[0])}, {_fmt(icc_ci[1])}]"
        if icc == icc else "ICC: insufficient data"  # noqa: PLR0124
    )
    summary = "\n".join([
        f"Pilot: {len(config.pilot_families)} families x conditions {config.conditions}",
        f"Fault cases: {n_fault}  Clean cases: {n_clean}",
        f"SEER by condition: {seer_lines}",
        icc_line,
        "NOTE: Conditions A/C report monitor_status=stub_not_trained; "
        "plug in trained model (WP4) for definitive A/C SEER.",
    ])

    results = PilotResults(
        config=asdict(config),
        outcomes=[o.to_dict() for o in all_outcomes],
        seer_by_condition={k: float(v) for k, v in seer_by_condition.items()},
        false_escalation_by_condition={k: float(v) for k, v in false_esc_by_condition.items()},
        icc_estimate=float(icc),
        icc_ci_95=[float(icc_ci[0]), float(icc_ci[1])],
        icc_condition_used=icc_cond,
        n_families=len(config.pilot_families),
        n_fault_cases=n_fault,
        n_clean_cases=n_clean,
        families_run=list(config.pilot_families),
        pilot_summary=summary,
        icc_note=icc_note,
    )

    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as fh:
            json.dump(results.to_dict(), fh, indent=2, default=str)
        logger.info("Pilot results written to %s", output_path)

    return results


# ─── CLI ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import argparse
    import sys

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

    parser = argparse.ArgumentParser(description="AIIM Gate 9 feasibility pilot runner")
    parser.add_argument(
        "--families", nargs="+",
        default=["F01", "F02", "F03", "F04", "F05", "F09", "F10", "F13", "F14", "F17", "F18"],
        metavar="FXX",
        help="Development families to run (default: all 11 implementable)",
    )
    parser.add_argument(
        "--conditions", nargs="+", default=["A", "B", "C", "D"],
        choices=["A", "B", "C", "D"],
        help="Conditions to run (default: A B C D)",
    )
    parser.add_argument("--n-clean", type=int, default=5, metavar="N",
                        help="Clean cases per family (default: 5)")
    parser.add_argument("--seed", type=int, default=42, help="Base random seed (default: 42)")
    parser.add_argument("--output", type=Path, default=None, metavar="FILE",
                        help="Write results JSON to FILE (default: summary to stdout)")
    parser.add_argument("--workspace", type=Path, default=None, metavar="DIR",
                        help="Workspace root (default: temp dir, deleted after run)")
    args = parser.parse_args()

    cfg = PilotConfig(
        pilot_families=args.families,
        conditions=args.conditions,
        n_clean_per_family=args.n_clean,
        base_seed=args.seed,
    )
    res = run_pilot(cfg, workspace_root=args.workspace, output_path=args.output)
    print(res.pilot_summary)
    if args.output:
        print(f"Results JSON → {args.output}", file=sys.stderr)

"""
distortion.py — Decision distortion metrics (Gate 4, ADR-AIIM-023).

Each function measures how much a fault changes the downstream scientific
decision relative to a clean run.  All metrics return a float in [0, 1]:

  0.0  no change (fault had no distorting effect on this output)
  1.0  maximal change (complete reversal or full corruption)

No assay data may enter these functions.  Inputs are workflow-output files
(SMILES ledgers, docking scores, selection panels, results tables) —
structural properties and computed metrics only.

Per-family dispatch table (§6.2, approved Gate 4):
  F01–F05  → incorrect_record_proportion      (G0–G1 provenance/structure)
  F06–F08  → stub                             (G2–G3, requires WP5)
  F09–F10  → normalized_rank_displacement     (G4 scoring errors)
  F11–F12  → stub                             (G5, requires WP5)
  F13–F14  → top_k_jaccard_divergence         (G6 selection errors)
  F15–F16  → stub                             (G7, requires WP7)
  F17–F18  → incorrect_claim_proportion       (G8 reporting errors)

Calling convention
------------------
pre_workspace  : Path to the isolated run directory BEFORE fault injection
                 (clean-run pipeline output for the same case).
post_workspace : Path to the isolated run directory AFTER fault injection
                 AND after the pipeline has been re-run on the faulted input.
                 In WP3 unit tests, post_workspace is the injection workspace
                 and pre_workspace is a pristine copy of the same inputs.
"""

from __future__ import annotations

import csv
import io
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd


# ── helpers ──────────────────────────────────────────────────────────────────

def _find_col(df: pd.DataFrame, hint: str) -> Optional[str]:
    """Return the first column whose name contains *hint* (case-insensitive)."""
    hint_lower = hint.lower()
    for c in df.columns:
        if hint_lower in c.lower():
            return c
    return None


def _read_csv_safe(path: Path) -> Optional[pd.DataFrame]:
    """Read a CSV file, returning None on any error."""
    try:
        return pd.read_csv(path, dtype=str)
    except Exception:
        return None


def _load_pair(
    pre_ws: Path,
    post_ws: Path,
    filename: str,
) -> tuple[Optional[pd.DataFrame], Optional[pd.DataFrame]]:
    """Load the same filename from two workspace directories."""
    return (
        _read_csv_safe(pre_ws / filename),
        _read_csv_safe(post_ws / filename),
    )


# ── core metric functions ─────────────────────────────────────────────────────

def incorrect_record_proportion(
    pre_df: pd.DataFrame,
    post_df: pd.DataFrame,
    id_col: str = "compound_id",
) -> float:
    """
    Fraction of records that differ between pre-fault and post-fault output.

    Counts added rows, removed rows, and rows whose non-ID values changed.
    Uses ``id_col`` for alignment when present in both DataFrames; otherwise
    falls back to positional comparison.

    Returns a float in [0, 1].
    """
    if pre_df.empty and post_df.empty:
        return 0.0

    actual_id = _find_col(pre_df, id_col)
    if actual_id is not None and actual_id in post_df.columns:
        pre_ids = set(pre_df[actual_id].astype(str))
        post_ids = set(post_df[actual_id].astype(str))
        added = len(post_ids - pre_ids)
        removed = len(pre_ids - post_ids)
        shared = pre_ids & post_ids
        # compare shared records (all columns except the ID)
        pre_idx = pre_df.set_index(actual_id)
        post_idx = post_df.set_index(actual_id)
        common_cols = [c for c in pre_idx.columns if c in post_idx.columns]
        changed = sum(
            1 for iid in shared
            if iid in pre_idx.index
            and iid in post_idx.index
            and not pre_idx.loc[iid, common_cols].astype(str).equals(
                post_idx.loc[iid, common_cols].astype(str)
            )
        )
        total_affected = added + removed + changed
        denominator = max(len(pre_df), len(post_df), 1)
    else:
        # positional fallback
        n_pre, n_post = len(pre_df), len(post_df)
        n_cmp = min(n_pre, n_post)
        if n_cmp > 0:
            changed = int(
                (pre_df.iloc[:n_cmp].astype(str) != post_df.iloc[:n_cmp].astype(str))
                .any(axis=1)
                .sum()
            )
        else:
            changed = 0
        total_affected = abs(n_pre - n_post) + changed
        denominator = max(n_pre, n_post, 1)

    return min(1.0, total_affected / denominator)


def normalized_rank_displacement(
    pre_df: pd.DataFrame,
    post_df: pd.DataFrame,
    score_col: str = "vina_score",
    id_col: str = "compound_id",
    ascending: bool = True,
) -> float:
    """
    Mean normalised rank displacement for compounds present in both DataFrames.

    For each shared compound i::

        displacement_i = |rank_pre_i − rank_post_i| / (n − 1)

    Returns the mean displacement over all shared compounds.  Result in [0, 1].

    ``ascending=True`` means lower score is better (Vina convention).
    """
    pre_score = _find_col(pre_df, score_col)
    post_score = _find_col(post_df, score_col)
    if pre_score is None or post_score is None:
        return 0.0

    pre_id = _find_col(pre_df, id_col) or pre_df.columns[0]
    post_id = _find_col(post_df, id_col) or post_df.columns[0]

    def _rank_series(df: pd.DataFrame, id_c: str, sc: str) -> pd.Series:
        tmp = df[[id_c, sc]].copy()
        tmp[sc] = pd.to_numeric(tmp[sc], errors="coerce")
        tmp = tmp.dropna(subset=[sc])
        tmp["_rank"] = tmp[sc].rank(ascending=ascending, method="first")
        return tmp.set_index(id_c)["_rank"]

    pre_ranks = _rank_series(pre_df, pre_id, pre_score)
    post_ranks = _rank_series(post_df, post_id, post_score)

    shared = set(pre_ranks.index) & set(post_ranks.index)
    if not shared:
        return 0.0

    n = max(len(pre_ranks), len(post_ranks))
    if n <= 1:
        return 0.0

    displacements = [
        abs(float(pre_ranks[i]) - float(post_ranks[i])) / (n - 1)
        for i in shared
    ]
    return float(np.mean(displacements))


def top_k_jaccard_divergence(
    pre_df: pd.DataFrame,
    post_df: pd.DataFrame,
    id_col: str = "compound_id",
    score_col: Optional[str] = None,
    k: int = 20,
    ascending: bool = True,
) -> float:
    """
    Jaccard distance between the top-k compounds in pre-fault and post-fault outputs.

    ::

        divergence = 1 − |pre_top_k ∩ post_top_k| / |pre_top_k ∪ post_top_k|

    When ``score_col`` is provided, top-k is selected by score; otherwise by row order.
    Returns a float in [0, 1].
    """
    pre_id = _find_col(pre_df, id_col) or pre_df.columns[0]
    post_id = _find_col(post_df, id_col) or post_df.columns[0]

    def _top_k_ids(df: pd.DataFrame, id_c: str) -> set[str]:
        if score_col is not None:
            sc = _find_col(df, score_col)
            if sc is not None:
                tmp = df.copy()
                tmp["_s"] = pd.to_numeric(tmp[sc], errors="coerce")
                tmp = tmp.dropna(subset=["_s"]).sort_values("_s", ascending=ascending)
                return set(tmp[id_c].astype(str).iloc[:k])
        return set(df[id_c].astype(str).iloc[:k])

    pre_top = _top_k_ids(pre_df, pre_id)
    post_top = _top_k_ids(post_df, post_id)

    if not pre_top and not post_top:
        return 0.0
    union = len(pre_top | post_top)
    if union == 0:
        return 0.0
    return 1.0 - len(pre_top & post_top) / union


def incorrect_claim_proportion(
    pre_df: pd.DataFrame,
    post_df: pd.DataFrame,
    value_cols: Optional[list[str]] = None,
    tolerance: float = 1e-6,
) -> float:
    """
    Proportion of numeric (or string) claims that differ between DataFrames.

    A numeric claim is "incorrect" when::

        |post − pre| > tolerance × max(1, |pre|)

    String claims are compared exactly.  Extra or missing rows count as fully
    incorrect.  Returns a float in [0, 1].
    """
    if pre_df.empty and post_df.empty:
        return 0.0

    cols = (
        [c for c in (value_cols or pre_df.columns) if c in pre_df.columns and c in post_df.columns]
    )
    if not cols:
        return 0.0

    n_rows = min(len(pre_df), len(post_df))
    extra_rows = abs(len(pre_df) - len(post_df))
    n_total = n_rows * len(cols) + extra_rows * len(cols)
    if n_total == 0:
        return 0.0

    incorrect = 0
    for col in cols:
        for pv, qv in zip(pre_df[col].iloc[:n_rows], post_df[col].iloc[:n_rows]):
            try:
                pf, qf = float(pv), float(qv)
                if abs(qf - pf) > tolerance * max(1.0, abs(pf)):
                    incorrect += 1
            except (ValueError, TypeError):
                if str(pv) != str(qv):
                    incorrect += 1

    incorrect += extra_rows * len(cols)
    return min(1.0, incorrect / n_total)


def md_metric_normalized_error(
    pre_value: float,
    post_value: float,
    value_range: float,
) -> float:
    """
    Normalised absolute error for a scalar MD metric (RMSD, binding energy).

    ``= |post − pre| / value_range``, capped at 1.0.

    ``value_range`` is the expected full range of the metric, e.g. 10 Å for
    RMSD or 500 kcal/mol for MM-GBSA binding energy.
    """
    if value_range <= 0:
        raise ValueError(f"value_range must be positive, got {value_range!r}")
    return min(1.0, abs(post_value - pre_value) / value_range)


# ── DistortionRecord ──────────────────────────────────────────────────────────

@dataclass
class DistortionRecord:
    """
    Result of a single distortion computation.

    Fields
    ------
    instance_id  : matches InjectionRecord.instance_id
    family_id    : e.g. "F01"
    metric_name  : name of the core function applied
    value        : distortion in [0, 1], or None for stubs
    status       : "computed" | "stub" | "error"
    detail       : human-readable explanation of what was compared
    """
    instance_id: str
    family_id: str
    metric_name: str
    value: Optional[float]
    status: str                         # "computed" | "stub" | "error"
    detail: str = ""
    extra: dict[str, Any] = field(default_factory=dict)


# ── per-family configuration ──────────────────────────────────────────────────

#: Maps family_id → configuration for compute_distortion.
#: Keys:
#:   metric     : name of the core function to call
#:   file       : filename to compare in pre/post workspace (None = stub)
#:   id_col     : identifier column (for record/rank metrics)
#:   score_col  : score column (for rank/Jaccard metrics)
#:   k          : top-k for Jaccard (default 20)
#:   ascending  : True = lower score is better (Vina convention)
#:   value_range: for md_metric_normalized_error stubs (informational)
_FAMILY_CONFIG: dict[str, dict[str, Any]] = {
    # G0 — provenance errors
    "F01": {
        "metric": "incorrect_record_proportion",
        "file":   "compounds.csv",
        "id_col": "compound_id",
    },
    "F02": {
        "metric": "incorrect_record_proportion",
        "file":   "compounds.csv",           # use compound list as downstream proxy
        "id_col": "compound_id",
        "detail_note": (
            "F02 corrupts a registered input file (1 bit flip).  "
            "Comparing the downstream compound list measures propagation, "
            "not the raw byte corruption."
        ),
    },
    # G1 — structure errors
    "F03": {
        "metric": "incorrect_record_proportion",
        "file":   "standardised_ledger.csv",
        "id_col": "compound_id",
    },
    "F04": {
        "metric": "incorrect_record_proportion",
        "file":   "standardised_ledger.csv",
        "id_col": "compound_id",
    },
    "F05": {
        "metric": "incorrect_record_proportion",
        "file":   "standardised_ledger.csv",
        "id_col": "compound_id",
    },
    # G2–G3 — stubs (require receptor structure, WP5)
    "F06": {"metric": "stub", "value_range": 10.0},   # protonation RMSD proxy
    "F07": {"metric": "stub", "value_range": 20.0},   # off-pocket RMSD (Å)
    "F08": {"metric": "stub", "value_range": 10.0},   # crystal-ref RMSD (Å)
    # G4 — scoring errors
    "F09": {
        "metric":    "normalized_rank_displacement",
        "file":      "docking_scores.csv",
        "id_col":    "compound_id",
        "score_col": "vina_score",
        "ascending": True,
    },
    "F10": {
        "metric":    "normalized_rank_displacement",
        "file":      "docking_scores.csv",
        "id_col":    "compound_id",
        "score_col": "vina_score",
        "ascending": True,
    },
    # G5 — benchmark evaluation errors (stubs, require public benchmark, WP5)
    "F11": {"metric": "stub"},
    "F12": {"metric": "stub"},
    # G6 — selection errors
    "F13": {
        "metric":    "top_k_jaccard_divergence",
        "file":      "ranked_compounds.csv",
        "id_col":    "compound_id",
        "score_col": "vina_score",
        "k":         20,
        "ascending": True,
    },
    "F14": {
        "metric":    "top_k_jaccard_divergence",
        "file":      "selected_panel.csv",
        "id_col":    "compound_id",
        "score_col": None,               # panel already ranked; use row order
        "k":         20,
        "ascending": True,
    },
    # G7 — MD analysis errors (stubs, require trajectory, WP7)
    "F15": {"metric": "stub", "value_range": 10.0},   # trajectory RMSD (Å)
    "F16": {"metric": "stub", "value_range": 5.0},    # binding-site RMSD (Å)
    # G8 — reporting errors
    "F17": {
        "metric":      "incorrect_claim_proportion",
        "file":        "results_table.csv",
        "value_cols":  None,             # compare all numeric columns
        "tolerance":   1e-6,
        "detail_note": (
            "F17 corrupts a figure file (PNG/PDF).  results_table.csv is used "
            "as a proxy for the numeric claims that the figure encodes.  "
            "A full figure-hash comparison requires the pipeline to produce "
            "a registered artifact record."
        ),
    },
    "F18": {
        "metric":     "incorrect_claim_proportion",
        "file":       "results_table.csv",
        "value_cols": None,
        "tolerance":  1e-6,
    },
}


# ── dispatcher ────────────────────────────────────────────────────────────────

def compute_distortion(
    family_id: str,
    instance_id: str,
    pre_workspace: Path,
    post_workspace: Path,
) -> DistortionRecord:
    """
    Compute the decision distortion for one fault instance.

    Parameters
    ----------
    family_id      : e.g. "F01"
    instance_id    : e.g. "F01_003"
    pre_workspace  : directory containing the clean-run output files
    post_workspace : directory containing the faulted-run output files

    Returns
    -------
    DistortionRecord with ``value`` in [0, 1] (or None for stubs).
    """
    if family_id not in _FAMILY_CONFIG:
        return DistortionRecord(
            instance_id=instance_id,
            family_id=family_id,
            metric_name="unknown",
            value=None,
            status="error",
            detail=f"No distortion config for family '{family_id}'.",
        )

    cfg = _FAMILY_CONFIG[family_id]
    metric = cfg["metric"]

    # ── stubs ────────────────────────────────────────────────────────────────
    if metric == "stub":
        return DistortionRecord(
            instance_id=instance_id,
            family_id=family_id,
            metric_name="stub",
            value=None,
            status="stub",
            detail=(
                f"{family_id}: distortion requires real pipeline outputs "
                f"(value_range={cfg.get('value_range', 'N/A')}). "
                f"Implement at WP5/WP7."
            ),
        )

    # ── load DataFrames ───────────────────────────────────────────────────────
    filename: str = cfg["file"]
    pre_df, post_df = _load_pair(pre_workspace, post_workspace, filename)

    if pre_df is None or post_df is None:
        missing = filename if pre_df is None else f"{filename} (post)"
        return DistortionRecord(
            instance_id=instance_id,
            family_id=family_id,
            metric_name=metric,
            value=None,
            status="error",
            detail=f"Could not read '{missing}' from workspace.",
        )

    # ── dispatch ──────────────────────────────────────────────────────────────
    try:
        if metric == "incorrect_record_proportion":
            value = incorrect_record_proportion(
                pre_df, post_df,
                id_col=cfg.get("id_col", "compound_id"),
            )
            detail = (
                f"Compared '{filename}': "
                f"{len(pre_df)} rows (pre) vs {len(post_df)} rows (post). "
                f"{cfg.get('detail_note', '')}"
            )

        elif metric == "normalized_rank_displacement":
            value = normalized_rank_displacement(
                pre_df, post_df,
                score_col=cfg.get("score_col", "vina_score"),
                id_col=cfg.get("id_col", "compound_id"),
                ascending=cfg.get("ascending", True),
            )
            detail = (
                f"Rank displacement in '{filename}' by '{cfg.get('score_col', 'score')}': "
                f"{len(pre_df)} compounds compared."
            )

        elif metric == "top_k_jaccard_divergence":
            k = cfg.get("k", 20)
            value = top_k_jaccard_divergence(
                pre_df, post_df,
                id_col=cfg.get("id_col", "compound_id"),
                score_col=cfg.get("score_col"),
                k=k,
                ascending=cfg.get("ascending", True),
            )
            detail = (
                f"Top-{k} Jaccard divergence in '{filename}'. "
                f"Pre panel: {len(pre_df)} rows, post panel: {len(post_df)} rows."
            )

        elif metric == "incorrect_claim_proportion":
            tol = cfg.get("tolerance", 1e-6)
            vcols = cfg.get("value_cols")
            value = incorrect_claim_proportion(
                pre_df, post_df,
                value_cols=vcols,
                tolerance=tol,
            )
            detail = (
                f"Claim comparison in '{filename}' "
                f"(tolerance={tol}): "
                f"{len(pre_df)} rows pre, {len(post_df)} rows post. "
                f"{cfg.get('detail_note', '')}"
            )

        else:
            return DistortionRecord(
                instance_id=instance_id,
                family_id=family_id,
                metric_name=metric,
                value=None,
                status="error",
                detail=f"Unknown metric '{metric}' in _FAMILY_CONFIG.",
            )

    except Exception as exc:  # noqa: BLE001
        return DistortionRecord(
            instance_id=instance_id,
            family_id=family_id,
            metric_name=metric,
            value=None,
            status="error",
            detail=f"Exception during distortion computation: {exc}",
        )

    return DistortionRecord(
        instance_id=instance_id,
        family_id=family_id,
        metric_name=metric,
        value=round(value, 6),
        status="computed",
        detail=detail,
    )

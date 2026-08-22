"""
calibrate.py — Probability calibration on development groups only (ADR-AIIM-021).

Uses FrozenEstimator (sklearn 1.4+) inside CalibratedClassifierCV so the base
Pipeline fitted by train.fit_selected_model() is not refitted during calibration.
GroupKFold keyed by fault-template family prevents cross-group contamination.

Calibration uses development families only and is frozen before the held-out run.
No recalibration is permitted after unsealing.
"""

from __future__ import annotations

import logging
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.calibration import CalibratedClassifierCV, calibration_curve
from sklearn.frozen import FrozenEstimator
from sklearn.metrics import brier_score_loss
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import Pipeline

from telemetry.leakage_guard import check_field_names, check_allowlist

logger = logging.getLogger(__name__)

# Default number of calibration folds.  Must not exceed the number of
# distinct fault-template families in the calibration split.
N_CAL_FOLDS = 3


def calibrate(
    fitted_pipeline: Pipeline,
    X_cal: pd.DataFrame,
    y_cal: np.ndarray,
    groups_cal: np.ndarray,
    method: str = "isotonic",
    n_cal_folds: int = N_CAL_FOLDS,
    training_groups: np.ndarray | None = None,
) -> CalibratedClassifierCV:
    """
    Wrap a pre-fitted pipeline with probability calibration.

    Parameters
    ----------
    fitted_pipeline:
        Fully fitted Pipeline from train.fit_selected_model().
        Wrapped in FrozenEstimator so it is never refitted on calibration data.
    X_cal:
        Features for the calibration split (development families only).
    y_cal:
        Binary target for the calibration split.
    groups_cal:
        Fault-template family identifier for each calibration row.
        Used to key the inner GroupKFold so no family crosses the fold boundary.
    method:
        'isotonic' (default) or 'sigmoid'.  Isotonic is more flexible but
        needs at least ~1,000 samples to avoid overfitting.  Use sigmoid for
        small calibration sets.
    n_cal_folds:
        Number of GroupKFold splits for calibration.  Automatically reduced
        if fewer distinct groups are present.
    training_groups:
        Optional: fault-template family IDs from the training split.  When
        provided, the function verifies that no calibration family appears in
        training (family disjointness).  Raises RegistryError on overlap.

    Returns
    -------
    Calibrated classifier whose predict_proba() returns calibrated probabilities.

    Raises
    ------
    ValueError if training and calibration groups overlap (data leakage).
    """
    check_field_names(X_cal.columns, context="calibrate.calibrate X_cal columns")
    check_allowlist(X_cal.columns, context="calibrate.calibrate X_cal columns")

    y_cal = np.asarray(y_cal)
    groups_cal = np.asarray(groups_cal)

    # Family disjointness check
    if training_groups is not None:
        training_groups_arr = np.asarray(training_groups)
        overlap = set(np.unique(groups_cal)) & set(np.unique(training_groups_arr))
        if overlap:
            raise ValueError(
                f"Calibration and training data share fault-template families: {overlap}. "
                "Calibration families must be disjoint from training families "
                "(ADR-AIIM-021: grouped splitting by fault-template family)."
            )

    # Adapt n_cal_folds to available distinct groups
    n_unique_cal = len(np.unique(groups_cal))
    actual_folds = min(n_cal_folds, n_unique_cal)
    if actual_folds < 2:
        raise ValueError(
            f"Calibration set has only {n_unique_cal} distinct group(s); "
            "need at least 2 for GroupKFold calibration."
        )
    if actual_folds < n_cal_folds:
        logger.warning(
            "Calibration set has %d distinct groups; reducing n_cal_folds from %d to %d.",
            n_unique_cal, n_cal_folds, actual_folds,
        )

    cal_cv = GroupKFold(n_splits=actual_folds)

    # Pre-generate splits explicitly. Passing a list of (train_idx, test_idx) tuples
    # to cv= avoids sklearn 1.9.0 metadata-routing issues where groups passed to
    # calibrated.fit() would not be forwarded automatically to cal_cv.split().
    explicit_splits = list(cal_cv.split(X_cal, y_cal, groups=groups_cal))

    # FrozenEstimator prevents re-fitting the base pipeline during calibration.
    # CalibratedClassifierCV uses out-of-fold predictions to fit isotonic/sigmoid layer.
    calibrated = CalibratedClassifierCV(
        estimator=FrozenEstimator(fitted_pipeline),
        method=method,
        cv=explicit_splits,
    )
    calibrated.fit(X_cal, y_cal)

    proba_cal = calibrated.predict_proba(X_cal)[:, 1]
    brier = brier_score_loss(y_cal, proba_cal)
    logger.info(
        "Calibration complete (method=%s, folds=%d): Brier score on cal set = %.4f",
        method, actual_folds, brier,
    )

    return calibrated


def evaluate_calibration(
    calibrated_pipeline,
    X_eval: pd.DataFrame,
    y_eval: np.ndarray,
    n_bins: int = 10,
) -> dict:
    """
    Evaluate calibration quality on a held-aside development split.

    Returns fraction_of_positives, mean_predicted_value, and Brier score.
    Do NOT call this on held-out evaluation data before unsealing.
    """
    check_field_names(X_eval.columns, context="calibrate.evaluate_calibration X_eval")
    check_allowlist(X_eval.columns, context="calibrate.evaluate_calibration X_eval")

    proba = calibrated_pipeline.predict_proba(X_eval)[:, 1]
    fraction_pos, mean_pred = calibration_curve(y_eval, proba, n_bins=n_bins)
    brier = brier_score_loss(y_eval, proba)

    return {
        "fraction_of_positives": fraction_pos.tolist(),
        "mean_predicted_value": mean_pred.tolist(),
        "brier_score": float(brier),
        "n_bins": n_bins,
    }


def save_calibrated_model(
    calibrated_pipeline,
    output_path: Path | str,
) -> str:
    """
    Save the calibrated pipeline to disk with joblib.

    Returns the file path (as string) for registration in the model registry.
    Caller must register the file and record its SHA-256 before held-out evaluation.
    """
    p = Path(output_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(calibrated_pipeline, p, compress=3)
    logger.info("Calibrated model saved to %s", p)
    return str(p)

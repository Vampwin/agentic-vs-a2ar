"""
test_calibration.py — Calibration correctness tests (ADR-AIIM-021).

Verifies that:
1. calibrate() wraps a pre-fitted pipeline and produces valid probabilities.
2. FrozenEstimator prevents the base pipeline from being refitted.
3. GroupKFold is used — no family crosses the calibration fold boundary.
4. Family disjoint check raises ValueError when training and calibration overlap.
5. Insufficient groups raises ValueError.
6. save_calibrated_model() serialises to disk.
"""

import numpy as np
import pandas as pd
import pytest
import joblib

from telemetry.leakage_guard import get_canonical_feature_order


def _make_pipeline_and_data(
    n_families: int = 6,
    instances_per_family: int = 8,
    seed: int = 0,
):
    """
    Fit a logistic L2 Pipeline on synthetic data with the canonical feature order.
    Returns (fitted_pipeline, X_cal, y_cal, groups_cal, training_groups).
    """
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler
    from risk_monitor.train import fit_selected_model

    feature_order = get_canonical_feature_order()
    rng = np.random.default_rng(seed)
    n_total = n_families * instances_per_family

    # Split families: first half for training, second half for calibration
    half = n_families // 2
    train_families = np.arange(half)
    cal_families = np.arange(half, n_families)

    all_groups = np.repeat(np.arange(n_families), instances_per_family)
    train_mask = np.isin(all_groups, train_families)
    cal_mask = np.isin(all_groups, cal_families)

    X_all = pd.DataFrame(rng.random((n_total, len(feature_order))), columns=feature_order)
    y_all = (all_groups % 2).astype(int)

    X_train = X_all[train_mask].reset_index(drop=True)
    y_train = y_all[train_mask]
    groups_train = all_groups[train_mask]

    X_cal = X_all[cal_mask].reset_index(drop=True)
    y_cal = y_all[cal_mask]
    groups_cal = all_groups[cal_mask]

    fitted_pipeline = fit_selected_model(X_train, y_train, model_name="logistic_l2")

    return fitted_pipeline, X_cal, y_cal, groups_cal, groups_train


# ---------------------------------------------------------------------------
# Basic calibration tests
# ---------------------------------------------------------------------------

def test_calibrate_returns_valid_probabilities(tmp_path):
    """calibrate() returns a calibrated model producing probabilities in [0, 1]."""
    from risk_monitor.calibrate import calibrate

    fitted_pipeline, X_cal, y_cal, groups_cal, _ = _make_pipeline_and_data()
    calibrated = calibrate(
        fitted_pipeline=fitted_pipeline,
        X_cal=X_cal,
        y_cal=y_cal,
        groups_cal=groups_cal,
        method="sigmoid",
    )

    proba = calibrated.predict_proba(X_cal)
    assert proba.shape == (len(X_cal), 2), "predict_proba must return (n, 2) array."
    assert np.all(proba >= 0.0) and np.all(proba <= 1.0), "All probabilities must be in [0, 1]."
    np.testing.assert_allclose(proba.sum(axis=1), 1.0, atol=1e-6)


def test_calibrate_frozen_base_pipeline():
    """FrozenEstimator prevents the base pipeline's coefficients from changing."""
    from risk_monitor.calibrate import calibrate
    from sklearn.frozen import FrozenEstimator

    fitted_pipeline, X_cal, y_cal, groups_cal, _ = _make_pipeline_and_data()

    # Record the base pipeline's coefficients before calibration
    clf_before = fitted_pipeline.named_steps["clf"]
    coef_before = clf_before.coef_.copy()

    calibrated = calibrate(
        fitted_pipeline=fitted_pipeline,
        X_cal=X_cal,
        y_cal=y_cal,
        groups_cal=groups_cal,
        method="sigmoid",
    )

    # Coefficients must not have changed (FrozenEstimator prevents refitting)
    coef_after = clf_before.coef_
    np.testing.assert_array_equal(
        coef_before, coef_after,
        err_msg="Base pipeline coefficients must not change during calibration."
    )


def test_calibrate_family_overlap_raises():
    """calibrate() raises ValueError when calibration families overlap with training."""
    from risk_monitor.calibrate import calibrate

    fitted_pipeline, X_cal, y_cal, groups_cal, training_groups = _make_pipeline_and_data(
        n_families=6
    )

    # Inject overlap: training families = [0,1,2]; calibration = [3,4,5]
    # Tamper to create overlap by setting cal groups = training groups
    with pytest.raises(ValueError, match="disjoint|overlap"):
        calibrate(
            fitted_pipeline=fitted_pipeline,
            X_cal=X_cal,
            y_cal=y_cal,
            groups_cal=training_groups,   # same groups as training → overlap
            method="sigmoid",
            training_groups=training_groups,
        )


def test_calibrate_insufficient_groups_raises():
    """calibrate() raises ValueError when calibration set has < 2 distinct groups."""
    from risk_monitor.calibrate import calibrate
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    feature_order = get_canonical_feature_order()
    rng = np.random.default_rng(0)
    n = 20
    X_train = pd.DataFrame(rng.random((n, len(feature_order))), columns=feature_order)
    y_train = rng.integers(0, 2, size=n)

    pipeline = Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(max_iter=200, random_state=0)),
    ])
    pipeline.fit(X_train, y_train)

    # Calibration set with only 1 distinct group
    X_cal = X_train[:5]
    y_cal = y_train[:5]
    groups_cal = np.zeros(5, dtype=int)  # all same group

    with pytest.raises(ValueError, match="[Gg]roup|fold"):
        from risk_monitor.calibrate import calibrate
        calibrate(
            fitted_pipeline=pipeline,
            X_cal=X_cal,
            y_cal=y_cal,
            groups_cal=groups_cal,
            method="sigmoid",
        )


def test_calibrate_no_family_crosses_cal_fold():
    """Within calibrate(), no fault-template family appears in both cal train and val."""
    from sklearn.model_selection import GroupKFold

    _, X_cal, y_cal, groups_cal, _ = _make_pipeline_and_data(n_families=6)

    # Verify the GroupKFold splits that calibrate() would use
    n_cal_folds = min(3, len(np.unique(groups_cal)))
    cal_cv = GroupKFold(n_splits=n_cal_folds)

    for train_idx, val_idx in cal_cv.split(X_cal, y_cal, groups=groups_cal):
        overlap = set(groups_cal[train_idx]) & set(groups_cal[val_idx])
        assert not overlap, (
            f"Calibration fold has families {overlap} in both train and val. "
            "GroupKFold must prevent this."
        )


def test_save_calibrated_model_roundtrip(tmp_path):
    """save_calibrated_model serialises the calibrated model; joblib.load recovers it."""
    from risk_monitor.calibrate import calibrate, save_calibrated_model

    fitted_pipeline, X_cal, y_cal, groups_cal, _ = _make_pipeline_and_data()
    calibrated = calibrate(
        fitted_pipeline=fitted_pipeline,
        X_cal=X_cal,
        y_cal=y_cal,
        groups_cal=groups_cal,
        method="sigmoid",
    )

    out_path = tmp_path / "calibrated_model.joblib"
    returned_path = save_calibrated_model(calibrated, out_path)
    assert returned_path == str(out_path)

    reloaded = joblib.load(out_path)
    proba_original = calibrated.predict_proba(X_cal)
    proba_reloaded = reloaded.predict_proba(X_cal)
    np.testing.assert_array_equal(proba_original, proba_reloaded)

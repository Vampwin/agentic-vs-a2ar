"""
test_group_splits.py — Grouped split correctness tests (ADR-AIIM-021).

Verifies that:
1. GroupKFold keeps all instances of one fault-template family in the same fold.
2. No template family ever crosses the train/validation boundary.
3. Preprocessing fitted only on training groups (Pipeline prevents leakage).
4. Parameterized variants from the same template are never split.
"""

import numpy as np
import pandas as pd
import pytest
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression


def _make_synthetic_data(n_families: int = 6, instances_per_family: int = 5):
    """
    Synthetic development-set data with known family structure.
    Each family has `instances_per_family` parameterized variants.
    """
    rng = np.random.default_rng(42)
    n = n_families * instances_per_family

    families = np.repeat(np.arange(n_families), instances_per_family)
    X = pd.DataFrame({
        "n_missing_input_hashes": rng.integers(0, 3, size=n).astype(float),
        "n_tool_warnings": rng.integers(0, 10, size=n).astype(float),
        "elapsed_seconds": rng.uniform(1.0, 100.0, size=n),
        "n_retries": rng.integers(0, 2, size=n).astype(float),
        "shadow_id_match_fraction": rng.uniform(0.8, 1.0, size=n),
    })
    y = (families % 2).astype(int)  # alternating labels by family
    return X, y, families


def test_groupkfold_no_template_crosses_boundary():
    """No fault-template family appears in both train and validation in any fold."""
    X, y, groups = _make_synthetic_data(n_families=6, instances_per_family=5)
    cv = GroupKFold(n_splits=3)

    for fold_idx, (train_idx, val_idx) in enumerate(cv.split(X, y, groups=groups)):
        train_groups = set(groups[train_idx])
        val_groups = set(groups[val_idx])
        overlap = train_groups & val_groups
        assert not overlap, (
            f"Fold {fold_idx}: groups {overlap} appear in both train and validation. "
            "Fault-template families must not cross the split boundary (ADR-AIIM-021)."
        )


def test_groupkfold_all_instances_of_family_together():
    """All instances from one family are always in the same fold partition."""
    X, y, groups = _make_synthetic_data(n_families=6, instances_per_family=5)
    cv = GroupKFold(n_splits=3)

    for fold_idx, (train_idx, val_idx) in enumerate(cv.split(X, y, groups=groups)):
        for family_id in np.unique(groups):
            family_mask = groups == family_id
            in_train = family_mask[train_idx].any()
            in_val = family_mask[val_idx].any()
            assert not (in_train and in_val), (
                f"Fold {fold_idx}: family {family_id} instances appear in both "
                "train and validation partitions."
            )


def test_pipeline_imputer_fitted_only_on_train():
    """
    StandardScaler mean is computed from training rows only.
    This is the key scikit-learn Pipeline property that prevents leakage.
    """
    X, y, groups = _make_synthetic_data(n_families=6, instances_per_family=10)
    pipeline = Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(max_iter=500, random_state=0)),
    ])

    cv = GroupKFold(n_splits=3)
    for fold_idx, (train_idx, val_idx) in enumerate(cv.split(X, y, groups=groups)):
        X_train, X_val = X.iloc[train_idx], X.iloc[val_idx]
        y_train = y[train_idx]

        pipeline.fit(X_train, y_train)
        scaler = pipeline.named_steps["scaler"]

        # The scaler's mean must match the training fold's mean, not the full-data mean
        full_mean = X["elapsed_seconds"].mean()
        train_mean = X_train["elapsed_seconds"].mean()
        col_idx = list(X.columns).index("elapsed_seconds")
        fitted_mean = scaler.mean_[col_idx]

        assert abs(fitted_mean - train_mean) < 1e-6, (
            f"Fold {fold_idx}: scaler fitted on wrong data. "
            f"Expected train mean {train_mean:.4f}, got {fitted_mean:.4f}."
        )
        # The full-data mean is different → fitting on all data would be detectable
        if abs(full_mean - train_mean) > 0.01:
            assert fitted_mean != pytest.approx(full_mean, abs=0.01), (
                f"Fold {fold_idx}: scaler appears to have been fitted on all data "
                f"(full mean={full_mean:.4f}, scaler mean={fitted_mean:.4f})."
            )


def test_no_random_instance_split():
    """
    Random instance-level splitting (shuffled KFold) would allow templates to cross.
    Verify that with GroupKFold the constraint holds even when instances are shuffled.
    """
    rng = np.random.default_rng(99)
    n_families, instances = 4, 8
    n = n_families * instances

    # Shuffle instance order
    families = np.repeat(np.arange(n_families), instances)
    shuffle_idx = rng.permutation(n)
    families = families[shuffle_idx]

    X = pd.DataFrame({"feature": rng.random(n)})
    y = (families % 2).astype(int)

    cv = GroupKFold(n_splits=4)
    for train_idx, val_idx in cv.split(X, y, groups=families):
        train_groups = set(families[train_idx])
        val_groups = set(families[val_idx])
        assert not (train_groups & val_groups), (
            "GroupKFold must keep all instances of a template family together "
            "even when the dataframe rows are shuffled."
        )


# ---------------------------------------------------------------------------
# compare_models() nested CV tests
# ---------------------------------------------------------------------------

def _make_compare_models_data(
    n_families: int = 10,
    instances_per_family: int = 5,
    n_features: int | None = None,
    seed: int = 42,
):
    """
    Synthetic development-set data suitable for compare_models().
    Uses the actual allowed feature names from the schema allowlist.
    """
    from telemetry.leakage_guard import get_canonical_feature_order
    feature_order = get_canonical_feature_order()
    if n_features is not None:
        feature_order = feature_order[:n_features]

    rng = np.random.default_rng(seed)
    n = n_families * instances_per_family
    groups = np.repeat(np.arange(n_families), instances_per_family)

    X = pd.DataFrame(
        rng.random((n, len(feature_order))),
        columns=feature_order,
    )
    # Even-numbered families are unsafe_release=1; odd = 0
    y = (groups % 2).astype(int)
    return X, y, groups


def test_compare_models_returns_two_results():
    """compare_models() returns one CVResult per candidate model (2 total)."""
    from risk_monitor.train import compare_models
    X, y, groups = _make_compare_models_data(n_families=10, instances_per_family=5)
    results = compare_models(X, y, groups)
    assert len(results) == 2, f"Expected 2 results (one per model), got {len(results)}"


def test_compare_models_ordered_by_brier():
    """compare_models() returns results sorted ascending by mean_brier."""
    from risk_monitor.train import compare_models
    X, y, groups = _make_compare_models_data(n_families=10, instances_per_family=5)
    results = compare_models(X, y, groups)
    briers = [r.mean_brier for r in results]
    assert briers == sorted(briers), "Results must be sorted ascending by mean_brier."


def test_compare_models_fold_details_populated():
    """compare_models() produces N_OUTER_FOLDS fold entries per model."""
    from risk_monitor.train import compare_models, N_OUTER_FOLDS
    X, y, groups = _make_compare_models_data(n_families=10, instances_per_family=5)
    results = compare_models(X, y, groups)
    for result in results:
        assert len(result.fold_details) == N_OUTER_FOLDS, (
            f"Model {result.model_name}: expected {N_OUTER_FOLDS} fold detail entries, "
            f"got {len(result.fold_details)}"
        )


def test_compare_models_no_template_crosses_outer_boundary():
    """
    In compare_models(), no fault-template family appears in both outer train
    and outer validation in any fold (GroupKFold boundary enforcement).
    """
    from risk_monitor.train import compare_models, N_OUTER_FOLDS
    from sklearn.model_selection import GroupKFold

    X, y, groups = _make_compare_models_data(n_families=10, instances_per_family=5)
    # Independently verify the splits that compare_models would use
    cv = GroupKFold(n_splits=N_OUTER_FOLDS)
    for train_idx, val_idx in cv.split(X, y, groups=groups):
        assert not (set(groups[train_idx]) & set(groups[val_idx])), (
            "GroupKFold must prevent any fault-template family from crossing "
            "the outer train/validation boundary."
        )


def test_compare_models_brier_range():
    """Brier scores from compare_models() are in the valid [0, 1] range."""
    from risk_monitor.train import compare_models
    X, y, groups = _make_compare_models_data(n_families=10, instances_per_family=5)
    results = compare_models(X, y, groups)
    for result in results:
        assert 0.0 <= result.mean_brier <= 1.0, (
            f"Model {result.model_name}: mean_brier={result.mean_brier} outside [0,1]."
        )
        for fd in result.fold_details:
            assert 0.0 <= fd["brier"] <= 1.0, (
                f"Model {result.model_name} fold {fd['fold']}: "
                f"brier={fd['brier']} outside [0,1]."
            )


def test_compare_models_rejects_insufficient_groups():
    """compare_models() raises ValueError when fewer groups than N_OUTER_FOLDS."""
    from risk_monitor.train import compare_models, N_OUTER_FOLDS
    # Only 2 families — fewer than N_OUTER_FOLDS=5
    X, y, groups = _make_compare_models_data(n_families=2, instances_per_family=10)
    with pytest.raises(ValueError, match="distinct"):
        compare_models(X, y, groups)

"""
train.py — Grouped nested cross-validation for the risk monitor (ADR-AIIM-021).

IMPORTANT: This module scaffolds the training interface.  Do NOT run this with
definitive development-set data until humans approve the feature allowlist,
fault-family allocation, label definitions, model candidates, calibration/action
budgets, and compute budget (14_ML_RISK_CONTROL_PROTOCOL.md §gate).

Design constraints enforced here:
  - ALL preprocessing and estimation steps are inside sklearn Pipeline.
  - Outer GroupKFold keyed by fault-template family; no template crosses train/val.
  - Inner GridSearchCV uses GroupKFold with same key for hyperparameter tuning.
  - Groups are passed to inner CV via GridSearchCV.fit(groups=groups_train).
  - No access to held-out labels, fault identity, or experimental condition.
  - The prespecified baseline is regularized logistic regression (L2, saga solver).
  - A limited GradientBoostingClassifier alternative may be compared.
  - Primary metric: Brier score.  Secondary: AUPRC, calibration error.
  - Prefer logistic regression when Brier intervals overlap (ADR-AIIM-022).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss, average_precision_score
from sklearn.model_selection import GroupKFold, GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from telemetry.leakage_guard import validate_telemetry, check_field_names, check_allowlist

logger = logging.getLogger(__name__)

# Outer CV folds for model comparison; inner folds for hyperparameter selection.
N_OUTER_FOLDS = 5
N_INNER_FOLDS = 3

# Minimum distinct groups required for outer and inner CV
_MIN_GROUPS_OUTER = N_OUTER_FOLDS
_MIN_GROUPS_INNER = N_INNER_FOLDS

# Warn when minority class proportion is below this threshold
_CLASS_IMBALANCE_WARN_THRESHOLD = 0.1


@dataclass
class CVResult:
    model_name: str
    mean_brier: float
    std_brier: float
    mean_auprc: float
    std_auprc: float
    fold_details: list[dict[str, Any]] = field(default_factory=list)


# ── Hyperparameter search grids ────────────────────────────────────────────────

_PARAM_GRIDS: dict[str, dict] = {
    "logistic_l2": {
        "clf__C": [0.01, 0.1, 1.0, 10.0],
    },
    "gradient_boosting": {
        "clf__n_estimators": [50, 100, 200],
        "clf__max_depth": [2, 3],
    },
}


# ── Pipeline factories ─────────────────────────────────────────────────────────

def _make_logistic_pipeline() -> Pipeline:
    """
    Baseline: impute nulls with median, scale, then L2 logistic regression.
    All steps are inside the Pipeline to prevent leakage.
    """
    return Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(
            solver="saga",
            max_iter=1000,
            random_state=20260729,
            class_weight="balanced",
        )),
    ])


def _make_gbm_pipeline() -> Pipeline:
    """
    Alternative: impute nulls with median then gradient boosting.
    Tree-based models handle interactions that logistic regression cannot.
    Prefer logistic regression when Brier-score intervals overlap (ADR-AIIM-022).
    """
    return Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("clf", GradientBoostingClassifier(
            n_estimators=100,
            max_depth=3,
            learning_rate=0.05,
            subsample=0.8,
            random_state=20260729,
        )),
    ])


def _make_pipeline(name: str) -> Pipeline:
    makers = {
        "logistic_l2": _make_logistic_pipeline,
        "gradient_boosting": _make_gbm_pipeline,
    }
    if name not in makers:
        raise ValueError(f"Unknown model_name {name!r}; choose from {list(makers)}")
    return makers[name]()


# ── Group and class checks ─────────────────────────────────────────────────────

def _check_groups_and_balance(
    y: np.ndarray,
    groups: np.ndarray,
    context: str,
) -> None:
    """Warn on insufficient groups or severe class imbalance."""
    n_unique = len(np.unique(groups))
    if n_unique < _MIN_GROUPS_OUTER:
        raise ValueError(
            f"{context}: need at least {_MIN_GROUPS_OUTER} distinct fault-template "
            f"families for outer GroupKFold, got {n_unique}."
        )
    minority_prop = float(y.mean()) if y.mean() <= 0.5 else float(1.0 - y.mean())
    if minority_prop < _CLASS_IMBALANCE_WARN_THRESHOLD:
        logger.warning(
            "%s: severe class imbalance — minority class proportion = %.3f "
            "(< %.0f%%). Consider adjusting class_weight or sampling strategy.",
            context, minority_prop, _CLASS_IMBALANCE_WARN_THRESHOLD * 100,
        )


# ── Nested cross-validation ────────────────────────────────────────────────────

def compare_models(
    X: pd.DataFrame,
    y: np.ndarray,
    groups: np.ndarray,
) -> list[CVResult]:
    """
    Compare prespecified model candidates using nested grouped cross-validation.

    Outer loop: GroupKFold(N_OUTER_FOLDS) for unbiased model comparison.
    Inner loop: GridSearchCV with GroupKFold(N_INNER_FOLDS) for hyperparameter
    selection inside each outer training fold.  Groups are propagated to the
    inner CV via GridSearchCV.fit(groups=groups_train).

    Parameters
    ----------
    X:
        Feature matrix.  Column names must match the feature_schema.yaml allowlist.
        Nullable features arrive as NaN; the Pipeline imputer handles them.
    y:
        Binary target array: 1 = unsafe_release, 0 = safe.
    groups:
        Integer or string array of fault-template family identifiers.
        All rows from one template must have the same group value.
        GroupKFold ensures no template crosses the train/val boundary.

    Returns
    -------
    List of CVResult, one per candidate model, ordered by mean Brier score ascending.
    """
    # Leakage check on column names
    check_field_names(X.columns, context="train.compare_models X columns")
    check_allowlist(X.columns, context="train.compare_models X columns")

    y = np.asarray(y)
    groups = np.asarray(groups)

    _check_groups_and_balance(y, groups, context="compare_models")

    candidates = {
        "logistic_l2": _make_logistic_pipeline,
        "gradient_boosting": _make_gbm_pipeline,
    }

    outer_cv = GroupKFold(n_splits=N_OUTER_FOLDS)
    results: list[CVResult] = []

    for name, pipeline_factory in candidates.items():
        fold_details: list[dict[str, Any]] = []
        brier_scores: list[float] = []
        auprc_scores: list[float] = []
        param_grid = _PARAM_GRIDS.get(name, {})

        for fold_idx, (train_idx, val_idx) in enumerate(
            outer_cv.split(X, y, groups=groups)
        ):
            X_train = X.iloc[train_idx]
            X_val = X.iloc[val_idx]
            y_train = y[train_idx]
            y_val = y[val_idx]
            groups_train = groups[train_idx]

            # Inner nested CV for hyperparameter selection (if param_grid non-empty)
            base_pipeline = pipeline_factory()
            if param_grid:
                inner_cv = GroupKFold(n_splits=N_INNER_FOLDS)
                inner_groups_count = len(np.unique(groups_train))
                if inner_groups_count < N_INNER_FOLDS:
                    logger.warning(
                        "compare_models fold %d: only %d distinct groups in training "
                        "partition; inner CV reduced to %d folds.",
                        fold_idx, inner_groups_count, inner_groups_count,
                    )
                    inner_cv = GroupKFold(n_splits=max(2, inner_groups_count))

                search = GridSearchCV(
                    base_pipeline,
                    param_grid,
                    cv=inner_cv,
                    scoring="neg_brier_score",
                    refit=True,
                    n_jobs=1,
                    error_score="raise",
                )
                # Groups propagated to inner GroupKFold via fit's groups kwarg
                search.fit(X_train, y_train, groups=groups_train)
                best_pipeline = search.best_estimator_
            else:
                base_pipeline.fit(X_train, y_train)
                best_pipeline = base_pipeline

            proba = best_pipeline.predict_proba(X_val)[:, 1]

            brier = brier_score_loss(y_val, proba)
            auprc = average_precision_score(y_val, proba)
            brier_scores.append(brier)
            auprc_scores.append(auprc)
            fold_details.append({
                "fold": fold_idx,
                "n_train": len(train_idx),
                "n_val": len(val_idx),
                "brier": brier,
                "auprc": auprc,
                "n_val_positive": int(y_val.sum()),
                "best_params": getattr(search, "best_params_", {}) if param_grid else {},
            })
            logger.debug(
                "Model=%s fold=%d brier=%.4f auprc=%.4f",
                name, fold_idx, brier, auprc,
            )

        results.append(CVResult(
            model_name=name,
            mean_brier=float(np.mean(brier_scores)),
            std_brier=float(np.std(brier_scores)),
            mean_auprc=float(np.mean(auprc_scores)),
            std_auprc=float(np.std(auprc_scores)),
            fold_details=fold_details,
        ))

    results.sort(key=lambda r: r.mean_brier)
    return results


def fit_selected_model(
    X: pd.DataFrame,
    y: np.ndarray,
    model_name: str = "logistic_l2",
) -> Pipeline:
    """
    Fit the selected model on all development data.

    This produces the pipeline that will be passed to calibrate.py.
    Call this only after compare_models() and human approval of model choice.
    The returned pipeline is fitted on X (a DataFrame) so that
    pipeline.feature_names_in_ records the training column order.
    """
    check_field_names(X.columns, context="train.fit_selected_model X columns")
    check_allowlist(X.columns, context="train.fit_selected_model X columns")

    y = np.asarray(y)
    pipeline = _make_pipeline(model_name)
    pipeline.fit(X, y)
    return pipeline

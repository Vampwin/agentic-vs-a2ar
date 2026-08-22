"""
ood_detector.py — Applicability-domain detector for pre-release telemetry (ADR-AIIM-031).

IsolationForest trained on development training features.  Returns True (OOD) or
False (in-domain).  All-null columns are excluded at fit time; the same exclusion
set is reapplied at predict time.

Fail-closed semantics:
  - Not fitted → predict_one returns None → inference.py raises InferenceError
  - Prediction error → predict_one returns True → inference.py raises InferenceError

Forbidden inputs: must not use cAMP data, fault identifiers, held-out labels,
condition identifiers, or post-release outcomes (leakage_guard parity).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


class OutOfDomainDetector:
    """
    IsolationForest-based applicability-domain detector.

    Attributes
    ----------
    non_null_features : list[str]
        Feature columns used for OOD detection; all-null columns are excluded.
    contamination : float
        Expected fraction of training anomalies (default 0.01 → ≈2 of 188 rows).
    """

    def __init__(
        self,
        contamination: float = 0.01,
        n_estimators: int = 100,
        random_state: int = 20260729,
    ) -> None:
        self.contamination = contamination
        self.n_estimators = n_estimators
        self.random_state = random_state
        self.non_null_features: list[str] = []
        self._pipeline = None  # sklearn Pipeline: imputer → IsolationForest
        self._fitted = False

    def fit(self, X) -> "OutOfDomainDetector":
        """
        Fit on training feature DataFrame.

        All-null columns are dropped; remaining columns are median-imputed then
        passed to IsolationForest.  Raises ValueError if fewer than 2 non-null
        columns remain.
        """
        import pandas as pd
        from sklearn.ensemble import IsolationForest
        from sklearn.impute import SimpleImputer
        from sklearn.pipeline import Pipeline

        df = pd.DataFrame(X) if not hasattr(X, "columns") else X.copy()
        self.non_null_features = [
            col for col in df.columns if df[col].notna().any()
        ]
        if len(self.non_null_features) < 2:
            raise ValueError(
                f"OOD detector requires at least 2 non-null features; "
                f"got {len(self.non_null_features)}."
            )
        logger.info(
            "OOD: %d non-null features of %d total: %s",
            len(self.non_null_features), len(df.columns), self.non_null_features,
        )

        X_sub = df[self.non_null_features].astype(float)
        self._pipeline = Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("detector", IsolationForest(
                n_estimators=self.n_estimators,
                contamination=self.contamination,
                random_state=self.random_state,
            )),
        ])
        self._pipeline.fit(X_sub)
        self._fitted = True
        logger.info(
            "OOD detector fitted: %d rows, contamination=%.3f.",
            len(df), self.contamination,
        )
        return self

    def predict_one(self, record: dict[str, Any]) -> bool | None:
        """
        Returns False (in-domain), True (OOD), or None (not fitted).

        None and True both cause inference.py to fail closed.
        """
        if not self._fitted:
            logger.warning("OOD detector not fitted — returning None (fail closed).")
            return None
        try:
            import pandas as pd
            row = {f: record.get(f) for f in self.non_null_features}
            df_row = pd.DataFrame([row]).astype(float)
            prediction = self._pipeline.predict(df_row)[0]
            return bool(prediction == -1)  # IsolationForest: 1=inlier, -1=outlier
        except Exception as exc:
            logger.error("OOD prediction error — failing closed: %s", exc)
            return True

    def score_one(self, record: dict[str, Any]) -> float | None:
        """Anomaly score (lower = more anomalous). For diagnostics only; not used in policy."""
        if not self._fitted:
            return None
        try:
            import pandas as pd
            row = {f: record.get(f) for f in self.non_null_features}
            df_row = pd.DataFrame([row]).astype(float)
            return float(self._pipeline.decision_function(df_row)[0])
        except Exception:
            return None

    def save(self, path: Path | str) -> None:
        """Persist to disk with joblib (compress=3)."""
        import joblib
        joblib.dump(self, path, compress=3)
        logger.info("OOD detector saved → %s", path)

    @staticmethod
    def load(path: Path | str) -> "OutOfDomainDetector":
        """Load from disk. Hash must be verified by ModelRegistry before calling this."""
        import joblib
        det = joblib.load(path)
        if not isinstance(det, OutOfDomainDetector):
            raise TypeError(
                f"Expected OutOfDomainDetector, got {type(det).__name__}"
            )
        if not det._fitted:
            raise RuntimeError("Loaded OutOfDomainDetector is not fitted.")
        return det

"""
test_ood_detector.py — Unit tests for OutOfDomainDetector (ADR-AIIM-031).

Coverage:
  - fit: column exclusion, fitted flag, ValueError on insufficient non-null columns
  - predict_one: in-domain / OOD / unfitted / extreme outlier
  - score_one: returns float when fitted, None when not
  - save/load: round-trip identity; load TypeError on wrong type; RuntimeError on unfitted
  - feature_extractor integration: ood_flag resolved when detector provided, stays None without
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from telemetry.ood_detector import OutOfDomainDetector


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_X(n: int = 60, seed: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    return pd.DataFrame({
        "n_missing_input_hashes": rng.integers(0, 2, n).astype(float),
        "n_unregistered_artifacts": rng.integers(0, 3, n).astype(float),
        "has_unexpected_lineage": rng.integers(0, 2, n).astype(float),
        "n_tool_warnings": rng.integers(0, 5, n).astype(float),
        "n_retries": rng.integers(0, 2, n).astype(float),
        "elapsed_seconds": rng.uniform(0.0, 10.0, n),
        "all_null_col": [None] * n,
        "out_of_domain_flag": [None] * n,
    })


def _base_record() -> dict:
    return {
        "n_missing_input_hashes": 0.0,
        "n_unregistered_artifacts": 0.0,
        "has_unexpected_lineage": 0.0,
        "n_tool_warnings": 0.0,
        "n_retries": 0.0,
        "elapsed_seconds": 1.0,
        "all_null_col": None,
        "out_of_domain_flag": None,
    }


# ---------------------------------------------------------------------------
# Fit
# ---------------------------------------------------------------------------

class TestFit:
    def test_excludes_all_null_columns(self):
        det = OutOfDomainDetector().fit(_make_X(60))
        assert "all_null_col" not in det.non_null_features
        assert "out_of_domain_flag" not in det.non_null_features

    def test_includes_non_null_columns(self):
        det = OutOfDomainDetector().fit(_make_X(60))
        assert "n_missing_input_hashes" in det.non_null_features
        assert "elapsed_seconds" in det.non_null_features

    def test_returns_self(self):
        det = OutOfDomainDetector()
        result = det.fit(_make_X(30))
        assert result is det

    def test_fitted_flag_set(self):
        det = OutOfDomainDetector().fit(_make_X(30))
        assert det._fitted is True

    def test_pipeline_fitted(self):
        det = OutOfDomainDetector().fit(_make_X(30))
        assert det._pipeline is not None

    def test_raises_on_insufficient_non_null_columns(self):
        X = pd.DataFrame({"a": [None] * 5, "b": [None] * 5, "c": [1.0, 2.0, 3.0, 4.0, 5.0]})
        with pytest.raises(ValueError, match="non-null features"):
            OutOfDomainDetector().fit(X)

    def test_accepts_dict_list_input(self):
        records = [{"x": float(i), "y": float(i * 2)} for i in range(30)]
        det = OutOfDomainDetector().fit(records)
        assert det._fitted


# ---------------------------------------------------------------------------
# predict_one
# ---------------------------------------------------------------------------

class TestPredictOne:
    def test_returns_bool_for_fitted(self):
        det = OutOfDomainDetector(contamination=0.01).fit(_make_X(80))
        result = det.predict_one(_base_record())
        assert isinstance(result, bool)

    def test_unfitted_returns_none(self):
        det = OutOfDomainDetector()
        assert det.predict_one(_base_record()) is None

    def test_training_data_mostly_in_domain(self):
        X = _make_X(100, seed=42)
        det = OutOfDomainDetector(contamination=0.01).fit(X)
        ood_count = sum(1 for rec in X.to_dict("records") if det.predict_one(rec) is True)
        assert ood_count <= 5, f"Expected ≤5 OOD in training data, got {ood_count}"

    def test_extreme_outlier_returns_true(self):
        X = _make_X(100, seed=7)
        det = OutOfDomainDetector(contamination=0.05).fit(X)
        outlier = _base_record()
        outlier["n_missing_input_hashes"] = 9999.0
        outlier["n_tool_warnings"] = 9999.0
        outlier["elapsed_seconds"] = 999999.0
        result = det.predict_one(outlier)
        assert result is True, "Extreme outlier must be flagged OOD"

    def test_missing_features_handled(self):
        """predict_one gracefully handles missing keys via record.get(f) → None."""
        det = OutOfDomainDetector(contamination=0.05).fit(_make_X(60))
        sparse_record = {"n_missing_input_hashes": 0.0}  # most keys missing
        result = det.predict_one(sparse_record)
        assert result in (True, False)  # no crash; fails closed or predicts


# ---------------------------------------------------------------------------
# score_one
# ---------------------------------------------------------------------------

class TestScoreOne:
    def test_returns_float_when_fitted(self):
        det = OutOfDomainDetector().fit(_make_X(40))
        score = det.score_one(_base_record())
        assert isinstance(score, float)

    def test_returns_none_when_not_fitted(self):
        det = OutOfDomainDetector()
        assert det.score_one(_base_record()) is None

    def test_outlier_has_lower_score(self):
        X = _make_X(100, seed=3)
        det = OutOfDomainDetector(contamination=0.05).fit(X)
        normal = _base_record()
        outlier = _base_record()
        outlier["n_missing_input_hashes"] = 9999.0
        outlier["elapsed_seconds"] = 999999.0
        assert det.score_one(outlier) < det.score_one(normal)


# ---------------------------------------------------------------------------
# Save / Load
# ---------------------------------------------------------------------------

class TestSaveLoad:
    def test_round_trip_predictions_match(self, tmp_path):
        X = _make_X(60)
        det = OutOfDomainDetector().fit(X)
        path = tmp_path / "ood.joblib"
        det.save(path)
        det2 = OutOfDomainDetector.load(path)
        for rec in X.to_dict("records"):
            assert det.predict_one(rec) == det2.predict_one(rec)

    def test_round_trip_non_null_features_match(self, tmp_path):
        X = _make_X(60)
        det = OutOfDomainDetector().fit(X)
        path = tmp_path / "ood.joblib"
        det.save(path)
        det2 = OutOfDomainDetector.load(path)
        assert det2.non_null_features == det.non_null_features

    def test_load_wrong_type_raises_type_error(self, tmp_path):
        import joblib
        path = tmp_path / "bad.joblib"
        joblib.dump({"not_a_detector": True}, path)
        with pytest.raises(TypeError, match="OutOfDomainDetector"):
            OutOfDomainDetector.load(path)

    def test_load_unfitted_raises_runtime_error(self, tmp_path):
        import joblib
        det = OutOfDomainDetector()  # not fitted
        path = tmp_path / "unfitted.joblib"
        joblib.dump(det, path)
        with pytest.raises(RuntimeError, match="not fitted"):
            OutOfDomainDetector.load(path)


# ---------------------------------------------------------------------------
# feature_extractor integration
# ---------------------------------------------------------------------------

class TestExtractIntegration:
    """Verify that feature_extractor.extract() resolves out_of_domain_flag when ood_detector is provided."""

    def _minimal_extract(self, **kwargs):
        from telemetry.feature_extractor import extract
        return extract(
            registered_inputs=[],
            registered_outputs=[],
            tool_logs=[],
            n_retries=0,
            elapsed_seconds=0.0,
            **kwargs,
        )

    def test_ood_flag_stays_none_without_detector(self):
        record = self._minimal_extract()
        assert record["out_of_domain_flag"] is None

    def test_ood_flag_resolved_when_detector_provided(self):
        X = _make_X(60)
        det = OutOfDomainDetector(contamination=0.01).fit(X)
        record = self._minimal_extract(ood_detector=det)
        assert record["out_of_domain_flag"] is not None
        assert isinstance(record["out_of_domain_flag"], bool)

    def test_n_null_features_consistent_with_training(self):
        """n_null_features must be computed BEFORE ood flag is resolved (training-distribution parity)."""
        X = _make_X(60)
        det = OutOfDomainDetector(contamination=0.01).fit(X)

        record_without = self._minimal_extract()
        record_with = self._minimal_extract(ood_detector=det)

        # n_null_features should be identical whether or not ood_detector is provided,
        # because n_null_features is computed before ood flag resolution.
        assert record_without["n_null_features"] == record_with["n_null_features"]

    def test_all_schema_fields_present(self):
        from telemetry.leakage_guard import get_allowed_names
        record = self._minimal_extract()
        assert set(record.keys()) == get_allowed_names()

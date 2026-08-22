"""
Tests for fault_injection/distortion.py.

Covers each core metric function independently, then the dispatcher
(compute_distortion) for every non-stub family.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from fault_injection.distortion import (
    DistortionRecord,
    _FAMILY_CONFIG,
    compute_distortion,
    incorrect_claim_proportion,
    incorrect_record_proportion,
    md_metric_normalized_error,
    normalized_rank_displacement,
    top_k_jaccard_divergence,
)


# ── helpers ───────────────────────────────────────────────────────────────────

def _df(rows: list[dict]) -> pd.DataFrame:
    return pd.DataFrame(rows).astype(str)


def _write_csv(path: Path, rows: list[dict]) -> None:
    df = pd.DataFrame(rows)
    df.to_csv(path, index=False)


_COMPOUNDS = [
    {"compound_id": f"CMP{i:05d}", "SMILES": "CCO", "name": f"c{i}"}
    for i in range(20)
]

_SCORES = [
    {"compound_id": f"CMP{i:05d}", "SMILES": "CCO", "vina_score": str(-5.0 - i * 0.1)}
    for i in range(20)
]

_LEDGER = [
    {"compound_id": f"CMP{i:05d}", "SMILES": "[C@@H](O)(F)Cl", "name": f"c{i}"}
    for i in range(20)
]

_RESULTS = [
    {"condition": "A", "seer_pct": "0.15", "brier_score": "0.12", "n_compounds": "172"},
    {"condition": "B", "seer_pct": "0.35", "brier_score": "0.28", "n_compounds": "172"},
]


# ═══════════════════════════════════════════════════════════════════════════════
# incorrect_record_proportion
# ═══════════════════════════════════════════════════════════════════════════════

class TestIncorrectRecordProportion:

    def test_identical_dfs_return_zero(self):
        df = _df(_COMPOUNDS)
        assert incorrect_record_proportion(df, df.copy()) == 0.0

    def test_all_rows_changed(self):
        pre = _df([{"compound_id": f"A{i}", "SMILES": "CCO"} for i in range(5)])
        post = _df([{"compound_id": f"B{i}", "SMILES": "CCC"} for i in range(5)])
        val = incorrect_record_proportion(pre, post)
        assert val == 1.0

    def test_one_row_added(self):
        pre = _df(_COMPOUNDS[:10])
        post = _df(_COMPOUNDS[:11])   # one extra row (same ID set + 1 new)
        val = incorrect_record_proportion(pre, post, id_col="compound_id")
        assert val > 0.0
        assert val <= 1.0

    def test_one_row_duplicated(self):
        rows = list(_COMPOUNDS)
        rows.append(_COMPOUNDS[0])    # duplicate first row
        pre = _df(_COMPOUNDS)
        post = _df(rows)
        val = incorrect_record_proportion(pre, post, id_col="compound_id")
        assert val > 0.0

    def test_value_changed_in_shared_row(self):
        pre = _df([{"compound_id": "A", "val": "1.0"}])
        post = _df([{"compound_id": "A", "val": "9.9"}])
        val = incorrect_record_proportion(pre, post)
        assert val == 1.0

    def test_both_empty_returns_zero(self):
        assert incorrect_record_proportion(pd.DataFrame(), pd.DataFrame()) == 0.0

    def test_positional_fallback_when_no_id_col(self):
        pre = pd.DataFrame({"x": ["a", "b", "c"]})
        post = pd.DataFrame({"x": ["a", "Z", "c"]})
        val = incorrect_record_proportion(pre, post, id_col="compound_id")
        assert 0.0 < val <= 1.0

    def test_result_capped_at_one(self):
        pre = _df([{"compound_id": "A", "v": "1"}])
        post = _df([
            {"compound_id": "B", "v": "2"},
            {"compound_id": "C", "v": "3"},
            {"compound_id": "D", "v": "4"},
        ])
        val = incorrect_record_proportion(pre, post)
        assert val <= 1.0


# ═══════════════════════════════════════════════════════════════════════════════
# normalized_rank_displacement
# ═══════════════════════════════════════════════════════════════════════════════

class TestNormalizedRankDisplacement:

    def _scores_df(self, scores: list[float], n: int = None) -> pd.DataFrame:
        n = n or len(scores)
        return pd.DataFrame({
            "compound_id": [f"CMP{i:05d}" for i in range(n)],
            "vina_score": scores,
        })

    def test_identical_ranking_returns_zero(self):
        scores = [-5.0, -4.5, -4.0, -3.5, -3.0]
        df = self._scores_df(scores)
        val = normalized_rank_displacement(df, df.copy())
        assert val == pytest.approx(0.0)

    def test_fully_inverted_ranking_near_theoretical_max(self):
        # For a complete rank reversal with n=10, the mean normalised displacement
        # is sum(|n+1-2i| for i=1..n) / (n*(n-1)) = 50/90 = 5/9 ≈ 0.556.
        # This is the theoretical maximum for a symmetric reversal; 1.0 is not
        # achievable when n > 1 (it would require all n compounds to simultaneously
        # hold rank 1 and rank n).
        n = 10
        pre_scores = [-float(i) for i in range(n)]  # 0.0, -1.0, ..., -9.0
        post_scores = list(reversed(pre_scores))      # -9.0, ..., 0.0 per compound
        pre = self._scores_df(pre_scores)
        post = self._scores_df(post_scores)
        val = normalized_rank_displacement(pre, post)
        assert val == pytest.approx(5 / 9, abs=0.01)

    def test_partial_inversion(self):
        pre_scores = [-5.0, -4.0, -3.0, -2.0, -1.0]
        # swap first and last only
        post_scores = [-1.0, -4.0, -3.0, -2.0, -5.0]
        pre = self._scores_df(pre_scores)
        post = self._scores_df(post_scores)
        val = normalized_rank_displacement(pre, post)
        assert 0.0 < val < 1.0

    def test_missing_score_col_returns_zero(self):
        df = pd.DataFrame({"compound_id": ["A"], "other": [1.0]})
        val = normalized_rank_displacement(df, df.copy(), score_col="vina_score")
        assert val == 0.0

    def test_no_shared_compounds_returns_zero(self):
        pre = pd.DataFrame({"compound_id": ["A", "B"], "vina_score": [-1.0, -2.0]})
        post = pd.DataFrame({"compound_id": ["C", "D"], "vina_score": [-3.0, -4.0]})
        val = normalized_rank_displacement(pre, post)
        assert val == 0.0

    def test_result_in_range(self):
        import numpy as np
        rng = np.random.default_rng(42)
        pre_scores = rng.uniform(-10, 0, 50).tolist()
        post_scores = rng.uniform(-10, 0, 50).tolist()
        pre = self._scores_df(pre_scores)
        post = self._scores_df(post_scores)
        val = normalized_rank_displacement(pre, post)
        assert 0.0 <= val <= 1.0


# ═══════════════════════════════════════════════════════════════════════════════
# top_k_jaccard_divergence
# ═══════════════════════════════════════════════════════════════════════════════

class TestTopKJaccardDivergence:

    def _panel(self, ids: list[str], scores: list[float] = None) -> pd.DataFrame:
        if scores is not None:
            return pd.DataFrame({"compound_id": ids, "vina_score": scores})
        return pd.DataFrame({"compound_id": ids})

    def test_identical_panels_return_zero(self):
        df = self._panel([f"CMP{i}" for i in range(20)])
        val = top_k_jaccard_divergence(df, df.copy(), k=10)
        assert val == pytest.approx(0.0)

    def test_completely_different_panels_return_one(self):
        pre = self._panel([f"A{i}" for i in range(10)])
        post = self._panel([f"B{i}" for i in range(10)])
        val = top_k_jaccard_divergence(pre, post, k=10)
        assert val == pytest.approx(1.0)

    def test_partial_overlap(self):
        # 5 shared, 5 different out of k=10
        pre_ids = [f"SHARED{i}" for i in range(5)] + [f"PRE{i}" for i in range(5)]
        post_ids = [f"SHARED{i}" for i in range(5)] + [f"POST{i}" for i in range(5)]
        pre = self._panel(pre_ids)
        post = self._panel(post_ids)
        val = top_k_jaccard_divergence(pre, post, k=10)
        # |A∩B|=5, |A∪B|=15 → Jaccard = 5/15 = 1/3 → distance = 2/3
        assert val == pytest.approx(2 / 3, abs=0.01)

    def test_with_score_sorting(self):
        # pre: A(-1) B(-2) C(-3) D(-4) E(-5) — top3 = A,B,C
        # post: E(-1) D(-2) C(-3) B(-4) A(-5) — top3 = E,D,C
        # |A∩B|={C}, |A∪B|={A,B,C,D,E} → 1/5 → distance = 4/5
        ids = list("ABCDE")
        pre_scores = [-1.0, -2.0, -3.0, -4.0, -5.0]
        post_scores = [-5.0, -4.0, -3.0, -2.0, -1.0]
        pre = self._panel(ids, pre_scores)
        post = self._panel(ids, post_scores)
        val = top_k_jaccard_divergence(pre, post, score_col="vina_score", k=3, ascending=True)
        assert val == pytest.approx(4 / 5, abs=0.01)

    def test_both_empty_returns_zero(self):
        val = top_k_jaccard_divergence(pd.DataFrame({"compound_id": []}),
                                       pd.DataFrame({"compound_id": []}), k=5)
        assert val == 0.0

    def test_result_in_range(self):
        import numpy as np
        rng = np.random.default_rng(7)
        ids_pre = [f"C{i}" for i in rng.choice(40, 20, replace=False)]
        ids_post = [f"C{i}" for i in rng.choice(40, 20, replace=False)]
        pre = self._panel(ids_pre)
        post = self._panel(ids_post)
        val = top_k_jaccard_divergence(pre, post, k=10)
        assert 0.0 <= val <= 1.0


# ═══════════════════════════════════════════════════════════════════════════════
# incorrect_claim_proportion
# ═══════════════════════════════════════════════════════════════════════════════

class TestIncorrectClaimProportion:

    def test_identical_tables_return_zero(self):
        df = pd.DataFrame(_RESULTS)
        assert incorrect_claim_proportion(df, df.copy()) == pytest.approx(0.0)

    def test_one_value_changed(self):
        pre = pd.DataFrame([{"metric": "auc", "value": "0.85"}])
        post = pd.DataFrame([{"metric": "auc", "value": "0.75"}])
        val = incorrect_claim_proportion(pre, post, value_cols=["value"])
        assert val == pytest.approx(1.0)

    def test_tolerance_allows_rounding(self):
        pre = pd.DataFrame([{"v": "0.850000"}])
        post = pd.DataFrame([{"v": "0.850001"}])
        val = incorrect_claim_proportion(pre, post, value_cols=["v"], tolerance=1e-4)
        assert val == pytest.approx(0.0)

    def test_all_values_wrong(self):
        pre = pd.DataFrame([{"a": "1.0", "b": "2.0"}, {"a": "3.0", "b": "4.0"}])
        post = pd.DataFrame([{"a": "9.0", "b": "8.0"}, {"a": "7.0", "b": "6.0"}])
        val = incorrect_claim_proportion(pre, post)
        assert val == pytest.approx(1.0)

    def test_extra_rows_counted_as_incorrect(self):
        pre = pd.DataFrame([{"v": "1.0"}])
        post = pd.DataFrame([{"v": "1.0"}, {"v": "2.0"}])  # one extra row
        val = incorrect_claim_proportion(pre, post, value_cols=["v"])
        assert val > 0.0

    def test_string_comparison(self):
        pre = pd.DataFrame([{"label": "hybrid"}])
        post = pd.DataFrame([{"label": "neither"}])
        val = incorrect_claim_proportion(pre, post, value_cols=["label"])
        assert val == pytest.approx(1.0)

    def test_both_empty_returns_zero(self):
        assert incorrect_claim_proportion(pd.DataFrame(), pd.DataFrame()) == 0.0


# ═══════════════════════════════════════════════════════════════════════════════
# md_metric_normalized_error
# ═══════════════════════════════════════════════════════════════════════════════

class TestMdMetricNormalizedError:

    def test_no_change(self):
        assert md_metric_normalized_error(2.5, 2.5, 10.0) == pytest.approx(0.0)

    def test_full_range_change(self):
        assert md_metric_normalized_error(0.0, 10.0, 10.0) == pytest.approx(1.0)

    def test_capped_at_one(self):
        assert md_metric_normalized_error(0.0, 100.0, 10.0) == pytest.approx(1.0)

    def test_symmetry(self):
        a = md_metric_normalized_error(2.0, 5.0, 10.0)
        b = md_metric_normalized_error(5.0, 2.0, 10.0)
        assert a == pytest.approx(b)

    def test_invalid_range_raises(self):
        with pytest.raises(ValueError, match="positive"):
            md_metric_normalized_error(1.0, 2.0, 0.0)
        with pytest.raises(ValueError, match="positive"):
            md_metric_normalized_error(1.0, 2.0, -1.0)


# ═══════════════════════════════════════════════════════════════════════════════
# compute_distortion dispatcher
# ═══════════════════════════════════════════════════════════════════════════════

_STUB_FAMILIES  = {"F06", "F07", "F08", "F11", "F12", "F15", "F16"}
_ACTIVE_FAMILIES = {fid for fid in _FAMILY_CONFIG if fid not in _STUB_FAMILIES}


class TestComputeDistortionStubs:

    @pytest.mark.parametrize("fid", sorted(_STUB_FAMILIES))
    def test_stub_returns_none_value(self, fid, tmp_path):
        rec = compute_distortion(fid, f"{fid}_001", tmp_path, tmp_path)
        assert rec.status == "stub"
        assert rec.value is None
        assert rec.family_id == fid

    @pytest.mark.parametrize("fid", sorted(_STUB_FAMILIES))
    def test_stub_metric_name(self, fid, tmp_path):
        rec = compute_distortion(fid, f"{fid}_001", tmp_path, tmp_path)
        assert rec.metric_name == "stub"


class TestComputeDistortionUnknownFamily:
    def test_unknown_family_returns_error(self, tmp_path):
        rec = compute_distortion("F99", "F99_001", tmp_path, tmp_path)
        assert rec.status == "error"
        assert rec.value is None


class TestComputeDistortionG0:
    """F01, F02 — incorrect_record_proportion on compounds.csv"""

    def _setup(self, tmp_path: Path, n_pre=20, n_post=20,
               extra_row=False, change_val=False) -> tuple[Path, Path]:
        pre_ws = tmp_path / "pre"; pre_ws.mkdir()
        post_ws = tmp_path / "post"; post_ws.mkdir()
        _write_csv(pre_ws / "compounds.csv", _COMPOUNDS[:n_pre])
        post_rows = list(_COMPOUNDS[:n_post])
        if extra_row:
            post_rows.append({"compound_id": "CMP99999", "SMILES": "CCC", "name": "extra"})
        if change_val:
            post_rows[0] = {**post_rows[0], "SMILES": "CCC"}
        _write_csv(post_ws / "compounds.csv", post_rows)
        return pre_ws, post_ws

    @pytest.mark.parametrize("fid", ["F01", "F02"])
    def test_identical_returns_zero(self, fid, tmp_path):
        pre, post = self._setup(tmp_path)
        rec = compute_distortion(fid, f"{fid}_001", pre, post)
        assert rec.status == "computed"
        assert rec.value == pytest.approx(0.0)

    @pytest.mark.parametrize("fid", ["F01", "F02"])
    def test_extra_row_gives_positive(self, fid, tmp_path):
        pre, post = self._setup(tmp_path, extra_row=True)
        rec = compute_distortion(fid, f"{fid}_001", pre, post)
        assert rec.status == "computed"
        assert rec.value > 0.0

    @pytest.mark.parametrize("fid", ["F01", "F02"])
    def test_missing_file_returns_error(self, fid, tmp_path):
        pre_ws = tmp_path / "pre"; pre_ws.mkdir()
        post_ws = tmp_path / "post"; post_ws.mkdir()
        rec = compute_distortion(fid, f"{fid}_001", pre_ws, post_ws)
        assert rec.status == "error"

    @pytest.mark.parametrize("fid", ["F01", "F02"])
    def test_value_in_range(self, fid, tmp_path):
        pre, post = self._setup(tmp_path, change_val=True)
        rec = compute_distortion(fid, f"{fid}_001", pre, post)
        assert 0.0 <= rec.value <= 1.0


class TestComputeDistortionG1:
    """F03, F04, F05 — incorrect_record_proportion on standardised_ledger.csv"""

    @pytest.mark.parametrize("fid", ["F03", "F04", "F05"])
    def test_identical_returns_zero(self, fid, tmp_path):
        pre_ws = tmp_path / "pre"; pre_ws.mkdir()
        post_ws = tmp_path / "post"; post_ws.mkdir()
        _write_csv(pre_ws / "standardised_ledger.csv", _LEDGER)
        _write_csv(post_ws / "standardised_ledger.csv", _LEDGER)
        rec = compute_distortion(fid, f"{fid}_001", pre_ws, post_ws)
        assert rec.status == "computed"
        assert rec.value == pytest.approx(0.0)

    @pytest.mark.parametrize("fid", ["F03", "F04", "F05"])
    def test_all_smiles_changed(self, fid, tmp_path):
        pre_ws = tmp_path / "pre"; pre_ws.mkdir()
        post_ws = tmp_path / "post"; post_ws.mkdir()
        _write_csv(pre_ws / "standardised_ledger.csv", _LEDGER)
        altered = [{**r, "SMILES": "CCC"} for r in _LEDGER]
        _write_csv(post_ws / "standardised_ledger.csv", altered)
        rec = compute_distortion(fid, f"{fid}_001", pre_ws, post_ws)
        assert rec.status == "computed"
        assert rec.value == pytest.approx(1.0)


class TestComputeDistortionG4:
    """F09, F10 — normalized_rank_displacement on docking_scores.csv"""

    @pytest.mark.parametrize("fid", ["F09", "F10"])
    def test_identical_returns_zero(self, fid, tmp_path):
        pre_ws = tmp_path / "pre"; pre_ws.mkdir()
        post_ws = tmp_path / "post"; post_ws.mkdir()
        _write_csv(pre_ws / "docking_scores.csv", _SCORES)
        _write_csv(post_ws / "docking_scores.csv", _SCORES)
        rec = compute_distortion(fid, f"{fid}_001", pre_ws, post_ws)
        assert rec.status == "computed"
        assert rec.value == pytest.approx(0.0)

    @pytest.mark.parametrize("fid", ["F09", "F10"])
    def test_negated_scores_give_high_displacement(self, fid, tmp_path):
        pre_ws = tmp_path / "pre"; pre_ws.mkdir()
        post_ws = tmp_path / "post"; post_ws.mkdir()
        _write_csv(pre_ws / "docking_scores.csv", _SCORES)
        negated = [{**r, "vina_score": str(-float(r["vina_score"]))} for r in _SCORES]
        _write_csv(post_ws / "docking_scores.csv", negated)
        rec = compute_distortion(fid, f"{fid}_001", pre_ws, post_ws)
        assert rec.status == "computed"
        assert rec.value > 0.5     # fully inverted → near 1.0

    @pytest.mark.parametrize("fid", ["F09", "F10"])
    def test_value_in_range(self, fid, tmp_path):
        pre_ws = tmp_path / "pre"; pre_ws.mkdir()
        post_ws = tmp_path / "post"; post_ws.mkdir()
        _write_csv(pre_ws / "docking_scores.csv", _SCORES)
        shuffled = list(reversed(_SCORES))
        _write_csv(post_ws / "docking_scores.csv", shuffled)
        rec = compute_distortion(fid, f"{fid}_001", pre_ws, post_ws)
        assert 0.0 <= rec.value <= 1.0


class TestComputeDistortionG6:
    """F13, F14 — top_k_jaccard_divergence"""

    _PANEL = [{"compound_id": f"CMP{i:05d}", "SMILES": "CCO", "vina_score": str(-5.0 - i * 0.1)}
              for i in range(30)]

    @pytest.mark.parametrize("fid,fname", [
        ("F13", "ranked_compounds.csv"),
        ("F14", "selected_panel.csv"),
    ])
    def test_identical_returns_zero(self, fid, fname, tmp_path):
        pre_ws = tmp_path / "pre"; pre_ws.mkdir()
        post_ws = tmp_path / "post"; post_ws.mkdir()
        _write_csv(pre_ws / fname, self._PANEL)
        _write_csv(post_ws / fname, self._PANEL)
        rec = compute_distortion(fid, f"{fid}_001", pre_ws, post_ws)
        assert rec.status == "computed"
        assert rec.value == pytest.approx(0.0)

    @pytest.mark.parametrize("fid,fname", [
        ("F13", "ranked_compounds.csv"),
        ("F14", "selected_panel.csv"),
    ])
    def test_different_panels_give_positive(self, fid, fname, tmp_path):
        pre_ws = tmp_path / "pre"; pre_ws.mkdir()
        post_ws = tmp_path / "post"; post_ws.mkdir()
        _write_csv(pre_ws / fname, self._PANEL[:20])
        # completely different IDs in top-20
        alt_panel = [{"compound_id": f"ALT{i:05d}", "SMILES": "CCC", "vina_score": str(-3.0 - i * 0.1)}
                     for i in range(20)]
        _write_csv(post_ws / fname, alt_panel)
        rec = compute_distortion(fid, f"{fid}_001", pre_ws, post_ws)
        assert rec.status == "computed"
        assert rec.value > 0.0


class TestComputeDistortionG8:
    """F17, F18 — incorrect_claim_proportion on results_table.csv"""

    @pytest.mark.parametrize("fid", ["F17", "F18"])
    def test_identical_returns_zero(self, fid, tmp_path):
        pre_ws = tmp_path / "pre"; pre_ws.mkdir()
        post_ws = tmp_path / "post"; post_ws.mkdir()
        _write_csv(pre_ws / "results_table.csv", _RESULTS)
        _write_csv(post_ws / "results_table.csv", _RESULTS)
        rec = compute_distortion(fid, f"{fid}_001", pre_ws, post_ws)
        assert rec.status == "computed"
        assert rec.value == pytest.approx(0.0)

    @pytest.mark.parametrize("fid", ["F17", "F18"])
    def test_changed_value_gives_positive(self, fid, tmp_path):
        pre_ws = tmp_path / "pre"; pre_ws.mkdir()
        post_ws = tmp_path / "post"; post_ws.mkdir()
        _write_csv(pre_ws / "results_table.csv", _RESULTS)
        altered = [dict(r) for r in _RESULTS]
        altered[0]["seer_pct"] = "0.99"   # big change
        _write_csv(post_ws / "results_table.csv", altered)
        rec = compute_distortion(fid, f"{fid}_001", pre_ws, post_ws)
        assert rec.status == "computed"
        assert rec.value > 0.0

    @pytest.mark.parametrize("fid", ["F17", "F18"])
    def test_missing_file_returns_error(self, fid, tmp_path):
        pre_ws = tmp_path / "pre"; pre_ws.mkdir()
        post_ws = tmp_path / "post"; post_ws.mkdir()
        rec = compute_distortion(fid, f"{fid}_001", pre_ws, post_ws)
        assert rec.status == "error"


class TestDistortionRecordFields:
    """DistortionRecord always has the right fields."""

    def test_computed_record_has_value(self, tmp_path):
        pre_ws = tmp_path / "pre"; pre_ws.mkdir()
        post_ws = tmp_path / "post"; post_ws.mkdir()
        _write_csv(pre_ws / "compounds.csv", _COMPOUNDS)
        _write_csv(post_ws / "compounds.csv", _COMPOUNDS)
        rec = compute_distortion("F01", "F01_001", pre_ws, post_ws)
        assert isinstance(rec, DistortionRecord)
        assert rec.instance_id == "F01_001"
        assert rec.family_id == "F01"
        assert isinstance(rec.detail, str)

    def test_stub_record_fields(self, tmp_path):
        rec = compute_distortion("F06", "F06_001", tmp_path, tmp_path)
        assert rec.metric_name == "stub"
        assert rec.status == "stub"
        assert "WP5" in rec.detail or "WP7" in rec.detail

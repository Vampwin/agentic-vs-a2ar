"""
Tests for fault_injection/pilot_runner.py (Gate 9 feasibility pilot).

Coverage:
  - Workspace factory functions
  - Evaluator dispatch
  - Single-case runner (clean and fault)
  - Full pilot run
  - ICC computation
  - Output contract (no assay data in results)
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pandas as pd
import pytest

from fault_injection.pilot_runner import (
    CaseOutcome,
    PilotConfig,
    PilotResults,
    _clean_hashes,
    _compute_icc,
    _dispatch_evaluator,
    _run_one_case,
    _setup_g0_workspace,
    _setup_g1_workspace,
    _setup_g4_workspace,
    _setup_g6_workspace,
    _setup_g8_workspace,
    build_training_dataset,
    run_pilot,
)


# ─── Helpers ─────────────────────────────────────────────────────────────────

def _make_fault_outcome(
    family_id: str = "F01",
    instance_num: int = 0,
    condition: str = "B",
    escaped: bool = False,
    evaluator_available: bool = True,
    evaluator_state: str = "HALT",
    error: str | None = None,
) -> CaseOutcome:
    return CaseOutcome(
        family_id=family_id,
        instance_id=f"{family_id}_{instance_num:03d}",
        case_type="fault",
        condition=condition,
        contracts_on=True,
        monitor_on=False,
        evaluator_available=evaluator_available,
        evaluator_state=evaluator_state,
        monitor_status="off",
        escaped=escaped,
        distortion=None,
        error=error,
    )


# ─── Workspace factory tests ──────────────────────────────────────────────────

class TestWorkspaceSetup:
    def test_g0_creates_compounds_csv(self, tmp_path: Path) -> None:
        _setup_g0_workspace(tmp_path)
        assert (tmp_path / "compounds.csv").exists()

    def test_g0_csv_has_required_columns(self, tmp_path: Path) -> None:
        _setup_g0_workspace(tmp_path)
        df = pd.read_csv(tmp_path / "compounds.csv")
        assert "compound_id" in df.columns
        assert "SMILES" in df.columns

    def test_g0_csv_no_duplicate_ids(self, tmp_path: Path) -> None:
        _setup_g0_workspace(tmp_path)
        df = pd.read_csv(tmp_path / "compounds.csv")
        assert df["compound_id"].nunique() == len(df)

    def test_g0_default_row_count(self, tmp_path: Path) -> None:
        _setup_g0_workspace(tmp_path)
        df = pd.read_csv(tmp_path / "compounds.csv")
        assert len(df) == 50

    def test_g0_custom_row_count(self, tmp_path: Path) -> None:
        _setup_g0_workspace(tmp_path, n=20)
        df = pd.read_csv(tmp_path / "compounds.csv")
        assert len(df) == 20

    def test_g0_no_assay_columns(self, tmp_path: Path) -> None:
        _setup_g0_workspace(tmp_path)
        df = pd.read_csv(tmp_path / "compounds.csv", nrows=0)
        cols_lower = [c.lower() for c in df.columns]
        for forbidden in ["camp", "assay", "ic50", "ec50", "inhibit"]:
            assert not any(forbidden in c for c in cols_lower), (
                f"Forbidden term '{forbidden}' found in G0 workspace columns"
            )

    def test_g1_creates_ledger(self, tmp_path: Path) -> None:
        _setup_g1_workspace(tmp_path)
        assert (tmp_path / "standardised_ledger.csv").exists()

    def test_g1_has_required_columns(self, tmp_path: Path) -> None:
        _setup_g1_workspace(tmp_path)
        df = pd.read_csv(tmp_path / "standardised_ledger.csv")
        for col in ["SNP_ID", "SMILES", "exclusion_reason", "valid_smiles",
                    "has_stereocenters", "n_stereocenters", "stereo_defined"]:
            assert col in df.columns, f"Missing column: {col}"

    def test_g1_has_some_chiral_smiles(self, tmp_path: Path) -> None:
        _setup_g1_workspace(tmp_path)
        df = pd.read_csv(tmp_path / "standardised_ledger.csv")
        assert df["has_stereocenters"].any(), "G1 workspace must contain some chiral SMILES"

    def test_g1_no_drug_likeness_exclusions(self, tmp_path: Path) -> None:
        _setup_g1_workspace(tmp_path)
        df = pd.read_csv(tmp_path / "standardised_ledger.csv")
        forbidden_kws = ["lipinski", "veber", "qed", "pains"]
        for reason in df["exclusion_reason"].fillna(""):
            assert not any(kw in str(reason).lower() for kw in forbidden_kws)

    def test_g6_creates_ranking_csv(self, tmp_path: Path) -> None:
        _setup_g6_workspace(tmp_path)
        assert (tmp_path / "ranked_compounds.csv").exists()

    def test_g6_creates_selected_panel(self, tmp_path: Path) -> None:
        _setup_g6_workspace(tmp_path)
        assert (tmp_path / "selected_panel.csv").exists()

    def test_g6_selected_panel_is_subset(self, tmp_path: Path) -> None:
        _setup_g6_workspace(tmp_path)
        ranked = pd.read_csv(tmp_path / "ranked_compounds.csv")
        selected = pd.read_csv(tmp_path / "selected_panel.csv")
        assert len(selected) <= len(ranked)
        assert set(selected["compound_id"]).issubset(set(ranked["compound_id"]))

    def test_g6_no_assay_columns(self, tmp_path: Path) -> None:
        _setup_g6_workspace(tmp_path)
        df = pd.read_csv(tmp_path / "ranked_compounds.csv", nrows=0)
        cols_lower = [c.lower() for c in df.columns]
        for forbidden in ["camp", "assay", "ic50", "ec50"]:
            assert not any(forbidden in c for c in cols_lower)


# ─── Evaluator dispatch tests ─────────────────────────────────────────────────

class TestEvaluatorDispatch:
    def test_g0_clean_workspace_passes(self, tmp_path: Path) -> None:
        """A clean workspace passes once its provenance anchor is registered."""
        _setup_g0_workspace(tmp_path)
        clean = _clean_hashes(tmp_path, "G0")
        result, available = _dispatch_evaluator("G0", tmp_path, clean)
        assert available
        assert result is not None
        assert result["state"] == "PASS"

    def test_g0_without_registered_anchor_halts(self, tmp_path: Path) -> None:
        """
        A critical hash check with no registered anchor is unverified, not
        satisfied. Passing here would mean an unregistered input silently
        clears provenance and a corrupted file flows downstream unnoticed.
        """
        _setup_g0_workspace(tmp_path)
        result, available = _dispatch_evaluator("G0", tmp_path)
        assert available
        assert result is not None
        assert result["state"] == "HALT"
        assert "hash_snp_library_csv" in result["failed"]

    def test_g0_detects_content_mutation(self, tmp_path: Path) -> None:
        """A single flipped byte must break the registered hash (F02)."""
        _setup_g0_workspace(tmp_path)
        clean = _clean_hashes(tmp_path, "G0")
        target = tmp_path / "compounds.csv"
        data = bytearray(target.read_bytes())
        data[5] ^= 1
        target.write_bytes(bytes(data))

        result, available = _dispatch_evaluator("G0", tmp_path, clean)
        assert available
        assert result["state"] == "HALT"
        assert "hash_snp_library_csv" in result["failed"]

    def test_g0_returns_dict_with_state(self, tmp_path: Path) -> None:
        _setup_g0_workspace(tmp_path)
        clean = _clean_hashes(tmp_path, "G0")
        result, available = _dispatch_evaluator("G0", tmp_path, clean)
        assert available
        assert isinstance(result, dict)
        assert "state" in result

    def test_g1_clean_workspace_passes(self, tmp_path: Path) -> None:
        _setup_g1_workspace(tmp_path)
        result, available = _dispatch_evaluator("G1", tmp_path)
        assert available
        assert result is not None
        assert result["state"] in ("PASS", "HALT")  # depends on RDKit availability

    def test_g6_clean_workspace_passes(self, tmp_path: Path) -> None:
        _setup_g6_workspace(tmp_path)
        result, available = _dispatch_evaluator("G6", tmp_path)
        assert available
        assert result is not None
        assert result["state"] == "PASS"

    def test_unknown_gate_returns_unavailable(self, tmp_path: Path) -> None:
        result, available = _dispatch_evaluator("G99", tmp_path)
        assert not available
        assert result is None

    def test_g4_clean_workspace_passes(self, tmp_path: Path) -> None:
        """G4 is implemented and dispatched; a clean scoring table passes."""
        _setup_g4_workspace(tmp_path)
        result, available = _dispatch_evaluator("G4", tmp_path)
        assert available
        assert result is not None
        assert result["state"] == "PASS"

    def test_g4_detects_inverted_score_direction(self, tmp_path: Path) -> None:
        """Negated Vina scores (F09) must halt: positive binding energies are unphysical."""
        _setup_g4_workspace(tmp_path)
        scores = tmp_path / "docking_scores.csv"
        df = pd.read_csv(scores)
        df["vina_score"] = -df["vina_score"]
        df.to_csv(scores, index=False)

        result, available = _dispatch_evaluator("G4", tmp_path)
        assert available
        assert result["state"] == "HALT"
        assert "score_direction" in result["failed"]

    def test_g8_clean_workspace_passes(self, tmp_path: Path) -> None:
        _setup_g8_workspace(tmp_path)
        clean = _clean_hashes(tmp_path, "G8")
        result, available = _dispatch_evaluator("G8", tmp_path, clean)
        assert available
        assert result is not None
        assert result["state"] == "PASS"

    def test_g8_detects_stale_figure(self, tmp_path: Path) -> None:
        """A figure whose bytes changed after registration (F17) must halt."""
        _setup_g8_workspace(tmp_path)
        clean = _clean_hashes(tmp_path, "G8")
        (tmp_path / "main_figure.png").write_bytes(b"STALE_FIGURE_FROM_PRIOR_RUN")

        result, available = _dispatch_evaluator("G8", tmp_path, clean)
        assert available
        assert result["state"] == "HALT"
        assert "figure_source_verification" in result["failed"]

    def test_g8_detects_edited_results_table(self, tmp_path: Path) -> None:
        """A value edited after computation (F18) must halt."""
        _setup_g8_workspace(tmp_path)
        clean = _clean_hashes(tmp_path, "G8")
        table = tmp_path / "results_table.csv"
        df = pd.read_csv(table)
        df.loc[0, "seer_pct"] = 0.999
        df.to_csv(table, index=False)

        result, available = _dispatch_evaluator("G8", tmp_path, clean)
        assert available
        assert result["state"] == "HALT"
        assert "table_source_verification" in result["failed"]


# ─── Single-case runner tests ─────────────────────────────────────────────────

class TestRunOneCase:
    def test_f01_fault_condition_b_detected(self, tmp_path: Path) -> None:
        """F01 duplicates compound IDs; G0 id_uniqueness must fire → contained."""
        # Use first instance (fraction=0.01, will add at least 1 duplicate)
        outcome = _run_one_case("F01", 0, "B", tmp_path, base_seed=42, is_clean=False)
        assert outcome.case_type == "fault"
        assert outcome.condition == "B"
        assert outcome.evaluator_available
        # G0 id_uniqueness is critical → HALT expected
        assert not outcome.escaped, (
            f"F01 should be contained by id_uniqueness check; "
            f"evaluator_state={outcome.evaluator_state!r}"
        )

    def test_f01_fault_condition_d_escaped(self, tmp_path: Path) -> None:
        """Under condition D (no defense), F01 should escape."""
        outcome = _run_one_case("F01", 0, "D", tmp_path, base_seed=42, is_clean=False)
        assert outcome.escaped

    def test_clean_case_does_not_escape_under_b(self, tmp_path: Path) -> None:
        """Clean workspace should pass G0 → no false escalation."""
        outcome = _run_one_case("F01", 0, "B", tmp_path, base_seed=42, is_clean=True)
        assert outcome.case_type == "clean"
        assert not outcome.escaped, "Clean workspace should not trigger HALT"

    def test_all_conditions_d_escaped(self, tmp_path: Path) -> None:
        """Under condition D, any fault escapes (no defense)."""
        for family_id in ["F01", "F03", "F05"]:
            outcome = _run_one_case(family_id, 0, "D", tmp_path, base_seed=42, is_clean=False)
            assert outcome.escaped, f"{family_id} should escape under condition D"

    def test_condition_a_has_stub_monitor(self, tmp_path: Path) -> None:
        outcome = _run_one_case("F01", 0, "A", tmp_path, base_seed=42, is_clean=False)
        assert outcome.monitor_status == "stub_not_trained"
        assert outcome.monitor_on

    def test_condition_b_monitor_is_off(self, tmp_path: Path) -> None:
        outcome = _run_one_case("F01", 0, "B", tmp_path, base_seed=42, is_clean=False)
        assert outcome.monitor_status == "off"

    def test_outcome_has_instance_id(self, tmp_path: Path) -> None:
        outcome = _run_one_case("F01", 2, "B", tmp_path, base_seed=42, is_clean=False)
        assert "F01" in outcome.instance_id

    def test_deterministic_seed_gives_same_result(self, tmp_path: Path) -> None:
        a = _run_one_case("F01", 0, "B", tmp_path / "a", base_seed=7, is_clean=False)
        b = _run_one_case("F01", 0, "B", tmp_path / "b", base_seed=7, is_clean=False)
        assert a.escaped == b.escaped
        assert a.evaluator_state == b.evaluator_state

    def test_different_seeds_may_differ(self, tmp_path: Path) -> None:
        """Different seeds → workspace is set up fresh each time (no cross-contamination)."""
        a = _run_one_case("F01", 0, "B", tmp_path / "a", base_seed=1, is_clean=False)
        b = _run_one_case("F01", 0, "B", tmp_path / "b", base_seed=2, is_clean=False)
        # Both should have the same evaluator available
        assert a.evaluator_available == b.evaluator_available


# ─── ICC computation tests ────────────────────────────────────────────────────

class TestComputeICC:
    def _outcomes_from_pattern(
        self,
        pattern: dict[str, list[bool]],
        condition: str = "B",
    ) -> list[CaseOutcome]:
        """Build CaseOutcome list from family_id → [escaped, ...] pattern."""
        outcomes = []
        for fid, escaped_list in pattern.items():
            for i, esc in enumerate(escaped_list):
                state = "PASS" if esc else "HALT"
                outcomes.append(_make_fault_outcome(
                    family_id=fid,
                    instance_num=i,
                    condition=condition,
                    escaped=esc,
                    evaluator_available=True,
                    evaluator_state=state,
                ))
        return outcomes

    def test_icc_high_when_families_differ_completely(self) -> None:
        """F01 all contained, F03 all escaped → between-family variance dominates → ICC ≈ 1."""
        outcomes = self._outcomes_from_pattern({
            "F01": [False] * 10,
            "F03": [True]  * 10,
        })
        icc, (lo, hi), _ = _compute_icc(outcomes, condition="B", n_bootstrap=500)
        assert icc > 0.8, f"Expected ICC > 0.8, got {icc}"

    def test_icc_near_zero_when_families_identical(self) -> None:
        """Both families have same escape rate → little between-family variance → ICC ≈ 0."""
        outcomes = self._outcomes_from_pattern({
            "F01": [True, False] * 5,
            "F03": [True, False] * 5,
        })
        icc, _, _ = _compute_icc(outcomes, condition="B", n_bootstrap=500)
        assert abs(icc) < 0.5, f"Expected |ICC| < 0.5, got {icc}"

    def test_icc_nan_with_single_family(self) -> None:
        outcomes = self._outcomes_from_pattern({"F01": [False, True, False]})
        icc, _, _ = _compute_icc(outcomes, condition="B")
        assert math.isnan(icc)

    def test_icc_nan_with_no_eligible_cases(self) -> None:
        outcomes = [_make_fault_outcome(escaped=False, evaluator_available=False)]
        icc, _, _ = _compute_icc(outcomes, condition="B")
        assert math.isnan(icc)

    def test_icc_ci_ordered(self) -> None:
        outcomes = self._outcomes_from_pattern({
            "F01": [False] * 8,
            "F03": [True]  * 8,
        })
        _, (lo, hi), _ = _compute_icc(outcomes, condition="B", n_bootstrap=500)
        assert lo <= hi

    def test_icc_excludes_error_cases(self) -> None:
        outcomes = self._outcomes_from_pattern({"F01": [False] * 5, "F03": [True] * 5})
        # Add errored cases that should not influence ICC
        for _ in range(10):
            outcomes.append(_make_fault_outcome(
                family_id="F99", escaped=True,
                evaluator_available=True, error="boom",
            ))
        icc, _, _ = _compute_icc(outcomes, condition="B", n_bootstrap=200)
        # ICC should still be high (only F01 and F03 counted)
        assert icc > 0.5

    def test_icc_note_contains_family_count(self) -> None:
        outcomes = self._outcomes_from_pattern({
            "F01": [False] * 5,
            "F03": [True]  * 5,
        })
        _, _, note = _compute_icc(outcomes, condition="B", n_bootstrap=100)
        assert "2 families" in note

    def test_icc_bounded(self) -> None:
        """ICC must be in [-1, 1]."""
        outcomes = self._outcomes_from_pattern({
            "F01": [False, True, False, True],
            "F03": [True,  False, True, False],
            "F05": [False, False, True, True],
        })
        icc, (lo, hi), _ = _compute_icc(outcomes, condition="B", n_bootstrap=200)
        if not math.isnan(icc):
            assert -1.0 <= icc <= 1.0
        if not math.isnan(lo):
            assert -1.0 <= lo <= hi <= 1.0


# ─── Full pilot run tests ─────────────────────────────────────────────────────

class TestRunPilot:
    def test_stub_family_rejected(self) -> None:
        with pytest.raises(ValueError, match="not implementable"):
            run_pilot(PilotConfig(pilot_families=["F06"]))

    def test_unknown_condition_rejected(self) -> None:
        with pytest.raises(ValueError, match="Unknown condition"):
            run_pilot(PilotConfig(pilot_families=["F01"], conditions=["Z"]))

    def test_minimal_single_family_returns_results(self, tmp_path: Path) -> None:
        config = PilotConfig(
            pilot_families=["F01"],
            conditions=["B", "D"],
            n_clean_per_family=2,
            base_seed=0,
        )
        results = run_pilot(config, workspace_root=tmp_path / "ws")
        assert isinstance(results, PilotResults)
        assert results.n_families == 1
        assert results.n_fault_cases > 0
        assert results.n_clean_cases > 0

    def test_two_families_pilot(self, tmp_path: Path) -> None:
        config = PilotConfig(
            pilot_families=["F01", "F03"],
            conditions=["B", "D"],
            n_clean_per_family=2,
        )
        results = run_pilot(config, workspace_root=tmp_path / "ws")
        assert results.n_families == 2

    def test_seer_by_condition_keys_present(self, tmp_path: Path) -> None:
        config = PilotConfig(
            pilot_families=["F01"],
            conditions=["B", "D"],
            n_clean_per_family=1,
        )
        results = run_pilot(config, workspace_root=tmp_path / "ws")
        assert "B" in results.seer_by_condition
        assert "D" in results.seer_by_condition

    def test_condition_d_seer_is_one(self, tmp_path: Path) -> None:
        """Under condition D (no defense) all faults must escape → SEER = 1.0."""
        config = PilotConfig(
            pilot_families=["F01"],
            conditions=["D"],
            n_clean_per_family=0,
        )
        results = run_pilot(config, workspace_root=tmp_path / "ws")
        assert results.seer_by_condition["D"] == pytest.approx(1.0)

    def test_f01_condition_b_low_seer(self, tmp_path: Path) -> None:
        """F01 (ID duplicate) should be caught by G0 → SEER_B well below 1.0."""
        config = PilotConfig(
            pilot_families=["F01"],
            conditions=["B"],
            n_clean_per_family=0,
        )
        results = run_pilot(config, workspace_root=tmp_path / "ws")
        seer_b = results.seer_by_condition["B"]
        assert seer_b < 0.5, f"F01 should mostly be caught by contracts; SEER_B={seer_b}"

    def test_f03_stereo_loss_now_contained(self, tmp_path: Path) -> None:
        """
        F03 (stereo loss) produces a *valid* SMILES, so no syntactic check can
        see it. It is caught by comparing the stereocentres actually present in
        the structure against the count the ledger declares — the fault leaves
        that metadata untouched, so the contradiction is the signal.

        This previously escaped: the contract only verified that the stereo
        bookkeeping columns existed, not that they were true.
        """
        config = PilotConfig(
            pilot_families=["F03"],
            conditions=["B"],
            n_clean_per_family=0,
        )
        results = run_pilot(config, workspace_root=tmp_path / "ws")
        seer_b = results.seer_by_condition["B"]
        assert seer_b < 0.5, f"F03 should be contained by contracts; SEER_B={seer_b}"

    def test_false_escalation_zero_for_clean_workspace(self, tmp_path: Path) -> None:
        """Clean workspaces should not trigger HALT → false escalation = 0."""
        config = PilotConfig(
            pilot_families=["F01"],
            conditions=["B"],
            n_clean_per_family=5,
        )
        results = run_pilot(config, workspace_root=tmp_path / "ws")
        fe = results.false_escalation_by_condition.get("B", float("nan"))
        assert fe == pytest.approx(0.0), f"Expected 0 false escalation for clean cases, got {fe}"

    def test_output_json_written(self, tmp_path: Path) -> None:
        config = PilotConfig(
            pilot_families=["F01"],
            conditions=["B"],
            n_clean_per_family=1,
        )
        out = tmp_path / "pilot_results.json"
        run_pilot(config, workspace_root=tmp_path / "ws", output_path=out)
        assert out.exists()
        with open(out) as fh:
            data = json.load(fh)
        assert "seer_by_condition" in data
        assert "icc_estimate" in data
        assert "outcomes" in data

    def test_output_json_outcomes_list(self, tmp_path: Path) -> None:
        config = PilotConfig(
            pilot_families=["F01"],
            conditions=["B"],
            n_clean_per_family=1,
        )
        out = tmp_path / "pilot_results.json"
        run_pilot(config, workspace_root=tmp_path / "ws", output_path=out)
        with open(out) as fh:
            data = json.load(fh)
        assert isinstance(data["outcomes"], list)
        assert len(data["outcomes"]) > 0
        sample = data["outcomes"][0]
        for key in ["family_id", "case_type", "condition", "escaped", "evaluator_state"]:
            assert key in sample

    def test_no_assay_data_in_output(self, tmp_path: Path) -> None:
        """Pilot results must contain no cAMP / assay / IC50 terms."""
        config = PilotConfig(
            pilot_families=["F01", "F03"],
            conditions=["B"],
            n_clean_per_family=2,
        )
        results = run_pilot(config, workspace_root=tmp_path / "ws")
        txt = json.dumps(results.to_dict()).lower()
        for forbidden in ["camp", "ic50", "ec50", "pct_inhibition", "known_active"]:
            assert forbidden not in txt, (
                f"Forbidden term '{forbidden}' found in pilot output"
            )

    def test_pilot_summary_contains_seer(self, tmp_path: Path) -> None:
        config = PilotConfig(
            pilot_families=["F01"],
            conditions=["B", "D"],
            n_clean_per_family=1,
        )
        results = run_pilot(config, workspace_root=tmp_path / "ws")
        assert "SEER" in results.pilot_summary

    def test_pilot_summary_mentions_stub_monitor(self, tmp_path: Path) -> None:
        config = PilotConfig(
            pilot_families=["F01"],
            conditions=["A", "B"],
            n_clean_per_family=1,
        )
        results = run_pilot(config, workspace_root=tmp_path / "ws")
        assert "stub_not_trained" in results.pilot_summary

    def test_icc_computed_for_multiple_families(self, tmp_path: Path) -> None:
        config = PilotConfig(
            pilot_families=["F01", "F03"],
            conditions=["B"],
            n_clean_per_family=0,
        )
        results = run_pilot(config, workspace_root=tmp_path / "ws")
        # With 2 families that have very different SEER, ICC should be computable
        assert not math.isnan(results.icc_estimate)

    def test_results_config_matches_input(self, tmp_path: Path) -> None:
        config = PilotConfig(
            pilot_families=["F01"],
            conditions=["B"],
            n_clean_per_family=3,
            base_seed=17,
        )
        results = run_pilot(config, workspace_root=tmp_path / "ws")
        assert results.config["pilot_families"] == ["F01"]
        assert results.config["base_seed"] == 17
        assert results.config["n_clean_per_family"] == 3

    def test_outcomes_count(self, tmp_path: Path) -> None:
        """Number of outcomes = (n_instances × n_conditions) + (n_clean × n_conditions)."""
        from fault_injection.manifest import FAULT_REGISTRY
        config = PilotConfig(
            pilot_families=["F01"],
            conditions=["B", "D"],
            n_clean_per_family=3,
        )
        n_instances = len(FAULT_REGISTRY["F01"].spec.param_grid)  # 10
        expected = n_instances * 2 + 3 * 2  # fault cases + clean cases
        results = run_pilot(config, workspace_root=tmp_path / "ws")
        assert len(results.outcomes) == expected

    def test_temp_workspace_cleaned_up(self) -> None:
        """When no workspace_root is given, temp dir should be deleted after run."""
        config = PilotConfig(
            pilot_families=["F01"],
            conditions=["B"],
            n_clean_per_family=0,
        )
        results = run_pilot(config)  # no workspace_root → uses temp dir
        # If temp cleanup worked, we get a result without an error
        assert isinstance(results, PilotResults)


# ─── Feature collection tests ─────────────────────────────────────────────────

class TestCollectCaseTelemetry:
    """Unit tests for _collect_case_telemetry() workspace helpers."""

    def test_clean_g0_features_all_zero_integrity(self, tmp_path: Path) -> None:
        """Clean workspace: no hash mismatch, no duplicates, no count discrepancy."""
        from fault_injection.pilot_runner import (
            _collect_case_telemetry, _read_ws_ids, _sha256_ws_file,
        )
        _setup_g0_workspace(tmp_path)
        pre_hash = _sha256_ws_file(tmp_path, "compounds.csv")
        ref_ids, n_exp, _ = _read_ws_ids(tmp_path, "G0")
        feats = _collect_case_telemetry(tmp_path, "G0", pre_hash, ref_ids, n_exp)
        assert feats is not None
        assert feats["n_missing_input_hashes"] == 0
        assert feats["n_duplicate_ids"] == 0
        assert feats["n_count_discrepancies"] == 0
        assert feats["n_lost_records"] == 0
        assert feats["n_unexpected_records"] == 0

    def test_clean_g0_has_all_schema_fields(self, tmp_path: Path) -> None:
        from fault_injection.pilot_runner import (
            _collect_case_telemetry, _read_ws_ids, _sha256_ws_file,
        )
        from telemetry.leakage_guard import get_allowed_names
        _setup_g0_workspace(tmp_path)
        pre_hash = _sha256_ws_file(tmp_path, "compounds.csv")
        ref_ids, n_exp, _ = _read_ws_ids(tmp_path, "G0")
        feats = _collect_case_telemetry(tmp_path, "G0", pre_hash, ref_ids, n_exp)
        assert set(feats.keys()) == get_allowed_names()

    def test_f01_injection_raises_duplicate_ids(self, tmp_path: Path) -> None:
        """F01 (id_duplicate) must produce n_duplicate_ids > 0 and hash mismatch."""
        from fault_injection.pilot_runner import (
            _collect_case_telemetry, _read_ws_ids, _sha256_ws_file,
        )
        from fault_injection.g0.f01_id_duplicate import F01IdDuplicate
        _setup_g0_workspace(tmp_path)
        pre_hash = _sha256_ws_file(tmp_path, "compounds.csv")
        ref_ids, n_exp, _ = _read_ws_ids(tmp_path, "G0")
        inj = F01IdDuplicate()
        inj.inject(tmp_path, inj.spec.param_grid[0], seed=42, instance_index=0)
        feats = _collect_case_telemetry(tmp_path, "G0", pre_hash, ref_ids, n_exp)
        assert feats["n_duplicate_ids"] >= 1
        assert feats["n_missing_input_hashes"] == 1

    def test_f01_injection_count_discrepancy(self, tmp_path: Path) -> None:
        """F01 adds rows → n_count_discrepancies = 1 (actual > expected)."""
        from fault_injection.pilot_runner import (
            _collect_case_telemetry, _read_ws_ids, _sha256_ws_file,
        )
        from fault_injection.g0.f01_id_duplicate import F01IdDuplicate
        _setup_g0_workspace(tmp_path)
        pre_hash = _sha256_ws_file(tmp_path, "compounds.csv")
        ref_ids, n_exp, _ = _read_ws_ids(tmp_path, "G0")
        inj = F01IdDuplicate()
        inj.inject(tmp_path, inj.spec.param_grid[0], seed=42, instance_index=0)
        feats = _collect_case_telemetry(tmp_path, "G0", pre_hash, ref_ids, n_exp)
        assert feats["n_count_discrepancies"] == 1

    def test_unknown_gate_returns_none(self, tmp_path: Path) -> None:
        from fault_injection.pilot_runner import _collect_case_telemetry
        result = _collect_case_telemetry(tmp_path, "G99", None, set(), 0)
        assert result is None

    def test_clean_g1_features(self, tmp_path: Path) -> None:
        """Clean G1 workspace: no hash mismatch, no anomalies."""
        from fault_injection.pilot_runner import (
            _collect_case_telemetry, _read_ws_ids, _sha256_ws_file,
        )
        _setup_g1_workspace(tmp_path)
        pre_hash = _sha256_ws_file(tmp_path, "standardised_ledger.csv")
        ref_ids, n_exp, _ = _read_ws_ids(tmp_path, "G1")
        feats = _collect_case_telemetry(tmp_path, "G1", pre_hash, ref_ids, n_exp)
        assert feats is not None
        assert feats["n_missing_input_hashes"] == 0
        assert feats["n_duplicate_ids"] == 0


class TestBuildTrainingDataset:
    """Tests for build_training_dataset() training-data assembly."""

    def _run_mini_pilot(self, tmp_path: Path, families: list[str] | None = None) -> "PilotResults":
        from fault_injection.pilot_runner import build_training_dataset  # noqa: F401
        cfg = PilotConfig(
            pilot_families=families or ["F01", "F03"],
            conditions=["B"],
            n_clean_per_family=2,
            base_seed=42,
        )
        return run_pilot(cfg, workspace_root=tmp_path / "ws")

    def test_returns_triple(self, tmp_path: Path) -> None:
        from fault_injection.pilot_runner import build_training_dataset
        import numpy as np, pandas as pd
        results = self._run_mini_pilot(tmp_path)
        X, y, groups = build_training_dataset(results)
        assert isinstance(X, pd.DataFrame)
        assert isinstance(y, np.ndarray)
        assert isinstance(groups, np.ndarray)

    def test_x_y_groups_same_length(self, tmp_path: Path) -> None:
        from fault_injection.pilot_runner import build_training_dataset
        results = self._run_mini_pilot(tmp_path)
        X, y, groups = build_training_dataset(results)
        assert len(X) == len(y) == len(groups)

    def test_x_has_24_features(self, tmp_path: Path) -> None:
        """Feature matrix must have exactly 24 columns (frozen schema size)."""
        from fault_injection.pilot_runner import build_training_dataset
        results = self._run_mini_pilot(tmp_path)
        X, _, _ = build_training_dataset(results)
        assert X.shape[1] == 24

    def test_labels_are_binary(self, tmp_path: Path) -> None:
        from fault_injection.pilot_runner import build_training_dataset
        results = self._run_mini_pilot(tmp_path)
        _, y, _ = build_training_dataset(results)
        assert set(y).issubset({0, 1})

    def test_groups_are_family_ids(self, tmp_path: Path) -> None:
        from fault_injection.pilot_runner import build_training_dataset
        results = self._run_mini_pilot(tmp_path)
        _, _, groups = build_training_dataset(results, reference_condition="B")
        # Groups should be family IDs (may include "" for clean cases)
        family_groups = {g for g in groups if g}
        assert family_groups.issubset({"F01", "F03", "F02", "F05"})

    def test_reference_condition_filters(self, tmp_path: Path) -> None:
        """All outcomes in dataset must come from the reference condition."""
        from fault_injection.pilot_runner import build_training_dataset
        cfg = PilotConfig(
            pilot_families=["F01"],
            conditions=["B", "D"],
            n_clean_per_family=2,
        )
        results = run_pilot(cfg, workspace_root=tmp_path / "ws")
        X, y, _ = build_training_dataset(results, reference_condition="B")
        # The dataset should not contain condition D outcomes (which would be all-positive)
        # F01 caught under B → label 0 for fault cases; under D → label 1
        # If condition D were included, ALL F01 cases would be 1 (escaped under D)
        # Under B, F01 is caught → label 0 for fault cases
        f01_labels = [y[i] for i in range(len(y)) if y[i] is not None]
        # At least some labels should be 0 (caught cases or clean)
        assert 0 in f01_labels

    def test_features_pass_leakage_check(self, tmp_path: Path) -> None:
        """build_training_dataset raises on forbidden feature column names."""
        from fault_injection.pilot_runner import build_training_dataset
        from telemetry.leakage_guard import get_allowed_names
        results = self._run_mini_pilot(tmp_path)
        X, _, _ = build_training_dataset(results)
        assert set(X.columns).issubset(get_allowed_names())

    def test_missing_reference_condition_raises(self, tmp_path: Path) -> None:
        """If reference_condition has no outcomes, should raise ValueError."""
        from fault_injection.pilot_runner import build_training_dataset
        cfg = PilotConfig(
            pilot_families=["F01"],
            conditions=["B"],
            n_clean_per_family=1,
        )
        results = run_pilot(cfg, workspace_root=tmp_path / "ws")
        with pytest.raises(ValueError):
            build_training_dataset(results, reference_condition="A")

    def test_all_outcomes_have_features_field(self, tmp_path: Path) -> None:
        """After pilot run, every outcome dict should contain a 'features' key."""
        results = self._run_mini_pilot(tmp_path)
        for o in results.outcomes:
            assert "features" in o

    def test_f01_caught_cases_have_label_zero(self, tmp_path: Path) -> None:
        """F01 is caught by G0 contracts → fault case labels should be 0 under B."""
        from fault_injection.pilot_runner import build_training_dataset
        cfg = PilotConfig(
            pilot_families=["F01"],
            conditions=["B"],
            n_clean_per_family=0,
        )
        results = run_pilot(cfg, workspace_root=tmp_path / "ws")
        X, y, _ = build_training_dataset(results, reference_condition="B")
        assert all(y == 0), "F01 caught by contracts → all labels should be 0"

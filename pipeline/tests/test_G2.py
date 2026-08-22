"""
Known-answer tests for the G2 evaluator — receptor/ligand preparation.

structural_waters_recorded (held-out finding H02, ADR-AIIM-060): this check
was declared in contracts/G2_preparation.yaml but had no implementation at
all in G2_evaluator.py's evaluate() loop, so it silently never ran. These
tests pin the real implementation added in response.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from evaluators.G2_evaluator import evaluate

MINIMAL_CONTRACT = {
    "contract_id": "G2_test",
    "checks": [{"id": "structural_waters_recorded", "severity": "high"}],
}

_PDB_3_WATERS = "\n".join([
    "ATOM      1  N   ALA A   1      11.104  13.207   2.001  1.00 20.00           N",
    "ATOM      2  CA  ALA A   1      12.560  13.207   2.001  1.00 20.00           C",
    "HETATM    3  O   HOH A 201      15.000  10.000   5.000  1.00 30.00           O",
    "HETATM    4  O   HOH A 202      16.000  11.000   6.000  1.00 30.00           O",
    "HETATM    5  O   HOH A 203      17.000  12.000   7.000  1.00 30.00           O",
]) + "\n"


def test_passes_when_actual_count_matches_recorded_decision(tmp_path):
    p = tmp_path / "receptor.pdb"
    p.write_text(_PDB_3_WATERS)
    inputs = {
        "prepared_receptor_4eiy_pdbqt": str(p),
        "structural_water_decision": {"cutoff": "5A", "reference_ligand": "ZM241385",
                                       "n_waters_retained": 3},
    }
    result = evaluate(MINIMAL_CONTRACT, inputs, None)
    assert "structural_waters_recorded" in result["passed"]


def test_fails_when_a_water_was_silently_removed(tmp_path):
    """
    Held-out case H02: a water is stripped from the prepared receptor after
    the retention decision was recorded, without updating the decision.
    Before the fix, this check simply did not exist in the evaluator and the
    fault passed by omission (ADR-AIIM-060). It must now fail.
    """
    lines = _PDB_3_WATERS.splitlines()
    removed = "\n".join(l for l in lines if "HOH A 201" not in l) + "\n"
    p = tmp_path / "receptor.pdb"
    p.write_text(removed)
    inputs = {
        "prepared_receptor_4eiy_pdbqt": str(p),
        "structural_water_decision": {"cutoff": "5A", "reference_ligand": "ZM241385",
                                       "n_waters_retained": 3},  # stale: actual is now 2
    }
    result = evaluate(MINIMAL_CONTRACT, inputs, None)
    assert "structural_waters_recorded" in result["failed"]
    assert result["evidence"]["structural_waters_recorded"]["actual_n_waters_in_receptor"] == 2


def test_fails_closed_when_no_decision_recorded(tmp_path):
    p = tmp_path / "receptor.pdb"
    p.write_text(_PDB_3_WATERS)
    contract = {"contract_id": "G2_test",
                "checks": [{"id": "structural_waters_recorded", "severity": "critical"}]}
    result = evaluate(contract, {"prepared_receptor_4eiy_pdbqt": str(p)}, None)
    assert "structural_waters_recorded" in result["failed"]
    assert result["state"] == "HALT"


# ── water_shell_selection_criterion_registered (ADR-AIIM-073, 2026-08-15) ────
#
# Earlier work chose a water-shell cutoff after observing it admitted more of
# the library, which the same work later flagged as a circularity risk. These
# tests pin the check that makes the underlying rule mechanically enforceable.

_SHELL_CONTRACT = {
    "contract_id": "G2_test",
    "checks": [{"id": "water_shell_selection_criterion_registered",
               "allowed_criterion_bases": ["external_recovery", "retrospective_enrichment",
                                          "redocking_accuracy_on_labeled_benchmark"],
               "severity": "high"}],
}


def test_single_shell_is_not_applicable():
    result = evaluate(_SHELL_CONTRACT, {"water_shell_candidates": ["5A"]}, None)
    assert "water_shell_selection_criterion_registered" in result["passed"]
    assert result["evidence"]["water_shell_selection_criterion_registered"]["note"].startswith(
        "not applicable")


def test_no_shells_evaluated_is_also_not_applicable():
    result = evaluate(_SHELL_CONTRACT, {}, None)
    assert "water_shell_selection_criterion_registered" in result["passed"]


def test_fails_when_no_criterion_recorded():
    result = evaluate(_SHELL_CONTRACT,
                      {"water_shell_candidates": ["3A", "4A", "5A"]}, None)
    assert "water_shell_selection_criterion_registered" in result["failed"]
    assert "no selection criterion recorded" in \
        result["evidence"]["water_shell_selection_criterion_registered"]["error"]


def test_fails_when_criterion_is_admission_rate():
    """The exact failure mode this check exists to catch."""
    result = evaluate(
        _SHELL_CONTRACT,
        {"water_shell_candidates": ["3A", "4A", "5A"],
         "water_shell_selection_criterion": {"basis": "library_admission_rate"},
         "water_shell_criterion_registered_before_selection": True},
        None,
    )
    assert "water_shell_selection_criterion_registered" in result["failed"]
    ev = result["evidence"]["water_shell_selection_criterion_registered"]
    assert "library admission rate must not drive the choice" in ev["error"]


def test_fails_when_criterion_registered_after_selection():
    result = evaluate(
        _SHELL_CONTRACT,
        {"water_shell_candidates": ["3A", "4A", "5A"],
         "water_shell_selection_criterion": {"basis": "external_recovery"},
         "water_shell_criterion_registered_before_selection": False},
        None,
    )
    assert "water_shell_selection_criterion_registered" in result["failed"]
    assert "not confirmed registered before" in \
        result["evidence"]["water_shell_selection_criterion_registered"]["error"]


def test_passes_with_a_pre_registered_external_recovery_criterion():
    result = evaluate(
        _SHELL_CONTRACT,
        {"water_shell_candidates": ["3A", "4A", "5A"],
         "water_shell_selection_criterion": {"basis": "retrospective_enrichment",
                                             "detail": "shell maximizing benchmark AUC"},
         "water_shell_criterion_registered_before_selection": True},
        None,
    )
    assert "water_shell_selection_criterion_registered" in result["passed"]

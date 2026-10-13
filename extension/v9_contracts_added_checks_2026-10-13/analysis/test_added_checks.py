"""Tests for the four added checks (ADR-AIIM-427 to 433) and the contracts and evaluator wrappers that carry them. The frozen evaluators are not
edited; these tests exercise the new module `evaluators.added_checks`, the three `*_added_checks` evaluators and the three `*_added_checks.yaml` contracts."""

from __future__ import annotations

import importlib
from pathlib import Path

import sys

import pandas as pd
import pytest
import yaml

TREE = Path(r"D:/Agentic VS Workflow/10_new_system")
sys.path.insert(0, str(TREE))
from evaluators import added_checks as ac      # noqa: E402

CONTRACTS = TREE / "contracts"


def _his_lines(state: str, resseq: int = 77) -> list[str]:
    """Atom lines of one histidine in PDBQT column layout; state is HID (HD1), HIE (HE2) or NONE (no ring hydrogen)."""
    atoms = [("N", "N"), ("CA", "C"), ("CG", "A"), ("ND1", "N" if state == "HID" else "NA"), ("CD2", "A"), ("CE1", "A"),
             ("NE2", "N" if state == "HIE" else "NA")]
    if state == "HID":
        atoms.append(("HD1", "HD"))
    if state == "HIE":
        atoms.append(("HE2", "HD"))
    out = []
    for i, (name, typ) in enumerate(atoms, 1):
        nm = f" {name:<3s}" if len(name) < 4 else name
        out.append(f"ATOM  {i:5d} {nm} HIS A{resseq:4d}    {1.0:8.3f}{2.0:8.3f}{3.0:8.3f}  1.00  0.00    +0.000 {typ:<2s}")
    return out


def test_protonation_matches_record(tmp_path):
    rec = tmp_path / "r.pdbqt"
    rec.write_text("\n".join(_his_lines("HID")))
    assert ac.protonation_matches_record(rec, {"A  77 ": "HID"})[0]
    assert not ac.protonation_matches_record(rec, {"A  77 ": "HIE"})[0]
    assert not ac.protonation_matches_record(rec, None)[0]                     # no record fails closed
    assert not ac.protonation_matches_record(rec, {"A  77 ": "HID", "A  99 ": "HID"})[0]   # a record naming a residue the file lacks
    assert not ac.protonation_matches_record(tmp_path / "missing.pdbqt", {"A  77 ": "HID"})[0]


def test_trajectory_covers_plan(tmp_path):
    rows = [{"compound_id": "C", "replicate": "r1", "time_ns": t} for t in (0.0, 25.0, 49.9)]
    p = tmp_path / "eq.csv"
    pd.DataFrame(rows).to_csv(p, index=False)
    assert ac.trajectory_covers_plan(p, 50.0)[0]
    pd.DataFrame(rows[:2]).to_csv(p, index=False)
    assert not ac.trajectory_covers_plan(p, 50.0)[0]
    assert not ac.trajectory_covers_plan(tmp_path / "missing.csv", 50.0)[0]


class _Registry:
    def __init__(self, log):
        self._log = log

    def get_file_access_log(self):
        return self._log


def test_benchmark_ids_not_exposed(tmp_path):
    bench = tmp_path / "bench.csv"
    pd.DataFrame({"compound_id": ["A1", "A2", "A3"], "label": [1, 0, 0]}).to_csv(bench, index=False)
    scores = tmp_path / "scores.csv"
    scores.write_text("compound_id,score\nA1,1\nA2,2\nA3,3\n")
    leak = tmp_path / "leak.txt"
    leak.write_text("A1\nA9")
    other = tmp_path / "other.txt"
    other.write_text("X1\nX2")
    log = lambda *paths: _Registry([{"reader": "run_ranking", "file_path": str(p), "timestamp": 1.0} for p in paths])
    assert ac.benchmark_ids_not_exposed(log(scores), bench, [str(scores)], 10.0)[0]
    assert ac.benchmark_ids_not_exposed(log(scores, other), bench, [str(scores)], 10.0)[0]          # unrelated IDs pass
    assert not ac.benchmark_ids_not_exposed(log(scores, leak), bench, [str(scores)], 10.0)[0]       # a benchmark ID in a file read before the gate
    assert not ac.benchmark_ids_not_exposed(log(scores, tmp_path / "gone.txt"), bench, [str(scores)], 10.0)[0]   # unreadable fails closed
    late = _Registry([{"reader": "G5b", "file_path": str(bench), "timestamp": 20.0}])
    assert ac.benchmark_ids_not_exposed(late, bench, [str(scores)], 10.0)[0]                       # reads at the gate are the gate's own
    assert not ac.benchmark_ids_not_exposed(None, bench, [], 10.0)[0]


def test_merge_the_worse_state_wins():
    base = {"state": "PASS", "passed": ["a"], "failed": [], "warned": [], "critical_failures": [], "evidence": {}}
    out = ac.merge(base, [("x", False, "critical", {})])
    assert out["state"] == "HALT" and out["failed"] == ["x"] and out["critical_failures"] == ["x"]
    assert ac.merge(base, [("x", True, "critical", {})])["state"] == "PASS"
    assert ac.merge({**base, "state": "ESCALATE"}, [("x", False, "critical", {})])["state"] == "HALT"
    assert ac.merge({**base, "state": "HALT"}, [("x", True, "critical", {})])["state"] == "HALT"


@pytest.mark.parametrize("stem,evaluator,added", [
    ("G2_preparation", "G2_evaluator", ["protonation_matches_record"]),
    ("G7_md", "G7_evaluator", ["trajectory_covers_plan", "binding_site_matches_receptor"]),
    ("G5a_validity", "G5a_evaluator", ["benchmark_ids_not_exposed"]),
])
def test_added_contract_extends_the_frozen_contract(stem, evaluator, added):
    frozen = yaml.safe_load((CONTRACTS / f"{stem}.yaml").read_text(encoding="utf-8"))
    new = yaml.safe_load((CONTRACTS / f"{stem}_added_checks.yaml").read_text(encoding="utf-8"))
    assert new["contract_id"] == frozen["contract_id"] + "_added_checks"
    assert [c["id"] for c in new["checks"]][: len(frozen["checks"])] == [c["id"] for c in frozen["checks"]]
    assert [c["id"] for c in new["checks"]][len(frozen["checks"]):] == added
    assert all(c["severity"] == "critical" for c in new["checks"][len(frozen["checks"]):])
    assert new["evaluator"] == f"evaluators.{evaluator}_added_checks" and frozen["evaluator"] == f"evaluators.{evaluator}"
    mod = importlib.import_module(new["evaluator"])
    assert callable(mod.evaluate)


def test_g7_skip_by_protocol_is_kept():
    from evaluators import G7_evaluator_added_checks as g7
    r = g7.evaluate({"checks": []}, {"md_decision_recorded": "no_MD"}, registry=None)
    assert r["state"] == "SKIPPED_BY_PROTOCOL"

"""ADR-AIIM-431 (pre-registration) / ADR-AIIM-432 (results): the 37 clean and the 160 faulted workspaces of the run v9_added_checks_2026-10-13, evaluated
with the contracts and evaluator modules that now stand in 10_new_system (`*_added_checks.yaml`, `evaluators.*_added_checks`) instead of the copies of
that run, directly and deployed. The harness, the fault injectors and the grading are those of that run; only the contract directory, the contract
file names and the evaluator modules are replaced. Output: results.json, anchor_clean.json and equivalence.json in this folder."""
from __future__ import annotations

import hashlib
import importlib
import json
import sys
from pathlib import Path

OUT = Path(__file__).resolve().parent.parent
PREV = OUT.parent / "v9_added_checks_2026-10-13"
sys.path.insert(0, str(PREV / "analysis"))
import harness as H      # noqa: E402

TREE = H.FROZEN
H.OUT = OUT                                           # results and anchor are written here, read from here
H.PLUS = TREE / "contracts"                           # the deployed path loads the contracts that stand in the tree
H.CONTRACT_FILE = {"G2": "G2_preparation_added_checks.yaml", "G7": "G7_md_added_checks.yaml", "G5a": "G5a_validity_added_checks.yaml",
                   "G5b": "G5b_performance.yaml"}
H.PLUS_MODULE = {"G2": "evaluators.G2_evaluator_added_checks", "G7": "evaluators.G7_evaluator_added_checks", "G5a": "evaluators.G5a_evaluator_added_checks"}
H.write_plus_contracts = lambda: None                 # nothing is written into the tree by this run
import v9ac      # noqa: E402
v9ac.G5a_plus = importlib.import_module(H.PLUS_MODULE["G5a"])
sys.modules["v9ac.G5a_plus"] = v9ac.G5a_plus

import anchor_clean      # noqa: E402
import run_added         # noqa: E402

anchor_clean.main()
run_added.main()

# ---- equivalence with the earlier run
def key(c):
    st = (lambda v: v["G5a"]["state"] if c["gate"] == "G5" else v["state"])
    fl = (lambda v: sorted(v["G5a"]["checks_failed"] if c["gate"] == "G5" else v["checks_failed"]))
    return (c.get("set"), c["family"], json.dumps(c["params"], sort_keys=True, default=str), st(c["frozen"]), fl(c["frozen"]), st(c["plus"]), fl(c["plus"]),
            st(c["deployed_plus"]), fl(c["deployed_plus"]))


old = json.loads((PREV / "results.json").read_text("utf-8"))
new = json.loads((OUT / "results.json").read_text("utf-8"))
same_cases = [key(a) for a in old["cases"]] == [key(b) for b in new["cases"]]
oa = json.loads((PREV / "anchor_clean.json").read_text("utf-8"))
na = json.loads((OUT / "anchor_clean.json").read_text("utf-8"))


def akey(r):
    st = (lambda v: v["G5a"]["state"] if r["gate"] == "G5" else v["state"])
    return (r["gate"], r["case"], st(r["frozen"]), st(r["plus"]), st(r["deployed_plus"]))


same_anchor = [akey(r) for r in oa["rows"]] == [akey(r) for r in na["rows"]]
before = json.loads((OUT / "tree_before.json").read_text("utf-8"))
after = {str(p.relative_to(TREE)).replace("\\", "/"): hashlib.sha256(p.read_bytes()).hexdigest()
         for p in sorted(TREE.rglob("*")) if p.is_file() and "__pycache__" not in p.parts}
changed = [k for k in before if after.get(k) != before[k]]
added = sorted(set(after) - set(before))
eq = {"n_cases": len(new["cases"]), "cases_equal_to_earlier_run": same_cases, "n_clean": len(na["rows"]), "clean_equal_to_earlier_run": same_anchor,
      "predictions_met_now": new["n_met"], "predictions_met_before": old["n_met"], "existing_files_in_tree": len(before),
      "existing_files_changed": changed, "files_added_to_tree": added}
(OUT / "equivalence.json").write_text(json.dumps(eq, indent=1), encoding="utf-8")
print(json.dumps(eq, indent=1))

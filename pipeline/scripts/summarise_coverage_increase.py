"""summarise_coverage_increase.py -- the RA1 increase, measured before and after.

ADR-AIIM-198 reported that `kdss_audit/contract_engine.py` and
`kdss_audit/orchestrator.py` executed no branch under the 682-case suite.
ADR-AIIM-199 records what happened when tests were written for them.

This script does not run the tests. It reads the two `coverage.py` JSON reports
produced by the commands recorded in the run directory's README and emits the
comparison, so that the reported deltas are derived from artifacts rather than
retyped from a terminal.

Both reports were produced with coverage.py 7.16.0 -- the same version the
external audit used -- so the two measurements are comparable to each other.
They are *not* directly comparable to the audit's aggregate, because the audit
ran against a clean checkout of the public repository and this runs against the
development tree, which carries modules the public checkout does not: 1,732
branches here against the audit's 1,578. The per-module figures for the two
modules under test are identical in both (0 of 12, 0 of 22), which is what the
claim rests on; the aggregates are reported on their own basis and labelled.

Output -> 20_new_runs/branch_coverage_increase_2026-09-21/coverage_increase.json
"""
import io
import json
import os
import sys
from datetime import datetime, timezone

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

HERE = os.path.dirname(os.path.abspath(__file__))
SYS = os.path.dirname(HERE)
ROOT = os.path.dirname(SYS)
OUT = os.path.join(ROOT, "20_new_runs", "branch_coverage_increase_2026-09-21")

BEFORE = os.path.join(OUT, "coverage_before.json")
AFTER = os.path.join(OUT, "coverage_after.json")

FOCUS = ["kdss_audit\\contract_engine.py", "kdss_audit\\orchestrator.py"]

# Branches that no input can reach, established in ADR-AIIM-199 and pinned by
# the two tests named here. Recorded so the shortfall is not mistaken for an
# untested branch that somebody could close later.
UNREACHABLE = {
    "kdss_audit/orchestrator.py": [{
        "arc": [169, -160],
        "line": 169,
        "source": "if result.state != AuditState.PASS:",
        "why": ("the enclosing branch is entered only when result.state is "
                "AUTO_REPAIR, so the comparison cannot evaluate false; the "
                "arc is a tautology's false exit"),
        "pinned_by": "tests/test_orchestrator.py::"
                     "test_auto_repair_branch_is_unreachable_from_the_engine",
    }],
}


def load(p):
    with io.open(p, encoding="utf-8") as f:
        return json.load(f)


def summarise(d, key):
    s = d["files"][key]["summary"]
    return {
        "statements": s["num_statements"],
        "percent_statements": round(s["percent_covered"], 1),
        "branches": s["num_branches"],
        "covered_branches": s["covered_branches"],
        "percent_branches": round(s["percent_branches_covered"], 1),
        "missing_branches": d["files"][key]["missing_branches"],
    }


def main():
    for p in (BEFORE, AFTER):
        if not os.path.exists(p):
            print(f"missing artifact: {p}")
            print("see the run directory's README.md for the commands that "
                  "produce it")
            raise SystemExit(1)

    b, a = load(BEFORE), load(AFTER)
    modules = {}
    for key in FOCUS:
        rel = key.replace("\\", "/")
        sb, sa = summarise(b, key), summarise(a, key)
        unreach = UNREACHABLE.get(rel, [])
        reachable = sa["branches"] - len(unreach)
        modules[rel] = {
            "before": sb,
            "after": sa,
            "branches_gained": sa["covered_branches"] - sb["covered_branches"],
            "unreachable_branches": unreach,
            "reachable_branches": reachable,
            "reachable_branch_coverage_percent": (
                round(100.0 * sa["covered_branches"] / reachable, 1)
                if reachable else None),
            "all_reachable_branches_covered":
                sa["covered_branches"] == reachable,
        }

    tb, ta = b["totals"], a["totals"]
    out = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "coverage_tool": f"coverage.py {a['meta']['version']}",
        "branch_coverage_enabled": a["meta"]["branch_coverage"],
        "basis": ("development tree, not the public clean checkout the external "
                  "audit measured; aggregates are not comparable to the audit's, "
                  "per-module figures are"),
        "tests_before": 682,
        "tests_after": 715,
        "tests_added": 33,
        "core_totals": {
            "branches": ta["num_branches"],
            "covered_before": tb["covered_branches"],
            "covered_after": ta["covered_branches"],
            "percent_branches_before": round(
                100.0 * tb["covered_branches"] / tb["num_branches"], 1),
            "percent_branches_after": round(
                100.0 * ta["covered_branches"] / ta["num_branches"], 1),
            "percent_combined_before": int(tb["percent_covered_display"]),
            "percent_combined_after": int(ta["percent_covered_display"]),
        },
        "modules": modules,
        "defects_found_by_writing_the_tests": [
            {
                "id": "D1",
                "module": "kdss_audit/orchestrator.py",
                "line": 166,
                "kind": "unreachable branch in production",
                "statement": ("`_handle_outcome`'s AUTO_REPAIR branch cannot be "
                              "reached through ContractEngine, which transitions "
                              "out of AUTO_REPAIR on every path before "
                              "returning"),
                "patched": False,
            },
            {
                "id": "D2",
                "module": "kdss_audit/orchestrator.py",
                "line": 169,
                "kind": "tautological guard",
                "statement": ("`if result.state != AuditState.PASS` inside the "
                              "AUTO_REPAIR branch can never be false"),
                "patched": False,
            },
            {
                "id": "D3",
                "module": "kdss_audit/orchestrator.py",
                "line": 162,
                "kind": "fail-open default",
                "statement": ("`_handle_outcome` is an if/elif chain with no "
                              "else: an unrecognised state sets neither "
                              "_halted nor _escalation_pending and logs "
                              "nothing, so the workflow continues as though "
                              "the gate had passed"),
                "patched": False,
            },
        ],
        "production_code_modified": False,
    }

    os.makedirs(OUT, exist_ok=True)
    with io.open(os.path.join(OUT, "coverage_increase.json"), "w",
                 encoding="utf-8", newline="\n") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
        f.write("\n")

    print("=" * 72)
    print("RA1 INCREASE -- BRANCH COVERAGE BEFORE AND AFTER")
    print("=" * 72)
    print(f"tool  : {out['coverage_tool']}")
    print(f"tests : {out['tests_before']} -> {out['tests_after']} "
          f"(+{out['tests_added']})")
    print()
    for rel, m in modules.items():
        print(f"{rel}")
        print(f"   branches   {m['before']['covered_branches']:>3}/"
              f"{m['before']['branches']:<3} "
              f"({m['before']['percent_branches']:5.1f} %)  ->  "
              f"{m['after']['covered_branches']:>3}/{m['after']['branches']:<3} "
              f"({m['after']['percent_branches']:5.1f} %)")
        print(f"   statements {m['before']['percent_statements']:5.1f} %"
              f"           ->  {m['after']['percent_statements']:5.1f} %")
        if m["unreachable_branches"]:
            for u in m["unreachable_branches"]:
                print(f"   unreachable: L{u['line']} {u['source']}")
                print(f"                {u['why']}")
        print(f"   all reachable branches covered: "
              f"{m['all_reachable_branches_covered']}")
        print()
    ct = out["core_totals"]
    print(f"core branch   : {ct['percent_branches_before']} % -> "
          f"{ct['percent_branches_after']} %  "
          f"({ct['covered_before']} -> {ct['covered_after']} of {ct['branches']})")
    print(f"core combined : {ct['percent_combined_before']} % -> "
          f"{ct['percent_combined_after']} %")
    print()
    print("defects found by writing the tests, none patched:")
    for d in out["defects_found_by_writing_the_tests"]:
        print(f"  {d['id']}  L{d['line']:<4} {d['kind']}")
    print()
    print(f"written -> {OUT}")


if __name__ == "__main__":
    main()

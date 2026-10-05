"""run_v1_v2_comparison.py -- the four defects closed, and what closing them changed.

Executes M4.2.

`10_new_system/` is untouched. Every number the manuscript reports is a
measurement of v1 and stays exactly as published; v2 is measured beside it. The
patched tree is `11_patched_pipeline_2026-09-22/`, built by copy with a
SHA-256 manifest of its source in `_PATCH_MANIFEST.json`.

What was closed:

  F1  the AUTO_REPAIR branch, unreachable from the engine -- removed, not
      repaired. Code that claims to handle a case it cannot receive reads as
      coverage.
  F2  the tautological guard inside it, which re-tested `result.state` against
      PASS after its own alias had been matched as AUTO_REPAIR.
  F3  the missing final `else`. A state outside the enumeration set neither
      flag and the workflow continued as though the gate had passed.
  F4  the source of that state: `StateMachine.current` is `AuditState | None`
      and reached a `ContractResult` field typed `AuditState`. mypy flags the
      line; the engine now refuses.
  F5  six G8 checks that passed by vacuous truth on empty input, while the
      contract's declared action on passing is to release the report.
  F6  `F821 undefined name wsl_root` -- a script that could not be imported.

The cheap comparisons run here. The expensive ones -- the test suite, branch
coverage, the seeded-defect benchmark -- are recorded with the command that
produced them, because re-running them inside this script would triple its cost
for no extra evidence.

Output -> 20_new_runs/patched_pipeline_v2_2026-09-22/
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kdss_paths import project_root, under  # noqa: E402

ROOT = project_root()
V1 = ROOT / "10_new_system"
V2 = ROOT / "11_patched_pipeline_2026-09-22"

# Measured outside this script; the command is recorded so each is reproducible.
RECORDED = {
    "test_suite": {
        "command": "python -m pytest -q   (run from the tree root)",
        "v1": {"passed": 715, "failed": 0},
        "v2": {"passed": 720, "failed": 0},
        "note": "two v1 tests were characterisation tests asserting the "
                "defective behaviour and necessarily fail under v2; they are "
                "replaced by their v2 counterparts. Five tests are new: four "
                "for the G8 fix and one for the stateless-result guard.",
    },
    "branch_coverage": {
        "command": "python -m coverage run --branch -m pytest -q "
                   "tests/test_orchestrator.py tests/test_contract_engine.py",
        "v1": {"orchestrator": "21/22", "contract_engine": "12/12"},
        "v2": {"orchestrator": "18/18", "contract_engine": "14/14"},
        "note": "the orchestrator's branch count falls from 22 to 18 because "
                "the dead branch and its tautology are gone, and every "
                "remaining branch is now covered. The single uncovered arc of "
                "v1 was the tautology's false arc, which no input could reach.",
    },
    "seeded_defect_benchmark": {
        "command": "python scripts/run_realism_benchmark.py   (from each tree)",
        "v1": {"conventional": "128/128, 55/55 FE", "tolerant": "128/128, 0/55 FE",
               "contracts": "124/128, 0/55 FE"},
        "v2": {"conventional": "128/128, 55/55 FE", "tolerant": "128/128, 0/55 FE",
               "contracts": "124/128, 0/55 FE"},
        "note": "IDENTICAL. Fixing four defects in the verification layer moves "
                "the benchmark not at all, because the harness calls the "
                "evaluators directly and re-implements the gating decision "
                "itself, so it never reaches the orchestrator (ADR-AIIM-198). "
                "The G8 fix does not move it either: the harness supplies every "
                "evidence key, so no check is asked to pass on absence. This is "
                "the benchmark's blind spot, measured rather than argued.",
    },
}


def probe_empty_gate(tree: Path) -> dict:
    """Ask each tree's G8 evaluator what it returns when handed nothing."""
    code = (
        "import io,json,sys,yaml;"
        f"sys.path.insert(0, r'{tree}');"
        "from evaluators import G8_evaluator as G;"
        f"c=yaml.safe_load(io.open(r'{tree / 'contracts' / 'G8_reporting.yaml'}',"
        "encoding='utf-8'));"
        "r=G.evaluate(contract=c, inputs={}, registry=None);"
        "print(json.dumps({'state': r['state'], 'n_passed': len(r['passed']),"
        "'n_failed': len(r['failed']), 'n_warned': len(r['warned']),"
        "'passed': r['passed'], 'failed': r['failed']}))"
    )
    out = subprocess.run([sys.executable, "-c", code], capture_output=True,
                         text=True, cwd=str(tree))
    try:
        return json.loads(out.stdout.strip().splitlines()[-1])
    except Exception:                                        # noqa: BLE001
        return {"error": (out.stderr or out.stdout)[-400:]}


def ruff_f821(tree: Path) -> int:
    out = subprocess.run(
        [sys.executable, "-m", "ruff", "check", "--no-cache", "--select",
         "F821", "--output-format", "json", str(tree)],
        capture_output=True, text=True)
    try:
        return len(json.loads(out.stdout or "[]"))
    except json.JSONDecodeError:
        return -1


def mypy_audit(tree: Path, cache: Path) -> int:
    out = subprocess.run(
        [sys.executable, "-m", "mypy", "--cache-dir", str(cache),
         "--ignore-missing-imports", "--no-error-summary",
         str(tree / "kdss_audit")], capture_output=True, text=True)
    return len([ln for ln in (out.stdout or "").splitlines() if ": error:" in ln])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(
        under("20_new_runs", "patched_pipeline_v2_2026-09-22")))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    if not V2.is_dir():
        print(f"ABORT: patched tree missing at {V2}")
        return 1

    cache = Path(tempfile.mkdtemp(prefix="v1v2_mypy_"))
    gate = {"v1": probe_empty_gate(V1), "v2": probe_empty_gate(V2)}
    f821 = {"v1": ruff_f821(V1), "v2": ruff_f821(V2)}
    mypy = {"v1": mypy_audit(V1, cache / "a"), "v2": mypy_audit(V2, cache / "b")}

    print("=" * 78)
    print("THE FOUR DEFECTS CLOSED, AND WHAT CLOSING THEM CHANGED   (M4.2)")
    print("=" * 78)
    print(f"  {'measurement':44}{'v1':>14}{'v2':>16}")
    print("  " + "-" * 74)
    print(f"  {'G8 handed no evidence: state':44}"
          f"{str(gate['v1'].get('state')):>14}{str(gate['v2'].get('state')):>16}")
    print(f"  {'   checks passing by vacuous truth':44}"
          f"{gate['v1'].get('n_passed', '?'):>14}{gate['v2'].get('n_passed', '?'):>16}")
    print(f"  {'   checks failing closed':44}"
          f"{gate['v1'].get('n_failed', '?'):>14}{gate['v2'].get('n_failed', '?'):>16}")
    print(f"  {'ruff F821 undefined names in the tree':44}"
          f"{f821['v1']:>14}{f821['v2']:>16}")
    print(f"  {'mypy errors in kdss_audit/':44}{mypy['v1']:>14}{mypy['v2']:>16}")
    for k, d in RECORDED.items():
        print(f"  {k.replace('_', ' '):44}"
              f"{str(list(d['v1'].values())[0]):>14}"
              f"{str(list(d['v2'].values())[0]):>16}")
    print()
    print("  The seeded-defect benchmark is unchanged, and that is a result:")
    print("  the harness never reaches the component the fixes are in.")

    summary = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "executes": "M4.2",
        "v1_tree": "10_new_system  (frozen; every published number is a "
                   "measurement of this tree and is unchanged)",
        "v2_tree": "11_patched_pipeline_2026-09-22  (copy; see "
                   "_PATCH_MANIFEST.json for the source SHA-256 of every file)",
        "measured_here": {
            "g8_empty_input_probe": gate,
            "ruff_F821_undefined_names": f821,
            "mypy_errors_in_kdss_audit": mypy,
        },
        "measured_separately": RECORDED,
        "what_this_does_not_establish": (
            "that the repair validates the v1 results retrospectively. The "
            "reported numbers are measurements of v1 and remain so. v2 "
            "establishes that the defects were real, that they are closeable, "
            "and what closing them costs -- for a next release with a re-run "
            "behind it."
        ),
        "not_fixed": [
            {"what": "40_deposit/pipeline/scripts/build_figures_1_2.py, which "
                     "raises NameError",
             "why": "the live tree is already correct; only the deposited copy "
                    "is stale, and the deposit must stay at the commit both "
                    "manuscripts cite (ADR-AIIM-190). A re-release fixes it; a "
                    "patch would invalidate the citation."},
        ],
    }
    with open(out / "v1_v2_comparison.json", "w", encoding="utf-8",
              newline="\n") as f:
        json.dump(summary, f, indent=1, ensure_ascii=False)
        f.write("\n")
    print(f"\nwritten -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

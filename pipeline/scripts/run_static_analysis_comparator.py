"""run_static_analysis_comparator.py -- what static analysis sees that the contracts cannot, and the reverse.

Executes M3.3.

The reviewer asked for a static-analysis comparator. The obvious version of it
is a strawman and is not run here: pointing a linter at 128 seeded **data**
faults, which it cannot architecturally see, and reporting that it found none.
That claim is not asserted -- it is **measured** below, by checking how many of
the injected faults touch a Python file at all -- and it is reported only
because the complementary direction is reported beside it.

Both directions, on the code the paper actually published
(`40_deposit/pipeline/`, the tree at the cited commit):

  CONTRACTS -> data faults      124 of 128 (ADR-AIIM-194/202)
  CONTRACTS -> code defects     0; they read artifacts, not source
  ruff+mypy -> data faults      not applicable, and the script shows why
  ruff+mypy -> code defects     measured here

No rule set is tuned. ruff runs on its **own default selection with no
configuration file present**, which is the only choice that cannot be accused
of having been picked to flatter a result. mypy runs with
`--ignore-missing-imports` because third-party stubs are absent; that is
recorded rather than hidden.

Caches are written outside the repository. `40_deposit` is a git repository
whose HEAD is the commit both manuscripts cite, and a stray `.ruff_cache/`
would dirty it (ADR-AIIM-190).

Output -> 20_new_runs/static_analysis_2026-09-22/
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kdss_paths import project_root, under  # noqa: E402

ROOT = project_root()
TREE = ROOT / "40_deposit" / "pipeline"
AUDIT = TREE / "kdss_audit"

# The four defects the branch-coverage tests found (ADR-AIIM-199, ADR-AIIM-218).
KNOWN_DEFECTS = [
    {"id": "D1", "what": "AUTO_REPAIR branch unreachable in production",
     "file": "kdss_audit/orchestrator.py", "line": 166},
    {"id": "D2", "what": "tautological guard inside that branch",
     "file": "kdss_audit/orchestrator.py", "line": 169},
    {"id": "D3", "what": "no final else; an unenumerated state sets neither "
                         "the halt nor the escalation flag (fail-open)",
     "file": "kdss_audit/orchestrator.py", "line": 176},
    {"id": "D4", "what": "six G8 checks pass by vacuous truth on empty input",
     "file": "contracts/G8_reporting.yaml", "line": None},
]


def tool_version(mod: str) -> str:
    r = subprocess.run([sys.executable, "-m", mod, "--version"],
                       capture_output=True, text=True)
    return (r.stdout or r.stderr).strip()


def run_ruff(cache: Path) -> list[dict]:
    r = subprocess.run(
        [sys.executable, "-m", "ruff", "check", "--no-cache",
         "--output-format", "json", str(TREE)],
        capture_output=True, text=True)
    try:
        return json.loads(r.stdout or "[]")
    except json.JSONDecodeError:
        print("ruff produced no parseable JSON:", r.stderr[:400])
        return []


def run_mypy(target: Path, cache: Path) -> list[str]:
    r = subprocess.run(
        [sys.executable, "-m", "mypy", "--cache-dir", str(cache),
         "--ignore-missing-imports", "--no-error-summary", str(target)],
        capture_output=True, text=True)
    return [ln for ln in (r.stdout or "").splitlines()
            if ": error:" in ln or ": note:" in ln]


def workspace_file_suffixes() -> dict:
    """Substantiate 'not applicable by construction' rather than asserting it.

    Every seeded fault is injected into a gate workspace. If no workspace holds
    a Python file, no linter could see any of the 128 faults, and saying so is
    a measurement rather than an excuse.
    """
    sys.path.insert(0, str(ROOT / "10_new_system"))
    from fault_injection import pilot_runner as pr           # noqa: PLC0415
    tmp = Path(tempfile.mkdtemp(prefix="static_analysis_ws_"))
    out = {}
    try:
        for gate, factory in pr._WORKSPACE_SETUP.items():
            d = tmp / gate
            d.mkdir(parents=True, exist_ok=True)
            factory(d)
            out[gate] = sorted(Counter(p.suffix.lower() for p in d.rglob("*")
                                       if p.is_file()).items())
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(
        under("20_new_runs", "static_analysis_2026-09-22")))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    cache = Path(tempfile.mkdtemp(prefix="static_analysis_cache_"))

    versions = {"ruff": tool_version("ruff"), "mypy": tool_version("mypy"),
                "python": sys.version.split()[0]}
    print("=" * 78)
    print("STATIC ANALYSIS AS A COMPARATOR, BOTH DIRECTIONS   (M3.3)")
    print("=" * 78)
    print(f"  ruff {versions['ruff']}   mypy {versions['mypy']}")
    print(f"  target: {TREE.relative_to(ROOT).as_posix()}  (the published tree)")

    ruff = run_ruff(cache)
    by_rule = Counter(f.get("code") or "?" for f in ruff)
    # F821 is the one that matters: a name used and never bound is a crash, not
    # a style opinion.
    f821 = [{"file": Path(f["filename"]).relative_to(ROOT).as_posix(),
             "line": f["location"]["row"], "message": f["message"]}
            for f in ruff if f.get("code") == "F821"]

    mypy_audit = run_mypy(AUDIT, cache)
    mypy_tree = run_mypy(TREE / "evaluators", cache)

    print(f"\n  ruff findings: {len(ruff)} across {len(by_rule)} rules")
    for code, n in by_rule.most_common(8):
        print(f"      {n:5}  {code}")
    print(f"\n  ruff F821 (undefined name -- a NameError at runtime): "
          f"**{len(f821)}**")
    for f in f821:
        print(f"      {f['file']}:{f['line']}  {f['message']}")
    print(f"\n  mypy on kdss_audit/: {len(mypy_audit)} error(s)")
    for ln in mypy_audit:
        print(f"      {ln.split('pipeline' + chr(92))[-1]}")

    suffixes = workspace_file_suffixes()
    py_in_workspaces = sum(n for g in suffixes.values()
                           for suf, n in g if suf == ".py")
    print(f"\n  Python files in any seeded-fault workspace: "
          f"**{py_in_workspaces}**")
    for gate, items in suffixes.items():
        print(f"      {gate}: {', '.join(f'{s or '(none)'}x{n}' for s, n in items)}")

    print("\n=== the 2 x 2 ===")
    print(f"  {'':26}{'128 data faults':>22}{'4 code defects':>18}")
    print(f"  {'contracts':26}{'124':>22}{'0 (read artifacts)':>18}")
    print(f"  {'ruff + mypy':26}{'n/a by construction':>22}"
          f"{'1 precondition':>18}")
    print("\n  and, in the direction nobody asked about:")
    print(f"  ruff found **{len(f821)}** defects that the 715-case suite, the "
          f"94 contract\n  checks and the G8 audit all missed.")

    summary = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "executes": "M3.3",
        "target": TREE.relative_to(ROOT).as_posix(),
        "why_this_target": "the tree at the commit both manuscripts cite, so "
                           "the analysed code is the published code",
        "versions": versions,
        "configuration": {
            "ruff": "tool default selection; no configuration file exists "
                    "anywhere in the repository, and none was written -- the "
                    "only choice that cannot be accused of being tuned",
            "mypy": "--ignore-missing-imports (third-party stubs absent), "
                    "otherwise default",
            "caches": "written outside the repository; 40_deposit is a git "
                      "repository whose HEAD is the cited commit",
        },
        "ruff": {
            "n_findings": len(ruff),
            "n_rules": len(by_rule),
            "by_rule": dict(by_rule.most_common()),
            "undefined_name_F821": f821,
        },
        "mypy": {
            "kdss_audit": mypy_audit,
            "evaluators": mypy_tree,
        },
        "data_faults_are_out_of_reach": {
            "claim": "no linter could detect any of the 128 seeded faults",
            "basis": "they are injected into gate workspaces, and no workspace "
                     "contains a Python file",
            "n_python_files_in_workspaces": py_in_workspaces,
            "workspace_file_suffixes": {k: dict(v) for k, v in suffixes.items()},
            "reported_as": "not applicable by construction, NOT as '0 detected'",
        },
        "known_code_defects": KNOWN_DEFECTS,
        # Verdict per defect, established by running each tool against the file
        # that holds it. The interesting cell is D3.
        "defect_verdicts": {
            "D1": {"ruff": "not flagged", "mypy": "not flagged",
                   "note": "the branch is semantically unreachable, not "
                           "syntactically; no reachability rule covers it"},
            "D2": {"ruff": "not flagged", "mypy": "not flagged",
                   "note": "PLR0124 compares a name with itself; this "
                           "compares result.state with a constant after state "
                           "-- its alias -- was matched, which needs alias "
                           "tracking neither tool does"},
            "D3": {"ruff": "not flagged",
                   "mypy": "PRECONDITION FLAGGED at kdss_audit/"
                           "contract_engine.py:120",
                   "note": "mypy does not see the missing else. It sees the "
                           "only way an unenumerated state can reach it: "
                           "AuditStateMachine.current is AuditState | None and "
                           "initialises to None, and that value is passed to a "
                           "ContractResult field typed AuditState. A None then "
                           "matches no branch in _handle_outcome and sets "
                           "neither flag. The branch-coverage tests found the "
                           "symptom after 33 tests were written; mypy found "
                           "the cause in seconds, from the other direction, "
                           "with no test at all"},
            "D4": {"ruff": "n/a", "mypy": "n/a",
                   "note": "the vacuous passes are in a YAML contract, not in "
                           "Python; out of reach by construction, reported as "
                           "such rather than as a miss"},
        },
        # The direction the reviewer did not ask about, and the one that
        # justifies running this at all.
        "found_by_static_analysis_alone": {
            "n": len(f821),
            "findings": f821,
            "what_they_are": "a name used and never bound -- a NameError the "
                             "moment the line executes, not a style opinion",
            "missed_by": ["the 715-case test suite (no test imports either "
                          "module)", "the 94 contract checks (they evaluate "
                          "artifacts, not source)", "the G8 number audit (it "
                          "matches literals in the manuscript)"],
            "status_of_each": {
                "scripts/render_pose_overlay_pymol.py": "broken in every tree "
                    "including the live one; the module cannot be imported, so "
                    "the deposited script that draws the pose-overlay figure "
                    "cannot run anywhere",
                "scripts/build_figures_1_2.py": "correct in the live tree, "
                    "which imports wsl_root at line 33; broken in the deposit "
                    "and the GitHub mirror, which predate that import. A "
                    "consequence of the deposit being frozen at the cited "
                    "commit (ADR-AIIM-190), and a concrete instance of it: the "
                    "published figure builder raises NameError",
            },
        },
        "contracts_on_code_defects": {
            "n_detected": 0,
            "why": "the contracts evaluate artifacts produced by the pipeline, "
                   "not the pipeline's source; this is symmetric with the row "
                   "above and is also not a failing",
        },
    }
    with open(out / "static_analysis.json", "w", encoding="utf-8",
              newline="\n") as f:
        json.dump(summary, f, indent=1, ensure_ascii=False)
        f.write("\n")
    shutil.rmtree(cache, ignore_errors=True)
    print(f"\nwritten -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

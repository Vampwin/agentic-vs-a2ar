"""run_naive_baseline.py -- the driver `sanity_check_baseline.py` never had.

`sanity_check_baseline.py` defines a naive-analyst comparator: existence and
parse checks, loose row-count plausibility, duplicate-ID hygiene, RDKit
parseability. Its scope was fixed before it was run against any fault family and
is documented in its own module docstring. It was run once, on 2026-08-10,
producing `20_new_runs/naive_baseline_2026-08-10/naive_baseline_outcomes.csv`
and a summary reporting SEER 0.477 -- and the code that drove it was not kept.
`evaluate()` has had no caller in this repository since.

That number now carries more weight than it did. After ADR-AIIM-131 corrected
the governed arm's denominator, Section 3.1 reports 0.031 against 0.477 as the
paper's like-for-like comparison, in preference to condition D's 1.000, which is
1.000 by construction. A comparator that load-bearing should be re-runnable, and
its value should be audited by G8 rather than transcribed. This script is that
driver (ADR-AIIM-131).

**The workspaces are the ones the contracts were scored on.** `_case_seed()`
mixes the condition label into the per-case seed, so running the baseline under
a condition label of its own would build different workspaces and compare two
evaluators on two different sets of files. It runs under condition **B** --
contracts on, monitor off -- so that every workspace is constructed and faulted
exactly as it was for the contract arm, and the only thing that differs between
the two results is which evaluator reads the workspace. The contracts' own
verdicts are discarded here; only the workspace is used.

Vacuity is recorded the same way `pilot_runner.py` records it, from the
injector's own `affected_items`, and excluded from the escape rate for the same
reason: an instance that changed no bytes carries no fault to escape.

Usage:
    python run_naive_baseline.py [--out DIR] [--keep-workspaces]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
import tempfile
from pathlib import Path
from kdss_paths import project_root, resolve, under, dev_dir

import pandas as pd

# --- pre-registered before the run, transcribed from ADR-AIIM-131 ------------
CONDITION_FOR_WORKSPACES = "B"   # contracts-on workspaces; see docstring
BASE_SEED = 42                   # PilotConfig default, unchanged
N_CLEAN_PER_FAMILY = 5           # PilotConfig default, unchanged
EXCLUDE_VACUOUS = True
EXCLUDE_ERROR = True
# A naive analyst has no escalation channel, so the baseline returns only
# PASS or HALT. A fault escapes when the baseline does not HALT.
DETECTED_STATES = {"HALT"}


def _root() -> Path:
    # Delegates to kdss_paths, which finds whichever layout this file is inside
    # -- the development tree or a deposit -- by walking up from its own
    # location. It previously tried two absolute paths and exited if neither
    # existed, which on the author's machine meant a copy of the deposit read
    # and wrote the development tree instead of itself.
    return project_root()


ROOT = _root()
RUNS = dev_dir("20_new_runs")
PREVIOUS = under("20_new_runs", "naive_baseline_2026-08-10")
sys.path.insert(0, str(ROOT / "10_new_system"))

from fault_injection import pilot_runner as pr           # noqa: E402
from fault_injection import sanity_check_baseline as nb  # noqa: E402
from fault_injection.manifest import FAULT_REGISTRY      # noqa: E402


def run_case(family_id: str, idx: int, is_clean: bool, ws_parent: Path) -> dict:
    injector = FAULT_REGISTRY[family_id]
    spec = injector.spec
    gate = spec.gate
    case_type = "clean" if is_clean else "fault"
    seed = pr._case_seed(BASE_SEED, family_id, idx, CONDITION_FOR_WORKSPACES,
                         case_type)
    instance_id = (f"{family_id}_clean_{idx:03d}" if is_clean
                   else injector.instance_id(idx))

    ws = ws_parent / f"{family_id}_{idx:03d}_{case_type}"
    ws.mkdir(parents=True, exist_ok=True)
    row = {"family_id": family_id, "instance_id": instance_id, "gate": gate,
           "case_type": case_type, "detected": False, "vacuous": False,
           "state": None, "findings": "", "error": None}
    try:
        if gate in pr._WORKSPACE_SETUP:
            pr._WORKSPACE_SETUP[gate](ws)
        if not is_clean:
            params = spec.param_grid[idx % len(spec.param_grid)]
            rec = injector.inject(ws, params, seed, idx)
            row["vacuous"] = not getattr(rec, "affected_items", None)
        out = nb.evaluate(gate, ws)
        row["state"] = out["state"]
        row["detected"] = out["state"] in DETECTED_STATES
        row["findings"] = "; ".join(out["findings"])
    except NotImplementedError as e:
        row["error"] = f"injector stub: {e}"
    except Exception as e:                       # noqa: BLE001
        row["error"] = f"{type(e).__name__}: {e}"
    return row


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(under("20_new_runs", "naive_baseline_2026-08-29")))
    ap.add_argument("--keep-workspaces", action="store_true")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    cfg = pr.PilotConfig()
    ws_root = Path(tempfile.mkdtemp(prefix="naive_baseline_"))
    rows: list[dict] = []
    try:
        for family_id in cfg.pilot_families:
            n = len(FAULT_REGISTRY[family_id].spec.param_grid)
            for idx in range(n):
                rows.append(run_case(family_id, idx, False, ws_root))
            for idx in range(N_CLEAN_PER_FAMILY):
                rows.append(run_case(family_id, idx, True, ws_root))
    finally:
        if not args.keep_workspaces:
            shutil.rmtree(ws_root, ignore_errors=True)

    df = pd.DataFrame(rows)
    faults = df[df.case_type == "fault"]
    counted = faults[~faults.vacuous & faults.error.isna()] if EXCLUDE_VACUOUS \
        else faults[faults.error.isna()]
    escaped = int((~counted.detected).sum())
    clean = df[df.case_type == "clean"]

    vac = sorted(faults.loc[faults.vacuous, "instance_id"])
    per_family = {
        f: {"n_injected": int((faults.family_id == f).sum()),
            "n_counted": int((counted.family_id == f).sum()),
            "n_detected": int(counted.loc[counted.family_id == f, "detected"].sum())}
        for f in cfg.pilot_families}

    summary = {
        "what_this_is": "The naive-analyst comparator of sanity_check_baseline.py, "
                        "re-run by a driver rather than transcribed. Section 3.1 "
                        "reports its SEER against the governed arm's.",
        "authority": "ADR-AIIM-131",
        "workspaces": f"Built under condition {CONDITION_FOR_WORKSPACES} seeds, so "
                      f"each workspace is byte-identical to the one the contract "
                      f"arm was scored on. Only the evaluator differs.",
        "base_seed": BASE_SEED,
        "seer_point": round(escaped / len(counted), 6),
        "n_families": len(cfg.pilot_families),
        "n_fault_cases_injected": len(faults),
        "n_fault_cases_counted": len(counted),
        "n_detected": int(counted.detected.sum()),
        "n_escaped": escaped,
        "n_clean_cases": len(clean),
        "false_escalation": round(float(clean.detected.sum()) / len(clean), 6)
                            if len(clean) else 0.0,
        "vacuous_excluded": {"n": len(vac), "instances": vac,
                             "all_in_family": sorted({i.split("_")[0] for i in vac})},
        "errors": int(faults.error.notna().sum()),
        "per_family": per_family,
    }

    prev_path = PREVIOUS / "naive_baseline_summary.json"
    if prev_path.is_file():
        prev = json.loads(prev_path.read_text(encoding="utf-8"))
        # Compare the unrounded rate against the stored one at the precision the
        # stored value carries. The first version compared a value already
        # rounded to six places against an unrounded 0.4765625 at a tolerance of
        # 1e-9 and reported a reproduction failure on a difference of 5e-7 that
        # the comparison itself had introduced.
        exact = escaped / len(counted)
        per_family_match = all(
            prev["per_family_n"][f] == per_family[f]["n_counted"]
            and prev["per_family_detected"][f] == per_family[f]["n_detected"]
            for f in prev["per_family_n"])
        summary["previous_run_2026-08-10"] = {
            "seer_point": prev["seer_point"],
            "n_fault_cases": prev["n_fault_cases"],
            "seer_matches": abs(prev["seer_point"] - exact) < 1e-9,
            "n_cases_matches": prev["n_fault_cases"] == len(counted),
            "per_family_matches": per_family_match,
            "reproduced": (abs(prev["seer_point"] - exact) < 1e-9
                           and prev["n_fault_cases"] == len(counted)
                           and per_family_match),
            "note": "That run had no driver script. Agreement here is what makes "
                    "its number reproducible rather than merely recorded; "
                    "disagreement would itself be the finding.",
        }

    df.to_csv(out / "naive_baseline_outcomes.csv", index=False, encoding="utf-8")
    (out / "naive_baseline_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")

    print(f"families            : {summary['n_families']}")
    print(f"fault cases injected: {summary['n_fault_cases_injected']}")
    print(f"  vacuous excluded  : {summary['vacuous_excluded']['n']}  "
          f"{', '.join(vac)}")
    print(f"  errors            : {summary['errors']}")
    print(f"fault cases counted : {summary['n_fault_cases_counted']}")
    print(f"  detected          : {summary['n_detected']}")
    print(f"  escaped           : {summary['n_escaped']}")
    print(f"SEER                : {summary['seer_point']:.6f}")
    print(f"false escalation    : {summary['false_escalation']:.6f} "
          f"on {summary['n_clean_cases']} clean cases")
    p = summary.get("previous_run_2026-08-10")
    if p:
        print(f"\n2026-08-10 run      : SEER {p['seer_point']:.6f} over "
              f"{p['n_fault_cases']} cases")
        print(f"reproduced          : {p['reproduced']}")
    print(f"\nwritten -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

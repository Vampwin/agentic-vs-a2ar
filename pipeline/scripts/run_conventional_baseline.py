"""run_conventional_baseline.py -- RA2 / major comment 2.

Runs the conventional-test comparator of `conventional_test_baseline.py` over
the same seeded faults the contract arm was scored on, and reports SEER the
same way `run_naive_baseline.py` does.

**The workspaces are the ones the contracts were scored on.** `_case_seed()`
mixes the condition label into the per-case seed, so running under a label of
its own would build different workspaces and compare two evaluators on two
different sets of files. It runs under condition **B**, exactly as the naive
baseline does, so the only thing differing between the three arms is which
evaluator reads the workspace.

Scope and the falsifiable prediction are pre-registered in ADR-AIIM-193, which
was written before this ran.

Usage:
    python run_conventional_baseline.py [--out DIR] [--keep-workspaces]
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
import tempfile
from collections import Counter
from pathlib import Path

import pandas as pd

from kdss_paths import project_root, under, dev_dir

# --- pre-registered before the run, transcribed from ADR-AIIM-193 -----------
CONDITION_FOR_WORKSPACES = "B"
BASE_SEED = 42
N_CLEAN_PER_FAMILY = 5
EXCLUDE_VACUOUS = True
EXCLUDE_ERROR = True
DETECTED_STATES = {"HALT"}

ROOT = project_root()
sys.path.insert(0, str(ROOT / "10_new_system"))

from fault_injection import pilot_runner as pr                    # noqa: E402
from fault_injection import conventional_test_baseline as ct      # noqa: E402
from fault_injection.manifest import FAULT_REGISTRY               # noqa: E402


def run_case(family_id: str, idx: int, is_clean: bool, ws_parent: Path,
             golden_root: Path) -> dict:
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
           "state": None, "findings": "", "practices": "", "error": None}
    try:
        if gate in pr._WORKSPACE_SETUP:
            pr._WORKSPACE_SETUP[gate](ws)
        if not is_clean:
            params = spec.param_grid[idx % len(spec.param_grid)]
            rec = injector.inject(ws, params, seed, idx)
            row["vacuous"] = not getattr(rec, "affected_items", None)
        out = ct.evaluate(gate, ws, golden_root / gate)
        row["state"] = out["state"]
        row["detected"] = out["state"] in DETECTED_STATES
        row["findings"] = "; ".join(out["findings"])
        row["practices"] = ",".join(sorted({f.split(":")[0]
                                            for f in out["findings"]}))
    except NotImplementedError as e:
        row["error"] = f"injector stub: {e}"
    except Exception as e:                       # noqa: BLE001
        row["error"] = f"{type(e).__name__}: {e}"
    return row


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(
        under("20_new_runs", "conventional_baseline_2026-09-20")))
    ap.add_argument("--keep-workspaces", action="store_true")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    cfg = pr.PilotConfig()
    ws_root = Path(tempfile.mkdtemp(prefix="conventional_baseline_"))
    golden_root = ws_root / "_reference"

    # practice 4: deterministic-rerun equality. Build the reference twice and
    # confirm it is stable, so a drifting reference cannot be mistaken for a
    # clean result. Recorded, not assumed.
    rerun: dict[str, bool] = {}
    for gate in pr._WORKSPACE_SETUP:
        a = ct.build_reference(gate, golden_root / gate, pr._WORKSPACE_SETUP)
        b = ct.build_reference(gate, ws_root / "_rerun" / gate, pr._WORKSPACE_SETUP)
        rerun[gate] = all(
            (a / f).is_file() and (b / f).is_file()
            and ct._sha(a / f) == ct._sha(b / f)
            for f in ct.GATE_FILES[gate])

    rows: list[dict] = []
    try:
        for family_id in cfg.pilot_families:
            n = len(FAULT_REGISTRY[family_id].spec.param_grid)
            for idx in range(n):
                rows.append(run_case(family_id, idx, False, ws_root, golden_root))
            for idx in range(N_CLEAN_PER_FAMILY):
                rows.append(run_case(family_id, idx, True, ws_root, golden_root))
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

    practice_counts = Counter()
    for p in counted.loc[counted.detected, "practices"]:
        for tag in (p.split(",") if p else []):
            practice_counts[tag] += 1

    summary = {
        "what_this_is": "A conventional software-testing comparator -- snapshot, "
                        "schema and referential-integrity checks -- over the same "
                        "seeded faults the contract arm was scored on. Answers "
                        "reviewer RA2 and the first half of major comment 2.",
        "authority": "ADR-AIIM-193 (pre-registration, written before this ran)",
        "workspaces": f"Built under condition {CONDITION_FOR_WORKSPACES} seeds, so "
                      f"each workspace is byte-identical to the one the contract "
                      f"arm and the naive comparator were scored on.",
        "base_seed": BASE_SEED,
        "seer_point": round(escaped / len(counted), 6) if len(counted) else None,
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
        "per_family": per_family,
        "detections_by_practice": dict(practice_counts),
        "deterministic_rerun_equality": rerun,
        "prediction_registered": "Golden-file testing should catch value-changing "
                                 "and stale-artifact faults and MISS statistical, "
                                 "graded ones. ADR-AIIM-193.",
        "confound_registered": "The harness builds workspaces deterministically, so "
                               "byte-exact snapshot comparison is far stronger here "
                               "than on a real pipeline whose outputs vary between "
                               "runs. Stated in ADR-AIIM-193 before the run.",
    }

    (out / "conventional_baseline_outcomes.csv").write_text(
        df.to_csv(index=False), encoding="utf-8")
    (out / "conventional_baseline_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")

    print(f"conventional baseline: SEER {summary['seer_point']}  "
          f"({summary['n_detected']}/{summary['n_fault_cases_counted']} detected, "
          f"{escaped} escaped)")
    print(f"  false escalation on {len(clean)} clean cases: "
          f"{summary['false_escalation']}")
    print(f"  detections by practice: {dict(practice_counts)}")
    print(f"  deterministic rerun equality: {rerun}")
    print(f"  vacuous excluded: {len(vac)}")
    print(f"written -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

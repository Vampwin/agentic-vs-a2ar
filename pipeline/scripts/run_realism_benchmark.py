"""run_realism_benchmark.py -- both detectors, same workspaces, realistic noise.

Executes the experiment pre-registered in ADR-AIIM-201. Every prediction and
falsification condition was fixed in the decision log before this ran.

Design note: the comparison is **paired**. Each case builds one workspace,
perturbs it with the neutral run-to-run variation of
`fault_injection.realistic_nondeterminism`, and then runs *both* detectors on
that same directory:

  * the contract evaluator, anchored on `_clean_hashes` captured from this
    workspace before injection;
  * the conventional golden-file comparator, anchored on a reference built by
    re-running the factory.

Running them on separately generated workspaces would confound the detector
with the draw. They see the same bytes.

Scoring follows `pilot_runner` exactly, so the numbers are comparable to
ADR-AIIM-194's:

  contracts   detected  = evaluator state in {HALT, ESCALATE}
  conventional detected = state == HALT   (it has no escalation channel)

and on clean cases the same rule is a false escalation.

Output -> 20_new_runs/realism_benchmark_2026-09-21/
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kdss_paths import project_root, under  # noqa: E402

BASE_SEED = 42
CONDITION_FOR_WORKSPACES = "B"
NOISE_SEED = 20260921

ROOT = project_root()
sys.path.insert(0, str(ROOT / "10_new_system"))

from fault_injection import pilot_runner as pr                    # noqa: E402
from fault_injection import conventional_test_baseline as ct      # noqa: E402
from fault_injection import realistic_nondeterminism as rn        # noqa: E402
from fault_injection.manifest import FAULT_REGISTRY               # noqa: E402


def run_case(family_id: str, idx: int, is_clean: bool, ws_parent: Path,
             golden_root: Path, setup_map: dict) -> dict:
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
    row = {
        "family_id": family_id, "instance_id": instance_id, "gate": gate,
        "case_type": case_type, "vacuous": False,
        "conv_state": None, "conv_detected": False, "conv_practices": "",
        "tol_state": None, "tol_detected": False,
        "contract_state": None, "contract_available": False,
        "contract_detected": False,
        "error": None,
    }
    try:
        if gate in setup_map:
            setup_map[gate](ws)

        # Contract anchor: captured from THIS workspace before injection.
        clean_hashes = pr._clean_hashes(ws, gate)

        if not is_clean:
            params = spec.param_grid[idx % len(spec.param_grid)]
            rec = injector.inject(ws, params, seed, idx)
            row["vacuous"] = not getattr(rec, "affected_items", None)

        # --- detector 1: conventional golden-file comparison ---------------
        out = ct.evaluate(gate, ws, golden_root / gate)
        row["conv_state"] = out["state"]
        row["conv_detected"] = out["state"] == "HALT"
        row["conv_practices"] = ",".join(sorted({f.split(":")[0]
                                                 for f in out["findings"]}))

        # --- detector 2: the strong-man, normalised before comparing -------
        tol = ct.evaluate_tolerant(gate, ws, golden_root / gate)
        row["tol_state"] = tol["state"]
        row["tol_detected"] = tol["state"] == "HALT"

        # --- detector 3: the contract evaluator ----------------------------
        ev, available = pr._dispatch_evaluator(gate, ws, clean_hashes)
        row["contract_available"] = bool(available)
        state = ev.get("state") if ev else None
        row["contract_state"] = state
        row["contract_detected"] = bool(available) and state in ("HALT", "ESCALATE")
    except NotImplementedError as e:
        row["error"] = f"injector stub: {e}"
    except Exception as e:                       # noqa: BLE001
        row["error"] = f"{type(e).__name__}: {e}"
    return row


def summarise(df: pd.DataFrame, arm: str) -> dict:
    fault = df[(df.case_type == "fault") & (~df.vacuous) & (df.error.isna())]
    clean = df[(df.case_type == "clean") & (df.error.isna())]
    det = int(fault[f"{arm}_detected"].sum())
    fe = int(clean[f"{arm}_detected"].sum())
    sens = det / len(fault) if len(fault) else None
    spec = 1.0 - fe / len(clean) if len(clean) else None
    return {
        "n_fault_cases_counted": int(len(fault)),
        "n_detected": det,
        "detection_rate": round(sens, 4) if sens is not None else None,
        "seer": round(1.0 - sens, 4) if sens is not None else None,
        "n_clean_cases": int(len(clean)),
        "n_false_escalations": fe,
        "false_escalation_rate": round(fe / len(clean), 4) if len(clean) else None,
        "specificity": round(spec, 4) if spec is not None else None,
        # A detector that fires on everything scores perfect detection and is
        # worthless, so ADR-AIIM-202 fixes this as the single summary figure.
        "balanced_accuracy": (round((sens + spec) / 2, 4)
                              if sens is not None and spec is not None else None),
        "youden_j": (round(sens + spec - 1.0, 4)
                     if sens is not None and spec is not None else None),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(
        under("20_new_runs", "realism_benchmark_2026-09-21")))
    ap.add_argument("--keep-workspaces", action="store_true")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    cfg = pr.PilotConfig()
    ws_root = Path(tempfile.mkdtemp(prefix="realism_benchmark_"))
    golden_root = ws_root / "_reference"

    rng = np.random.default_rng(NOISE_SEED)
    setup_map = rn.make_setup_map(pr._WORKSPACE_SETUP, rng)

    # Validity gate, from ADR-AIIM-201: two builds must differ in bytes and
    # agree in science. If either fails the experiment is invalid and no
    # result is written.
    print("=== validity gate (ADR-AIIM-201) ===")
    import hashlib
    validity = {}
    ok_all = True
    for gate in pr._WORKSPACE_SETUP:
        a = ws_root / "_v" / f"{gate}_a"
        b = ws_root / "_v" / f"{gate}_b"
        a.mkdir(parents=True, exist_ok=True)
        b.mkdir(parents=True, exist_ok=True)
        setup_map[gate](a)
        setup_map[gate](b)

        def h(d: Path) -> str:
            return hashlib.sha256(b"".join(
                sorted(p.read_bytes() for p in d.rglob("*") if p.is_file()))
            ).hexdigest()

        bytes_differ = h(a) != h(b)
        same_science, why = rn.semantic_equal(a, b)
        validity[gate] = {"bytes_differ": bytes_differ,
                          "semantically_identical": same_science, "note": why}
        flag = "ok " if (bytes_differ and same_science) else "FAIL"
        print(f"  [{flag}] {gate}: bytes differ={bytes_differ} "
              f"science identical={same_science}")
        if not (bytes_differ and same_science):
            ok_all = False
    if not ok_all:
        print("\nABORT: perturbation is not neutral, or not reaching the files.")
        print("ADR-AIIM-201 calls this invalid rather than negative.")
        return 2

    # The reference the conventional comparator is anchored on, built by
    # re-running the factory exactly as build_reference does in production use.
    for gate in pr._WORKSPACE_SETUP:
        ct.build_reference(gate, golden_root / gate, setup_map)

    print("\n=== running both detectors over the same workspaces ===")
    rows = []
    for family_id in cfg.pilot_families:
        spec = FAULT_REGISTRY[family_id].spec
        n = len(spec.param_grid)
        for idx in range(n):
            rows.append(run_case(family_id, idx, False, ws_root, golden_root,
                                 setup_map))
        for idx in range(cfg.n_clean_per_family):
            rows.append(run_case(family_id, idx, True, ws_root, golden_root,
                                 setup_map))
        done = [r for r in rows if r["family_id"] == family_id]
        print(f"  {family_id}: {len(done)} cases")

    df = pd.DataFrame(rows)
    df.to_csv(out / "realism_benchmark_outcomes.csv", index=False)

    conv = summarise(df, "conv")
    tol = summarise(df, "tol")
    contract = summarise(df, "contract")

    summary = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "pre_registration": "ADR-AIIM-201",
        "design": "paired; both detectors run on the same perturbed workspace",
        "noise_seed": NOISE_SEED,
        "base_seed": BASE_SEED,
        "workspace_condition": CONDITION_FOR_WORKSPACES,
        "validity_gate": validity,
        "arms": {"conventional_golden_file": conv,
                 "tolerant_golden_file": tol,
                 "contracts": contract},
        "baseline_for_comparison": {
            "source": "ADR-AIIM-194, deterministic workspaces",
            "conventional": {"n_detected": 128, "n_fault_cases_counted": 128,
                             "false_escalation_rate": 0.0},
            "contracts": {"n_detected": 124, "n_fault_cases_counted": 128,
                          "false_escalation_rate": 0.0},
        },
        "n_errors": int(df.error.notna().sum()),
    }
    with open(out / "realism_benchmark_summary.json", "w", encoding="utf-8",
              newline="\n") as f:
        json.dump(summary, f, indent=1, ensure_ascii=False)
        f.write("\n")

    print("\n" + "=" * 72)
    print("RESULT -- reported as (detection, false escalation), per ADR-AIIM-201")
    print("=" * 72)
    print(f"{'arm':30} {'detection':>14} {'false escalation':>18} {'bal.acc':>9}")
    for name, s in (("conventional golden-file", conv),
                    ("tolerant golden-file", tol),
                    ("contracts", contract)):
        print(f"{name:30} {s['n_detected']:>6}/{s['n_fault_cases_counted']:<7} "
              f"{s['n_false_escalations']:>9}/{s['n_clean_cases']:<8} "
              f"{s['balanced_accuracy']:>9}")
    print()
    print("under deterministic workspaces (ADR-AIIM-194):")
    print(f"{'conventional golden-file':34} {128:>7}/{128:<7} {0:>10}/{55:<8}")
    print(f"{'contracts':34} {124:>7}/{128:<7} {0:>10}/{55:<8}")
    print()
    print(f"errors: {summary['n_errors']}")
    print(f"written -> {out}")

    if not args.keep_workspaces:
        shutil.rmtree(ws_root, ignore_errors=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

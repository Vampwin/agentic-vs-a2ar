"""run_contaminated_reference_benchmark.py -- what each detector does when its reference is defective.

Executes M3.1, pre-registered in ADR-AIIM-223 before any case ran.

ADR-AIIM-202 measured both detector families against a **trustworthy**
reference and found the tolerant golden-file comparator ahead of the contracts.
That result stands. It assumed something an agent-written pipeline does not get
to assume: that a known-good earlier run exists to compare against. On a first
run there is none -- the first output *is* the baseline -- and the test-oracle
literature names this as the standing limit of comparison-based oracles.

**Symmetric contamination is the whole design.** The contract arm is anchored on
`pilot_runner._clean_hashes`, which is a reference too. Contaminating the golden
directory while leaving that anchor clean would invent an advantage for the
contracts and the result would be worthless. Both references are therefore built
from a workspace carrying the same fault, produced by a separate execution with
a different seed -- the defect is in the code, so every run carries it,
including the one whose output became the baseline.

Paired at the case level. Each case is evaluated twice against **identical
workspace bytes**:

    CLEAN         reference built from a clean workspace   (reproduces ADR-AIIM-202)
    CONTAMINATED  reference built from a faulted workspace

The pre-registration withdrew the fraction-p curve of the approved plan: with
independent per-case contamination the rate at any mixture is the linear
combination of these two endpoints, so the intermediate points interpolate a
straight line at four times the cost.

Clean cases matter as much as fault cases here. A contaminated reference makes a
*correct* run look wrong, so contamination should cost specificity as well as
sensitivity, and both are reported.

Output -> 20_new_runs/contaminated_reference_2026-09-22/
"""
from __future__ import annotations

import argparse
import hashlib
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

BASE_SEED = 42            # carried over from ADR-AIIM-201, unchanged
REF_BASE_SEED = 43        # the reference run is a *different* execution
CONDITION_FOR_WORKSPACES = "B"
NOISE_SEED = 20260921     # carried over unchanged
REF_NOISE_SEED = 20260922  # a separate stream; see run_case

# The honesty gate of ADR-AIIM-223: arm CLEAN must return exactly this.
PUBLISHED = {
    "conv": {"n_detected": 128, "n_false_escalations": 55},
    "tol": {"n_detected": 128, "n_false_escalations": 0},
    "contract": {"n_detected": 124, "n_false_escalations": 0},
}
PUBLISHED_N_FAULT = 128
PUBLISHED_N_CLEAN = 55

ROOT = project_root()
sys.path.insert(0, str(ROOT / "10_new_system"))

from fault_injection import pilot_runner as pr                    # noqa: E402
from fault_injection import conventional_test_baseline as ct      # noqa: E402
from fault_injection import realistic_nondeterminism as rn        # noqa: E402
from fault_injection.manifest import FAULT_REGISTRY               # noqa: E402

ARMS = ("conv", "tol", "contract")
REFS = ("clean", "contam")


def run_case(family_id: str, idx: int, is_clean: bool, ws_parent: Path,
             clean_golden_root: Path, setup_map: dict,
             ref_setup_map: dict) -> dict:
    """One case, three detectors, two references.

    `setup_map` and `ref_setup_map` draw their neutral noise from **separate**
    generators, and that separation is load-bearing. `realistic_nondeterminism`
    closes each factory over one shared `Generator`, so every workspace build
    advances a single stream. Building the contaminated reference from the same
    stream shifts the noise every later test workspace receives, the run stops
    being byte-comparable with ADR-AIIM-202, and the honesty gate fails by one
    case. It did, on the first execution of this script, which is what the gate
    is for.
    """
    injector = FAULT_REGISTRY[family_id]
    spec = injector.spec
    gate = spec.gate
    case_type = "clean" if is_clean else "fault"
    params = spec.param_grid[idx % len(spec.param_grid)]
    instance_id = (f"{family_id}_clean_{idx:03d}" if is_clean
                   else injector.instance_id(idx))

    row = {
        "family_id": family_id, "instance_id": instance_id, "gate": gate,
        "case_type": case_type, "vacuous": False, "ref_vacuous": False,
        "error": None,
    }
    for arm in ARMS:
        for ref in REFS:
            row[f"{arm}_{ref}_state"] = None
            row[f"{arm}_{ref}_detected"] = False
    row["contract_available"] = False

    try:
        # ---- the workspace under test ---------------------------------
        ws = ws_parent / f"{family_id}_{idx:03d}_{case_type}"
        ws.mkdir(parents=True, exist_ok=True)
        if gate in setup_map:
            setup_map[gate](ws)

        # the contract anchor as the existing harness captures it: from this
        # workspace, before injection
        hashes_clean = pr._clean_hashes(ws, gate)

        if not is_clean:
            seed = pr._case_seed(BASE_SEED, family_id, idx,
                                 CONDITION_FOR_WORKSPACES, case_type)
            rec = injector.inject(ws, params, seed, idx)
            row["vacuous"] = not getattr(rec, "affected_items", None)

        # ---- the contaminated reference -------------------------------
        # A separate execution of the same defective code: same family, same
        # parameters, different seed. Injected even for clean cases, because
        # the scenario is a defective *baseline* against which a correct run is
        # later compared.
        ref_ws = ws_parent / f"{family_id}_{idx:03d}_{case_type}_REF"
        ref_ws.mkdir(parents=True, exist_ok=True)
        ct.build_reference(gate, ref_ws, ref_setup_map)
        ref_seed = pr._case_seed(REF_BASE_SEED, family_id, idx,
                                 CONDITION_FOR_WORKSPACES, case_type)
        ref_rec = injector.inject(ref_ws, params, ref_seed, idx)
        row["ref_vacuous"] = not getattr(ref_rec, "affected_items", None)

        # the contract anchor as it would have been registered from that
        # defective run: captured AFTER the fault is in place
        hashes_contam = pr._clean_hashes(ref_ws, gate)

        goldens = {"clean": clean_golden_root / gate, "contam": ref_ws}
        anchors = {"clean": hashes_clean, "contam": hashes_contam}

        for ref in REFS:
            out = ct.evaluate(gate, ws, goldens[ref])
            row[f"conv_{ref}_state"] = out["state"]
            row[f"conv_{ref}_detected"] = out["state"] == "HALT"

            tol = ct.evaluate_tolerant(gate, ws, goldens[ref])
            row[f"tol_{ref}_state"] = tol["state"]
            row[f"tol_{ref}_detected"] = tol["state"] == "HALT"

            ev, available = pr._dispatch_evaluator(gate, ws, anchors[ref])
            row["contract_available"] = bool(available)
            state = ev.get("state") if ev else None
            row[f"contract_{ref}_state"] = state
            row[f"contract_{ref}_detected"] = (bool(available)
                                               and state in ("HALT", "ESCALATE"))
    except NotImplementedError as e:
        row["error"] = f"injector stub: {e}"
    except Exception as e:                                   # noqa: BLE001
        row["error"] = f"{type(e).__name__}: {e}"
    return row


def summarise(df: pd.DataFrame, arm: str, ref: str) -> dict:
    col = f"{arm}_{ref}_detected"
    fault = df[(df.case_type == "fault") & (~df.vacuous) & (df.error.isna())]
    clean = df[(df.case_type == "clean") & (df.error.isna())]
    det = int(fault[col].sum())
    fe = int(clean[col].sum())
    sens = det / len(fault) if len(fault) else None
    spec = 1.0 - fe / len(clean) if len(clean) else None
    return {
        "n_fault_cases_counted": int(len(fault)),
        "n_detected": det,
        "detection_rate": round(sens, 4) if sens is not None else None,
        "seer": round(1.0 - sens, 4) if sens is not None else None,
        "n_clean_cases": int(len(clean)),
        "n_false_escalations": fe,
        "false_escalation_rate": (round(fe / len(clean), 4)
                                  if len(clean) else None),
        "specificity": round(spec, 4) if spec is not None else None,
        "balanced_accuracy": (round((sens + spec) / 2, 4)
                              if sens is not None and spec is not None else None),
    }


def per_gate(df: pd.DataFrame, arm: str) -> dict:
    """The decomposition the prediction turns on: which gates carry a hash anchor."""
    out = {}
    fault = df[(df.case_type == "fault") & (~df.vacuous) & (df.error.isna())]
    for gate in sorted(fault.gate.unique()):
        sub = fault[fault.gate == gate]
        c = int(sub[f"{arm}_clean_detected"].sum())
        k = int(sub[f"{arm}_contam_detected"].sum())
        out[gate] = {
            "n": int(len(sub)), "detected_clean_ref": c,
            "detected_contaminated_ref": k, "lost": c - k,
            "has_hash_anchor": gate in ("G0", "G8"),
        }
    return out


def per_gate_fe(df: pd.DataFrame, arm: str) -> dict:
    """False escalations by gate -- where a contaminated anchor makes a correct
    run look wrong. For the contract arm this is the whole mechanism."""
    out = {}
    clean = df[(df.case_type == "clean") & (df.error.isna())]
    for gate in sorted(clean.gate.unique()):
        sub = clean[clean.gate == gate]
        out[gate] = {
            "n_clean": int(len(sub)),
            "fe_clean_ref": int(sub[f"{arm}_clean_detected"].sum()),
            "fe_contaminated_ref": int(sub[f"{arm}_contam_detected"].sum()),
            "has_hash_anchor": gate in ("G0", "G8"),
        }
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(
        under("20_new_runs", "contaminated_reference_2026-09-22")))
    ap.add_argument("--keep-workspaces", action="store_true")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    cfg = pr.PilotConfig()
    ws_root = Path(tempfile.mkdtemp(prefix="contaminated_reference_"))
    clean_golden_root = ws_root / "_reference_clean"

    rng = np.random.default_rng(NOISE_SEED)
    setup_map = rn.make_setup_map(pr._WORKSPACE_SETUP, rng)
    # Separate generator for reference builds, so the draw sequence the test
    # workspaces see is identical to ADR-AIIM-202's.
    ref_rng = np.random.default_rng(REF_NOISE_SEED)
    ref_setup_map = rn.make_setup_map(pr._WORKSPACE_SETUP, ref_rng)

    # The validity gate of ADR-AIIM-201, unchanged: the perturbation must move
    # bytes and preserve the science, or the harness is measuring the wrong
    # thing and nothing is written.
    print("=== validity gate (ADR-AIIM-201) ===")
    validity, ok_all = {}, True
    for gate in pr._WORKSPACE_SETUP:
        a, b = ws_root / "_v" / f"{gate}_a", ws_root / "_v" / f"{gate}_b"
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
        ok_all = ok_all and bytes_differ and same_science
    if not ok_all:
        print("\nABORT: perturbation is not neutral. Invalid, not negative.")
        return 2

    for gate in pr._WORKSPACE_SETUP:
        ct.build_reference(gate, clean_golden_root / gate, setup_map)

    print("\n=== three detectors x two references, identical bytes ===")
    rows = []
    for family_id in cfg.pilot_families:
        spec = FAULT_REGISTRY[family_id].spec
        for idx in range(len(spec.param_grid)):
            rows.append(run_case(family_id, idx, False, ws_root,
                                 clean_golden_root, setup_map,
                                 ref_setup_map))
        for idx in range(cfg.n_clean_per_family):
            rows.append(run_case(family_id, idx, True, ws_root,
                                 clean_golden_root, setup_map,
                                 ref_setup_map))
        print(f"  {family_id}: {len([r for r in rows if r['family_id'] == family_id])} cases")

    df = pd.DataFrame(rows)
    df.to_csv(out / "contaminated_reference_outcomes.csv", index=False)

    arms = {a: {r: summarise(df, a, r) for r in REFS} for a in ARMS}

    # --- the honesty gate --------------------------------------------------
    gate_ok, gate_detail = True, {}
    for a, want in PUBLISHED.items():
        got = arms[a]["clean"]
        ok = (got["n_detected"] == want["n_detected"]
              and got["n_false_escalations"] == want["n_false_escalations"])
        gate_detail[a] = {"expected": want,
                          "got": {"n_detected": got["n_detected"],
                                  "n_false_escalations":
                                      got["n_false_escalations"]},
                          "reproduced": ok}
        gate_ok = gate_ok and ok
    n_err = int(df.error.notna().sum())
    counts_ok = (arms["contract"]["clean"]["n_fault_cases_counted"]
                 == PUBLISHED_N_FAULT
                 and arms["contract"]["clean"]["n_clean_cases"]
                 == PUBLISHED_N_CLEAN)

    print("\n=== honesty gate: arm CLEAN must reproduce ADR-AIIM-202 ===")
    for a, d in gate_detail.items():
        print(f"  [{'ok ' if d['reproduced'] else 'FAIL'}] {a:9} "
              f"detected {d['got']['n_detected']} (want "
              f"{d['expected']['n_detected']}), false escalations "
              f"{d['got']['n_false_escalations']} (want "
              f"{d['expected']['n_false_escalations']})")
    print(f"  [{'ok ' if counts_ok else 'FAIL'}] case counts "
          f"{arms['contract']['clean']['n_fault_cases_counted']}/"
          f"{arms['contract']['clean']['n_clean_cases']} "
          f"(want {PUBLISHED_N_FAULT}/{PUBLISHED_N_CLEAN})")
    print(f"  [{'ok ' if n_err == 0 else 'FAIL'}] errors {n_err}")

    print("\n=== detection: clean reference -> contaminated reference ===")
    print(f"  {'arm':26}{'clean':>10}{'contaminated':>15}{'bal.acc c->k':>16}")
    print("  " + "-" * 66)
    for a, label in (("conv", "strict golden file"),
                     ("tol", "tolerant golden file"),
                     ("contract", "contracts")):
        c, k = arms[a]["clean"], arms[a]["contam"]
        print(f"  {label:26}{c['n_detected']:>10}{k['n_detected']:>15}"
              f"{str(c['balanced_accuracy']) + ' -> ' + str(k['balanced_accuracy']):>16}")
    print("\n  false escalations on clean runs (contamination costs specificity too):")
    for a, label in (("conv", "strict golden file"),
                     ("tol", "tolerant golden file"),
                     ("contract", "contracts")):
        c, k = arms[a]["clean"], arms[a]["contam"]
        print(f"    {label:26}{c['n_false_escalations']:>4} -> "
              f"{k['n_false_escalations']:>4}  of {c['n_clean_cases']}")

    pg = per_gate(df, "contract")
    fe_pg = per_gate_fe(df, "contract")
    print("\n=== contracts, per gate: prediction 2 lives or dies here ===")
    print(f"  {'gate':6}{'n':>5}{'clean':>8}{'contam':>9}{'lost':>7}   hash anchor")
    for gate, d in pg.items():
        print(f"  {gate:6}{d['n']:>5}{d['detected_clean_ref']:>8}"
              f"{d['detected_contaminated_ref']:>9}{d['lost']:>7}   "
              f"{'yes' if d['has_hash_anchor'] else 'no'}")

    print()
    print("=== contracts, false escalations per gate: the mechanism ===")
    print(f"  {'gate':6}{'clean n':>9}{'FE clean':>10}{'FE contam':>11}   hash anchor")
    for gate, d in fe_pg.items():
        print(f"  {gate:6}{d['n_clean']:>9}{d['fe_clean_ref']:>10}"
              f"{d['fe_contaminated_ref']:>11}   "
              f"{'yes' if d['has_hash_anchor'] else 'no'}")

    summary = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "executes": "M3.1",
        "pre_registration": "ADR-AIIM-223",
        "design": "paired at the case level; three detectors evaluated against "
                  "two references over identical workspace bytes; both "
                  "references contaminated symmetrically",
        "base_seed": BASE_SEED,
        "reference_base_seed": REF_BASE_SEED,
        "noise_seed": NOISE_SEED,
        "reference_noise_seed": REF_NOISE_SEED,
        "workspace_condition": CONDITION_FOR_WORKSPACES,
        "validity_gate": validity,
        "honesty_gate": {
            "requirement": "arm CLEAN reproduces ADR-AIIM-202 exactly",
            "per_arm": gate_detail,
            "case_counts_match": counts_ok,
            "n_errors": n_err,
            "passed": bool(gate_ok and counts_ok and n_err == 0),
        },
        "arms": arms,
        "contracts_per_gate": pg,
        "false_escalations_per_gate": {a: per_gate_fe(df, a) for a in ARMS},
        "tolerant_golden_per_gate": per_gate(df, "tol"),
        "n_reference_injections_vacuous": int(df.ref_vacuous.sum()),
        "mixture_note": "Any contamination fraction p is the linear "
                        "combination p*contaminated + (1-p)*clean; the "
                        "intermediate points of the approved plan were "
                        "withdrawn in ADR-AIIM-223 because they interpolate a "
                        "straight line.",
    }
    with open(out / "contaminated_reference.json", "w", encoding="utf-8",
              newline="\n") as f:
        json.dump(summary, f, indent=1, ensure_ascii=False)
        f.write("\n")

    if not summary["honesty_gate"]["passed"]:
        print("\nHONESTY GATE FAILED (ADR-AIIM-223): the contaminated result is")
        print("written for inspection but MUST NOT be reported.")
    print(f"\nwritten -> {out}")
    if not args.keep_workspaces:
        shutil.rmtree(ws_root, ignore_errors=True)
    return 0 if summary["honesty_gate"]["passed"] else 3


if __name__ == "__main__":
    raise SystemExit(main())

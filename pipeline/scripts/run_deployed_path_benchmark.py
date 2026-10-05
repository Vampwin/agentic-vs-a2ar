"""run_deployed_path_benchmark.py -- E1: does the benchmark measure the system,
or only its evaluators? Pre-registered in ADR-AIIM-236.

ADR-AIIM-226 closed six defects and the seeded-defect benchmark did not move,
because the harness calls the evaluators directly and never reaches the
orchestrator, where three of the four defects lived. Detection 96.9 % and SEER
0.031 are therefore properties of the *evaluators*. This measures the system.

Four arms, paired at the case level, same cases, same evaluator functions:

    A  DIRECT        harness contract dict (21 checks), evaluator called direct
    B  DIRECT-YAML   frozen YAML contract (55 checks), evaluator called direct
    C  DEPLOYED-v1   frozen YAML, through orchestrator + engine + report lock
    D  DEPLOYED-v2   the same, on the patched tree

A->B isolates the contract definition, B->C the execution path, C->D the repair.

### The case generation is driven, not reimplemented

Every arm calls `pilot_runner.run_pilot` with the frozen `PilotConfig` and swaps
only `_dispatch_evaluator`, which is the single seam between the harness and the
detector. Workspace construction, fault injection, `_case_seed`, the vacuity
rule and the telemetry are therefore *the published code*, unmodified. That is
what makes the honesty gate meaningful: arm A is not a re-creation of the
published run, it **is** the published path.

Rewriting the case loop instead would have made the gate decorative, and the
project has a precedent for exactly that failure -- a v2 run that silently
executed v1 and reported "no change" because only one of two path idioms was
rewired.

### Containment, not `evaluator_state`

The harness records what the evaluator returned. ADR-AIIM-236 asks the stricter
question:

    contained  iff  the workflow halted  OR  the artifact is not reportable
    escaped    otherwise

`ReportLock.can_report` requires a PASS contract record or a human adjudication
record, so an ESCALATE artifact is contained by the lock although nothing
halted. That matters because `_escalation_pending` is set by the orchestrator
and read nowhere in the tree. Both definitions are recorded per case; neither is
privileged, which is the convention `_run_one_case` already follows.

### Arms C and D run in separate processes, deliberately

Module resolution decides which tree arm D actually executes, and an in-process
`sys.path` edit leaves already-imported modules cached -- the exact way a v2 run
silently becomes a v1 run. Each arm is therefore a subprocess with `PYTHONPATH`
set to its tree, and **each asserts at startup that the orchestrator it imported
lives under the tree it was asked for**, aborting if not.

Nothing here is analysed beyond the pre-registered gradings. It records per-case
outcomes and stops.

    python run_deployed_path_benchmark.py --self-test
    python run_deployed_path_benchmark.py --coverage-only
    python run_deployed_path_benchmark.py --out 20_new_runs/deployed_path_benchmark_<date>
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml

ROOT = Path(__file__).resolve().parents[2]
FROZEN_TREE = ROOT / "10_new_system"
PATCHED_TREE = ROOT / "11_patched_pipeline_2026-09-22"
PATCHED_V3 = ROOT / "12_patched_pipeline_2026-09-23"

GATES = ["G0", "G1", "G4", "G6", "G8"]
CONTRACT_FILE = {"G0": "G0_input.yaml", "G1": "G1_standardization.yaml",
                 "G4": "G4_scoring.yaml", "G6": "G6_selection.yaml",
                 "G8": "G8_reporting.yaml"}

ARMS = {
    "A": {"name": "DIRECT",      "contract": "harness", "path": "direct",   "tree": "frozen"},
    "B": {"name": "DIRECT-YAML", "contract": "yaml",    "path": "direct",   "tree": "frozen"},
    "C": {"name": "DEPLOYED-v1", "contract": "yaml",    "path": "deployed", "tree": "frozen"},
    "D": {"name": "DEPLOYED-v2", "contract": "yaml",    "path": "deployed", "tree": "patched"},
    "E": {"name": "DEPLOYED-v3", "contract": "yaml",    "path": "deployed", "tree": "patched_v3"},
}

# ADR-AIIM-238: v3 adds eligible_set_completeness and scaffold_uniqueness to the
# deployed G6 contract, which the evaluator already implemented and no contract
# named. P5 predicts arm E reaches 128 of 128; P6 that clean-case false
# escalation stays at 35 of 55.
PREDICTED_E = {"detected": 128, "false_escalation": 35}

# What the published path produces, taken from the artifact rather than from
# prose: seeded_defect_recount_2026-08-29/mutation_score_by_gate_recounted.json
# states "124 cases are detected before and after", and 124/128 is the reported
# detection of 96.9 % with SEER 4/128 = 0.031.
#
# An earlier revision of this file expected 128 detected, read out of
# ADR-AIIM-224's "128/55, 128/0, 124/0". That line reports TWO detectors -- the
# tolerant golden-file comparator at 128/0 and the contracts at 124/0 -- and the
# comparator's number was taken for the contracts'. The honesty gate caught it
# on the first run, which is the only reason it is a footnote and not a result.
PUBLISHED = {"n_fault_nonvacuous": 128, "detected": 124, "strict": 124,
             "n_clean": 55, "false_escalation": 0}

sys.stdout.reconfigure(encoding="utf-8", errors="replace")


# ── the static half: certain, and needs nothing executed ─────────────────────

def check_coverage(tree: Path = FROZEN_TREE) -> dict:
    """How many of each gate's deployed checks the benchmark actually exercises.

    Read rather than executed, so no later step can affect it. This is the
    result ADR-AIIM-236's amendment records.
    """
    import ast
    src = (tree / "fault_injection" / "pilot_runner.py").read_text(encoding="utf-8")
    pilot = None
    for node in ast.parse(src).body:
        if isinstance(node, ast.AnnAssign) and getattr(node.target, "id", "") == "_PILOT_CONTRACT":
            pilot = ast.literal_eval(ast.unparse(node.value).replace("_N_COMPOUNDS", "236"))
    if pilot is None:
        raise SystemExit("ABORT: _PILOT_CONTRACT not found in pilot_runner.py")

    out: dict[str, Any] = {"by_gate": {}, "harness_total": 0, "deployed_total": 0}
    for g in GATES:
        y = yaml.safe_load((tree / "contracts" / CONTRACT_FILE[g]).read_text(encoding="utf-8"))
        yids = [c.get("id") for c in y.get("checks", [])]
        hids = [c.get("id") for c in pilot[g]["checks"]]
        out["by_gate"][g] = {
            "harness_checks": sorted(hids), "deployed_checks": sorted(yids),
            "n_harness": len(hids), "n_deployed": len(yids),
            "never_exercised": sorted(set(yids) - set(hids)),
            "harness_only": sorted(set(hids) - set(yids)),
        }
        out["harness_total"] += len(hids)
        out["deployed_total"] += len(yids)
    out["never_exercised_total"] = sum(len(v["never_exercised"]) for v in out["by_gate"].values())
    out["harness_only_total"] = sum(len(v["harness_only"]) for v in out["by_gate"].values())
    heldout = {"H01": ("G0", "name_structure_mapping"),
               "H05": ("G8", "no_unregistered_claims")}
    out["heldout_targets_in_benchmark"] = {
        case: {"gate": g, "check": c,
               "in_deployed_contract": c in out["by_gate"][g]["deployed_checks"],
               "exercised_by_benchmark": c in out["by_gate"][g]["harness_checks"]}
        for case, (g, c) in heldout.items()}
    out["note"] = ("H02 (G2) and H03 (G3) target gates the benchmark does not "
                   "cover at all, so they have no harness contract to compare.")
    return out


# ── containment, per ADR-AIIM-236 ────────────────────────────────────────────

def contained(halted: bool, reportable: bool | None) -> bool:
    """Halted, or the artifact cannot be reported."""
    return bool(halted or (reportable is False))


# ── one arm, in its own process ──────────────────────────────────────────────

def _assert_tree(tree: Path) -> None:
    """Refuse to run if module resolution did not land in the intended tree."""
    import kdss_audit.orchestrator as orch
    got = Path(orch.__file__).resolve()
    want = tree.resolve()
    if want not in got.parents:
        raise SystemExit(
            f"ABORT: asked for {want.name} but imported orchestrator from {got}. "
            f"A v2 run that silently executes v1 is the failure this check exists "
            f"to prevent.")
    print(f"tree verified: orchestrator imported from {got.parent.parent.name}")


def run_arm(arm: str, out: Path) -> int:
    cfg = ARMS[arm]
    tree = {"patched": PATCHED_TREE, "patched_v3": PATCHED_V3}.get(
        cfg["tree"], FROZEN_TREE)
    if cfg["path"] == "deployed":
        _assert_tree(tree)

    from fault_injection import pilot_runner as PR

    original_dispatch = PR._dispatch_evaluator
    side: dict[str, dict] = {}
    contracts_dir = tree / "contracts"
    yaml_cache: dict[str, dict] = {}

    def yaml_contract(gate: str) -> dict:
        if gate not in yaml_cache:
            yaml_cache[gate] = yaml.safe_load(
                (contracts_dir / CONTRACT_FILE[gate]).read_text(encoding="utf-8"))
        return yaml_cache[gate]

    def dispatch(gate, ws, clean_hashes=None):
        key = Path(ws).name
        if arm == "A":
            res, avail = original_dispatch(gate, ws, clean_hashes)
            # The evaluators return "failed"/"passed"; only ContractResult uses
            # checks_failed. Reading the wrong key gave an empty list and made a
            # per-check claim unsupportable on the first run of this script.
            side[key] = {"gate": gate, "state": (res or {}).get("state"),
                         "checks_failed": list((res or {}).get("failed", [])),
                         "checks_passed": list((res or {}).get("passed", [])),
                         "halted": None, "reportable": None, "contained": None}
            return res, avail

        inputs = PR._get_evaluator_inputs(gate, ws, clean_hashes)

        if arm == "B":
            import importlib
            mod = importlib.import_module(f"evaluators.{gate}_evaluator")
            try:
                res = mod.evaluate(contract=yaml_contract(gate), inputs=inputs,
                                   registry=None)
            except NotImplementedError:
                side[key] = {"gate": gate, "error": "NotImplementedError"}
                return None, False
            except Exception as exc:                          # noqa: BLE001
                side[key] = {"gate": gate, "error": f"{type(exc).__name__}: {exc}"}
                return {"state": "ERROR", "error": str(exc)}, False
            side[key] = {"gate": gate, "state": res.get("state"),
                         "checks_failed": list(res.get("failed", [])),
                         "checks_passed": list(res.get("passed", [])),
                         "halted": None, "reportable": None, "contained": None}
            return res, True

        # arms C and D: the deployed path
        from kdss_audit.artifact_registry import ArtifactRegistry
        from kdss_audit.orchestrator import WorkflowOrchestrator
        run_dir = Path(ws) / "_deployed"
        run_dir.mkdir(parents=True, exist_ok=True)
        registry = ArtifactRegistry(run_dir / "artifact_registry.jsonl")
        orch = WorkflowOrchestrator(run_dir=run_dir, contracts_dir=contracts_dir,
                                    config_path=run_dir / "config.yaml",
                                    registry=registry)
        try:
            sr = orch.run_stage(stage_name=gate, contract_file=CONTRACT_FILE[gate],
                                tool_fn=None, tool_outputs=inputs)
        except Exception as exc:                              # noqa: BLE001
            side[key] = {"gate": gate, "error": f"{type(exc).__name__}: {exc}",
                         "halted": None, "reportable": None, "contained": None}
            return {"state": "ERROR", "error": str(exc)}, False

        state = getattr(sr.contract_result.state, "value", str(sr.contract_result.state))
        lock = orch.build_report_lock(run_dir / "report.lock")
        reportable = (all(lock.can_report(a) for a in sr.artifacts_produced)
                      if sr.artifacts_produced else None)
        side[key] = {"gate": gate, "state": state,
                     "checks_failed": list(sr.contract_result.checks_failed),
                     "halted": bool(orch._halted), "reportable": reportable,
                     "contained": contained(bool(orch._halted), reportable),
                     "n_artifacts": len(sr.artifacts_produced)}
        return {"state": state,
                "checks_failed": list(sr.contract_result.checks_failed)}, True

    PR._dispatch_evaluator = dispatch
    try:
        config = PR.PilotConfig(conditions=["B"])
        with tempfile.TemporaryDirectory(prefix=f"e1_{arm}_") as td:
            results = PR.run_pilot(config=config, workspace_root=Path(td))
    finally:
        PR._dispatch_evaluator = original_dispatch

    # Joining the side channel back to the outcomes has to be exact. Parsing the
    # index out of `instance_id` is not: each injector formats its own, and for
    # F17/F18 it carries no index at all, which silently left 11 cases unjoined
    # and 9 more mismatched on the first run of this script.
    #
    # `run_pilot` appends outcomes in a fixed order -- per family, fault indices
    # 0..n-1, then clean indices 0..n_clean-1, conditions innermost -- and
    # `_run_one_case` names the workspace
    # f"{family}_{index:03d}_{condition}_{case_type}". With one condition the
    # ordinal position within (family, case_type) IS the index, so the key is
    # reconstructed from position rather than from a parsed string.
    rows = []
    seen: dict[tuple[str, str], int] = {}
    for o in results.outcomes:          # run_pilot returns dicts, not dataclasses
        k = (o["family_id"], o["case_type"])
        idx = seen.get(k, 0)
        seen[k] = idx + 1
        ws_name = f"{o['family_id']}_{idx:03d}_{o['condition']}_{o['case_type']}"
        s = side.get(ws_name, {})
        rows.append({
            "arm": arm, "arm_name": cfg["name"], "family_id": o["family_id"],
            "instance_id": o["instance_id"], "case_type": o["case_type"],
            "condition": o["condition"], "vacuous": bool(o.get("vacuous", False)),
            "evaluator_available": o["evaluator_available"],
            "evaluator_state": o["evaluator_state"],
            "escaped_strict": o["escaped"], "escaped_detected": o["escaped_detected"],
            "error": o.get("error"),
            "checks_failed": s.get("checks_failed", []),
            "checks_passed": s.get("checks_passed", []),
            "halted": s.get("halted"), "reportable": s.get("reportable"),
            "contained": s.get("contained"),
        })

    unjoined = [r["instance_id"] for r in rows
                if cfg["path"] == "deployed" and r["halted"] is None
                and not r["error"]]
    if unjoined:
        raise SystemExit(
            f"ABORT: {len(unjoined)} cases did not join to the side channel, "
            f"e.g. {unjoined[:3]}. A partial join silently understates "
            f"containment; it is not reported.")

    out.mkdir(parents=True, exist_ok=True)
    (out / f"arm_{arm}_cases.json").write_text(
        json.dumps({"arm": arm, **cfg, "tree": str(tree),
                    "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                    "n_cases": len(rows), "cases": rows}, indent=2), encoding="utf-8")
    print(f"arm {arm} ({cfg['name']}): {len(rows)} cases -> arm_{arm}_cases.json")
    return 0


# ── summary and the pre-registered gradings ──────────────────────────────────

def summarise(out: Path) -> dict:
    arms = {}
    for a in ARMS:
        p = out / f"arm_{a}_cases.json"
        if p.exists():
            arms[a] = json.loads(p.read_text(encoding="utf-8"))
    s: dict[str, Any] = {"by_arm": {}}
    for a, d in arms.items():
        faults = [c for c in d["cases"] if c["case_type"] == "fault" and not c["vacuous"]]
        clean = [c for c in d["cases"] if c["case_type"] == "clean"]
        s["by_arm"][a] = {
            "arm_name": d["name"],
            "n_fault_nonvacuous": len(faults), "n_clean": len(clean),
            "detected": sum(1 for c in faults if not c["escaped_detected"]),
            "strict": sum(1 for c in faults if not c["escaped_strict"]),
            "false_escalation": sum(1 for c in clean if c["escaped_detected"]),
            "contained": sum(1 for c in faults if c["contained"] is True),
            "clean_contained": sum(1 for c in clean if c["contained"] is True),
            "errors": sum(1 for c in d["cases"] if c["error"]),
        }
    a = s["by_arm"].get("A", {})
    s["honesty_gate"] = {
        "expected": PUBLISHED,
        "observed": {k: a.get(k) for k in PUBLISHED},
        "pass": all(a.get(k) == v for k, v in PUBLISHED.items()) if a else False,
    }
    if "C" in s["by_arm"]:
        c, b = s["by_arm"]["C"], s["by_arm"].get("B", {})
        s["P1_deployed_contains_at_least_direct"] = (
            c["contained"] >= a.get("detected", 0) if a else None)
        s["P3_lock_adds_containment"] = c["contained"] > c["strict"]
    if "C" in s["by_arm"] and "D" in s["by_arm"]:
        s["P2_repair_changes_nothing"] = (
            s["by_arm"]["C"]["contained"] == s["by_arm"]["D"]["contained"]
            and s["by_arm"]["C"]["detected"] == s["by_arm"]["D"]["detected"])
    if "E" in s["by_arm"]:
        e = s["by_arm"]["E"]
        s["P5_v3_recovers_every_case"] = e["detected"] == PREDICTED_E["detected"]
        s["P6_false_escalation_unchanged"] = (
            e["false_escalation"] == PREDICTED_E["false_escalation"])
    if "B" in s["by_arm"] and "A" in s["by_arm"]:
        s["P4_frozen_contract_raises_clean_false_escalation"] = (
            s["by_arm"]["B"]["false_escalation"] > s["by_arm"]["A"]["false_escalation"])
    return s


# ── self-test ────────────────────────────────────────────────────────────────

def self_test() -> int:
    ok = []

    def check(name, cond):
        ok.append((name, bool(cond)))
        print(("  PASS  " if cond else "  FAIL  ") + name)

    cov = check_coverage()
    check("the 21/55 split is what ADR-AIIM-236 recorded",
          cov["harness_total"] == 21 and cov["deployed_total"] == 55)
    check("two harness checks exist in no deployed contract",
          cov["harness_only_total"] == 2)
    h = cov["heldout_targets_in_benchmark"]
    check("H01's target check is deployed but never exercised",
          h["H01"]["in_deployed_contract"] and not h["H01"]["exercised_by_benchmark"])
    check("H05's target check is deployed but never exercised",
          h["H05"]["in_deployed_contract"] and not h["H05"]["exercised_by_benchmark"])

    check("containment: halted counts", contained(True, True))
    check("containment: unreportable counts", contained(False, False))
    check("containment: reportable and unhalted escapes", not contained(False, True))
    check("containment: unknown reportability is not containment",
          not contained(False, None))

    check("all five arms are defined", set(ARMS) == {"A", "B", "C", "D", "E"})
    check("every mapped contract file exists in the frozen tree",
          all((FROZEN_TREE / "contracts" / CONTRACT_FILE[g]).is_file() for g in GATES))
    check("every mapped contract file exists in the patched tree",
          all((PATCHED_TREE / "contracts" / CONTRACT_FILE[g]).is_file() for g in GATES))
    check("the v3 tree exists for arm E", PATCHED_V3.is_dir())
    import yaml as _y
    _v3 = _y.safe_load((PATCHED_V3 / "contracts" / "G6_selection.yaml").read_text(encoding="utf-8"))
    _ids = [c["id"] for c in _v3["checks"]]
    check("v3 G6 names both previously unnamed checks",
          "eligible_set_completeness" in _ids and "scaffold_uniqueness" in _ids)
    _v1 = _y.safe_load((FROZEN_TREE / "contracts" / "G6_selection.yaml").read_text(encoding="utf-8"))
    check("the frozen G6 contract is untouched",
          "eligible_set_completeness" not in [c["id"] for c in _v1["checks"]]
          and _v1["version"] == "1.3.0")

    sys.path.insert(0, str(FROZEN_TREE))
    from fault_injection import pilot_runner as PR
    check("the seam this harness swaps still exists",
          callable(getattr(PR, "_dispatch_evaluator", None)))
    check("the input builder it reuses still exists",
          callable(getattr(PR, "_get_evaluator_inputs", None)))
    cfg = PR.PilotConfig(conditions=["B"])
    n_inst = sum(len(PR.FAULT_REGISTRY[f].spec.param_grid) for f in cfg.pilot_families)
    check("the frozen config still yields 133 fault instances", n_inst == 133)
    check("the frozen config still yields 55 clean cases",
          cfg.n_clean_per_family * len(cfg.pilot_families) == 55)

    # workspace-name reconstruction must match what _run_one_case builds
    check("PUBLISHED matches the recounted artifact, not ADR prose",
          PUBLISHED["detected"] == 124 and PUBLISHED["strict"] == 124)

    print(f"\n{sum(v for _, v in ok)} of {len(ok)} checks pass")
    return 0 if all(v for _, v in ok) else 1


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "20_new_runs" /
                                         "deployed_path_benchmark_2026-09-23"))
    ap.add_argument("--arm", choices=sorted(ARMS),
                    help="run one arm in this process (the driver sets PYTHONPATH)")
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--coverage-only", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        return self_test()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    if args.arm:
        return run_arm(args.arm, out)

    cov = check_coverage()
    (out / "check_coverage.json").write_text(
        json.dumps({"generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                    "adr": "ADR-AIIM-236", **cov}, indent=2), encoding="utf-8")
    print(f"static: the benchmark exercises {cov['harness_total']} of "
          f"{cov['deployed_total']} checks these five gates define; "
          f"{cov['never_exercised_total']} never exercised, "
          f"{cov['harness_only_total']} exist only in the harness")
    if args.coverage_only:
        return 0

    for a in sorted(ARMS):
        tree = {"patched": PATCHED_TREE, "patched_v3": PATCHED_V3}.get(
            ARMS[a]["tree"], FROZEN_TREE)
        env = dict(os.environ, PYTHONPATH=str(tree))
        print(f"\n--- arm {a} ({ARMS[a]['name']}), PYTHONPATH={tree.name} ---")
        rc = subprocess.run([sys.executable, str(Path(__file__).resolve()),
                             "--arm", a, "--out", str(out)], env=env).returncode
        if rc != 0:
            print(f"arm {a} exited {rc}; stopping. Nothing is summarised from a "
                  f"partial set of arms.")
            return rc

    s = summarise(out)
    (out / "summary.json").write_text(json.dumps(
        {"generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
         "adr": "ADR-AIIM-236", **s}, indent=2), encoding="utf-8")
    g = s["honesty_gate"]
    print("\nhonesty gate:", "PASS" if g["pass"] else "FAIL")
    print("  expected", g["expected"])
    print("  observed", g["observed"])
    if not g["pass"]:
        print("\nADR-AIIM-236: arm A did not reproduce the published counts. "
              "THE RUN IS VOID and no other number from it may be reported.")
        return 1
    for k in ("P1_deployed_contains_at_least_direct", "P2_repair_changes_nothing",
              "P3_lock_adds_containment",
              "P4_frozen_contract_raises_clean_false_escalation",
              "P5_v3_recovers_every_case", "P6_false_escalation_unchanged"):
        if k in s:
            print(f"  {k}: {s[k]}")
    print(f"\nwrote {out / 'summary.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

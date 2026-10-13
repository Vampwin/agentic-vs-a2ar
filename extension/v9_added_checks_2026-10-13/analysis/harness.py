"""ADR-AIIM-427/428: workspaces, held-out fault sets and the evaluation paths for the four added checks. Nothing in 10_new_system is modified.

Paths: F = frozen evaluator, frozen contract, called directly (arm B of ADR-AIIM-236); P = the same with the added checks (wrapper evaluator
`v9ac.*_plus`, direct); D = the deployed path (orchestrator, contract engine, report lock) on the copied contracts in `contracts_plus/`.
The workspaces of G2 and G7 are those of ADR-AIIM-416 (`ext.py`) and those of G5 are those of ADR-AIIM-420 (`ext2.py`); both are imported."""
from __future__ import annotations

import hashlib
import importlib
import json
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

HERE = Path(__file__).parent
OUT = HERE.parent
RUNS = OUT.parent
FROZEN = RUNS.parent / "10_new_system"
for p in (FROZEN, HERE, RUNS / "v9_g2_g7_extension_2026-10-12" / "analysis", RUNS / "v9_g3_g5_extension_2026-10-12" / "analysis"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))
import ext      # noqa: E402  (G2 and G7 workspaces and injectors)
import ext2     # noqa: E402  (G5 workspaces)
from v9ac import checks      # noqa: E402

PLUS = OUT / "contracts_plus"
CONTRACT_FILE = {"G2": "G2_preparation.yaml", "G7": "G7_md.yaml", "G5a": "G5a_validity.yaml", "G5b": "G5b_performance.yaml"}
PLUS_MODULE = {"G2": "v9ac.G2_plus", "G7": "v9ac.G7_plus", "G5a": "v9ac.G5a_plus"}
ADDED = {"G2": ["protonation_matches_record"], "G7": ["trajectory_covers_plan", "binding_site_matches_receptor"],
         "G5a": ["benchmark_ids_not_exposed"]}
ADDED_DESC = {
    "protonation_matches_record": "Every histidine of the prepared receptor is in the protonation state (the ring hydrogen it carries) that the "
                                  "preparation step recorded. A missing record fails.",
    "trajectory_covers_plan": "Every per-frame series reaches at least 95% of the planned production length.",
    "binding_site_matches_receptor": "The residue indices that define the binding site in the results table equal the receptor C-alpha residues within "
                                     "10 A of the reference ligand, recomputed from the prepared receptor and the reference ligand.",
    "benchmark_ids_not_exposed": "No file read before the gate, other than the registered ranking inputs, contains a compound ID of the benchmark.",
}


def write_plus_contracts() -> None:
    PLUS.mkdir(exist_ok=True)
    for gate, fn in CONTRACT_FILE.items():
        c = yaml.safe_load((FROZEN / "contracts" / fn).read_text(encoding="utf-8"))
        if gate in PLUS_MODULE:
            c["evaluator"] = PLUS_MODULE[gate]
            for cid in ADDED[gate]:
                c["checks"].append({"id": cid, "severity": "critical", "description": ADDED_DESC[cid]})
            c["version"] = str(c.get("version", "1.0.0")) + "+v9ac"
        (PLUS / fn).write_text(yaml.safe_dump(c, sort_keys=False, allow_unicode=True), encoding="utf-8")


# ------------------------------------------------------------------------------------------------ evaluation paths
def run_frozen(gate: str, inputs: dict) -> dict:
    return ext.run_direct(gate, inputs)


def run_plus(gate: str, inputs: dict) -> dict:
    contract = yaml.safe_load((FROZEN / "contracts" / CONTRACT_FILE[gate]).read_text(encoding="utf-8"))
    mod = importlib.import_module(PLUS_MODULE[gate])
    try:
        r = mod.evaluate(contract=contract, inputs=inputs, registry=None)
    except Exception as exc:                                     # noqa: BLE001
        return {"state": "ERROR", "error": f"{type(exc).__name__}: {exc}", "checks_failed": [], "checks_passed": []}
    return {"state": r.get("state"), "checks_failed": list(r.get("failed", [])), "checks_passed": list(r.get("passed", []))}


def run_deployed_plus(gate: str, inputs: dict, ws: Path) -> dict:
    from kdss_audit.artifact_registry import ArtifactRegistry
    from kdss_audit.orchestrator import WorkflowOrchestrator
    run_dir = ws / "_deployed_plus"
    run_dir.mkdir(parents=True, exist_ok=True)
    registry = ArtifactRegistry(run_dir / "artifact_registry.jsonl")
    orch = WorkflowOrchestrator(run_dir=run_dir, contracts_dir=PLUS, config_path=run_dir / "config.yaml", registry=registry)
    try:
        sr = orch.run_stage(stage_name=gate, contract_file=CONTRACT_FILE[gate], tool_fn=None, tool_outputs=inputs)
    except Exception as exc:                                     # noqa: BLE001
        return {"state": "ERROR", "error": f"{type(exc).__name__}: {exc}", "halted": None, "reportable": None, "contained": None, "checks_failed": []}
    state = getattr(sr.contract_result.state, "value", str(sr.contract_result.state))
    lock = orch.build_report_lock(run_dir / "report.lock")
    reportable = (all(lock.can_report(a) for a in sr.artifacts_produced) if sr.artifacts_produced else None)
    halted = bool(orch._halted)
    return {"state": state, "checks_failed": list(sr.contract_result.checks_failed), "halted": halted, "reportable": reportable,
            "contained": bool(halted or (reportable is False))}


# ------------------------------------------------------------------------------------------------ G2 and G7 workspaces
def build_g2_plus(ws: Path, name: str, receptor: Path, cutoff: float) -> dict:
    inp = ext.build_g2(ws, name, receptor, cutoff)
    st = checks.his_states(Path(inp["prepared_receptor_4eiy_pdbqt"]))
    inp["protonation_decision"] = {k: v["state"] for k, v in st.items()}       # what a compliant preparation step records
    return inp


def build_g7_plus(ws: Path, sid: str, reps: list[str]) -> dict:
    inp = ext.build_g7(ws, sid, reps)
    shutil.copyfile(ext.G2_DIR / "4EIY_prepared.pdbqt", ws / "receptor_4EIY_prepared.pdbqt")
    shutil.copyfile(ext.G2_DIR / "ZMA_reference_4EIY.sdf", ws / "reference_ligand.sdf")
    inp["receptor_pdbqt"] = str(ws / "receptor_4EIY_prepared.pdbqt")
    inp["reference_ligand_sdf"] = str(ws / "reference_ligand.sdf")
    return inp


# ------------------------------------------------------------------------------------------------ G5 workspace with leak modes
def build_g5_plus(ws: Path, variant: str, leak_ids: list[str] | None = None, leak_mode: str = "logged", context_ids: list[str] | None = None):
    """Like ext2.build_g5 (copied, one parameter added). leak_mode: 'logged' (the leak file is in the access log), 'unlogged' (the file exists and
    shapes the ranking but no access record names it), 'benchmark' (the ranking step read the benchmark file). context_ids: a file of IDs that are
    not benchmark compounds, read by the ranking step (a benign control). The registered ranking input is the scores file."""
    ws.mkdir(parents=True, exist_ok=True)
    name, src, col, sign = next(v for v in ext2.G5_VARIANTS if v[0] == variant)
    df = pd.read_csv(src)
    bench = ext2.bench_frame().copy()
    bench_path, ranked_path = ws / "chembl_benchmark.csv", ws / "ranked_scores.csv"
    id_hash = ext2.hash_ids(bench["compound_id"].tolist())
    scores = (sign * df.set_index("compound_id").loc[bench["compound_id"], col].to_numpy(float)).copy()
    if leak_ids:
        scores[bench["compound_id"].isin(leak_ids).to_numpy()] = scores.max() + 1.0
    bench.to_csv(bench_path, index=False)
    pd.DataFrame({"compound_id": bench["compound_id"], "score": scores}).to_csv(ranked_path, index=False)
    m = ext2.metrics(scores, bench["label"].to_numpy())
    reg = ext2.new_registry(ws / "_registry" / "artifact_registry.jsonl")
    reg.record_benchmark_role_declaration("Option_A_development", "Adisorn Sirichotejirakul", "the set was used during legacy method development")
    ext2.log_read(reg, "run_ranking", ranked_path, "ranking of the benchmark compounds by score")
    if leak_ids:
        leak_path = ws / "ranking_context_known_ids.txt"
        leak_path.write_text("\n".join(leak_ids))
        if leak_mode == "logged":
            ext2.log_read(reg, "run_ranking", leak_path, "ranking step read the known IDs")
        elif leak_mode == "benchmark":
            ext2.log_read(reg, "run_ranking", bench_path, "ranking step read the benchmark labels")
    if context_ids:
        cp = ws / "ranking_context_unrelated_ids.txt"
        cp.write_text("\n".join(context_ids))
        ext2.log_read(reg, "run_ranking", cp, "ranking step read a list of compounds that are not in the benchmark")
    time.sleep(0.02)
    gate_start = time.time()
    time.sleep(0.02)
    ext2.log_read(reg, "G5b", bench_path, "G5b reads the labels at the gate")
    g5a = {"chembl_benchmark_csv": str(bench_path), "benchmark_role_declaration": True,
           "g5a_evaluator_source_hash": hashlib.sha256((FROZEN / "evaluators" / "G5a_evaluator.py").read_bytes()).hexdigest().upper(),
           "kdss_audit_module_hashes": ext2.kdss_hashes(), "benchmark_id_list_hash_registered": id_hash,
           "bootstrap_known_answer_test": {"pairs_resampled_together": True}, "bootstrap_n_resamples": 2000,
           "ef1pct_count": int(bench["label"].to_numpy()[np.argsort(-scores, kind="stable")[:max(1, int(np.ceil(0.01 * len(bench))))]].sum()),
           "_gate_start_time": gate_start, "registered_ranking_inputs": [str(ranked_path)]}
    g5b = {"auc_roc": m.get("auc_roc"), "ef_1pct": m.get("ef_1pct"), "bedroc20": m.get("bedroc20"), "registered_auc_min_threshold": 0.60,
           "benchmark_role": "Option_A_development"}
    return {"G5a": g5a, "G5b": g5b, "metrics": m}, reg


def _g5_chain(run_a, run_b, g5):
    a = run_a(g5["G5a"])
    b = run_b(g5["G5b"]) if a["state"] == "PASS" else None
    return {"G5a": a, "G5b": b}


def run_frozen_g5(g5: dict, reg) -> dict:
    return ext2.run_direct_g5(g5, reg)


def run_plus_g5(g5: dict, reg) -> dict:
    import evaluators.G5b_evaluator as eb
    import v9ac.G5a_plus as ap
    ca = yaml.safe_load((FROZEN / "contracts" / CONTRACT_FILE["G5a"]).read_text(encoding="utf-8"))
    try:
        ra = ap.evaluate(contract=ca, inputs=g5["G5a"], registry=reg)
    except Exception as exc:                                     # noqa: BLE001
        return {"G5a": {"state": "ERROR", "error": f"{type(exc).__name__}: {exc}", "checks_failed": []}, "G5b": None}
    out = {"G5a": {"state": ra["state"], "checks_failed": list(ra["failed"])}, "G5b": None}
    if ra["state"] == "PASS":
        cb = yaml.safe_load((FROZEN / "contracts" / CONTRACT_FILE["G5b"]).read_text(encoding="utf-8"))
        rb = eb.evaluate(contract=cb, inputs=g5["G5b"], registry=reg)
        out["G5b"] = {"state": rb["state"], "checks_failed": list(rb["failed"])}
    return out


def run_deployed_plus_g5(ws: Path, g5: dict) -> dict:
    from kdss_audit.artifact_registry import ArtifactRegistry
    from kdss_audit.orchestrator import WorkflowOrchestrator
    run_dir = ws / "_registry"
    reg = ArtifactRegistry(run_dir / "artifact_registry.jsonl")
    orch = WorkflowOrchestrator(run_dir=run_dir, contracts_dir=PLUS, config_path=run_dir / "config.yaml", registry=reg)

    def stage(name, inputs):
        try:
            sr = orch.run_stage(stage_name=name, contract_file=CONTRACT_FILE[name], tool_fn=None, tool_outputs=inputs)
        except Exception as exc:                                 # noqa: BLE001
            return {"state": "ERROR", "error": f"{type(exc).__name__}: {exc}", "checks_failed": []}
        return {"state": getattr(sr.contract_result.state, "value", str(sr.contract_result.state)),
                "checks_failed": list(sr.contract_result.checks_failed), "halted": bool(orch._halted)}
    a = stage("G5a", g5["G5a"])
    return {"G5a": a, "G5b": stage("G5b", g5["G5b"]) if a["state"] == "PASS" else None}


def frozen_tree_hashes() -> dict:
    return {str(p.relative_to(FROZEN)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(FROZEN.rglob("*")) if p.is_file() and "__pycache__" not in p.parts and p.suffix in (".py", ".yaml", ".yml", ".json")}

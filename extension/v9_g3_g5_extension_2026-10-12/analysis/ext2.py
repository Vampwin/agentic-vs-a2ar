"""ADR-AIIM-420: workspaces, injectors and the two evaluation paths for G3 and G5. Nothing in 10_new_system is modified; the frozen
evaluators, contracts, orchestrator and docking script are imported and read. The injectors implement the parameter grids of the stubs
F07, F08, F11 and F12 of `10_new_system/fault_injection/` outside that tree.

G3 workspaces: the real redocking artifacts of this project (poses, atom map, reference ligand, conformers, receptor and box) and,
where a fault needs it, a real redocking with the frozen multi-conformer Vina script of `run_snp_screen.py`.
G5 workspaces: the real ChEMBL A2AR benchmark of this project (245 compounds, 50 active and 195 inactive) with its real scores, the
metrics computed by the project's own function, and a registry that carries the benchmark-role declaration and the file-access log."""
from __future__ import annotations

import hashlib
import importlib.util
import json
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(r"D:\Agentic VS Workflow")
RUNS = ROOT / "20_new_runs"
FROZEN = ROOT / "10_new_system"
OUT = RUNS / "v9_g3_g5_extension_2026-10-12"
sys.path.insert(0, str(FROZEN))
sys.path.insert(0, str(FROZEN / "scripts"))

SNP = RUNS / "snp_screen_2026-08-04"
G2D, G3D = SNP / "G2_preparation", SNP / "G3_redocking"
CRYSTAL = ROOT / "04_legacy_evidence" / "redocking" / "ZMA_crystal_4EIY.sdf"
CFG = yaml.safe_load((FROZEN / "configs" / "a2a_frozen_config.yaml").read_text(encoding="utf-8"))
SEED, EXH, NPOSES = CFG["vina"]["seed"], CFG["vina"]["exhaustiveness"], CFG["vina"]["n_poses"]
TEMPLATE = CFG["g3_redocking"]["reference_ligand_smiles"]
CONTRACT_FILE = {"G3": "G3_pose.yaml", "G5a": "G5a_validity.yaml", "G5b": "G5b_performance.yaml"}
BOX = json.loads((G2D / "docking_box.json").read_text())
CENTER, SIZE = BOX["centers_by_receptor"]["4EIY"], BOX["size"]


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def snp_screen():
    if not hasattr(snp_screen, "m"):
        snp_screen.m = _load("_run_snp_screen", FROZEN / "scripts" / "run_snp_screen.py")
    return snp_screen.m


def ifp_module():
    if not hasattr(ifp_module, "m"):
        ifp_module.m = _load("_ifp", FROZEN / "scripts" / "ifp_rescore_benchmark.py")
    return ifp_module.m


# ------------------------------------------------------------------------------------------------ registry helpers
def new_registry(path: Path):
    from kdss_audit.artifact_registry import ArtifactRegistry
    if path.exists():
        path.unlink()
    return ArtifactRegistry(path)


def register(reg, artifact_id: str, path: Path, stage: str, desc: str):
    from kdss_audit.artifact_registry import ArtifactEntry, sha256_file
    reg.register_artifact(ArtifactEntry(artifact_id=artifact_id, path=str(path), sha256=sha256_file(path), stage=stage,
                                        description=desc, tool="v9_g3_g5_extension", tool_version="", conda_env="", command="",
                                        random_seed=None))


def log_read(reg, reader: str, path: Path, purpose: str):
    from kdss_audit.artifact_registry import FileAccessEntry
    reg.record_file_access(FileAccessEntry(reader=reader, file_path=str(path), purpose=purpose))


# ------------------------------------------------------------------------------------------------ G3
_LIB = None


def library_smiles() -> list[str]:
    global _LIB
    if _LIB is None:
        led = pd.read_csv(SNP / "G1_standardization" / "snp_source_ledger.csv")
        _LIB = [s for s in led["smiles_canonical"].tolist() if isinstance(s, str)]
    return _LIB


def retro_enrichment() -> dict:
    p = RUNS / "chembl245_ifp_rescoring_2026-08-15" / "ifp_rescoring_summary.json"
    vo = json.loads(p.read_text())["methods"]["vina_only"]
    return {"auc_roc": vo["auc_roc"], "bedroc20": vo["bedroc20"], "protocol_matches_screen": True, "n": vo["n"], "source": str(p)}


def redock(out_sdf: Path, center, seed: int, tag: str) -> Path:
    """The frozen multi-conformer redocking of run_snp_screen.py (same script, receptor, conformers, box size, exhaustiveness, poses),
    with the box centre and the seed as the only parameters. The script's fixed /tmp file names are made unique per call."""
    m = snp_screen()
    script = m._MULTICONF_DOCK_SCRIPT.replace("/tmp/_g3_lig.pdbqt", f"/tmp/_g3_lig_{tag}.pdbqt").replace("/tmp/_g3_out.pdbqt", f"/tmp/_g3_out_{tag}.pdbqt")
    energy = out_sdf.with_suffix(".energies.json")
    out_sdf.parent.mkdir(parents=True, exist_ok=True)
    args = [m._win_to_wsl(G2D / "4EIY_prepared.pdbqt"), m._win_to_wsl(G3D / "ZM241385_conformers.sdf"), m._win_to_wsl(out_sdf),
            m._win_to_wsl(energy), TEMPLATE, str(center[0]), str(center[1]), str(center[2]), str(SIZE[0]), str(SIZE[1]), str(SIZE[2]),
            str(seed), str(EXH), str(NPOSES)]
    proc = m._wsl_python("import sys; sys.argv = ['x', " + ", ".join(repr(x) for x in args) + "]\n" + script, purpose=f"redock {tag}")
    if proc.returncode != 0 or not out_sdf.exists():
        raise RuntimeError(f"redock {tag} failed:\n{proc.stdout}\n{proc.stderr}")
    return out_sdf


def build_g3(ws: Path, poses_src: Path, reference_src: Path = CRYSTAL, atom_map_src: Path = G3D / "atom_map_4eiy.json"):
    """A G3 workspace and its registry: the atom map is registered first, the poses after it, as in the production run."""
    ws.mkdir(parents=True, exist_ok=True)
    poses, ref, amap = ws / "ZMA_redocked_4EIY.sdf", ws / "ZMA_crystal_4EIY.sdf", ws / "atom_map_4eiy.json"
    shutil.copyfile(poses_src, poses)
    shutil.copyfile(reference_src, ref)
    shutil.copyfile(atom_map_src, amap)
    shutil.copyfile(G3D / "ZM241385_conformers.sdf", ws / "ZM241385_conformers.sdf")
    reg = new_registry(ws / "_registry" / "artifact_registry.jsonl")
    register(reg, "G3::atom_map_4eiy", amap, "G3", "pre-docking MCS atom map")
    time.sleep(0.02)
    register(reg, "G3::redocked_poses_4eiy", poses, "G3", "pooled Vina poses")
    inputs = {"reference_ligand_4eiy_sdf": str(ref), "redocked_poses_4eiy": str(poses), "atom_map_4eiy_json": str(amap),
              "docking_input_conformer_sdf": str(ws / "ZM241385_conformers.sdf"), "control_ligand_smiles": [TEMPLATE],
              "screening_library_smiles": library_smiles(), "retrospective_enrichment": retro_enrichment(), "frame_check_4eiy": True}
    return inputs, reg


F07_GRID = [{"offset_angstrom": o, "axis": a} for o in (5, 10, 15, 20, 25) for a in ("x", "y", "z")]
F08_GRID = [{"rotation_degrees": d, "rotation_axis": a} for d in (45, 90, 120, 135, 150, 165, 180) for a in ("x", "y", "z")][:15]


def f07_center(p: dict):
    c = list(CENTER)
    c["xyz".index(p["axis"])] += p["offset_angstrom"]
    return c


def rotate_reference(src: Path, dst: Path, degrees: float, axis: str) -> dict:
    """F08: the reference ligand coordinates rotated by `degrees` about `axis` through the ligand centroid (the public crystal ligand)."""
    from rdkit import Chem
    mol = Chem.MolFromMolFile(str(src), removeHs=False)
    conf = mol.GetConformer()
    xyz = conf.GetPositions()
    c = xyz.mean(axis=0)
    th = np.deg2rad(degrees)
    cth, sth = np.cos(th), np.sin(th)
    R = {"x": np.array([[1, 0, 0], [0, cth, -sth], [0, sth, cth]]), "y": np.array([[cth, 0, sth], [0, 1, 0], [-sth, 0, cth]]),
         "z": np.array([[cth, -sth, 0], [sth, cth, 0], [0, 0, 1]])}[axis]
    new = (xyz - c) @ R.T + c
    for i, p in enumerate(new):
        conf.SetAtomPosition(i, p.tolist())
    w = Chem.SDWriter(str(dst))
    w.write(mol)
    w.close()
    return {"max_atom_shift_A": float(np.max(np.linalg.norm(new - xyz, axis=1)))}


def run_direct_g3(inputs: dict, reg) -> dict:
    import evaluators.G3_evaluator as ev
    contract = yaml.safe_load((FROZEN / "contracts" / CONTRACT_FILE["G3"]).read_text(encoding="utf-8"))
    try:
        r = ev.evaluate(contract=contract, inputs=inputs, registry=reg)
    except Exception as exc:                                     # noqa: BLE001
        return {"state": "ERROR", "error": f"{type(exc).__name__}: {exc}", "checks_failed": []}
    ev_r = r.get("evidence", {}).get("rmsd_rank1_4eiy", {})
    return {"state": r["state"], "checks_failed": list(r["failed"]), "critical": list(r["critical_failures"]),
            "warned": list(r["warned"]), "rank1_rmsd": ev_r.get("rank1_rmsd")}


def _orch(run_dir: Path, reg):
    from kdss_audit.orchestrator import WorkflowOrchestrator
    return WorkflowOrchestrator(run_dir=run_dir, contracts_dir=FROZEN / "contracts", config_path=run_dir / "config.yaml", registry=reg)


def run_deployed_stage(orch, stage: str, inputs: dict) -> dict:
    try:
        sr = orch.run_stage(stage_name=stage, contract_file=CONTRACT_FILE[stage], tool_fn=None, tool_outputs=inputs)
    except Exception as exc:                                     # noqa: BLE001
        return {"state": "ERROR", "error": f"{type(exc).__name__}: {exc}", "checks_failed": []}
    state = getattr(sr.contract_result.state, "value", str(sr.contract_result.state))
    return {"state": state, "checks_failed": list(sr.contract_result.checks_failed), "halted": bool(orch._halted)}


def run_deployed_g3(ws: Path, inputs: dict) -> dict:
    """Deployed path: the orchestrator is created on the registry that holds the registrations made in this workspace."""
    from kdss_audit.artifact_registry import ArtifactRegistry
    run_dir = ws / "_registry"
    reg = ArtifactRegistry(run_dir / "artifact_registry.jsonl")
    return run_deployed_stage(_orch(run_dir, reg), "G3", inputs)


# ------------------------------------------------------------------------------------------------ G5
PER = RUNS / "chembl245_ifp_rescoring_2026-08-15" / "ifp_rescoring_per_compound.csv"
MULTI = RUNS / "chembl245_ifp_rescoring_2026-08-15" / "multi_reference_ifp_per_compound.csv"
# (name, source, column, sign): higher "score" is better after the sign, as in ifp_rescore_benchmark.py
G5_VARIANTS = [("vina_only", PER, "vina_score", -1), ("vina_at_best_ifp_pose", PER, "vina_at_best_ifp_pose", -1),
               ("max_ifp_tanimoto", PER, "max_ifp_tanimoto", 1), ("asn253_filtered", PER, "asn253_filtered_score", 1),
               ("combined_zscore", PER, "combined_zscore", 1), ("z_ifp", PER, "z_ifp", 1), ("z_vina", PER, "z_vina", -1),
               ("ifp_zm_only", MULTI, "ifp_zm_only", 1), ("ifp_9xt_only", MULTI, "ifp_9xt_only", 1),
               ("ifp_union_ref", MULTI, "ifp_union_ref", 1), ("ifp_max_of_two", MULTI, "ifp_max_of_two", 1)]
MODULE_HASHES = None


def kdss_hashes() -> dict:
    global MODULE_HASHES
    if MODULE_HASHES is None:
        MODULE_HASHES = {p.name: hashlib.sha256(p.read_bytes()).hexdigest().upper() for p in (FROZEN / "kdss_audit").glob("*.py")}
    return MODULE_HASHES


def bench_frame() -> pd.DataFrame:
    return pd.read_csv(PER)[["compound_id", "label"]]


def metrics(scores: np.ndarray, labels: np.ndarray) -> dict:
    return ifp_module().compute_metrics([(float(s), int(l)) for s, l in zip(scores, labels)])


def hash_ids(ids) -> str:
    return hashlib.sha256("\n".join(sorted(map(str, ids))).encode()).hexdigest().upper()


def build_g5(ws: Path, variant: str, leak_ids: list[str] | None = None, leak_from_benchmark: bool = False, label_swap: dict | None = None):
    """A G5 workspace: the benchmark file, the ranking scores, the registry (declaration, accesses) and the inputs of G5a and G5b.
    leak_ids: public active IDs exposed to the ranking step (F11); the leaked actives are ranked first. label_swap: F12."""
    ws.mkdir(parents=True, exist_ok=True)
    name, src, col, sign = next(v for v in G5_VARIANTS if v[0] == variant)
    df = pd.read_csv(src)
    bench = bench_frame().copy()
    bench_path, ranked_path = ws / "chembl_benchmark.csv", ws / "ranked_scores.csv"
    original_ids = bench["compound_id"].tolist()
    id_hash = hash_ids(original_ids)                       # registered at the start of the run, before any fault
    scores = (sign * df.set_index("compound_id").loc[bench["compound_id"], col].to_numpy(float)).copy()
    if label_swap:
        bench["label"] = label_swap["labels"]
    if leak_ids:
        pos = bench["compound_id"].isin(leak_ids).to_numpy()
        scores[pos] = scores.max() + 1.0                   # the leaked actives end up on top of the ranking
    bench.to_csv(bench_path, index=False)
    pd.DataFrame({"compound_id": bench["compound_id"], "score": scores}).to_csv(ranked_path, index=False)
    m = metrics(scores, bench["label"].to_numpy())
    reg = new_registry(ws / "_registry" / "artifact_registry.jsonl")
    reg.record_benchmark_role_declaration("Option_A_development", "Adisorn Sirichotejirakul", "the set was used during legacy method development")
    log_read(reg, "run_ranking", ranked_path, "ranking of the benchmark compounds by score")
    if leak_ids:
        leak_path = ws / "ranking_context_known_actives.txt"
        leak_path.write_text("\n".join(leak_ids))
        log_read(reg, "run_ranking", bench_path if leak_from_benchmark else leak_path,
                 "ranking step read the known actives" if not leak_from_benchmark else "ranking step read the benchmark labels")
    time.sleep(0.02)
    gate_start = time.time()
    time.sleep(0.02)
    log_read(reg, "G5b", bench_path, "G5b reads the labels at the gate")        # after the gate start: allowed
    g5a = {"chembl_benchmark_csv": str(bench_path), "benchmark_role_declaration": True, "g5a_evaluator_source_hash":
           hashlib.sha256((FROZEN / "evaluators" / "G5a_evaluator.py").read_bytes()).hexdigest().upper(),
           "kdss_audit_module_hashes": kdss_hashes(), "benchmark_id_list_hash_registered": id_hash,
           "bootstrap_known_answer_test": {"pairs_resampled_together": True}, "bootstrap_n_resamples": 2000,
           "ef1pct_count": int(bench["label"].to_numpy()[np.argsort(-scores, kind="stable")[:max(1, int(np.ceil(0.01 * len(bench))))]].sum()),
           "_gate_start_time": gate_start}
    g5b = {"auc_roc": m.get("auc_roc"), "ef_1pct": m.get("ef_1pct"), "bedroc20": m.get("bedroc20"), "registered_auc_min_threshold": 0.60,
           "benchmark_role": "Option_A_development"}
    return {"G5a": g5a, "G5b": g5b, "metrics": m}, reg


def run_direct_g5(g5: dict, reg) -> dict:
    import evaluators.G5a_evaluator as ea
    import evaluators.G5b_evaluator as eb
    out = {}
    ca = yaml.safe_load((FROZEN / "contracts" / CONTRACT_FILE["G5a"]).read_text(encoding="utf-8"))
    try:
        ra = ea.evaluate(contract=ca, inputs=g5["G5a"], registry=reg)
    except Exception as exc:                                     # noqa: BLE001
        return {"G5a": {"state": "ERROR", "error": f"{type(exc).__name__}: {exc}", "checks_failed": []}, "G5b": None}
    out["G5a"] = {"state": ra["state"], "checks_failed": list(ra["failed"]), "warned": list(ra["warned"])}
    if ra["state"] == "PASS":                              # G5b is conditional on G5a passing
        cb = yaml.safe_load((FROZEN / "contracts" / CONTRACT_FILE["G5b"]).read_text(encoding="utf-8"))
        rb = eb.evaluate(contract=cb, inputs=g5["G5b"], registry=reg)
        out["G5b"] = {"state": rb["state"], "checks_failed": list(rb["failed"]), "warned": list(rb["warned"])}
    else:
        out["G5b"] = None
    return out


def run_deployed_g5(ws: Path, g5: dict) -> dict:
    from kdss_audit.artifact_registry import ArtifactRegistry
    run_dir = ws / "_registry"
    reg = ArtifactRegistry(run_dir / "artifact_registry.jsonl")
    orch = _orch(run_dir, reg)
    out = {"G5a": run_deployed_stage(orch, "G5a", g5["G5a"])}
    out["G5b"] = run_deployed_stage(orch, "G5b", g5["G5b"]) if out["G5a"]["state"] == "PASS" else None
    return out


F11_GRID = [{"n_leaked_ids": n} for n in (5, 10, 15, 20, 25, 30, 40, 50, 75, 100)]
F12_GRID = [{"swap_mode": m} for m in ("full_swap", "partial_50pct", "partial_25pct")]


def swapped_labels(mode: str, seed: int = 7) -> np.ndarray:
    lab = bench_frame()["label"].to_numpy().copy()
    if mode == "full_swap":
        return 1 - lab
    frac = 0.5 if mode == "partial_50pct" else 0.25
    rng = np.random.default_rng(seed)
    act, ina = np.where(lab == 1)[0], np.where(lab == 0)[0]
    k = int(round(frac * len(act)))
    a, i = rng.choice(act, k, replace=False), rng.choice(ina, k, replace=False)
    lab[a], lab[i] = 0, 1
    return lab


def leaked_ids(n: int, seed: int = 7) -> list[str]:
    b = bench_frame()
    act = b.loc[b["label"] == 1, "compound_id"].tolist()
    return list(np.random.default_rng(seed).permutation(act)[:n])

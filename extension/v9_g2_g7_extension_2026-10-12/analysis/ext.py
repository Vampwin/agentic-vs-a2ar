"""ADR-AIIM-416: workspaces, injectors and the two evaluation paths for G2 and G7. Nothing in 10_new_system is modified; the frozen
evaluators and contracts are imported and read. The injectors implement the parameter grids of the stubs F06, F15 and F16 of
`10_new_system/fault_injection/` outside that tree.

G2 workspaces: the real prepared receptors of 4EIY of this project, the real reference ligand and docking box.
G7 workspaces: the real ligand-RMSD series of the molecular-dynamics runs of this project (GROMACS xvg files), their seeds read from
the run parameter files, and the real selected panel."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Agentic VS Workflow")
RUNS = ROOT / "20_new_runs"
FROZEN = ROOT / "10_new_system"
OUT = RUNS / "v9_g2_g7_extension_2026-10-12"
sys.path.insert(0, str(FROZEN))

# ------------------------------------------------------------------------------------------------ G2
G2_DIR = RUNS / "snp_screen_2026-08-04" / "G2_preparation"
G2_RECEPTORS = [   # (name, path, water cutoff in A): the ten distinct prepared 4EIY receptors of this project
    ("snp_screen", G2_DIR / "4EIY_prepared.pdbqt", 5.0),
    ("snp_screen_repro1", RUNS / "snp_screen_repro_run1_2026-08-15" / "G2_preparation" / "4EIY_prepared.pdbqt", 5.0),
    ("snp_screen_repro2", RUNS / "snp_screen_repro_run2_2026-08-15" / "G2_preparation" / "4EIY_prepared.pdbqt", 5.0),
    ("shell_cut3", RUNS / "routeA_water_shell_2026-08-09" / "4EIY_cut3_prepared.pdbqt", 3.0),
    ("shell_cut4", RUNS / "routeA_water_shell_2026-08-09" / "4EIY_cut4_prepared.pdbqt", 4.0),
    ("shell_cut5", RUNS / "routeA_water_shell_2026-08-09" / "4EIY_cut5_prepared.pdbqt", 5.0),
    ("clinical_reference", RUNS / "clinical_reference_2026-09-03" / "receptor" / "4EIY_prepared.pdbqt", 5.0),
    ("crossdock_zma", RUNS / "crossdock_zma_2026-09-03" / "receptors" / "4EIY_prepared.pdbqt", 5.0),
    ("crossdock_validation", RUNS / "crossdock_validation_2026-09-03" / "receptors" / "4EIY_prepared.pdbqt", 0.0),
    ("redock_9xt", RUNS / "redock_9xt_2026-09-05" / "receptors" / "4EIY_prepared.pdbqt", 5.0),
]
N_VALID = 223            # distinct valid structures of the 236-entry library (G1 ledger); all embedded in the compliant clean run


def _water_count(path: Path) -> int:
    seen = set()
    for ln in path.read_text(errors="replace").splitlines():
        if ln.startswith(("ATOM  ", "HETATM")) and ln[17:20].strip() in ("HOH", "WAT"):
            seen.add((ln[21:22], ln[22:26].strip()))
    return len(seen)


def build_g2(ws: Path, name: str, receptor: Path, cutoff: float) -> dict:
    """A compliant clean G2 workspace: the receptor with the retained-water decision that the file itself bears, the real reference
    ligand and box, the cofactor decisions and the embedding counts. Returns the evaluator inputs."""
    ws.mkdir(parents=True, exist_ok=True)
    rec = ws / "4EIY_prepared.pdbqt"
    shutil.copyfile(receptor, rec)
    shutil.copyfile(G2_DIR / "ZMA_reference_4EIY.sdf", ws / "ZMA_reference_4EIY.sdf")
    shutil.copyfile(G2_DIR / "docking_box.json", ws / "docking_box.json")
    n_w = _water_count(rec)
    return {
        "prepared_receptor_4eiy_pdbqt": str(rec),
        "reference_ligand_4eiy_sdf": str(ws / "ZMA_reference_4EIY.sdf"),
        "docking_box_json": str(ws / "docking_box.json"),
        "structural_water_decision": {"cutoff": cutoff, "reference_ligand": "ZM241385", "n_waters_retained": n_w},
        "cofactor_decisions": {"ZMA": "removed: reference ligand extracted and kept as a separate file",
                               "HOH": f"retained: {n_w} waters within {cutoff} A of ZM241385" if n_w else "removed: no waters retained"},
        "n_embedded_successfully": N_VALID, "count_valid_for_docking": N_VALID, "embedding_failures": {},
        "water_shell_candidates": [],
    }


def _his_residues(lines: list[str]) -> dict[str, list[int]]:
    """Residue key -> indices of its ATOM lines, for every HIS residue."""
    out: dict[str, list[int]] = {}
    for i, ln in enumerate(lines):
        if ln.startswith("ATOM  ") and ln[17:20] == "HIS":
            out.setdefault(ln[21:27], []).append(i)
    return out


def _atom(lines, idxs, name):
    for i in idxs:
        if lines[i][12:16].strip() == name:
            return i
    return None


def _xyz(ln):
    return np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])


def _fmt(ln, name=None, xyz=None, typ=None):
    s = ln
    if name is not None:
        s = s[:12] + f" {name:<3s}" + s[16:] if len(name) < 4 else s[:12] + name + s[16:]
    if xyz is not None:
        s = s[:30] + f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}" + s[54:]
    if typ is not None:
        s = s[:77] + f"{typ:<2s}" + s[79:] if len(s) >= 79 else s[:77] + typ
    return s


def flip_his_tautomer(lines: list[str], key: str, idxs: list[int]) -> bool:
    """Move the ring hydrogen of one histidine from ND1 to NE2 (or back), the protonation-state flip HID <-> HIE: the hydrogen record is
    renamed in place and placed 1.01 A from its new nitrogen away from the ring centre, and the two nitrogens swap AutoDock type (the
    protonated one is N, the other NA). Line count, serial numbers and residue keys are unchanged. Returns False when the residue
    holds neither or both of HD1 and HE2 (nothing to flip)."""
    hd1, he2 = _atom(lines, idxs, "HD1"), _atom(lines, idxs, "HE2")
    nd1, ne2 = _atom(lines, idxs, "ND1"), _atom(lines, idxs, "NE2")
    ring = [_atom(lines, idxs, n) for n in ("CG", "ND1", "CD2", "CE1", "NE2")]
    if (hd1 is None) == (he2 is None) or nd1 is None or ne2 is None or any(r is None for r in ring):
        return False
    centre = np.mean([_xyz(lines[r]) for r in ring], axis=0)

    def place(n_idx):
        n = _xyz(lines[n_idx])
        v = n - centre
        return n + 1.01 * v / np.linalg.norm(v)

    if hd1 is not None:                      # HID -> HIE
        lines[hd1] = _fmt(lines[hd1], name="HE2", xyz=place(ne2))
        lines[nd1] = _fmt(lines[nd1], typ="NA")
        lines[ne2] = _fmt(lines[ne2], typ="N")
    else:                                    # HIE -> HID
        lines[he2] = _fmt(lines[he2], name="HD1", xyz=place(nd1))
        lines[nd1] = _fmt(lines[nd1], typ="N")
        lines[ne2] = _fmt(lines[ne2], typ="NA")
    return True


F06_GRID = [{"n_his_flipped": n, "residue_selection": sel} for n in (1, 2, 3) for sel in ("first", "random", "all")][:10]


def inject_f06(ws: Path, params: dict, seed: int) -> dict:
    """F06 WRONG_PROTONATION on the receptor of the workspace: the retained-water decision, the cofactor decisions and every other
    input stay as the clean run recorded them."""
    rec = ws / "4EIY_prepared.pdbqt"
    lines = rec.read_text(errors="replace").split("\n")
    his = _his_residues(lines)
    keys = sorted(his, key=lambda k: int(k[1:5]))
    sel = params["residue_selection"]
    n = params["n_his_flipped"]
    rng = np.random.default_rng(seed)
    if sel == "first":
        chosen = keys[:n]
    elif sel == "random":
        chosen = list(rng.choice(keys, size=min(n, len(keys)), replace=False))
    else:                                    # "all": every histidine, whatever n is (the grid pairs it with 1, 2 and 3)
        chosen = keys
    flipped = [k for k in chosen if flip_his_tautomer(lines, k, his[k])]
    rec.write_text("\n".join(lines))
    return {"flipped_residues": [k.strip() for k in flipped], "n_flipped": len(flipped), "vacuous": len(flipped) == 0}


# ------------------------------------------------------------------------------------------------ G7
CAND = RUNS / "md_candidates_2026-08-22"
SYSTEMS = {}
for sid in ("NPLIB01185", "NPLIB01643", "NPLIB03056", "SNP00157", "SNP00163"):
    SYSTEMS[sid] = [(f"r{r}", CAND / sid / "04_run" / f"analysis_prod_r{r}" / "ligand_rmsd.xvg",
                     CAND / sid / "04_run" / f"prod_r{r}.mdp") for r in (1, 2, 3)]
CTRL_REP = RUNS / "md_control_replicates_2026-10-06"
SYSTEMS["ZM241385"] = [
    ("r1", RUNS / "md_control_anchor_2026-10-07" / "analysis_prod" / "ligand_rmsd.xvg", RUNS / "md_production_2026-08-21" / "04_run" / "nvt.mdp"),   # the velocity-generating step of the continued production
    ("r2", CTRL_REP / "analysis_prod_r2" / "ligand_rmsd.xvg", CTRL_REP / "prod_r2.mdp"),
    ("r3", CTRL_REP / "analysis_prod_r3" / "ligand_rmsd.xvg", CTRL_REP / "prod_r3.mdp"),
]
PLANNED_NS = 50.0
CLAIM = "MD results are binding-stability annotations only and do not constitute experimental validation."
SHORTLIST = ["SNP00157", "SNP00163", "SNP00154", "SNP00097", "SNP00052", "SNP00033", "SNP00142", "SNP00006"]


def read_xvg(p: Path) -> tuple[np.ndarray, np.ndarray]:
    t, y = [], []
    for ln in p.read_text().splitlines():
        if ln and ln[0] not in "#@":
            a = ln.split()
            t.append(float(a[0]))
            y.append(float(a[1]))
    return np.array(t), np.array(y)


def read_seed(mdp: Path):
    if not mdp.exists():
        return None
    for ln in mdp.read_text(errors="replace").splitlines():
        m = re.match(r"\s*gen[-_]seed\s*=\s*(-?\d+)", ln)
        if m:
            return int(m.group(1))
    return None


def pocket_residues(radius: float = 10.0) -> list[int]:
    """Receptor C-alpha residues within 10 A of the crystal ligand (the definition of the binding site used for the distance series)."""
    from rdkit import Chem
    mol = Chem.MolFromMolFile(str(G2_DIR / "ZMA_reference_4EIY.sdf"), removeHs=False)
    lig = mol.GetConformer().GetPositions()
    out = []
    for ln in (G2_DIR / "4EIY_prepared.pdbqt").read_text(errors="replace").splitlines():
        if ln.startswith("ATOM  ") and ln[12:16].strip() == "CA":
            if np.min(np.linalg.norm(lig - _xyz(ln), axis=1)) <= radius:
                out.append(int(ln[22:26]))
    return sorted(set(out))


def clean_g7_cases() -> list[tuple[str, str, list[str]]]:
    """Eleven clean systems: the six systems with all their replicates, and the five compounds with their first two replicates."""
    cases = []
    for sid, runs in SYSTEMS.items():
        cases.append((f"{sid}_all", sid, [r[0] for r in runs]))
    for sid in SYSTEMS:
        if sid != "ZM241385":
            cases.append((f"{sid}_r1r2", sid, ["r1", "r2"]))
    return cases


def build_g7(ws: Path, sid: str, reps: list[str]) -> dict:
    ws.mkdir(parents=True, exist_ok=True)
    runs = {r[0]: r for r in SYSTEMS[sid] if r[0] in reps}
    comp, eq, res, seeds = [], [], [], {}
    pocket = pocket_residues()
    for rep, xvg, mdp in runs.values():
        t, y = read_xvg(xvg)
        seed = read_seed(mdp)
        if seed is None:
            raise RuntimeError(f"no gen-seed read for {sid} {rep} ({mdp})")
        seeds.setdefault(sid, seed)
        comp.append({"compound_id": sid, "replicate": rep, "status": "completed", "prod_ns": round(float(t[-1]), 3), "seed": seed})
        w = t <= PLANNED_NS + 1e-6
        for ti, yi in zip(t[w], y[w]):
            eq.append({"compound_id": sid, "replicate": rep, "time_ns": round(float(ti), 4), "ligand_rmsd_nm": round(float(yi), 6)})
        res.append({"compound_id": sid, "replicate": rep, "rmsd_mean_nm": round(float(y[w].mean()), 5),
                    "rmsd_max_nm": round(float(y[w].max()), 5), "n_frames": int(w.sum()),
                    "binding_site_residues": " ".join(map(str, pocket))})
    for name, rows in (("md_completion_log.csv", comp), ("equilibration_check.csv", eq), ("md_results.csv", res)):
        with open(ws / name, "w", newline="", encoding="utf-8") as f:
            wr = csv.DictWriter(f, fieldnames=list(rows[0]))
            wr.writeheader()
            wr.writerows(rows)
    with open(ws / "selected_candidates.csv", "w", newline="", encoding="utf-8") as f:
        wr = csv.writer(f)
        wr.writerow(["compound_id", "shortlist_position"])
        for i, c in enumerate(SHORTLIST, 1):
            wr.writerow([c, i])
    return {"md_decision_recorded": "MD_run", "md_completion_log": str(ws / "md_completion_log.csv"),
            "equilibration_check_csv": str(ws / "equilibration_check.csv"), "md_results_csv": str(ws / "md_results.csv"),
            "planned_production_ns": PLANNED_NS, "required_claim_boundary": CLAIM,
            "selected_candidates_csv": str(ws / "selected_candidates.csv"), "selection_reranked_by_md": False,
            "md_seeds_by_compound": seeds}


def _rewrite(path: Path, fn):
    import pandas as pd
    df = pd.read_csv(path)
    df = fn(df)
    df.to_csv(path, index=False)
    return df


F15_GRID = [{"fraction_kept": f} for f in (0.05, 0.10, 0.15, 0.20, 0.25, 0.33, 0.40, 0.50, 0.67, 0.80)]
F16_GRID = [{"residue_offset": o} for o in (-50, -30, -20, -10, -5, 5, 10, 20, 30, 50)]


def inject_f15(ws: Path, params: dict, log_updated: bool) -> dict:
    """F15 TRUNCATED_TRAJ: the per-frame series is cut to `fraction_kept` of the planned production window before analysis. In F15
    the completion log stays as the run wrote it (the truncation is silent, like the other families); in F15L the log records the
    shorter production length."""
    f = params["fraction_kept"]
    cut = f * PLANNED_NS
    d = _rewrite(ws / "equilibration_check.csv", lambda df: df[df["time_ns"] <= cut + 1e-9])
    if log_updated:
        _rewrite(ws / "md_completion_log.csv", lambda df: df.assign(prod_ns=np.minimum(df["prod_ns"], round(cut, 3))))
    return {"fraction_kept": f, "frames_left": int(len(d)), "vacuous": False}


def inject_f16(ws: Path, params: dict) -> dict:
    """F16 RESIDUE_INDEX_ERROR: the residue indices that define the binding site used for the RMSD are shifted by a fixed offset."""
    o = params["residue_offset"]

    def shift(df):
        df["binding_site_residues"] = df["binding_site_residues"].map(lambda s: " ".join(str(int(x) + o) for x in str(s).split()))
        return df
    _rewrite(ws / "md_results.csv", shift)
    return {"residue_offset": o, "vacuous": False}


# ------------------------------------------------------------------------------------------------ evaluation paths
CONTRACT_FILE = {"G2": "G2_preparation.yaml", "G7": "G7_md.yaml"}


def run_direct(gate: str, inputs: dict) -> dict:
    """Arm B: the frozen evaluator called directly with the frozen YAML contract (as arm B of ADR-AIIM-236)."""
    import importlib
    import yaml
    contract = yaml.safe_load((FROZEN / "contracts" / CONTRACT_FILE[gate]).read_text(encoding="utf-8"))
    mod = importlib.import_module(f"evaluators.{gate}_evaluator")
    try:
        r = mod.evaluate(contract=contract, inputs=inputs, registry=None)
    except Exception as exc:                                     # noqa: BLE001
        return {"state": "ERROR", "error": f"{type(exc).__name__}: {exc}", "checks_failed": [], "checks_passed": []}
    return {"state": r.get("state"), "checks_failed": list(r.get("failed", [])), "checks_passed": list(r.get("passed", [])),
            "warned": list(r.get("warned", []))}


def run_deployed(gate: str, inputs: dict, ws: Path) -> dict:
    """Arm C: the deployed path of the frozen tree: orchestrator, contract engine and report lock, as arm C of ADR-AIIM-236."""
    from kdss_audit.artifact_registry import ArtifactRegistry
    from kdss_audit.orchestrator import WorkflowOrchestrator
    run_dir = ws / "_deployed"
    run_dir.mkdir(parents=True, exist_ok=True)
    registry = ArtifactRegistry(run_dir / "artifact_registry.jsonl")
    orch = WorkflowOrchestrator(run_dir=run_dir, contracts_dir=FROZEN / "contracts", config_path=run_dir / "config.yaml", registry=registry)
    try:
        sr = orch.run_stage(stage_name=gate, contract_file=CONTRACT_FILE[gate], tool_fn=None, tool_outputs=inputs)
    except Exception as exc:                                     # noqa: BLE001
        return {"state": "ERROR", "error": f"{type(exc).__name__}: {exc}", "halted": None, "reportable": None, "contained": None,
                "checks_failed": []}
    state = getattr(sr.contract_result.state, "value", str(sr.contract_result.state))
    lock = orch.build_report_lock(run_dir / "report.lock")
    reportable = (all(lock.can_report(a) for a in sr.artifacts_produced) if sr.artifacts_produced else None)
    halted = bool(orch._halted)
    return {"state": state, "checks_failed": list(sr.contract_result.checks_failed), "halted": halted, "reportable": reportable,
            "contained": bool(halted or (reportable is False))}

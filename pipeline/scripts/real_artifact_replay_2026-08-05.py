"""
real_artifact_replay_2026-08-05.py

Injects real, historically-motivated faults directly into copies of the
actual A2AR/SNP-Lib production artifacts and runs the real, current
evaluators against them -- not the synthetic fault_injection/ toy harness.
Answers a narrower but stronger question than the fault-injection pilot:
does the audit layer catch faults on the files this project actually
produced, not just on templated development data.

Four replays are attempted, chosen because a real artifact exists to fault:

  G0 -- byte-corrupt one character of the real source library CSV
  G1 -- strip stereochemistry from one real compound's SMILES in the real
        236-row standardization ledger, leaving its declared stereo columns
        unchanged (the mismatch pattern stereochemistry_preserved exists for)
  G3 -- promote a genuine decoy pose -- already present among the 80 real
        pooled poses from the actual production redocking -- to rank 1, in
        place of the real near-native pose. Same class of defect as the
        historical Kabsch/MCS-submatch RMSD bug (a genuinely wrong pose
        scored as acceptable), reproduced with real geometry rather than a
        synthetic rotation.
  G4 -- negate the real Vina scores from the 13-compound on-pocket
        calibration study (g3_validation_study_2026-08-04) -- the actual
        historical sign-convention bug class (fault family F09 in the
        synthetic harness, replayed here on real docking output)

G6 and G8 are NOT attempted: G4-G8 were never formally run for this
application (ADR-AIIM-037, negative result -- no ranking was released), so
there is no real selection or report artifact to fault. Fabricating one to
complete the set would be exactly the kind of manufactured-looking-real data
this project's own audit philosophy exists to catch.

Each replay runs a clean control (the real, unmodified artifact) alongside
the faulted case, so a catch is attributable to the fault, not to a check
that always fires regardless of input.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
from rdkit import Chem

from evaluators import G0_evaluator, G1_evaluator, G3_evaluator, G4_evaluator

FROZEN = Path(r"D:/Agentic VS Workflow/02_frozen_inputs")
PROD_RUN = Path(r"D:/Agentic VS Workflow/20_new_runs/snp_screen_2026-08-04")
CALIB_DIR = Path(r"D:/Agentic VS Workflow/20_new_runs/g3_validation_study_2026-08-04")
OUT_DIR = Path(r"D:/Agentic VS Workflow/20_new_runs/real_artifact_replay_2026-08-05")
CONTRACTS = Path(r"D:/Agentic VS Workflow/10_new_system/contracts")

import yaml


def load_contract(name: str) -> dict:
    return yaml.safe_load((CONTRACTS / name).read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest().upper()


def log(msg: str) -> None:
    print(f"[replay] {msg}")


results: list[dict] = []


def record(gate: str, check_id: str, description: str,
           clean_outcome: dict, fault_outcome: dict, caught: bool) -> None:
    results.append({
        "gate": gate, "check_id": check_id, "description": description,
        "clean_outcome": clean_outcome, "fault_outcome": fault_outcome,
        "caught": caught,
    })
    status = "CAUGHT" if caught else "MISSED"
    log(f"{gate} [{check_id}] {status} -- clean={clean_outcome.get('state')} "
        f"fault={fault_outcome.get('state')}")


def state_of(report: dict) -> str:
    return report.get("state", "?")


# ─────────────────────────────────────────────────────────────────────────────
# G0 -- byte-corrupt the real source library CSV
# ─────────────────────────────────────────────────────────────────────────────

def replay_g0() -> None:
    ws = OUT_DIR / "G0"
    ws.mkdir(parents=True, exist_ok=True)
    real_csv = FROZEN / "snp_library" / "01_library_validated.csv"
    real_sha = sha256(real_csv)

    clean_copy = ws / "01_library_validated_clean.csv"
    shutil.copyfile(real_csv, clean_copy)

    raw = real_csv.read_bytes()
    # Flip one character inside the file body (not the header line) --
    # the same class of corruption F02 (hash_mismatch) models in the
    # synthetic harness, applied to the real 236-compound source file.
    header_end = raw.index(b"\n") + 1
    target = header_end + 200  # well inside the data region
    corrupted = bytearray(raw)
    corrupted[target] = (corrupted[target] + 1) % 256
    fault_copy = ws / "01_library_validated_faulted.csv"
    fault_copy.write_bytes(bytes(corrupted))

    def run(csv_path: Path) -> dict:
        inputs = {
            "snp_library_csv": str(csv_path),
            "snp_library_sdf": str(FROZEN / "snp_library" / "01_library_validated.sdf"),
            "snp_library_raw_xlsx": str(FROZEN / "snp_library" /
                "Sritularak Natural Products Library (SNP-Lib) - Prof. Boonchoo Sritularak_update.xlsx"),
            "receptor_pdb_4eiy": str(FROZEN / "a2ar_structures" / "4EIY.pdb"),
            "receptor_pdb_3eml": str(FROZEN / "a2ar_structures" / "3EML.pdb"),
            "chembl_benchmark_csv": str(FROZEN / "public_benchmarks" / "chembl_a2ar" /
                "chembl_a2ar_validation_set.csv"),
            "registered_input_hashes": {"hash_snp_library_csv": real_sha},
        }
        contract = load_contract("G0_input.yaml")
        return G0_evaluator.evaluate(contract, inputs, registry=None)

    clean_report = run(clean_copy)
    fault_report = run(fault_copy)
    cid = "hash_snp_library_csv"
    caught = cid in clean_report["passed"] and cid in fault_report["failed"]
    record(
        "G0", cid,
        "One byte flipped inside the real 236-row source library CSV; real "
        "SHA-256 anchor supplied for both runs.",
        {"state": clean_report["state"], cid: clean_report["evidence"].get(cid)},
        {"state": fault_report["state"], cid: fault_report["evidence"].get(cid)},
        caught,
    )


# ─────────────────────────────────────────────────────────────────────────────
# G1 -- strip stereochemistry from one real compound in the real ledger
# ─────────────────────────────────────────────────────────────────────────────

def replay_g1() -> None:
    ws = OUT_DIR / "G1"
    ws.mkdir(parents=True, exist_ok=True)
    real_ledger = PROD_RUN / "G1_standardization" / "snp_source_ledger.csv"

    df = pd.read_csv(real_ledger)
    candidates = df[(df["has_stereocenters"] == True) & (df["stereo_defined"] == True) &
                     (df["n_stereocenters"] > 0)]
    if candidates.empty:
        log("G1: no real compound with declared+defined stereocentres found; skipping")
        return
    row_idx = candidates.index[0]
    snp_id = df.loc[row_idx, "SNP_ID"]
    original_smiles = df.loc[row_idx, "smiles_canonical"]

    mol = Chem.MolFromSmiles(original_smiles)
    Chem.RemoveStereochemistry(mol)
    stripped_smiles = Chem.MolToSmiles(mol)
    log(f"G1: {snp_id} stereo stripped -- '{original_smiles}' -> '{stripped_smiles}' "
        f"(declared n_stereocenters={df.loc[row_idx, 'n_stereocenters']} left unchanged)")

    clean_copy = ws / "snp_source_ledger_clean.csv"
    fault_copy = ws / "snp_source_ledger_faulted.csv"
    df.to_csv(clean_copy, index=False)
    df_fault = df.copy()
    df_fault.loc[row_idx, "smiles_canonical"] = stripped_smiles
    df_fault.to_csv(fault_copy, index=False)

    def run(csv_path: Path) -> dict:
        inputs = {"snp_source_ledger_csv": str(csv_path), "standardized_smiles_sdf": ""}
        contract = load_contract("G1_standardization.yaml")
        return G1_evaluator.evaluate(contract, inputs, registry=None)

    clean_report = run(clean_copy)
    fault_report = run(fault_copy)
    cid = "stereochemistry_preserved"
    caught = cid in clean_report["passed"] and cid in fault_report["failed"]
    record(
        "G1", cid,
        f"Stereo markers stripped from {snp_id}'s SMILES in the real 236-row "
        "ledger; declared stereocentre metadata left saying stereo is present.",
        {"state": clean_report["state"], cid: clean_report["evidence"].get(cid)},
        {"state": fault_report["state"], cid: fault_report["evidence"].get(cid)},
        caught,
    )


# ─────────────────────────────────────────────────────────────────────────────
# G3 -- promote a genuine decoy (from the real 80-pose pool) to rank 1
# ─────────────────────────────────────────────────────────────────────────────

def replay_g3() -> None:
    ws = OUT_DIR / "G3"
    ws.mkdir(parents=True, exist_ok=True)
    g3_dir = PROD_RUN / "G3_redocking"
    real_poses_sdf = g3_dir / "ZMA_redocked_4EIY.sdf"
    reference_sdf = Path(r"D:/Agentic VS Workflow/04_legacy_evidence/redocking/ZMA_crystal_4EIY.sdf")
    if not reference_sdf.is_file():
        log("G3: could not locate ZMA_crystal_4EIY.sdf; skipping")
        return

    supplier = Chem.SDMolSupplier(str(real_poses_sdf), removeHs=False)
    poses = [m for m in supplier if m is not None]
    if len(poses) < 2:
        log("G3: fewer than 2 real poses parsed; skipping")
        return

    ref_supplier = Chem.SDMolSupplier(str(reference_sdf), removeHs=False)
    ref_mol = next((m for m in ref_supplier if m is not None), None)
    if ref_mol is None:
        log("G3: could not parse reference ligand; skipping")
        return
    ref_conf = ref_mol.GetConformer()
    ref_centroid = [sum(ref_conf.GetAtomPosition(i)[k] for i in range(ref_mol.GetNumAtoms())) / ref_mol.GetNumAtoms()
                    for k in range(3)]

    def centroid_dist(mol) -> float:
        conf = mol.GetConformer()
        n = mol.GetNumAtoms()
        c = [sum(conf.GetAtomPosition(i)[k] for i in range(n)) / n for k in range(3)]
        return sum((c[k] - ref_centroid[k]) ** 2 for k in range(3)) ** 0.5

    dists = [(i, centroid_dist(m)) for i, m in enumerate(poses)]
    dists_sorted = sorted(dists, key=lambda t: -t[1])
    decoy_idx, decoy_dist = dists_sorted[0]
    rank1_dist = dists[0][1]
    log(f"G3: real rank-1 pose centroid distance={rank1_dist:.2f} A; "
        f"real decoy pose (index {decoy_idx} of {len(poses)}) centroid distance={decoy_dist:.2f} A")

    clean_copy = ws / "ZMA_redocked_4EIY_clean.sdf"
    shutil.copyfile(real_poses_sdf, clean_copy)

    reordered = [poses[decoy_idx]] + [m for i, m in enumerate(poses) if i != decoy_idx]
    fault_copy = ws / "ZMA_redocked_4EIY_faulted.sdf"
    writer = Chem.SDWriter(str(fault_copy))
    for m in reordered:
        writer.write(m)
    writer.close()

    def run(poses_path: Path) -> dict:
        inputs = {
            "reference_ligand_4eiy_sdf": str(reference_sdf),
            "redocked_poses_4eiy": str(poses_path),
            "atom_map_4eiy_json": str(g3_dir / "atom_map_4eiy.json"),
            "docking_input_conformer_sdf": str(g3_dir / "ZM241385_conformers.sdf"),
            "frame_check_4eiy": True,
        }
        contract = load_contract("G3_pose.yaml")
        return G3_evaluator.evaluate(contract, inputs, registry=None)

    clean_report = run(clean_copy)
    fault_report = run(fault_copy)
    cid = "rmsd_rank1_4eiy"
    caught = cid in clean_report["passed"] and cid in fault_report["failed"]
    record(
        "G3", cid,
        f"Rank-1 pose swapped for a genuine decoy already present among the real "
        f"80 pooled poses (centroid distance {decoy_dist:.2f} A vs real rank-1's "
        f"{rank1_dist:.2f} A from the reference) -- real geometry, reordered, no "
        "synthetic coordinates.",
        {"state": clean_report["state"], cid: clean_report["evidence"].get(cid)},
        {"state": fault_report["state"], cid: fault_report["evidence"].get(cid)},
        caught,
    )


# ─────────────────────────────────────────────────────────────────────────────
# G4 -- negate the real Vina scores from the calibration study
# ─────────────────────────────────────────────────────────────────────────────

def replay_g4() -> None:
    ws = OUT_DIR / "G4"
    ws.mkdir(parents=True, exist_ok=True)
    frames = []
    for fname in ("_calib.csv", "_calib_strat.csv"):
        p = CALIB_DIR / fname
        if p.is_file():
            frames.append(pd.read_csv(p))
    if not frames:
        log("G4: no real calibration Vina-score CSVs found; skipping")
        return
    df = pd.concat(frames, ignore_index=True)
    df = df[df["status"] == "ok"][["compound_id", "vina_score"]].dropna()
    if df.empty:
        log("G4: no valid real Vina scores found; skipping")
        return
    log(f"G4: {len(df)} real Vina scores from the on-pocket calibration study "
        f"(range {df['vina_score'].min():.2f} to {df['vina_score'].max():.2f} kcal/mol)")

    clean_copy = ws / "vina_scores_clean.csv"
    fault_copy = ws / "vina_scores_faulted.csv"
    df.to_csv(clean_copy, index=False)
    df_fault = df.copy()
    df_fault["vina_score"] = -df_fault["vina_score"]
    df_fault.to_csv(fault_copy, index=False)

    def run(csv_path: Path) -> dict:
        inputs = {"vina_score_csv": str(csv_path), "count_for_docking": len(df)}
        contract = load_contract("G4_scoring.yaml")
        return G4_evaluator.evaluate(contract, inputs, registry=None)

    clean_report = run(clean_copy)
    fault_report = run(fault_copy)
    cid = "score_direction"
    caught = cid in clean_report["passed"] and cid in fault_report["failed"]
    record(
        "G4", cid,
        f"Sign of all {len(df)} real Vina binding-energy scores from the "
        "13-compound on-pocket calibration study inverted (the historical "
        "sign-convention bug class).",
        {"state": clean_report["state"], cid: clean_report["evidence"].get(cid)},
        {"state": fault_report["state"], cid: fault_report["evidence"].get(cid)},
        caught,
    )


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    replay_g0()
    replay_g1()
    replay_g3()
    replay_g4()

    n_caught = sum(1 for r in results if r["caught"])
    summary = {
        "run_date": "2026-08-05",
        "n_attempted": len(results),
        "n_caught": n_caught,
        "detection_rate": n_caught / len(results) if results else None,
        "not_attempted": {
            "G6": "no real selection artifact exists (G4-G8 not run for A2AR, ADR-AIIM-037)",
            "G8": "no real report artifact exists (G4-G8 not run for A2AR, ADR-AIIM-037)",
        },
        "results": results,
    }
    out_path = OUT_DIR / "real_artifact_replay_results.json"
    out_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    log(f"\n{n_caught}/{len(results)} real-artifact faults caught")
    log(f"wrote {out_path}")


if __name__ == "__main__":
    main()

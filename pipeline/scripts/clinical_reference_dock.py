"""clinical_reference_dock.py -- what does the score say about compounds we know?

Criteria are module constants transcribed from ADR-AIIM-163 **before any docking
was run**, the convention R2, R3, R4, ADR-AIIM-161 and ADR-AIIM-162 followed.

§3.4 states the docking score's discriminative power as an AUC over 245
compounds. This restates the same limitation on a scale a reader can picture:
ten A2A ligands with measured affinities spanning pChEMBL 6.82 to 9.17 -- one
approved drug, four that reached clinical trials, three reference antagonists,
two agonists -- docked under the protocol that produced every reported score,
beside the three candidates.

**No new input is registered.** Every SMILES and every affinity is a row of
`chembl_a2ar_validation_set.csv` and `chembl_a2a_raw_v2.csv`, both hash-anchored
at G0 already.

**The three rules from crossdock_zma.py are enforced here too.** A run directory
is written once; every external input is copied in and hashed before anything is
computed; and the control gates the result -- ZM241385 is in the set, and if it
does not redock within 1.0 A nothing else from this run is reported.

Usage:
    python clinical_reference_dock.py --out <a new directory>
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import shlex
import shutil
import subprocess
from pathlib import Path

import numpy as np

# --- pre-registered, ADR-AIIM-163 -------------------------------------------
CONTROL = "ZM241385"
CONTROL_MAX_RMSD_A = 1.0
HBOND_CUTOFF_A = 3.5          # matches ADR-AIIM-161 and the §2.5 IFP definition
WATER_CUTOFF_A = 5.0
PROTONATION_PH = 7.4
N_CONFORMERS, EXHAUSTIVENESS, N_POSES, SEED = 8, 32, 10, 20260729
RECEPTOR_PDB = "4EIY"
LIGAND_RESN = "ZMA"
# name -> ChEMBL id. Ten ligands with measured affinity, already frozen.
PANEL = {
    "istradefylline": "CHEMBL431770",   # approved (Nourianz)
    "preladenant": "CHEMBL240624",      # phase III
    "tozadenant": "CHEMBL2105747",      # phase III
    "vipadenant": "CHEMBL447664",       # phase II
    "imaradenant": "CHEMBL4594442",     # phase II
    "SCH-58261": "CHEMBL17127",         # reference antagonist
    "CGS-15943": "CHEMBL16687",         # reference antagonist
    "ZM241385": "CHEMBL113142",         # the co-crystal ligand, and the control
    "NECA": "CHEMBL464859",             # agonist
    "adenosine": "CHEMBL477",           # endogenous agonist
}
CANDIDATES = ["NPLIB01185", "NPLIB01643", "NPLIB03056"]


def _root() -> Path:
    for c in (Path(r"D:\Agentic VS Workflow"), Path("/mnt/d/Agentic VS Workflow")):
        if (c / "00_project_control/08_DECISION_LOG.md").is_file():
            return c
    raise SystemExit("project root not found")


ROOT = _root()
SCRIPTS = ROOT / "10_new_system/scripts"
STRUCT = ROOT / "02_frozen_inputs/a2ar_structures"
BENCH = ROOT / "02_frozen_inputs/public_benchmarks/chembl_a2ar"
VALIDATION_SET = BENCH / "chembl_a2ar_validation_set.csv"
RAW_SET = BENCH / "chembl_a2a_raw_v2.csv"
BOX_SRC = ROOT / "20_new_runs/snp_screen_2026-08-04/G2_preparation/docking_box.json"
TEMPLATE_SRC = ROOT / "20_new_runs/snp_screen_2026-08-04/G3_redocking/reference_template.json"
NP_LEDGER = ROOT / "20_new_runs/np_library_funnel_2026-08-15/e8_np_library_prioritization.csv"
PDBFIXER_ENV = "mdenv"
WSL_DISTRO, CONDA_ENV = "Ubuntu-22.04", "diffdock"


def _wsl_path(p) -> str:
    s = str(p).replace(chr(92), "/")
    return f"/mnt/{s[0].lower()}{s[2:]}" if s[1:3] == ":/" else s


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


CD = _load("crossdock_zma", SCRIPTS / "crossdock_zma.py")
RSS = _load("run_snp_screen", SCRIPTS / "run_snp_screen.py")


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def wsl(cmd, timeout=14400):
    full = (f"source ~/miniconda3/etc/profile.d/conda.sh && conda activate {CONDA_ENV} "
            f"&& cd '/mnt/d/Agentic VS Workflow' && {cmd}")
    return subprocess.run(["wsl", "-d", WSL_DISTRO, "--", "bash", "-lc", full],
                          capture_output=True, text=True, timeout=timeout)


def asn253_contact(sdf: Path, receptor_pdb: Path) -> dict:
    """Min distance from any ligand N/O to the Asn253 sidechain, rank-1 pose.

    The asparagine is located the way ADR-AIIM-161 located it -- by its position
    in the TM6 run TRP/LEU/HIS/ASN at spacings 3-1-3 -- not by residue number,
    which is not unique in these files.
    """
    from rdkit import Chem
    atoms = CD.MR.read_pdb(receptor_pdb)
    inst, prev = [], None
    for a in atoms:
        key = (a["resi"], a["resn"])
        if key != prev:
            inst.append({"resn": a["resn"], "atoms": []})
            prev = key
        inst[-1]["atoms"].append(a)
    pos = None
    for i in range(7, len(inst)):
        if (inst[i]["resn"] == "ASN" and inst[i - 7]["resn"] == "TRP"
                and inst[i - 4]["resn"] == "LEU"
                and inst[i - 3]["resn"] in ("HIS", "HIE", "HID", "HIP")):
            pos = i
            break
    if pos is None:
        return {"available": False, "reason": "TM6 spacing 3-1-3 not found"}
    sc = np.array([[a["x"], a["y"], a["z"]] for a in inst[pos]["atoms"]
                   if a["name"] in ("OD1", "ND2")])
    ms = [m for m in Chem.SDMolSupplier(str(sdf), removeHs=True) if m]
    if not ms or sc.size == 0:
        return {"available": False, "reason": "no pose or no sidechain atoms"}
    conf = ms[0].GetConformer()
    lig = np.array([list(conf.GetAtomPosition(i))
                    for i, at in enumerate(ms[0].GetAtoms())
                    if at.GetSymbol() in ("N", "O")])
    if lig.size == 0:
        return {"available": True, "min_distance_A": None,
                "contact": False, "reason": "ligand has no N or O"}
    d = float(np.linalg.norm(lig[:, None, :] - sc[None, :, :], axis=2).min())
    return {"available": True, "min_distance_A": round(d, 2),
            "contact": d <= HBOND_CUTOFF_A}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out = Path(args.out).resolve()
    if out.exists() and any(out.iterdir()):
        raise SystemExit(f"REFUSED: {out} exists and is not empty. "
                         "A run directory is written once (ADR-AIIM-092).")
    (out / "inputs").mkdir(parents=True, exist_ok=True)

    import pandas as pd
    manifest = {"declared_inputs": {}, "authority": "ADR-AIIM-163",
                "recipe": {"receptor": RECEPTOR_PDB,
                           "preparation": "chain A + structural waters within "
                                          f"{WATER_CUTOFF_A} A -> pdbfixer "
                                          f"(pH {PROTONATION_PH}, keepWater) -> "
                                          "obabel -xr",
                           "docking": f"{N_CONFORMERS} conformers, "
                                      f"exhaustiveness {EXHAUSTIVENESS}, "
                                      f"{N_POSES} poses, seed {SEED}"}}
    for src in (VALIDATION_SET, RAW_SET, BOX_SRC, TEMPLATE_SRC,
                STRUCT / f"{RECEPTOR_PDB}.pdb", NP_LEDGER):
        shutil.copy2(src, out / "inputs" / src.name)
        manifest["declared_inputs"][src.name] = {
            "from": src.relative_to(ROOT).as_posix(), "sha256": sha256(src)}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2),
                                       encoding="utf-8")

    # --- panel, from frozen rows only ---------------------------------------
    v = pd.read_csv(VALIDATION_SET).drop_duplicates("chembl_id").set_index("chembl_id")
    med = pd.read_csv(RAW_SET).groupby("molecule_chembl_id")["pchembl_value"].median()
    npl = pd.read_csv(NP_LEDGER).set_index("compound_id")
    rows = []
    for name, cid in PANEL.items():
        rows.append({"compound_id": name, "chembl_id": cid,
                     "SMILES": v.loc[cid, "SMILES"],
                     "pchembl_median": round(float(med.get(cid, float("nan"))), 2),
                     "class": "reference ligand"})
    for cid in CANDIDATES:
        rows.append({"compound_id": cid, "chembl_id": "", "SMILES": npl.loc[cid, "SMILES"],
                     "pchembl_median": None, "class": "candidate, never assayed"})
    panel = pd.DataFrame(rows)
    ledger = out / "panel_ledger.csv"
    panel.to_csv(ledger, index=False, lineterminator="\n")
    print(f"panel: {len(panel)} compounds "
          f"({len(PANEL)} with measured affinity, {len(CANDIDATES)} candidates)")

    # --- receptor, the validated recipe -------------------------------------
    rec_dir = out / "receptor"
    rec_dir.mkdir(exist_ok=True)
    src = STRUCT / f"{RECEPTOR_PDB}.pdb"
    filt = rec_dir / f"{RECEPTOR_PDB}_protein_plus_waters.pdb"
    wat = RSS._write_receptor_with_waters(src, filt, chain="A", cutoff=WATER_CUTOFF_A)
    fixed = rec_dir / f"{RECEPTOR_PDB}_fixed.pdb"
    info = CD.pdbfix(filt, fixed)
    pdbqt = rec_dir / f"{RECEPTOR_PDB}_prepared.pdbqt"
    wsl(f'obabel -ipdb "{fixed.relative_to(ROOT).as_posix()}" -opdbqt '
        f'-O "{pdbqt.relative_to(ROOT).as_posix()}" -xr')
    print(f"receptor: {wat.get('n_waters_retained')} waters retained, {info}")
    if not (pdbqt.is_file() and pdbqt.stat().st_size > 0):
        raise SystemExit("receptor preparation failed")

    # --- dock ---------------------------------------------------------------
    box = out / "box.json"
    box.write_text(BOX_SRC.read_text(encoding="utf-8"), encoding="utf-8", newline="\n")
    scores = out / "scores.csv"
    poses = out / "poses"
    print("\ndocking the panel ...", flush=True)
    r = wsl(f'python 10_new_system/scripts/g4_dock_library.py '
            f'--ledger "{ledger.relative_to(ROOT).as_posix()}" '
            f'--receptor "{pdbqt.relative_to(ROOT).as_posix()}" '
            f'--box "{box.relative_to(ROOT).as_posix()}" '
            f'--out "{scores.relative_to(ROOT).as_posix()}" '
            f'--poses-dir "{poses.relative_to(ROOT).as_posix()}" '
            f'--id-col compound_id --smiles-col SMILES '
            f'--n-conformers {N_CONFORMERS} --exhaustiveness {EXHAUSTIVENESS} '
            f'--n-poses {N_POSES} --seed {SEED} --cpu 6')
    if not scores.is_file():
        raise SystemExit(f"docking produced nothing:\n{(r.stderr or r.stdout)[-500:]}")

    # --- score, with the control as a gate ----------------------------------
    import warnings; warnings.filterwarnings("ignore")
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdMolAlign
    from posebusters import PoseBusters
    from scipy.stats import spearmanr
    RDLogger.DisableLog("rdApp.*")

    smiles = json.loads(TEMPLATE_SRC.read_text(encoding="utf-8"))["smiles"]
    ref = CD.crystal_ligand(src, rec_dir / "ZMA_crystal.pdb", smiles)
    sc = pd.read_csv(scores).set_index("compound_id")
    pb = PoseBusters("dock")
    LOAD = {"mol_pred_loaded", "mol_cond_loaded", "sanitization", "inchi_convertible"}

    res = {}
    for r_ in panel.itertuples():
        cid = r_.compound_id
        sdf = poses / f"{cid}.sdf"
        entry = {"chembl_id": r_.chembl_id, "class": r_._5,
                 "pchembl_median": r_.pchembl_median,
                 "vina_score": (float(sc.loc[cid, "vina_score"])
                                if cid in sc.index and pd.notna(sc.loc[cid, "vina_score"])
                                else None)}
        if sdf.is_file():
            entry["asn253"] = asn253_contact(sdf, filt)
            try:
                df = pb.bust([sdf], None, filt)
                checks = [c for c in df.columns if c not in LOAD]
                entry["posebusters_rank1_valid"] = bool(df.iloc[0][checks].all())
            except Exception as exc:
                entry["posebusters_rank1_valid"] = None
            if cid == CONTROL:
                ms = [m for m in Chem.SDMolSupplier(str(sdf), removeHs=True) if m]
                entry["rmsd_to_crystal_A"] = round(rdMolAlign.CalcRMS(ms[0], ref), 3)
        res[cid] = entry

    ctrl = res.get(CONTROL, {})
    ctrl_ok = (ctrl.get("rmsd_to_crystal_A") is not None
               and ctrl["rmsd_to_crystal_A"] <= CONTROL_MAX_RMSD_A)

    known = [(v_["pchembl_median"], v_["vina_score"]) for v_ in res.values()
             if v_["class"] == "reference ligand"
             and v_["pchembl_median"] is not None and v_["vina_score"] is not None]
    rho = pv = None
    if len(known) >= 3:
        rho, pv = spearmanr([k[0] for k in known], [k[1] for k in known])

    report = {
        "question": "What does the docking score say about compounds whose "
                    "affinity is already known?",
        "authority": "ADR-AIIM-163, pre-registered before any docking was run",
        "recipe": manifest["recipe"],
        "control": {"compound": CONTROL, "max_rmsd_A": CONTROL_MAX_RMSD_A,
                    "observed_rmsd_A": ctrl.get("rmsd_to_crystal_A"),
                    "passed": bool(ctrl_ok)},
        "by_compound": res,
    }
    if not ctrl_ok:
        report["outcome"] = "CONTROL_FAILED"
        report["results_withheld"] = (
            f"{CONTROL} redocked at {ctrl.get('rmsd_to_crystal_A')} A against a "
            f"ceiling of {CONTROL_MAX_RMSD_A} A. No other number here is reported.")
        print(f"\nCONTROL FAILED -- results withheld")
    else:
        vs = [v_["vina_score"] for v_ in res.values() if v_["vina_score"] is not None]
        kn = [v_ for v_ in res.values() if v_["class"] == "reference ligand"
              and v_["vina_score"] is not None]
        cd = [v_ for v_ in res.values() if v_["class"].startswith("candidate")
              and v_["vina_score"] is not None]
        report["outcome"] = "MEASURED"
        report["summary"] = {
            "n_known": len(kn), "n_candidates": len(cd),
            "pchembl_range": [min(k["pchembl_median"] for k in kn),
                              max(k["pchembl_median"] for k in kn)],
            "vina_range_known": [round(min(k["vina_score"] for k in kn), 3),
                                 round(max(k["vina_score"] for k in kn), 3)],
            "vina_range_candidates": [round(min(c["vina_score"] for c in cd), 3),
                                      round(max(c["vina_score"] for c in cd), 3)],
            "spearman_score_vs_pchembl": (round(float(rho), 3) if rho is not None else None),
            "spearman_n": len(known),
            "no_significance_claim": "n is 10; the correlation is descriptive "
                                     "and is reported with its n",
        }
        print(f"\ncontrol {ctrl['rmsd_to_crystal_A']} A -- passed")
        print(f"\n{'compound':<16}{'pChEMBL':>9}{'vina':>9}{'Asn253':>9}  valid")
        for cid, v_ in sorted(res.items(),
                              key=lambda kv: (kv[1]["vina_score"] is None,
                                              kv[1]["vina_score"] or 0)):
            a = v_.get("asn253", {})
            print(f"{cid:<16}"
                  f"{(v_['pchembl_median'] if v_['pchembl_median'] is not None else float('nan')):>9.2f}"
                  f"{(v_['vina_score'] if v_['vina_score'] is not None else float('nan')):>9.3f}"
                  f"{(a.get('min_distance_A') if a.get('min_distance_A') is not None else float('nan')):>9.2f}"
                  f"  {v_.get('posebusters_rank1_valid')}")
        s = report["summary"]
        print(f"\npChEMBL {s['pchembl_range'][0]}-{s['pchembl_range'][1]} "
              f"spans the known set; their Vina scores span "
              f"{s['vina_range_known'][0]} to {s['vina_range_known'][1]}")
        print(f"candidates: {s['vina_range_candidates'][0]} to "
              f"{s['vina_range_candidates'][1]}")
        print(f"Spearman score vs pChEMBL: {s['spearman_score_vs_pchembl']} "
              f"(n = {s['spearman_n']}, descriptive)")

    (out / "clinical_reference.json").write_text(json.dumps(report, indent=2),
                                                 encoding="utf-8")
    print(f"\nwritten -> {out / 'clinical_reference.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

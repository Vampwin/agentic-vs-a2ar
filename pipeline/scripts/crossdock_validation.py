"""crossdock_validation.py -- does pose accuracy survive a change of receptor?

Criteria are module constants transcribed from ADR-AIIM-162 **before any docking
was run**, the convention R2, R3, R4 and ADR-AIIM-161 followed.

The pipeline's whole geometry validation is one number on one molecule:
ZM241385 redocked into 4EIY at 0.379 A. A ligand redocked into the structure
solved with it is the easiest case that exists. This measures the next case out
-- same ligand, same protocol, a receptor conformation solved with something
else.

**Superposition, reused not rewritten.** The five frozen structures are not in a
common frame (pairwise CA RMSD without superposition runs 0.64 to 101 A). Each
is carried onto 4EIY by `kabsch()` from `multi_reference_ifp_variant.py`, which
this project already wrote and used for its second IFP reference.

**All five receptors are prepared identically, by obabel alone.** `pdbfixer`,
which the production G2 step used, is installed in no environment here. What has
to be constant for this comparison is the preparation across the five, not
agreement with production -- so 4EIY is re-docked under this same preparation as
the within-study reference, and the production 0.379 A is not treated as
comparable.

Usage:
    python crossdock_validation.py [--stage all|prep|dock|score] [--out DIR]
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import subprocess
import sys
from pathlib import Path
from kdss_paths import project_root, resolve, under, wsl_root

import numpy as np

# --- pre-registered, ADR-AIIM-162 -------------------------------------------
RMSD_PASS_A = 2.0            # this project's frozen G3 criterion
N_CONFORMERS = 8
EXHAUSTIVENESS = 32
N_POSES = 10
SEED = 20260729
RECEPTORS = ["3EML", "4EIY", "5OLG", "5OLH", "6MH8"]
REFERENCE_FRAME = "4EIY"
LIGAND_RESN = "ZMA"
# 5OLH was solved with 9XT, so it has no ZM241385 of its own. Its reference is
# 4EIY's crystal ZMA carried into 5OLH's frame -- constructed, not observed, and
# named here because it is the one comparison of which that is true.
NO_NATIVE_LIGAND = {"5OLH"}
# Every one of these structures is a fusion construct, and the fusion is numbered
# from 1000 in the same chain as the receptor: T4 lysozyme 1002-1161 in 3EML,
# BRIL 1001-1106 in 4EIY. Matching CA atoms by residue number therefore pairs
# T4 lysozyme with BRIL -- two unrelated proteins -- and the first run of this
# script superposed 3EML onto 4EIY at 8.94 A because of it. The receptor proper
# ends before 1000 in all five.
RECEPTOR_RESI_MAX = 999


def _root() -> Path:
    # Delegates to kdss_paths, which finds whichever layout this file is inside
    # -- the development tree or a deposit -- by walking up from its own
    # location. It previously tried two absolute paths and exited if neither
    # existed, which on the author's machine meant a copy of the deposit read
    # and wrote the development tree instead of itself.
    return project_root()


ROOT = _root()
SCRIPTS = resolve("10_new_system/scripts")
STRUCT = resolve("02_frozen_inputs/a2ar_structures")
PROD_BOX = resolve("20_new_runs/snp_screen_2026-08-04/G2_preparation/docking_box.json")
ZMA_TEMPLATE = resolve("20_new_runs/snp_screen_2026-08-04/G3_redocking/reference_template.json")
WSL_DISTRO = "Ubuntu-22.04"
CONDA_ENV = "diffdock"


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


MR = _load("multi_reference_ifp_variant", under("10_new_system/scripts", "multi_reference_ifp_variant.py"))


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def wsl(cmd: str, timeout=7200):
    full = (f"source ~/miniconda3/etc/profile.d/conda.sh && conda activate {CONDA_ENV} "
            f"&& cd {wsl_root()!r} && {cmd}")
    return subprocess.run(["wsl", "-d", WSL_DISTRO, "--", "bash", "-lc", full],
                          capture_output=True, text=True, timeout=timeout)


def transform_pdb(src: Path, dst: Path, R, t, mode: str):
    """Write src with every coordinate mapped by (R, t).

    mode "protein"  -- ATOM records only, for the receptor
    mode "ligand"   -- the LIGAND_RESN HETATM block only

    The first version took a keep_het flag and, when it was set, wrote the
    protein *and* the ligand into the file called `..._crystal.pdb`. RDKit then
    read a 3,037-atom molecule as the reference. Substructure matching still
    found the right 25 atoms for four receptors and returned None for the fifth,
    which is how it surfaced.
    """
    if mode not in ("protein", "ligand"):
        raise ValueError(mode)
    out = []
    for line in src.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            resn = line[17:20].strip()
            if mode == "protein" and not line.startswith("ATOM"):
                continue
            if mode == "ligand" and not (line.startswith("HETATM")
                                         and resn == LIGAND_RESN):
                continue
            # Alternate conformations. 4EIY carries 93 atoms as altloc A and the
            # same 93 as B; keeping both puts two copies of a side chain a
            # fraction of an angstrom apart, and the receptor then clashes with
            # itself. This is what broke the first run: ZM241385 redocked into
            # its own structure at 8.691 A instead of 0.398 A. pdbfixer, which
            # the production step used and which is not installed here, resolves
            # altlocs as a side effect; doing it explicitly is the part of that
            # preparation this comparison actually needed.
            if line[16] not in (" ", "A"):
                continue
            v = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            v = R @ v + t
            line = f"{line[:30]}{v[0]:8.3f}{v[1]:8.3f}{v[2]:8.3f}{line[54:]}"
            out.append(line)
        elif line.startswith("TER") and mode == "protein":
            out.append(line)
    dst.write_text("\n".join(out) + "\nEND\n", encoding="utf-8", newline="\n")


def stage_prep(out: Path) -> dict:
    """Superpose every receptor onto 4EIY, strip to protein, convert to pdbqt."""
    ref_atoms = MR.read_pdb(STRUCT / f"{REFERENCE_FRAME}.pdb")
    ref_ca = MR.ca_by_resi(ref_atoms)
    prep = out / "receptors"
    prep.mkdir(parents=True, exist_ok=True)
    info = {}
    for pdb in RECEPTORS:
        src = STRUCT / f"{pdb}.pdb"
        atoms = MR.read_pdb(src)
        ca = MR.ca_by_resi(atoms)
        shared = sorted((r for r in set(ca) & set(ref_ca)
                         if int(r) <= RECEPTOR_RESI_MAX), key=lambda r: int(r))
        if len(shared) < 100:
            raise SystemExit(f"{pdb}: only {len(shared)} shared CA with {REFERENCE_FRAME}")
        P = np.array([ca[r] for r in shared])
        Q = np.array([ref_ca[r] for r in shared])
        R, t = MR.kabsch(P, Q) if pdb != REFERENCE_FRAME else (np.eye(3), np.zeros(3))
        resid = np.sqrt((((P @ R.T + t) - Q) ** 2).sum(1).mean())

        prot = prep / f"{pdb}_super.pdb"
        transform_pdb(src, prot, R, t, mode="protein")
        lig = prep / f"{pdb}_{LIGAND_RESN}_crystal.pdb"
        has_lig = any(a["resn"] == LIGAND_RESN
                      for a in MR.read_pdb(src, hetres=LIGAND_RESN))
        if has_lig:
            transform_pdb(src, lig, R, t, mode="ligand")

        pdbqt = prep / f"{pdb}_prepared.pdbqt"
        r = wsl(f'obabel -ipdb "{prot.relative_to(ROOT).as_posix()}" -opdbqt '
                f'-O "{pdbqt.relative_to(ROOT).as_posix()}" -xr -p 7.4')
        ok = pdbqt.is_file() and pdbqt.stat().st_size > 0
        info[pdb] = {
            "source_sha256": sha256(src),
            "n_shared_CA": len(shared),
            "superposition_rmsd_A": round(float(resid), 3),
            "has_native_ZMA": has_lig,
            "receptor_pdbqt": str(pdbqt.relative_to(ROOT).as_posix()) if ok else None,
            "obabel_stderr": (r.stderr or "").strip()[-200:] if not ok else "",
        }
        print(f"  {pdb}: {len(shared)} shared CA, superposition {resid:.3f} A, "
              f"native ZMA {'yes' if has_lig else 'no'}, pdbqt {'ok' if ok else 'FAILED'}")
    (out / "prep.json").write_text(json.dumps(info, indent=2), encoding="utf-8")
    return info


def stage_dock(out: Path, info: dict):
    """Dock ZM241385 into each superposed receptor, production protocol and box."""
    smiles = json.loads(ZMA_TEMPLATE.read_text(encoding="utf-8"))["smiles"]
    ledger = out / "ligand_ledger.csv"
    ledger.write_text(f"compound_id,SMILES\nZM241385,{smiles}\n",
                      encoding="utf-8", newline="\n")
    box = out / "box.json"
    box.write_text(PROD_BOX.read_text(encoding="utf-8"), encoding="utf-8", newline="\n")
    for pdb in RECEPTORS:
        rec = info[pdb]["receptor_pdbqt"]
        if not rec:
            print(f"  {pdb}: no receptor, skipped")
            continue
        dest = out / f"scores_{pdb}.csv"
        poses = out / f"poses_{pdb}"
        if dest.is_file() and dest.stat().st_size > 0:
            print(f"  {pdb}: already docked")
            continue
        poses.mkdir(parents=True, exist_ok=True)
        cmd = (f'python 10_new_system/scripts/g4_dock_library.py '
               f'--ledger "{ledger.relative_to(ROOT).as_posix()}" '
               f'--receptor "{rec}" --box "{box.relative_to(ROOT).as_posix()}" '
               f'--out "{dest.relative_to(ROOT).as_posix()}" '
               f'--poses-dir "{poses.relative_to(ROOT).as_posix()}" '
               f'--id-col compound_id --smiles-col SMILES '
               f'--n-conformers {N_CONFORMERS} --exhaustiveness {EXHAUSTIVENESS} '
               f'--n-poses {N_POSES} --seed {SEED} --cpu 6')
        print(f"  docking into {pdb} ...", flush=True)
        r = wsl(cmd)
        if not dest.is_file():
            print(f"    FAILED: {(r.stderr or r.stdout)[-300:]}")


def _crystal_mol(pdb_path: Path, smiles: str):
    """The crystal ligand with bond orders assigned from the frozen template.
    A PDB HETATM block carries no bond orders; inferring them instead of taking
    them from the registered SMILES is guessing at chemistry."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    raw = Chem.MolFromPDBFile(str(pdb_path), removeHs=True, sanitize=False)
    if raw is None:
        return None
    tmpl = Chem.MolFromSmiles(smiles)
    try:
        return AllChem.AssignBondOrdersFromTemplate(tmpl, raw)
    except Exception:
        return None


def stage_score(out: Path, info: dict) -> dict:
    """RMSD of the rank-1 pose to that receptor's own crystal ligand."""
    import warnings
    warnings.filterwarnings("ignore")
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdMolAlign
    from posebusters import PoseBusters
    RDLogger.DisableLog("rdApp.*")

    smiles = json.loads(ZMA_TEMPLATE.read_text(encoding="utf-8"))["smiles"]
    prep = out / "receptors"
    ref_for = {}
    for pdb in RECEPTORS:
        p = prep / f"{pdb}_{LIGAND_RESN}_crystal.pdb"
        if p.is_file():
            ref_for[pdb] = _crystal_mol(p, smiles)
    # 5OLH has no ZMA of its own; its reference is 4EIY's, already in the shared
    # frame. Named in ADR-AIIM-162 as the one constructed reference.
    for pdb in NO_NATIVE_LIGAND:
        ref_for.setdefault(pdb, ref_for.get(REFERENCE_FRAME))

    pb = PoseBusters("dock")
    LOAD = {"mol_pred_loaded", "mol_cond_loaded", "sanitization", "inchi_convertible"}
    results = {}
    for pdb in RECEPTORS:
        sdf = out / f"poses_{pdb}" / "ZM241385.sdf"
        rec = info.get(pdb, {}).get("receptor_pdbqt")
        if not sdf.is_file():
            results[pdb] = {"available": False, "reason": "no poses written"}
            continue
        mols = [m for m in Chem.SDMolSupplier(str(sdf), removeHs=True) if m]
        ref = ref_for.get(pdb)
        rmsds = []
        for m in mols:
            try:
                rmsds.append(round(rdMolAlign.CalcRMS(m, ref), 3))
            except Exception:
                rmsds.append(None)
        good = [r for r in rmsds if r is not None]
        prot = prep / f"{pdb}_super.pdb"
        try:
            df = pb.bust([sdf], None, prot)
            checks = [c for c in df.columns if c not in LOAD]
            valid = [bool(df.iloc[i][checks].all()) for i in range(len(df))]
        except Exception as exc:
            valid = [None] * len(mols)
        results[pdb] = {
            "available": True,
            "n_poses": len(mols),
            "reference": ("own crystal ZM241385" if pdb not in NO_NATIVE_LIGAND
                          else "4EIY crystal ZM241385 carried into this frame"),
            "superposition_rmsd_A": info.get(pdb, {}).get("superposition_rmsd_A"),
            "rank1_rmsd_A": rmsds[0] if rmsds else None,
            "best_of_10_rmsd_A": min(good) if good else None,
            "all_rmsd_A": rmsds,
            "rank1_within_2A": (rmsds[0] <= RMSD_PASS_A) if rmsds and rmsds[0] else False,
            "posebusters_valid_rank1": valid[0] if valid else None,
            "posebusters_n_valid": sum(1 for v in valid if v),
        }
        r = results[pdb]
        print(f"  {pdb}: rank-1 RMSD {r['rank1_rmsd_A']} A  "
              f"(best of {r['n_poses']}: {r['best_of_10_rmsd_A']})  "
              f"{'PASS' if r['rank1_within_2A'] else 'fail'} at {RMSD_PASS_A} A  "
              f"valid {r['posebusters_n_valid']}/{r['n_poses']}")

    passing = [p for p in RECEPTORS
               if results.get(p, {}).get("rank1_within_2A")]
    n = len(passing)
    outcome = ("HOLDS" if n == 5 else "FAILS" if n <= 1 else "DEGRADES")
    report = {
        "question": "Does pose accuracy survive a change of receptor conformation?",
        "authority": "ADR-AIIM-162, pre-registered before any docking was run",
        "criteria": {"rmsd_pass_A": RMSD_PASS_A,
                     "source_of_threshold": "this project's frozen G3 criterion",
                     "HOLDS": "5 of 5", "DEGRADES": "2 to 4 of 5",
                     "FAILS": "0 or 1 of 5"},
        "preparation": "all five receptors by obabel alone; pdbfixer, used by the "
                       "production G2 step, is installed in no environment here. "
                       "4EIY is re-docked under this same preparation as the "
                       "within-study reference and the production 0.379 A is not "
                       "treated as comparable.",
        "superposition": "kabsch() from multi_reference_ifp_variant.py, matched "
                         f"chain-A CA with residue number <= {RECEPTOR_RESI_MAX}. "
                         "The fusion partners (T4 lysozyme in 3EML, BRIL in 4EIY) "
                         "are numbered from 1000 in the same chain and matching "
                         "across them superposed 3EML at 8.94 A.",
        "by_receptor": results,
        "n_passing": n,
        "outcome": outcome,
    }
    (out / "crossdock_results.json").write_text(json.dumps(report, indent=2),
                                               encoding="utf-8")
    print(f"\n{n} of 5 rank-1 poses within {RMSD_PASS_A} A  ->  {outcome}")
    return report


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all",
                    choices=["all", "prep", "dock", "score"])
    ap.add_argument("--out", default=str(resolve("20_new_runs/crossdock_validation_2026-09-03")))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    info = None
    if args.stage in ("all", "prep"):
        print("stage: prep (superpose onto 4EIY, strip to protein, obabel -> pdbqt)")
        info = stage_prep(out)
    if info is None and (out / "prep.json").is_file():
        info = json.loads((out / "prep.json").read_text(encoding="utf-8"))
    if args.stage in ("all", "dock"):
        print("\nstage: dock (ZM241385 into each, production protocol and box)")
        stage_dock(out, info)
    if args.stage in ("all", "score"):
        print("\nstage: score (RMSD to crystal, PoseBusters on every pose)")
        stage_score(out, info)
    print(f"\nartifacts -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

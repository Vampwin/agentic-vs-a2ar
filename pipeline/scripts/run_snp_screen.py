"""
run_snp_screen.py — the runner that drives the audited SNP-Lib screen.

This is the piece that was missing. `WorkflowOrchestrator.run_stage` expects a
`tool_fn(run_dir, config_path, registry) -> dict`, and nothing implemented that
signature, so no stage had ever been executed under contract and the append-only
artifact registry had never written a record.

Each `tool_*` adapter below runs one deterministic domain step and returns the
paths and values its gate contract declares as inputs. The orchestrator hashes
every returned path into the registry, then runs the gate. The adapters never
decide whether a stage passed — that is the evaluator's job, and the evaluator
cannot be instructed by the code that produced the artifacts.

Two rules the adapters must not break:

1. **Every file read goes through the registry access log.** G5a's label
   isolation check reads that log as its source of truth, so a read that
   bypasses it is invisible to the audit and silently weakens the assay-blind
   claim.

2. **The G3 atom map is registered before docking runs.** The map defines which
   reference atom corresponds to which pose atom. Computing it from the docked
   coordinates would let the correspondence be chosen to flatter the RMSD, so it
   is derived from the crystal ligand and the *input* conformer, written, and
   registered before any pose exists.

Usage
-----
    python scripts/run_snp_screen.py --approver "Full Name" [--stages G0,G1]
    python scripts/run_snp_screen.py --approver "Full Name" --run-id my_run

Stages G4-G8 are not implemented here yet (Week 3); requesting them raises.
"""

from __future__ import annotations

import argparse
import json
import logging
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml

# Make the package importable when run as a script from anywhere.
_SYSTEM_ROOT = Path(__file__).resolve().parent.parent
if str(_SYSTEM_ROOT) not in sys.path:
    sys.path.insert(0, str(_SYSTEM_ROOT))

from kdss_audit.artifact_registry import (  # noqa: E402
    ArtifactEntry,
    ArtifactRegistry,
    FileAccessEntry,
    sha256_file,
)
from kdss_audit.orchestrator import WorkflowOrchestrator  # noqa: E402
from kdss_audit.state_machine import AuditState  # noqa: E402

logger = logging.getLogger("run_snp_screen")

ROOT = _SYSTEM_ROOT.parent                       # D:\Agentic VS Workflow
FROZEN = ROOT / "02_frozen_inputs"
LEGACY = ROOT / "04_legacy_evidence"
CONTRACTS_DIR = _SYSTEM_ROOT / "contracts"
CONFIG_PATH = _SYSTEM_ROOT / "configs" / "a2a_frozen_config.yaml"

# WSL relay: Vina, RDKit and meeko live in the `diffdock` conda env on the
# Ubuntu-22.04 distro. Windows Python has RDKit but no Vina.
WSL_DISTRO = "Ubuntu-22.04"
WSL_CONDA_ENV = "diffdock"


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _record_read(registry: ArtifactRegistry, path: Path, purpose: str,
                 reader: str = "run_snp_screen") -> Path:
    """Log a file read. Returns the path so it can be used inline."""
    registry.record_file_access(
        FileAccessEntry(reader=reader, file_path=str(path), purpose=purpose)
    )
    return path


def _register(registry: ArtifactRegistry, artifact_id: str, path: Path, stage: str,
              description: str, tool: str, tool_version: str = "",
              command: str = "", seed: int | None = None,
              conda_env: str = "", decision: str = "") -> str:
    """Register an artifact explicitly, when ordering or metadata matters."""
    registry.register_artifact(ArtifactEntry(
        artifact_id=artifact_id,
        path=str(path),
        sha256=sha256_file(path),
        stage=stage,
        description=description,
        tool=tool,
        tool_version=tool_version,
        conda_env=conda_env,
        command=command,
        random_seed=seed,
        agent_decision_summary=decision,
    ))
    return artifact_id


def _win_to_wsl(p: Path) -> str:
    """C:\\x -> /mnt/c/x, preserving spaces (caller must quote)."""
    s = str(Path(p).resolve())
    drive, rest = s[0].lower(), s[2:].replace("\\", "/")
    return f"/mnt/{drive}{rest}"


def _wsl_python(script: str, purpose: str,
                env: str = WSL_CONDA_ENV) -> subprocess.CompletedProcess:
    """Run a Python snippet inside a named conda env on WSL."""
    cmd = (
        f"source ~/miniconda3/etc/profile.d/conda.sh && "
        f"conda activate {env} && python -"
    )
    logger.info("WSL(%s): %s", env, purpose)
    return subprocess.run(
        ["wsl", "-d", WSL_DISTRO, "--", "bash", "-lc", cmd],
        input=script, text=True, capture_output=True, timeout=1800,
    )


def _load_config(config_path: Path) -> dict[str, Any]:
    with open(config_path, encoding="utf-8") as f:
        return yaml.safe_load(f)


# ─────────────────────────────────────────────────────────────────────────────
# G0 — input provenance
# ─────────────────────────────────────────────────────────────────────────────

def tool_g0_verify_inputs(run_dir: Path, config_path: Path,
                          registry: ArtifactRegistry) -> dict[str, Any]:
    """
    Resolve the frozen inputs. No transformation happens here: G0 exists to
    establish that what the pipeline is about to consume is byte-identical to
    what was registered, so this adapter only locates files and logs the reads.
    """
    inputs = {
        "snp_library_csv": FROZEN / "snp_library" / "01_library_validated.csv",
        "snp_library_sdf": FROZEN / "snp_library" / "01_library_validated.sdf",
        "snp_library_raw_xlsx": FROZEN / "snp_library" /
            "Sritularak Natural Products Library (SNP-Lib) - Prof. Boonchoo Sritularak_update.xlsx",
        "receptor_pdb_4eiy": FROZEN / "a2ar_structures" / "4EIY.pdb",
        "receptor_pdb_3eml": FROZEN / "a2ar_structures" / "3EML.pdb",
        "chembl_benchmark_csv": FROZEN / "public_benchmarks" / "chembl_a2ar" /
            "chembl_a2ar_validation_set.csv",
        "chembl_a2ar_raw_csv": FROZEN / "public_benchmarks" / "chembl_a2ar" /
            "chembl_a2a_raw_v2.csv",
    }
    for key, path in inputs.items():
        _record_read(registry, path, purpose=f"G0 provenance verification: {key}")
    return {k: str(v) for k, v in inputs.items()}


# ─────────────────────────────────────────────────────────────────────────────
# G1 — standardization
# ─────────────────────────────────────────────────────────────────────────────

def tool_g1_standardize(run_dir: Path, config_path: Path,
                        registry: ArtifactRegistry) -> dict[str, Any]:
    """
    Produce the 236-record accountability ledger.

    Reuses `scripts/01_standardize_library.py` verbatim rather than
    reimplementing salt stripping and stereo perception, rebinding only its
    output paths so results land in this run's directory instead of the
    module's hard-coded WP1 folder.
    """
    import importlib.util

    src = _SYSTEM_ROOT / "scripts" / "01_standardize_library.py"
    _record_read(registry, src, purpose="G1 standardization implementation")
    _record_read(registry, FROZEN / "snp_library" / "01_library_validated.csv",
                 purpose="G1 standardization input")

    spec = importlib.util.spec_from_file_location("_standardize_library", src)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)

    out_dir = run_dir / "G1_standardization"
    out_dir.mkdir(parents=True, exist_ok=True)
    mod.OUTPUT_DIR = out_dir
    mod.OUTPUT_LEDGER = out_dir / "snp_source_ledger.csv"
    mod.OUTPUT_META = out_dir / "snp_source_ledger_meta.json"

    try:
        mod.run()
    except SystemExit as exc:
        # The module exits non-zero when the ledger row count is wrong. That is
        # a finding for G1 to report, not a reason to abort the runner.
        if exc.code not in (0, None):
            logger.warning("standardization exited with code %s", exc.code)

    return {
        "snp_source_ledger_csv": str(mod.OUTPUT_LEDGER),
        "snp_source_ledger_meta_json": str(mod.OUTPUT_META),
    }


# ─────────────────────────────────────────────────────────────────────────────
# G2 — receptor / ligand preparation
# ─────────────────────────────────────────────────────────────────────────────

def _write_receptor_with_waters(
    pdb_path: Path, out_pdb: Path, chain: str, cutoff: float,
) -> dict[str, Any]:
    """
    Protein chain + the ordered waters that sit within `cutoff` of the
    co-crystal ligand.

    Retaining these is not a detail. ZM241385 binds A2AR through a bridging
    water network, and removing it moved rank-1 redocking from 0.28 A to
    8.67 A with every other variable held constant. Waters beyond the cutoff
    are dropped so the receptor is not padded with bulk solvent.
    """
    import math

    lines = pdb_path.read_text(encoding="utf-8", errors="replace").splitlines()

    def xyz(l: str) -> tuple[float, float, float]:
        return (float(l[30:38]), float(l[38:46]), float(l[46:54]))

    def usable(l: str) -> bool:
        return l[16] in (" ", "A")          # primary altloc only

    def is_hydrogen(l: str) -> bool:
        return (l[76:78].strip() or l[12:16].strip()[:1]).upper().startswith("H")

    ligand = [xyz(l) for l in lines
              if l.startswith(("ATOM", "HETATM")) and l[17:20].strip() == "ZMA"
              and usable(l) and not is_hydrogen(l)]
    if not ligand:
        raise RuntimeError(f"No ZMA heavy atoms in {pdb_path}")

    groups: dict[tuple[str, str], list[Any]] = {}
    for l in lines:
        if not l.startswith("HETATM") or l[17:20].strip() != "HOH" or not usable(l):
            continue
        d = min(math.dist(xyz(l), z) for z in ligand)
        key = (l[21], l[22:27])
        entry = groups.setdefault(key, [d, []])
        entry[0] = min(entry[0], d)
        entry[1].append(l)

    kept_keys = [k for k, (d, _) in groups.items() if d <= cutoff]
    kept = [l for k in kept_keys for l in groups[k][1]]
    protein = [l for l in lines
               if l.startswith("ATOM") and l[21] == chain and usable(l)]
    out_pdb.write_text("\n".join(protein + kept) + "\nEND\n", encoding="utf-8")

    return {
        "cutoff_angstrom": cutoff,
        "reference_ligand_resname": "ZMA",
        "n_waters_in_file": len(groups),
        "n_waters_retained": len(kept_keys),
        "n_water_atoms_retained": len(kept),
        "n_protein_atoms": len(protein),
        "chain": chain,
    }


def _zma_centroid(pdb_path: Path) -> tuple[list[float], int]:
    """Heavy-atom centroid of the co-crystal ZMA ligand, in the receptor frame."""
    coords: list[tuple[float, float, float]] = []
    with open(pdb_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if line[17:20].strip() != "ZMA":
                continue
            element = line[76:78].strip() or line[12:16].strip()[:1]
            if element.upper().startswith("H"):
                continue
            coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    if not coords:
        raise RuntimeError(f"No ZMA heavy atoms found in {pdb_path}")
    n = len(coords)
    centre = [round(sum(c[i] for c in coords) / n, 3) for i in range(3)]
    return centre, n


def tool_g2_prepare(run_dir: Path, config_path: Path,
                    registry: ArtifactRegistry) -> dict[str, Any]:
    """
    Prepare receptors and define the docking box.

    The box centre is computed here from the crystallographic ZMA centroid
    rather than copied from the config, so the value the contract sees is
    derived from structure geometry and is reproducible from the PDB alone.
    That is precisely what `no_assay_constants_in_box` is protecting: a box
    quietly recentred on a compound of interest is undetectable downstream.
    """
    cfg = _load_config(config_path)
    out_dir = run_dir / "G2_preparation"
    out_dir.mkdir(parents=True, exist_ok=True)

    # A box centre is only meaningful in its own structure's coordinate frame,
    # so it is derived from that structure's own co-crystal ligand rather than
    # carried over from another crystal.
    centres: dict[str, list[float]] = {}
    n_ref_atoms: dict[str, int] = {}
    for label in ("4EIY",):
        pdb = _record_read(registry, FROZEN / "a2ar_structures" / f"{label}.pdb",
                           purpose=f"G2 box derivation from crystal ZMA ({label})")
        centres[label], n_ref_atoms[label] = _zma_centroid(pdb)

    size = [cfg["vina"]["box_size_x"], cfg["vina"]["box_size_y"],
            cfg["vina"]["box_size_z"]]
    box = {
        "center": centres["4EIY"],          # primary structure
        "size": size,
        "centers_by_receptor": centres,
        "n_reference_atoms_by_receptor": n_ref_atoms,
        "pdb_used": "4EIY (primary and only docking structure)",
        "centroid_reference": "ZMA heavy-atom centroid (co-crystal ligand)",
        "derived_from": "crystallographic geometry only",
        "config_center_for_comparison": [
            cfg["vina"]["box_center_x"], cfg["vina"]["box_center_y"],
            cfg["vina"]["box_center_z"],
        ],
    }
    box_path = out_dir / "docking_box.json"
    box_path.write_text(json.dumps(box, indent=2), encoding="utf-8")

    # Receptors are prepared here rather than inherited. The legacy PDBQT files
    # were produced in 2026-06 outside this protocol and were never validated
    # under it; carrying them forward would leave the provenance of every
    # docking score resting on an unregistered conversion.
    # Only the primary structure is prepared. 3EML was excluded on 2026-08-04
    # (ADR-AIIM-036) after failing redocking at every water shell tested; a
    # receptor that cannot reproduce a known co-crystal pose cannot be used to
    # rank unknown compounds.
    water_cfg = cfg["receptor"]["structural_waters"]
    chain = cfg["receptor"]["primary"]["chain"]
    prepared: dict[str, Path] = {}
    prep_reports: dict[str, dict[str, Any]] = {}
    water_reports: dict[str, dict[str, Any]] = {}

    for label in ("4EIY",):
        src_pdb = _record_read(registry, FROZEN / "a2ar_structures" / f"{label}.pdb",
                               purpose=f"G2 receptor preparation input ({label})")
        filtered_pdb = out_dir / f"{label}_protein_plus_waters.pdb"
        fixed_pdb = out_dir / f"{label}_fixed.pdb"
        pdbqt = out_dir / f"{label}_prepared.pdbqt"

        water_reports[label] = _write_receptor_with_waters(
            src_pdb, filtered_pdb, chain=chain,
            cutoff=float(water_cfg["cutoff_angstrom"]),
        )
        logger.info("  %s retained %d/%d waters within %s A",
                    label, water_reports[label]["n_waters_retained"],
                    water_reports[label]["n_waters_in_file"],
                    water_cfg["cutoff_angstrom"])

        # PDBFixer lives in mdenv; Vina tooling lives in diffdock. The input is
        # already filtered to chain + retained waters, so keepWater=True keeps
        # exactly those and nothing else.
        fix = _wsl_python(
            "from pdbfixer import PDBFixer\n"
            "from openmm.app import PDBFile\n"
            f"f = PDBFixer(filename={_win_to_wsl(filtered_pdb)!r})\n"
            "f.findMissingResidues(); f.missingResidues = {}\n"
            "f.findNonstandardResidues(); f.replaceNonstandardResidues()\n"
            "f.removeHeterogens(keepWater=True)\n"
            "f.findMissingAtoms(); f.addMissingAtoms()\n"
            f"f.addMissingHydrogens({cfg['receptor']['protonation_ph']})\n"
            f"PDBFile.writeFile(f.topology, f.positions, open({_win_to_wsl(fixed_pdb)!r}, 'w'))\n"
            "nw = sum(1 for r in f.topology.residues() if r.name in ('HOH','WAT'))\n"
            "n_res = sum(1 for _ in f.topology.residues())\n"
            "n_at = sum(1 for _ in f.topology.atoms())\n"
            "print('residues:', n_res, 'waters:', nw, 'atoms:', n_at)\n",
            purpose=f"PDBFixer receptor preparation ({label}, waters retained)",
            env="mdenv",
        )
        if fix.returncode != 0 or not fixed_pdb.exists():
            raise RuntimeError(f"receptor fixing failed for {label}:\n{fix.stdout}\n{fix.stderr}")
        logger.info("  %s %s", label, fix.stdout.strip().splitlines()[-1])

        conv = subprocess.run(
            ["wsl", "-d", WSL_DISTRO, "--", "bash", "-lc",
             "source ~/miniconda3/etc/profile.d/conda.sh && "
             f"conda activate {WSL_CONDA_ENV} && "
             f'obabel -ipdb "{_win_to_wsl(fixed_pdb)}" -opdbqt '
             f'-O "{_win_to_wsl(pdbqt)}" -xr'],
            capture_output=True, text=True, timeout=900,
        )
        if not pdbqt.exists() or pdbqt.stat().st_size == 0:
            raise RuntimeError(f"receptor PDBQT conversion failed for {label}:\n"
                               f"{conv.stdout}\n{conv.stderr}")

        types: dict[str, int] = {}
        for line in pdbqt.read_text(encoding="utf-8", errors="replace").splitlines():
            if line.startswith(("ATOM", "HETATM")):
                t = line[77:79].strip()
                types[t] = types.get(t, 0) + 1
        prep_reports[label] = {"atom_type_counts": types,
                               "n_atoms": sum(types.values()),
                               "polar_hydrogens": types.get("HD", 0),
                               "structural_waters": water_reports[label]}
        prepared[label] = pdbqt
        _register(
            registry, artifact_id=f"G2::prepared_receptor_{label.lower()}",
            path=pdbqt, stage="G2",
            description=f"receptor prepared in-run from frozen {label}.pdb",
            tool="PDBFixer + Open Babel 3.1.0",
            command=f"PDBFixer(chain A, removeHeterogens, addH pH "
                    f"{cfg['receptor']['protonation_ph']}) | obabel -ipdb -opdbqt -xr",
            decision="prepared here rather than inherited so every docking score "
                     "traces to a registered conversion",
        )

    rec_4eiy = prepared["4EIY"]

    ref_lig_src = _record_read(registry, LEGACY / "redocking" / "ZMA_crystal_4EIY.sdf",
                               purpose="G2 reference ligand extracted from 4EIY")
    ref_lig = out_dir / "ZMA_reference_4EIY.sdf"
    shutil.copy2(ref_lig_src, ref_lig)

    cofactors = {
        "HOH": ("retained within %s A of the co-crystal ligand; bulk solvent "
                "beyond that cutoff removed" % water_cfg["cutoff_angstrom"]),
        "ZMA": "removed: co-crystal ligand removed from receptor; retained separately "
               "as the redocking reference",
        "other_heterogens": "removed: PDBFixer removeHeterogens(keepWater=False)",
    }

    return {
        "prepared_receptor_4eiy_pdbqt": str(rec_4eiy),
        "reference_ligand_4eiy_sdf": str(ref_lig),
        "docking_box_json": str(box_path),
        "cofactor_decisions": cofactors,
        "receptor_preparation_report": prep_reports,
        "structural_waters_recorded": water_reports,
        "receptor_preparation_provenance": (
            "prepared in this run from the frozen PDB files: PDBFixer "
            "(chain A, heterogens removed, hydrogens added at the configured pH) "
            "followed by Open Babel PDBQT conversion"
        ),
    }


# ─────────────────────────────────────────────────────────────────────────────
# G3 — redocking validation
# ─────────────────────────────────────────────────────────────────────────────

_ATOM_MAP_SCRIPT = r'''
import json, sys
from rdkit import Chem
from rdkit.Chem import rdFMCS

ref_path, probe_path, out_path, label, template_smiles = sys.argv[1:6]


def read_like_evaluator(path):
    """
    Read exactly as G3_evaluator does: SDMolSupplier(removeHs=True).

    The atom map is an index correspondence, so it is only valid under the same
    reader. A different parse (or a different removeHs outcome) renumbers the
    atoms and the map silently points at the wrong ones.
    """
    mol = next((m for m in Chem.SDMolSupplier(path, removeHs=True) if m is not None), None)
    if mol is None:
        # Crystal ligands extracted from PDB carry no bond orders and can fail
        # sanitization; fall back to an unsanitized read, then drop hydrogens
        # explicitly so the numbering still matches.
        mol = next((m for m in Chem.SDMolSupplier(path, removeHs=False, sanitize=False)
                    if m is not None), None)
        if mol is not None:
            try:
                mol = Chem.RemoveHs(mol, sanitize=False)
            except Exception:
                pass
    return mol


ref = read_like_evaluator(ref_path)
probe = read_like_evaluator(probe_path)
if ref is None or probe is None:
    raise SystemExit(f"could not read reference or probe for {label}")

res = rdFMCS.FindMCS(
    [ref, probe],
    ringMatchesRingOnly=True,
    completeRingsOnly=True,
    atomCompare=rdFMCS.AtomCompare.CompareElements,
    bondCompare=rdFMCS.BondCompare.CompareAny,
    timeout=60,
)
patt = Chem.MolFromSmarts(res.smartsString)
ref_match = ref.GetSubstructMatch(patt)
probe_match = probe.GetSubstructMatch(patt)
if not ref_match or not probe_match or len(ref_match) != len(probe_match):
    raise SystemExit(f"MCS produced no usable correspondence for {label}")

if max(probe_match) >= probe.GetNumAtoms() or max(ref_match) >= ref.GetNumAtoms():
    raise SystemExit(
        f"{label}: atom map indices out of range "
        f"(ref {ref.GetNumAtoms()} atoms, probe {probe.GetNumAtoms()} atoms)"
    )

payload = {
    "method": "rdkit_MCS",
    "label": label,
    "template_smiles": template_smiles,
    "n_ref_atoms": ref.GetNumAtoms(),
    "n_probe_atoms": probe.GetNumAtoms(),
    "smarts": res.smartsString,
    "n_mapped_atoms": len(ref_match),
    "parameters": {
        "ringMatchesRingOnly": True,
        "completeRingsOnly": True,
        "atomCompare": "CompareElements",
        "bondCompare": "CompareAny",
    },
    "atom_map": {str(r): int(p) for r, p in zip(ref_match, probe_match)},
    "derivation": (
        "Computed from the crystal reference ligand and the docking INPUT "
        "conformer, before any pose was generated. Vina preserves input atom "
        "ordering, so this correspondence holds for the docked poses."
    ),
}
with open(out_path, "w") as fh:
    json.dump(payload, fh, indent=2)
print(f"{label}: mapped {len(ref_match)} atoms")
'''

_MULTICONF_DOCK_SCRIPT = r'''
import json, sys
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from meeko import MoleculePreparation, PDBQTMolecule
from vina import Vina

(receptor, confs_sdf, out_sdf, energy_json, template_smiles,
 cx, cy, cz, sx, sy, sz, seed, exhaustiveness, n_poses) = sys.argv[1:15]

tmpl = Chem.MolFromSmiles(template_smiles)
tmpl_formula = rdMolDescriptors.CalcMolFormula(tmpl)
confs = [m for m in Chem.SDMolSupplier(confs_sdf, removeHs=False) if m is not None]
if not confs:
    raise SystemExit("no conformers to dock")

pool = []
for ci, conf in enumerate(confs):
    lig = Chem.AddHs(conf, addCoords=True)
    prep = MoleculePreparation()
    prep.prepare(lig)
    tmp_lig = "/tmp/_g3_lig.pdbqt"
    with open(tmp_lig, "w") as fh:
        fh.write(prep.write_pdbqt_string())

    v = Vina(sf_name="vina", seed=int(seed), verbosity=0)
    v.set_receptor(receptor)
    v.set_ligand_from_file(tmp_lig)
    v.compute_vina_maps(center=[float(cx), float(cy), float(cz)],
                        box_size=[float(sx), float(sy), float(sz)])
    v.dock(exhaustiveness=int(exhaustiveness), n_poses=int(n_poses))
    tmp_out = "/tmp/_g3_out.pdbqt"
    v.write_poses(tmp_out, n_poses=int(n_poses), overwrite=True)
    es = [float(e[0]) for e in v.energies()]

    for k, pose in enumerate(PDBQTMolecule.from_file(tmp_out, skip_typing=True)):
        mols = pose.export_rdkit_mol()
        m = mols[0] if isinstance(mols, (list, tuple)) else mols
        m = Chem.RemoveHs(m)
        if rdMolDescriptors.CalcMolFormula(m) != tmpl_formula:
            raise SystemExit("pose formula %s != template %s"
                             % (rdMolDescriptors.CalcMolFormula(m), tmpl_formula))
        pool.append((es[k] if k < len(es) else 0.0, ci, m))

# Energy-ranked across all conformers: pose 1 of this file is the pooled rank-1.
pool.sort(key=lambda t: t[0])
w = Chem.SDWriter(out_sdf)
for e, ci, m in pool:
    m.SetProp("vina_energy", str(e))
    m.SetProp("conformer_index", str(ci))
    w.write(m)
w.close()

with open(energy_json, "w") as fh:
    json.dump({"energies": [e for e, _, _ in pool],
               "conformer_index": [c for _, c, _ in pool],
               "n_conformers": len(confs),
               "n_pooled_poses": len(pool)}, fh, indent=2)
print("conformers: %d  pooled poses: %d  best energy: %.3f"
      % (len(confs), len(pool), pool[0][0]))
'''

_REDOCK_SCRIPT = r'''
import sys, json
from vina import Vina

receptor_pdbqt, ligand_pdbqt, out_pdbqt, cx, cy, cz, sx, sy, sz, seed, exhaustiveness, n_poses = sys.argv[1:13]
v = Vina(sf_name="vina", seed=int(seed), verbosity=0)
v.set_receptor(receptor_pdbqt)
v.set_ligand_from_file(ligand_pdbqt)
v.compute_vina_maps(center=[float(cx), float(cy), float(cz)],
                    box_size=[float(sx), float(sy), float(sz)])
v.dock(exhaustiveness=int(exhaustiveness), n_poses=int(n_poses))
v.write_poses(out_pdbqt, n_poses=int(n_poses), overwrite=True)
print(json.dumps({"energies": [float(e[0]) for e in v.energies()]}))
'''


def tool_g3_redock(run_dir: Path, config_path: Path,
                   registry: ArtifactRegistry) -> dict[str, Any]:
    """
    Redock the co-crystal ligand and validate the pose against the reference.

    Ordering is the point of this stage. The MCS atom map is derived from the
    crystal ligand and the docking input conformer, written, and registered
    *before* Vina is invoked, so the atom correspondence cannot have been
    influenced by the poses it will later be used to score.
    """
    cfg = _load_config(config_path)
    out_dir = run_dir / "G3_redocking"
    out_dir.mkdir(parents=True, exist_ok=True)

    crystal_4eiy = _record_read(registry, LEGACY / "redocking" / "ZMA_crystal_4EIY.sdf",
                                purpose="G3 reference ligand 4EIY")

    # The standardized library, for the representativeness assessment.
    ledger_csv = run_dir / "G1_standardization" / "snp_source_ledger.csv"
    library_smiles: list[str] = []
    if ledger_csv.is_file():
        import pandas as pd
        _record_read(registry, ledger_csv,
                     purpose="G3 control-representativeness: screening library")
        ledger = pd.read_csv(ledger_csv)
        col = next((c for c in ("smiles_canonical", "SMILES_std", "SMILES")
                    if c in ledger.columns), None)
        if col:
            library_smiles = [s for s in ledger[col].tolist() if isinstance(s, str)]
    # Recorded as inspected-and-rejected, not used: this file encodes C16N6O3.
    _record_read(registry, LEGACY / "redocking" / "ZM241385_conformer.sdf",
                 purpose="G3 legacy conformer inspected; rejected as wrong molecule")

    # ── Step 0: prepare the ligand, then round-trip it through the same
    # PDBQT->SDF conversion the poses will take.
    #
    # The atom map is an index correspondence, so it is only valid against
    # molecules in the same atom order. Vina preserves the order of its input
    # PDBQT, but the obabel conversion applied to the poses may not preserve the
    # order of the original conformer SDF. Deriving the map against a probe that
    # has been through the identical conversion keeps the indices aligned, and
    # still uses no pose information.
    # The conformer SDF carries 3D coordinates but no reliable bond orders, so
    # RDKit re-perceives its chemistry on every read. Left uncorrected, meeko
    # types an N-oxide where ZM241385 has an NH linker and the molecule that
    # gets docked is C16N6O3 rather than C16N7O2 — a different compound. The
    # template makes reference, docking input and poses one molecule.
    template_smiles = cfg["g3_redocking"]["reference_ligand_smiles"]
    template_path = out_dir / "reference_template.json"
    template_path.write_text(json.dumps({
        "name": cfg["g3_redocking"]["reference_ligand_name"],
        "smiles": template_smiles,
        "formula": cfg["g3_redocking"]["reference_ligand_formula"],
        "purpose": ("bond-order template for the crystal reference, the docking "
                    "input, and the docked poses"),
        "source": "frozen config g3_redocking.reference_ligand_smiles",
    }, indent=2), encoding="utf-8")
    _register(
        registry, artifact_id="G3::reference_template", path=template_path,
        stage="G3", description="pre-registered ZM241385 bond-order template",
        tool="frozen config",
        decision="registered before docking; fixes the chemistry of every "
                 "molecule G3 compares so RMSD is computed on one compound",
    )

    # The docking input is built from the declared SMILES rather than from the
    # legacy conformer file. That file contains C16N6O3 — one nitrogen is
    # recorded as oxygen — so it is not ZM241385, and anything docked from it
    # is a different compound. Generating the conformer here makes the docked
    # species verifiable against a declared structure instead of inherited from
    # an unchecked artifact, and it is what a redocking control should do
    # anyway: start from a generated geometry, not the answer.
    probe_sdf = out_dir / "ZM241385_conformers.sdf"
    seed = cfg["vina"]["seed"]
    n_conf = int(cfg["vina"]["n_conformers"])
    prune = float(cfg["vina"]["conformer_prune_rms"])
    prep = _wsl_python(
        "from rdkit import Chem\n"
        "from rdkit.Chem import AllChem, rdMolDescriptors\n"
        f"m = Chem.AddHs(Chem.MolFromSmiles({template_smiles!r}))\n"
        "ps = AllChem.ETKDGv3()\n"
        f"ps.randomSeed = {seed}\n"
        f"ps.pruneRmsThresh = {prune}\n"
        f"cids = AllChem.EmbedMultipleConfs(m, numConfs={n_conf}, params=ps)\n"
        "if not len(cids):\n"
        "    raise SystemExit('conformer embedding failed')\n"
        "AllChem.MMFFOptimizeMoleculeConfs(m)\n"
        "heavy = Chem.RemoveHs(m)\n"
        f"w = Chem.SDWriter({_win_to_wsl(probe_sdf)!r})\n"
        "for ci in cids: w.write(heavy, confId=ci)\n"
        "w.close()\n"
        "print('conformers:', len(cids), rdMolDescriptors.CalcMolFormula(heavy))\n",
        purpose=f"generate {n_conf} ZM241385 conformers from the declared SMILES",
    )
    if prep.returncode != 0 or not probe_sdf.exists():
        raise RuntimeError(f"conformer generation failed:\n{prep.stdout}\n{prep.stderr}")
    logger.info("  %s", prep.stdout.strip().splitlines()[-1] if prep.stdout.strip() else "")
    _register(
        registry, artifact_id="G3::docking_input_conformers", path=probe_sdf,
        stage="G3",
        description=f"{n_conf} ZM241385 conformers generated from the declared SMILES",
        tool="RDKit ETKDGv3 + MMFF94", conda_env=WSL_CONDA_ENV, seed=seed,
        decision="generated rather than reusing the legacy conformer (which encodes "
                 "C16N6O3, not ZM241385); multiple conformers because Vina holds "
                 "bond geometry rigid, making single-conformer results bimodal",
    )

    # ── Step 1: atom maps, registered BEFORE docking ─────────────────────────
    atom_maps: dict[str, Path] = {}
    for label, crystal in (("4EIY", crystal_4eiy),):
        dst = out_dir / f"atom_map_{label.lower()}.json"
        proc = _wsl_python(
            f"import sys; sys.argv = ['x', {_win_to_wsl(crystal)!r}, "
            f"{_win_to_wsl(probe_sdf)!r}, {_win_to_wsl(dst)!r}, {label!r}, "
            f"{template_smiles!r}]\n"
            + _ATOM_MAP_SCRIPT,
            purpose=f"MCS atom map for {label}",
        )
        if proc.returncode != 0 or not dst.exists():
            raise RuntimeError(
                f"atom map generation failed for {label}:\n{proc.stdout}\n{proc.stderr}"
            )
        logger.info("  %s", proc.stdout.strip())
        atom_maps[label] = dst
        _register(
            registry,
            artifact_id=f"G3::atom_map_{label.lower()}",
            path=dst, stage="G3", description=f"pre-docking MCS atom map ({label})",
            tool="rdkit.Chem.rdFMCS.FindMCS", conda_env=WSL_CONDA_ENV,
            decision="registered before docking so the atom correspondence "
                     "cannot be chosen from the poses it will score",
        )

    # ── Step 2: dock every conformer, pool poses by energy ───────────────────
    #
    # Vina samples torsions but keeps bond lengths, angles and ring geometry
    # fixed at whatever the input conformer had, so a single starting geometry
    # can cap the score the native binding mode is able to reach. Docking each
    # conformer separately and ranking the combined pose set by energy removes
    # that dependence; pose 1 of the pooled SDF is the pooled rank-1 and stays
    # the sole pass criterion.
    box = json.loads((run_dir / "G2_preparation" / "docking_box.json").read_text())
    sx, sy, sz = box["size"]
    exh = cfg["vina"]["exhaustiveness"]
    n_poses = cfg["vina"]["n_poses"]

    poses: dict[str, Path] = {}
    energies: dict[str, list[float]] = {}
    for label in ("4EIY",):
        receptor = run_dir / "G2_preparation" / f"{label}_prepared.pdbqt"
        _record_read(registry, receptor, purpose=f"G3 redocking receptor {label}")
        cx, cy, cz = box["centers_by_receptor"][label]
        out_sdf = out_dir / f"ZMA_redocked_{label}.sdf"
        energy_json = out_dir / f"ZMA_redocked_{label}_energies.json"

        proc = _wsl_python(
            "import sys; sys.argv = ['x', "
            + ", ".join(repr(x) for x in [
                _win_to_wsl(receptor), _win_to_wsl(probe_sdf),
                _win_to_wsl(out_sdf), _win_to_wsl(energy_json), template_smiles,
                str(cx), str(cy), str(cz), str(sx), str(sy), str(sz),
                str(seed), str(exh), str(n_poses),
            ])
            + "]\n" + _MULTICONF_DOCK_SCRIPT,
            purpose=f"multi-conformer Vina redocking into {label}",
        )
        if proc.returncode != 0 or not out_sdf.exists():
            raise RuntimeError(f"redocking failed for {label}:\n{proc.stdout}\n{proc.stderr}")
        logger.info("  %s", proc.stdout.strip().splitlines()[-1] if proc.stdout.strip() else "")
        try:
            energies[label] = json.loads(energy_json.read_text())["energies"]
        except Exception:
            energies[label] = []

        poses[label] = out_sdf
        _register(
            registry, artifact_id=f"G3::redocked_poses_{label.lower()}",
            path=out_sdf, stage="G3",
            description=f"Vina poses from {n_conf} conformers, pooled and ranked by energy ({label})",
            tool="AutoDock Vina (multi-conformer) + meeko", conda_env=WSL_CONDA_ENV,
            seed=seed,
            command=f"vina dock exhaustiveness={exh} n_poses={n_poses} "
                    f"center={[cx, cy, cz]} size={box['size']} "
                    f"x {n_conf} ETKDGv3 conformers, poses pooled by energy",
        )
        _register(
            registry, artifact_id=f"G3::pose_energies_{label.lower()}",
            path=energy_json, stage="G3",
            description=f"pooled pose energies and source conformer index ({label})",
            tool="AutoDock Vina", conda_env=WSL_CONDA_ENV, seed=seed,
        )

    # Structure-based transferability evidence (ADR-AIIM-073/074): the
    # ChEMBL A2AR benchmark (n=245) docked under this exact receptor
    # preparation -- box and every heavy-atom coordinate, including the
    # retained waters, verified byte-identical to this run's own G2 output
    # (ADR-AIIM-080) -- gives retrospective_enrichment a factual basis for
    # protocol_matches_screen rather than an assumed one. This is read-only
    # evidence already computed this session; it is not recomputed here.
    retro_path = (ROOT / "20_new_runs" / "chembl245_ifp_rescoring_2026-08-15"
                 / "ifp_rescoring_summary.json")
    retrospective_enrichment = None
    if retro_path.is_file():
        _record_read(registry, retro_path,
                     purpose="G3 retrospective_enrichment: ChEMBL A2AR benchmark evidence")
        vina_only = json.loads(retro_path.read_text())["methods"]["vina_only"]
        retrospective_enrichment = {
            "auc_roc": vina_only["auc_roc"],
            "bedroc20": vina_only["bedroc20"],
            "protocol_matches_screen": True,
            "n": vina_only["n"],
            "source": str(retro_path),
        }

    return {
        "reference_ligand_4eiy_sdf": str(crystal_4eiy),
        "redocked_poses_4eiy": str(poses["4EIY"]),
        "atom_map_4eiy_json": str(atom_maps["4EIY"]),
        "docking_input_conformer_sdf": str(probe_sdf),
        # Supplied so G3 can test whether a passing pose-validation result says
        # anything about the compounds this protocol would go on to rank.
        "control_ligand_smiles": [template_smiles],
        "screening_library_smiles": library_smiles,
        "retrospective_enrichment": retrospective_enrichment,
        "vina_energies": energies,
        # Both ligand files are expressed in their receptor's crystallographic
        # frame and no superposition is applied at any point in this stage.
        "frame_check_4eiy": True,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Stage table
# ─────────────────────────────────────────────────────────────────────────────

STAGES: list[tuple[str, str, Any]] = [
    ("G0", "G0_input.yaml", tool_g0_verify_inputs),
    ("G1", "G1_standardization.yaml", tool_g1_standardize),
    ("G2", "G2_preparation.yaml", tool_g2_prepare),
    ("G3", "G3_pose.yaml", tool_g3_redock),
]
_IMPLEMENTED = {name for name, _, _ in STAGES}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--approver", required=True,
                    help="Name of the human approving this scientific run. "
                         "Recorded in the registry; the frozen config requires "
                         "a human signature before any run.")
    ap.add_argument("--stages", default="G0,G1,G2,G3",
                    help="Comma-separated stages to run (default: G0,G1,G2,G3)")
    ap.add_argument("--run-id", default=None, help="Run directory name")
    ap.add_argument("--benchmark-role", default="Option_A_development",
                    choices=["Option_A_development", "Option_B_heldout"],
                    help="Declared role of the ChEMBL benchmark. Option A is the "
                         "truthful declaration: the set was used during legacy "
                         "method development, so it is not held out.")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)-7s %(message)s",
    )

    if not args.approver.strip():
        print("ERROR: --approver must be a real human name.", file=sys.stderr)
        return 2

    requested = [s.strip().upper() for s in args.stages.split(",") if s.strip()]
    unknown = [s for s in requested if s not in _IMPLEMENTED]
    if unknown:
        print(f"ERROR: stages not implemented in this runner: {unknown}. "
              f"G4-G8 arrive in Week 3.", file=sys.stderr)
        return 2

    run_id = args.run_id or f"snp_screen_{datetime.now(timezone.utc):%Y-%m-%dT%H%M%SZ}"
    run_dir = ROOT / "20_new_runs" / run_id
    run_dir.mkdir(parents=True, exist_ok=True)

    registry = ArtifactRegistry(run_dir / "artifact_registry.jsonl")
    registry.record_benchmark_role_declaration(
        role=args.benchmark_role,
        declared_by=args.approver,
        rationale=(
            "The ChEMBL A2AR set informed threshold selection during legacy "
            "method development, so it cannot be described as held out. "
            "Declaring Option A forces the development-benchmark caveat onto "
            "every metric G5b reports."
        ),
    )
    _record_read(registry, CONFIG_PATH, purpose="frozen configuration for this run")

    orch = WorkflowOrchestrator(
        run_dir=run_dir, contracts_dir=CONTRACTS_DIR,
        config_path=CONFIG_PATH, registry=registry,
    )

    print(f"\nRun directory : {run_dir}")
    print(f"Approver      : {args.approver}")
    print(f"Benchmark role: {args.benchmark_role}")
    print(f"Stages        : {', '.join(requested)}\n")

    outcomes: list[tuple[str, str, str]] = []
    halted = False
    for name, contract_file, tool_fn in STAGES:
        if name not in requested:
            continue
        print(f"-- {name} " + "-" * 56)
        try:
            result = orch.run_stage(name, contract_file, tool_fn)
        except Exception as exc:  # a tool failure is a run failure, not a pass
            logger.error("%s tool raised: %s", name, exc)
            outcomes.append((name, "TOOL_ERROR", str(exc)[:200]))
            halted = True
            break

        cr = result.contract_result
        state = cr.state.value if isinstance(cr.state, AuditState) else str(cr.state)
        summary = (f"{len(cr.checks_passed)} passed, {len(cr.checks_failed)} failed, "
                   f"{len(cr.checks_warned)} warned")
        outcomes.append((name, state, summary))
        print(f"   state   : {state}")
        print(f"   checks  : {summary}")
        if cr.checks_failed:
            print(f"   FAILED  : {cr.checks_failed}")
            for cid in cr.checks_failed[:5]:
                print(f"       {cid}: {str(cr.evidence.get(cid, ''))[:180]}")
        if state == AuditState.HALT.value:
            halted = True
            break

    print("\n" + "=" * 64)
    print(f"{'stage':8}{'state':14}detail")
    for name, state, detail in outcomes:
        print(f"{name:8}{state:14}{detail}")

    reg_path = run_dir / "artifact_registry.jsonl"
    n_records = sum(1 for _ in open(reg_path, encoding="utf-8")) if reg_path.exists() else 0
    print(f"\nartifact_registry.jsonl: {n_records} records -> {reg_path}")
    if halted:
        print("\nWorkflow halted. The gate that halted it is shown above; a human "
              "must adjudicate before this run can proceed.")
    return 1 if halted else 0


if __name__ == "__main__":
    raise SystemExit(main())

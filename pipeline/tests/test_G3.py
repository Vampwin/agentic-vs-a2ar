"""
Known-answer tests for the G3 evaluator — direct RMSD computation.

These tests use synthetic SDF files and pre-registered atom maps to verify
that the RMSD calculation is geometrically correct and uses the direct
(non-Kabsch) method.

Note: these tests require RDKit. Skip gracefully if not installed.
"""

import json
import math
import sys
import tempfile
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

rdkit = pytest.importorskip("rdkit", reason="RDKit not installed")
from rdkit import Chem
from rdkit.Chem import AllChem


def write_mol_sdf(mol: "Chem.Mol", path: Path) -> None:
    writer = Chem.SDWriter(str(path))
    writer.write(mol)
    writer.close()


def make_atom_map_json(path: Path, n_atoms: int) -> None:
    # Identity map: ref atom i → pose atom i
    atom_map = {str(i): i for i in range(n_atoms)}
    data = {"method": "rdkit_MCS", "mcs_smarts": "CCCC", "n_atoms": n_atoms, "atom_map": atom_map}
    with open(path, "w") as f:
        json.dump(data, f)


def mol_at_origin() -> "Chem.Mol":
    """Ethane with atoms at specific coordinates."""
    mol = Chem.MolFromSmiles("CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    mol = Chem.RemoveHs(mol)
    return mol


# ── Import the evaluator ──────────────────────────────────────────────────────
from evaluators.G3_evaluator import _compute_direct_rmsd


class TestDirectRMSD:

    def test_zero_rmsd_identical_positions(self, tmp_path):
        """Pose identical to reference must give RMSD = 0."""
        mol = mol_at_origin()
        ref_sdf = tmp_path / "ref.sdf"
        pose_sdf = tmp_path / "pose.sdf"
        atom_map_json = tmp_path / "atom_map.json"
        write_mol_sdf(mol, ref_sdf)
        write_mol_sdf(mol, pose_sdf)
        n_heavy = mol.GetNumHeavyAtoms()
        make_atom_map_json(atom_map_json, n_heavy)

        result = _compute_direct_rmsd(ref_sdf, pose_sdf, atom_map_json, "TEST")
        assert result.get("error") is None
        assert math.isclose(result["rank1_rmsd"], 0.0, abs_tol=1e-6)

    def test_known_displacement(self, tmp_path):
        """Shift each atom by exactly 1 Å along X. Expected RMSD = 1.0 Å."""
        from rdkit.Geometry import Point3D

        mol = mol_at_origin()
        ref_sdf = tmp_path / "ref.sdf"

        # Build shifted pose
        shifted_mol = Chem.RWMol(mol)
        conf = shifted_mol.GetConformer()
        for i in range(shifted_mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, Point3D(pos.x + 1.0, pos.y, pos.z))
        pose_sdf = tmp_path / "pose.sdf"
        write_mol_sdf(ref_sdf=ref_sdf, mol=mol, path=ref_sdf)  # wrong arg order
        write_mol_sdf(shifted_mol.GetMol(), pose_sdf)

        n_heavy = mol.GetNumHeavyAtoms()
        atom_map_json = tmp_path / "atom_map.json"
        make_atom_map_json(atom_map_json, n_heavy)

        result = _compute_direct_rmsd(ref_sdf, pose_sdf, atom_map_json, "TEST")
        assert result.get("error") is None
        assert math.isclose(result["rank1_rmsd"], 1.0, abs_tol=0.01)

    def test_rank1_is_first_pose(self, tmp_path):
        """Rank-1 RMSD must always use the first molecule in the SDF, not the best."""
        from rdkit.Geometry import Point3D

        mol = mol_at_origin()
        ref_sdf = tmp_path / "ref.sdf"
        write_mol_sdf(mol, ref_sdf)

        # First pose: shifted 3 Å; second pose: 0.5 Å (better but not rank-1)
        writer = Chem.SDWriter(str(tmp_path / "poses.sdf"))
        for shift in [3.0, 0.5]:
            m = Chem.RWMol(mol)
            conf = m.GetConformer()
            for i in range(m.GetNumAtoms()):
                pos = conf.GetAtomPosition(i)
                conf.SetAtomPosition(i, Point3D(pos.x + shift, pos.y, pos.z))
            writer.write(m.GetMol())
        writer.close()

        n_heavy = mol.GetNumHeavyAtoms()
        atom_map_json = tmp_path / "atom_map.json"
        make_atom_map_json(atom_map_json, n_heavy)

        result = _compute_direct_rmsd(ref_sdf, tmp_path / "poses.sdf", atom_map_json, "TEST")
        # rank-1 RMSD must be ~3.0, NOT ~0.5
        assert result["rank1_rmsd"] > 2.0
        assert result["best_rmsd"] < 1.0
        assert result["rank1_rmsd"] > result["best_rmsd"]

    def test_per_atom_distances_count(self, tmp_path):
        """Per-atom distances list must have length = number of atom pairs in map."""
        mol = mol_at_origin()
        ref_sdf = tmp_path / "ref.sdf"
        pose_sdf = tmp_path / "pose.sdf"
        write_mol_sdf(mol, ref_sdf)
        write_mol_sdf(mol, pose_sdf)
        n_heavy = mol.GetNumHeavyAtoms()
        atom_map_json = tmp_path / "atom_map.json"
        make_atom_map_json(atom_map_json, n_heavy)

        result = _compute_direct_rmsd(ref_sdf, pose_sdf, atom_map_json, "TEST")
        assert len(result["per_atom_distances_rank1"]) == n_heavy

    def test_missing_ref_returns_error(self, tmp_path):
        result = _compute_direct_rmsd(
            tmp_path / "nonexistent_ref.sdf",
            tmp_path / "pose.sdf",
            tmp_path / "map.json",
            "TEST",
        )
        assert result.get("error") is not None


# Helper overload fix (write_mol_sdf signature fix above)
def write_mol_sdf(mol, path, ref_sdf=None):  # noqa: F811
    if ref_sdf is not None:  # called with wrong keyword; ignore ref_sdf
        path = ref_sdf
    writer = Chem.SDWriter(str(path))
    writer.write(mol)
    writer.close()


# ── atom_map_pre_registration (held-out finding H03, ADR-AIIM-060) ─────────

import hashlib as _hashlib
import time as _time

from evaluators.G3_evaluator import evaluate as g3_evaluate
from kdss_audit.artifact_registry import ArtifactEntry, ArtifactRegistry


def _make_atom_map_and_pose(tmp_path):
    amap_path = tmp_path / "atom_map.json"
    amap_path.write_text(json.dumps({"atom_map": {"0": 0, "1": 1}, "method": "rdkit_MCS"}))
    pose_path = tmp_path / "pose.sdf"
    pose_path.write_text("mock pose\n$$$$\n")
    return amap_path, pose_path


def _registered_registry(tmp_path, amap_path, pose_path):
    registry = ArtifactRegistry(tmp_path / "registry.jsonl")
    registry.register_artifact(ArtifactEntry(
        artifact_id="G3::atom_map_4eiy", path=str(amap_path),
        sha256=_hashlib.sha256(amap_path.read_bytes()).hexdigest().upper(),
        stage="G3", description="atom map", tool="rdkit_MCS", tool_version="",
        conda_env="", command="", random_seed=None))
    _time.sleep(0.02)
    registry.register_artifact(ArtifactEntry(
        artifact_id="G3::redocked_poses_4eiy", path=str(pose_path),
        sha256=_hashlib.sha256(pose_path.read_bytes()).hexdigest().upper(),
        stage="G3", description="redocked poses", tool="vina", tool_version="",
        conda_env="", command="", random_seed=None))
    return registry


class TestAtomMapPreRegistration:

    def test_passes_when_registered_before_docking_and_unmodified(self, tmp_path):
        amap_path, pose_path = _make_atom_map_and_pose(tmp_path)
        registry = _registered_registry(tmp_path, amap_path, pose_path)
        contract = {"checks": [{"id": "atom_map_pre_registration", "severity": "critical"}]}
        result = g3_evaluate(contract, {"atom_map_4eiy_json": str(amap_path)}, registry)
        assert "atom_map_pre_registration" in result["passed"]

    def test_fails_on_posthoc_rewrite_after_registration(self, tmp_path):
        """
        Held-out case H03: the atom map is rewritten after both it and the
        docking result are registered, keeping valid JSON schema throughout.
        Before the fix, only file existence and schema were checked -- the
        rewrite passed silently (ADR-AIIM-060). It must now fail: the file's
        hash no longer matches what was registered.
        """
        amap_path, pose_path = _make_atom_map_and_pose(tmp_path)
        registry = _registered_registry(tmp_path, amap_path, pose_path)

        data = json.loads(amap_path.read_text())
        data["atom_map"] = {"0": 1, "1": 0}  # post-hoc remap, same valid schema
        amap_path.write_text(json.dumps(data))

        contract = {"checks": [{"id": "atom_map_pre_registration", "severity": "critical"}]}
        result = g3_evaluate(contract, {"atom_map_4eiy_json": str(amap_path)}, registry)
        assert "atom_map_pre_registration" in result["failed"]
        assert "atom_map_pre_registration" in result["critical_failures"]

    def test_fails_when_registered_after_docking(self, tmp_path):
        """Even an unmodified file should fail if it was registered too late."""
        amap_path, pose_path = _make_atom_map_and_pose(tmp_path)
        registry = ArtifactRegistry(tmp_path / "registry.jsonl")
        # Docking registered FIRST, atom map registered after -- wrong order.
        registry.register_artifact(ArtifactEntry(
            artifact_id="G3::redocked_poses_4eiy", path=str(pose_path),
            sha256=_hashlib.sha256(pose_path.read_bytes()).hexdigest().upper(),
            stage="G3", description="redocked poses", tool="vina", tool_version="",
            conda_env="", command="", random_seed=None))
        _time.sleep(0.02)
        registry.register_artifact(ArtifactEntry(
            artifact_id="G3::atom_map_4eiy", path=str(amap_path),
            sha256=_hashlib.sha256(amap_path.read_bytes()).hexdigest().upper(),
            stage="G3", description="atom map", tool="rdkit_MCS", tool_version="",
            conda_env="", command="", random_seed=None))

        contract = {"checks": [{"id": "atom_map_pre_registration", "severity": "critical"}]}
        result = g3_evaluate(contract, {"atom_map_4eiy_json": str(amap_path)}, registry)
        assert "atom_map_pre_registration" in result["failed"]

    def test_passes_when_no_registry_available(self, tmp_path):
        """No registry supplied (e.g. a lightweight unit test elsewhere in this
        suite): falls back to schema-only validation rather than failing
        closed on an ordering check it has no way to evaluate."""
        amap_path, _pose_path = _make_atom_map_and_pose(tmp_path)
        contract = {"checks": [{"id": "atom_map_pre_registration", "severity": "critical"}]}
        result = g3_evaluate(contract, {"atom_map_4eiy_json": str(amap_path)}, None)
        assert "atom_map_pre_registration" in result["passed"]

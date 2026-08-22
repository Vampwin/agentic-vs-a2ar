"""
Known-answer tests for the G0 evaluator.
"""

import hashlib
import sys
import tempfile
from pathlib import Path

import pytest
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))
from evaluators.G0_evaluator import evaluate
# The assay-blindness blocklist now lives in one shared module so G0, G4 and G6
# cannot drift apart; _has_blocked_column was folded into blocked_columns_in.
from evaluators._assay_blindness import blocked_columns_in


MINIMAL_CONTRACT = {
    "contract_id": "G0_test",
    "checks": [
        {"id": "compound_count", "expected_count": 5, "severity": "critical"},
        {"id": "id_uniqueness", "severity": "critical"},
        {"id": "smiles_column_present", "severity": "critical"},
        {"id": "no_assay_columns", "severity": "critical"},
    ],
}


def make_csv(tmp_path: Path, n: int = 5, extra_col: str | None = None) -> Path:
    data = {
        "SNP_ID": [f"SNP{i:05d}" for i in range(n)],
        "SMILES": ["CCO"] * n,
        "name": [f"compound_{i}" for i in range(n)],
    }
    if extra_col:
        data[extra_col] = ["x"] * n
    df = pd.DataFrame(data)
    p = tmp_path / "test_lib.csv"
    df.to_csv(p, index=False)
    return p


def test_pass_clean_csv(tmp_path):
    csv = make_csv(tmp_path, n=5)
    result = evaluate(
        contract=MINIMAL_CONTRACT,
        inputs={"snp_library_csv": str(csv)},
        registry=None,
    )
    assert "compound_count" in result["passed"]
    assert "id_uniqueness" in result["passed"]
    assert "smiles_column_present" in result["passed"]
    assert result["state"] == "PASS"


def test_fail_wrong_count(tmp_path):
    csv = make_csv(tmp_path, n=3)
    result = evaluate(
        contract=MINIMAL_CONTRACT,
        inputs={"snp_library_csv": str(csv)},
        registry=None,
    )
    assert "compound_count" in result["failed"]
    assert result["state"] == "HALT"


def test_fail_duplicate_ids(tmp_path):
    data = {
        "SNP_ID": ["SNP00001", "SNP00001", "SNP00002", "SNP00003", "SNP00004"],
        "SMILES": ["CCO"] * 5,
    }
    df = pd.DataFrame(data)
    p = tmp_path / "dupes.csv"
    df.to_csv(p, index=False)
    result = evaluate(
        contract=MINIMAL_CONTRACT,
        inputs={"snp_library_csv": str(p)},
        registry=None,
    )
    assert "id_uniqueness" in result["failed"]


def test_fail_blocked_column(tmp_path):
    csv = make_csv(tmp_path, n=5, extra_col="known_active")
    result = evaluate(
        contract=MINIMAL_CONTRACT,
        inputs={"snp_library_csv": str(csv)},
        registry=None,
    )
    assert "no_assay_columns" in result["failed"]
    assert result["state"] == "HALT"


def test_fail_camp_substring(tmp_path):
    csv = make_csv(tmp_path, n=5, extra_col="cAMP_score")
    result = evaluate(
        contract=MINIMAL_CONTRACT,
        inputs={"snp_library_csv": str(csv)},
        registry=None,
    )
    assert "no_assay_columns" in result["failed"]


def test_hash_check_correct(tmp_path):
    csv = make_csv(tmp_path, n=5)
    h = hashlib.sha256(csv.read_bytes()).hexdigest().upper()
    contract = {
        "contract_id": "G0_test",
        "checks": [
            {"id": "hash_snp_library_csv", "expected_sha256": h, "severity": "critical"},
        ],
    }
    result = evaluate(contract=contract, inputs={"snp_library_csv": str(csv)}, registry=None)
    assert "hash_snp_library_csv" in result["passed"]


def test_hash_check_mismatch(tmp_path):
    csv = make_csv(tmp_path, n=5)
    wrong_hash = "A" * 64
    contract = {
        "contract_id": "G0_test",
        "checks": [
            {"id": "hash_snp_library_csv", "expected_sha256": wrong_hash, "severity": "critical"},
        ],
    }
    result = evaluate(contract=contract, inputs={"snp_library_csv": str(csv)}, registry=None)
    assert "hash_snp_library_csv" in result["failed"]
    assert result["state"] == "HALT"


def test_file_not_found_halts(tmp_path):
    result = evaluate(
        contract=MINIMAL_CONTRACT,
        inputs={"snp_library_csv": str(tmp_path / "nonexistent.csv")},
        registry=None,
    )
    assert result["state"] == "HALT"


# ── chembl_a2ar_raw_csv (ADR-AIIM-049) ──────────────────────────────────────

def test_chembl_a2ar_raw_present_and_hash_pass(tmp_path):
    p = tmp_path / "raw.csv"
    p.write_text("molecule_chembl_id,pchembl_value\nCHEMBL1,7.0\n", encoding="utf-8")
    expected = hashlib.sha256(p.read_bytes()).hexdigest().upper()
    contract = {
        "contract_id": "G0_test",
        "checks": [
            {"id": "chembl_a2ar_raw_present", "severity": "high"},
            {"id": "hash_chembl_a2ar_raw_csv", "expected_sha256": expected, "severity": "high"},
        ],
    }
    result = evaluate(contract=contract, inputs={"chembl_a2ar_raw_csv": str(p)}, registry=None)
    assert "chembl_a2ar_raw_present" in result["passed"]
    assert "hash_chembl_a2ar_raw_csv" in result["passed"]


def test_chembl_a2ar_raw_hash_mismatch_fails(tmp_path):
    p = tmp_path / "raw.csv"
    p.write_text("molecule_chembl_id,pchembl_value\nCHEMBL1,7.0\n", encoding="utf-8")
    contract = {
        "contract_id": "G0_test",
        "checks": [{"id": "hash_chembl_a2ar_raw_csv", "expected_sha256": "A" * 64, "severity": "high"}],
    }
    result = evaluate(contract=contract, inputs={"chembl_a2ar_raw_csv": str(p)}, registry=None)
    assert "hash_chembl_a2ar_raw_csv" in result["failed"]


def test_chembl_benchmark_hash_no_longer_vacuously_passes(tmp_path):
    """
    Before ADR-AIIM-049, chembl_benchmark_hash had no expected_sha256 and
    passed regardless of file content. It must now behave like every other
    anchored hash check: unverifiable is a warning (severity=high, not
    critical), never a silent pass.
    """
    p = tmp_path / "chembl.csv"
    p.write_text("x\n", encoding="utf-8")
    contract = {
        "contract_id": "G0_test",
        "checks": [{"id": "chembl_benchmark_hash", "severity": "high"}],  # no expected_sha256
    }
    result = evaluate(contract=contract, inputs={"chembl_benchmark_csv": str(p)}, registry=None)
    assert "chembl_benchmark_hash" not in result["passed"]
    assert "chembl_benchmark_hash" in result["warned"]


# ── name_structure_mapping (held-out finding H01, ADR-AIIM-060) ────────────

def _write_sdf(path, ids, smiles_list):
    rdkit = pytest.importorskip("rdkit", reason="RDKit not installed")
    from rdkit import Chem
    writer = Chem.SDWriter(str(path))
    for sid, smi in zip(ids, smiles_list):
        mol = Chem.MolFromSmiles(smi)
        mol.SetProp("_Name", sid)
        writer.write(mol)
    writer.close()


def test_name_structure_mapping_passes_when_csv_matches_sdf(tmp_path):
    ids = [f"SNP{i:03d}" for i in range(5)]
    smiles = ["CCO", "c1ccccc1", "CC(=O)O", "CCN", "CCC"]
    csv = tmp_path / "lib.csv"
    pd.DataFrame({"SNP_ID": ids, "name": ids, "SMILES": smiles}).to_csv(csv, index=False)
    sdf = tmp_path / "lib.sdf"
    _write_sdf(sdf, ids, smiles)
    contract = {"contract_id": "G0_test",
                "checks": [{"id": "name_structure_mapping", "severity": "medium"}]}
    result = evaluate(contract=contract,
                       inputs={"snp_library_csv": str(csv), "snp_library_sdf": str(sdf)},
                       registry=None)
    assert "name_structure_mapping" in result["passed"]


def test_name_structure_mapping_fails_on_row_shifted_structures(tmp_path):
    """
    Held-out case H01: compound IDs and the 'name' column are untouched, but
    the structure column is shifted by one row relative to an independent SDF
    source. Before the fix, this check only verified a 'name' column existed
    and passed silently (ADR-AIIM-060). It must now fail.
    """
    ids = [f"SNP{i:03d}" for i in range(5)]
    smiles = ["CCO", "c1ccccc1", "CC(=O)O", "CCN", "CCC"]
    shifted = smiles[1:] + smiles[:1]
    csv = tmp_path / "lib.csv"
    pd.DataFrame({"SNP_ID": ids, "name": ids, "SMILES": shifted}).to_csv(csv, index=False)
    sdf = tmp_path / "lib.sdf"
    _write_sdf(sdf, ids, smiles)  # independent source, unshifted
    contract = {"contract_id": "G0_test",
                "checks": [{"id": "name_structure_mapping", "severity": "medium"}]}
    result = evaluate(contract=contract,
                       inputs={"snp_library_csv": str(csv), "snp_library_sdf": str(sdf)},
                       registry=None)
    assert "name_structure_mapping" in result["failed"]


def test_name_structure_mapping_falls_back_when_no_sdf(tmp_path):
    """No independent source to cross-check against: falls back to the
    column-presence check rather than failing or crashing."""
    ids = [f"SNP{i:03d}" for i in range(5)]
    csv = tmp_path / "lib.csv"
    pd.DataFrame({"SNP_ID": ids, "name": ids, "SMILES": ["CCO"] * 5}).to_csv(csv, index=False)
    contract = {"contract_id": "G0_test",
                "checks": [{"id": "name_structure_mapping", "severity": "medium"}]}
    result = evaluate(contract=contract, inputs={"snp_library_csv": str(csv)}, registry=None)
    assert "name_structure_mapping" in result["passed"]
    assert result["evidence"]["name_structure_mapping"]["cross_check"] == "not_run"

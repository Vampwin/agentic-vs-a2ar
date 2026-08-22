"""
WP1 standardization — SNP-Lib 236-compound ledger production.

Reads:  02_frozen_inputs/snp_library/01_library_validated.csv
Writes: 20_new_runs/WP0_WP1_audit_2026-07-29/snp_source_ledger.csv

Contract: G1_standardization_v1
Skills:   rdkit/SKILL.md (salt stripping, canonical SMILES, InChIKey, annotations)

Rules enforced:
- ADR-AIIM-006: drug-likeness is annotation-only; never drives exclusion
- Only valid exclusion_reason values: standardization_failure, mixture_excluded
  (embedding_failure and docking_failure are set by later stages)
- All 236 compounds must appear in the output ledger
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path
from datetime import datetime, timezone

import pandas as pd

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors, inchi as rd_inchi
    from rdkit.Chem.MolStandardize import rdMolStandardize
    RDKIT_OK = True
except ImportError:
    RDKIT_OK = False
    print("ERROR: RDKit not available. Install via: pip install rdkit", file=sys.stderr)
    sys.exit(1)


# ── Paths ──────────────────────────────────────────────────────────────────────
ROOT = Path(__file__).parent.parent.parent          # D:\Agentic VS Workflow
INPUT_CSV = ROOT / "02_frozen_inputs" / "snp_library" / "01_library_validated.csv"
OUTPUT_DIR = ROOT / "20_new_runs" / "WP0_WP1_audit_2026-07-29"
OUTPUT_LEDGER = OUTPUT_DIR / "snp_source_ledger.csv"
OUTPUT_META = OUTPUT_DIR / "snp_source_ledger_meta.json"

EXPECTED_COUNT = 236


# ── Drug-likeness exclusion keyword guard ──────────────────────────────────────
_DL_FORBIDDEN_REASONS = {
    "lipinski", "veber", "qed", "pains", "drug_like",
    "ro5", "hbd", "hba", "mw", "logp", "tpsa",
}


def _assert_no_dl_exclusion(reason: str) -> None:
    """Raise if any DL keyword appears in an exclusion_reason."""
    reason_lower = reason.lower()
    for kw in _DL_FORBIDDEN_REASONS:
        if kw in reason_lower:
            raise ValueError(
                f"ADR-AIIM-006 VIOLATION: DL keyword '{kw}' found in exclusion_reason '{reason}'. "
                "Drug-likeness must not drive exclusions."
            )


# ── RDKit helpers ──────────────────────────────────────────────────────────────
def _canonical(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def _inchikey(mol: Chem.Mol) -> str:
    try:
        inchi_str = rd_inchi.MolToInchi(mol)
        if inchi_str is None:
            return ""
        return rd_inchi.InchiToInchiKey(inchi_str) or ""
    except Exception:
        return ""


def _strip_salts(mol: Chem.Mol) -> tuple[Chem.Mol, bool, str]:
    """
    Return (parent_mol, had_salts, removed_fragment_smiles).
    Uses LargestFragmentChooser — keeps the biggest fragment.
    """
    chooser = rdMolStandardize.LargestFragmentChooser()
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) <= 1:
        return mol, False, ""
    parent = chooser.choose(mol)
    removed = [_canonical(f) for f in frags if _canonical(f) != _canonical(parent)]
    return parent, True, ".".join(removed)


def _is_mixture(mol: Chem.Mol) -> bool:
    """True if the mol still has multiple fragments after salt stripping."""
    frags = Chem.GetMolFrags(mol, asMols=True)
    return len(frags) > 1


def _stereo_info(mol: Chem.Mol) -> tuple[bool, int, bool]:
    """Return (has_stereocenters, n_stereocenters, all_defined)."""
    centers = Chem.FindMolChiralCenters(
        mol, includeUnassigned=True, useLegacyImplementation=False
    )
    n = len(centers)
    if n == 0:
        return False, 0, True
    undefined = sum(1 for _, chirality in centers if chirality == "?")
    return True, n, (undefined == 0)


def _descriptors(mol: Chem.Mol) -> dict:
    return {
        "mol_wt": round(Descriptors.MolWt(mol), 4),
        "logp": round(Descriptors.MolLogP(mol), 4),
        "hbd": rdMolDescriptors.CalcNumHBD(mol),
        "hba": rdMolDescriptors.CalcNumHBA(mol),
        "tpsa": round(rdMolDescriptors.CalcTPSA(mol), 4),
        "n_rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "n_rings": rdMolDescriptors.CalcNumRings(mol),
        "n_aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
    }


def _lipinski_ro5_pass(desc: dict) -> bool:
    return (
        desc["mol_wt"] <= 500
        and desc["logp"] <= 5
        and desc["hbd"] <= 5
        and desc["hba"] <= 10
    )


def _qed(mol: Chem.Mol) -> float:
    try:
        from rdkit.Chem import QED
        return round(QED.qed(mol), 4)
    except Exception:
        return float("nan")


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest().upper()


# ── Main standardization logic ─────────────────────────────────────────────────
def standardize_row(row: pd.Series) -> dict:
    snp_id = str(row["Compound"]).strip()
    name = str(row.get("Compound name", "")).strip()
    smiles_input = str(row.get("SMILES", "")).strip()
    pubchem_cid = row.get("Pubchem Compound ID", "")
    groups = str(row.get("Groups", "")).strip()
    source_no = int(row.get("No", 0))

    record: dict = {
        "SNP_ID": snp_id,
        "name": name,
        "source_row_no": source_no,
        "groups": groups,
        "pubchem_cid": pubchem_cid if pd.notna(pubchem_cid) else "",
        "smiles_input": smiles_input,
        # Populated below
        "valid_smiles": False,
        "standardization_failure_reason": "",
        "has_salts": False,
        "removed_fragment": "",
        "parent_smiles": "",
        "is_mixture": False,
        "smiles_canonical": "",
        "inchikey": "",
        "has_stereocenters": False,
        "n_stereocenters": 0,
        "stereo_defined": True,
        # Drug-likeness annotations (ADR-AIIM-006: annotation only)
        "mol_wt": float("nan"),
        "logp": float("nan"),
        "hbd": float("nan"),
        "hba": float("nan"),
        "tpsa": float("nan"),
        "n_rotatable_bonds": float("nan"),
        "n_rings": float("nan"),
        "n_aromatic_rings": float("nan"),
        "qed": float("nan"),
        "lipinski_ro5_pass": None,
        # Exclusion (only non-DL reasons permitted)
        "exclusion_reason": "",
        # Filled in post-processing
        "duplicate_inchikey_flag": False,
    }

    if not smiles_input:
        record["standardization_failure_reason"] = "empty_smiles"
        record["exclusion_reason"] = "standardization_failure"
        return record

    mol = Chem.MolFromSmiles(smiles_input)
    if mol is None:
        record["standardization_failure_reason"] = "rdkit_parse_fail"
        record["exclusion_reason"] = "standardization_failure"
        return record

    record["valid_smiles"] = True

    # Salt stripping
    try:
        parent, had_salts, removed_frag = _strip_salts(mol)
        record["has_salts"] = had_salts
        record["removed_fragment"] = removed_frag
    except Exception as exc:
        record["standardization_failure_reason"] = f"salt_strip_error:{exc}"
        record["exclusion_reason"] = "standardization_failure"
        record["valid_smiles"] = False
        return record

    # Mixture check (after stripping)
    if _is_mixture(parent):
        record["is_mixture"] = True
        record["exclusion_reason"] = "mixture_excluded"
        # Still compute canonical SMILES for the record
        try:
            record["parent_smiles"] = _canonical(parent)
        except Exception:
            pass
        return record

    # Canonical SMILES
    try:
        record["parent_smiles"] = _canonical(parent)
        record["smiles_canonical"] = record["parent_smiles"]
    except Exception as exc:
        record["standardization_failure_reason"] = f"canonical_smiles_error:{exc}"
        record["exclusion_reason"] = "standardization_failure"
        record["valid_smiles"] = False
        return record

    # InChIKey
    record["inchikey"] = _inchikey(parent)

    # Stereocenters
    has_sc, n_sc, all_def = _stereo_info(parent)
    record["has_stereocenters"] = has_sc
    record["n_stereocenters"] = n_sc
    record["stereo_defined"] = all_def

    # Descriptors (annotation only — ADR-AIIM-006)
    desc = _descriptors(parent)
    record.update(desc)
    record["qed"] = _qed(parent)
    record["lipinski_ro5_pass"] = _lipinski_ro5_pass(desc)

    # Guard: no DL keywords in exclusion_reason
    _assert_no_dl_exclusion(record["exclusion_reason"])

    return record


def run() -> None:
    if not INPUT_CSV.exists():
        print(f"ERROR: Input CSV not found: {INPUT_CSV}", file=sys.stderr)
        sys.exit(1)

    print(f"Reading: {INPUT_CSV}")
    src_df = pd.read_csv(INPUT_CSV, encoding="utf-8-sig")
    print(f"  Source rows: {len(src_df)}")

    if len(src_df) != EXPECTED_COUNT:
        print(
            f"WARNING: Expected {EXPECTED_COUNT} rows, got {len(src_df)}. "
            "Proceeding — G1 will flag this.",
            file=sys.stderr,
        )

    records = []
    for _, row in src_df.iterrows():
        rec = standardize_row(row)
        records.append(rec)
        status = "OK" if rec["valid_smiles"] and not rec["is_mixture"] else f"SKIP({rec['exclusion_reason']})"
        print(f"  {rec['SNP_ID']}  {status}")

    ledger = pd.DataFrame(records)

    # Post-processing: flag duplicate InChIKeys
    inchikey_counts = ledger[ledger["inchikey"] != ""]["inchikey"].value_counts()
    dup_inchikeys = set(inchikey_counts[inchikey_counts > 1].index)
    ledger["duplicate_inchikey_flag"] = ledger["inchikey"].isin(dup_inchikeys)

    # Column ordering for legibility
    col_order = [
        "SNP_ID", "name", "source_row_no", "groups", "pubchem_cid",
        "smiles_input", "valid_smiles", "standardization_failure_reason",
        "has_salts", "removed_fragment", "parent_smiles", "smiles_canonical",
        "is_mixture", "inchikey", "duplicate_inchikey_flag",
        "has_stereocenters", "n_stereocenters", "stereo_defined",
        "mol_wt", "logp", "hbd", "hba", "tpsa", "n_rotatable_bonds",
        "n_rings", "n_aromatic_rings", "qed", "lipinski_ro5_pass",
        "exclusion_reason",
    ]
    ledger = ledger[col_order]

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    ledger.to_csv(OUTPUT_LEDGER, index=False, encoding="utf-8")
    print(f"\nLedger written: {OUTPUT_LEDGER}")

    # Summary stats
    n_valid = int(ledger["valid_smiles"].sum())
    n_mixture = int(ledger["is_mixture"].sum())
    n_std_fail = int((ledger["exclusion_reason"] == "standardization_failure").sum())
    n_dockable = int(
        (ledger["valid_smiles"] & ~ledger["is_mixture"] & (ledger["exclusion_reason"] == "")).sum()
    )
    n_with_salts = int(ledger["has_salts"].sum())
    n_dup_ik = int(ledger["duplicate_inchikey_flag"].sum())
    n_stereo = int(ledger["has_stereocenters"].sum())
    n_undef_stereo = int((ledger["has_stereocenters"] & ~ledger["stereo_defined"]).sum())
    n_lipinski_pass = int(ledger["lipinski_ro5_pass"].sum()) if "lipinski_ro5_pass" in ledger else -1

    summary = {
        "script": "01_standardize_library.py",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "input_csv": str(INPUT_CSV),
        "input_sha256": _sha256(INPUT_CSV),
        "output_ledger": str(OUTPUT_LEDGER),
        "total_rows": len(ledger),
        "valid_smiles": n_valid,
        "standardization_failures": n_std_fail,
        "mixtures_excluded": n_mixture,
        "dockable_compounds": n_dockable,
        "compounds_with_salts": n_with_salts,
        "duplicate_inchikey_pairs": n_dup_ik,
        "compounds_with_stereocenters": n_stereo,
        "compounds_unassigned_stereo": n_undef_stereo,
        "lipinski_ro5_pass_annotation": n_lipinski_pass,
        "adraim006_dl_exclusions": 0,
        "note": (
            "lipinski_ro5_pass is an annotation column only (ADR-AIIM-006). "
            "No compound was excluded by drug-likeness criteria."
        ),
    }

    with open(OUTPUT_META, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print(f"Metadata written: {OUTPUT_META}")

    print("\n-- Summary --------------------------------------------------")
    print(f"  Total rows:                  {len(ledger)}")
    print(f"  Valid SMILES:                {n_valid}")
    print(f"  Standardization failures:    {n_std_fail}")
    print(f"  Mixtures excluded:           {n_mixture}")
    print(f"  Dockable (clean):            {n_dockable}")
    print(f"  Had counterions stripped:    {n_with_salts}")
    print(f"  Duplicate InChIKey pairs:    {n_dup_ik}")
    print(f"  Has stereocenters:           {n_stereo}")
    print(f"  Unassigned stereocenters:    {n_undef_stereo}")
    print(f"  Lipinski pass (annotation):  {n_lipinski_pass}")
    print(f"  DL exclusions:               0  (ADR-AIIM-006 enforced)")
    print("-------------------------------------------------------------")

    if len(ledger) != EXPECTED_COUNT:
        print(f"\nWARNING: Ledger has {len(ledger)} rows, expected {EXPECTED_COUNT}.", file=sys.stderr)
        sys.exit(2)

    print("\nStandardization complete. G1 contract check is next.")


if __name__ == "__main__":
    run()

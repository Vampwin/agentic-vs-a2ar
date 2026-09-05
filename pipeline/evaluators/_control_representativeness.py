"""
Validation representativeness — is the positive control like what we screen?

A redocking control answers a narrow question: does this protocol reproduce
*this ligand's* crystallographic pose. It is routinely read as the much larger
claim that the protocol is valid. Where the binding site's water structure or
induced fit is ligand-dependent, the first can hold while the second fails, and
nothing downstream detects the gap because the scores, poses and rankings are
all internally consistent.

This module measures the distance between the control ligands and the screening
library, so that a passing redocking control cannot be silently generalised to a
chemically unrelated set.

Three complementary signals are reported rather than one composite, because they
fail in different ways and a reader should see which one is responsible:

  scaffold coverage    fraction of the library sharing a Bemis-Murcko framework
                       with some control — asks whether the control belongs to
                       any chemotype actually present in the library
  similarity coverage  fraction of the library within a fingerprint-similarity
                       cutoff of its nearest control — asks how much of the
                       library the control speaks for at all
  descriptor position  how many physicochemical descriptors place the control
                       inside the library's central range — asks whether the
                       control is typical of the library or an outlier beside it

The thresholds are conventions, not derived constants, so they are contract-
configurable and the underlying numbers are always reported. The check is meant
to escalate to a human with the evidence, not to silently decide.
"""

from __future__ import annotations

from typing import Any, Iterable

# ECFP4 Tanimoto below roughly 0.3-0.4 is conventionally treated as dissimilar.
DEFAULT_SIMILARITY_CUTOFF = 0.35
# At least half the library should be at least weakly similar to some control
# before a control-based validation is assumed to transfer to it.
DEFAULT_MIN_SIMILARITY_COVERAGE = 0.50
DEFAULT_MIN_SCAFFOLD_COVERAGE = 0.10
# Percentile band defining "the library's central range" for descriptors.
DEFAULT_PERCENTILES = (5.0, 95.0)

_DESCRIPTORS = ("MW", "cLogP", "TPSA", "HBD", "HBA", "RotB", "AromaticRings", "HeavyAtoms")


def _descriptor_vector(mol: Any) -> dict[str, float] | None:
    try:
        from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors as rdMD
    except ImportError:
        return None
    if mol is None:
        return None
    try:
        return {
            "MW": float(Descriptors.MolWt(mol)),
            "cLogP": float(Crippen.MolLogP(mol)),
            "TPSA": float(rdMD.CalcTPSA(mol)),
            "HBD": float(rdMD.CalcNumHBD(mol)),
            "HBA": float(rdMD.CalcNumHBA(mol)),
            "RotB": float(rdMD.CalcNumRotatableBonds(mol)),
            "AromaticRings": float(rdMD.CalcNumAromaticRings(mol)),
            "HeavyAtoms": float(mol.GetNumHeavyAtoms()),
        }
    except Exception:
        return None


def _framework(mol: Any) -> str | None:
    """Bemis-Murcko framework (ring systems plus linkers, atom types generalised)."""
    try:
        from rdkit import Chem
        from rdkit.Chem.Scaffolds import MurckoScaffold
    except ImportError:
        return None
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        if scaffold is None or scaffold.GetNumAtoms() == 0:
            return ""      # acyclic: no ring framework, still a valid answer
        generic = MurckoScaffold.MakeScaffoldGeneric(scaffold)
        return Chem.MolToSmiles(generic)
    except Exception:
        return None


def assess(
    control_smiles: Iterable[str],
    library_smiles: Iterable[str],
    similarity_cutoff: float = DEFAULT_SIMILARITY_CUTOFF,
    percentiles: tuple[float, float] = DEFAULT_PERCENTILES,
) -> dict[str, Any]:
    """
    Compare control ligands with a screening library.

    Returns a report dict; the caller applies thresholds and decides the state.
    `available` is False when RDKit is missing, so a caller can warn rather than
    assert representativeness it could not actually measure.
    """
    try:
        import numpy as np
        from rdkit import Chem, DataStructs
        from rdkit.Chem import rdFingerprintGenerator
    except ImportError:
        return {"available": False, "reason": "RDKit/numpy not available"}

    controls = [(s, Chem.MolFromSmiles(s)) for s in control_smiles if s]
    controls = [(s, m) for s, m in controls if m is not None]
    library = [(s, Chem.MolFromSmiles(s)) for s in library_smiles if isinstance(s, str) and s]
    library = [(s, m) for s, m in library if m is not None]
    if not controls or not library:
        return {"available": False,
                "reason": f"parsed {len(controls)} controls and {len(library)} library molecules"}

    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    ctrl_fps = [gen.GetFingerprint(m) for _, m in controls]
    lib_fps = [gen.GetFingerprint(m) for _, m in library]

    # For each library compound, similarity to its nearest control.
    nearest = [max(DataStructs.BulkTanimotoSimilarity(fp, ctrl_fps)) for fp in lib_fps]
    nearest_arr = np.asarray(nearest, dtype=float)
    sim_cov = float((nearest_arr >= similarity_cutoff).mean())

    ctrl_frames = {f for f in (_framework(m) for _, m in controls) if f}
    lib_frames = [_framework(m) for _, m in library]
    scaffold_cov = (
        float(sum(1 for f in lib_frames if f and f in ctrl_frames) / len(lib_frames))
        if lib_frames else 0.0
    )

    lib_desc = [d for d in (_descriptor_vector(m) for _, m in library) if d]
    ctrl_desc = [(s, d) for (s, _), d in
                 zip(controls, (_descriptor_vector(m) for _, m in controls)) if d]
    desc_report: dict[str, Any] = {}
    inside_counts: dict[str, int] = {}
    if lib_desc and ctrl_desc:
        lo_p, hi_p = percentiles
        bands = {
            k: (float(np.percentile([d[k] for d in lib_desc], lo_p)),
                float(np.percentile([d[k] for d in lib_desc], hi_p)))
            for k in _DESCRIPTORS
        }
        desc_report["library_band"] = {k: [round(v[0], 2), round(v[1], 2)]
                                       for k, v in bands.items()}
        for smi, d in ctrl_desc:
            inside = [k for k in _DESCRIPTORS if bands[k][0] <= d[k] <= bands[k][1]]
            outside = [k for k in _DESCRIPTORS if k not in inside]
            inside_counts[smi] = len(inside)
            desc_report[smi] = {
                "values": {k: round(d[k], 2) for k in _DESCRIPTORS},
                "inside_library_band": inside,
                "outside_library_band": outside,
                "n_inside": len(inside),
                "n_descriptors": len(_DESCRIPTORS),
            }

    return {
        "available": True,
        "n_controls": len(controls),
        "n_library": len(library),
        "similarity_cutoff": similarity_cutoff,
        "similarity_coverage": round(sim_cov, 4),
        "nearest_similarity_max": round(float(nearest_arr.max()), 4),
        "nearest_similarity_median": round(float(np.median(nearest_arr)), 4),
        "nearest_similarity_p95": round(float(np.percentile(nearest_arr, 95)), 4),
        "scaffold_coverage": round(scaffold_cov, 4),
        "n_control_frameworks": len(ctrl_frames),
        "n_library_frameworks": len({f for f in lib_frames if f}),
        "descriptor_percentiles": list(percentiles),
        "descriptors": desc_report,
        "control_descriptors_inside": inside_counts,
    }

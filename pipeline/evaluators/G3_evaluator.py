"""
G3 evaluator — pose validation by direct heavy-atom RMSD.

Specification (frozen 2026-07-30, per ADR-AIIM revision):
  1. Atom map established by MCS BEFORE docking, serialized as JSON, registered.
  2. Same coordinate frame (receptor crystallographic); no Kabsch / superimposition.
  3. Heavy atoms only.
  4. Per-atom-pair distances recorded.
  5. PASS criterion: rank-1 pose RMSD ≤ 2.0 Å.
  6. Best-of-N is diagnostic only; selecting a non-rank-1 pose is forbidden.

Requires RDKit for MCS and RMSD computation.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

from ._control_representativeness import (
    DEFAULT_MIN_SCAFFOLD_COVERAGE,
    DEFAULT_MIN_SIMILARITY_COVERAGE,
    DEFAULT_SIMILARITY_CUTOFF,
    assess,
)


def _registry_artifact_record(registry: Any, artifact_id: str) -> dict[str, Any] | None:
    """The registered record for an artifact_id (timestamp + sha256), or None."""
    if registry is None or artifact_id is None:
        return None
    records: list[dict[str, Any]] | None = None
    if isinstance(registry, list):
        records = registry
    elif hasattr(registry, "get_all_records"):
        try:
            records = registry.get_all_records()
        except Exception:
            return None
    if records is None:
        return None
    for rec in records:
        if rec.get("record_type") == "artifact" and str(rec.get("record_id")) == str(artifact_id):
            return rec
    return None


def _sha256(path: Path) -> str:
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest().upper()


def evaluate(
    contract: dict[str, Any],
    inputs: dict[str, Any],
    registry: Any,
) -> dict[str, Any]:
    passed, failed, warned, evidence = [], [], [], {}
    critical_failures = []

    atom_map_4eiy = Path(inputs.get("atom_map_4eiy_json", ""))
    # 3EML is excluded from docking (ADR-AIIM-036); these stay defined so the
    # evaluator still works against a contract that declares them.
    atom_map_3eml = Path(inputs.get("atom_map_3eml_json", ""))
    poses_4eiy = Path(inputs.get("redocked_poses_4eiy", ""))
    poses_3eml = Path(inputs.get("redocked_poses_3eml", ""))
    ref_4eiy = Path(inputs.get("reference_ligand_4eiy_sdf", ""))
    ref_3eml = Path(inputs.get("reference_ligand_3eml_sdf", ""))

    # An evaluator that runs checks the contract does not declare is not
    # contract-driven: the contract is supposed to be the single statement of
    # what this gate verifies.
    declared = {c["id"] for c in contract.get("checks", [])}
    want_3eml = "rmsd_rank1_3eml" in declared

    # ── Consolidated atom map pre-registration check ──────────────────────────
    cid = "atom_map_pre_registration"
    missing = []
    if cid in declared:
        if not atom_map_4eiy.exists():
            missing.append("atom_map_4eiy_json")
        if want_3eml and not atom_map_3eml.exists():
            missing.append("atom_map_3eml_json")
    if cid not in declared:
        pass                      # not declared by this contract
    elif missing:
        failed.append(cid)
        critical_failures.append(cid)
        evidence[cid] = {"error": f"Atom map files not found: {missing}"}
    else:
        # Validate JSON structure
        map_ok = True
        maps = [(atom_map_4eiy, "4EIY")] + ([(atom_map_3eml, "3EML")] if want_3eml else [])
        for mp, label in maps:
            with open(mp) as f:
                m = json.load(f)
            if not isinstance(m, dict) or "atom_map" not in m or "method" not in m:
                map_ok = False
                evidence[cid] = {
                    "error": f"{label} atom map JSON missing 'atom_map' or 'method' key"
                }
            elif m.get("method") != "rdkit_MCS":
                map_ok = False
                evidence[cid] = {
                    "error": f"{label} atom map method is {m.get('method')!r}; expected 'rdkit_MCS'"
                }
        if map_ok:
            # Held-out finding H03 (2026-08-08): file existence + JSON schema
            # validity say nothing about *when* the map was written relative to
            # the docking run, or whether it was silently rewritten after being
            # registered -- a post-hoc remap that keeps the same schema passed
            # silently. The contract's own precondition
            # (atom_map_..._registered_before_docking) promises the contract
            # engine verifies registration order; nothing anywhere in
            # kdss_audit/ actually did. This closes the gap with two
            # complementary registry-backed signals: (1) hash match -- was the
            # file modified at all since it was registered? A post-hoc rewrite
            # changes the file's content and therefore its SHA-256, even if the
            # rewrite itself is never re-registered, which is exactly the H03
            # scenario. (2) Timestamp ordering -- if the map genuinely was
            # re-registered, was that registration still before docking?
            order_ok = True
            order_detail: dict[str, Any] = {}
            targets = [("4eiy", atom_map_4eiy)] + ([("3eml", atom_map_3eml)] if want_3eml else [])
            for label, mp in targets:
                amap_artifact_id = inputs.get(f"atom_map_artifact_id_{label}", f"G3::atom_map_{label}")
                pose_artifact_id = inputs.get(f"docking_artifact_id_{label}", f"G3::redocked_poses_{label}")
                amap_rec = _registry_artifact_record(registry, amap_artifact_id)
                pose_rec = _registry_artifact_record(registry, pose_artifact_id)
                entry: dict[str, Any] = {"atom_map_artifact_id": amap_artifact_id,
                                          "docking_artifact_id": pose_artifact_id}
                if amap_rec is None or pose_rec is None:
                    entry["checked"] = False
                    entry["reason"] = "one or both artifacts not found in registry"
                else:
                    entry["checked"] = True
                    registered_hash = str(amap_rec.get("sha256", "")).upper()
                    actual_hash = _sha256(mp)
                    entry["registered_sha256"] = registered_hash
                    entry["actual_sha256"] = actual_hash
                    entry["hash_match"] = actual_hash == registered_hash
                    amap_ts, pose_ts = amap_rec.get("timestamp"), pose_rec.get("timestamp")
                    entry["atom_map_registered_at"] = amap_ts
                    entry["docking_registered_at"] = pose_ts
                    entry["order_ok"] = amap_ts is not None and pose_ts is not None and amap_ts < pose_ts
                    if not entry["hash_match"]:
                        entry["error"] = ("atom map file content no longer matches what was "
                                           "registered -- modified after registration")
                        order_ok = False
                    elif not entry["order_ok"]:
                        entry["error"] = "atom map registered at or after the docking run"
                        order_ok = False
                order_detail[label] = entry

            if order_ok:
                passed.append(cid)
                evidence[cid] = {"4EIY_map_valid": True,
                                 "3EML_map_valid": True if want_3eml else "not_applicable",
                                 "registration_order": order_detail}
            else:
                failed.append(cid)
                critical_failures.append(cid)
                evidence[cid] = {"error": "atom map was not registered before the docking run",
                                  "registration_order": order_detail}

    # ── Coordinate frame checks ───────────────────────────────────────────────
    frame_targets = [t for t in [("coordinate_frame_4eiy", "4EIY", poses_4eiy)]
                     if t[0] in declared]
    if want_3eml:
        frame_targets.append(("coordinate_frame_3eml", "3EML", poses_3eml))
    for cid, receptor_label, poses_path in frame_targets:
        # Verifying coordinate frame requires reading receptor CA atoms from both
        # the reference ligand's associated receptor file and the docked pose file.
        # Delegated to the domain tool; we check that the tool produced a frame_check artifact.
        frame_key = f"frame_check_{receptor_label.lower()}"
        frame_result = inputs.get(frame_key)
        if frame_result is None:
            warned.append(cid)
            evidence[cid] = {"warning": f"frame_check_{receptor_label.lower()} not in inputs; cannot verify"}
        elif frame_result is True or frame_result == "PASS":
            passed.append(cid)
            evidence[cid] = {"frame_match": True}
        else:
            failed.append(cid)
            critical_failures.append(cid)
            evidence[cid] = {"frame_match": False, "detail": str(frame_result)}

    # ── Direct RMSD computation ───────────────────────────────────────────────
    pass_threshold = 2.0
    rmsd_targets = [t for t in
                    [("rmsd_rank1_4eiy", ref_4eiy, poses_4eiy, atom_map_4eiy, "4EIY")]
                    if t[0] in declared]
    if want_3eml:
        rmsd_targets.append(("rmsd_rank1_3eml", ref_3eml, poses_3eml, atom_map_3eml, "3EML"))
    for cid, ref_sdf, poses_sdf, atom_map_json, label in rmsd_targets:
        result = _compute_direct_rmsd(ref_sdf, poses_sdf, atom_map_json, label)
        evidence[cid] = result
        if result.get("error"):
            failed.append(cid)
            critical_failures.append(cid)
        else:
            rank1_rmsd = result.get("rank1_rmsd")
            ok = rank1_rmsd is not None and rank1_rmsd <= pass_threshold
            (passed if ok else failed).append(cid)
            if not ok:
                critical_failures.append(cid)

    # ── Per-atom distances recorded check ────────────────────────────────────
    cid = "per_atom_distances_recorded"
    for label in ([] if cid not in declared else (["4EIY", "3EML"] if want_3eml else ["4EIY"])):
        key = f"per_atom_distances_{label.lower()}"
        if key in inputs and inputs[key]:
            passed.append(f"{cid}_{label}")
            evidence[f"{cid}_{label}"] = {"recorded": True}
        else:
            warned.append(f"{cid}_{label}")
            evidence[f"{cid}_{label}"] = {"recorded": False}

    # ── Best-of-N diagnostic ─────────────────────────────────────────────────
    cid = "best_of_n_diagnostic"
    for label in ([] if cid not in declared else (["4EIY", "3EML"] if want_3eml else ["4EIY"])):
        key = f"all_pose_rmsds_{label.lower()}"
        if key in inputs:
            evidence[f"best_of_n_{label}"] = {
                "all_rmsds": inputs[key],
                "best_rmsd": min(inputs[key]) if inputs[key] else None,
                "note": "diagnostic only — rank-1 is the sole PASS criterion",
            }
    if cid in declared:
        passed.append(cid)  # This check always passes (it only records)

    # ── Reproducibility check for 4EIY ───────────────────────────────────────
    cid = "reproducibility_check_4eiy"
    rmsd_4eiy_evidence = evidence.get("rmsd_rank1_4eiy", {})
    rank1_rmsd_4eiy = rmsd_4eiy_evidence.get("rank1_rmsd")
    if cid not in declared:
        pass
    elif rank1_rmsd_4eiy is not None:
        expected = 1.909
        band = [1.80, 2.00]
        in_band = band[0] <= rank1_rmsd_4eiy <= band[1]
        evidence[cid] = {
            "rank1_rmsd": rank1_rmsd_4eiy,
            "legacy_reference": expected,
            "tolerance_band": band,
            "in_band": in_band,
        }
        if in_band:
            passed.append(cid)
        else:
            warned.append(cid)
            evidence[cid]["note"] = "Outside reproducibility band — investigate before accepting"
    elif cid in declared:
        warned.append(cid)
        evidence[cid] = {"note": "Rank-1 RMSD for 4EIY not available"}

    # ── Validation representativeness ────────────────────────────────────────
    #
    # Everything above establishes that the protocol reproduces the *control
    # ligand's* pose. That is not the same as establishing it is valid for the
    # compounds about to be ranked, and where pocket hydration or induced fit is
    # ligand-dependent the two come apart without any downstream symptom. This
    # check makes the gap visible instead of leaving it assumed.
    cid = "control_representativeness"
    if any(c["id"] == cid for c in contract.get("checks", [])):
        check = next(c for c in contract["checks"] if c["id"] == cid)
        sev = check.get("severity", "high")
        controls = inputs.get("control_ligand_smiles") or []
        if isinstance(controls, str):
            controls = [controls]
        library = inputs.get("screening_library_smiles") or []

        if not controls or not library:
            warned.append(cid)
            evidence[cid] = {"warning": "control or library SMILES not supplied; "
                                        "representativeness not assessed"}
        else:
            report = assess(
                controls, library,
                similarity_cutoff=float(check.get("similarity_cutoff",
                                                  DEFAULT_SIMILARITY_CUTOFF)),
            )
            if not report.get("available"):
                warned.append(cid)
                evidence[cid] = {"warning": f"not assessable: {report.get('reason')}"}
            else:
                min_sim = float(check.get("min_similarity_coverage",
                                          DEFAULT_MIN_SIMILARITY_COVERAGE))
                min_scaf = float(check.get("min_scaffold_coverage",
                                           DEFAULT_MIN_SCAFFOLD_COVERAGE))
                sim_ok = report["similarity_coverage"] >= min_sim
                scaf_ok = report["scaffold_coverage"] >= min_scaf
                report.update({
                    "min_similarity_coverage": min_sim,
                    "min_scaffold_coverage": min_scaf,
                    "similarity_ok": sim_ok,
                    "scaffold_ok": scaf_ok,
                })
                if sim_ok or scaf_ok:
                    passed.append(cid)
                else:
                    # Deliberately not a HALT: the science is sound for the
                    # control, and whether the gap matters is a human call. It
                    # must be made explicitly and on the record.
                    report["interpretation"] = (
                        "The redocking control is chemically unrelated to the "
                        "screening library, so a passing pose-validation result "
                        "does not demonstrate the protocol transfers to it. "
                        "Ranking the library on this basis requires explicit "
                        "human justification."
                    )
                    failed.append(cid)
                    if sev == "critical":
                        critical_failures.append(cid)
                evidence[cid] = report

    # ── Retrospective structure-based enrichment ─────────────────────────────
    #
    # Direct evidence that this protocol separates actives from inactives,
    # independent of whether the redocking control resembles the library.
    # Unlike control_representativeness (a similarity proxy), this cannot be
    # defeated by a structurally novel library — it measures the protocol's
    # actual behaviour on labeled compounds.
    cid = "retrospective_enrichment"
    if any(c["id"] == cid for c in contract.get("checks", [])):
        check = next(c for c in contract["checks"] if c["id"] == cid)
        sev = check.get("severity", "critical")
        re_input = inputs.get("retrospective_enrichment")
        if not re_input:
            warned.append(cid)
            evidence[cid] = {"warning": "no retrospective-enrichment evidence supplied"}
        else:
            auc_min = float(check.get("threshold_auc_min", 0.60))
            bedroc_min = float(check.get("threshold_bedroc20_min", 0.40))
            auc = re_input.get("auc_roc")
            bedroc = re_input.get("bedroc20")
            protocol_match = re_input.get("protocol_matches_screen", True)
            detail: dict[str, Any] = {
                "auc_roc": auc, "bedroc20": bedroc,
                "threshold_auc_min": auc_min, "threshold_bedroc20_min": bedroc_min,
                "protocol_matches_screen": protocol_match,
            }
            if not protocol_match:
                # Evidence generated under a different protocol says nothing
                # about the one that will actually screen the library.
                detail["error"] = ("benchmark was not docked under the identical protocol "
                                    "about to be used for the screen -- not valid evidence "
                                    "of transferability")
                failed.append(cid)
                if sev == "critical":
                    critical_failures.append(cid)
            elif auc is None:
                warned.append(cid)
                detail["warning"] = "auc_roc not provided"
            else:
                ok = auc >= auc_min
                detail["auc_pass"] = ok
                if bedroc is not None:
                    detail["bedroc20_pass"] = bedroc >= bedroc_min
                (passed if ok else failed).append(cid)
                if not ok and sev == "critical":
                    critical_failures.append(cid)
            evidence[cid] = detail

    # ── Structure-based evidence can supersede a failed similarity proxy ──────
    #
    # control_representativeness only shows whether the control resembles the
    # library; it does not show whether the protocol separates actives from
    # inactives for a library the control does not resemble. When the proxy
    # fails, passing retrospective_enrichment evidence closes that gap
    # directly. Contracts that do not declare retrospective_enrichment behave
    # exactly as before this amendment (ADR-AIIM-073).
    if "control_representativeness" in failed:
        re_cid = "retrospective_enrichment"
        re_declared = any(c["id"] == re_cid for c in contract.get("checks", []))
        if re_cid in passed:
            evidence.setdefault(re_cid, {})["supersedes_control_representativeness"] = True
            evidence.setdefault("control_representativeness", {})["superseded_by"] = re_cid
        elif re_declared:
            # retrospective_enrichment was required (declared by this contract)
            # but is missing or itself failed: the similarity proxy failed and
            # nothing structure-based corroborates transferability. That
            # combination is a critical gap, not an open ESCALATE.
            if re_cid not in failed:
                failed.append(re_cid)
            if re_cid not in critical_failures:
                critical_failures.append(re_cid)
            evidence.setdefault(re_cid, {})["error"] = (
                "control_representativeness failed and no passing retrospective_enrichment "
                "evidence corroborates transferability to this library"
            )

    # ── Informational checks (not gating) ────────────────────────────────────
    for cid in [c for c in ["pocket_centroid_check", "clash_count",
                            "failed_pose_accounting"] if c in declared]:
        info = inputs.get(cid)
        evidence[cid] = {"value": info, "role": "informational_only"}
        passed.append(cid)

    # Critical breaches halt; non-critical failures still require a human
    # decision. Returning PASS while `failed` is non-empty would let a
    # high-severity finding through silently.
    if critical_failures:
        state = "HALT"
    elif failed:
        state = "ESCALATE"
    else:
        state = "PASS"
    return {
        "state": state,
        "passed": passed,
        "failed": failed,
        "warned": warned,
        "critical_failures": critical_failures,
        "evidence": evidence,
        "summary": (
            f"{len(passed)} checks passed/informational, {len(failed)} failed "
            f"({len(critical_failures)} critical), {len(warned)} warned"
        ),
    }


def template_smiles_for(atom_map_data: dict[str, Any]) -> str | None:
    """
    Bond-order template for the reference chemistry.

    Read from the pre-registered atom-map artifact so the template is fixed
    before docking rather than chosen when the RMSD is computed.
    """
    return (atom_map_data.get("template_smiles")
            or atom_map_data.get("reference_ligand_smiles"))


def _compute_direct_rmsd(
    ref_sdf: Path,
    poses_sdf: Path,
    atom_map_json: Path,
    label: str,
) -> dict[str, Any]:
    """
    Direct heavy-atom RMSD of rank-1 docked pose vs. reference ligand.
    Uses the pre-registered atom map (ref_idx → pose_idx mapping).
    No Kabsch rotation; atoms must already be in the same coordinate frame.
    """
    if not ref_sdf.exists():
        return {"error": f"Reference SDF not found: {ref_sdf}"}
    if not poses_sdf.exists():
        return {"error": f"Poses SDF not found: {poses_sdf}"}
    if not atom_map_json.exists():
        return {"error": f"Atom map JSON not found: {atom_map_json}"}

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, rdMolAlign
    except ImportError:
        return {"error": "RDKit not installed"}

    with open(atom_map_json) as f:
        atom_map_data = json.load(f)
    atom_map: dict[str, int] = atom_map_data.get("atom_map", {})

    # Crystal ligands extracted from a PDB carry no bond orders, so each reader
    # perceives their chemistry differently. Applying the pre-registered
    # template makes the reference and the poses the same compound; without it
    # any graph-based comparison either fails outright or silently compares
    # different molecules.
    template_smiles = template_smiles_for(atom_map_data)
    template = Chem.MolFromSmiles(template_smiles) if template_smiles else None

    def _load(path: Path) -> list[Any]:
        mols = [m for m in Chem.SDMolSupplier(str(path), removeHs=False, sanitize=False)
                if m is not None]
        out = []
        for m in mols:
            try:
                m = Chem.RemoveHs(m, sanitize=False)
                if template is not None:
                    m = AllChem.AssignBondOrdersFromTemplate(template, m)
                out.append(m)
            except Exception:
                continue
        return out

    ref_list = _load(ref_sdf)
    ref_mol = ref_list[0] if ref_list else None
    if ref_mol is None:
        return {"error": f"Could not parse reference SDF: {ref_sdf}"}

    all_poses = _load(poses_sdf)
    if not all_poses:
        return {"error": f"No valid poses in {poses_sdf}"}

    # Symmetry-corrected RMSD in place. CalcRMS enumerates the symmetry-
    # equivalent atom correspondences and takes the minimum, and — unlike
    # GetBestRMS — performs no superposition, which is what the contract
    # requires. A single MCS submatch cannot do this: on a substructure with
    # 2-fold symmetry it resolves the correspondence arbitrarily and can report
    # >10 A for a pose sitting 4 A from the reference.
    if template is not None:
        sym_rmsds: list[float] = []
        for pose_mol in all_poses:
            try:
                sym_rmsds.append(round(float(rdMolAlign.CalcRMS(pose_mol, ref_mol)), 4))
            except Exception as exc:
                return {"error": f"symmetry-corrected RMSD failed for {label}: {exc}"}
        return {
            "rank1_rmsd": sym_rmsds[0],
            "n_poses": len(sym_rmsds),
            "all_pose_rmsds": sym_rmsds,
            "best_rmsd": min(sym_rmsds),
            "method": "symmetry-corrected in-place RMSD (rdMolAlign.CalcRMS), "
                      "all heavy atoms, no superposition",
            "n_heavy_atoms": ref_mol.GetNumAtoms(),
            "template_smiles": template_smiles,
            "mcs_atom_map_n": len(atom_map),
            "mcs_atom_map_role": "diagnostic only; not the RMSD basis",
        }

    ref_conf = ref_mol.GetConformer()
    per_atom_distances: list[float] = []
    all_pose_rmsds: list[float] = []

    # An index that does not exist in the molecule means the map was built
    # against a different parse or a different atom ordering, so it is not a
    # valid correspondence. Report that as a finding; letting RDKit raise would
    # surface it as an evaluator crash and lose the reason.
    n_ref, n_pose = ref_mol.GetNumAtoms(), all_poses[0].GetNumAtoms()
    bad = [
        (r, p) for r, p in atom_map.items()
        if int(r) >= n_ref or int(p) >= n_pose
    ]
    if bad:
        return {
            "error": (
                f"atom map indices out of range for {label}: reference has "
                f"{n_ref} atoms, pose has {n_pose}; offending pairs {bad[:10]}. "
                "The map was built against a different parse or atom ordering."
            ),
            "n_mapped_atoms": len(atom_map),
        }

    for pose_mol in all_poses:
        pose_conf = pose_mol.GetConformer()
        sq_dists = []
        for ref_idx_str, pose_idx in atom_map.items():
            ref_idx = int(ref_idx_str)
            ref_pos = ref_conf.GetAtomPosition(ref_idx)
            pose_pos = pose_conf.GetAtomPosition(pose_idx)
            d2 = (ref_pos.x - pose_pos.x)**2 + (ref_pos.y - pose_pos.y)**2 + (ref_pos.z - pose_pos.z)**2
            sq_dists.append(d2)
        rmsd = math.sqrt(sum(sq_dists) / len(sq_dists)) if sq_dists else float("inf")
        all_pose_rmsds.append(rmsd)
        if len(all_pose_rmsds) == 1:
            # Rank-1 pose: record per-atom distances
            per_atom_distances = [math.sqrt(d) for d in sq_dists]

    rank1_rmsd = all_pose_rmsds[0] if all_pose_rmsds else None

    return {
        "rank1_rmsd": round(rank1_rmsd, 4) if rank1_rmsd is not None else None,
        "n_poses": len(all_pose_rmsds),
        "all_pose_rmsds": [round(r, 4) for r in all_pose_rmsds],
        "best_rmsd": round(min(all_pose_rmsds), 4) if all_pose_rmsds else None,
        "per_atom_distances_rank1": [round(d, 4) for d in per_atom_distances],
        "atom_pairs_used": len(atom_map),
        "label": label,
    }


def build_atom_map(ref_sdf: Path, output_json: Path) -> dict[str, Any]:
    """
    Pre-docking utility: compute MCS between reference ligand and a template
    pose (typically the crystal ligand itself), then serialize the atom map.

    Call this BEFORE any docking run. Register the output JSON in the artifact
    registry with a timestamp before the docking job starts.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import rdFMCS
    except ImportError:
        return {"error": "RDKit not installed"}

    supplier = Chem.SDMolSupplier(str(ref_sdf), removeHs=True)
    mols = [m for m in supplier if m is not None]
    if len(mols) < 1:
        return {"error": "No valid molecules in reference SDF"}

    ref_mol = mols[0]

    # MCS of the molecule with itself → complete atom map (identity)
    mcs = rdFMCS.FindMCS(
        [ref_mol, ref_mol],
        ringMatchesRingOnly=True,
        completeRingsOnly=True,
    )
    query = Chem.MolFromSmarts(mcs.smartsString)
    matches = ref_mol.GetSubstructMatches(query)
    if not matches:
        return {"error": "MCS self-match returned no matches"}

    atom_map = {str(ref_idx): pose_idx for ref_idx, pose_idx in zip(matches[0], matches[0])}
    data = {
        "method": "rdkit_MCS",
        "mcs_smarts": mcs.smartsString,
        "n_atoms": len(atom_map),
        "atom_map": atom_map,
        "source_sdf": str(ref_sdf),
        "note": (
            "Identity map for crystal redocking. "
            "For cross-docking, replace ref_mol with the docking query mol."
        ),
    }
    output_json.parent.mkdir(parents=True, exist_ok=True)
    with open(output_json, "w") as f:
        json.dump(data, f, indent=2)
    return {"success": True, "n_atoms": len(atom_map), "output": str(output_json)}

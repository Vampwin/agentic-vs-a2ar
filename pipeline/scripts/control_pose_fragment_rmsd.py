"""control_pose_fragment_rmsd.py -- m11. Where the control's 1.92 A sits.

Figure 5A reports that DiffDock reproduced the crystallographic ZM241385 pose
to 1.92 A, inside the declared 2.0 A gate by 0.08 A, and Section 3.6 rests a
conclusion on that gate having passed. A single whole-molecule RMSD does not
say whether the deviation is spread evenly or concentrated in one flexible
substituent -- and the two readings license different amounts of confidence.
The rendered panel shows the furanyl end pointing differently between the two
poses, which is the visual reason a reviewer would ask.

So: split the ligand at its rotatable bonds into a rigid core and the flexible
tail, and report each. RMSD is in place and symmetry-corrected, matching the
convention used everywhere else in this project -- never superposed, because a
superposed value would flatter a pose that is correctly shaped and wrongly
placed.

Usage:
    python control_pose_fragment_rmsd.py [--out DIR]
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolAlign

RDLogger.DisableLog("rdApp.*")


def _root() -> Path:
    marker = "00_project_control/08_DECISION_LOG.md"
    for cand in (Path(r"D:\Agentic VS Workflow"), Path("/mnt/d/Agentic VS Workflow")):
        if (cand / marker).is_file():
            return cand
    raise SystemExit("project root not found")


ROOT = _root()
RUNS = ROOT / "20_new_runs"
CRYSTAL = RUNS / "snp_screen_2026-08-04/G2_preparation/ZMA_reference_4EIY.sdf"
DIFFDOCK = RUNS / "diffdock_orthogonal_2026-08-23/results/ZM241385_control/rank1.sdf"
VINA = RUNS / "snp_screen_2026-08-04/G3_redocking/ZMA_redocked_4EIY.sdf"

# Three fragments, because the visual objection is specifically about the
# furanyl end. The two substituents are matched by SMARTS on the reference,
# so the split is defined by chemistry rather than by atom index order --
# which differs between the crystal file and the docked poses. The fused
# bicyclic core is then whatever neither substituent claims: writing SMARTS
# for a triazolo[1,5-a]triazine is easy to get subtly wrong, and a silent
# zero-atom match is exactly the failure this project keeps finding.
TAIL_SMARTS = "NCCc1ccc(O)cc1"   # 2-(4-hydroxyphenyl)ethylamino tail
FURAN_SMARTS = "c1ccoc1"         # the furanyl substituent


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


# The crystal SDF carries no bond orders -- RDKit reads it as a saturated
# ring system, not as ZM241385 -- so it must be templated before anything is
# compared. This is the same correction diffdock_orthogonal_check.py applies
# before computing the 1.92 A the manuscript quotes, reused rather than
# reinvented so the fragment split cannot disagree with the whole-molecule
# number it decomposes.
ZMA_TEMPLATE = "Nc1nc(NCCc2ccc(O)cc2)nc2nc(-c3ccco3)nn12"


def load(p: Path, template: bool = False) -> Chem.Mol:
    m = Chem.MolFromMolFile(str(p), removeHs=True)
    if m is None:
        raise SystemExit(f"could not read {p}")
    if template:
        m = AllChem.AssignBondOrdersFromTemplate(
            Chem.MolFromSmiles(ZMA_TEMPLATE), m)
    return m


def subset_rmsd(pose: Chem.Mol, ref: Chem.Mol, ref_idx: list[int],
                amap: tuple[int, ...]) -> float:
    """In-place RMSD over a subset of REFERENCE atoms, no superposition.

    amap[i] is the pose atom matching reference atom i. Atom order differs
    between the crystal file and the docked poses, so the pairing cannot be
    by index."""
    cp, cr = pose.GetConformer(), ref.GetConformer()
    d = [np.linalg.norm(np.array(cp.GetAtomPosition(amap[i])) -
                        np.array(cr.GetAtomPosition(i))) for i in ref_idx]
    return float(np.sqrt(np.mean(np.square(d))))


def best_map(pose: Chem.Mol, ref: Chem.Mol) -> tuple[int, ...]:
    """The reference-to-pose atom pairing that minimises in-place RMSD.

    ZM241385's 4-hydroxyphenyl ring is two-fold symmetric, so several
    substructure matches are chemically valid and they do not all describe
    the same geometry. Taking the first one gave a tail RMSD of 1.61 A
    against a whole-molecule 0.38 A -- arithmetically impossible, since ten
    of twenty-five atoms at 1.61 A force the whole to at least 1.02 A. That
    impossibility is the check: the mapping must be the one CalcRMS itself
    would choose, and the assertion below enforces it."""
    cp, cr = pose.GetConformer(), ref.GetConformer()
    best, best_sq = None, None
    for m in pose.GetSubstructMatches(ref, uniquify=False, maxMatches=10000):
        sq = sum(float(np.sum((np.array(cp.GetAtomPosition(m[i])) -
                               np.array(cr.GetAtomPosition(i))) ** 2))
                 for i in range(ref.GetNumAtoms()))
        if best_sq is None or sq < best_sq:
            best, best_sq = m, sq
    if best is None:
        raise SystemExit("no substructure match between pose and reference")
    rms = float(np.sqrt(best_sq / ref.GetNumAtoms()))
    ref_rms = float(rdMolAlign.CalcRMS(pose, ref))
    if abs(rms - ref_rms) > 5e-3:
        raise SystemExit(
            f"chosen mapping gives {rms:.4f} A but CalcRMS gives "
            f"{ref_rms:.4f} A -- the fragment split would not decompose the "
            f"number it claims to decompose")
    return best


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(RUNS / "control_pose_fragments_2026-08-28"))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    ref = load(CRYSTAL, template=True)
    report = {"reference": {"path": str(CRYSTAL.relative_to(ROOT).as_posix()),
                            "sha256": sha256(CRYSTAL),
                            "n_heavy_atoms": ref.GetNumHeavyAtoms()},
              "note": "In-place, symmetry-corrected where the whole molecule is "
                      "compared; per-fragment values are in-place over the matched "
                      "atom subset. No superposition anywhere.",
              "poses": {}}

    tail = list(ref.GetSubstructMatch(Chem.MolFromSmarts(TAIL_SMARTS)))
    furan = list(ref.GetSubstructMatch(Chem.MolFromSmarts(FURAN_SMARTS)))
    core = [i for i in range(ref.GetNumAtoms())
            if i not in set(tail) | set(furan)]
    if not tail or not furan:
        raise SystemExit("a substituent SMARTS matched nothing -- refusing to "
                         "report a fragment split that is not the molecule")
    report["fragments"] = {
        "tail_smarts": TAIL_SMARTS, "n_tail_atoms": len(tail),
        "furan_smarts": FURAN_SMARTS, "n_furan_atoms": len(furan),
        "n_core_atoms": len(core),
        "core_definition": "every heavy atom claimed by neither substituent",
        "total_check": len(tail) + len(furan) + len(core)}

    for name, path in (("vina", VINA), ("diffdock", DIFFDOCK)):
        mol = load(path)
        entry = {"path": str(path.relative_to(ROOT).as_posix()),
                 "sha256": sha256(path),
                 "whole_molecule_rmsd": round(
                     float(rdMolAlign.CalcRMS(mol, ref)), 3)}
        amap = best_map(mol, ref)
        if len(amap) == ref.GetNumAtoms():
            entry["core_rmsd"] = round(subset_rmsd(mol, ref, core, amap), 3)
            entry["furan_rmsd"] = round(subset_rmsd(mol, ref, furan, amap), 3)
            entry["tail_rmsd"] = round(subset_rmsd(mol, ref, tail, amap), 3)
        report["poses"][name] = entry

    (out / "control_pose_fragments.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")

    print(f"reference: {CRYSTAL.name}  ({ref.GetNumHeavyAtoms()} heavy atoms)")
    print(f"core {len(core)}  furan {len(furan)}  tail {len(tail)}  "
          f"= {len(core) + len(furan) + len(tail)} of {ref.GetNumAtoms()}")
    for name, e in report["poses"].items():
        extra = (f"   core {e['core_rmsd']:.2f}   furan {e['furan_rmsd']:.2f}"
                 f"   tail {e['tail_rmsd']:.2f}"
                 if "core_rmsd" in e else "   (fragment split unavailable)")
        print(f"  {name:>9}: whole {e['whole_molecule_rmsd']:.2f} A{extra}")
    print(f"\nwritten -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

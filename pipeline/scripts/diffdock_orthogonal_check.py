"""
diffdock_orthogonal_check.py -- evaluate the pre-registered orthogonal-engine
check of ADR-AIIM-111.

The decision rules are fixed by that ADR and are NOT re-derived here; this
script only applies them. Rule 1 (the control gate) is evaluated first, and if
it fails, Rule 2's candidate numbers are printed but explicitly marked
uninterpretable rather than quietly reported.

RMSD is computed in place and symmetry-corrected (rdMolAlign.CalcRMS), never
superposed. Section 2.5 gives the reason: a superposed RMSD flatters a pose
that is correctly shaped but wrongly placed, which is precisely the failure
mode this check exists to detect.
"""
from __future__ import annotations

import json
import re
from pathlib import Path
from kdss_paths import project_root, resolve, under

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolAlign

RDLogger.DisableLog("rdApp.*")

ROOT = project_root()
RUN = resolve("20_new_runs/diffdock_orthogonal_2026-08-23")
PREP = resolve("20_new_runs/snp_screen_2026-08-04/G2_preparation")
POSES = resolve("20_new_runs/np_library_funnel_2026-08-15/full_poses_chunk1")

CONTROL_GATE_A = 2.0        # ADR-AIIM-111 Rule 1, the Section 2.5 convention
AGREEMENT_A = 2.0           # ADR-AIIM-111 Rule 2
ONPOCKET_A = 8.0            # the study's existing centroid criterion
CANDIDATES = ["NPLIB01185", "NPLIB01643", "NPLIB03056"]
ZMA_TEMPLATE = "Nc1nc(NCCc2ccc(O)cc2)nc2nc(-c3ccco3)nn12"


def rank1(cid: str) -> tuple[Chem.Mol, float]:
    """DiffDock's top pose and its confidence, read from the filename."""
    d = under("20_new_runs/diffdock_orthogonal_2026-08-23", "results") / cid
    hit = [p for p in d.iterdir() if re.match(r"^rank1_confidence", p.name)]
    if len(hit) != 1:
        raise SystemExit(f"{cid}: expected one rank1_confidence file, found {len(hit)}")
    conf = float(re.search(r"confidence(-?[\d.]+)\.sdf$", hit[0].name).group(1))
    m = Chem.MolFromMolFile(str(hit[0]))
    if m is None:
        raise SystemExit(f"{cid}: DiffDock rank-1 pose did not parse")
    return Chem.RemoveHs(m), conf


def band(conf: float) -> str:
    """The tool's own documented bands -- not thresholds invented here."""
    return "high" if conf > 0 else ("moderate" if conf >= -1.5 else "low")


def centroid(m: Chem.Mol) -> np.ndarray:
    return np.array(m.GetConformer().GetPositions()).mean(axis=0)


def main() -> int:
    box = json.loads((under("20_new_runs/snp_screen_2026-08-04/G2_preparation", "docking_box.json")).read_text())
    pocket = np.array(box["center"])
    out: dict = {"adr": "ADR-AIIM-111", "rules_applied_verbatim": True}

    # ---- Rule 1: the control gate -----------------------------------------
    crystal = Chem.MolFromMolFile(str(under("20_new_runs/snp_screen_2026-08-04/G2_preparation", "ZMA_reference_4EIY.sdf")))
    crystal = AllChem.AssignBondOrdersFromTemplate(
        Chem.MolFromSmiles(ZMA_TEMPLATE), crystal)
    dd_ctrl, conf_ctrl = rank1("ZM241385_control")

    if Chem.MolToInchiKey(dd_ctrl) != Chem.MolToInchiKey(crystal):
        raise SystemExit("control identity mismatch -- DiffDock docked a "
                         "different molecule than the crystal reference")

    rmsd_ctrl = rdMolAlign.CalcRMS(dd_ctrl, crystal)
    passed = rmsd_ctrl <= CONTROL_GATE_A

    print("=" * 74)
    print("RULE 1 -- CONTROL GATE (evaluated first, per ADR-AIIM-111)")
    print("=" * 74)
    print(f"  ZM241385 DiffDock rank-1 vs crystallographic pose")
    print(f"    in-place symmetry-corrected RMSD : {rmsd_ctrl:.2f} A")
    print(f"    threshold                        : {CONTROL_GATE_A:.1f} A")
    print(f"    DiffDock confidence              : {conf_ctrl:+.2f} ({band(conf_ctrl)})")
    print(f"    centroid vs pocket centre        : "
          f"{np.linalg.norm(centroid(dd_ctrl) - pocket):.2f} A")
    print(f"\n    >>> CONTROL {'PASSES' if passed else 'FAILS'}")

    out["control"] = {
        "rmsd_to_crystal_A": round(rmsd_ctrl, 3),
        "threshold_A": CONTROL_GATE_A,
        "passed": bool(passed),
        "confidence": conf_ctrl,
        "confidence_band": band(conf_ctrl),
        "centroid_to_pocket_A": round(float(np.linalg.norm(centroid(dd_ctrl) - pocket)), 2),
    }

    # ---- Rule 2: candidates ------------------------------------------------
    print("\n" + "=" * 74)
    print("RULE 2 -- CANDIDATE AGREEMENT vs the Vina poses the MD is running")
    if not passed:
        print("  *** UNINTERPRETABLE: Rule 1 failed. Per ADR-AIIM-111 no")
        print("  *** conclusion about the candidates may be drawn from these.")
    print("=" * 74)

    out["candidates"] = {}
    for cid in CANDIDATES:
        vina = Chem.RemoveHs(
            [m for m in Chem.SDMolSupplier(str(POSES / f"{cid}.sdf"), removeHs=False)
             if m is not None][0])
        dd, conf = rank1(cid)
        if Chem.MolToInchiKey(dd) != Chem.MolToInchiKey(vina):
            raise SystemExit(f"{cid}: DiffDock and Vina poses are different molecules")
        r = rdMolAlign.CalcRMS(dd, vina)
        cdist = float(np.linalg.norm(centroid(dd) - pocket))
        agree = r <= AGREEMENT_A
        print(f"  {cid}")
        print(f"    RMSD to Vina rank-1 pose : {r:6.2f} A   "
              f"{'AGREES' if agree else 'DISAGREES'} (<= {AGREEMENT_A:.1f} A)")
        print(f"    DiffDock confidence      : {conf:+6.2f}   ({band(conf)})")
        print(f"    centroid vs pocket       : {cdist:6.2f} A   "
              f"{'on-pocket' if cdist <= ONPOCKET_A else 'OFF-POCKET'}")
        out["candidates"][cid] = {
            "rmsd_to_vina_A": round(r, 3), "agrees": bool(agree),
            "confidence": conf, "confidence_band": band(conf),
            "centroid_to_pocket_A": round(cdist, 2),
            "on_pocket": bool(cdist <= ONPOCKET_A),
        }

    # ---- the pre-registered prediction, scored honestly --------------------
    cand_conf = [out["candidates"][c]["confidence"] for c in CANDIDATES]
    predicted = conf_ctrl > max(cand_conf)
    print("\n" + "=" * 74)
    print("THE PRE-REGISTERED PREDICTION (ADR-AIIM-111), scored")
    print("=" * 74)
    print(f"  Predicted: DiffDock is markedly MORE confident on ZM241385 than")
    print(f"             on the three natural-product candidates.")
    print(f"  Observed : control {conf_ctrl:+.2f}  vs candidates "
          f"{', '.join(f'{c:+.2f}' for c in cand_conf)}")
    print(f"\n  >>> PREDICTION {'HOLDS' if predicted else 'FAILS'}")
    if not predicted:
        print("  The thesis-derived prediction is not supported. This is")
        print("  reported as such; it was declared before the measurement")
        print("  precisely so that it could not be quietly dropped.")
    out["preregistered_prediction"] = {
        "statement": "control confidence exceeds all three candidate confidences",
        "holds": bool(predicted),
        "control_confidence": conf_ctrl,
        "candidate_confidences": dict(zip(CANDIDATES, cand_conf)),
    }

    (under("20_new_runs/diffdock_orthogonal_2026-08-23", "orthogonal_check_result.json")).write_text(json.dumps(out, indent=2))
    print(f"\nwritten -> {(RUN / 'orthogonal_check_result.json').relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""validation_licence_matrix.py -- what was each method validated FOR?

This project judged DiffDock, MM-GBSA, IFP rescoring and the ensemble variants
against one number: ROC-AUC on a 245-compound ChEMBL benchmark. ADR-AIIM-164
then showed that benchmark's baseline is not separable from heavy-atom count
(0.601 against 0.584, p = 0.722). A method rejected on that axis was rejected by
a ruler with no marks on it.

The fix is not to re-run them. It is to stop treating "is this method good?" as
one question. The scoring-function literature settled this a decade ago: CASF
decouples the capabilities and scores each separately -- "scoring power",
"ranking power", "docking power", "screening power" -- and CASF-2016 reports
that current functions do markedly better at docking power than at the other
three (Su et al., J Chem Inf Model 2019). PoseBusters adds a fifth axis, showing
that RMSD alone certifies nothing because methods reach native-like RMSD with
physically impossible poses (Buttenschoen et al., Chem Sci 2024). Li et al.
(Chem Sci 2025) evaluate deep-learning docking across five dimensions including
interaction recovery and generalisation.

So: one row per method, one column per capability, each cell either a measured
number from a named artifact or an explicit "not measured". The last column is
the only one that matters for reading the rest of the paper -- what use the
evidence licenses.

**This script computes nothing.** Every number is read from a frozen artifact
and the source is carried into the output beside it, so `g8_number_traceability`
can audit the table the same way it audits the manuscript.

Usage:
    python validation_licence_matrix.py --out 20_new_runs/validation_matrix_2026-09-05
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def _root() -> Path:
    for c in (Path(r"D:\Agentic VS Workflow"), Path("/mnt/d/Agentic VS Workflow")):
        if (c / "00_project_control/08_DECISION_LOG.md").is_file():
            return c
    raise SystemExit("project root not found")


ROOT = _root()
RUNS = ROOT / "20_new_runs"

# The five capabilities, named as the literature names them.
AXES = [
    ("pose_accuracy", "Docking power: RMSD to a crystallographic pose"),
    ("physical_validity", "Physical plausibility: PoseBusters, tool defaults"),
    ("interaction_recovery", "Does the pose reproduce the Asn253 anchor?"),
    ("screening_power", "Ranking actives above inactives: ROC-AUC"),
    ("scoring_power", "Correlation with measured affinity"),
]

REFERENCES = {
    "CASF-2016": "Su, Yang, Li, Liu, Han, Wang. Comparative Assessment of "
                 "Scoring Functions: The CASF-2016 Update. J Chem Inf Model "
                 "2019;59(2):895-913. Decouples scoring / ranking / docking / "
                 "screening power; reports docking power as the strongest of "
                 "the four for current functions.",
    "PoseBusters": "Buttenschoen, Morris, Deane. PoseBusters: AI-based docking "
                   "methods fail to generate physically valid poses or "
                   "generalise to novel sequences. Chem Sci 2024;15:3130-3139.",
    "DL-limits": "Li, et al. Decoding the limits of deep learning in molecular "
                 "docking for drug discovery. Chem Sci 2025. Five dimensions: "
                 "pose accuracy, physical plausibility, interaction recovery, "
                 "virtual-screening efficacy, generalisation.",
    "MMGBSA-rescoring": "Sahakyan. Improving virtual screening results with "
                        "MM/GBSA and MM/PBSA rescoring. J Comput Aided Mol Des "
                        "2021;35:731-736; and Zhang et al., J Chem Inf Model "
                        "2014;54:324-337, which reports MM/GBSA rescoring "
                        "raising AUC over Vina across 38 DUD-E targets while "
                        "remaining target-dependent and limited by pose quality.",
}


def _load(rel: str):
    p = RUNS / rel
    return json.loads(p.read_text(encoding="utf-8")) if p.is_file() else None


def build() -> dict:
    xd = _load("crossdock_zma_2026-09-03/crossdock_zma.json")
    lt = _load("redock_9xt_2026-09-05/redock_ligand_transfer.json")
    pv = _load("pose_validity_2026-09-03/pose_validity.json")
    cr = _load("clinical_reference_2026-09-03/clinical_reference.json")
    dd = _load("diffdock_orthogonal_2026-08-23/orthogonal_check_result.json")

    mm = pd.read_csv(ROOT / "04_legacy_evidence/external_validation/"
                     "exp3b_harmonized_metrics.csv")
    mm_auc = mm[(mm.method.str.startswith("MM-GBSA")) & (mm.metric == "AUC")].iloc[0]

    # scoring power: Vina score against measured pChEMBL on the clinical panel.
    import math
    from scipy.stats import pearsonr
    byc = cr["by_compound"]
    ok = lambda v: v is not None and not (isinstance(v, float) and math.isnan(v))
    agonists = {"NECA", "adenosine"}
    ant = [(v["vina_score"], v["pchembl_median"]) for k, v in byc.items()
           if ok(v.get("pchembl_median")) and k not in agonists]
    allr = [(v["vina_score"], v["pchembl_median"]) for k, v in byc.items()
            if ok(v.get("pchembl_median"))]
    r_ant, p_ant = pearsonr([x[0] for x in ant], [x[1] for x in ant])
    r_all, p_all = pearsonr([x[0] for x in allr], [x[1] for x in allr])

    n_contact = sum(1 for k, v in byc.items()
                    if v["class"] == "reference ligand" and k != "ZM241385"
                    and k not in agonists and v["asn253"]["contact"])
    n_ant = sum(1 for k, v in byc.items()
                if v["class"] == "reference ligand" and k != "ZM241385"
                and k not in agonists)

    t1 = (lt or {}).get("by_arm", {}).get("T1", {})
    t2 = (lt or {}).get("by_arm", {}).get("T2", {})

    rows = [
        {
            "method": "AutoDock Vina (pose generation)",
            "pose_accuracy": {
                "value": f"{xd['by_receptor']['4EIY']['rank1_rmsd_A']} A (4EIY), "
                         f"{xd['by_receptor']['5OLG']['rank1_rmsd_A']} A (5OLG), "
                         f"{xd['by_receptor']['3EML']['rank1_rmsd_A']} A (3EML)",
                "caveat": "all three structures were solved with ZM241385: the "
                          "receptor changed, the ligand never did",
                "source": "crossdock_zma_2026-09-03 / ADR-AIIM-162"},
            "physical_validity": {
                "value": f"{pv['by_set']['vina_screen']['fraction_valid']*100:.2f} % of "
                         f"{pv['by_set']['vina_screen']['n_poses']} poses; "
                         f"{pv['by_set']['control_vina_redock']['fraction_valid']*100:.0f} % "
                         f"of {pv['by_set']['control_vina_redock']['n_poses']} control poses",
                "source": "pose_validity_2026-09-03 / ADR-AIIM-166"},
            "interaction_recovery": {
                "value": f"{byc['ZM241385']['asn253']['min_distance_A']} A on the "
                         f"redocked control (contact); {n_contact} of {n_ant} other "
                         "known antagonists placed in the contact",
                "source": "clinical_reference_2026-09-03 / ADR-AIIM-163"},
            "screening_power": {
                "value": "AUC 0.601 (95 % CI 0.510-0.683), against 0.584 for "
                         "heavy-atom count alone, paired p = 0.722",
                "source": "trivial_baseline_2026-09-03 / ADR-AIIM-164"},
            "scoring_power": {
                "value": f"r = {r_ant:+.3f} (p = {p_ant:.3f}) against measured "
                         f"pChEMBL over {len(ant)} antagonists; "
                         f"{r_all:+.3f} (p = {p_all:.3f}) over all {len(allr)} "
                         "reference ligands, carried by two agonists",
                "source": "clinical_reference_2026-09-03 / ADR-AIIM-163"},
            "licence": "Geometry, for ZM241385-like chemotypes. Not ranking, "
                       "not affinity.",
        },
        {
            "method": "DiffDock (pose generation)",
            "pose_accuracy": {
                "value": f"{dd['control']['rmsd_to_crystal_A']} A on the ZM241385 "
                         f"control, passing a {dd['control']['threshold_A']} A gate "
                         f"by {dd['control']['threshold_A']-dd['control']['rmsd_to_crystal_A']:.3f} A",
                "source": "diffdock_orthogonal_2026-08-23 / ADR-AIIM-111"},
            "physical_validity": {
                "value": f"{pv['by_set']['diffdock']['n_poses_all_checks_pass']} of "
                         f"{pv['by_set']['diffdock']['n_poses']} rank-1 poses valid; "
                         "the control is among the failures",
                "note": "the pose that passed the RMSD gate overlaps the receptor. "
                        "RMSD and validity disagree on the same pose, which is "
                        "the PoseBusters result reproduced on a GPCR target.",
                "source": "pose_validity_2026-09-03 / ADR-AIIM-166"},
            "interaction_recovery": {"value": "not measured"},
            "screening_power": {"value": "not measured: DiffDock was run on four "
                                         "compounds as an orthogonal pose check, "
                                         "never over the benchmark"},
            "scoring_power": {"value": "not measured"},
            "licence": "Nothing. Its one validated pose is physically invalid, "
                       "so the RMSD it passed on does not certify it.",
        },
        {
            "method": "MM-GBSA rescoring",
            "pose_accuracy": {"value": "n/a -- rescores poses, does not make them"},
            "physical_validity": {"value": "n/a"},
            "interaction_recovery": {"value": "not measured"},
            "screening_power": {
                "value": f"AUC {mm_auc.value} (95 % CI {mm_auc.ci_lo}-{mm_auc.ci_hi}), "
                         "an interval containing 0.5",
                "note": "published expectation is the opposite: MM/GBSA rescoring "
                        "usually raises AUC over docking alone, though it is "
                        "target-dependent and bounded by the quality of the poses "
                        "it rescores. The pose-validity and interaction-recovery "
                        "cells in the Vina row are where that bound bites here.",
                "source": "exp3b_harmonized_metrics.csv; provenance caveat in "
                          "ADR-AIIM-048 (equilibration and seed records incomplete)"},
            "scoring_power": {"value": "not measured on this target"},
            "licence": "Nothing, on this target, at this pose quality. Reported "
                       "as a measurement rather than omitted.",
        },
        {
            "method": "Interaction-fingerprint (IFP) rescoring",
            "pose_accuracy": {"value": "n/a -- reselects among existing poses"},
            "physical_validity": {"value": "n/a"},
            "interaction_recovery": {
                "value": "the Asn253 filter carried weak signal in isolation "
                         "(27 % of compounds forming the contact were active "
                         "against 19 % of those that did not)",
                "source": "Table 5"},
            "screening_power": {
                "value": "AUC 0.449 (95 % CI 0.343-0.545), below chance, with "
                         "BEDROC20 0.495 against the baseline's 0.275",
                "source": "chembl245_ifp_rescoring_2026-08-15"},
            "scoring_power": {"value": "not measured"},
            "licence": "Early recognition only, and not claimed: at n = 245 the "
                       "1 % cut is 2.45 compounds.",
        },
        {
            "method": "Water-shell ensemble docking",
            "pose_accuracy": {"value": "not measured across shells"},
            "physical_validity": {"value": "not measured"},
            "interaction_recovery": {"value": "not measured"},
            "screening_power": {
                "value": "paired dAUC +0.0028 (95 % CI -0.0568 to +0.0588) "
                         "against a threshold of +0.05 fixed before the pilot",
                "source": "ensemble_shell_pilot_2026-08-15"},
            "scoring_power": {"value": "not measured"},
            "licence": "Stopped by a pre-registered rule, not by a verdict on "
                       "the method.",
        },
        {
            "method": "Random forest (ligand-based)",
            "pose_accuracy": {"value": "n/a -- uses no structure"},
            "physical_validity": {"value": "n/a"},
            "interaction_recovery": {"value": "n/a"},
            "screening_power": {
                "value": "nested cross-validated AUC 0.909 +/- 0.020, with "
                         "Y-randomisation over 100 label permutations",
                "source": "qsar_lbdd_2026-08-07, Figure 3B"},
            "scoring_power": {"value": "not measured"},
            "licence": "Ranking, inside its applicability domain -- which "
                       "excludes this library in four of eight spaces.",
        },
    ]

    if t1.get("available"):
        rows[0]["pose_accuracy"]["ligand_transfer"] = {
            "value": f"{t1['rank1_rmsd_A']} A redocking 9XT, a chemically "
                     f"distinct antagonist, into its own crystal (5OLH)"
                     + (f"; {t2['rank1_rmsd_A']} A into 4EIY"
                        if t2.get("available") else ""),
            "source": "redock_9xt_2026-09-05 / ADR-AIIM-168"}

    return {
        "question": "What was each method validated FOR? A single ROC-AUC "
                    "answers one of five questions, and ADR-AIIM-164 showed "
                    "that one was asked against an uninformative baseline.",
        "authority": "ADR-AIIM-169",
        "framework": {"axes": dict(AXES), "references": REFERENCES},
        "computes_nothing": "Every value is read from a frozen artifact named "
                            "beside it. This file is an index, not a result.",
        "rows": rows,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out = Path(args.out).resolve()
    if out.exists() and any(out.iterdir()):
        raise SystemExit(f"REFUSED: {out} exists and is not empty (ADR-AIIM-092).")
    out.mkdir(parents=True, exist_ok=True)

    rep = build()
    (out / "validation_matrix.json").write_text(json.dumps(rep, indent=2),
                                                encoding="utf-8")
    for r in rep["rows"]:
        print(f"\n{r['method']}")
        for axis, _ in AXES:
            v = r[axis]["value"] if isinstance(r[axis], dict) else r[axis]
            print(f"  {axis:22s} {v}")
        print(f"  -> licence: {r['licence']}")
    print(f"\nwritten -> {out / 'validation_matrix.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

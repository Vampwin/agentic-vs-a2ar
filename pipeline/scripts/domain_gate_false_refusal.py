"""domain_gate_false_refusal.py -- what the applicability-domain gate does to
known actives it has never seen.

The G6 contract names a false-refusal measurement on a panel of known actives.
Until 2026-08-29 nothing implemented it, and no such number appeared anywhere in
the manuscript: `_domain_gate_bias_probe.assess()` checked only how many
descriptor spaces a verdict rested on. Review comment M3 recorded that as an
instance of this paper's own defect class -- a contract asserting more than its
evaluator performs -- and this script is the other half of closing it
(ADR-AIIM-147).

**The panel and the reference set must be disjoint, and getting that wrong is
the whole trap.** The gate defines its domain from
`chembl_a2a_training_v3.csv`, all 5,970 compounds. Running any subset of those
through the gate returns a refusal rate near zero by construction, and the check
would report a pass while measuring nothing. Two candidate panels were rejected
for exactly this reason before the run:

  * `chembl_a2ar_validation_set.csv` shares 3,307 of its 3,372 compounds with
    the training file.
  * The clinical A2A antagonists -- istradefylline, preladenant, tozadenant,
    vipadenant -- are all in the training file. They would have made the most
    legible illustration and they are all inside the domain by construction.

What is used instead is the scaffold-disjoint split the QSAR pipeline already
performs: `scaffold_train_test_split` in `run_qsar_lbdd_2026-08-07.py`, which is
deterministic and reproduced here rather than re-derived. It yields 4,753
training and 1,217 held-out compounds sharing **no Murcko scaffold**. The
reference set is the 4,753; the panel is the 1,096 actives among the 1,217.

That reduced reference set is not the one production used, and that is the
point: measuring what a gate does to actives it has never seen requires holding
some out. Every other element of the rule -- five nearest neighbours, the tenth
percentile of reference self-similarity, the descriptor spaces -- is taken from
`ad_descriptor_sweep.py` unchanged.

Usage:
    python domain_gate_false_refusal.py [--out DIR]
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


def _root() -> Path:
    marker = "00_project_control/08_DECISION_LOG.md"
    for cand in (Path(r"D:\Agentic VS Workflow"), Path("/mnt/d/Agentic VS Workflow")):
        if (cand / marker).is_file():
            return cand
    raise SystemExit("project root not found")


ROOT = _root()
SCRIPTS = ROOT / "10_new_system/scripts"
TRAIN = ROOT / "02_frozen_inputs/public_benchmarks/chembl_a2ar/chembl_a2a_training_v3.csv"
sys.path.insert(0, str(ROOT / "10_new_system"))

from evaluators import _domain_gate_bias_probe as probe   # noqa: E402


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(mod)
    except SystemExit:
        pass
    return mod


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "20_new_runs" /
                                         "domain_gate_false_refusal_2026-08-29"))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    sweep = _load("ad_sweep", SCRIPTS / "ad_descriptor_sweep.py")
    qsar = _load("qsar", SCRIPTS / "run_qsar_lbdd_2026-08-07.py")

    df = pd.read_csv(TRAIN)
    smiles = df["canonical_smiles"].astype(str).tolist()
    labels = df["label"].to_numpy()

    tr_idx, te_idx = qsar.scaffold_train_test_split(smiles, labels)
    panel_idx = [i for i in te_idx if labels[i] == 1]

    def parse(idx):
        mols, keep = [], []
        for i in idx:
            m = Chem.MolFromSmiles(smiles[i])
            if m is not None:
                mols.append(m)
                keep.append(i)
        return mols, keep

    ref_mols, ref_keep = parse(tr_idx)
    panel_mols, panel_keep = parse(panel_idx)

    print(f"reference set (domain-defining) : {len(ref_mols)} parsed of {len(tr_idx)}")
    print(f"known-active panel              : {len(panel_mols)} parsed of {len(panel_idx)}")
    print(f"rule: in-domain if mean Tanimoto to {sweep.K}-NN >= {sweep.PCTL}th "
          f"percentile of reference self-similarity (frozen G6 rule, unchanged)\n")

    binary = ["ECFP4", "FCFP4", "MACCS", "AtomPair", "TopTorsion", "Gobbi_Pharm2D"]
    real = ["ErG", "PhysChem"]

    in_domain: dict[str, list[bool]] = {}
    detail = {}
    for space in binary + real:
        if space in binary:
            R = sweep.bitvec_matrix(ref_mols, space)
            P = sweep.bitvec_matrix(panel_mols, space)
            self_sim = sweep.knn_mean(sweep.tanimoto_block(R, R), sweep.K, True)
            panel_sim = sweep.knn_mean(sweep.tanimoto_block(P, R), sweep.K, False)
        else:
            R = sweep.realvec_matrix(ref_mols, space)
            P = sweep.realvec_matrix(panel_mols, space)
            self_sim = sweep.knn_mean(sweep.cosine_block(R, R), sweep.K, True)
            panel_sim = sweep.knn_mean(sweep.cosine_block(P, R), sweep.K, False)

        thr = float(np.percentile(self_sim, sweep.PCTL))
        inside = panel_sim >= thr
        in_domain[space] = [bool(x) for x in inside]
        detail[space] = {
            "threshold": round(thr, 6),
            "panel_mean_similarity": round(float(panel_sim.mean()), 6),
            "reference_self_similarity_median": round(float(np.median(self_sim)), 6),
        }
        print(f"  {space:<14} threshold {thr:.4f}   refused "
              f"{int((~inside).sum()):>4} / {len(inside)}  "
              f"({(~inside).mean():.1%})")

    # The natural-product library through the IDENTICAL reduced reference set.
    # Without this the paper would set a refusal rate measured against 4,753
    # reference compounds beside a coverage figure measured against 5,970, and
    # call them the same gate. They are the same rule; they were not the same
    # reference. Scoring both against 4,753 makes the comparison exact instead
    # of caveated.
    nplib = pd.read_csv(ROOT / "20_new_runs/np_library_funnel_2026-08-15/"
                               "e8_np_library_prioritization.csv", low_memory=False)
    np_mols = []
    for s in nplib["SMILES"].astype(str):
        m = Chem.MolFromSmiles(s)
        if m is not None:
            np_mols.append(m)
    print(f"\nnatural-product library through the same reduced reference: "
          f"{len(np_mols)} parsed of {len(nplib)}")
    np_refused = {}
    for space in binary + real:
        if space in binary:
            R = sweep.bitvec_matrix(ref_mols, space)
            N = sweep.bitvec_matrix(np_mols, space)
            self_sim = sweep.knn_mean(sweep.tanimoto_block(R, R), sweep.K, True)
            np_sim = sweep.knn_mean(sweep.tanimoto_block(N, R), sweep.K, False)
        else:
            R = sweep.realvec_matrix(ref_mols, space)
            N = sweep.realvec_matrix(np_mols, space)
            self_sim = sweep.knn_mean(sweep.cosine_block(R, R), sweep.K, True)
            np_sim = sweep.knn_mean(sweep.cosine_block(N, R), sweep.K, False)
        thr = float(np.percentile(self_sim, sweep.PCTL))
        inside = np_sim >= thr
        np_refused[space] = {
            "n": int(len(inside)),
            "n_refused": int((~inside).sum()),
            "refusal_rate": round(float((~inside).mean()), 4),
        }
        print(f"  {space:<14} refused {int((~inside).sum()):>4} / {len(inside)}  "
              f"({(~inside).mean():.1%})")

    report_probe = probe.assess_false_refusal(
        in_domain_by_space=in_domain,
        panel_disjoint_from_reference=(
            "The panel and the reference set are the two halves of "
            "scaffold_train_test_split() in run_qsar_lbdd_2026-08-07.py, applied "
            "to chembl_a2a_training_v3.csv. They share no Murcko scaffold: 4,753 "
            "compounds define the domain, 1,217 are held out, and the panel is "
            "the 1,096 actives among the held-out half. The gate's production "
            "reference set is all 5,970, which contains the panel; measuring a "
            "refusal rate against that would return zero by construction."),
        panel_description=(
            "1,096 ChEMBL A2A receptor actives (pChEMBL >= 6.0) whose Murcko "
            "scaffolds are absent from the compounds defining the domain."),
    )

    report = {
        "question": "How many compounds known to bind A2A does the "
                    "applicability-domain gate refuse, in each descriptor space, "
                    "when it has never seen their scaffolds?",
        "authority": "ADR-AIIM-147",
        "closes": "M3 -- the mechanism the G6 contract names and no code performed",
        "rule": {"k": sweep.K, "percentile": sweep.PCTL,
                 "definition": "in-domain if mean Tanimoto (or cosine) to k "
                               "nearest reference neighbours >= the given "
                               "percentile of reference self-similarity",
                 "source": "ad_descriptor_sweep.py, unchanged"},
        "reference_set": {"n": len(ref_mols),
                          "source": str(TRAIN.relative_to(ROOT).as_posix()),
                          "sha256": sha256(TRAIN),
                          "split": "scaffold_train_test_split, training half"},
        "panels_rejected_before_the_run": {
            "chembl_a2ar_validation_set.csv":
                "3,307 of its 3,372 compounds are in the training file",
            "clinical A2A antagonists":
                "istradefylline, preladenant, tozadenant and vipadenant are all "
                "in the training file and therefore in-domain by construction",
        },
        "per_space_detail": detail,
        "false_refusal": report_probe,
        "natural_product_library_same_reference": {
            "n": len(np_mols),
            "source": "20_new_runs/np_library_funnel_2026-08-15/"
                      "e8_np_library_prioritization.csv",
            "why": "Scored against the same 4,753-compound reference as the "
                   "known-active panel, so the two refusal rates are "
                   "comparable. The published coverage figures use the full "
                   "5,970-compound reference and are not repointed.",
            "by_space": np_refused,
        },
    }

    dest = out / "domain_gate_false_refusal.json"
    dest.write_text(json.dumps(report, indent=2), encoding="utf-8")

    fr = report_probe
    print(f"\nworst space : {fr['worst_space']}  {fr['worst_rate']:.1%}")
    print(f"best space  : {fr['best_space']}  {fr['best_rate']:.1%}")
    print(f"spread      : {fr['spread']:.1%}")
    if "structural_vs_physicochemical" in fr:
        sp = fr["structural_vs_physicochemical"]
        print(f"structural worst {sp['structural_worst']:.1%} vs "
              f"physicochemical worst {sp['physicochemical_worst']:.1%}")
    print(f"\n{fr['interpretation']}")
    print(f"\nwritten -> {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""run_clinical_reference_ad_gate.py -- do real A2A drugs pass the domain gate?

The question, put by the author on 2026-09-22: compounds that are actual drugs
or clinical candidates for this target should pass the applicability-domain
gate. If they do not, the gate is broken.

The first attempt at this control (`confidence_tier_positive_control_2026-08-08`)
used seven small purines -- xanthine, theobromine, theophylline and relatives.
All seven were refused, and that was read as a gate failure. Checking the
artifacts shows the control was mis-specified rather than the gate broken: none
of the seven appears in the 6,565-compound A2A pChEMBL pool at all, the only
small purine-core compound that does appear sits at pChEMBL 5.50 -- below the
model's own 6.0 active threshold -- and they carry 11 to 18 heavy atoms against
a training median of 28. A model fitted on sub-micromolar ligands should not be
confident about caffeine metabolites.

This is the properly specified version. The panel is the ten A2A reference
ligands of `clinical_reference_2026-09-03`, spanning pChEMBL 6.82 to 9.17, and
including one approved drug and four compounds that reached clinical trials.
These are potent, measured, and in the range the model was fitted for.

METHOD, identical to production
-------------------------------
ECFP4 (Morgan radius 2, 2048 bits), 5-nearest-neighbour mean Tanimoto against
the domain-defining set, compared with the frozen threshold 0.6810766.

**Self-exclusion is mandatory.** These are ChEMBL A2A ligands, so most are
likely inside the training set; scoring a compound against a set containing
itself returns a pass by construction. Each query is therefore scored against
the training set with its own record removed, and membership is reported per
compound so the reader can see which comparisons were self-excluded.

PREDICTION, recorded before the run
-----------------------------------
Unlike the seven purines, these sit in the model's own potency and size range,
so I expect **most or all ten to pass**. Any that fail are the informative
cases, and would be evidence that the gate refuses genuinely potent ligands
rather than only weak ones.

If a majority fail, the gate is refusing real A2A drugs and the author's
criterion is not met -- that result would be reported as the headline.

Output -> 20_new_runs/clinical_reference_ad_gate_2026-09-22/
"""
from __future__ import annotations

import argparse
import io
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Descriptors

RDLogger.DisableLog("rdApp.*")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kdss_paths import under  # noqa: E402

AD_K = 5
AD_THRESHOLD = 0.6810766499582288          # frozen, from run_qsar_lbdd_2026-08-07
TRAIN = under("02_frozen_inputs/public_benchmarks/chembl_a2ar",
              "chembl_a2a_training_v3.csv")
RAW = under("02_frozen_inputs/public_benchmarks/chembl_a2ar",
            "chembl_a2a_raw_v2.csv")
REF = under("20_new_runs", "clinical_reference_2026-09-03/clinical_reference.json")


def fp(smi: str):
    m = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
    return AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048) if m else None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(
        under("20_new_runs", "clinical_reference_ad_gate_2026-09-22")))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    tr = pd.read_csv(TRAIN)
    raw = pd.read_csv(RAW)
    smi_by_id = raw.drop_duplicates("molecule_chembl_id").set_index(
        "molecule_chembl_id")["canonical_smiles"]
    pch_by_id = raw.groupby("molecule_chembl_id")["pchembl_value"].median()

    ref = json.load(io.open(REF, encoding="utf-8"))
    panel = [(cid, v) for cid, v in ref["by_compound"].items()
             if v.get("chembl_id")] if isinstance(ref["by_compound"], dict) else []
    # by_compound is keyed by name; the ChEMBL id is inside the record
    rows = []
    for name, v in ref["by_compound"].items():
        cid = v.get("chembl_id") or ""
        if not cid:
            continue
        rows.append({"name": name, "chembl_id": cid,
                     "class": v.get("class"),
                     "pchembl_median": v.get("pchembl_median")})
    panel = pd.DataFrame(rows)
    panel["smiles"] = panel.chembl_id.map(smi_by_id)
    panel["pchembl_pool"] = panel.chembl_id.map(pch_by_id)

    train_ids = tr.molecule_chembl_id.tolist()
    train_smi = tr.canonical_smiles.tolist()
    train_fps = [fp(s) for s in train_smi]
    keep = [i for i, f in enumerate(train_fps) if f is not None]
    train_ids = [train_ids[i] for i in keep]
    train_fps = [train_fps[i] for i in keep]

    results = []
    for _, r in panel.iterrows():
        q = fp(r.smiles)
        if q is None:
            results.append({**r.to_dict(), "resolved": False})
            continue
        in_train = r.chembl_id in set(train_ids)
        idx = [i for i, cid in enumerate(train_ids) if cid != r.chembl_id]
        sims = DataStructs.BulkTanimotoSimilarity(q, [train_fps[i] for i in idx])
        score = float(np.mean(sorted(sims, reverse=True)[:AD_K]))
        m = Chem.MolFromSmiles(r.smiles)
        results.append({
            "name": r["name"], "chembl_id": r.chembl_id, "class": r["class"],
            "pchembl": r.pchembl_median, "heavy_atoms": m.GetNumHeavyAtoms(),
            "mw": round(Descriptors.MolWt(m), 1),
            "in_training_set": bool(in_train),
            "ad_5nn_mean_tanimoto_self_excluded": round(score, 4),
            "threshold": round(AD_THRESHOLD, 4),
            "passes_gate": bool(score >= AD_THRESHOLD),
            "resolved": True,
        })

    df = pd.DataFrame(results)
    df.to_csv(out / "clinical_reference_ad_gate.csv", index=False)
    ok = df[df.resolved]
    n_pass = int(ok.passes_gate.sum())

    summary = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "question": "Do potent, measured A2A reference ligands pass the "
                    "applicability-domain gate?",
        "method": "ECFP4 Morgan r=2 2048 bits, 5-NN mean Tanimoto against the "
                  "domain-defining training set with the query's own record "
                  "removed; frozen threshold 0.6810766",
        "self_exclusion_applied": True,
        "panel_source": "20_new_runs/clinical_reference_2026-09-03",
        "n_panel": int(len(ok)),
        "n_in_training_set": int(ok.in_training_set.sum()),
        "n_pass": n_pass,
        "n_fail": int(len(ok) - n_pass),
        "pass_rate": round(n_pass / len(ok), 4) if len(ok) else None,
        "pchembl_range": [float(ok.pchembl.min()), float(ok.pchembl.max())],
        "comparison": {
            "seven_small_purines_2026-08-08": "0 of 7 passed; mis-specified "
                                              "control, none in the A2A pChEMBL pool",
            "held_out_actives_probe": "1,096 potent scaffold-disjoint actives, "
                                      "false refusal 5.84 % (PhysChem) to "
                                      "14.96 % (MACCS)",
        },
        "per_compound": results,
    }
    with open(out / "clinical_reference_ad_gate.json", "w", encoding="utf-8",
              newline="\n") as f:
        json.dump(summary, f, indent=1, ensure_ascii=False)
        f.write("\n")

    print("=" * 78)
    print("DO REAL A2A DRUGS PASS THE APPLICABILITY-DOMAIN GATE?")
    print("=" * 78)
    print(f"threshold {AD_THRESHOLD:.4f} | 5-NN mean Tanimoto, ECFP4, self-excluded\n")
    print(f"{'compound':26} {'pChEMBL':>8} {'heavy':>6} {'inTrain':>8} "
          f"{'AD score':>9}  gate")
    for r in sorted(results, key=lambda x: -(x.get("ad_5nn_mean_tanimoto_self_excluded") or 0)):
        if not r.get("resolved"):
            print(f"{r['name'][:26]:26} {'-':>8} {'-':>6} {'-':>8} "
                  f"{'unresolved':>9}")
            continue
        print(f"{r['name'][:26]:26} {r['pchembl']:>8} {r['heavy_atoms']:>6} "
              f"{str(r['in_training_set']):>8} "
              f"{r['ad_5nn_mean_tanimoto_self_excluded']:>9.4f}  "
              f"{'PASS' if r['passes_gate'] else 'REFUSED'}")
    print()
    print(f"  {n_pass} of {len(ok)} pass  ({summary['pass_rate']:.0%})")
    print(f"  for comparison: 0 of 7 small purines passed on 2026-08-08")
    print()
    print(f"written -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

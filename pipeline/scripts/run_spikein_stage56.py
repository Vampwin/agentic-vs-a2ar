"""run_spikein_stage56.py -- M2.1 stages 5 and 6, which grade P1 and P3.

ADR-AIIM-256 completed stage 4 and graded **P2 MET**: of the 89 spike-ins the
funnel recovered, only 32 reach the production top-500 cut. **P1 and P3 were
emitted as PENDING rather than guessed**, because stages 5 and 6 -- the
applicability-domain gate and the ligand-based arm -- had not been run against
the spiked shortlist. ADR-AIIM-259 recorded that this is the one thing standing
between O2 and closure. This runs them.

### The predictions, from ADR-AIIM-229, unchanged

| P1 | the AD gate and the ligand-based arm admit most surviving spike-ins |
| P3 | if P1 and P2 hold, the loss is localised to docking rank and 3.5's undifferentiated negative becomes a diagnosis |

Falsification, also from ADR-AIIM-229: **if the AD gate rather than docking is
the largest loss, P1 is wrong** and ADR-AIIM-210's 9-of-10 does not generalise
past clinical-stage chemistry.

### What ADR-AIIM-229 left underspecified, fixed HERE and BEFORE any result

"Admit most surviving spike-ins" does not say which spike-ins, which rule, or
what "most" means. Fixed now, in the same way ADR-AIIM-246 fixed the stage-4
cut-off before the docking finished:

**Population.** The **primary** population is the **32 spike-ins inside the
production top-500 cut**, because those are what a deployed run would actually
hand to stages 5 and 6 -- the tiering operates on the redocked 500. The
**secondary** population is **all 89 docked spike-ins**, reported beside it,
which answers the different question of whether the gate would have lost them
regardless of rank. Neither is chosen after seeing the other: both are written
here, and both are reported whatever they say.

**Rules, both taken from production and neither re-tuned.**

* *AD gate* -- ECFP4 (Morgan radius 2, 2048 bits), 5-nearest-neighbour mean
  Tanimoto against the frozen ChEMBL A2AR training set, compared with the
  frozen threshold **0.6810766499582288**. Identical to
  `run_clinical_reference_ad_gate.py`, constants imported from it rather than
  restated, so the two cannot drift.
* *Ligand-based arm* -- the frozen Mondrian conformal predictor, loaded and not
  refit, `prediction_set` at alpha = 0.05. A compound is **admitted** when that
  set is not `{0}`, i.e. the arm does **not** confidently exclude activity.
  This is the dual-criterion rule's own threshold, not a new one.

**Grading.**

    P1a  AD gate admits            > 50 % of the primary population
    P1b  ligand-based arm admits   > 50 % of the primary population
    P1   MET if and only if P1a AND P1b
    P3   MET if and only if P1 AND P2, where P2 is read from the stage-4
         artifact rather than recomputed

Reporting both halves separately matters: P1 is a conjunction, and a conjunction
that fails tells you nothing about *which* conjunct failed unless both are
published.

### Self-exclusion is not needed here, and that is checked rather than assumed

`run_clinical_reference_ad_gate.py` had to exclude each query's own record
because its panel were ChEMBL A2A ligands inside the training set. **These
spike-ins were selected precisely because they are not** (ADR-AIIM-222: the
panel is drawn from the 169 DUD-E actives whose InChIKey connectivity block is
absent from the 6,565-compound A2A set). The script verifies that and
**refuses to run if any spike-in is found in the domain-defining set**, because
a spike-in scored against a set containing itself would pass by construction and
inflate P1.

### What this does not do

It does not apply the dual-criterion tier rule, and says so rather than
pretending to: Tier A requires an **on-pocket pose**, and the spike-in docking
recorded scores, not poses. The tier rule is therefore reported as *not
evaluable from this run*. P1 does not depend on it -- it names the AD gate and
the ligand-based arm -- and P3 is a statement about where the loss sits, which
the two stage rates answer.

Nothing is docked, nothing is refit, no threshold is moved.

Usage:
    python run_spikein_stage56.py --self-test     # falsifiability, no real data
    python run_spikein_stage56.py

Output -> 20_new_runs/spike_in_dock_2026-09-23/analysis/
"""
from __future__ import annotations

import argparse
import io
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kdss_paths import under  # noqa: E402

# Imported, never restated, so the gate here and the gate in the clinical
# control cannot drift apart.
from run_clinical_reference_ad_gate import AD_K, AD_THRESHOLD, TRAIN, fp  # noqa: E402

from rdkit import DataStructs, RDLogger  # noqa: E402

RDLogger.DisableLog("rdApp.*")
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

RUN = Path(under("20_new_runs", "spike_in_dock_2026-09-23"))
SHORTLIST = Path(under("20_new_runs", "spike_in_recovery_2026-09-23",
                       "np_library_docking_shortlist.csv"))
STAGE4 = RUN / "analysis" / "spikein_recovery.json"
OUT = RUN / "analysis"
ALPHA = "0.05"
PRIMARY_CUTOFF = 500


def wilson(k: int, n: int) -> list[float]:
    if n == 0:
        return [0.0, 0.0]
    z, p = 1.959963984540054, k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * ((p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5) / d
    return [round(max(0.0, c - h), 4), round(min(1.0, c + h), 4)]


def rate(k: int, n: int) -> dict:
    return {"k": int(k), "n": int(n),
            "rate": round(k / n, 4) if n else None, "wilson95": wilson(k, n)}


def ad_admits(smiles: list[str], train_fps: list) -> np.ndarray:
    """Production AD gate: 5-NN mean ECFP4 Tanimoto >= the frozen threshold."""
    out = np.zeros(len(smiles), dtype=bool)
    for i, s in enumerate(smiles):
        f = fp(s)
        if f is None:
            continue
        sims = np.array(DataStructs.BulkTanimotoSimilarity(f, train_fps))
        out[i] = float(np.sort(sims)[-AD_K:].mean()) >= AD_THRESHOLD
    return out


def grade(ad: np.ndarray, lb: np.ndarray, p2_met: bool) -> dict:
    n = len(ad)
    p1a = bool(ad.sum() > n / 2)
    p1b = bool(lb.sum() > n / 2)
    p1 = bool(p1a and p1b)
    return {"P1a_ad_gate_admits_most": p1a, "P1b_ligand_arm_admits_most": p1b,
            "P1_met": p1, "P2_met": bool(p2_met), "P3_met": bool(p1 and p2_met)}


def self_test() -> int:
    """Both predictions must be reachable, or the grading proves nothing."""
    print("self-test: grading is falsifiable in both directions")
    ok = True
    cases = [
        ("both arms admit almost everything", np.ones(32, bool), np.ones(32, bool),
         True, True),
        ("AD gate refuses almost everything", np.zeros(32, bool), np.ones(32, bool),
         True, False),
        ("ligand arm refuses almost everything", np.ones(32, bool), np.zeros(32, bool),
         True, False),
        ("both admit but P2 had missed", np.ones(32, bool), np.ones(32, bool),
         False, False),
    ]
    for name, ad, lb, p2, want_p3 in cases:
        g = grade(ad, lb, p2)
        good = g["P3_met"] == want_p3
        ok &= good
        print(f"  {'ok ' if good else 'FAIL'} {name:38s} "
              f"P1={g['P1_met']} P2={g['P2_met']} P3={g['P3_met']}")
    # a knife-edge case: exactly half is NOT "most"
    g = grade(np.array([True] * 16 + [False] * 16), np.ones(32, bool), True)
    good = g["P1a_ad_gate_admits_most"] is False
    ok &= good
    print(f"  {'ok ' if good else 'FAIL'} exactly half is not 'most'            "
          f"P1a={g['P1a_ad_gate_admits_most']}")
    print("self-test:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args()
    if args.self_test:
        return self_test()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    stage4 = json.load(io.open(STAGE4, encoding="utf-8"))
    if stage4.get("provisional"):
        print("REFUSING: the stage-4 artifact is provisional.")
        return 2
    p2_met = bool(stage4["predictions"]["P2"]["met"])

    scores = pd.read_csv(RUN / "scores.csv")
    short = pd.read_csv(SHORTLIST)
    d = scores.merge(short[["structure_inchikey", "smiles_standardized"]],
                     on="structure_inchikey", how="left")
    if d.smiles_standardized.isna().any():
        print(f"REFUSING: {d.smiles_standardized.isna().sum()} rows without SMILES.")
        return 2
    d = d[d.vina_score.notna()].copy()
    d["rank"] = d.vina_score.rank(method="first").astype(int)
    d = d.sort_values("vina_score").reset_index(drop=True)
    d["rank"] = np.arange(1, len(d) + 1)

    spikes = d[d.spike_in.astype(bool)].copy()
    print(f"{len(d)} scored, {len(spikes)} spike-ins")

    # ---- the ligand-based arm, frozen predictor, run as its own script -----
    pred_csv = out / "spikein_conformal_predictions.csv"
    if not pred_csv.exists():
        q = out / "_conformal_query.csv"
        d[["structure_inchikey", "smiles_standardized"]].to_csv(q, index=False)
        cmd = [sys.executable, str(Path(__file__).parent / "conformal_qsar.py"),
               "--query-csv", str(q), "--query-smiles-col", "smiles_standardized",
               "--query-id-col", "structure_inchikey",
               "--out", str(out), "--out-name", pred_csv.name]
        print("running the frozen conformal predictor ...")
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            print(r.stdout[-2000:]); print(r.stderr[-2000:])
            print("REFUSING: the conformal predictor failed.")
            return 2
    pred = pd.read_csv(pred_csv)
    setcol = f"prediction_set_sig{ALPHA}"
    if setcol not in pred.columns:
        print(f"REFUSING: {setcol} missing from the predictor output.")
        return 2
    lb_map = dict(zip(pred.structure_inchikey,
                      pred[setcol].astype(str).str.strip() != "{0}"))

    # ---- the AD gate, production rule, imported constants ------------------
    tr = pd.read_csv(TRAIN)
    train_fps = [f for f in (fp(s) for s in tr.canonical_smiles) if f is not None]
    print(f"domain-defining set: {len(train_fps)} fingerprints, "
          f"k={AD_K}, threshold={AD_THRESHOLD:.7f}")

    # membership control -- a spike-in inside the domain set would pass by
    # construction, so the run refuses rather than reporting an inflated P1.
    tr_keys = set(tr.canonical_smiles.astype(str))
    inside = spikes[spikes.smiles_standardized.astype(str).isin(tr_keys)]
    if len(inside):
        print(f"REFUSING: {len(inside)} spike-in(s) appear in the "
              f"domain-defining set; self-exclusion would be required.")
        return 2
    print("membership control: 0 spike-ins inside the domain-defining set")

    spikes["ad_admits"] = ad_admits(spikes.smiles_standardized.tolist(), train_fps)
    spikes["lb_admits"] = spikes.structure_inchikey.map(lb_map)
    if spikes.lb_admits.isna().any():
        print("REFUSING: some spike-ins carry no conformal prediction.")
        return 2
    spikes["lb_admits"] = spikes.lb_admits.astype(bool)

    primary = spikes[spikes["rank"] <= PRIMARY_CUTOFF]
    pops = {"primary_top500": primary, "secondary_all_docked": spikes}
    report_pops = {}
    for name, g in pops.items():
        report_pops[name] = {
            "n": int(len(g)),
            "ad_gate_admits": rate(int(g.ad_admits.sum()), len(g)),
            "ligand_arm_admits": rate(int(g.lb_admits.sum()), len(g)),
            "both_admit": rate(int((g.ad_admits & g.lb_admits).sum()), len(g)),
        }

    g = grade(primary.ad_admits.to_numpy(), primary.lb_admits.to_numpy(), p2_met)

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "executes": "M2.1 stages 5 and 6",
        "pre_registration": "ADR-AIIM-229; operationalised in this script's "
                            "docstring before any result was read",
        "grades": "P1 and P3; P2 is read from the stage-4 artifact, not recomputed",
        "rules": {
            "ad_gate": f"ECFP4 r=2 2048 bits, {AD_K}-NN mean Tanimoto vs the "
                       f"frozen ChEMBL A2AR training set, threshold "
                       f"{AD_THRESHOLD}; constants imported from "
                       f"run_clinical_reference_ad_gate.py",
            "ligand_arm": f"frozen Mondrian conformal predictor, loaded not "
                          f"refit; admitted when the alpha={ALPHA} prediction "
                          f"set is not {{0}}",
            "primary_population": "the spike-ins inside the production top-500 "
                                  "cut, because that is what a deployed run "
                                  "hands to stages 5 and 6",
        },
        "membership_control": {
            "spike_ins_in_domain_defining_set": 0,
            "why": "ADR-AIIM-222 drew the panel from DUD-E actives absent from "
                   "the A2A set; scoring a compound against a set containing "
                   "itself would pass by construction",
        },
        "populations": report_pops,
        "predictions": {
            "P1": {"statement": "the AD gate and the ligand-based arm admit "
                                "most surviving spike-ins", **g},
            "P2": {"statement": "the structure-based rank step loses the majority",
                   "met": p2_met, "source": "stage-4 artifact, not recomputed"},
            "P3": {"statement": "if P1 and P2 hold the loss is localised to "
                                "docking rank",
                   "met": g["P3_met"]},
        },
        "not_evaluable": [
            "the dual-criterion tier rule: Tier A requires an on-pocket pose "
            "and the spike-in docking recorded scores, not poses. P1 names the "
            "AD gate and the ligand-based arm and does not depend on it.",
        ],
        "not_claimed": [
            "that any spike-in is active in any assay run by this project",
            "that admitting a known active makes the gate correct in general; "
            "this is one panel, one target, one representation",
        ],
    }
    (out / "spikein_stage56.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")
    spikes.to_csv(out / "spikein_stage56_per_compound.csv", index=False)

    print("=" * 74)
    print("M2.1 STAGES 5 AND 6")
    print("=" * 74)
    for name, r in report_pops.items():
        print(f"  {name}  (n = {r['n']})")
        for k in ("ad_gate_admits", "ligand_arm_admits", "both_admit"):
            v = r[k]
            print(f"    {k:20s} {v['k']:3d} / {v['n']:3d}  {v['rate']:.3f}  "
                  f"[{v['wilson95'][0]:.3f}, {v['wilson95'][1]:.3f}]")
    print(f"\n  P1a AD gate admits most     : {g['P1a_ad_gate_admits_most']}")
    print(f"  P1b ligand arm admits most  : {g['P1b_ligand_arm_admits_most']}")
    print(f"  P1  : {'MET' if g['P1_met'] else 'MISSED'}")
    print(f"  P2  : {'MET' if p2_met else 'MISSED'}  (from stage 4)")
    print(f"  P3  : {'MET' if g['P3_met'] else 'MISSED'}")
    print(f"\nwritten -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

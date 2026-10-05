"""analyse_dude_enrichment.py -- the M5.1 analysis, written before the data exists.

Executes the analysis pre-registered in ADR-AIIM-214 and amended by
ADR-AIIM-215. **This script was written while the docking run was at roughly
half its compounds and no score had been looked at.** That is the point: an
analysis written after the results are visible is not the analysis that was
pre-registered, whatever it claims in its header.

Every metric is **imported** from `run_matched_screening_benchmark.py`, which
holds this project's frozen conventions -- `ceil(n*f)` enrichment, mergesort
tie handling, 1,000-sample bootstrap at seed 42, BEDROC alpha = 20. Reusing
them is what makes the DUD-E number comparable with the 245-compound benchmark
rather than merely similar to it.

The pre-registration fixes, and this script implements:

  PRIMARY    paired delta-AUC, Vina minus the **better-performing** trivial
             null, with a bootstrap interval. "Better-performing" is a rule
             applied to the data, not a choice made after seeing it: of the two
             declared nulls, the one with the higher AUC on the identical rows.

  NULLS      heavy-atom count and cLogP, on the identical rows, **both reported
             whichever way they fall**. cLogP was the stronger on the full set
             (0.594 against 0.540) and reporting heavy atoms alone would
             understate the baseline (ADR-AIIM-212).

  SECONDARY  whether Vina's own AUC interval excludes 0.5; BEDROC20 and EF at
             1 / 5 / 10 % against the same nulls.

Both nulls are scored in their raw direction, higher value ranking higher,
because that is the direction the pre-check measured them in and both came out
above 0.5 there. Vina is negated, its score being negative-is-better.

The three predictions of ADR-AIIM-214 are evaluated **mechanically** and their
verdicts written to the artifact, so that a disappointing result cannot be
reworded into a qualified success afterwards.

Refuses to run on a partial file. Analysing a half-finished benchmark and
calling the number provisional is how a provisional number ends up in a table.

Output -> 20_new_runs/dude_enrichment_2026-09-22/analysis/
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kdss_paths import under  # noqa: E402

# One implementation of every metric, shared with the 245-compound benchmark.
from run_matched_screening_benchmark import (  # noqa: E402
    BEDROC_ALPHA,
    BOOTSTRAP_N,
    BOOTSTRAP_SEED,
    EF_FRACTIONS,
    bedroc,
    boot_ci,
    enrichment,
    paired_delta_ci,
    roc_auc,
)

RDLogger.DisableLog("rdApp.*")

RUN = Path(under("20_new_runs", "dude_enrichment_2026-09-22"))
SCORES = RUN / "scores_shard0.csv"
SET = Path(under("02_frozen_inputs",
                 "public_benchmarks/dude_a2ar/dude_a2ar_docking_set.csv"))
N_EXPECTED = 2402

# ADR-AIIM-214, section 6. Recorded here so the script grades itself.
PREDICTIONS = {
    "P1": {"statement": "Vina ROC-AUC lands in 0.62 - 0.78",
           "lo": 0.62, "hi": 0.78},
    "P2": {"statement": "the paired delta-AUC against the better trivial null "
                        "is positive and its 95 % interval excludes zero"},
    "P3": {"statement": "BEDROC20 exceeds the corresponding null"},
}


def descriptors(smiles: str) -> tuple[float, float] | None:
    m = Chem.MolFromSmiles(str(smiles))
    if m is None:
        return None
    return float(m.GetNumHeavyAtoms()), float(Crippen.MolLogP(m))


def block(y, s, name: str) -> dict:
    lo, hi = boot_ci(y, s, roc_auc)
    return {
        "name": name,
        "roc_auc": round(float(roc_auc(y, s)), 4),
        "roc_auc_ci95": [round(lo, 4), round(hi, 4)],
        "auc_interval_excludes_0.5": bool(lo > 0.5 or hi < 0.5),
        "bedroc20": round(float(bedroc(y, s, BEDROC_ALPHA)), 4),
        "ef": {f"{int(f*100)}%": round(float(enrichment(y, s, f)), 4)
               for f in EF_FRACTIONS},
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(RUN / "analysis"))
    ap.add_argument("--scores", default=str(SCORES),
                    help="override the score file; used only by the "
                         "synthetic smoke test, which exercises this code "
                         "without anyone reading the real scores")
    ap.add_argument("--allow-partial", action="store_true",
                    help="analyse an unfinished run; the output is marked "
                         "PROVISIONAL and must not be reported")
    args = ap.parse_args()
    out = Path(args.out)

    d = pd.read_csv(args.scores)
    n_rows = len(d)
    scored = d[d.vina_score.notna() & (d.vina_score != "")].copy()
    scored["vina_score"] = scored.vina_score.astype(float)

    if n_rows < N_EXPECTED and not args.allow_partial:
        print(f"REFUSING: {n_rows} of {N_EXPECTED} compounds are scored.")
        print("A half-finished benchmark analysed 'provisionally' is how a")
        print("provisional number reaches a table. Re-run when the docking")
        print("completes, or pass --allow-partial and accept the label.")
        return 1

    # The protocol gate: one conformer docked, or this is not comparable with
    # the ChEMBL 245 re-run.
    bad = scored[scored.n_conformers_docked != 1]
    if len(bad):
        print(f"ABORT: {len(bad)} rows were docked with a conformer count "
              f"other than 1. This run cannot be compared across protocols.")
        return 2

    # The frozen set has 2,402 rows but 2,401 distinct compound_id: DUD-E gave
    # C19324335 to both oxime isomers. Joining on that key without collapsing it
    # turns two scored rows into four, which is how a 2,402-row benchmark
    # reported 2,404 scored and 1,922 decoys (ADR-AIIM-233). Applied before any
    # real score was read; the docking was unfinished when this was written.
    # compound_id is NOT a key of the frozen set: it has 2,402 rows and 2,401
    # distinct ids, because DUD-E issued C19324335 to both oxime isomers. The
    # docking walks the set by row and scores both, so the score file carries
    # two rows under that id as well -- two different molecules, correctly.
    #
    # Joining on the id alone turned 2,402 rows into 2,404 (ADR-AIIM-233).
    # De-duplicating the reference, which that entry proposed, fixes the count
    # but hands both scored rows the SAME SMILES, so one of the two isomers is
    # described by the other's structure. The pairing has to survive instead:
    # each occurrence of an id joins to the corresponding occurrence in the set,
    # which is well defined because the docking preserves set order.
    ref = pd.read_csv(SET)[["compound_id", "SMILES"]]
    ref["compound_id"] = ref.compound_id.astype(str)
    scored["compound_id"] = scored.compound_id.astype(str)
    ref["_occ"] = ref.groupby("compound_id").cumcount()
    scored["_occ"] = scored.groupby("compound_id").cumcount()

    m = scored.merge(ref, on=["compound_id", "_occ"], how="left")
    if len(m) != len(scored):
        print(f"ABORT: the SMILES join changed the row count, "
              f"{len(scored)} -> {len(m)}.")
        return 5
    if m.SMILES.isna().any():
        print(f"ABORT: {int(m.SMILES.isna().sum())} scored compounds have no "
              f"SMILES in the frozen set.")
        return 3

    desc = m.SMILES.map(descriptors)
    if desc.isna().any():
        print(f"ABORT: {int(desc.isna().sum())} SMILES did not parse.")
        return 4
    m["heavy_atoms"] = [x[0] for x in desc]
    m["clogp"] = [x[1] for x in desc]

    y = m.label.to_numpy(int)
    s_vina = -m.vina_score.to_numpy(float)      # negative-is-better
    s_heavy = m.heavy_atoms.to_numpy(float)     # raw direction, per the pre-check
    s_clogp = m.clogp.to_numpy(float)

    arms = {"vina": block(y, s_vina, "Vina, one minimised conformer"),
            "heavy_atoms": block(y, s_heavy, "heavy-atom count"),
            "clogp": block(y, s_clogp, "cLogP")}

    # The rule, applied rather than chosen: the stronger of the two declared
    # nulls on these rows.
    better = max(("heavy_atoms", "clogp"), key=lambda k: arms[k]["roc_auc"])
    s_best = {"heavy_atoms": s_heavy, "clogp": s_clogp}[better]

    lo, hi, p_le0 = paired_delta_ci(y, s_vina, s_best)
    delta = float(roc_auc(y, s_vina) - roc_auc(y, s_best))
    primary = {
        "comparator": better,
        "comparator_auc": arms[better]["roc_auc"],
        "delta_auc": round(delta, 4),
        "delta_ci95": [round(lo, 4), round(hi, 4)],
        "interval_excludes_zero": bool(lo > 0 or hi < 0),
        "bootstrap_fraction_delta_le_zero": round(p_le0, 4),
        "vina_beats_the_better_null": bool(delta > 0 and lo > 0),
    }
    # Reported for completeness; the primary is against the better one.
    other = "clogp" if better == "heavy_atoms" else "heavy_atoms"
    s_other = {"heavy_atoms": s_heavy, "clogp": s_clogp}[other]
    lo2, hi2, _ = paired_delta_ci(y, s_vina, s_other)
    secondary_delta = {
        "comparator": other,
        "delta_auc": round(float(roc_auc(y, s_vina) - roc_auc(y, s_other)), 4),
        "delta_ci95": [round(lo2, 4), round(hi2, 4)],
    }

    verdicts = {
        "P1": {**PREDICTIONS["P1"],
               "met": bool(PREDICTIONS["P1"]["lo"] <= arms["vina"]["roc_auc"]
                           <= PREDICTIONS["P1"]["hi"]),
               "observed": arms["vina"]["roc_auc"]},
        "P2": {**PREDICTIONS["P2"],
               "met": primary["vina_beats_the_better_null"],
               "observed": {"delta": primary["delta_auc"],
                            "ci95": primary["delta_ci95"]}},
        "P3": {**PREDICTIONS["P3"],
               "met": bool(arms["vina"]["bedroc20"] > arms[better]["bedroc20"]),
               "observed": {"vina": arms["vina"]["bedroc20"],
                            better: arms[better]["bedroc20"]}},
    }

    # ADR-AIIM-214, section 7, written before the run and applied without
    # rewording.
    if primary["vina_beats_the_better_null"]:
        licensed = ("The protocol has demonstrated screening power on "
                    "chemotypes resembling the co-crystal ligand, and the "
                    "0.601 of section 3.4 is a benchmark artifact rather than "
                    "a protocol failure. It does NOT license ranking the "
                    "natural-product library: ADR-AIIM-214 section 8 records "
                    "that in advance.")
    else:
        licensed = ("The protocol has no demonstrated screening power on a "
                    "benchmark constructed to be favourable to it. O5 closes. "
                    "This is the headline and the screening arm is presented "
                    "as a cautionary methodological result only.")

    out.mkdir(parents=True, exist_ok=True)
    m.to_csv(out / "dude_enrichment_per_compound.csv", index=False)
    summary = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "executes": "M5.1",
        "pre_registration": "ADR-AIIM-214, amended by ADR-AIIM-215",
        "analysis_written": "before the docking run finished and before any "
                            "score was inspected",
        "provisional": bool(n_rows < N_EXPECTED),
        "n_rows_in_file": int(n_rows),
        "n_scored": int(len(m)),
        "n_failed_to_dock": int(len(d) - len(scored)),
        "n_actives": int((m.label == 1).sum()),
        "n_decoys": int((m.label == 0).sum()),
        "protocol": {"conformers_generated": 8, "conformers_docked": 1,
                     "exhaustiveness": 32, "poses": 10, "seed": 20260729},
        "metrics_from": "run_matched_screening_benchmark.py (frozen "
                        "conventions, shared with the 245-compound benchmark)",
        "bootstrap": {"n": BOOTSTRAP_N, "seed": BOOTSTRAP_SEED},
        "arms": arms,
        "primary_outcome": primary,
        "secondary_delta_against_the_weaker_null": secondary_delta,
        "predictions": verdicts,
        "what_this_licenses": licensed,
        "not_claimed": [
            "that any compound binds; nothing here has been assayed",
            "that a favourable result licenses ranking the natural-product "
            "library (ADR-AIIM-214 section 8)",
            "anything about the ligand-based arm: 65 % of DUD-E actives are in "
            "its training set (ADR-AIIM-222)",
        ],
    }
    with open(out / "dude_enrichment_analysis.json", "w", encoding="utf-8",
              newline="\n") as f:
        json.dump(summary, f, indent=1, ensure_ascii=False)
        f.write("\n")

    print("=" * 78)
    print("DUD-E A2AR ENRICHMENT   (M5.1)"
          + ("   ** PROVISIONAL **" if summary["provisional"] else ""))
    print("=" * 78)
    print(f"  {len(m)} scored: {summary['n_actives']} active, "
          f"{summary['n_decoys']} decoy; {summary['n_failed_to_dock']} failed")
    print()
    print(f"  {'arm':32}{'ROC-AUC':>10}{'95 % CI':>20}{'BEDROC20':>11}")
    print("  " + "-" * 73)
    for k in ("vina", "heavy_atoms", "clogp"):
        a = arms[k]
        ci = f"[{a['roc_auc_ci95'][0]:.3f}, {a['roc_auc_ci95'][1]:.3f}]"
        print(f"  {a['name']:32}{a['roc_auc']:>10.4f}{ci:>20}"
              f"{a['bedroc20']:>11.4f}")
    print()
    print(f"  PRIMARY  Vina minus {better} (the better null):")
    print(f"           delta = {primary['delta_auc']:+.4f}  "
          f"95 % CI [{primary['delta_ci95'][0]:+.4f}, "
          f"{primary['delta_ci95'][1]:+.4f}]")
    print(f"           interval excludes zero: "
          f"{primary['interval_excludes_zero']}")
    print()
    print("  predictions recorded before the run:")
    for k, v in verdicts.items():
        print(f"    [{'MET ' if v['met'] else 'MISS'}] {k}: {v['statement']}")
    print()
    print(f"  {licensed}")
    print(f"\nwritten -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

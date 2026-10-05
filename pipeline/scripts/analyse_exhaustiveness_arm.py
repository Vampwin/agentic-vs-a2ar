"""analyse_exhaustiveness_arm.py -- the M2.2 analysis, written before any score
of the arm was read.

Executes the analysis pre-registered in ADR-AIIM-231. The exhaustiveness-8
docking finished at 06:29 on 2026-09-23; **this script was written after that and
before a single `vina_score` of it was inspected.** Only the row count, the
conformer column, the error count and the elapsed times had been looked at, none
of which is a result. That is a weaker guarantee than `analyse_dude_enrichment.py`
had -- which was written while its run was still going -- and it is stated here
rather than glossed, because the difference is real.

### The question

The frozen protocol docks at exhaustiveness 32. Vina's documented default is 8.
Does quadrupling the Monte-Carlo search budget change screening power on the
identical rows?

### Paired, on the rows both arms scored

Every compound is compared with itself. Rows that failed to dock in either arm
are dropped from both, so the comparison never rests on a set difference. Both
arms dock one minimised conformer of eight generated, the same receptor, box,
seed and pose count; **exhaustiveness is the only thing that differs.**

### Metrics are imported, not reimplemented

From `run_matched_screening_benchmark.py`, exactly as the M5.1 analysis does:
`ceil(n*f)` enrichment, mergesort tie handling, 1,000-sample bootstrap at seed
42, BEDROC alpha 20. That is what makes this arm comparable with M5.1 and with
the 245-compound benchmark rather than merely similar to them.

### The three predictions, graded mechanically

ADR-AIIM-231 fixed them before the run and this script writes its verdicts to the
artifact, so a disappointing result cannot be reworded afterwards.

  P1  no difference: |delta AUC| < 0.05 and the paired bootstrap interval
      includes zero. The basis was that exhaustiveness sets restarts, and a
      method already near chance at 32 has little room to lose.
  P2  individual scores move even if the ranking does not. This is the arm's
      own validity check: if per-compound scores were identical the flag did
      not take effect and the run is void, whatever the AUC says.
  P3  wall-clock falls roughly four-fold, tracking the budget.

**P2 is checked first and is fatal.** The other two are not reported if it fails.

### Refusals

A partial file in either arm, and any scored row not docked with exactly one
conformer, abort. A mixed-protocol comparison measures the protocol.

Output -> 20_new_runs/dude_exhaustiveness8_2026-09-23/analysis/
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
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ARM32 = Path(under("20_new_runs", "dude_enrichment_2026-09-22",
                   "scores_shard0.csv"))
ARM8 = Path(under("20_new_runs", "dude_exhaustiveness8_2026-09-23",
                  "scores_shard0.csv"))
SET = Path(under("02_frozen_inputs",
                 "public_benchmarks/dude_a2ar/dude_a2ar_docking_set.csv"))
OUT = Path(under("20_new_runs", "dude_exhaustiveness8_2026-09-23", "analysis"))
N_EXPECTED = 2402

PREDICTIONS = {
    "P1": {"statement": "no difference: |delta AUC| < 0.05 and the paired "
                        "bootstrap 95 % interval includes zero",
           "max_abs_delta": 0.05},
    "P2": {"statement": "per-compound scores move between the two budgets; "
                        "identical scores would mean the flag did not take "
                        "effect and the run is void"},
    "P3": {"statement": "wall-clock per compound falls roughly four-fold",
           "lo": 2.0, "hi": 8.0},
}


def descriptors(smiles: str) -> tuple[float, float] | None:
    m = Chem.MolFromSmiles(str(smiles))
    if m is None:
        return None
    return float(m.GetNumHeavyAtoms()), float(Crippen.MolLogP(m))


def load_arm(path: Path, name: str, allow_partial: bool) -> pd.DataFrame:
    d = pd.read_csv(path)
    if len(d) < N_EXPECTED and not allow_partial:
        raise SystemExit(
            f"REFUSING: {name} has {len(d)} of {N_EXPECTED} rows. A half-finished "
            f"benchmark analysed 'provisionally' is how a provisional number "
            f"reaches a table.")
    d["compound_id"] = d.compound_id.astype(str)
    scored = d[d.vina_score.notna() & (d.vina_score != "")].copy()
    scored["vina_score"] = scored.vina_score.astype(float)
    bad = scored[scored.n_conformers_docked != 1]
    if len(bad):
        raise SystemExit(
            f"ABORT: {len(bad)} rows in {name} were docked with a conformer count "
            f"other than 1. This cannot be compared across protocols.")
    # compound_id is not a key of the frozen set -- C19324335 labels two oxime
    # isomers -- so occurrences are paired positionally, as ADR-AIIM-233 fixed
    # in the M5.1 analysis.
    scored["_occ"] = scored.groupby("compound_id").cumcount()
    return scored


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--arm32", default=str(ARM32))
    ap.add_argument("--arm8", default=str(ARM8))
    ap.add_argument("--allow-partial", action="store_true",
                    help="used only by the synthetic smoke test")
    args = ap.parse_args()
    out = Path(args.out)

    a32 = load_arm(Path(args.arm32), "exhaustiveness 32", args.allow_partial)
    a8 = load_arm(Path(args.arm8), "exhaustiveness 8", args.allow_partial)

    m = a32.merge(a8, on=["compound_id", "_occ"], suffixes=("_32", "_8"))
    if (m.label_32 != m.label_8).any():
        raise SystemExit("ABORT: labels disagree between arms on shared rows.")
    y = m.label_32.to_numpy().astype(int)
    s32 = -m.vina_score_32.to_numpy()          # negative-is-better -> negate
    s8 = -m.vina_score_8.to_numpy()

    # ---- P2 first: did the flag take effect at all? ------------------------
    diff = m.vina_score_8.to_numpy() - m.vina_score_32.to_numpy()
    n_moved = int((np.abs(diff) > 1e-9).sum())
    p2 = n_moved > 0
    verdicts = {"P2": {"statement": PREDICTIONS["P2"]["statement"], "met": p2,
                       "observed": {"n_rows_compared": len(m),
                                    "n_scores_that_moved": n_moved,
                                    "fraction_moved": round(n_moved / len(m), 4)}}}
    if not p2:
        (out).mkdir(parents=True, exist_ok=True)
        (out / "exhaustiveness_arm_analysis.json").write_text(json.dumps(
            {"void": True, "reason": "every per-compound score is identical "
             "between the two budgets; the exhaustiveness flag did not take "
             "effect and no other number from this run may be reported.",
             "predictions": verdicts}, indent=2), encoding="utf-8")
        print("VOID: per-compound scores are identical between the arms.")
        return 3

    # ---- arms and nulls on the identical rows ------------------------------
    ref = pd.read_csv(SET)[["compound_id", "SMILES"]]
    ref["compound_id"] = ref.compound_id.astype(str)
    ref["_occ"] = ref.groupby("compound_id").cumcount()
    m = m.merge(ref, on=["compound_id", "_occ"], how="left")
    if m.SMILES.isna().any():
        raise SystemExit(f"ABORT: {int(m.SMILES.isna().sum())} rows have no SMILES.")
    desc = m.SMILES.map(descriptors)
    if desc.isna().any():
        raise SystemExit(f"ABORT: {int(desc.isna().sum())} SMILES did not parse.")
    heavy = np.array([d[0] for d in desc])
    clogp = np.array([d[1] for d in desc])

    def arm(name, s):
        lo, hi = boot_ci(y, s, roc_auc)
        return {"name": name, "roc_auc": round(float(roc_auc(y, s)), 4),
                "roc_auc_ci95": [round(lo, 4), round(hi, 4)],
                "bedroc20": round(float(bedroc(y, s, BEDROC_ALPHA)), 4),
                "ef": {f"{int(f*100)}%": round(float(enrichment(y, s, f)), 4)
                       for f in EF_FRACTIONS}}

    arms = {"exhaustiveness_32": arm("Vina, exhaustiveness 32", s32),
            "exhaustiveness_8": arm("Vina, exhaustiveness 8", s8),
            "heavy_atoms": arm("heavy-atom count", heavy),
            "clogp": arm("cLogP", clogp)}

    d_lo, d_hi, frac_le0 = paired_delta_ci(y, s32, s8)
    delta = arms["exhaustiveness_32"]["roc_auc"] - arms["exhaustiveness_8"]["roc_auc"]

    verdicts["P1"] = {
        "statement": PREDICTIONS["P1"]["statement"],
        "met": bool(abs(delta) < PREDICTIONS["P1"]["max_abs_delta"]
                    and d_lo <= 0.0 <= d_hi),
        "observed": {"delta_auc_32_minus_8": round(delta, 4),
                     "delta_ci95": [round(d_lo, 4), round(d_hi, 4)],
                     "interval_includes_zero": bool(d_lo <= 0.0 <= d_hi)}}

    t32 = float(a32.elapsed_s.mean())
    t8 = float(a8.elapsed_s.mean())
    ratio = t32 / t8 if t8 else float("nan")
    verdicts["P3"] = {
        "statement": PREDICTIONS["P3"]["statement"],
        "met": bool(PREDICTIONS["P3"]["lo"] <= ratio <= PREDICTIONS["P3"]["hi"]),
        "observed": {"mean_s_exh32": round(t32, 2), "mean_s_exh8": round(t8, 2),
                     "speedup": round(ratio, 2)}}

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "executes": "M2.2",
        "pre_registration": "ADR-AIIM-231, amended before the run",
        "analysis_written": "after the arm finished docking and before any "
                            "vina_score of it was read",
        "n_rows_compared": len(m),
        "n_actives": int((y == 1).sum()), "n_decoys": int((y == 0).sum()),
        "metrics_from": "run_matched_screening_benchmark.py (frozen conventions, "
                        "shared with M5.1 and the 245-compound benchmark)",
        "bootstrap": {"n": BOOTSTRAP_N, "seed": BOOTSTRAP_SEED},
        "arms": arms,
        "primary_outcome": {
            "comparison": "exhaustiveness 32 minus exhaustiveness 8, paired",
            "delta_auc": round(delta, 4),
            "delta_ci95": [round(d_lo, 4), round(d_hi, 4)],
            "interval_includes_zero": bool(d_lo <= 0.0 <= d_hi),
            "bootstrap_fraction_delta_le_zero": round(float(frac_le0), 4)},
        "per_compound_score_difference": {
            "n_moved": n_moved, "fraction_moved": round(n_moved / len(m), 4),
            "mean_abs_kcal_mol": round(float(np.abs(diff).mean()), 4),
            "median_abs_kcal_mol": round(float(np.median(np.abs(diff))), 4),
            "max_abs_kcal_mol": round(float(np.abs(diff).max()), 4)},
        "predictions": verdicts,
        "not_claimed": [
            "that either setting is correct for any target other than this one",
            "that a change here transfers to the natural-product library, whose "
            "background is not property-matched",
            "anything about biological activity; no compound here has been assayed",
        ],
    }
    out.mkdir(parents=True, exist_ok=True)
    (out / "exhaustiveness_arm_analysis.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")
    m[["compound_id", "label_32", "vina_score_32", "vina_score_8"]].to_csv(
        out / "paired_scores.csv", index=False)

    print("=" * 78)
    print("DUD-E A2AR, EXHAUSTIVENESS 32 vs 8   (M2.2)")
    print("=" * 78)
    print(f"  {len(m)} paired rows: {report['n_actives']} active, "
          f"{report['n_decoys']} decoy\n")
    print(f"  {'arm':32s} {'ROC-AUC':>9s} {'95 % CI':>20s} {'BEDROC20':>9s}")
    print("  " + "-" * 74)
    for k in ("exhaustiveness_32", "exhaustiveness_8", "heavy_atoms", "clogp"):
        a = arms[k]
        print(f"  {a['name']:32s} {a['roc_auc']:9.4f} "
              f"[{a['roc_auc_ci95'][0]:.3f}, {a['roc_auc_ci95'][1]:.3f}]".rjust(21)
              + f" {a['bedroc20']:9.4f}")
    print(f"\n  PRIMARY  exhaustiveness 32 minus 8, paired:")
    print(f"           delta = {delta:+.4f}  95 % CI "
          f"[{d_lo:+.4f}, {d_hi:+.4f}]")
    print(f"           interval includes zero: {d_lo <= 0.0 <= d_hi}")
    print(f"\n  predictions recorded before the run:")
    for k in ("P1", "P2", "P3"):
        print(f"    [{'MET ' if verdicts[k]['met'] else 'MISS'}] {k}: "
              f"{verdicts[k]['statement']}")
    print(f"\nwritten -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

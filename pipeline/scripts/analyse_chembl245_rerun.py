"""analyse_chembl245_rerun.py -- the bridge between benchmark and protocol.

§3.4 reports Vina at ROC-AUC 0.601 on the 245-compound ChEMBL A2AR benchmark,
indistinguishable from heavy-atom count. ADR-AIIM-240 reports 0.6689 on DUD-E
A2AR. **Two things differ between those figures — the decoy set and the
conformer protocol — so neither can be blamed on the evidence so far.**

ADR-AIIM-227 re-docked the same 245 compounds under the single-conformer
protocol precisely to remove one of the two. This analyses that run.

    published 245, eight pooled conformers   -> AUC 0.601   (§3.4)
    re-run    245, one minimised conformer   -> measured here
    DUD-E,        one minimised conformer    -> AUC 0.6689  (ADR-AIIM-240)

With the protocol held fixed across the last two, any remaining gap is the
benchmark. With the benchmark held fixed across the first two, any gap is the
protocol. The design is a 2x2 with one cell missing, and the missing cell --
DUD-E under eight pooled conformers -- was priced at 121 hours and declined
(ADR-AIIM-215), so it is named rather than quietly ignored.

### Prediction, recorded before any re-run score was read

**The re-run lands close to the published 0.601, and the paired interval against
it includes zero.** Basis: ADR-AIIM-215 argued the conformer count is largely
redundant for a flexible-ligand docking engine that samples torsions internally,
and ADR-AIIM-241 has since measured that a four-fold change in the *search*
budget moves AUC by 0.0034. If that holds, the 0.601-to-0.6689 gap is the
benchmark, and §3.4's reading stands as written.

**If instead the re-run moves materially toward 0.6689, the protocol is doing
the work**, §3.4's baseline is an artifact of the eight-pooled protocol, and the
manuscript's central negative about this benchmark needs re-examination. That is
the more disruptive outcome and it is written here first so it cannot be
softened afterwards.

Metrics are imported from `run_matched_screening_benchmark.py`, the same
conventions as M5.1 and the matched-screening analysis.

Output -> 20_new_runs/chembl245_rerun_2026-09-22/analysis/
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kdss_paths import under  # noqa: E402

from run_matched_screening_benchmark import (  # noqa: E402
    BEDROC_ALPHA, BOOTSTRAP_N, BOOTSTRAP_SEED, EF_FRACTIONS,
    bedroc, boot_ci, enrichment, paired_delta_ci, roc_auc,
)

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

RERUN = Path(under("20_new_runs", "chembl245_rerun_2026-09-22/scores.csv"))
PUBLISHED = Path(under("20_new_runs", "chembl245_ifp_rescoring_2026-08-15/"
                       "ifp_rescoring_per_compound.csv"))
OUT = Path(under("20_new_runs", "chembl245_rerun_2026-09-22", "analysis"))
N_EXPECTED = 245
PUBLISHED_AUC = 0.601          # §3.4, eight pooled conformers
DUDE_AUC = 0.6689              # ADR-AIIM-240, one minimised conformer


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rerun", default=str(RERUN))
    ap.add_argument("--published", default=str(PUBLISHED))
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--allow-partial", action="store_true")
    args = ap.parse_args()
    out = Path(args.out)

    r = pd.read_csv(args.rerun)
    if len(r) < N_EXPECTED and not args.allow_partial:
        raise SystemExit(f"REFUSING: re-run has {len(r)} of {N_EXPECTED} rows.")
    r["compound_id"] = r.compound_id.astype(str)
    rs = r[r.vina_score.notna()].copy()
    rs["vina_score"] = rs.vina_score.astype(float)
    bad = rs[rs.n_conformers_docked != 1]
    if len(bad):
        raise SystemExit(f"ABORT: {len(bad)} re-run rows not docked with one "
                         f"conformer; the comparison would measure the protocol.")

    p = pd.read_csv(args.published)
    p["compound_id"] = p.compound_id.astype(str)
    ps = p[p.vina_score.notna()][["compound_id", "label", "vina_score"]].copy()

    m = rs.merge(ps, on="compound_id", suffixes=("_new", "_old"))
    if (m.label_new != m.label_old).any():
        raise SystemExit("ABORT: labels disagree between the two runs.")
    if len(m) < 2:
        raise SystemExit("ABORT: too few shared compounds.")
    y = m.label_new.to_numpy().astype(int)
    s_new = -m.vina_score_new.to_numpy()
    s_old = -m.vina_score_old.to_numpy()

    def arm(name, s):
        lo, hi = boot_ci(y, s, roc_auc)
        return {"name": name, "roc_auc": round(float(roc_auc(y, s)), 4),
                "roc_auc_ci95": [round(lo, 4), round(hi, 4)],
                "bedroc20": round(float(bedroc(y, s, BEDROC_ALPHA)), 4),
                "ef": {f"{int(f*100)}%": round(float(enrichment(y, s, f)), 4)
                       for f in EF_FRACTIONS}}

    arms = {"eight_pooled_published": arm("245, eight pooled conformers", s_old),
            "one_conformer_rerun": arm("245, one minimised conformer", s_new)}
    d_lo, d_hi, frac_le0 = paired_delta_ci(y, s_new, s_old)
    delta = arms["one_conformer_rerun"]["roc_auc"] - arms["eight_pooled_published"]["roc_auc"]
    diff = m.vina_score_new.to_numpy() - m.vina_score_old.to_numpy()

    protocol_effect = abs(delta)
    benchmark_gap = abs(DUDE_AUC - arms["one_conformer_rerun"]["roc_auc"])
    verdict = ("benchmark" if benchmark_gap > protocol_effect else "protocol")

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "executes": "the bridge ADR-AIIM-227 created, analysed under ADR-AIIM-242",
        "analysis_written": "before any re-run score was read",
        "n_shared_compounds": len(m),
        "n_actives": int((y == 1).sum()), "n_inactives": int((y == 0).sum()),
        "metrics_from": "run_matched_screening_benchmark.py",
        "bootstrap": {"n": BOOTSTRAP_N, "seed": BOOTSTRAP_SEED},
        "arms": arms,
        "protocol_effect_same_benchmark": {
            "delta_auc_new_minus_old": round(delta, 4),
            "delta_ci95": [round(d_lo, 4), round(d_hi, 4)],
            "interval_includes_zero": bool(d_lo <= 0.0 <= d_hi),
            "bootstrap_fraction_delta_le_zero": round(float(frac_le0), 4)},
        "per_compound_score_difference": {
            "n_moved": int((np.abs(diff) > 1e-9).sum()),
            "mean_abs_kcal_mol": round(float(np.abs(diff).mean()), 4),
            "median_abs_kcal_mol": round(float(np.median(np.abs(diff))), 4),
            "max_abs_kcal_mol": round(float(np.abs(diff).max()), 4)},
        "decomposition": {
            "published_245_eight_pooled": PUBLISHED_AUC,
            "rerun_245_one_conformer": arms["one_conformer_rerun"]["roc_auc"],
            "dude_one_conformer": DUDE_AUC,
            "protocol_effect_abs": round(protocol_effect, 4),
            "benchmark_gap_abs_same_protocol": round(benchmark_gap, 4),
            "dominant_explanation": verdict,
            "missing_cell": "DUD-E under eight pooled conformers, priced at 121 h "
                            "and declined (ADR-AIIM-215)"},
        "not_claimed": [
            "that either protocol is correct in general",
            "anything about biological activity; no compound here has been assayed",
        ],
    }
    out.mkdir(parents=True, exist_ok=True)
    (out / "chembl245_rerun_analysis.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")
    m[["compound_id", "label_new", "vina_score_old", "vina_score_new"]].to_csv(
        out / "paired_scores.csv", index=False)

    print("=" * 74)
    print("ChEMBL 245: EIGHT POOLED vs ONE MINIMISED CONFORMER")
    print("=" * 74)
    print(f"  {len(m)} shared compounds: {report['n_actives']} active, "
          f"{report['n_inactives']} inactive\n")
    for k in ("eight_pooled_published", "one_conformer_rerun"):
        a = arms[k]
        print(f"  {a['name']:34s} AUC {a['roc_auc']:.4f} "
              f"[{a['roc_auc_ci95'][0]:.3f}, {a['roc_auc_ci95'][1]:.3f}]  "
              f"BEDROC20 {a['bedroc20']:.4f}")
    print(f"\n  protocol effect, same benchmark: delta = {delta:+.4f} "
          f"[{d_lo:+.4f}, {d_hi:+.4f}]   includes zero: {d_lo <= 0.0 <= d_hi}")
    print(f"  benchmark gap, same protocol   : "
          f"|{DUDE_AUC} - {arms['one_conformer_rerun']['roc_auc']}| "
          f"= {benchmark_gap:.4f}")
    print(f"\n  dominant explanation for the 0.601 -> 0.6689 difference: "
          f"**{verdict}**")
    print(f"\nwritten -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

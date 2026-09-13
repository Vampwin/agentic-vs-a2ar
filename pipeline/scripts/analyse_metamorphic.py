"""
analyse_metamorphic.py -- evaluate the two metamorphic invariants.

Both variants were docked with the production command; only the input ledger
differed, and it differed in a way that cannot change the answer. Any
deviation is therefore a defect in the pipeline, not chemistry.

The comparison baseline is `rerun_passA_cpu6.csv`, the same 25 compounds
already shown bit-identical to production and to the 20-thread rerun, so a
deviation here cannot be blamed on run-to-run noise: there is none.

The script reports the verdict either way. A metamorphic test that is only
believed when it passes is not a test.

Usage:
    python analyse_metamorphic.py
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from kdss_paths import project_root, resolve, under

import pandas as pd
from scipy import stats


def _root() -> Path:
    # Delegates to kdss_paths, which finds whichever layout this file is inside
    # -- the development tree or a deposit -- by walking up from its own
    # location. It previously tried two absolute paths and exited if neither
    # existed, which on the author's machine meant a copy of the deposit read
    # and wrote the development tree instead of itself.
    return project_root()


ROOT = _root()
G4 = resolve("20_new_runs/g4_verification_2026-08-21")
RUN = resolve("20_new_runs/metamorphic_2026-08-21")
TOL = 0.01          # kcal/mol; the tolerance G4's seed_repeatability check uses

INVARIANT = {
    "M1": "docking score must not depend on the order of rows in the input ledger",
    "M2": "docking score must not depend on which valid SMILES string is used "
          "to write the same molecule",
}


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> int:
    base = pd.read_csv(under("20_new_runs/g4_verification_2026-08-21", "rerun_passA_cpu6.csv")).set_index("compound_id")["vina_score"]
    prod = pd.read_csv(under("20_new_runs/g4_verification_2026-08-21", "full_scores_MERGED.csv")).set_index("compound_id")["vina_score"]
    led = pd.read_csv(under("20_new_runs/g4_verification_2026-08-21", "g4_subset_ledger.csv")).set_index("compound_id")

    out: dict[str, object] = {
        "tolerance_kcal_mol": TOL,
        "baseline": {
            "file": str(under("20_new_runs/g4_verification_2026-08-21", "rerun_passA_cpu6.csv")),
            "sha256": sha256(under("20_new_runs/g4_verification_2026-08-21", "rerun_passA_cpu6.csv")),
            "note": "already established identical to production and to the "
                    "20-thread rerun, so run-to-run noise is zero by measurement",
        },
        "variants": {},
    }

    for v in ("M1", "M2"):
        s = pd.read_csv(RUN / f"scores_{v}.csv").set_index("compound_id")["vina_score"]
        c = s.index.intersection(base.index)
        d = (s.loc[c] - base.loc[c])
        ad = d.abs()
        n_ident = int((ad == 0).sum())
        passed = bool(ad.max() <= TOL)

        rec: dict[str, object] = {
            "invariant": INVARIANT[v],
            "ledger_sha256": sha256(RUN / (
                "ledger_M1_rowshuffle.csv" if v == "M1" else "ledger_M2_smiles_rewrite.csv")),
            "scores_sha256": sha256(RUN / f"scores_{v}.csv"),
            "n_compared": int(len(c)),
            "n_identical": n_ident,
            "max_abs_deviation": round(float(ad.max()), 6),
            "mean_abs_deviation": round(float(ad.mean()), 6),
            "verdict": "HOLDS" if passed else "VIOLATED",
        }

        if not passed:
            r1, r2 = base.loc[c].rank(), s.loc[c].rank()
            shift = (r1 - r2).abs()
            rho_mw, p_mw = stats.spearmanr(led.loc[c, "mw"], ad)
            rho_rb, p_rb = stats.spearmanr(led.loc[c, "rotatable_bonds"], ad)
            spread = float(base.loc[c].max() - base.loc[c].min())
            rec.update({
                "signed_mean_shift": round(float(d.mean()), 6),
                "n_more_negative": int((d < 0).sum()),
                "n_less_negative": int((d > 0).sum()),
                "spearman_rank_correlation_with_baseline":
                    round(float(base.loc[c].corr(s.loc[c], method="spearman")), 6),
                "score_range_of_subset": round(spread, 4),
                "mean_abs_deviation_as_fraction_of_range": round(float(ad.mean()) / spread, 4),
                "rank_movement": {"mean_places": round(float(shift.mean()), 2),
                                  "max_places": int(shift.max()),
                                  "n_moved_3_or_more": int((shift >= 3).sum())},
                "deviation_vs_size": {
                    "spearman_rho_molecular_weight": round(float(rho_mw), 4),
                    "p_molecular_weight": round(float(p_mw), 6),
                    "spearman_rho_rotatable_bonds": round(float(rho_rb), 4),
                    "p_rotatable_bonds": round(float(p_rb), 6)},
                "largest_deviations": [
                    {"compound_id": cid,
                     "baseline": round(float(base[cid]), 3),
                     "variant": round(float(s[cid]), 3),
                     "deviation": round(float(d[cid]), 3),
                     "mw": round(float(led.loc[cid, "mw"]), 1),
                     "rotatable_bonds": int(led.loc[cid, "rotatable_bonds"])}
                    for cid in ad.sort_values(ascending=False).head(6).index],
                "interpretation":
                    "The deviation grows with molecular size and flexibility, "
                    "which is the signature of under-sampled conformer "
                    "generation: ETKDGv3 embedding depends on atom ordering, "
                    "atom ordering is fixed by the SMILES traversal, and eight "
                    "conformers do not sample a large flexible molecule densely "
                    "enough for the choice to wash out. The reported ranking is "
                    "still reproducible -- production used one frozen, "
                    "standardised SMILES per compound and a fixed seed -- but "
                    "the score is a property of (molecule, input "
                    "representation, conformer count), not of the molecule "
                    "alone.",
            })
        out["variants"][v] = rec  # type: ignore[index]

        print(f"--- {v} : {rec['verdict']} ---")
        print(f"  {INVARIANT[v]}")
        print(f"  identical {n_ident}/{len(c)}   max|dev| {ad.max():.4f}   "
              f"mean|dev| {ad.mean():.4f} kcal/mol")
        if not passed:
            print(f"  Spearman vs baseline {rec['spearman_rank_correlation_with_baseline']}")
            print(f"  |dev| vs MW  rho={rec['deviation_vs_size']['spearman_rho_molecular_weight']}"
                  f"  p={rec['deviation_vs_size']['p_molecular_weight']}")
        print()

    # A deviation against production as well as against the rerun, so no one has
    # to take the baseline's equivalence on trust.
    out["cross_check_vs_production"] = {
        v: round(float((pd.read_csv(RUN / f"scores_{v}.csv")
                        .set_index("compound_id")["vina_score"] - prod).abs().max()), 6)
        for v in ("M1", "M2")}

    p = under("20_new_runs/metamorphic_2026-08-21", "metamorphic_results.json")
    p.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"written -> {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

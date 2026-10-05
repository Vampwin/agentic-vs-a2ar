"""tierA_conditioned_null.py -- what the chemical-coherence p-value is testing.

The manuscript reports that the three dual-criterion candidates are all *Pongamia
pinnata* flavones and gives p = 0.002 for that coincidence. The number comes
from `tierA_convergence_analysis.json`, which has no producer script anywhere in
this repository: the reference set it draws from -- 23 compounds -- could not be
matched to any deposited table until this script was written. It is
`conformal_p_active > 0.05` over the 5,000-compound docked library, and this
script reconstructs it and reproduces both frozen values before doing anything
else. If that reconstruction ever stops matching, the assertion at the top of
`main()` stops the run.

Review comment M11 raised two objections. Both are answered here by measurement.

**"The null ignores that the three were selected by docking score."** It does,
and this script shows what difference that makes: none. Of the 23 flagged
compounds, exactly three reached the full protocol, and they are the three best
fast-pass docking scores in the set -- -9.409, -9.391, -9.207, then a gap to
-8.251. So the selection is entirely by docking. But under the null being tested
-- that chemical identity is independent of docking rank -- permuting the labels
across the 23 with the scores held in place leaves the selected positions fixed,
and the permutation distribution is the hypergeometric the frozen artifact
already used. The conditioned p equals the uniform p. That is worth reporting
because it was not obvious in advance, not because it rescues the number.

**"Multiplicity over possible groupings is undeclared."** This is the objection
that bites. The p-value asks whether all three are *Pongamia*, a label chosen
after seeing them. The question a reader hears is whether the three are
chemically coherent at all, and coherence could have been noticed at any of
several levels. Asking "do the three share a label" rather than "are the three
Pongamia" -- summing over every label with at least three members in the
reference set -- is the test that matches the question, and it is computed here
at species, genus and structural-superclass level.

The uniform-null value is kept and labelled, as ADR-AIIM-133 requires. The
conditioned and multiplicity-aware values are added beside it. G8's existing
`tierA_pong` check is not repointed.

Usage:
    python tierA_conditioned_null.py [--out DIR] [--n-permutations N]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from itertools import combinations
from pathlib import Path
from kdss_paths import project_root, resolve, under, dev_dir

import numpy as np
import pandas as pd

# --- pre-registered before the run, transcribed from ADR-AIIM-145 ------------
FLAG_RULE = ("conformal_p_active > 0.05", 0.05)
N_SELECTED = 3                 # the Tier-A count, fixed by the funnel
N_PERMUTATIONS = 200_000
PERMUTATION_SEED = 20260829
# Levels at which "the three are chemically coherent" could have been noticed.
# Fixed here before any of them was computed.
LABEL_LEVELS = ("species", "genus", "superclass")


def _root() -> Path:
    # Delegates to kdss_paths, which finds whichever layout this file is inside
    # -- the development tree or a deposit -- by walking up from its own
    # location. It previously tried two absolute paths and exited if neither
    # existed, which on the author's machine meant a copy of the deposit read
    # and wrote the development tree instead of itself.
    return project_root()


ROOT = _root()
RUNS = dev_dir("20_new_runs")
FUNNEL = under("20_new_runs", "np_library_funnel_2026-08-15")
FROZEN = FUNNEL / "tierA_convergence_analysis.json"


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def load_flagged() -> pd.DataFrame:
    conf = pd.read_csv(FUNNEL / "nplibrary_conformal_predictions.csv", low_memory=False)
    tax = pd.read_csv(FUNNEL / "np_library_dock_ledger.csv", low_memory=False)
    fast = pd.read_csv(FUNNEL / "fastpass_scores_MERGED.csv", low_memory=False)
    e8 = pd.read_csv(FUNNEL / "e8_np_library_prioritization.csv", low_memory=False)

    ids = conf.loc[conf["conformal_p_active"] > FLAG_RULE[1], "compound_id"]
    d = (tax[tax.compound_id.isin(set(ids))]
         [["compound_id", "organism_name",
           "structure_taxonomy_npclassifier_02superclass"]]
         .merge(conf[["compound_id", "conformal_p_active"]], on="compound_id")
         .merge(fast[["compound_id", "vina_score", "on_pocket"]], on="compound_id",
                how="left")
         .merge(e8[["compound_id", "tier"]], on="compound_id", how="left"))
    d = d.rename(columns={"organism_name": "species",
                          "structure_taxonomy_npclassifier_02superclass": "superclass",
                          "vina_score": "fastpass_score"})
    d["genus"] = d["species"].astype(str).str.split().str[0]
    d["selected"] = d["tier"] == "A_converging_sbdd_lbdd"
    return d.sort_values("fastpass_score").reset_index(drop=True)


def hypergeom_all_same(counts: dict[str, int], n: int, k: int) -> float:
    """P(a uniform k-subset of n falls entirely inside one label class)."""
    total = math.comb(n, k)
    return sum(math.comb(c, k) for c in counts.values() if c >= k) / total


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(under("20_new_runs", "tierA_conditioned_null_2026-08-29")))
    ap.add_argument("--n-permutations", type=int, default=N_PERMUTATIONS)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    d = load_flagged()
    n = len(d)
    frozen = json.loads(FROZEN.read_text(encoding="utf-8"))["question_2_is_the_IDENTITY_random"]

    # --- the reconstruction has to match the frozen artifact, or stop ---------
    n_pong = int(d["genus"].eq("Pongamia").sum())
    n_flav = int(d["superclass"].eq("Flavonoids").sum())
    assert n == frozen["library_flagged_total"], (
        f"reconstructed {n} flagged compounds, artifact says "
        f"{frozen['library_flagged_total']}")
    assert n_pong == frozen["pongamia_among_flagged"], (
        f"reconstructed {n_pong} Pongamia, artifact says "
        f"{frozen['pongamia_among_flagged']}")
    assert n_flav == frozen["flavones_among_flagged"], (
        f"reconstructed {n_flav} Flavonoids, artifact says "
        f"{frozen['flavones_among_flagged']}")

    uniform_pong = math.comb(n_pong, N_SELECTED) / math.comb(n, N_SELECTED)
    uniform_flav = math.comb(n_flav, N_SELECTED) / math.comb(n, N_SELECTED)
    assert abs(uniform_pong - frozen["p_all_three_pongamia"]) < 5e-5
    assert abs(uniform_flav - frozen["p_all_three_flavones"]) < 5e-5

    # --- who was selected, and by what ---------------------------------------
    sel = d[d.selected]
    by_score = d.nsmallest(N_SELECTED, "fastpass_score")
    selection_is_top_by_docking = set(sel.compound_id) == set(by_score.compound_id)
    gap = float(d.nsmallest(N_SELECTED + 1, "fastpass_score")["fastpass_score"]
                .iloc[N_SELECTED] - by_score["fastpass_score"].max())

    # --- 1. the docking-conditioned permutation ------------------------------
    # Scores stay with their positions; the chemical labels are permuted. The
    # actual selection rule -- take the N_SELECTED best fast-pass scores -- is
    # re-applied to every permutation.
    rng = np.random.default_rng(PERMUTATION_SEED)
    order = np.argsort(d["fastpass_score"].to_numpy())
    top_positions = order[:N_SELECTED]
    genus = d["genus"].to_numpy()
    hits = 0
    for _ in range(args.n_permutations):
        perm = rng.permutation(genus)
        if perm[top_positions[0]] == "Pongamia" and \
           all(perm[p] == "Pongamia" for p in top_positions[1:]):
            hits += 1
    conditioned_pong = hits / args.n_permutations

    # --- 2. multiplicity: any label, at three levels -------------------------
    multiplicity = {}
    for level in LABEL_LEVELS:
        counts = d[level].value_counts().to_dict()
        eligible = {k: v for k, v in counts.items() if v >= N_SELECTED}
        multiplicity[level] = {
            "n_distinct_labels": len(counts),
            "labels_with_at_least_three_members": eligible,
            "observed_label_of_the_selected": (
                sel[level].iloc[0] if sel[level].nunique() == 1 else None),
            "p_all_three_share_the_observed_label": round(
                math.comb(counts.get(sel[level].iloc[0], 0), N_SELECTED)
                / math.comb(n, N_SELECTED), 6)
            if sel[level].nunique() == 1 else None,
            "p_all_three_share_ANY_label": round(
                hypergeom_all_same(counts, n, N_SELECTED), 6),
        }

    # --- 3. does chemistry actually predict docking here? --------------------
    pong_scores = d.loc[d["genus"] == "Pongamia", "fastpass_score"]
    rest_scores = d.loc[d["genus"] != "Pongamia", "fastpass_score"]
    ranks = d["fastpass_score"].rank().to_numpy()
    pong_ranks = sorted(int(r) for r in ranks[d["genus"].eq("Pongamia").to_numpy()])

    report = {
        "question": "The three convergent candidates are all Pongamia flavones. "
                    "How surprising is that, once it is asked the way a reader "
                    "hears it?",
        "authority": "ADR-AIIM-145",
        "reference_set": {
            "rule": FLAG_RULE[0],
            "over": "the 5,000-compound docked library",
            "n": n,
            "note": "Reconstructed here for the first time. "
                    "tierA_convergence_analysis.json has no producer script and "
                    "its reference set was not identifiable from any deposited "
                    "table until this reconstruction. Both frozen p-values are "
                    "reproduced exactly before anything else runs.",
            "reproduces_frozen_artifact": {
                "file": str(FROZEN.relative_to(ROOT).as_posix()),
                "sha256": sha256(FROZEN),
                "n_flagged": n, "n_pongamia": n_pong, "n_flavonoids": n_flav,
                "p_all_three_pongamia": round(uniform_pong, 6),
                "p_all_three_flavones": round(uniform_flav, 6),
            },
        },
        "how_the_three_were_selected": {
            "n_of_the_flagged_that_reached_full_protocol": int(d.tier.notna().sum()),
            "selected_are_the_top_three_by_fastpass_docking":
                bool(selection_is_top_by_docking),
            "selected_scores": [round(float(x), 3)
                                for x in sorted(sel["fastpass_score"])],
            "gap_to_the_fourth_kcal_per_mol": round(gap, 3),
            "reading": "Selection among the flagged set was by docking score "
                       "alone: only three of the flagged compounds reached the "
                       "full protocol at all, and they are the three best "
                       "fast-pass scores, ahead of the fourth by a clear margin.",
        },
        "uniform_null_as_published": {
            "p_all_three_pongamia": round(uniform_pong, 6),
            "note": "Kept and reported as published (ADR-AIIM-133). It draws "
                    "three uniformly from the flagged set and ignores that the "
                    "three were chosen by docking score.",
        },
        "docking_conditioned_null": {
            "design": "Fast-pass scores stay with their compounds; the chemical "
                      "labels are permuted across the flagged set; the actual "
                      "selection rule -- the three best scores -- is re-applied "
                      "to each permutation.",
            "n_permutations": args.n_permutations,
            "seed": PERMUTATION_SEED,
            "p_all_three_pongamia": round(conditioned_pong, 6),
            "reading": "Indistinguishable from the uniform value. Conditioning "
                       "on docking changes nothing, because under the null being "
                       "tested -- chemical identity independent of docking rank "
                       "-- the selected positions are fixed and permuting labels "
                       "over them reproduces the hypergeometric. The objection "
                       "is answered, and the answer is that it does not matter.",
        },
        "multiplicity": {
            "why": "p = 0.002 tests a label chosen after seeing the compounds. "
                   "A reader hears the broader claim that the three are "
                   "chemically coherent, which could have been noticed at any "
                   "level. Summing over every label with at least three members "
                   "is the test that matches that claim.",
            "levels": multiplicity,
        },
        "does_chemistry_predict_docking_here": {
            "pongamia_fastpass_ranks_within_the_flagged_set": pong_ranks,
            "n_flagged": n,
            "pongamia_median": round(float(pong_scores.median()), 3),
            "other_median": round(float(rest_scores.median()), 3),
            "reading": "The fourth Pongamia compound sits mid-set, so the genus "
                       "is not uniformly better docking; three of its four "
                       "members simply occupy the top three places.",
        },
        "flagged_set": d[["compound_id", "species", "genus", "superclass",
                          "conformal_p_active", "fastpass_score", "on_pocket",
                          "selected"]].to_dict(orient="records"),
    }

    dest = out / "tierA_conditioned_null.json"
    dest.write_text(json.dumps(report, indent=2, default=str), encoding="utf-8")

    print(f"reference set reconstructed: n={n}, Pongamia={n_pong}, "
          f"Flavonoids={n_flav}  (matches the frozen artifact)")
    print(f"selected are the top {N_SELECTED} by fast-pass docking: "
          f"{selection_is_top_by_docking}, gap to fourth {gap:.3f} kcal/mol")
    print()
    print(f"uniform null      p(all three Pongamia) = {uniform_pong:.6f}")
    print(f"conditioned null  p(all three Pongamia) = {conditioned_pong:.6f}  "
          f"({args.n_permutations:,} permutations)")
    print()
    print(f"{'level':<12}{'observed label':<22}{'p(that label)':>14}{'p(ANY label)':>14}")
    for level, m in multiplicity.items():
        po = m["p_all_three_share_the_observed_label"]
        lab = m["observed_label_of_the_selected"]
        print(f"{level:<12}{(lab if lab else 'NOT SHARED'):<22}"
              f"{('--' if po is None else f'{po:.6f}'):>14}"
              f"{m['p_all_three_share_ANY_label']:>14.6f}")
    print(f"\nPongamia fast-pass ranks within the {n}: {pong_ranks}")
    print(f"\nwritten -> {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

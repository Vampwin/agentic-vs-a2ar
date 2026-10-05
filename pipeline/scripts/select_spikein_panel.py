"""select_spikein_panel.py -- freeze the M2.1 spike-in panel before anything runs.

M2.1 asks whether the expansion funnel returns compounds known to be active when
they are placed in the library, and if not, **which stage loses them**. This
script selects the compounds and stops. It docks nothing, screens nothing and
computes no recovery rate; the recovery rate is the run's job, and the panel has
to exist before the run so that the panel cannot be chosen to suit it.

Three decisions are made here, and each is made because a measurement forced it.

1.  **The source is the 169 unseen actives, not all 482.**  ADR-AIIM-222
    measured that 313 of 482 DUD-E A2AR actives already share their InChIKey
    connectivity block with the ChEMBL set the ligand-based arm was fitted on.
    Recovering those would measure memory, not the pipeline.

2.  **The independent unit is a fingerprint cluster, not a Murcko scaffold.**
    The obvious key is the scaffold -- it is the key this project's own
    scaffold-disjoint splits use, and ADR-AIIM-212 reports 482 distinct
    scaffolds for 482 actives. On the 169 it gives 169 clusters of size one,
    which is not a measurement of independence but the absence of one: the key
    cannot see that two compounds differing by a methylene, or by a tautomer,
    are the same series. Morgan/Butina at the standard 0.70 Tanimoto cutoff
    finds 105 clusters, and 97 of the 169 have a nearest neighbour at or above
    that similarity. Sampling 50 compounds as if they were 50 independent tests
    when they are drawn from ~105 series is pseudoreplication, and it would show
    up as a Wilson interval narrower than the evidence supports.

3.  **The panel is not filtered by the pipeline's own criteria.**  Four of the
    169 exceed the funnel's rotatable-bond limit and would be dropped at stage 1.
    They stay in the draw. A positive control pre-filtered by the rule it is
    meant to test cannot test that rule; dropping them would be selecting on the
    outcome, and the stage-1 loss is itself a reportable property of the funnel.

What is frozen is a **seeded random order over the clusters**, not a set. The
primary panel is the first `--n-primary` of that order and the spiked pool holds
the first `--n-spike`. Both are fixed here, before any result exists, so that an
extension from 50 to the full 105 is a pre-registered replication rather than a
choice made after seeing the first 50.

Stages 1 and 2 of the funnel are deterministic functions of frozen inputs, so
what they will do to this panel is **derived here rather than predicted** -- see
`derived_expectations` in the provenance file. That is not a result; it is a
check the run has to reproduce, and a run that disagrees with it has an
implementation fault rather than a finding. The scientific content of M2.1 is at
stages 3 and later, which no derivation reaches.

Output -> 20_new_runs/spikein_panel_2026-09-22/

    spikein_cluster_order_105.csv    the frozen order, all clusters
    spikein_panel_primary_50.csv     the primary panel
    spiked_pool_<n>.csv              the pool the funnel is pointed at
    derived_stage12_expectation.csv  what stages 1-2 must do to the panel
    spikein_panel_provenance.json    seeds, hashes, derived expectations

The first three reproduce bit for bit. The fourth does not, and is kept
separate for that reason: it carries `pharmacophore_2d_score`, and for some
molecules that value is stable within a process and different in the next one,
whatever `PYTHONHASHSEED` is set to (ADR-AIIM-230). One of these 105 is such a
molecule; it sits in the replication arm far below the stage-2 cut, and no
count derived from the column moves. The seed is recorded in the provenance for
the record, not because pinning it makes the column reproducible.

Nothing here may start a run. ADR-AIIM-229 is the pre-registration and the run
waits on its approval.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Descriptors, rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.ML.Cluster import Butina

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kdss_paths import under  # noqa: E402

RDLogger.DisableLog("rdApp.*")
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

SEED = 42

# Morgan radius 2 / 2048 bits is this project's fingerprint everywhere it uses
# one (evaluators/_control_representativeness.py:116). Butina at distance 0.30
# is Tanimoto 0.70, the conventional series cutoff.
FP_RADIUS = 2
FP_BITS = 2048
BUTINA_DIST = 0.30

UNSEEN = Path(under("20_new_runs", "benchmark_overlap_2026-09-22",
                    "dude_actives_unseen_by_training_set.csv"))
LOTUS = Path(under("02_frozen_inputs", "lotus_np_library",
                   "lotus_pool_biosynthetic_25604.csv"))
FUNNEL_META = Path(under("20_new_runs", "np_library_funnel_2026-08-15",
                         "funnel_metadata.json"))
FUNNEL_STAGE2 = Path(under("20_new_runs", "np_library_funnel_2026-08-15",
                           "np_library_stage2_2d_ranked.csv"))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cluster(smiles: list[str]) -> list[tuple[int, ...]]:
    """Butina clusters over Morgan fingerprints; element 0 of each is its centroid."""
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=FP_RADIUS, fpSize=FP_BITS)
    fps = [gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in smiles]
    dists: list[float] = []
    for i in range(len(fps)):
        dists.extend(1.0 - s for s in DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i]))
    return list(Butina.ClusterData(dists, len(fps), distThresh=BUTINA_DIST,
                                   isDistData=True))


def nearest_neighbour_similarity(smiles: list[str]) -> list[float]:
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=FP_RADIUS, fpSize=FP_BITS)
    fps = [gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in smiles]
    out = []
    for i in range(len(fps)):
        others = fps[:i] + fps[i + 1:]
        out.append(max(DataStructs.BulkTanimotoSimilarity(fps[i], others)))
    return out


def derive_stages_1_2(panel: pd.DataFrame) -> dict:
    """What the funnel's first two stages must do to this panel.

    Imported from np_library_prefilter rather than restated, for the same reason
    run_chembl245_rerun.py imports its protocol from run_dude_enrichment.py
    (ADR-AIIM-227): a second copy of a rule is a rule that will drift.
    """
    import np_library_prefilter as P

    meta = json.loads(FUNNEL_META.read_text(encoding="utf-8"))
    P._init_2d_worker(meta["reference_smiles"])

    std, mw, rotb, s2score, standardise_failed = [], [], [], [], []
    for smi in panel.SMILES:
        s, m = P.standardize(smi)
        if s is None:
            standardise_failed.append(smi)
            std.append(None); mw.append(np.nan); rotb.append(np.nan)
            s2score.append(np.nan)
            continue
        std.append(s)
        mw.append(Descriptors.MolWt(m))
        rotb.append(rdMolDescriptors.CalcNumRotatableBonds(m))
        s2score.append(round(P._score_one_2d(s), 4))

    panel = panel.assign(smiles_standardized=std, mw=mw, rotatable_bonds=rotb,
                         pharmacophore_2d_score=s2score)
    survives_1 = (panel.mw <= P.MW_MAX) & (panel.rotatable_bonds <= P.ROTB_MAX)

    # Stage 2 keeps the top 8,000 by score. Every background compound outside
    # the kept 8,000 scores at or below the kept minimum, so the merged top
    # 8,000 is a subset of (kept background + panel) and the cut is exact.
    bg = pd.read_csv(FUNNEL_STAGE2, usecols=["pharmacophore_2d_score"])
    bg_scores = bg.pharmacophore_2d_score.to_numpy()
    n_keep = len(bg_scores)
    merged = np.sort(np.concatenate(
        [bg_scores, panel.loc[survives_1, "pharmacophore_2d_score"].dropna().to_numpy()]
    ))[::-1]
    cut = float(merged[n_keep - 1])
    survives_2 = survives_1 & (panel.pharmacophore_2d_score >= cut)
    ties_at_cut = int((merged[:n_keep] == cut).sum() +
                      (merged[n_keep:] == cut).sum()) if len(merged) > n_keep else 0

    return {
        "panel": panel,
        "summary": {
            "standardisation_failures": standardise_failed,
            "stage1_rule": f"MW <= {P.MW_MAX} and rotatable bonds <= {P.ROTB_MAX}",
            "stage1_survivors": int(survives_1.sum()),
            "stage1_lost": [str(c) for c in panel.loc[~survives_1, "compound_id"]],
            "stage2_background_kept_min_unspiked": float(bg_scores.min()),
            "stage2_cut_with_panel_spiked": cut,
            "stage2_ties_at_cut": ties_at_cut,
            "stage2_survivors": int(survives_2.sum()),
            "panel_score_median": float(panel.pharmacophore_2d_score.median()),
            "background_kept_score_max": float(bg_scores.max()),
            "panel_outranking_entire_background": int(
                (panel.pharmacophore_2d_score > bg_scores.max()).sum()),
            "note": "Derived from frozen inputs, not predicted. Exact up to "
                    "ties at the stage-2 cut, which pandas resolves arbitrarily "
                    "(sort_values uses quicksort).",
        },
    }


def build_spiked_pool(panel: pd.DataFrame) -> pd.DataFrame:
    """LOTUS schema plus a spike_in flag. The funnel selects its columns by name
    and carries the rest through untouched, so no edit to it is required."""
    lot = pd.read_csv(LOTUS)
    lot["spike_in"] = False
    rows = []
    for cid, smi in zip(panel.compound_id.astype(str), panel.SMILES):
        m = Chem.MolFromSmiles(smi)
        rows.append({
            "structure_inchikey": Chem.MolToInchiKey(m),
            "structure_smiles": smi,
            "structure_molecular_formula": rdMolDescriptors.CalcMolFormula(m),
            "structure_exact_mass": Descriptors.ExactMolWt(m),
            # LOTUS carries XLogP; RDKit's Crippen MolLogP substitutes for it.
            # The funnel does not read this column -- it is carried for schema
            # completeness only, and the substitution is recorded rather than
            # hidden.
            "structure_xlogp": round(Descriptors.MolLogP(m), 2),
            "structure_taxonomy_npclassifier_01pathway": "SPIKE_IN",
            "structure_taxonomy_npclassifier_02superclass": "SPIKE_IN",
            "structure_taxonomy_npclassifier_03class": "SPIKE_IN",
            "organism_name": f"SPIKE_IN DUD-E A2AR active {cid} (not a natural product)",
            "organism_taxonomy_06family": "SPIKE_IN",
            "is_orchidaceae": False,
            "spike_in": True,
        })
    spikes = pd.DataFrame(rows)[lot.columns]
    return pd.concat([lot, spikes], ignore_index=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-spike", type=int, default=0,
                    help="clusters to place in the spiked pool; 0 = all")
    ap.add_argument("--n-primary", type=int, default=50,
                    help="primary panel size, the first n of the frozen order")
    ap.add_argument("--no-derive", action="store_true",
                    help="skip the stage-1/2 derivation")
    ap.add_argument("--out", default=str(under("20_new_runs",
                                               "spikein_panel_2026-09-22")))
    args = ap.parse_args()

    out = Path(args.out)
    if out.exists() and any(out.iterdir()):
        print(f"refusing to write into non-empty {out} (20_new_runs is immutable, "
              f"ADR-AIIM-092)")
        return 2
    out.mkdir(parents=True, exist_ok=True)

    src = pd.read_csv(UNSEEN)
    src["compound_id"] = src.compound_id.astype(str)
    mols = [Chem.MolFromSmiles(s) for s in src.SMILES]
    if any(m is None for m in mols):
        print("unparsed SMILES in the source file; refusing to proceed")
        return 2
    src["skeleton_inchikey"] = [Chem.MolToInchiKey(m).split("-")[0] for m in mols]
    src["murcko_scaffold"] = [Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m))
                              for m in mols]
    src["nn_tanimoto"] = np.round(nearest_neighbour_similarity(src.SMILES.tolist()), 4)

    clusters = cluster(src.SMILES.tolist())
    sizes = [len(c) for c in clusters]
    print(f"{len(src)} unseen actives -> {src.murcko_scaffold.nunique()} Murcko "
          f"scaffolds, {len(clusters)} Butina clusters at Tanimoto "
          f"{1 - BUTINA_DIST:.2f} (largest {max(sizes)}, "
          f"{sum(1 for s in sizes if s == 1)} singletons)")

    # One representative per cluster: Butina puts the centroid first.
    rep_idx = [c[0] for c in clusters]
    order = np.random.default_rng(SEED).permutation(len(clusters))
    frozen = src.iloc[[rep_idx[i] for i in order]].copy().reset_index(drop=True)
    frozen.insert(0, "panel_order", range(1, len(frozen) + 1))
    frozen["cluster_size"] = [sizes[i] for i in order]
    frozen["cluster_members"] = [
        ";".join(src.compound_id.iloc[list(clusters[i])]) for i in order]
    frozen["arm"] = ["primary" if i < args.n_primary else "replication"
                     for i in range(len(frozen))]

    n_spike = len(frozen) if args.n_spike in (0, None) else args.n_spike
    panel_all = frozen.head(n_spike).copy()
    panel_primary = frozen.head(args.n_primary).copy()

    derived = None
    if not args.no_derive:
        d = derive_stages_1_2(panel_all)
        derived = d["summary"]
        panel_all = d["panel"]
        cut = derived["stage2_cut_with_panel_spiked"]

        # The deterministic part of the derivation -- the standardised SMILES,
        # molecular weight and rotatable-bond count -- goes onto the frozen
        # selection files. `pharmacophore_2d_score` does NOT: for some molecules
        # it is stable within a process and different in the next one, and
        # pinning PYTHONHASHSEED does not fix it (ADR-AIIM-230). One of these
        # 105 is such a molecule. Putting the column on the frozen artifact
        # would make the artifact irreproducible for a value the selection does
        # not depend on, so it is written to its own labelled file instead.
        frozen = frozen.merge(
            panel_all[["panel_order", "smiles_standardized", "mw",
                       "rotatable_bonds"]], on="panel_order", how="left")
        frozen["derived_survives_stage1"] = (
            (frozen.mw <= 700.0) & (frozen.rotatable_bonds <= 10))
        panel_primary = frozen.head(args.n_primary).copy()

        exp = panel_all[["panel_order", "compound_id", "arm", "mw",
                         "rotatable_bonds", "pharmacophore_2d_score"]].copy()
        exp["survives_stage1"] = (exp.mw <= 700.0) & (exp.rotatable_bonds <= 10)
        exp["survives_stage2"] = (exp.survives_stage1
                                  & (exp.pharmacophore_2d_score >= cut))
        exp_primary = exp.head(args.n_primary)
        derived["primary_arm"] = {
            "n": int(args.n_primary),
            "stage1_survivors": int(exp_primary.survives_stage1.sum()),
            "stage2_survivors": int(exp_primary.survives_stage2.sum()),
        }
        exp.to_csv(out / "derived_stage12_expectation.csv", index=False)
        print(f"derived: stage 1 keeps {derived['stage1_survivors']}/{n_spike}, "
              f"stage 2 keeps {derived['stage2_survivors']}/{n_spike} "
              f"(cut {cut:.4f}); "
              f"{derived['panel_outranking_entire_background']} outrank the "
              f"entire kept background. Primary arm: "
              f"{derived['primary_arm']['stage1_survivors']}/{args.n_primary} and "
              f"{derived['primary_arm']['stage2_survivors']}/{args.n_primary}")

    order_path = out / f"spikein_cluster_order_{len(frozen)}.csv"
    primary_path = out / f"spikein_panel_primary_{args.n_primary}.csv"
    frozen.to_csv(order_path, index=False)
    panel_primary.to_csv(primary_path, index=False)

    pool = build_spiked_pool(panel_all)
    pool_path = out / f"spiked_pool_{len(pool)}.csv"
    pool.to_csv(pool_path, index=False)
    print(f"spiked pool: {len(pool)} rows = {len(pool) - n_spike} LOTUS + "
          f"{n_spike} spike-ins -> {pool_path.name}")

    prov = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "adr": "ADR-AIIM-229",
        "executes": "M2.1 of 29_OBJECTIVES_AND_METHOD_DESIGN_v6.md",
        "source": {"path": str(UNSEEN), "sha256": sha256(UNSEEN),
                   "rows": len(src), "provenance": "ADR-AIIM-222"},
        "background_pool": {"path": str(LOTUS), "sha256": sha256(LOTUS)},
        "selection": {
            "seed": SEED, "rng": "numpy.random.default_rng",
            "python_hash_seed": os.environ.get("PYTHONHASHSEED"),
            "fingerprint": f"Morgan radius {FP_RADIUS}, {FP_BITS} bits",
            "clustering": f"Butina distThresh {BUTINA_DIST} "
                          f"(Tanimoto >= {1 - BUTINA_DIST:.2f})",
            "representative": "Butina cluster centroid (element 0)",
            "n_clusters": len(clusters),
            "n_murcko_scaffolds": int(src.murcko_scaffold.nunique()),
            "cluster_size_max": max(sizes),
            "n_singleton_clusters": sum(1 for s in sizes if s == 1),
            "nn_tanimoto_ge_0.70": int((src.nn_tanimoto >= 0.70).sum()),
            "n_spiked": n_spike,
            "n_primary": args.n_primary,
            "filtered_by_pipeline_criteria": False,
        },
        "outputs": {p.name: sha256(p) for p in
                    (order_path, primary_path, pool_path)},
        "outputs_not_reproducible": {
            p.name: sha256(p) for p in [out / "derived_stage12_expectation.csv"]
            if p.exists()},
        "derived_expectations": derived,
        "not_claimed": [
            "that the panel compounds are active in any assay run by this project",
            "that recovery through stages 1-2 is evidence the pipeline works; "
            "the stage-2 screen ranks by similarity to known A2A ligands and the "
            "panel are known A2A ligands",
            "that the LOTUS background is a property-matched decoy set; this is a "
            "recovery control, not an enrichment benchmark (enrichment is M5.1)",
        ],
    }
    prov_path = out / "spikein_panel_provenance.json"
    prov_path.write_text(json.dumps(prov, indent=2, default=str), encoding="utf-8")
    print(f"wrote {prov_path.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

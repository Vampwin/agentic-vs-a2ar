"""measure_benchmark_training_overlap.py -- how much of DUD-E is already in the training set.

Executes item 1 of the O2-O4 plan, and it exists because of a near-miss.

The spike-in recovery control (M2.1) seeds compounds of known activity into the
screening library and asks whether the pipeline returns them. DUD-E A2AR is the
obvious source of such compounds. It is also the wrong one if its actives are
the same molecules the ligand-based model was fitted on: "recovery" would then
be a memory test, and a favourable result would be an artefact of leakage rather
than evidence about the pipeline.

So the overlap is measured before the control is designed rather than after it
runs.

Three keys, from loosest to strictest:

  SKELETON   InChIKey first block -- the connectivity hash. Two molecules match
             if they are the same constitution, whatever their stereochemistry,
             salt form or charge state. This is the key that matters for a
             fingerprint model: ECFP4 does not see stereochemistry.

  FULL       the complete InChIKey, which separates stereoisomers.

  SCAFFOLD   Bemis-Murcko framework, which matches whole chemical series rather
             than individual compounds, and is the key the project's own
             scaffold-disjoint splits use.

What this measurement does and does not bear on is stated in the artifact, so
that a later reader cannot take it further than it goes:

  * It does NOT affect the M5.1 docking benchmark. AutoDock Vina has no
    training set -- an empirical scoring function fitted once by its authors on
    an unrelated complex set -- so a compound appearing in this project's ChEMBL
    extract confers no advantage on it (ADR-AIIM-212 makes this argument for
    the decoy side; it holds identically here).

  * It DOES rule out using DUD-E A2AR to validate the ligand-based arm, and it
    constrains the spike-in panel to the unseen subset.

Output -> 20_new_runs/benchmark_overlap_2026-09-22/
"""
from __future__ import annotations

import argparse
import io
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kdss_paths import under  # noqa: E402

RDLogger.DisableLog("rdApp.*")

DUDE = Path(under("02_frozen_inputs", "public_benchmarks/dude_a2ar/"
                  "dude_a2ar_docking_set.csv"))
CHEMBL = Path(under("02_frozen_inputs", "public_benchmarks/chembl_a2ar/"
                    "chembl_a2a_raw_v2.csv"))


def keys_for(smiles: str) -> tuple[str, str, str] | None:
    """(skeleton InChIKey, full InChIKey, Murcko scaffold SMILES)."""
    m = Chem.MolFromSmiles(str(smiles))
    if m is None:
        return None
    try:
        ik = Chem.MolToInchiKey(m)
    except Exception:                                        # noqa: BLE001
        return None
    if not ik:
        return None
    try:
        scaf = Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m))
    except Exception:                                        # noqa: BLE001
        scaf = ""
    return ik.split("-")[0], ik, scaf


def index(rows: list[tuple[str, str]]) -> dict:
    """rows = [(identifier, smiles)] -> per-key lookup plus parse failures."""
    out = {"skeleton": {}, "full": {}, "scaffold": {}, "unparsed": []}
    for ident, smi in rows:
        k = keys_for(smi)
        if k is None:
            out["unparsed"].append(ident)
            continue
        sk, fu, sc = k
        out["skeleton"].setdefault(sk, []).append(ident)
        out["full"].setdefault(fu, []).append(ident)
        if sc:
            out["scaffold"].setdefault(sc, []).append(ident)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(
        under("20_new_runs", "benchmark_overlap_2026-09-22")))
    args = ap.parse_args()
    out = Path(args.out)

    d = pd.read_csv(DUDE)
    act = d[d.label == 1]
    a_rows = list(zip(act.compound_id.astype(str), act.SMILES))

    c = pd.read_csv(CHEMBL).drop_duplicates("molecule_chembl_id")
    c_rows = list(zip(c.molecule_chembl_id.astype(str), c.canonical_smiles))

    A, C = index(a_rows), index(c_rows)

    shared = {k: sorted(set(A[k]) & set(C[k])) for k in
              ("skeleton", "full", "scaffold")}

    # A DUD-E active is "seen" if its connectivity appears in the training set.
    seen = {i for k in shared["skeleton"] for i in A["skeleton"][k]}
    all_ids = {i for ids in A["skeleton"].values() for i in ids}
    unseen = sorted(all_ids - seen)

    # For the spike-in panel we also want scaffold diversity among the unseen.
    unseen_set = set(unseen)
    unseen_scaffolds = sorted({sc for sc, ids in A["scaffold"].items()
                               if unseen_set & set(ids)})

    n_act = len(all_ids)
    n_seen = len(seen)

    print("=" * 78)
    print("IS THE BENCHMARK ALREADY IN THE TRAINING SET?")
    print("=" * 78)
    print(f"  DUD-E A2AR actives parsed : {n_act}"
          f"   (unparsed {len(A['unparsed'])})")
    print(f"  ChEMBL A2A molecules      : "
          f"{len({i for ids in C['skeleton'].values() for i in ids})}"
          f"   (unparsed {len(C['unparsed'])})")
    print()
    for k in ("skeleton", "full", "scaffold"):
        n = len(shared[k])
        print(f"  shared {k:9}: {n:5}")
    print()
    print(f"  DUD-E actives whose connectivity is in the training set: "
          f"**{n_seen} of {n_act}**  ({100*n_seen/n_act:.1f} %)")
    print(f"  usable as unseen spike-ins: **{len(unseen)}**, spanning "
          f"{len(unseen_scaffolds)} Murcko scaffolds")
    print()
    print("  Bearing:")
    print("    M5.1 docking benchmark  -- UNAFFECTED. Vina has no training set.")
    print("    validating the LBDD arm -- RULED OUT on this benchmark.")
    print("    M2.1 spike-in panel     -- must be drawn from the unseen subset.")

    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"compound_id": unseen}).merge(
        act[["compound_id", "SMILES"]].astype({"compound_id": str}),
        on="compound_id", how="left").to_csv(
        out / "dude_actives_unseen_by_training_set.csv", index=False)

    summary = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "executes": "item 1 of the O2-O4 plan; prerequisite for M2.1",
        "question": "How many DUD-E A2AR actives are already in the ChEMBL A2A "
                    "set the ligand-based model was fitted on?",
        "sources": {"benchmark": str(DUDE.relative_to(Path(under('', '')))),
                    "training": str(CHEMBL.relative_to(Path(under('', ''))))},
        "n_dude_actives": n_act,
        "n_chembl_molecules": len({i for ids in C["skeleton"].values()
                                   for i in ids}),
        "n_unparsed": {"dude": len(A["unparsed"]), "chembl": len(C["unparsed"])},
        "shared_counts": {k: len(v) for k, v in shared.items()},
        "n_dude_actives_seen_in_training": n_seen,
        "fraction_seen": round(n_seen / n_act, 4),
        "n_unseen": len(unseen),
        "n_unseen_scaffolds": len(unseen_scaffolds),
        "unseen_written_to": "dude_actives_unseen_by_training_set.csv",
        "bearing": {
            "m5_1_docking_benchmark": "unaffected -- AutoDock Vina is an "
                                      "empirical scoring function with no "
                                      "training set, so a compound's presence "
                                      "in this project's ChEMBL extract confers "
                                      "no advantage (same argument as "
                                      "ADR-AIIM-212 makes for the decoys)",
            "validating_the_ligand_based_arm": "ruled out on this benchmark",
            "m2_1_spike_in_panel": "must be drawn from the unseen subset",
        },
        "method": "InChIKey first block (connectivity, stereo-insensitive, "
                  "which is what an ECFP4 model sees), full InChIKey, and "
                  "Bemis-Murcko framework",
    }
    with open(out / "benchmark_overlap.json", "w", encoding="utf-8",
              newline="\n") as f:
        json.dump(summary, f, indent=1, ensure_ascii=False)
        f.write("\n")
    print(f"\nwritten -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

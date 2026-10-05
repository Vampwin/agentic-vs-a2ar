"""verify_spikein_panel.py -- check the frozen M2.1 panel, and measure the
between-process variation that forced one derived column out of it.

Three things are checked, and the script exits non-zero if any fails.

**Structure.** That the frozen order and the primary panel are what they claim:
105 rows and 50 rows, the primary the head of the order, one representative per
Butina cluster with the 169 source compounds partitioned between them, and a
spiked pool that is the LOTUS pool byte-for-byte plus the panel and nothing else.

**Reproducibility.** That re-running `select_spikein_panel.py` under the recorded
`PYTHONHASHSEED` reproduces the frozen selection bit for bit --
`derived_stage12_expectation.csv` excepted, which is excluded by name because it
is the one file carrying the irreproducible column.

**Between-process variation.** `np_library_prefilter._score_one_2d` returns
different values for the same molecule in different processes, and pinning
`PYTHONHASHSEED` does not stop it (ADR-AIIM-230). This re-runs the selection in
several fresh processes, reports which compounds move and by how much, and
confirms that the panel's derived stage-1 and stage-2 counts do not move with
them. That measurement is the evidence for freezing the panel despite the
defect. It is not a repair: the frozen funnel is untouched.

    PYTHONHASHSEED=0 python 10_new_system/scripts/verify_spikein_panel.py

Writes `spikein_panel_verification.json` beside the panel.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kdss_paths import under  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

PANEL_DIR = Path(under("20_new_runs", "spikein_panel_2026-09-22"))
UNSEEN = Path(under("20_new_runs", "benchmark_overlap_2026-09-22",
                    "dude_actives_unseen_by_training_set.csv"))
LOTUS = Path(under("02_frozen_inputs", "lotus_np_library",
                   "lotus_pool_biosynthetic_25604.csv"))
SELECT = Path(__file__).resolve().parent / "select_spikein_panel.py"
SEEDS = ["0", "1", "7", "99", "12345"]

# The seeds below are a way of getting fresh processes, not a controlled
# variable: the variation of ADR-AIIM-230 is between processes, not between
# seeds. This file is excluded from the bit-exactness check by design -- it
# carries the one column the frozen selection does not depend on and cannot
# reproduce.
NOT_REPRO = "derived_stage12_expectation.csv"


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def structural_checks(d: Path) -> list[tuple[str, bool, str]]:
    prov = json.loads((d / "spikein_panel_provenance.json").read_text(encoding="utf-8"))
    n_primary = prov["selection"]["n_primary"]
    order = pd.read_csv(d / f"spikein_cluster_order_{prov['selection']['n_clusters']}.csv")
    prim = pd.read_csv(d / f"spikein_panel_primary_{n_primary}.csv")
    pool = pd.read_csv(d / f"spiked_pool_{25604 + prov['selection']['n_spiked']}.csv")
    src = pd.read_csv(UNSEEN)
    lot = pd.read_csv(LOTUS)

    members = [m for row in order.cluster_members for m in str(row).split(";")]
    out = [
        ("order holds one row per cluster",
         len(order) == prov["selection"]["n_clusters"], f"{len(order)}"),
        ("primary panel is the head of the frozen order",
         list(prim.compound_id) == list(order.compound_id[:n_primary]), ""),
        ("panel compounds are distinct",
         order.compound_id.nunique() == len(order), ""),
        ("every panel compound comes from the unseen set",
         set(order.compound_id.astype(str)) <= set(src.compound_id.astype(str)), ""),
        ("clusters partition the source set",
         sorted(members) == sorted(src.compound_id.astype(str)), f"{len(members)}"),
        ("arms split primary/replication as recorded",
         dict(order.arm.value_counts()).get("primary") == n_primary, ""),
        ("pool = LOTUS + spike-ins, no other row",
         len(pool) == len(lot) + prov["selection"]["n_spiked"], f"{len(pool)}"),
        ("pool schema is LOTUS plus spike_in",
         list(pool.columns) == list(lot.columns) + ["spike_in"], ""),
        ("LOTUS rows carried through unchanged",
         pool[~pool.spike_in][lot.columns].reset_index(drop=True).equals(lot), ""),
        ("spike-ins collide with no LOTUS InChIKey",
         not (set(pool[pool.spike_in].structure_inchikey)
              & set(lot.structure_inchikey)), ""),
        ("spiked SMILES are the panel SMILES, in order",
         list(pool[pool.spike_in].structure_smiles) == list(order.SMILES), ""),
        ("panel was not filtered by the pipeline's own criteria",
         prov["selection"]["filtered_by_pipeline_criteria"] is False, ""),
        ("a hash seed is recorded",
         prov["selection"].get("python_hash_seed") is not None, ""),
    ]
    return out


def rerun(seed: str, tmp: Path) -> Path:
    dest = tmp / f"seed_{seed}"
    env = dict(os.environ, PYTHONHASHSEED=seed)
    r = subprocess.run([sys.executable, str(SELECT), "--out", str(dest)],
                       env=env, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"selection failed at seed {seed}: {r.stdout}{r.stderr}")
    return dest


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel-dir", default=str(PANEL_DIR))
    args = ap.parse_args()
    d = Path(args.panel_dir)
    prov = json.loads((d / "spikein_panel_provenance.json").read_text(encoding="utf-8"))
    pinned = str(prov["selection"]["python_hash_seed"])

    checks = structural_checks(d)
    for name, ok, note in checks:
        print(("  PASS  " if ok else "  FAIL  ") + name + (f"  [{note}]" if note else ""))

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        print(f"\nre-running the selection under hash seeds {', '.join(SEEDS)} ...")
        repro, derived = {}, {}
        for seed in SEEDS:
            dest = rerun(seed, tmp)
            files = {p.name: sha256(p) for p in sorted(dest.iterdir())
                     if p.suffix != ".json" and p.name != NOT_REPRO}
            repro[seed] = files
            derived[seed] = json.loads(
                (dest / "spikein_panel_provenance.json").read_text(encoding="utf-8")
            )["derived_expectations"]

        frozen = {p.name: sha256(p) for p in sorted(d.iterdir())
                  if p.suffix != ".json" and p.name != NOT_REPRO}
        bitexact = repro[pinned] == frozen
        print(("  PASS  " if bitexact else "  FAIL  ")
              + f"frozen outputs reproduce bit for bit at PYTHONHASHSEED={pinned}")

        seed_invariant_files = all(repro[s] == frozen for s in SEEDS)
        keys = ("stage1_survivors", "stage2_survivors", "stage2_cut_with_panel_spiked",
                "panel_outranking_entire_background")
        seed_invariant_counts = all(
            all(derived[s][k] == derived[SEEDS[0]][k] for k in keys) for s in SEEDS)
        seed_invariant_primary = all(
            derived[s]["primary_arm"] == derived[SEEDS[0]]["primary_arm"] for s in SEEDS)
        print(("  PASS  " if seed_invariant_counts and seed_invariant_primary
               else "  FAIL  ")
              + "derived stage-1/stage-2 counts are the same in every process")
        print(f"  note  frozen selection identical under every seed: "
              f"{seed_invariant_files} ({NOT_REPRO} excluded by name -- it is "
              f"where the irreproducible column lives, ADR-AIIM-230)")

        # Where the seed does bite: re-score the 169 candidates under each seed.
        moved = {}
        for seed in SEEDS:
            dest = tmp / f"seed_{seed}"
            o = pd.read_csv(dest / NOT_REPRO)
            moved[seed] = dict(zip(o.compound_id.astype(str),
                                   o.pharmacophore_2d_score))
        ids = list(moved[SEEDS[0]])
        spreads = {i: max(moved[s][i] for s in SEEDS) - min(moved[s][i] for s in SEEDS)
                   for i in ids}
        sensitive = {i: round(v, 6) for i, v in spreads.items() if v > 1e-9}
        print(f"\nbetween-process variation over the {len(ids)} panel compounds: "
              f"{len(sensitive)} move, max spread "
              f"{max(spreads.values()):.4f}")
        for i, v in sorted(sensitive.items(), key=lambda kv: -kv[1]):
            print(f"   compound {i}: spread {v:.4f}")

    ok = all(c[1] for c in checks) and bitexact and seed_invariant_counts \
        and seed_invariant_primary
    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "adr": ["ADR-AIIM-229", "ADR-AIIM-230"],
        "panel_dir": str(d),
        "pinned_python_hash_seed": pinned,
        "structural_checks": [{"check": n, "pass": bool(v)} for n, v, _ in checks],
        "bit_exact_at_pinned_seed": bool(bitexact),
        "frozen_selection_identical_in_every_process": bool(seed_invariant_files),
        "derived_counts_invariant": bool(seed_invariant_counts
                                         and seed_invariant_primary),
        "seeds_tested": SEEDS,
        "compounds_varying_between_processes": sensitive,
        "max_between_process_spread": round(max(spreads.values()), 6),
        "overall": "PASS" if ok else "FAIL",
    }
    (d / "spikein_panel_verification.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")
    print(f"\noverall: {report['overall']}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())

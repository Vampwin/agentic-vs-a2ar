"""run_spikein_dock.py -- dock the M2.1 spiked shortlist under the unified protocol.

M2.1 (ADR-AIIM-229) spikes 105 known A2A actives into the natural-product pool,
runs the expansion funnel unchanged, and asks which stage loses them. The funnel
produces a shortlist; this docks it. **It does not analyse anything** -- no
recovery rate, no enrichment, no ranking claim. It docks and records, for the
reason ADR-AIIM-227 gives: an analysis that lives inside the run that produced
its data cannot be shown to have been written first.

### The protocol is imported, not restated

Every setting, and `dock_one` itself, comes from `run_dude_enrichment.py`:

    8 ETKDGv3 conformers generated, all MMFF94-minimised,
    the LOWEST-ENERGY one docked, exhaustiveness 32, 10 poses, seed 20260729,
    receptor 4EIY_prepared.pdbqt, box centre [-0.424, 8.527, 17.127], 22 A cube

That box and receptor are byte-identical to `snp_screen_2026-08-04/
G2_preparation/docking_box.json`, the production geometry section 3.5 and
section 3.6 were produced with -- checked, not assumed. So one import statement
makes this run, the DUD-E benchmark, the ChEMBL 245 re-run and the M2.2 arm share
a protocol, and a reviewer verifies that by reading the import rather than by
diffing four files.

`g4_dock_library.py` is deliberately NOT reused. It embeds N conformers, docks
**every one** and pools by energy -- the eight-conformer production protocol --
and `--n-conformers 1` would give a third protocol again, a single unminimised
conformer. Mixing any of them measures the protocol, not recovery.

### Which environment

WSL `diffdock`, because vina lives only there and because the DUD-E run, the
ChEMBL 245 re-run and the M2.2 arm all ran there -- the comparators this run has
to match. The *funnel* is pinned the other way, to Windows rdkit 2026.03.3, whose
stage-3 O3A reproduces the published artifact 40/40 against the diffdock
environment's 1/40 (ADR-AIIM-232). Ranking on Windows, docking in WSL, each
consistent with its own comparator, and both facts written down because neither
is guessable.

### What it refuses

* a shortlist whose key is not unique -- the lesson of ADR-AIIM-233, where a
  benchmark with 2,402 rows and 2,401 distinct ids silently became 2,404
* a shortlist missing a column it needs, rather than a KeyError halfway through
* spike-ins in the shortlist that are not in the frozen panel of ADR-AIIM-229
* on completion, any scored row not docked with exactly one conformer

### Usage

    # measure the rate first -- natural products are larger than DUD-E ligands
    # and the 2026-08-15 fast pass measured 40.0 s/compound under the old
    # protocol. A seeded SAMPLE, not the head: the head is the top of the funnel
    # ranking. The probe writes into the real output directory, so its 30 docks
    # are the first 30 of the run rather than work thrown away.
    python run_spikein_dock.py --sample 30 --out 20_new_runs/spike_in_dock_<date> --cpu 20

    # the run, resumable: the same command picks up where it stopped
    python run_spikein_dock.py --out 20_new_runs/spike_in_dock_<date> --cpu 20

    # exercise every code path without vina and without docking anything
    python run_spikein_dock.py --self-test
"""
from __future__ import annotations

import argparse
import csv
import io
import os
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kdss_paths import under  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

FIELDS = ["inchikey_standardized", "structure_inchikey", "spike_in", "funnel_rank",
          "vina_score", "n_conformers_docked", "elapsed_s", "error"]
NEEDED = ["inchikey_standardized", "structure_inchikey", "smiles_standardized",
          "funnel_rank"]

DEFAULT_SHORTLIST = under("20_new_runs", "spike_in_recovery_2026-09-23",
                          "np_library_docking_shortlist.csv")
DEFAULT_POOL = under("20_new_runs", "spikein_panel_2026-09-22",
                     "spiked_pool_25709.csv")
DEFAULT_PANEL = under("20_new_runs", "spikein_panel_2026-09-22",
                      "spikein_cluster_order_105.csv")


def done_keys(path: Path) -> set[str]:
    if not path.exists():
        return set()
    try:
        return set(pd.read_csv(path).inchikey_standardized.astype(str))
    except Exception:                                        # noqa: BLE001
        return set()


def load_inputs(shortlist: Path, pool: Path, panel: Path) -> pd.DataFrame:
    """Read the shortlist, mark the spike-ins, and refuse anything malformed."""
    d = pd.read_csv(shortlist)
    missing = [c for c in NEEDED if c not in d.columns]
    if missing:
        raise SystemExit(f"REFUSING: shortlist lacks {missing}")

    d["inchikey_standardized"] = d.inchikey_standardized.astype(str)
    d["structure_inchikey"] = d.structure_inchikey.astype(str)
    dup = d.inchikey_standardized[d.inchikey_standardized.duplicated()].unique()
    if len(dup):
        raise SystemExit(
            f"REFUSING: {len(dup)} shortlist keys repeat, e.g. {list(dup[:3])}. "
            f"A resumable run keyed on a non-unique column loses compounds "
            f"silently (ADR-AIIM-233).")

    p = pd.read_csv(pool, usecols=["structure_inchikey", "spike_in"])
    spikes = set(p.loc[p.spike_in == True, "structure_inchikey"].astype(str))  # noqa: E712
    d["spike_in"] = d.structure_inchikey.isin(spikes)

    # Every spike-in that survived the funnel must be one this project froze.
    frozen = pd.read_csv(panel)
    n_frozen = len(frozen)
    found = int(d.spike_in.sum())
    if found > n_frozen:
        raise SystemExit(f"REFUSING: {found} spike-ins in the shortlist but the "
                         f"frozen panel holds {n_frozen} (ADR-AIIM-229)")
    print(f"shortlist {len(d)} compounds; {found} of the frozen {n_frozen} "
          f"spike-ins survived the funnel to this stage")
    return d


def run(d: pd.DataFrame, out: Path, cpu: int, limit: int, dock,
        sample: int = 0, seed: int = 42) -> int:
    out.mkdir(parents=True, exist_ok=True)
    scores = out / "scores.csv"
    if sample:
        # For a rate probe. The head of the shortlist is the top of the funnel
        # ranking, not a random draw from it, and docking cost tracks size and
        # flexibility -- so timing the first N would estimate the rate of the
        # most A2A-like end of the library rather than of the library. Seeded,
        # so the probe is reproducible, and ordered by rank afterwards so the
        # file stays readable.
        d = d.sample(min(sample, len(d)), random_state=seed).sort_values("funnel_rank")
    elif limit:
        d = d.head(limit)
    already = done_keys(scores)
    todo = d[~d.inchikey_standardized.isin(already)]
    print(f"assigned {len(d)}  done {len(already)}  to do {len(todo)}  cpu={cpu}",
          flush=True)
    if todo.empty:
        print("nothing to do", flush=True)
        return 0

    v = dock.setup(cpu)
    new = not scores.exists()
    fh = io.open(scores, "a", newline="", encoding="utf-8")
    w = csv.DictWriter(fh, fieldnames=FIELDS)
    if new:
        w.writeheader()
        fh.flush()

    t0 = time.time()
    for i, (_, r) in enumerate(todo.iterrows(), 1):
        c0 = time.time()
        score, n, err = dock.one(v, r.smiles_standardized)
        w.writerow({"inchikey_standardized": r.inchikey_standardized,
                    "structure_inchikey": r.structure_inchikey,
                    "spike_in": bool(r.spike_in),
                    "funnel_rank": int(r.funnel_rank),
                    "vina_score": "" if score is None else round(score, 4),
                    "n_conformers_docked": n,
                    "elapsed_s": round(time.time() - c0, 2),
                    "error": err or ""})
        fh.flush()
        os.fsync(fh.fileno())
        if i % 25 == 0 or i == len(todo):
            rate = (time.time() - t0) / i
            print(f"  {i}/{len(todo)}  {rate:.1f}s/cpd  "
                  f"eta {(len(todo)-i)*rate/3600:.1f}h", flush=True)
    fh.close()

    # The protocol gate, after the fact rather than as a promise.
    s = pd.read_csv(scores)
    scored = s[s.vina_score.notna() & (s.vina_score != "")]
    bad = scored[scored.n_conformers_docked != 1]
    if len(bad):
        print(f"FAIL: {len(bad)} rows docked with a conformer count other than 1; "
              f"this run cannot be compared with any other protocol.", flush=True)
        return 2
    print(f"complete: {len(scored)} scored, {len(s) - len(scored)} failed, "
          f"{int(scored.spike_in.sum())} of them spike-ins, in "
          f"{(time.time()-t0)/3600:.2f}h", flush=True)
    return 0


class RealDock:
    """The protocol, imported. Nothing here restates a setting."""

    def __init__(self) -> None:
        import run_dude_enrichment as D
        self.D = D
        print(f"protocol imported from run_dude_enrichment: "
              f"{D.N_CONFORMERS_GENERATED} conformers generated, "
              f"{D.N_CONFORMERS_DOCKED} docked, exhaustiveness {D.EXHAUSTIVENESS}, "
              f"{D.N_POSES} poses, seed {D.SEED}")
        print(f"receptor {D.RECEPTOR.name}, box centre {D.CENTER}, size {D.SIZE}")

    def setup(self, cpu: int):
        from vina import Vina
        v = Vina(sf_name="vina", cpu=cpu, seed=self.D.SEED, verbosity=0)
        v.set_receptor(str(self.D.RECEPTOR))
        t = time.time()
        v.compute_vina_maps(center=self.D.CENTER, box_size=self.D.SIZE)
        print(f"maps computed in {time.time()-t:.1f}s", flush=True)
        return v

    def one(self, v, smiles: str):
        return self.D.dock_one(v, smiles)


class FakeDock:
    """Deterministic stand-in for the self-test. Docks nothing."""

    def setup(self, cpu: int):
        return None

    def one(self, v, smiles: str):
        if smiles == "BAD":
            return None, 0, "unparseable SMILES"
        if smiles == "TWO":
            return -7.0, 2, None                 # wrong protocol, must be caught
        return -6.0 - (len(smiles) % 7) / 10.0, 1, None


def self_test() -> int:
    """Exercise every path without vina, rdkit or a single dock."""
    import json
    import tempfile
    ok = []

    def check(name, cond):
        ok.append((name, bool(cond)))
        print(("  PASS  " if cond else "  FAIL  ") + name)

    with tempfile.TemporaryDirectory() as td:
        t = Path(td)
        n = 40
        short = pd.DataFrame({
            "inchikey_standardized": [f"KEY{i:04d}-AAAAAA-N" for i in range(n)],
            "structure_inchikey": [f"SRC{i:04d}-AAAAAA-N" for i in range(n)],
            "smiles_standardized": [f"C{'C'*(i%9)}O" for i in range(n)],
            "funnel_rank": range(1, n + 1)})
        short.to_csv(t / "short.csv", index=False)
        pool = pd.DataFrame({
            "structure_inchikey": [f"SRC{i:04d}-AAAAAA-N" for i in range(n)],
            "spike_in": [i % 10 == 0 for i in range(n)]})
        pool.to_csv(t / "pool.csv", index=False)
        pd.DataFrame({"compound_id": range(105)}).to_csv(t / "panel.csv", index=False)

        d = load_inputs(t / "short.csv", t / "pool.csv", t / "panel.csv")
        check("shortlist loads and flags spike-ins", int(d.spike_in.sum()) == 4)

        rc = run(d, t / "run", cpu=1, limit=20, dock=FakeDock())
        s1 = pd.read_csv(t / "run" / "scores.csv")
        check("limit honoured, rows written", rc == 0 and len(s1) == 20)
        check("schema is exactly FIELDS", list(s1.columns) == FIELDS)

        rc = run(d, t / "run", cpu=1, limit=0, dock=FakeDock())
        s2 = pd.read_csv(t / "run" / "scores.csv")
        check("resume completes without duplicating",
              rc == 0 and len(s2) == n and s2.inchikey_standardized.is_unique)
        check("the first 20 rows were not re-docked",
              list(s2.inchikey_standardized[:20]) == list(s1.inchikey_standardized))
        check("spike-ins carried through", int(s2.spike_in.sum()) == 4)

        # protocol gate
        bad = short.copy()
        bad["smiles_standardized"] = "TWO"
        bad.to_csv(t / "bad.csv", index=False)
        db = load_inputs(t / "bad.csv", t / "pool.csv", t / "panel.csv")
        rc = run(db, t / "run_bad", cpu=1, limit=5, dock=FakeDock())
        check("two-conformer rows are refused", rc == 2)

        # failed docks are recorded, not fatal
        f = short.copy()
        f.loc[0, "smiles_standardized"] = "BAD"
        f.to_csv(t / "f.csv", index=False)
        df = load_inputs(t / "f.csv", t / "pool.csv", t / "panel.csv")
        rc = run(df, t / "run_f", cpu=1, limit=5, dock=FakeDock())
        sf = pd.read_csv(t / "run_f" / "scores.csv")
        check("a failed dock is recorded and does not stop the run",
              rc == 0 and sf.vina_score.isna().sum() == 1)

        # refusals
        dup = pd.concat([short, short.iloc[[0]]], ignore_index=True)
        dup.to_csv(t / "dup.csv", index=False)
        try:
            load_inputs(t / "dup.csv", t / "pool.csv", t / "panel.csv")
            check("duplicate shortlist key refused", False)
        except SystemExit:
            check("duplicate shortlist key refused", True)

        short.drop(columns=["funnel_rank"]).to_csv(t / "miss.csv", index=False)
        try:
            load_inputs(t / "miss.csv", t / "pool.csv", t / "panel.csv")
            check("missing column refused", False)
        except SystemExit:
            check("missing column refused", True)

        pd.DataFrame({"compound_id": range(2)}).to_csv(t / "small.csv", index=False)
        try:
            load_inputs(t / "short.csv", t / "pool.csv", t / "small.csv")
            check("more spike-ins than the frozen panel refused", False)
        except SystemExit:
            check("more spike-ins than the frozen panel refused", True)

    print(f"\n{sum(v for _, v in ok)} of {len(ok)} checks pass")
    return 0 if all(v for _, v in ok) else 1


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--shortlist", default=str(DEFAULT_SHORTLIST))
    ap.add_argument("--pool", default=str(DEFAULT_POOL))
    ap.add_argument("--panel", default=str(DEFAULT_PANEL))
    ap.add_argument("--out", default=str(under("20_new_runs",
                                               "spike_in_dock_2026-09-23")))
    ap.add_argument("--cpu", type=int, default=20)
    ap.add_argument("--limit", type=int, default=0,
                    help="dock only the first N by funnel rank")
    ap.add_argument("--sample", type=int, default=0,
                    help="dock a seeded random sample of N; this is the rate "
                         "probe, and it is a sample rather than a head because "
                         "docking cost tracks molecular size")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        return self_test()

    d = load_inputs(Path(args.shortlist), Path(args.pool), Path(args.panel))
    return run(d, Path(args.out), args.cpu, args.limit, RealDock(),
               sample=args.sample, seed=args.seed)


if __name__ == "__main__":
    raise SystemExit(main())

"""measure_funnel_environment_sensitivity.py -- the artifact ADR-AIIM-232 lacked.

ADR-AIIM-232 recorded that the expansion funnel's 3D stage reproduces the
published artifact under one interpreter and almost never under the other, while
its 2D stage reproduces under both. ADR-AIIM-249 then kept that finding out of
the manuscript for a reason it states plainly: **the numbers existed only inside
the ADR, with no artifact to anchor against**, and writing an unanchored number
into a paper arguing that numbers must be traceable would contradict the paper
in its own section 3.7.

This builds the artifact. It re-scores compounds drawn from the **published**
`np_library_stage2_2d_ranked.csv` with the funnel's own scoring functions and
compares against the values that file already stores.

    stage 2   `_score_one_2d`    Gobbi 2D pharmacophore + ErG, topological
    stage 3   `_score_one_o3a`   3D embed, MMFF94 minimise, O3A alignment

Run once per environment; the results merge on the `environment` key:

    # Windows, rdkit 2026.03.3
    python measure_funnel_environment_sensitivity.py --tag windows

    # WSL diffdock, rdkit 2025.09.5
    python measure_funnel_environment_sensitivity.py --tag wsl_diffdock

    # then
    python measure_funnel_environment_sensitivity.py --merge

### Why this is cheap and was mis-costed

ADR-AIIM-249 deferred it "only because the CPU is committed". O3A was measured
at **0.34 s per molecule** on one core (ADR-AIIM-232's own cost probe), so 40
compounds in two environments is about **thirty seconds of a single core** and
competes with nothing. The deferral was a misjudgement of the cost, not of the
value, and this file exists to correct it.

### What it does not do

It scores. It draws no conclusion about which interpreter is correct — both
reproduce stage 2, so neither is broken, and the published artifact was produced
by whichever one produced it. The finding is that the **3D stage is not portable
between them**, which is a statement about the pipeline's reproducibility
envelope and not about rdkit.

Output -> 20_new_runs/funnel_environment_sensitivity_2026-09-24/
"""
from __future__ import annotations

import argparse
import io
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kdss_paths import under  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

PUBLISHED = Path(under("20_new_runs", "np_library_funnel_2026-08-15",
                       "np_library_stage2_2d_ranked.csv"))
FUNNEL_META = Path(under("20_new_runs", "np_library_funnel_2026-08-15",
                         "funnel_metadata.json"))
OUT = Path(under("20_new_runs", "funnel_environment_sensitivity_2026-09-24"))

# Fixed so every environment scores the identical compounds, and so the numbers
# are the ones ADR-AIIM-232 reported rather than a fresh draw.
N_STAGE2, SEED_STAGE2 = 300, 0
N_STAGE3, SEED_STAGE3 = 40, 1
EXACT = 1e-9


def measure(tag: str) -> dict:
    import rdkit
    import np_library_prefilter as P
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")

    meta = json.loads(FUNNEL_META.read_text(encoding="utf-8"))
    refs = meta["reference_smiles"]
    pub = pd.read_csv(PUBLISHED, usecols=["smiles_standardized",
                                          "pharmacophore_2d_score",
                                          "shape_pharm_3d_score"])

    # ---- stage 2 -------------------------------------------------------
    P._init_2d_worker(refs)
    s2 = pub.sample(N_STAGE2, random_state=SEED_STAGE2).reset_index(drop=True)
    got2 = np.round([P._score_one_2d(s) for s in s2.smiles_standardized], 4)
    d2 = np.abs(got2 - s2.pharmacophore_2d_score.to_numpy())

    # ---- stage 3 -------------------------------------------------------
    refs_3d = [P.ZM241385_SMILES, P.get_9xt_smiles()] + refs[2:12]
    P._init_o3a_worker(refs_3d, P.SEED)
    s3 = pub.sample(N_STAGE3, random_state=SEED_STAGE3).reset_index(drop=True)
    got3 = np.round([P._score_one_o3a((s, P.SEED)) for s in s3.smiles_standardized], 4)
    d3 = np.abs(got3 - s3.shape_pharm_3d_score.to_numpy())

    return {
        "environment": tag,
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "interpreter": {"python": platform.python_version(),
                        "rdkit": rdkit.__version__,
                        "platform": platform.system()},
        "published_artifact": str(PUBLISHED),
        "stage2_2d": {
            "n": int(len(s2)), "sample_seed": SEED_STAGE2,
            "reproduced_exactly": int((d2 < EXACT).sum()),
            "max_abs_difference": round(float(d2.max()), 6),
            "mean_abs_difference": round(float(d2.mean()), 6)},
        "stage3_o3a": {
            "n": int(len(s3)), "sample_seed": SEED_STAGE3,
            "reproduced_exactly": int((d3 < EXACT).sum()),
            "max_abs_difference": round(float(d3.max()), 4),
            "mean_abs_difference": round(float(d3.mean()), 4)},
    }


def merge() -> int:
    files = sorted(OUT.glob("environment_*.json"))
    if len(files) < 2:
        print(f"REFUSING to merge: {len(files)} environment file(s) present, "
              f"need at least 2. Run this script in each environment first.")
        return 2
    envs = [json.loads(p.read_text(encoding="utf-8")) for p in files]
    rows = {e["environment"]: e for e in envs}

    s3 = {k: v["stage3_o3a"]["reproduced_exactly"] for k, v in rows.items()}
    s2 = {k: v["stage2_2d"]["reproduced_exactly"] for k, v in rows.items()}
    n3 = envs[0]["stage3_o3a"]["n"]
    n2 = envs[0]["stage2_2d"]["n"]

    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "adr": ["ADR-AIIM-232", "ADR-AIIM-249", "ADR-AIIM-250"],
        "question": "Does the expansion funnel reproduce its published scores "
                    "under a different interpreter?",
        "environments": {k: v["interpreter"] for k, v in rows.items()},
        "stage2_2d_reproduced": {"n": n2, "by_environment": s2,
                                 "portable": all(v == n2 for v in s2.values())},
        "stage3_o3a_reproduced": {"n": n3, "by_environment": s3,
                                  "portable": all(v == n3 for v in s3.values())},
        "per_environment": rows,
        "not_claimed": [
            "that either interpreter is wrong; both reproduce the topological "
            "stage, so neither is broken",
            "anything about rdkit as a library; this measures this pipeline's "
            "reproducibility envelope",
        ],
    }
    (OUT / "funnel_environment_sensitivity.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")

    print("=" * 70)
    print("FUNNEL REPRODUCIBILITY ACROSS INTERPRETERS")
    print("=" * 70)
    for k, v in rows.items():
        i = v["interpreter"]
        print(f"  {k:16s} python {i['python']:9s} rdkit {i['rdkit']}")
    print()
    print(f"  stage 2 (2D pharmacophore/ErG), n = {n2}")
    for k, v in s2.items():
        print(f"    {k:16s} reproduced exactly {v:4d} / {n2}")
    print(f"  stage 3 (3D O3A), n = {n3}")
    for k, v in s3.items():
        print(f"    {k:16s} reproduced exactly {v:4d} / {n3}   "
              f"max |diff| {rows[k]['stage3_o3a']['max_abs_difference']}")
    print(f"\n  stage 2 portable: {report['stage2_2d_reproduced']['portable']}   "
          f"stage 3 portable: {report['stage3_o3a_reproduced']['portable']}")
    print(f"\nwrote {OUT / 'funnel_environment_sensitivity.json'}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", help="environment label, e.g. windows or wsl_diffdock")
    ap.add_argument("--merge", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    if args.merge:
        return merge()
    if not args.tag:
        print("give --tag <name> or --merge")
        return 2

    r = measure(args.tag)
    (OUT / f"environment_{args.tag}.json").write_text(
        json.dumps(r, indent=2), encoding="utf-8")
    i = r["interpreter"]
    print(f"{args.tag}: python {i['python']} rdkit {i['rdkit']}")
    print(f"  stage 2 reproduced {r['stage2_2d']['reproduced_exactly']}"
          f"/{r['stage2_2d']['n']}")
    print(f"  stage 3 reproduced {r['stage3_o3a']['reproduced_exactly']}"
          f"/{r['stage3_o3a']['n']}  max |diff| "
          f"{r['stage3_o3a']['max_abs_difference']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

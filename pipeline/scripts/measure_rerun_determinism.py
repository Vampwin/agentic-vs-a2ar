"""measure_rerun_determinism.py -- is the docking stage reproducible between runs?

§3.1 justifies contracts over snapshot testing by saying the production
pipeline's outputs vary between runs. That is checkable, and this project had
already produced the artifacts to check it without anyone computing the
comparison: `20_new_runs/g4_verification_2026-08-21/` holds the merged
production scores and two independent reruns of the same gate at **different
core counts**, which is the setting where a non-deterministic search or a
thread-order-dependent reduction would show itself.

This reads those three files and reports the agreement. Nothing is re-docked.

The answer matters in both directions. If the scores drift, snapshot testing is
unusable on this pipeline and the manuscript's justification stands as written.
If they do not drift, the justification needs to name the artifacts that do vary
instead of implying the numbers do -- and the pipeline has a reproducibility
property worth reporting in its own right, since bitwise agreement across core
counts is not something most computational workflows can show.

Output -> 20_new_runs/rerun_determinism_2026-09-22/
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kdss_paths import under  # noqa: E402

SRC = under("20_new_runs", "g4_verification_2026-08-21")
FILES = {
    "production_merged": "full_scores_MERGED.csv",
    "rerun_passA_cpu6": "rerun_passA_cpu6.csv",
    "rerun_passB_cpu20": "rerun_passB_cpu20.csv",
}
# Agreement is reported at the precision the pipeline writes, and also exactly.
REPORTED_DECIMALS = 3


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(under("20_new_runs",
                                               "rerun_determinism_2026-09-22")))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    series, meta = {}, {}
    for tag, fn in FILES.items():
        p = Path(SRC) / fn
        d = pd.read_csv(p)
        series[tag] = d.set_index("compound_id")["vina_score"]
        meta[tag] = {"file": fn, "rows": int(len(d)), "sha256": sha256(p)}

    pairs = []
    for a, b in combinations(series, 2):
        sa, sb = series[a].align(series[b], join="inner")
        diff = (sa - sb).abs()
        pairs.append({
            "a": a, "b": b,
            "n_common_compounds": int(len(sa)),
            "n_values_differing_exactly": int((diff > 0).sum()),
            "max_absolute_difference": float(diff.max()) if len(diff) else None,
            "mean_absolute_difference": float(diff.mean()) if len(diff) else None,
            "identical_at_reported_precision": bool(
                (diff.round(REPORTED_DECIMALS) == 0).all()),
            "identical_exactly": bool((diff == 0).all()),
        })

    all_exact = all(p["identical_exactly"] for p in pairs)
    n_common = min(p["n_common_compounds"] for p in pairs)

    summary = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "question": "Do the docking scores of this pipeline change between "
                    "independent runs at different core counts?",
        "nothing_redocked": True,
        "source_run": "20_new_runs/g4_verification_2026-08-21",
        "inputs": meta,
        "reported_decimals": REPORTED_DECIMALS,
        "pairwise": pairs,
        "n_compounds_common_to_all": int(n_common),
        "all_pairs_bitwise_identical": bool(all_exact),
        "max_absolute_difference_over_all_pairs": max(
            (p["max_absolute_difference"] or 0.0) for p in pairs),
        "interpretation": (
            "The docking scores are reproducible; the variation that does exist "
            "between runs is in other artifacts, which a separate comparison "
            "measures." if all_exact else
            "The docking scores vary between runs, so snapshot comparison is "
            "unusable on this pipeline as the manuscript states."),
    }
    with open(out / "rerun_determinism.json", "w", encoding="utf-8",
              newline="\n") as f:
        json.dump(summary, f, indent=1, ensure_ascii=False)
        f.write("\n")

    print("=" * 70)
    print("DOCKING SCORE AGREEMENT BETWEEN INDEPENDENT RUNS")
    print("=" * 70)
    for tag, m in meta.items():
        print(f"  {tag:20} {m['rows']:>4} rows  {m['file']}")
    print()
    for p in pairs:
        print(f"  {p['a']:20} vs {p['b']:20} n={p['n_common_compounds']:>3}  "
              f"differing={p['n_values_differing_exactly']:>3}  "
              f"max|d|={p['max_absolute_difference']:.6f}")
    print()
    print(f"  compounds common to all three: {n_common}")
    print(f"  ALL PAIRS BITWISE IDENTICAL  : {all_exact}")
    print()
    print(f"written -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

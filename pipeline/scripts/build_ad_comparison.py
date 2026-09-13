"""build_ad_comparison.py -- producer for comparison_vs_snplib.json.

This file had no producer. It was read by figure_ad_multispace.py, hashed into
that figure's provenance record, copied into Supplementary S10 and published in
the public deposit -- and nothing in the repository wrote it. It had been
assembled by hand, which is exactly the class of artifact this project's
contract layer exists to prevent: a number reaching a released figure without a
traceable path from its inputs. Found by the internal review of 2026-08-25
(ADR-AIIM-136).

What it holds is a side-by-side of two applicability-domain sweeps run under
the identical frozen rule:

  * SNP-Lib, 236 Thai orchid natural products, the library the study started
    from;
  * the LOTUS biosynthetic expansion, narrowed by 2D-pharmacophore and 3D-shape
    similarity to 52 known A2AR ligands and then by docking score to the 497
    compounds that completed full-protocol redocking.

These are two different libraries, not one library before and after a filter.
The second was anchored to the biosynthetic superclasses that dominate the
first and then deliberately selected toward known A2AR chemistry, which is what
makes the comparison worth drawing: if similarity selection could move a
structural applicability domain, it would have moved it here.

Usage:
    python build_ad_comparison.py [--out DIR]
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from kdss_paths import project_root, resolve, under

import pandas as pd


def _root() -> Path:
    # Delegates to kdss_paths, which finds whichever layout this file is inside
    # -- the development tree or a deposit -- by walking up from its own
    # location. It previously tried two absolute paths and exited if neither
    # existed, which on the author's machine meant a copy of the deposit read
    # and wrote the development tree instead of itself.
    return project_root()


ROOT = _root()
SNPLIB_SWEEP = resolve("20_new_runs/ad_descriptor_sweep_2026-08-15")
EXPANDED_SWEEP = resolve("20_new_runs/ad_sweep_selfsim_2026-08-25")

ORDER = ["ECFP4", "FCFP4", "MACCS", "AtomPair", "TopTorsion",
         "Gobbi_Pharm2D", "ErG", "PhysChem"]
STRUCTURAL = [s for s in ORDER if s != "PhysChem"]


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(EXPANDED_SWEEP))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    a = pd.read_csv(under("20_new_runs/ad_descriptor_sweep_2026-08-15", "ad_coverage_by_descriptor.csv")).set_index("descriptor_space")
    b = pd.read_csv(under("20_new_runs/ad_sweep_selfsim_2026-08-25", "ad_coverage_by_descriptor.csv")).set_index("descriptor_space")
    pct = pd.read_csv(under("20_new_runs/ad_sweep_selfsim_2026-08-25", "ad_percentile_position.csv"))

    comparison = []
    for s in ORDER:
        comparison.append(dict(
            descriptor_space=s,
            n_inside_snplib=int(a.loc[s, "n_inside"]),
            n_total_snplib=int(a.loc[s, "n_total"]),
            coverage_fraction_snplib=round(float(a.loc[s, "coverage_fraction"]), 4),
            n_inside_expanded=int(b.loc[s, "n_inside"]),
            n_total_expanded=int(b.loc[s, "n_total"]),
            coverage_fraction_expanded=round(float(b.loc[s, "coverage_fraction"]), 4),
            median_percentile_expanded=round(float(pct[f"{s}_pctile"].median()), 2),
        ))

    med = {s: round(float(pct[f"{s}_pctile"].median()), 2) for s in STRUCTURAL}
    M = pct[[f"{s}_pctile" for s in STRUCTURAL]].to_numpy()
    below_all = int((M < 10).all(axis=1).sum())

    payload = {
        "comparison": comparison,
        "zero_coverage_spaces_expanded": [
            s for s in ORDER if float(b.loc[s, "coverage_fraction"]) == 0.0],
        "structural_median_percentile_range": [min(med.values()), max(med.values())],
        "n_below_10th_percentile_in_all_structural_spaces": below_all,
        "n_expanded": int(len(pct)),
        "key_finding":
            "Two natural-product libraries under one frozen rule. The expanded "
            "library was selected by 2D pharmacophore and 3D shape similarity to "
            "52 known A2AR ligands and then by docking score; it is still 0% "
            "in-domain in ECFP4, FCFP4, topological torsion and Gobbi 2D "
            "pharmacophore space. Deliberate similarity selection toward known "
            "A2AR chemistry did not move the structural applicability domain. "
            "These are two libraries, not one library before and after a filter.",
        "secondary":
            "Read as ranks rather than as admission rates, the seven structural "
            "and pharmacophoric spaces agree: the median expanded-library "
            f"compound sits between percentile {min(med.values())} and "
            f"{max(med.values())} of the training set's own self-similarity, and "
            f"{below_all} of {len(pct)} compounds fall below the 10th percentile "
            "in all seven at once. The spread in admission rates (0% against "
            "9.3% in MACCS) is where each space's threshold happens to fall in a "
            "differently shaped distribution, not a disagreement about the "
            "compounds. The one real division is structural against "
            "physicochemical, and that is what an applicability-domain verdict "
            "must name to be interpretable (ADR-AIIM-128).",
        "provenance": {
            "script": "build_ad_comparison.py",
            "note": "Written 2026-08-25. Before this, no script produced this "
                    "file; it was hand-assembled and then deposited.",
            "inputs": {
                "snplib_sweep": {
                    "path": str(SNPLIB_SWEEP.relative_to(ROOT).as_posix()),
                    "ad_coverage_by_descriptor.csv":
                        sha256(under("20_new_runs/ad_descriptor_sweep_2026-08-15", "ad_coverage_by_descriptor.csv"))},
                "expanded_sweep": {
                    "path": str(EXPANDED_SWEEP.relative_to(ROOT).as_posix()),
                    "ad_coverage_by_descriptor.csv":
                        sha256(under("20_new_runs/ad_sweep_selfsim_2026-08-25", "ad_coverage_by_descriptor.csv")),
                    "ad_percentile_position.csv":
                        sha256(under("20_new_runs/ad_sweep_selfsim_2026-08-25", "ad_percentile_position.csv"))},
            },
        },
    }

    dest = out / "comparison_vs_snplib.json"
    dest.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"structural median percentiles: {med}")
    print(f"below 10th percentile in all seven structural spaces: "
          f"{below_all} of {len(pct)}")
    print(f"written -> {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

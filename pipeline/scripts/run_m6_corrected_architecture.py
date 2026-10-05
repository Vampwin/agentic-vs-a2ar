"""run_m6_corrected_architecture.py -- M6, pre-registered in ADR-AIIM-263.

Would the corrected architecture have selected different compounds? Two parts,
both computable from frozen artifacts. **Nothing is docked, scored or fitted.**

M6.1  The precondition nobody measured
      `control_representativeness` requires similarity coverage >= 0.50 and
      scaffold coverage >= 0.10 between the redocking control and the screening
      library. It was measured at 0.000/0.000 for the original 236-compound
      library and **never recomputed after the library was replaced**
      (ADR-AIIM-211). This computes it, then re-runs the frozen confidence-tier
      evaluator over the same nine configurations with the measured value in
      place of the missing one.

M6.2  Selection under the completed gate
      G6 v1.4.0's `scaffold_uniqueness` HALTs the released panel because two of
      its three compounds share a Bemis-Murcko framework (ADR-AIIM-239). That
      was applied retrospectively. This applies it at **selection** time.

### Predictions, fixed in ADR-AIIM-263 before this file was run

    P1  control_representativeness still fails, so the consensus remains
        unlicensed in all nine configurations -- the omission was real but not
        decisive. FALSIFIED if it passes, which would mean the project refused
        a ranking it was entitled to make.
    P2  the Tier-A set changes: at least one of the two scaffold-identical
        compounds is replaced.
    P3  the replacement comes from the existing eligible pool, so the panel is
        not left short. FALSIFIED if no eligible compound of a distinct
        scaffold exists -- the completed gate would then shrink the panel
        rather than reorder it.

### Populations, fixed here before any number was read

M6.1's **primary** library is the **497 compounds that received the full
protocol and were tiered** -- that is the library the confidence-tier assessment
is about. The **5,000-compound shortlist** is reported beside it, because a
reader may reasonably ask whether the answer is an artefact of the top cut.

### Everything is imported, nothing restated

The coverage rule comes from `evaluators/_control_representativeness.assess`,
the tier rule from `evaluators/_confidence_tier`, the scaffold definition from
the same `_framework` the G6 evaluator uses, and the scorecard's own constants
from `run_selection_scorecard`. No threshold is redefined here.

Output -> 20_new_runs/m6_corrected_architecture_2026-09-24/
"""
from __future__ import annotations

import argparse
import io
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "10_new_system"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from kdss_paths import under  # noqa: E402

from evaluators._control_representativeness import assess, _framework  # noqa: E402
from np_library_prefilter import ZM241385_SMILES  # noqa: E402
import run_selection_scorecard as SC  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

E8 = Path(under("20_new_runs", "clinical_reference_2026-09-03/inputs/"
                "e8_np_library_prioritization.csv"))
SHORTLIST = Path(under("20_new_runs", "np_library_funnel_2026-08-15/"
                       "np_library_docking_shortlist.csv"))
OUT = Path(under("20_new_runs", "m6_corrected_architecture_2026-09-24"))

TIER_A_SET = "{0,1}"          # the conformal set that does not exclude activity
MIN_SIM_COV = SC.assess_sbdd.__globals__.get("DEFAULT_MIN_SIMILARITY_COVERAGE", 0.50)
MIN_SCAF_COV = SC.assess_sbdd.__globals__.get("DEFAULT_MIN_SCAFFOLD_COVERAGE", 0.10)


def coverage(library_smiles) -> dict:
    return assess([ZM241385_SMILES], list(library_smiles))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    e8 = pd.read_csv(E8)
    scored = e8[e8.vina_score.notna()].copy()
    print(f"E8 table: {len(e8)} rows, {len(scored)} with a full-protocol score")

    # ---- M6.1 -----------------------------------------------------------
    primary = coverage(scored.SMILES)
    if not primary.get("available"):
        print(f"REFUSING: coverage not computable -- {primary.get('reason')}")
        return 2
    secondary = None
    if SHORTLIST.is_file():
        sl = pd.read_csv(SHORTLIST)
        col = next((c for c in ("smiles_standardized", "SMILES") if c in sl), None)
        if col:
            secondary = coverage(sl[col])

    passes = (primary["similarity_coverage"] >= MIN_SIM_COV
              and primary["scaffold_coverage"] >= MIN_SCAF_COV)
    ctrl_measured = {"available": True,
                     "similarity_coverage": primary["similarity_coverage"],
                     "scaffold_coverage": primary["scaffold_coverage"]}

    # Re-run the frozen scorecard, changing only the one input that was absent.
    rows_before = [SC.score(r["case"], r["ad_coverage_used"],
                            SC.CTRL_236 if r["case"].startswith("A") else SC.CTRL_EXPANSION,
                            r["note"])
                   for r in json.load(io.open(
                       under("20_new_runs", "selection_scorecard_2026-09-22/"
                             "selection_scorecard.json"), encoding="utf-8"))["rows"]]
    rows_after = []
    for r in rows_before:
        ctrl = SC.CTRL_236 if r["case"].startswith("A") else ctrl_measured
        rows_after.append(SC.score(r["case"], r["ad_coverage_used"], ctrl, r["note"]))

    lic_before = sum(1 for r in rows_before if r["consensus_licensed"])
    lic_after = sum(1 for r in rows_after if r["consensus_licensed"])
    p1_met = bool(lic_after == 0)

    # ---- M6.2 -----------------------------------------------------------
    eligible = scored[(scored.on_pocket.astype(str).str.lower().isin(["true", "1", "yes"]))
                      & (scored["prediction_set_sig0.05"].astype(str).str.strip()
                         == TIER_A_SET)].copy()
    eligible["framework"] = [_framework_safe(s) for s in eligible.SMILES]
    eligible = eligible.sort_values("vina_score")
    released = sorted(scored[scored.tier == "A_converging_sbdd_lbdd"].compound_id)

    kept, seen = [], set()
    for _, row in eligible.iterrows():
        f = row["framework"]
        if f in seen:
            continue
        seen.add(f)
        kept.append(row)
    new_panel = sorted(r["compound_id"] for r in kept)

    p2_met = bool(set(new_panel) != set(released))
    p3_met = bool(len(new_panel) >= len(released))

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "executes": "M6, pre-registered in ADR-AIIM-263",
        "nothing_docked_scored_or_fitted": True,
        "m6_1_control_representativeness": {
            "control": "ZM241385",
            "thresholds": {"min_similarity_coverage": MIN_SIM_COV,
                           "min_scaffold_coverage": MIN_SCAF_COV},
            "primary_library": {
                "what": "the 497 compounds that received the full protocol",
                "n": int(len(scored)),
                "similarity_coverage": primary["similarity_coverage"],
                "scaffold_coverage": primary["scaffold_coverage"],
                "passes": bool(passes)},
            "secondary_library": (
                {"what": "the 5,000-compound docking shortlist",
                 "similarity_coverage": secondary["similarity_coverage"],
                 "scaffold_coverage": secondary["scaffold_coverage"]}
                if secondary and secondary.get("available") else None),
            "consensus_licensed_before": lic_before,
            "consensus_licensed_after": lic_after,
            "n_configurations": len(rows_after),
        },
        "m6_2_selection_under_completed_gate": {
            "rule": "tier-A criteria imported unchanged, then scaffold "
                    "uniqueness applied as a selection constraint",
            "n_eligible": int(len(eligible)),
            "n_distinct_frameworks_among_eligible": int(eligible.framework.nunique()),
            "released_panel": released,
            "panel_under_completed_gate": new_panel,
            "dropped": sorted(set(released) - set(new_panel)),
            "added": sorted(set(new_panel) - set(released)),
        },
        "predictions": {
            "P1": {"statement": "control_representativeness still fails, so the "
                                "consensus remains unlicensed in all nine",
                   "met": p1_met,
                   "observed": {"passes": bool(passes),
                                "consensus_licensed_after": lic_after}},
            "P2": {"statement": "the Tier-A set changes", "met": p2_met},
            "P3": {"statement": "the replacement comes from the existing "
                                "eligible pool, so the panel is not left short",
                   "met": p3_met,
                   "observed": {"released": len(released),
                                "under_completed_gate": len(new_panel)}},
        },
        "not_claimed": [
            "that any compound here is better than any other; no compound in "
            "this project has been assayed",
            "that this replaces the released selection, which stays in the "
            "record together with the fact that G6 v1.4.0 halts it",
        ],
    }
    (out / "m6_corrected_architecture.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")
    pd.DataFrame(kept).to_csv(out / "panel_under_completed_gate.csv", index=False)
    eligible.to_csv(out / "eligible_with_frameworks.csv", index=False)

    m1 = report["m6_1_control_representativeness"]
    print("=" * 74)
    print("M6.1  THE PRECONDITION NOBODY MEASURED")
    print("=" * 74)
    print(f"  control ZM241385 vs the {m1['primary_library']['n']} full-protocol compounds")
    print(f"    similarity coverage {m1['primary_library']['similarity_coverage']:.4f} "
          f"(needs >= {MIN_SIM_COV})")
    print(f"    scaffold coverage   {m1['primary_library']['scaffold_coverage']:.4f} "
          f"(needs >= {MIN_SCAF_COV})")
    print(f"    passes: {passes}")
    if m1["secondary_library"]:
        s = m1["secondary_library"]
        print(f"  the 5,000 shortlist: similarity {s['similarity_coverage']:.4f}, "
              f"scaffold {s['scaffold_coverage']:.4f}")
    print(f"\n  consensus licensed: {lic_before} of {len(rows_after)} before, "
          f"{lic_after} after supplying the measured value")
    print(f"  P1: {'MET' if p1_met else 'MISSED'}")

    m2 = report["m6_2_selection_under_completed_gate"]
    print("\n" + "=" * 74)
    print("M6.2  SELECTION UNDER THE COMPLETED GATE")
    print("=" * 74)
    print(f"  eligible for tier A: {m2['n_eligible']}, distinct frameworks "
          f"{m2['n_distinct_frameworks_among_eligible']}")
    print(f"  released          : {', '.join(m2['released_panel'])}")
    print(f"  completed gate    : {', '.join(m2['panel_under_completed_gate'])}")
    print(f"  dropped {m2['dropped'] or 'none'}   added {m2['added'] or 'none'}")
    print(f"  P2: {'MET' if p2_met else 'MISSED'}   P3: {'MET' if p3_met else 'MISSED'}")
    print(f"\nwritten -> {out}")
    return 0


def _framework_safe(smiles: str) -> str:
    try:
        from rdkit import Chem
        m = Chem.MolFromSmiles(smiles)
        return _framework(m) or "<unparsed>"
    except Exception:
        return "<unparsed>"


if __name__ == "__main__":
    raise SystemExit(main())

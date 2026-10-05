"""run_g6_retrospective.py -- what the completed G6 contract says about the
selection this project actually released.

ADR-AIIM-237 found that `eligible_set_completeness` and `scaffold_uniqueness`
are implemented in `evaluators/G6_evaluator.py` and named by no deployed
contract. ADR-AIIM-238 added them to a copy of the contract and measured, on the
seeded benchmark, that they recover 27 of 128 faults.

This asks the question the benchmark cannot: **what do they say about the real
Tier-A panel?**

It re-evaluates the published selection under both contract versions, with the
same inputs, by importing `run_g6_on_np_library` and swapping only the contract
file. Nothing is re-selected, re-ranked or re-docked; the selection is read from
`g6_selection_2026-09-02` exactly as it was released.

    v1.3.0   9 checks, the frozen contract, what was published
    v1.4.0  11 checks, the same plus the two the evaluator already implemented

The scaffold table is computed independently of the evaluator -- straight from
the published `selected_tierA.csv` with RDKit -- so the finding does not rest on
the check whose addition is under discussion.

Output -> 20_new_runs/g6_retrospective_<date>/
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
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(ROOT / "10_new_system"))

from kdss_paths import under  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

FROZEN_CONTRACT = ROOT / "10_new_system/contracts/G6_selection.yaml"
AMENDED_CONTRACT = ROOT / "12_patched_pipeline_2026-09-23/contracts/G6_selection.yaml"
PUBLISHED_PANEL = Path(under("20_new_runs", "g6_selection_2026-09-02",
                             "selected_tierA.csv"))
PUBLISHED_OUTCOME = Path(under("20_new_runs", "g6_selection_2026-09-02",
                               "G6_contract_outcome_tierA.json"))


def scaffold_table(panel_csv: Path, funnel_dir: Path) -> dict:
    """Murcko scaffold per published Tier-A member, computed from source.

    Structures come from the published panel file itself, which carries a SMILES
    column, so the table rests on the released artifact and not on a
    re-derivation. The funnel-ledger join is the fallback for a panel file
    written without one.
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem.Scaffolds import MurckoScaffold
    RDLogger.DisableLog("rdApp.*")

    panel = pd.read_csv(panel_csv)
    if "SMILES" in panel.columns:
        m, source = panel, "the published panel file itself"
    else:
        ledger = pd.read_csv(funnel_dir / "np_library_dock_ledger.csv",
                             usecols=["compound_id", "SMILES"])
        m = panel.merge(ledger, on="compound_id", how="left")
        source = "the funnel ledger the selection ranked"
    if m.SMILES.isna().any():
        raise SystemExit(f"ABORT: {int(m.SMILES.isna().sum())} panel members have "
                         f"no SMILES.")
    print(f"structures read from {source}")

    rows = []
    for cid, smi in zip(m.compound_id, m.SMILES):
        mol = Chem.MolFromSmiles(str(smi))
        scaf = Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol)) if mol else None
        rows.append({"compound_id": cid, "smiles": smi, "murcko_scaffold": scaf})
    by_scaf: dict[str, list[str]] = {}
    for r in rows:
        by_scaf.setdefault(r["murcko_scaffold"], []).append(r["compound_id"])
    dupes = {s: ids for s, ids in by_scaf.items() if len(ids) > 1}
    return {"panel": rows, "n_panel": len(rows), "n_distinct_scaffolds": len(by_scaf),
            "duplicate_scaffolds": dupes, "structures_from": source}


def evaluate_under(contract: Path, run_dir: Path) -> dict:
    """Run the released selection through G6 under one contract definition."""
    import run_g6_on_np_library as R
    saved_contract, saved_argv = R.CONTRACT, sys.argv
    try:
        R.CONTRACT = contract.resolve()
        sys.argv = ["run_g6_on_np_library.py", "--run-dir", str(run_dir)]
        R.main()
    finally:
        R.CONTRACT, sys.argv = saved_contract, saved_argv
    return json.loads((run_dir / "G6_contract_outcome_tierA.json")
                      .read_text(encoding="utf-8"))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(under("20_new_runs",
                                               "g6_retrospective_2026-09-23")))
    ap.add_argument("--funnel-dir", default=str(under(
        "20_new_runs", "np_library_funnel_2026-08-15")))
    args = ap.parse_args()
    out = Path(args.out)
    if out.exists() and any(p for p in out.iterdir()):
        print(f"refusing: {out.name} is not empty (ADR-AIIM-092)")
        return 2
    out.mkdir(parents=True, exist_ok=True)

    scaf = scaffold_table(PUBLISHED_PANEL, Path(args.funnel_dir))
    print(f"published Tier-A panel: {scaf['n_panel']} compounds, "
          f"{scaf['n_distinct_scaffolds']} distinct Murcko scaffolds")
    for s, ids in scaf["duplicate_scaffolds"].items():
        print(f"  shared scaffold {s}\n    {ids}")

    frozen = evaluate_under(FROZEN_CONTRACT, out / "under_v1.3.0")
    amended = evaluate_under(AMENDED_CONTRACT, out / "under_v1.4.0")
    published = json.loads(PUBLISHED_OUTCOME.read_text(encoding="utf-8"))

    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "adr": ["ADR-AIIM-237", "ADR-AIIM-238", "ADR-AIIM-239"],
        "what_this_is": "The released Tier-A selection re-evaluated under both "
                        "contract versions. Nothing re-selected, re-ranked or docked.",
        "published_state": published.get("state"),
        "under_frozen_v1_3_0": {"version": frozen.get("contract_version"),
                                "state": frozen.get("state"),
                                "failed": frozen.get("failed"),
                                "warned": frozen.get("warned")},
        "under_amended_v1_4_0": {"version": amended.get("contract_version"),
                                 "state": amended.get("state"),
                                 "failed": amended.get("failed"),
                                 "warned": amended.get("warned")},
        "reproduces_published_state": frozen.get("state") == published.get("state"),
        "scaffolds": scaf,
        "not_claimed": [
            "that the panel members are inactive, or active; no compound here has "
            "been assayed",
            "that the selection rule was applied incorrectly -- it was applied as "
            "frozen, and the gate that would have questioned the result did not "
            "exist in any contract",
        ],
    }
    (out / "g6_retrospective.json").write_text(json.dumps(report, indent=2),
                                               encoding="utf-8")
    print(f"\npublished state   : {published.get('state')}")
    print(f"under v1.3.0      : {frozen.get('state')}  "
          f"(reproduces published: {report['reproduces_published_state']})")
    print(f"under v1.4.0      : {amended.get('state')}  failed={amended.get('failed')}")
    print(f"\nwrote {out / 'g6_retrospective.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

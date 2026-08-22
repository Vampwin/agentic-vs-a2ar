"""
run_g6_on_np_library.py — execute the G6 selection contract over the E8
prioritisation.

G6 carries `applicability_domain_multi_space`, the corrected gate this project
introduces. Until now it had never been applied to the result the study
actually releases (see 23_VS_PIPELINE_STEPS.md). This runs it.

Two of G6's checks encode the *SNP-Lib* selection procedure — a Pareto/diversity
rule with Butina clustering, specified in the assay-blind protocol. The E8
tiering is a different rule (on-pocket + conformal convergence, plus a
descriptive score quartile) and performs no diversity clustering. Those inputs
are therefore genuinely absent rather than withheld, and the contract is
expected to say so. That outcome is reported, not engineered around: a gate
that passed by being fed substitutes would be worthless as evidence.

Usage:
    python run_g6_on_np_library.py
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd
import yaml


def _root() -> Path:
    marker = "00_project_control/08_DECISION_LOG.md"
    for cand in (Path(r"D:\Agentic VS Workflow"), Path("/mnt/d/Agentic VS Workflow")):
        if (cand / marker).is_file():
            return cand
    raise SystemExit("project root not found")


ROOT = _root()
sys.path.insert(0, str(ROOT / "10_new_system"))

from evaluators import G6_evaluator  # noqa: E402
from kdss_audit.artifact_registry import (  # noqa: E402
    ArtifactEntry, ArtifactRegistry, ContractOutcomeEntry, FileAccessEntry, sha256_file,
)

FUNNEL = ROOT / "20_new_runs/np_library_funnel_2026-08-15"
ADSWEEP = ROOT / "20_new_runs/ad_sweep_expanded_2026-08-21"
RUN = ROOT / "20_new_runs/g6_selection_2026-08-21"
CONTRACT = ROOT / "10_new_system/contracts/G6_selection.yaml"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tier", default="A", choices=["A", "AB"],
                    help="which tier(s) constitute the released selection")
    args = ap.parse_args()
    RUN.mkdir(parents=True, exist_ok=True)
    registry = ArtifactRegistry(RUN / "artifact_registry.jsonl")

    contract = yaml.safe_load(CONTRACT.read_text(encoding="utf-8"))
    registry.record_file_access(FileAccessEntry(
        reader="run_g6_on_np_library", file_path=str(CONTRACT),
        purpose="G6 contract definition"))

    ranking = FUNNEL / "e8_np_library_prioritization.csv"
    df = pd.read_csv(ranking)
    sel_tiers = (["A_converging_sbdd_lbdd"] if args.tier == "A"
                 else ["A_converging_sbdd_lbdd", "B_strong_sbdd_lbdd_null"])
    selected = df[df.tier.isin(sel_tiers)].copy()
    sel_path = RUN / f"selected_tier{args.tier}.csv"
    selected.to_csv(sel_path, index=False)

    registry.register_artifact(ArtifactEntry(
        artifact_id="G6::np_library_prioritisation",
        path=str(ranking), sha256=sha256_file(ranking), stage="G6",
        description="E8 confidence-tiered prioritisation of the top-500 NP-library "
                    "shortlist; the result this study releases.",
        tool="e8_prioritization.py", tool_version="2026-08-19",
        conda_env="a2ar", command="e8_prioritization.py", random_seed=None,
        extra={"n_rows": int(len(df)), "tiers_selected": sel_tiers,
               "n_selected": int(len(selected))}))

    # Multi-space applicability, computed on the redocked compounds themselves
    # (ad_sweep_expanded_2026-08-21), not inherited from the 236-compound run.
    cov = pd.read_csv(ADSWEEP / "ad_coverage_by_descriptor.csv")
    # _domain_gate_bias_probe expects one dict per space, not a bare float --
    # it is written against the row shape ad_descriptor_sweep.py emits.
    by_space = {
        r.descriptor_space: {
            "coverage_fraction": float(r.coverage_fraction),
            "n_inside": int(r.n_inside), "n_total": int(r.n_total),
            "threshold": float(r.threshold),
        }
        for r in cov.itertuples()
    }

    inputs: dict[str, object] = {
        "final_ranking_csv": str(ranking),
        "selected_candidates_csv": str(sel_path),
        "column_scan_log": [{"file": ranking.name, "columns": list(df.columns)},
                            {"file": sel_path.name, "columns": list(selected.columns)}],
        "applicability_domain_annotation_present": True,
        "applicability_domain_coverage_fraction": by_space["ECFP4"]["coverage_fraction"],
        "applicability_domain_by_space": by_space,
        "selection_count_report": {
            "total_selected": int(len(selected)),
            "by_tier": selected.tier.value_counts().to_dict(),
            "evaluated_on_pocket": int(((df.status == "ok") & (df.on_pocket == True)).sum()),  # noqa: E712
        },
        "selection_script_hash": sha256_file(ROOT / "10_new_system/scripts/e8_prioritization.py"),
        # Deliberately NOT supplied, because they do not exist for this selection:
        #   registered_selection_script_hash  - e8_prioritization.py was not
        #       hash-registered before it ran
        #   diversity_seed                    - no Butina clustering was performed
    }

    result = G6_evaluator.evaluate(contract, inputs, registry)

    registry.record_contract_outcome(ContractOutcomeEntry(
        contract_id=contract.get("contract_id", "G6_selection_v1"), gate="G6",
        artifact_ids=["G6::np_library_prioritisation"], outcome=result["state"],
        checks_passed=result["passed"], checks_failed=result["failed"],
        checks_warned=result["warned"],
        state_history=[{"state": result["state"], "reason": result["summary"]}]))

    out = RUN / f"G6_contract_outcome_tier{args.tier}.json"
    out.write_text(json.dumps({
        "contract": contract.get("contract_id"),
        "contract_version": contract.get("version"),
        "ranking_sha256": sha256_file(ranking),
        "ad_source": str(ADSWEEP / "ad_coverage_by_descriptor.csv"),
        "ad_by_space": by_space,
        "inputs_deliberately_absent": {
            "registered_selection_script_hash":
                "e8_prioritization.py was not hash-registered before execution",
            "diversity_seed":
                "the E8 tiering performs no Butina diversity clustering",
        },
        **result,
    }, indent=2, default=str), encoding="utf-8")

    print(f"selected (tier {args.tier}) : {len(selected)}")
    print(f"G6 state : {result['state']}")
    print(f"summary  : {result['summary']}\n")
    for cid in result["passed"]:
        print(f"  PASS  {cid}")
    for cid in result["warned"]:
        print(f"  WARN  {cid}: {json.dumps(result['evidence'].get(cid, {}), default=str)[:150]}")
    for cid in result["failed"]:
        print(f"  FAIL  {cid}: {json.dumps(result['evidence'].get(cid, {}), default=str)[:200]}")
    print(f"\nwritten -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""
run_g4_on_np_library.py — execute the G4 scoring contract over the
natural-product library full-protocol redock.

Why this script exists: the NP-library arm produced the manuscript's headline
result but was run by invoking `g4_dock_library.py` directly, outside the gate
architecture. No G4 contract outcome and no artifact registry exist for it (see
`00_project_control/23_VS_PIPELINE_STEPS.md`). This script applies the existing,
already-frozen G4 contract and evaluator to the artifacts that arm produced, so
the governance claim covers the result actually reported.

It does not re-dock the library. The only new computation is the
`seed_repeatability` rerun the contract itself mandates (5% seeded subset),
which is performed separately by `run_g4_verification.sh` and read in here.

Usage:
    python run_g4_on_np_library.py [--rerun-csv <path>] [--label <tag>]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(r"D:\Agentic VS Workflow")
sys.path.insert(0, str(ROOT / "10_new_system"))

from evaluators import G4_evaluator  # noqa: E402
from kdss_audit.artifact_registry import (  # noqa: E402
    ArtifactEntry,
    ArtifactRegistry,
    ContractOutcomeEntry,
    FileAccessEntry,
    sha256_file,
)

FUNNEL = ROOT / "20_new_runs/np_library_funnel_2026-08-15"
RUN = ROOT / "20_new_runs/g4_verification_2026-08-21"
CONTRACT = ROOT / "10_new_system/contracts/G4_scoring.yaml"


def _read(registry: ArtifactRegistry, path: Path, purpose: str) -> Path:
    registry.record_file_access(
        FileAccessEntry(reader="run_g4_on_np_library", file_path=str(path), purpose=purpose)
    )
    return path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rerun-csv", default=str(RUN / "rerun_passA_cpu6.csv"),
                    help="seed_repeatability rerun output; if absent the check "
                         "fails closed, which is the contract's intent")
    ap.add_argument("--label", default="passA_cpu6")
    args = ap.parse_args()

    RUN.mkdir(parents=True, exist_ok=True)
    registry = ArtifactRegistry(RUN / "artifact_registry.jsonl")

    contract = yaml.safe_load(CONTRACT.read_text(encoding="utf-8"))
    _read(registry, CONTRACT, "G4 contract definition")

    merged = _read(registry, RUN / "full_scores_MERGED.csv",
                   "G4 input: merged full-protocol docking scores")
    if not merged.is_file():
        print(f"ERROR: {merged} not found", file=sys.stderr)
        return 2

    scores = pd.read_csv(merged)
    registry.register_artifact(ArtifactEntry(
        artifact_id="G4::np_library_full_scores",
        path=str(merged), sha256=sha256_file(merged), stage="G4",
        description="Full-protocol Vina redock of the top-500 NP-library shortlist "
                    "(497 ok, 3 documented manual skips), merged from 3 parallel "
                    "chunks and left-joined with molecular weight and SMILES from "
                    "the funnel ledger so mw_bias_check has an MW column to work "
                    "with. Scores themselves are unmodified.",
        tool="g4_dock_library.py", tool_version="2026-08-18",
        conda_env="diffdock",
        command="g4_dock_library.py --n-conformers 8 --exhaustiveness 32 --n-poses 10",
        random_seed=20260729,
        extra={"n_rows": int(len(scores)),
               "n_ok": int((scores.status == "ok").sum()),
               "n_failed": int((scores.status == "failed").sum())},
    ))

    # Expected compound set and the documented-failure ledger, so
    # missing_output_accounting can distinguish "scored" from "explained".
    ledger = pd.concat(
        [pd.read_csv(_read(registry, FUNNEL / f"full_ledger_chunk{c}.csv",
                           "G4 input: expected compound list")) for c in (1, 2, 3)],
        ignore_index=True,
    )
    failure_ledger = {
        str(r.compound_id): str(r.failure_reason)
        for r in scores.itertuples() if str(r.status) == "failed"
    }

    inputs: dict[str, object] = {
        "vina_score_csv": str(merged),
        "count_for_docking": int(len(ledger)),
        "expected_compound_ids": ledger["compound_id"].astype(str).tolist(),
        "failure_ledger": failure_ledger,
    }

    rerun_path = Path(args.rerun_csv)
    if rerun_path.is_file():
        _read(registry, rerun_path, "G4 input: seed_repeatability rerun")
        rr = pd.read_csv(rerun_path)
        rr_ok = rr[rr.status == "ok"]
        inputs["seed_rerun_scores"] = dict(
            zip(rr_ok["compound_id"].astype(str), rr_ok["vina_score"].astype(float))
        )
        print(f"seed rerun loaded: {len(rr_ok)} compounds from {rerun_path.name}")
    else:
        print(f"NOTE: rerun not found at {rerun_path}; seed_repeatability will "
              f"fail closed, which is the contract's designed behaviour.")

    result = G4_evaluator.evaluate(contract, inputs, registry)

    registry.record_contract_outcome(ContractOutcomeEntry(
        contract_id=contract.get("contract_id", "G4_scoring_v1"),
        gate="G4",
        artifact_ids=["G4::np_library_full_scores"],
        outcome=result["state"],
        checks_passed=result["passed"],
        checks_failed=result["failed"],
        checks_warned=result["warned"],
        state_history=[{"state": result["state"], "reason": result["summary"]}],
    ))

    out = RUN / f"G4_contract_outcome_{args.label}.json"
    out.write_text(json.dumps({
        "contract": contract.get("contract_id"),
        "contract_version": contract.get("version"),
        "target_artifact": str(merged),
        "target_sha256": sha256_file(merged),
        "rerun_source": str(rerun_path) if rerun_path.is_file() else None,
        **result,
    }, indent=2, default=str), encoding="utf-8")

    print(f"\nG4 state : {result['state']}")
    print(f"summary  : {result['summary']}")
    for cid in result["passed"]:
        print(f"  PASS  {cid}")
    for cid in result["warned"]:
        print(f"  WARN  {cid}: {json.dumps(result['evidence'].get(cid, {}), default=str)[:160]}")
    for cid in result["failed"]:
        print(f"  FAIL  {cid}: {json.dumps(result['evidence'].get(cid, {}), default=str)[:220]}")
    print(f"\nwritten -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

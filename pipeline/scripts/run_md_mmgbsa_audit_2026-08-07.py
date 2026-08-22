"""
run_md_mmgbsa_audit_2026-08-07.py — Phase B: audit the existing MD/MM-GBSA
results under G7, do not recompute them.

Human decision (this session): rerunning the 50 ns MD trajectories is not
worth the wall-clock/GPU cost, and was already explicitly ruled out in this
project's own earlier plan ("MD 50 ns: reuse 15 systems, do not rerun").
What this script does instead is register the real, existing legacy MD/
MM-GBSA output files as artifacts with real hashes, and run the newly
implemented G7_evaluator (previously a stub) against them for the first
time -- auditing provenance completeness, not recomputing science.

Source files (real, unmodified, from the legacy project):
  D:\\kdss_SNP v2\\...\\data\\output\\md\\12b_md_final_panel_summary.csv
    10 real 50 ns trajectories, all status=completed
  D:\\kdss_SNP v2\\...\\data\\output\\mmgbsa\\15_mmgbsa_scores.csv
    172 real MM-GBSA interaction-energy rescores

NOT used: 16_final_mmgbsa_ranked.csv, which carries a `known_active_cAMP`
column -- forbidden under ADR-AIIM-004 and this project's assay-blindness
rules throughout. G7_evaluator itself also refuses to read any file with a
blocked column (test_refuses_to_read_a_file_with_assay_columns), so this is
enforced twice, not just by this script's own care.
"""

from __future__ import annotations

import sys
from pathlib import Path

_SYSTEM_ROOT = Path(__file__).resolve().parent.parent
if str(_SYSTEM_ROOT) not in sys.path:
    sys.path.insert(0, str(_SYSTEM_ROOT))

import json
import yaml

from kdss_audit.artifact_registry import ArtifactEntry, ArtifactRegistry, FileAccessEntry, sha256_file
from evaluators import G7_evaluator

ROOT = _SYSTEM_ROOT.parent
LEGACY_MD = Path(r"D:/kdss_SNP v2/KDSS_A2AR_Code_Data_Package_2026-08-03/data/output/md/12b_md_final_panel_summary.csv")
LEGACY_MMGBSA = Path(r"D:/kdss_SNP v2/KDSS_A2AR_Code_Data_Package_2026-08-03/data/output/mmgbsa/15_mmgbsa_scores.csv")
OUT_DIR = ROOT / "20_new_runs" / "md_mmgbsa_audit_2026-08-07"
APPROVED_BY = "adisorn.siric@gmail.com"

CLAIM_BOUNDARY = (
    "MD and MM-GBSA results are binding-stability and interaction-energy "
    "annotations only. They do not constitute experimental validation, "
    "were not used to alter or re-rank any compound selection, and (per "
    "prior external-benchmark evaluation, MM-GBSA AUC=0.533 against ChEMBL) "
    "are not used to drive ranking."
)


def log(msg: str) -> None:
    print(f"[md_audit] {msg}")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    registry = ArtifactRegistry(OUT_DIR / "artifact_registry.jsonl")

    def record_read(path: Path, purpose: str) -> None:
        registry.record_file_access(FileAccessEntry(reader="run_md_mmgbsa_audit",
                                                      file_path=str(path), purpose=purpose))

    def register(artifact_id: str, path: Path, description: str, decision: str = "") -> Path:
        local = OUT_DIR / path.name
        local.write_bytes(path.read_bytes())
        record_read(path, description)
        registry.register_artifact(ArtifactEntry(
            artifact_id=artifact_id, path=str(local), sha256=sha256_file(local), stage="G7_MD",
            description=description, tool="legacy MD/MM-GBSA pipeline (unchanged, not rerun)",
            tool_version="", conda_env="", command="", random_seed=None,
            agent_decision_summary=decision,
        ))
        return local

    if not LEGACY_MD.is_file() or not LEGACY_MMGBSA.is_file():
        raise FileNotFoundError("Legacy MD/MM-GBSA source files not found")

    md_local = register("md_trajectory_summary", LEGACY_MD,
                         "10 real 50ns MD trajectory summaries (legacy, not rerun)",
                         "Phase B of the QSAR/MD/MM-GBSA rerun decision (2026-08-07): "
                         "audit existing results under G7, do not recompute.")
    register("mmgbsa_scores", LEGACY_MMGBSA,
              "172 real MM-GBSA interaction-energy rescores (legacy, not rerun). "
              "16_final_mmgbsa_ranked.csv NOT used -- contains a known_active_cAMP "
              "column forbidden under ADR-AIIM-004.")

    contract = yaml.safe_load((_SYSTEM_ROOT / "contracts" / "G7_md.yaml").read_text(encoding="utf-8"))

    inputs = {
        "md_decision_recorded": "MD_run",
        "md_completion_log": str(md_local),
        "planned_production_ns": 50.0,
        # Genuinely absent from the legacy run -- not supplied, so the
        # equilibration and seed checks fail closed rather than pass on
        # summary statistics that are not the evidence the check asks for.
        # "equilibration_check_csv": not available,
        # "md_seeds_by_compound": not recorded in the legacy logs,
        "required_claim_boundary": CLAIM_BOUNDARY,
        # No SNP-Lib selection/ranking exists (ADR-AIIM-037, negative
        # result) for MD to have re-ranked.
    }

    result = G7_evaluator.evaluate(contract, inputs, registry=registry)
    log(f"G7 state: {result['state']}")
    log(f"passed: {result['passed']}")
    log(f"failed: {result['failed']}")
    log(f"warned: {result['warned']}")
    log(f"critical_failures: {result['critical_failures']}")
    for cid, ev in result["evidence"].items():
        log(f"  {cid}: {ev}")

    (OUT_DIR / "g7_audit_result.json").write_text(json.dumps(result, indent=2, default=str), encoding="utf-8")
    log(f"wrote {OUT_DIR / 'g7_audit_result.json'}")
    log("DONE")


if __name__ == "__main__":
    main()

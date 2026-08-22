"""
G7 evaluator — MD/energy refinement.

Conditional gate: if `md_decision_recorded` is "no_MD", every check is
skipped (SKIPPED_BY_PROTOCOL). If MD was run, this evaluator audits what the
trajectories actually demonstrate — completion, equilibration, replicate
consistency, an explicit annotation-only claim boundary, that MD did not
re-rank a selection, and that a seed was recorded. It does not run or
re-run any simulation.

A skipped or unverifiable critical check must never report PASS: if the
evidence needed for a check (e.g. per-frame RMSD-vs-time for the
equilibration check) was not captured by the run being audited, that is a
real provenance gap, not a pass.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from ._assay_blindness import blocked_columns_in


def evaluate(contract: dict[str, Any], inputs: dict[str, Any], registry: Any) -> dict[str, Any]:
    md_decision = inputs.get("md_decision_recorded", "no_MD")
    if md_decision == "no_MD":
        return {
            "state": "SKIPPED_BY_PROTOCOL",
            "passed": [], "failed": [], "warned": [], "critical_failures": [],
            "evidence": {"md_decision": "no_MD"},
            "summary": "MD not run; G7 skipped by protocol",
        }

    checks = contract.get("checks", [])
    declared = {c["id"] for c in checks}
    passed, failed, warned, evidence = [], [], [], {}
    critical_failures: list[str] = []

    def fail(cid: str, sev: str, detail: dict[str, Any]) -> None:
        failed.append(cid)
        if sev == "critical":
            critical_failures.append(cid)
        evidence[cid] = detail

    def sev_of(cid: str) -> str:
        return next((c.get("severity", "medium") for c in checks if c["id"] == cid), "medium")

    def read_csv(key: str) -> tuple[pd.DataFrame | None, str | None]:
        path_str = inputs.get(key)
        if not path_str:
            return None, f"{key} not supplied"
        path = Path(str(path_str))
        if not path.is_file():
            return None, f"file not found: {path}"
        try:
            df = pd.read_csv(path)
        except Exception as exc:
            return None, f"could not parse {path.name}: {exc}"
        blocked = blocked_columns_in(df.columns, strict=True)
        if blocked:
            return None, f"refused to read {path.name}: forbidden columns present ({blocked})"
        return df, None

    traj_df, traj_err = read_csv("md_completion_log")

    if "trajectory_completion" in declared:
        cid = "trajectory_completion"
        if traj_df is None:
            fail(cid, sev_of(cid), {"error": traj_err})
        else:
            planned_ns = float(inputs.get("planned_production_ns", 50.0))
            status_col = next((c for c in ("status", "completion_status") if c in traj_df.columns), None)
            ns_col = next((c for c in ("prod_ns", "sim_ns", "production_ns") if c in traj_df.columns), None)
            if status_col is None or ns_col is None:
                fail(cid, sev_of(cid), {"error": "status or production-length column not found",
                                        "columns": list(traj_df.columns)})
            else:
                incomplete = traj_df[(traj_df[status_col] != "completed")
                                      | (traj_df[ns_col] < planned_ns)]
                if len(incomplete):
                    fail(cid, sev_of(cid), {"n_incomplete": len(incomplete),
                                            "incomplete": incomplete.to_dict(orient="records")[:20]})
                else:
                    passed.append(cid)
                    evidence[cid] = {"n_trajectories": len(traj_df), "planned_production_ns": planned_ns,
                                     "all_completed": True}

    if "equilibration_verified" in declared:
        cid = "equilibration_verified"
        eq_path = inputs.get("equilibration_check_csv")
        if not eq_path or not Path(str(eq_path)).is_file():
            # A summary RMSD statistic is not the same evidence as a
            # per-frame RMSD-vs-time series showing a production-phase
            # plateau. Reporting PASS from a summary alone would be
            # exactly the "check ran, evidence didn't support it" gap
            # this project's audit layer exists to catch.
            fail(cid, sev_of(cid), {
                "error": "no per-frame RMSD-vs-time data available to verify a production-phase "
                         "plateau; only trajectory-level summary statistics were captured by this run",
            })
        else:
            eq_df, eq_err = read_csv("equilibration_check_csv")
            if eq_df is None:
                fail(cid, sev_of(cid), {"error": eq_err})
            else:
                passed.append(cid)
                evidence[cid] = {"n_trajectories_checked": eq_df["compound_id"].nunique()
                                 if "compound_id" in eq_df.columns else None}

    if "replicate_consistency" in declared:
        cid = "replicate_consistency"
        if traj_df is None:
            fail(cid, sev_of(cid), {"error": traj_err})
        else:
            id_col = next((c for c in ("compound_id",) if c in traj_df.columns), None)
            n_replicates = traj_df[id_col].value_counts() if id_col else pd.Series(dtype=int)
            has_replicates = bool((n_replicates > 1).any())
            if not has_replicates:
                warned.append(cid)
                evidence[cid] = {"warning": "no replicate trajectories were run for any compound; "
                                            "consistency across replicates cannot be assessed",
                                 "n_compounds": len(n_replicates)}
            else:
                passed.append(cid)
                evidence[cid] = {"n_compounds_with_replicates": int((n_replicates > 1).sum())}

    if "claim_boundary" in declared:
        cid = "claim_boundary"
        boundary = inputs.get("required_claim_boundary")
        if not boundary or not str(boundary).strip():
            fail(cid, sev_of(cid), {
                "error": "no required_claim_boundary field recorded for this MD result",
            })
        else:
            passed.append(cid)
            evidence[cid] = {"required_claim_boundary": str(boundary)}

    if "md_not_used_for_ranking" in declared:
        cid = "md_not_used_for_ranking"
        selected_csv = inputs.get("selected_candidates_csv")
        if not selected_csv:
            # No selection exists for this MD result to have re-ranked --
            # vacuously satisfied, not silently assumed.
            passed.append(cid)
            evidence[cid] = {"note": "no selected_candidates_csv exists for this application; "
                                     "nothing for MD output to have re-ranked"}
        else:
            reranked = inputs.get("selection_reranked_by_md", False)
            if reranked:
                fail(cid, sev_of(cid), {"error": "selection was re-ranked using MD output"})
            else:
                passed.append(cid)
                evidence[cid] = {"selected_candidates_csv": str(selected_csv), "reranked": False}

    if "seed_recorded" in declared:
        cid = "seed_recorded"
        seeds = inputs.get("md_seeds_by_compound")
        if traj_df is not None and "compound_id" in traj_df.columns:
            compounds = set(traj_df["compound_id"])
        else:
            compounds = set()
        missing = sorted(compounds - set((seeds or {}).keys())) if compounds else []
        if not seeds or missing:
            fail(cid, sev_of(cid), {"error": "MD random seed not recorded for one or more trajectories",
                                    "n_missing": len(missing) or (compounds and len(compounds)) or None,
                                    "missing_compounds": missing[:20] if missing else None})
        else:
            passed.append(cid)
            evidence[cid] = {"n_seeds_recorded": len(seeds)}

    if critical_failures:
        state = "HALT"
    elif failed:
        state = "ESCALATE"
    else:
        state = "PASS"

    return {
        "state": state, "passed": passed, "failed": failed, "warned": warned,
        "critical_failures": critical_failures, "evidence": evidence,
        "summary": (f"{len(passed)} checks passed, {len(failed)} failed "
                    f"({len(critical_failures)} critical), {len(warned)} warned"),
    }

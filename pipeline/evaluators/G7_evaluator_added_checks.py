"""G7 evaluator with the added checks `trajectory_covers_plan` and `binding_site_matches_receptor` (ADR-AIIM-431): the frozen `G7_evaluator` is
called unchanged and the checks are appended; a run that the frozen evaluator skips by protocol stays skipped."""
from evaluators import G7_evaluator as frozen
from evaluators import added_checks as checks


def evaluate(contract, inputs, registry):
    r = frozen.evaluate(contract=contract, inputs=inputs, registry=registry)
    if r.get("state") == "SKIPPED_BY_PROTOCOL":
        return r
    ok1, d1 = checks.trajectory_covers_plan(inputs.get("equilibration_check_csv"), float(inputs.get("planned_production_ns", 50.0)))
    ok2, d2 = checks.binding_site_matches_receptor(inputs.get("md_results_csv"), inputs.get("receptor_pdbqt"), inputs.get("reference_ligand_sdf"))
    return checks.merge(r, [("trajectory_covers_plan", ok1, "critical", d1), ("binding_site_matches_receptor", ok2, "critical", d2)])

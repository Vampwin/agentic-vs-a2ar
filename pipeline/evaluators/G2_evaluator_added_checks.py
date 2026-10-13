"""G2 evaluator with the added check `protonation_matches_record` (ADR-AIIM-431): the frozen `G2_evaluator` is called unchanged and the check is appended."""
from evaluators import G2_evaluator as frozen
from evaluators import added_checks as checks


def evaluate(contract, inputs, registry):
    r = frozen.evaluate(contract=contract, inputs=inputs, registry=registry)
    ok, det = checks.protonation_matches_record(inputs.get("prepared_receptor_4eiy_pdbqt"), inputs.get("protonation_decision"))
    return checks.merge(r, [("protonation_matches_record", ok, "critical", det)])

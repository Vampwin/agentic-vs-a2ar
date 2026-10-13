import evaluators.G2_evaluator as frozen

from v9ac import checks


def evaluate(contract, inputs, registry):
    r = frozen.evaluate(contract=contract, inputs=inputs, registry=registry)
    ok, det = checks.protonation_matches_record(inputs.get("prepared_receptor_4eiy_pdbqt"), inputs.get("protonation_decision"))
    return checks.merge(r, [("protonation_matches_record", ok, "critical", det)])

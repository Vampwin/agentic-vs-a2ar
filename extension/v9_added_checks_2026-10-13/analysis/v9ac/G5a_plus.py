import evaluators.G5a_evaluator as frozen

from v9ac import checks


def evaluate(contract, inputs, registry):
    r = frozen.evaluate(contract=contract, inputs=inputs, registry=registry)
    ok, det = checks.benchmark_ids_not_exposed(registry, inputs.get("chembl_benchmark_csv"), inputs.get("registered_ranking_inputs", []),
                                               inputs.get("_gate_start_time"))
    return checks.merge(r, [("benchmark_ids_not_exposed", ok, "critical", det)])

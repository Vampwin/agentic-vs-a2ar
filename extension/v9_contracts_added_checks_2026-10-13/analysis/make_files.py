"""ADR-AIIM-431: write the seven NEW files that carry the four added checks into 10_new_system (approved by the author as additive: no existing
file of the tree is opened for writing). Refuses to overwrite. Line ends are CRLF like the rest of the tree."""
import hashlib
import json
import re
from pathlib import Path

ROOT = Path(r"D:\Agentic VS Workflow")
TREE = ROOT / "10_new_system"
SRC = ROOT / "20_new_runs" / "v9_added_checks_2026-10-13" / "analysis" / "v9ac" / "checks.py"
OUT = Path(__file__).resolve().parent.parent
DATE = "2026-10-13"


def put(path: Path, text: str) -> str:
    assert not path.exists(), f"refusing to overwrite {path}"
    data = text.replace("\r\n", "\n").replace("\n", "\r\n").encode("utf-8")
    path.write_bytes(data)
    return hashlib.sha256(data).hexdigest()


written = {}

# ---- the checks module: the functions of the run v9_added_checks_2026-10-13, unchanged below the header
body = SRC.read_text(encoding="utf-8")
body = body.split('"""', 2)[2].lstrip("\n")
header = '''"""Four added checks (ADR-AIIM-427 to 431), kept in a module of their own so that no frozen evaluator is edited. Each returns (ok, detail) and
fails closed: a missing record or an unreadable file is a failure, not a pass. The evaluators `G2_evaluator_added_checks`,
`G7_evaluator_added_checks` and `G5a_evaluator_added_checks` call the frozen evaluator unchanged and append these checks to its output; the
contracts `*_added_checks.yaml` declare them. The frozen contracts and evaluators remain the ones the seeded benchmark of the article measured."""
'''
written["evaluators/added_checks.py"] = put(TREE / "evaluators" / "added_checks.py", header + body)

WRAP = {
    "G2": ("G2_evaluator", "protonation_matches_record", 'checks.protonation_matches_record(inputs.get("prepared_receptor_4eiy_pdbqt"), inputs.get("protonation_decision"))'),
    "G5a": ("G5a_evaluator", "benchmark_ids_not_exposed", 'checks.benchmark_ids_not_exposed(registry, inputs.get("chembl_benchmark_csv"), inputs.get("registered_ranking_inputs", []), inputs.get("_gate_start_time"))'),
}
for g, (mod, cid, call) in WRAP.items():
    text = f'''"""{g} evaluator with the added check `{cid}` (ADR-AIIM-431): the frozen `{mod}` is called unchanged and the check is appended."""
from evaluators import {mod} as frozen
from evaluators import added_checks as checks


def evaluate(contract, inputs, registry):
    r = frozen.evaluate(contract=contract, inputs=inputs, registry=registry)
    ok, det = {call}
    return checks.merge(r, [("{cid}", ok, "critical", det)])
'''
    written[f"evaluators/{mod}_added_checks.py"] = put(TREE / "evaluators" / f"{mod}_added_checks.py", text)

g7 = '''"""G7 evaluator with the added checks `trajectory_covers_plan` and `binding_site_matches_receptor` (ADR-AIIM-431): the frozen `G7_evaluator` is
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
'''
written["evaluators/G7_evaluator_added_checks.py"] = put(TREE / "evaluators" / "G7_evaluator_added_checks.py", g7)

# ---- contracts: the frozen file as text, with the identifier, the evaluator and the version changed and the check entries inserted
NEW_CHECKS = {
    "G2_preparation": [("protonation_matches_record", "critical",
                        "Every histidine of the prepared receptor carries the ring hydrogen (HD1, HE2 or both) that the preparation step recorded in protonation_decision; "
                        "a histidine absent from the record, a record naming a histidine absent from the file, and a missing record all fail.")],
    "G7_md": [("trajectory_covers_plan", "critical",
               "Every per-frame series (compound and replicate) reaches at least 0.95 of planned_production_ns."),
              ("binding_site_matches_receptor", "critical",
               "The residue indices that define the binding site in the results table equal, for every row, the receptor C-alpha residues within 10 A of the reference ligand, recomputed from the prepared receptor and the reference ligand.")],
    "G5a_validity": [("benchmark_ids_not_exposed", "critical",
                      "No file read before the gate by a step other than G5a and G5b, and not listed in registered_ranking_inputs, contains a compound ID of the benchmark (active or inactive) or cannot be opened.")],
}
EXTRA_INPUTS = {
    "G2_preparation": ["protonation_decision"],
    "G7_md": ["receptor_pdbqt", "reference_ligand_sdf"],
    "G5a_validity": ["registered_ranking_inputs"],
}
for stem, evaluator in (("G2_preparation", "G2_evaluator"), ("G7_md", "G7_evaluator"), ("G5a_validity", "G5a_evaluator")):
    t = (TREE / "contracts" / f"{stem}.yaml").read_bytes().decode("utf-8").replace("\r\n", "\n")
    cid = re.match(r"contract_id: (\S+)", t).group(1)
    assert t.count(f"evaluator: evaluators.{evaluator}\n") == 1 and t.count("\nthresholds:\n") == 1
    t = t.replace(f"contract_id: {cid}\n", f"contract_id: {cid}_added_checks\n", 1)
    t = t.replace(f"evaluator: evaluators.{evaluator}\n", f"evaluator: evaluators.{evaluator}_added_checks\n")
    block = "".join(f"  - id: {c}\n    description: >\n      {d}\n    severity: {s}\n" for c, s, d in NEW_CHECKS[stem])
    t = t.replace("\nthresholds:\n", "\n" + block + "\nthresholds:\n")
    m = re.search(r"\nversion: \"([^\"]+)\"\n", t)
    t = t.replace(m.group(0), f'\nversion: "{m.group(1)}+added-checks-{DATE}"\n')
    t = t.rstrip("\n") + (
        f'\nadded_checks:\n  - date: "{DATE}"\n    checks: [{", ".join(c for c, _, _ in NEW_CHECKS[stem])}]\n'
        f'    extra_inputs: [{", ".join(EXTRA_INPUTS[stem])}]\n'
        f'    note: "A separate file: the frozen {stem}.yaml and its evaluator are unchanged. This contract is evaluated by evaluators.{evaluator}_added_checks, which calls the frozen evaluator and appends the checks."\n')
    written[f"contracts/{stem}_added_checks.yaml"] = put(TREE / "contracts" / f"{stem}_added_checks.yaml", t)

(OUT / "written_files.json").write_text(json.dumps(written, indent=1), encoding="utf-8")
print(json.dumps(written, indent=1))

"""A3 and A5 of ADR-AIIM-399.

A3: a stratified 95% interval for the estimate of 1,251 discarded compounds reaching the fast-pass cut-off.
A5: a string-level map from every check id declared in the ten contracts to the source of its evaluator, with the checks the
    seeded benchmark exercised. Reads 10_new_system/ only; writes into this run folder.
"""
import json
import math
import re
from pathlib import Path

import yaml

ROOT = Path(r"D:\Agentic VS Workflow")
OUT = ROOT / "20_new_runs" / "v9_sensitivity_2026-10-11"

# ------------------------------------------------------------------ A3
src = json.loads((ROOT / "20_new_runs/funnel_recall_reweighted_2026-08-31/funnel_recall_reweighted.json").read_text("utf-8"))
tot, var, strata = 0.0, 0.0, {}
for name, s in src["by_stratum"].items():
    N, n, x = s["population"], s["n_docked_ok"], s["n_reaching_cutoff"]
    p = x / n
    tot += N * p
    v = N * N * (1 - n / N) * p * (1 - p) / (n - 1)
    var += v
    strata[name] = {"N": N, "n": n, "hits": x, "rate": round(p, 4), "sampling_fraction": round(n / N, 4)}
Ntot = sum(s["N"] for s in strata.values())
se = math.sqrt(var)
a3 = {"adr": "ADR-AIIM-399", "strata": strata, "total_discarded": Ntot,
      "estimate_count": round(tot, 1), "estimate_rate": round(tot / Ntot, 4),
      "se_count": round(se, 1), "ci95_count": [round(tot - 1.96 * se), round(tot + 1.96 * se)],
      "ci95_rate": [round((tot - 1.96 * se) / Ntot, 4), round((tot + 1.96 * se) / Ntot, 4)],
      "anchor_reproduces_1251_and_7.41pct": round(tot) == src["population_estimate"]["implied_count"]
                                             and round(tot / Ntot, 4) == src["population_estimate"]["weighted_rate"],
      "method": "stratified estimator of a total, finite-population correction, normal approximation"}
(OUT / "a3_stratified_interval.json").write_text(json.dumps(a3, indent=1), encoding="utf-8")
print("A3", json.dumps(a3, indent=1))

# ------------------------------------------------------------------ A5
cov = json.loads((ROOT / "20_new_runs/deployed_path_benchmark_2026-09-23/check_coverage.json").read_text("utf-8"))
exercised = {}
for gate, g in cov["by_gate"].items():
    never = set(g.get("never_exercised", []))
    exercised[gate] = {c: (c not in never) for c in g["deployed_checks"]}
other_src = ""
for d in ("evaluators", "kdss_audit"):
    for f in sorted((ROOT / "10_new_system" / d).rglob("*.py")):
        other_src += f.read_text(encoding="utf-8", errors="ignore") + "\n"
rows, per_gate = [], {}
for f in sorted((ROOT / "10_new_system" / "contracts").glob("*.yaml")):
    c = yaml.safe_load(f.read_text(encoding="utf-8"))
    gate = f.stem.split("_")[0]
    mod = c.get("evaluator", "")
    ev = ROOT / "10_new_system" / Path(*mod.split(".")).with_suffix(".py") if mod else None
    ev_src = ev.read_text(encoding="utf-8", errors="ignore") if ev and ev.exists() else ""
    n = n_own = n_any = 0
    for chk in c.get("checks", []):
        cid = chk["id"]
        lit = re.compile(r"""["']%s["']""" % re.escape(cid))
        own, anyw = bool(lit.search(ev_src)), bool(lit.search(other_src))
        g5 = gate[:2]
        ex = exercised.get(gate, exercised.get(g5, {})).get(cid)
        rows.append({"contract": f.name, "gate": gate, "check": cid, "severity": chk.get("severity"),
                     "evaluator": mod, "in_own_evaluator_source": own, "in_any_evaluator_source": anyw,
                     "exercised_by_seeded_benchmark": ex})
        n += 1; n_own += own; n_any += anyw
    per_gate[f.name] = {"declared": n, "in_own_evaluator": n_own, "in_any_source": n_any,
                        "evaluator_file_exists": bool(ev and ev.exists())}
missing_own = [r for r in rows if not r["in_own_evaluator_source"]]
a5 = {"adr": "ADR-AIIM-399", "tree": "10_new_system as on 2026-10-11 (includes the repairs made after the held-out evaluation)",
      "n_declared": len(rows), "n_in_own_evaluator": sum(r["in_own_evaluator_source"] for r in rows),
      "n_in_any_source": sum(r["in_any_evaluator_source"] for r in rows),
      "not_in_own_evaluator": [{"contract": r["contract"], "check": r["check"], "severity": r["severity"],
                                "in_any_source": r["in_any_evaluator_source"]} for r in missing_own],
      "benchmarked_gates": sorted(exercised),
      "n_declared_in_benchmarked_gates": sum(r["exercised_by_seeded_benchmark"] is not None for r in rows),
      "n_exercised_by_seeded_benchmark": sum(bool(r["exercised_by_seeded_benchmark"]) for r in rows),
      "per_contract": per_gate,
      "P11_declares_94": len(rows) == 94, "P12_at_most_three_absent": len(missing_own) <= 3,
      "limits": "string literals only; an id assembled at run time is missed; presence in the source is not a test of correctness"}
(OUT / "a5_check_map.json").write_text(json.dumps(a5, indent=1), encoding="utf-8")
import csv
with open(OUT / "a5_check_map.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("A5", json.dumps({k: v for k, v in a5.items()}, indent=1))

"""ADR-AIIM-447: how each of the fifteen selected SNP-Lib compounds persists across the analytic choices that were tested. A descriptive tabulation of
results that were already registered and reported (ADR-AIIM-312, 399, 400, 412, 413); nothing is fitted, no new threshold or rule is introduced. It is made
after those results were seen, and every analysis that bears on the fifteen is included: none was left out because of what it showed.

Columns: conventions (of 8) and distinct splits (of 5) that selected the compound (the registered vote); whether it is still selected, within the pool,
at the threshold that admits the same share of the expansion library (0.269) and at 0.30 (Supplementary Section 17); whether the ridge Tobit model selected
it. A compound 'persists under every choice tested' when it is selected by all eight conventions, in all five distinct splits, at both thresholds and by
the ridge Tobit model."""
import csv
import json
from pathlib import Path

RUNS = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parents[1]

surv = list(csv.DictReader(open(RUNS / "em_survivors_2026-10-06" / "survivors.csv", encoding="utf-8")))
a2 = json.load(open(RUNS / "v9_sensitivity_2026-10-11" / "a1_a2_vote_and_threshold.json", encoding="utf-8"))["A2"]["thresholds"]
tobit = set(json.load(open(RUNS / "v9_censored_sensitivity_2026-10-12" / "results.json", encoding="utf-8"))["C4_ridge_tobit"]["tobit_survivors"])
t_star = set(a2["t_star"]["vote_5_of_8"]["fifteen_still_selected"])
t_030 = set(a2["0.30"]["vote_5_of_8"]["fifteen_still_selected"])

rows = []
for r in surv:
    uid = r["SNP_ID"]
    rows.append({
        "compound": uid, "name": r["name"], "conventions_of_8": int(r["conventions_of_8"]), "distinct_splits_of_5": int(r["distinct_splits_max"]),
        "selected_at_threshold_0.269": uid in t_star, "selected_at_threshold_0.30": uid in t_030, "selected_by_ridge_tobit": uid in tobit,
        "shortlist_rank": int(r["shortlist_rank"]) if r["shortlist_rank"] else None, "dropped_by_scaffold_uniqueness": r["dropped_by_scaffold_uniqueness"] == "True",
        "murcko_framework": r["murcko_framework"],
    })
for r in rows:
    r["persists_under_every_choice_tested"] = (r["conventions_of_8"] == 8 and r["distinct_splits_of_5"] == 5 and r["selected_at_threshold_0.269"]
                                              and r["selected_at_threshold_0.30"] and r["selected_by_ridge_tobit"])
persist = [r["compound"] for r in rows if r["persists_under_every_choice_tested"]]
res = {"adr": "ADR-AIIM-447", "kind": "descriptive tabulation of registered results, made after they were seen", "n_compounds": len(rows),
       "persist_under_every_choice_tested": persist, "n_selected_at_0.269": len(t_star), "n_selected_at_0.30": len(t_030), "n_selected_by_ridge_tobit_among_fifteen": len(tobit & {r["compound"] for r in rows}),
       "rows": rows}
(OUT / "persistence.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
with open(OUT / "persistence.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
for r in rows:
    print(f"{r['compound']}  {r['conventions_of_8']}/8  {r['distinct_splits_of_5']}/5  0.269:{int(r['selected_at_threshold_0.269'])}  0.30:{int(r['selected_at_threshold_0.30'])}  tobit:{int(r['selected_by_ridge_tobit'])}  short:{r['shortlist_rank']}  persist:{r['persists_under_every_choice_tested']}")
print("persist:", persist)

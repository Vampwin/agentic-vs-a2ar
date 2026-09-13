"""recount_seeded_defects.py -- the seeded-defect figures, over the right denominator.

`pilot_runner.py` defines a vacuous case and says what to do with it: an instance
whose injection parameters match nothing in the workspace changes no bytes, so
there is no fault present to escape, and scoring the resulting PASS as an escape
would penalise the contracts for a fault that was never injected. Its own
`seer_by_condition` excludes such cases.

The pilot honoured that. `pilot_postfix_2026-08-04.json` records five vacuous
instances in F13 -- F13_007, F13_008, F13_012, F13_016, F13_020 -- and reports
SEER_B = 4/128 = 0.03125.

The next day, `wp4_conditions_ac_results_postfix_2026-08-05.json` re-exported the
per-case rows to add conditions A and C. Its row schema is family_id,
instance_id, case_type, condition, action, risk, reason, escaped,
escaped_detected. The `vacuous` field is absent, and everything downstream
divided by 133: SEER_A = 9/133 = 0.067669. Every one of the five F13 "escapes"
is a case where no fault was injected.

The defect is the re-export, not the analysis that trusted it. This script is the
durable fix: it takes the conditions from the file that has all four and the
vacuity from the file that recorded it, joins them on instance_id, and recomputes
over the non-vacuous set. It emits the previously reported figures beside the
corrected ones, because a correction that improves the authors' own headline
number should be checkable against what it replaces without hunting for the
older file.

The numerator does not move. 124 cases are detected either way; only the
denominator does.

Usage:
    python recount_seeded_defects.py [--out DIR]
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from kdss_paths import project_root, resolve, under, dev_dir

# --- the rule, transcribed from pilot_runner.py before this script was run ---
# A fault case counts towards SEER when a fault was actually injected into the
# workspace: not vacuous, and not an error. This is pilot_runner.py's own rule,
# not a new one chosen here, and it is applied identically to all four
# conditions so that no arm is scored on a denominator another arm does not use.
EXCLUDE_VACUOUS = True
EXCLUDE_ERROR = True


def _root() -> Path:
    # Delegates to kdss_paths, which finds whichever layout this file is inside
    # -- the development tree or a deposit -- by walking up from its own
    # location. It previously tried two absolute paths and exited if neither
    # existed, which on the author's machine meant a copy of the deposit read
    # and wrote the development tree instead of itself.
    return project_root()


ROOT = _root()
RUNS = dev_dir("20_new_runs")
PILOT = under("20_new_runs", "WP2.5_WP3A_readiness_2026-07-30/pilot_postfix_2026-08-04.json")
WP4 = RUNS / ("WP4_retrain_postfix_2026-08-05/"
              "wp4_conditions_ac_results_postfix_2026-08-05.json")
PREVIOUS = under("20_new_runs", "verification_summary_2026-08-21/mutation_score_by_gate.json")

CONDITION_NAMES = {
    "A": "A_contracts_and_monitor",
    "B": "B_contracts_only",
    "C": "C_monitor_only",
    "D": "D_ungoverned",
}


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def vacuity_map(pilot: dict) -> dict[str, bool]:
    """instance_id -> vacuous. Verified identical across the conditions the
    pilot ran, as it must be: vacuity is a property of the injection, and the
    conditions differ only in what is switched on to catch it."""
    per_condition: dict[str, set[str]] = {}
    for o in pilot["outcomes"]:
        if o["case_type"] != "fault":
            continue
        per_condition.setdefault(o["condition"], set())
        if o.get("vacuous"):
            per_condition[o["condition"]].add(o["instance_id"])
    sets = list(per_condition.values())
    if not all(s == sets[0] for s in sets):
        raise SystemExit("vacuity differs between conditions; the join key is "
                         "not sound and this script must not guess")
    return {o["instance_id"]: bool(o.get("vacuous"))
            for o in pilot["outcomes"] if o["case_type"] == "fault"}


def counted(row: dict, vac: dict[str, bool]) -> bool:
    if row["case_type"] != "fault":
        return False
    if EXCLUDE_VACUOUS and vac.get(row["instance_id"], False):
        return False
    if EXCLUDE_ERROR and row.get("error"):
        return False
    return True


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(under("20_new_runs", "seeded_defect_recount_2026-08-29")))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    pilot = json.loads(PILOT.read_text(encoding="utf-8"))
    wp4 = json.loads(WP4.read_text(encoding="utf-8"))
    prev = json.loads(PREVIOUS.read_text(encoding="utf-8"))

    vac = vacuity_map(pilot)
    vacuous_ids = sorted(i for i, v in vac.items() if v)

    # Conditions A and C exist only in the WP4 re-export; B and D only in the
    # pilot. Both are needed, and the vacuity flags come from the pilot for all
    # four -- which is the join this correction is.
    rows: dict[str, list[dict]] = {}
    for cond in ("A", "C"):
        rows[cond] = [r for r in wp4[f"condition_{cond}_rows"]]
    for cond in ("B", "D"):
        rows[cond] = [o for o in pilot["outcomes"] if o["condition"] == cond]

    # Every arm must see the same case set, or the ablation compares denominators
    # rather than conditions -- the exact failure this script exists to undo.
    fault_ids = {c: sorted(r["instance_id"] for r in rs if r["case_type"] == "fault")
                 for c, rs in rows.items()}
    ref = fault_ids["A"]
    for c, ids in fault_ids.items():
        if ids != ref:
            raise SystemExit(f"condition {c} does not cover the same fault cases as A")

    gate_of_family = {f: v["gate"] for f, v in prev["by_family"].items()}

    seer, n_counted = {}, {}
    for cond, rs in rows.items():
        fc = [r for r in rs if counted(r, vac)]
        n_counted[cond] = len(fc)
        seer[CONDITION_NAMES[cond]] = round(
            sum(1 for r in fc if r["escaped"]) / len(fc), 6)

    A = [r for r in rows["A"] if counted(r, vac)]
    n_esc = sum(1 for r in A if r["escaped"])
    clean_A = [r for r in rows["A"] if r["case_type"] == "clean"]

    by_family, by_gate = {}, {}
    for r in A:
        fam = r["family_id"]
        f = by_family.setdefault(fam, {
            "gate": gate_of_family.get(fam),
            "description": prev["by_family"].get(fam, {}).get("description"),
            "n_fault_cases": 0, "n_escaped_condition_A": 0, "n_vacuous_excluded": 0})
        f["n_fault_cases"] += 1
        f["n_escaped_condition_A"] += int(bool(r["escaped"]))
    for fam in by_family:
        by_family[fam]["n_vacuous_excluded"] = sum(
            1 for i in vacuous_ids if i.startswith(fam + "_"))
        n, e = by_family[fam]["n_fault_cases"], by_family[fam]["n_escaped_condition_A"]
        by_family[fam]["detection_rate"] = round((n - e) / n, 4)
    for fam, f in by_family.items():
        g = by_gate.setdefault(f["gate"], {"n_fault_cases": 0, "n_escaped": 0})
        g["n_fault_cases"] += f["n_fault_cases"]
        g["n_escaped"] += f["n_escaped_condition_A"]
    for g in by_gate.values():
        g["detection_rate"] = round((g["n_fault_cases"] - g["n_escaped"])
                                    / g["n_fault_cases"], 4)

    report = {
        "what_changed": "The denominator, and only the denominator. 124 cases "
                        "are detected before and after. Five F13 instances were "
                        "counted as escapes by an export that had dropped the "
                        "vacuous flag; in each of them the injector matched "
                        "nothing and the workspace was left byte-unchanged, so "
                        "there was no fault to escape.",
        "rule": "pilot_runner.py's own: a fault case counts when a fault was "
                "actually injected. Applied identically to all four conditions.",
        "authority": "ADR-AIIM-131",
        "sources": {
            "conditions_A_and_C": {"file": str(WP4.relative_to(ROOT).as_posix()),
                                   "sha256": sha256(WP4)},
            "conditions_B_and_D_and_the_vacuity_flags": {
                "file": str(PILOT.relative_to(ROOT).as_posix()),
                "sha256": sha256(PILOT)},
            "superseded": {"file": str(PREVIOUS.relative_to(ROOT).as_posix()),
                           "sha256": sha256(PREVIOUS)},
        },
        "gate_assignment_source": prev["gate_assignment_source"],
        "vacuous_excluded": {
            "n": len(vacuous_ids),
            "instances": vacuous_ids,
            "all_in_family": sorted({i.split("_")[0] for i in vacuous_ids}),
            "meaning": "Injection parameters matched nothing in the workspace; "
                       "no bytes changed. Recorded by pilot_runner.py at "
                       "injection time and preserved in "
                       "pilot_postfix_2026-08-04.json.",
        },
        "overall": {
            "n_fault_cases_injected": len(ref),
            "n_fault_cases_counted": len(A),
            "n_escaped": n_esc,
            "detection_rate": round((len(A) - n_esc) / len(A), 4),
            "SEER_condition_A": seer["A_contracts_and_monitor"],
            "false_escalation_rate_clean_runs": round(
                sum(1 for r in clean_A if r.get("action") in ("ESCALATE", "HALT"))
                / len(clean_A), 4) if clean_A else 0.0,
            "n_clean_cases": len(clean_A),
        },
        "ablation_SEER": seer,
        "n_counted_by_condition": n_counted,
        "by_gate": by_gate,
        "by_family": by_family,
        "coverage_limitation": prev["coverage_limitation"],
        "as_previously_reported": {
            "overall": prev["overall"],
            "ablation_SEER": prev["ablation_SEER"],
            "by_gate": prev["by_gate"],
            "note": "Published in the manuscript until 2026-08-29. Divided by "
                    "133 rather than 128 because the 2026-08-05 re-export had "
                    "dropped the vacuous field.",
        },
    }

    dest = out / "mutation_score_by_gate_recounted.json"
    dest.write_text(json.dumps(report, indent=2), encoding="utf-8")

    o, p = report["overall"], prev["overall"]
    print(f"vacuous excluded : {len(vacuous_ids)}  {', '.join(vacuous_ids)}")
    print(f"cases counted    : {o['n_fault_cases_counted']} of "
          f"{o['n_fault_cases_injected']} injected")
    print()
    print(f"{'':<22}{'reported':>10}{'corrected':>12}")
    print(f"{'detection rate':<22}{p['detection_rate'] * 100:>9.1f}%"
          f"{o['detection_rate'] * 100:>11.1f}%")
    print(f"{'escaped':<22}{p['n_escaped']:>10}{o['n_escaped']:>12}")
    for k in ("A_contracts_and_monitor", "B_contracts_only",
              "C_monitor_only", "D_ungoverned"):
        print(f"{'SEER ' + k.split('_')[0]:<22}"
              f"{prev['ablation_SEER'][k]:>10.4f}{seer[k]:>12.4f}")
    print()
    print(f"{'gate':<8}{'reported':>22}{'corrected':>24}")
    for g in sorted(by_gate):
        q = prev["by_gate"][g]
        print(f"{g:<8}{q['n_fault_cases']:>6} cases {q['n_escaped']:>2} esc "
              f"{q['detection_rate'] * 100:>6.1f}%"
              f"{by_gate[g]['n_fault_cases']:>8} cases {by_gate[g]['n_escaped']:>2} esc "
              f"{by_gate[g]['detection_rate'] * 100:>6.1f}%")
    print(f"\nwritten -> {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

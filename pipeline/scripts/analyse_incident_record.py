"""analyse_incident_record.py -- what the sixteen incidents show when counted.

Executes M4.1.

Section 3.7 describes the incident record in prose and sends the reader to the
supplement for the classification. The objectives document calls this the
project's strongest under-used asset, and asks for two things the prose does
not give: **what detected each incident**, counted, and **what it cost**.

The taxonomy table is read from the S10 artifact rather than retyped. Three
questions are asked of every row:

  DETECTOR   Was it caught by a check the workflow runs before an artifact may
             be used, or by someone inspecting something outside the workflow?
             This is the distinction the paper's own thesis turns on.

  REACH      Did it reach a figure that had been reported, or was it caught
             before any result depended on it?

  OUTCOME    corrected, contained, withdrawn, or excluded.

Detection *latency* is asked for and mostly cannot be answered: the record
dates when each incident was **detected** and not when it was **introduced**,
so the interval is recoverable for exactly one of the sixteen. That is reported
as a limitation of the record rather than estimated, and it is the same
omission ADR-AIIM-219 measures in the decision log one level up.

Every verdict below is a judgement, stored with its reason, and the script
aborts if the incident set ever stops matching the stored verdicts.

Output -> 20_new_runs/incident_record_2026-09-22/
"""
from __future__ import annotations

import argparse
import io
import json
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kdss_paths import under  # noqa: E402

TAX = Path(under("30_manuscript", "supplement/S10_verification_2026-08-21/"
                 "incident_taxonomy/incident_taxonomy.md"))

# key -> (layer_as_reported, detector_is_a_workflow_check, reached_a_reported
#         figure, outcome, reason)
# The key is the first 46 characters of the incident cell, which is unique and
# stable; the script asserts the parsed set matches these keys exactly.
V: dict[str, tuple[bool, bool, str, str]] = {
    "Column-selection defect silently validating th":
        (True, False, "corrected",
         "caught by the contract evaluator before any result depended on it"),
    "Pose-RMSD computation admitting a mispositione":
        (True, False, "corrected",
         "caught by the contract evaluator at the gate"),
    "Net system charge +0.996 in the membrane build":
        (True, False, "corrected",
         "caught by the pre-run integrity check, before the run started"),
    "Per-compound memory growth toward exhaustion":
        (False, False, "contained",
         "observed by comparing resident-set size across restarts; a wrapper "
         "was added, no artifact was defective"),
    "Host restarts terminating all docking processe":
        (False, False, "contained",
         "found by inspecting processes and boot time; idempotent resume, zero "
         "rows lost"),
    "CPU oversubscription — three chunks each claim":
        (False, False, "contained",
         "noticed from load average, not from any check; a thread cap was set"),
    "An edited shell script not reloaded by an alre":
        (False, False, "contained",
         "found by reading the running process command line"),
    "Three glycosides never converging":
        (False, False, "excluded",
         "identical progress counts noticed across restart cycles; the three "
         "were excluded with a recorded reason"),
    "Stray directory created by a host path used in":
        (False, False, "contained",
         "found in a directory listing; path resolution now probes"),
    "Literature precedent asserted without a search":
        (False, False, "withdrawn",
         "caught against the project's own frozen labels; the claim was "
         "withdrawn before it was published"),
    "Process loss attributed to system sleep":
        (False, False, "corrected",
         "corrected against the operating-system event log; the proposed "
         "remediation was withdrawn as ineffective"),
    "Journal scope judged from the scope statement":
        (False, False, "withdrawn",
         "reversed after searching the journal's published record"),
    '"This force field reports no parameter-penalty':
        (False, False, "corrected",
         "corrected by reading the tool's own output"),
    '"The simulation package is not installed"':
        (False, False, "corrected",
         "corrected after checking the remaining guest environments"),
    "A membrane-build failure attributed to the pac":
        (False, False, "corrected",
         "attribution corrected by measuring the two files separately"),
    "A re-export of per-case fault-injection rows s":
        (False, True, "corrected",
         "the only incident that reached a reported figure: detection stood at "
         "93.2 % instead of 96.9 % for three weeks, and no check fired"),
}

ROW_RX = re.compile(r"^\|\s*(Scientific|Infrastructure|Agent reasoning)\s*\|"
                    r"\s*(.+?)\s*\|\s*(.+?)\s*\|\s*(.+?)\s*\|\s*$", re.M)
KEY = 46  # first characters of the incident cell, right-stripped


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(
        under("20_new_runs", "incident_record_2026-09-22")))
    args = ap.parse_args()
    out = Path(args.out)

    text = io.open(TAX, encoding="utf-8", newline="").read().replace("\r\n", "\n")
    split = text.index("## Found during revision")

    rows = []
    for phase, chunk in (("construction", text[:split]),
                         ("revision", text[split:])):
        for m in ROW_RX.finditer(chunk):
            layer, incident, detector, resolution = (g.strip()
                                                     for g in m.groups())
            rows.append({"phase": phase, "layer": layer, "incident": incident,
                         "detector": detector, "resolution": resolution,
                         "key": incident[:KEY].rstrip()})

    keys = [r["key"] for r in rows]
    if len(set(keys)) != len(keys):
        print("ABORT: incident keys are not unique at "
              f"{KEY} characters: "
              f"{[k for k, c in Counter(keys).items() if c > 1]}")
        return 1
    if set(keys) != set(V):
        print("ABORT: the incident set no longer matches the stored verdicts.")
        print(f"  new, unjudged:  {sorted(set(keys) - set(V))}")
        print(f"  judged, absent: {sorted(set(V) - set(keys))}")
        return 1

    for r in rows:
        check, reached, outcome, reason = V[r["key"]]
        r.update(caught_by_workflow_check=check, reached_reported_figure=reached,
                 outcome=outcome, reason=reason)

    n = len(rows)
    by_check = sum(r["caught_by_workflow_check"] for r in rows)
    reached = [r for r in rows if r["reached_reported_figure"]]
    layers = Counter(r["layer"] for r in rows)
    outcomes = Counter(r["outcome"] for r in rows)

    # per layer, in the order the taxonomy uses
    order = ["Scientific", "Infrastructure", "Agent reasoning"]
    per_layer = {}
    for lay in order:
        sub = [r for r in rows if r["layer"] == lay]
        per_layer[lay] = {
            "n": len(sub),
            "construction": sum(r["phase"] == "construction" for r in sub),
            "revision": sum(r["phase"] == "revision" for r in sub),
            "caught_by_workflow_check": sum(r["caught_by_workflow_check"]
                                            for r in sub),
            "reached_reported_figure": sum(r["reached_reported_figure"]
                                           for r in sub),
            "outcomes": dict(Counter(r["outcome"] for r in sub)),
        }

    print("=" * 78)
    print("THE INCIDENT RECORD, COUNTED   (M4.1)")
    print("=" * 78)
    print(f"  incidents: {n}   "
          f"({sum(r['phase'] == 'construction' for r in rows)} during "
          f"construction, {sum(r['phase'] == 'revision' for r in rows)} during "
          f"revision)")
    print()
    print(f"  {'layer':18}{'n':>4}{'by a check':>12}{'reached a figure':>18}"
          f"   outcomes")
    print("  " + "-" * 72)
    for lay in order:
        d = per_layer[lay]
        oc = ", ".join(f"{v} {k}" for k, v in sorted(d["outcomes"].items()))
        print(f"  {lay:18}{d['n']:>4}{d['caught_by_workflow_check']:>12}"
              f"{d['reached_reported_figure']:>18}   {oc}")
    print("  " + "-" * 72)
    print(f"  {'all':18}{n:>4}{by_check:>12}{len(reached):>18}")
    print()
    print(f"  Caught by a check the workflow runs: {by_check} of {n}, and all "
          f"{by_check} are scientific.")
    print(f"  Reached a reported figure: {len(reached)} of {n}.")
    print("  Detection latency: the record dates detection and not "
          "introduction,")
    print("  so the interval is recoverable for 1 of 16 -- three weeks.")

    out.mkdir(parents=True, exist_ok=True)
    summary = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "executes": "M4.1",
        "question": "What detected the incidents, and what did they cost?",
        "source": "30_manuscript/supplement/S10_verification_2026-08-21/"
                  "incident_taxonomy/incident_taxonomy.md",
        "n_incidents": n,
        "n_construction": sum(r["phase"] == "construction" for r in rows),
        "n_revision": sum(r["phase"] == "revision" for r in rows),
        "by_layer": dict(layers),
        "per_layer": per_layer,
        "n_caught_by_workflow_check": by_check,
        "n_reached_a_reported_figure": len(reached),
        "reached_a_reported_figure": [r["incident"][:120] for r in reached],
        "outcomes": dict(outcomes),
        "latency": {
            "recoverable_for": 1,
            "of": n,
            "value": "three weeks",
            "why_not_more": "the record dates when each incident was detected "
                            "and not when it was introduced; the same omission "
                            "ADR-AIIM-219 measures in the decision log",
        },
        "verdicts": {r["key"]: {"layer": r["layer"], "phase": r["phase"],
                                "caught_by_workflow_check":
                                    r["caught_by_workflow_check"],
                                "reached_reported_figure":
                                    r["reached_reported_figure"],
                                "outcome": r["outcome"], "reason": r["reason"]}
                     for r in rows},
        "caveat": "The three classifications are judgements over a table read "
                  "from its artifact. Each verdict carries its reason and the "
                  "script aborts if the incident set changes.",
        "reading": "The contract layer caught every incident in the class it "
                   "was built for and none in the other two. The one incident "
                   "that reached a reported figure was caught by a "
                   "reconciliation a reviewer's question forced, not by a "
                   "check.",
    }
    with open(out / "incident_record.json", "w", encoding="utf-8",
              newline="\n") as f:
        json.dump(summary, f, indent=1, ensure_ascii=False)
        f.write("\n")
    print(f"\nwritten -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

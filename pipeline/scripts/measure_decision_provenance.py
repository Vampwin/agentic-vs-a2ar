"""measure_decision_provenance.py -- what the decision record says about what was not chosen.

Executes M1.5.

The five-objective design document asserted that the decision log "contains the
word `Alternatives` zero times". That assertion was written from a single grep
and is **false**: the word occurs three times, all inside ADR-AIIM-082, which
did screen an algorithm against three alternatives on identical nested
cross-validation. The assertion is corrected here by measurement, and the
finding that survives is narrower and better supported than the one claimed.

Three counts, in decreasing order of how mechanical they are:

  STRUCTURED  Does an ADR carry a *field* naming the options available and not
              taken -- a label whose absence a check could notice, the way
              `Status:` and `Decision:` can be checked for? Purely mechanical.

  GENEROUS    Does the ADR's prose mention a rejected option anywhere, under a
              deliberately wide keyword net? Mechanical, and it over-counts on
              purpose.

  READ        Of the generous hits, which name an alternative *course of
              action*, rather than a rejected null hypothesis, a discarded
              intermediate file, or compounds rejected by a filter? This one is
              a judgement. Every verdict is stored with its reason, and the
              script fails loudly if the set of generous hits ever stops
              matching the set of stored verdicts, so the judgement cannot go
              stale silently.

Output -> 20_new_runs/decision_provenance_2026-09-22/
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

LOG = Path(under("00_project_control", "08_DECISION_LOG.md"))

# A structured field is a bolded label opening a line.
FIELD_RX = re.compile(r"^\*\*([A-Z][A-Za-z /-]{2,40})[:*]", re.M)
HEAD_RX = re.compile(r"^#{2,4}\s*(.+?)\s*$", re.M)
ADR_RX = re.compile(r"^## (ADR-AIIM-(\d+)[^\n]*)$", re.M)

ALT_WORDS = ("alternativ", "option", "rejected", "considered",
             "discarded", "counterfactual", "not chosen", "ruled out")

# Deliberately wide. It over-counts, and the READ pass below says by how much.
GENEROUS = [
    r"\balternativ\w*", r"\bconsidered\b", r"\brejected\b", r"\bdiscarded\b",
    r"\bruled out\b", r"\bin preference to\b", r"\bchosen over\b",
    r"\bover the\b.{0,20}\boption", r"\bthe other option\b",
    r"\boptions?\s+(were|was|are|considered)\b", r"\bwhy not\b",
]
GEN_RX = [re.compile(p, re.I) for p in GENEROUS]

# The hand reading. True = the ADR names an alternative course of action.
# Recorded per ADR with its reason so the judgement is auditable, and asserted
# against the generous set so it cannot drift.
READ: dict[str, tuple[bool, str]] = {
    "023": (True, "records the permitted model-candidate set, gradient boosting as the single admitted alternative"),
    "032": (False, "'H0 rejected at J=11' -- a statistical rejection, not a design choice"),
    "038": (True, "carries an explicit 'What was explicitly rejected, and why' list"),
    "041": (True, "names a plausible alternative reading of the evidence and why the prespecified rule was kept"),
    "048": (True, "redoing ten trajectories for provenance explicitly ruled out on cost"),
    "049": (False, "'rejected-design framing' describes historical text, not a choice made here"),
    "073": (False, "'discarded poses' -- an intermediate file, not an option"),
    "075": (False, "states variants are 'not ruled out', i.e. records the absence of a decision"),
    "082": (True, "RandomForest screened against three named alternatives on identical nested CV"),
    "083": (True, "an alternative label threshold tested against the frozen one and compared"),
    "089": (False, "'ruled out, with evidence' concerns causes of a WSL2 crash, not a design option"),
    "090": (True, "splicing chosen over generating afresh, with the cost of the rejected route stated"),
    "091": (True, "the venue compared against the alternatives surveyed"),
    "096": (True, "flagging the missing citation weighed against guessing it"),
    "101": (True, "a single-replicate MD result 'considered and rejected'"),
    "105": (True, "records which highlights were omitted deliberately rather than missed"),
    "107": (True, "'why this one, and why not the alternative the author also liked'"),
    "111": (False, "a 0.0028 improvement rejected against a pre-set threshold -- a result, not an option"),
    "115": (False, "'a considered decision at the time' with no option named"),
    "118": (False, "'nothing already computed is discarded' -- about data retention"),
    "124": (True, "the alternative directory layout named, considered and rejected with reasons"),
    "125": (False, "'computed and discarded' at a source line -- an intermediate value"),
    "126": (False, "same discarded intermediate value, now retained"),
    "134": (False, "'nothing was recomputed or discarded' -- about run resumption"),
    "138": (True, "states the alternative to the recommendation and why it is worse"),
    "140": (False, "'before applicability is considered' -- ordinary prose use"),
    "146": (True, "the alternative of reporting one candidate at half length, named and rejected"),
    "147": (True, "a section headed 'the two panels rejected because of it', each with a reason"),
    "150": (False, "compounds rejected by a screening filter, not options"),
    "152": (True, "a section headed 'Why not', answering a naming alternative the author raised"),
    "158": (True, "'a tempting alternative was rejected before it was written'"),
    "161": (True, "names which parts of a reference protocol were adopted and which were not used"),
    "167": (False, "'self_test() rejected a category' -- a guard firing, not a choice"),
    "169": (True, "re-running the nulls named as a candidate fix and rejected with reasons"),
    "170": (True, "revisits techniques previously discarded and states what each was judged on"),
    "194": (True, "records that the reporting options were put to the author"),
    "198": (True, "a section headed 'Three innocent explanations, checked and ruled out'"),
    "200": (False, "'reviewed and rejected' describes the manuscript's fate, not an option"),
}


def split_adrs(text: str, through: int) -> list[tuple[str, str, str]]:
    """Every ADR up to `through`, so the count does not move as the log grows.

    Without the cut this measurement changes every time a decision is recorded
    -- including the one that records this measurement -- and a manuscript
    number computed from it would go stale on the next ADR. The cut makes the
    population "every decision written before this analysis", which is both
    stable and the population the claim is actually about.
    """
    marks = [(m.start(), m.group(1), m.group(2)) for m in ADR_RX.finditer(text)]
    out = []
    for i, (pos, title, num) in enumerate(marks):
        end = marks[i + 1][0] if i + 1 < len(marks) else len(text)
        if int(num) <= through:
            out.append((num, title, text[pos:end]))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(
        under("20_new_runs", "decision_provenance_2026-09-22")))
    ap.add_argument("--through", type=int, default=218,
                    help="last ADR number in the measured population; 218 is "
                         "the last recorded before this analysis was written")
    args = ap.parse_args()
    out = Path(args.out)

    text = io.open(LOG, encoding="utf-8", newline="").read().replace("\r\n", "\n")
    adrs = split_adrs(text, args.through)
    n = len(adrs)

    fields: Counter[str] = Counter()
    for _num, _title, body in adrs:
        for m in FIELD_RX.finditer(body):
            fields[m.group(1).strip()] += 1
    alt_fields = {k: v for k, v in fields.items()
                  if any(w in k.lower() for w in ALT_WORDS)}

    alt_headings = []
    for num, _title, body in adrs:
        for m in HEAD_RX.finditer(body):
            h = m.group(1).lower()
            if h.startswith("adr-aiim"):
                continue
            if any(w in h for w in ALT_WORDS):
                alt_headings.append([num, m.group(1)])

    generous = []
    for num, title, body in adrs:
        quotes = []
        for rx in GEN_RX:
            for m in rx.finditer(body):
                s, e = max(0, m.start() - 90), min(len(body), m.end() + 90)
                quotes.append(re.sub(r"\s+", " ", body[s:e]).strip())
        if quotes:
            generous.append((num, title, quotes))

    gen_ids = {g[0] for g in generous}
    if gen_ids != set(READ):
        print("ABORT: the generous hit set no longer matches the stored reading.")
        print(f"  new, unread:   {sorted(gen_ids - set(READ))}")
        print(f"  read, now absent: {sorted(set(READ) - gen_ids)}")
        return 1

    named = sorted(k for k, (v, _r) in READ.items() if v)
    incidental = sorted(k for k, (v, _r) in READ.items() if not v)

    print("=" * 78)
    print("WHAT THE DECISION RECORD SAYS ABOUT WHAT WAS NOT CHOSEN   (M1.5)")
    print("=" * 78)
    print(f"  decision records: {n}")
    print("  the fields they use:")
    for k, v in fields.most_common(8):
        print(f"      {v:4}  **{k}**")
    print()
    print(f"  STRUCTURED  fields naming options not taken: {len(alt_fields)}")
    print(f"              sections so headed:              {len(alt_headings)}"
          f"  {[h[0] for h in alt_headings]}")
    print(f"  GENEROUS    ADRs mentioning a rejected option anywhere: "
          f"{len(generous)} of {n}  ({100*len(generous)/n:.1f} %)")
    print(f"  READ        of those, naming an alternative course of action: "
          f"{len(named)}  ({100*len(named)/n:.1f} % of all)")
    print(f"              incidental use of the same words: {len(incidental)}")
    print()
    print("  Selection quality is auditable through outcomes, not through")
    print("  process, and the omission is undetectable because no field")
    print("  exists whose absence a check could notice.")

    out.mkdir(parents=True, exist_ok=True)
    summary = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "executes": "M1.5",
        "question": "Does the decision record say what was decided against?",
        "source": "00_project_control/08_DECISION_LOG.md",
        "population": f"every ADR up to and including ADR-AIIM-{args.through}, "
                      f"i.e. every decision recorded before this analysis; the "
                      f"cut keeps the count from moving as the log grows",
        "through_adr": args.through,
        "corrects": "the five-objective design document's claim that the log "
                    "contains the word 'Alternatives' zero times; it occurs "
                    "three times, all in ADR-AIIM-082",
        "n_decision_records": n,
        "field_inventory": dict(fields.most_common()),
        "n_structured_alternatives_fields": len(alt_fields),
        "structured_alternatives_fields": alt_fields,
        "sections_headed_for_rejected_options": alt_headings,
        "generous": {
            "patterns": GENEROUS,
            "n": len(generous),
            "fraction": round(len(generous) / n, 4),
            "ids": sorted(gen_ids),
            "quotes": {g[0]: g[2][:4] for g in generous},
        },
        "read": {
            "n_naming_an_alternative_course_of_action": len(named),
            "n_incidental": len(incidental),
            "ids_naming_an_alternative": named,
            "ids_incidental": incidental,
            "verdicts": {k: {"names_an_alternative": v, "reason": r}
                         for k, (v, r) in sorted(READ.items())},
            "caveat": "This pass is a judgement, not a measurement. Every "
                      "verdict carries its reason and the script aborts if "
                      "the generous hit set changes.",
        },
        "reading": "The record has fields for status, date, approver, decision "
                   "and consequence, and none for the options available and "
                   "not taken. A record with no field for something cannot "
                   "report that it is missing.",
    }
    with open(out / "decision_provenance.json", "w", encoding="utf-8",
              newline="\n") as f:
        json.dump(summary, f, indent=1, ensure_ascii=False)
        f.write("\n")
    print(f"\nwritten -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

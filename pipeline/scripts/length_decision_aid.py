"""length_decision_aid.py -- what each section costs, measured not remembered.

`LENGTH_AUDIT_2026-09-22.md` was written by hand and was wrong within days: its
route A claimed 3,700 relocatable words, which re-measurement put at 1,798, then
1,530, then 1,169 once figure captions were separated from prose and sections
carrying load-bearing Methods were excluded (ADR-AIIM-251, -252). A hand-written
audit of a moving manuscript goes stale silently.

This generates the same decision aid from the manuscript itself, so it cannot.
Re-run it after any edit.

For every section it reports:

    words       total, and split into prose and figure/table material, because
                a figure caption cannot be relocated without its figure
    anchors     audited G8 literals living there, split the same way. This is
                the real cost of moving a section: each anchor that leaves the
                article must be declared in `IN_SUPPLEMENT` and verified in the
                built supplement (ADR-AIIM-253), and each must be checked absent
                from the article AND present in the supplement (ADR-AIIM-254).
    relocatable prose words with no anchor in them -- the only words that move
                for free

It also prints the journal's own distribution, measured in ADR-AIIM-251, so a
target can be chosen against what AIIM publishes rather than against a
percentage of whatever the manuscript happened to weigh.

**It recommends nothing.** Which sections may leave the article is a scientific
judgement for the author; this supplies the arithmetic that judgement needs.

Output -> 30_manuscript/CURRENT_v5/LENGTH_DECISION_AID.md
"""
from __future__ import annotations

import io
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kdss_paths import under  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

CLEAN = Path(under("30_manuscript", "CURRENT_v5/manuscript_AIIM_v6_CLEAN.md"))
G8_REPORT = Path(under("20_new_runs", "verification_summary_2026-08-21/"
                       "g8_number_traceability.json"))
AIIM = Path(under("20_new_runs", "aiim_length_reference_2026-09-24/"
                  "aiim_length_reference.json"))
OUT = Path(under("30_manuscript", "CURRENT_v5/LENGTH_DECISION_AID.md"))

SENTINEL = ("figcopies_identical", "figpng_not_stale", "supp_sections_resolve",
            "supp_reference_lists_resolve", "deposit_url_named")


def classify(block: list[str]) -> tuple[str, str]:
    """(figure/table text, prose text) for one section body."""
    fig, pro, in_cap = [], [], False
    for l in block:
        s = l.strip()
        if s.startswith("!["):
            fig.append(l); in_cap = True; continue
        if s.startswith("|"):
            fig.append(l); continue
        if re.match(r"^\*\*(Figure|Table)\s", s):
            fig.append(l); in_cap = True; continue
        if in_cap:
            if not s:
                in_cap = False
            fig.append(l); continue
        pro.append(l)
    return "\n".join(fig), "\n".join(pro)


def main() -> int:
    text = io.open(CLEAN, encoding="utf-8", newline="").read().replace("\r\n", "\n")
    lines = text.split("\n")

    rep = json.load(io.open(G8_REPORT, encoding="utf-8"))
    lits = [(c["check"], c["literal"]) for c in rep["checks"]
            if isinstance(c.get("literal"), str) and c["literal"].strip()
            and not c["check"].startswith(("companion:", "figprov:"))
            and c["literal"] not in SENTINEL and c["check"] not in SENTINEL]
    in_supp = sum(1 for c in rep["checks"] if c.get("located_in") == "supplement")

    heads = [(i, l) for i, l in enumerate(lines) if re.match(r"^#{1,4}\s", l)]
    rows = []
    for k, (i, l) in enumerate(heads):
        j = heads[k + 1][0] if k + 1 < len(heads) else len(lines)
        body = lines[i + 1:j]
        ft, pt = classify(body)
        title = re.sub(r"^#+\s*", "", l).strip()
        # prose paragraphs holding no audited literal move for free
        free = sum(len(p.split()) for p in pt.split("\n\n")
                   if p.strip() and not any(re.sub(r"\s+", " ", x) in
                                            re.sub(r"\s+", " ", p) for _, x in lits))
        rows.append({
            "title": title, "depth": len(re.match(r"^(#+)", l).group(1)),
            "total": len(" ".join(body).split()) + len(l.split()),
            "prose": len(pt.split()), "figtab": len(ft.split()),
            "anchors_prose": sum(1 for _, x in lits
                                 if re.sub(r"\s+", " ", x) in re.sub(r"\s+", " ", pt)),
            "anchors_figtab": sum(1 for _, x in lits
                                  if re.sub(r"\s+", " ", x) in re.sub(r"\s+", " ", ft)),
            "relocatable": free,
        })

    body_rows = [r for r in rows if not r["title"].lower().startswith(
        ("reference", "declaration", "title", "authors"))]
    refs = next((r["total"] for r in rows
                 if r["title"].lower().startswith("reference")), 0)
    total = len(text.split())

    a = json.load(io.open(AIIM, encoding="utf-8")) if AIIM.is_file() else None

    L = ["# Length decision aid", "",
         f"*Generated {datetime.now(timezone.utc).isoformat(timespec='seconds')} "
         f"by `length_decision_aid.py` from the manuscript itself. "
         f"**Re-run it after any edit** — the hand-written "
         f"`LENGTH_AUDIT_2026-09-22.md` went stale within days and its route-A "
         f"figure was wrong by a factor of three (ADR-AIIM-251, -252).*", "",
         "## Where the paper stands", "",
         f"- **{total:,}** words total",
         f"- **{total - refs:,}** excluding the reference list ({refs:,})",
         f"- **{sum(r['figtab'] for r in body_rows):,}** of the body is figure "
         f"captions and tables, which cannot move without their figure", ""]

    if a:
        d = a["reference_distribution_words_excluding_refs"]
        L += ["## What this journal publishes", "",
              f"Measured from **{a['n_measured']}** AIIM papers, excluding "
              f"references (ADR-AIIM-251). **AIIM states no word or page limit "
              f"on the article body**; the only stated cap is 250 words on the "
              f"abstract.", "",
              "| | words |", "|---|---|",
              f"| median | {d['median']:,} |",
              f"| 75th percentile | {d['p75']:,} |",
              f"| 90th percentile | {d['p90']:,} |",
              f"| **maximum observed** | **{d['max']:,}** |",
              f"| **this manuscript** | **{total - refs:,}** |", "",
              f"A target below **{d['max']:,}** is the first one that puts the "
              f"paper inside the range this journal has actually printed.", ""]

    L += ["## Cost of cutting, by section", "",
          "`relocatable` is prose carrying **no** audited number: it moves for "
          "free. Everything else costs anchor work — each audited literal that "
          "leaves the article must be declared in `IN_SUPPLEMENT`, verified "
          "present in the built supplement, and verified absent from the "
          "article (ADR-AIIM-253, -254).", "",
          "| section | total | prose | fig/tab | anchors | relocatable |",
          "|---|---|---|---|---|---|"]
    for r in body_rows:
        ind = "&nbsp;&nbsp;" * max(0, r["depth"] - 2)
        anc = r["anchors_prose"] + r["anchors_figtab"]
        L.append(f"| {ind}{r['title'][:56]} | {r['total']:,} | {r['prose']:,} | "
                 f"{r['figtab']:,} | {anc} | {r['relocatable']:,} |")
    L += ["",
          f"**Total relocatable at zero anchor cost: "
          f"{sum(r['relocatable'] for r in body_rows):,} words.** Anything "
          f"beyond that moves audited numbers with it.", "",
          f"**{in_supp}** checks already live in the supplement from earlier "
          f"relocations, so the mechanism is in place and exercised.", "",
          "## What this does not say", "",
          "Which sections may leave the article. That is a scientific "
          "judgement — an exploratory analysis and a load-bearing control cost "
          "the same number of words and are not interchangeable. This file "
          "supplies arithmetic, not priorities.", ""]

    OUT.write_text("\n".join(L), encoding="utf-8")
    print(f"{total:,} words, {total - refs:,} excluding references")
    print(f"relocatable at zero anchor cost: "
          f"{sum(r['relocatable'] for r in body_rows):,}")
    print(f"already in the supplement: {in_supp} checks")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

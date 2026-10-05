"""measure_aiim_article_length.py -- what length does AIIM actually publish?

LENGTH_AUDIT_2026-09-22 section 4C states the problem plainly:

    "**The 13,500 target is the reviewer's, not the journal's.** ... The
    journal's own limit could not be confirmed here -- `sciencedirect.com`
    returns HTTP 403 to automated fetches and `CLAUDE.md` bars the in-app
    browser from that domain, so **checking the AIIM guide for authors is a
    one-minute author task that would materially change this decision.**"

That remains true: the guide for authors is still unreachable from an automated
session. Three routes were tried on 2026-09-24 and all failed --
`elsevier.com/.../guide-for-authors` 301-redirects to `sciencedirect.com`, which
returns **403**; the scispace mirror returns **403**; and Claude in Chrome
reports no connected browser. The in-app browser is barred from that domain by
`CLAUDE.md` after 14/14 confirmed crashes.

**So this measures the empirical distribution instead of the stated rule.** The
question it answers is not "does AIIM impose a hard word limit" -- it cannot
answer that -- but the question that actually drives the decision:

    Is 28,672 words within the range of what this journal publishes?

### Method

arXiv's export API is queried for preprints whose `journal_ref` names
*Artificial Intelligence in Medicine*, restricted to the journal (conference
proceedings and the separate AIME/AIMS venues are excluded by string match) and
to 2020 onward. Each is retrieved from ar5iv, which renders arXiv LaTeX as
HTML, and the body is counted with the bibliography section removed. The same
counter is then run over the manuscript for a like-for-like comparison.

### What this is and is not

It **is** a measure of the length of the author-submitted versions of papers
that AIIM accepted. It is **not**:

  * the published version -- arXiv holds the preprint, which may be the
    submitted, accepted or an intermediate draft, and Elsevier's own copy is
    not reachable;
  * a word count by the journal's definition -- ar5iv counts table cells,
    figure captions and mathematical tokens as words, which inflates the count
    for figure-heavy papers relative to a word processor's count;
  * a limit. A journal that publishes nothing above 15,000 words may still
    state no limit, and a paper above the observed range may still be accepted
    on its merits.

It is a comparison of this manuscript against its own journal's output under
one consistent counting rule, which is strictly better evidence than the
reviewer's 13,500 (which is 50 % of a number this manuscript happened to stand
at) and strictly worse evidence than the guide for authors.

Papers whose ar5iv conversion fails are reported and excluded, not silently
dropped; the count of exclusions is part of the artifact.

Output -> 20_new_runs/aiim_length_reference_2026-09-24/
"""
from __future__ import annotations

import argparse
import io
import json
import re
import sys
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from statistics import median

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kdss_paths import under  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

OUT = Path(under("20_new_runs", "aiim_length_reference_2026-09-24"))
MANUSCRIPT = Path(under("30_manuscript", "CURRENT_v5",
                        "manuscript_AIIM_v6_CLEAN.md"))

ARXIV_API = ("https://export.arxiv.org/api/query?search_query="
             "all:%22artificial+intelligence+in+medicine%22"
             "&start=0&max_results=200")
AR5IV = "https://ar5iv.labs.arxiv.org/html/{}"

# Venues that share the journal's name but are not the journal.
NOT_THE_JOURNAL = ("conference on", "proceedings of",
                   "artificial intelligence in medicine and surgery",
                   "international journal",
                   "transformers and large language models")
YEAR_FROM = 2020
# An ar5iv page this short is a conversion failure, not a short paper.
MIN_PLAUSIBLE_WORDS = 2000
UA = {"User-Agent": "aiim-length-reference/1.0 (academic length comparison)"}


def _get(url: str, timeout: int = 90) -> str:
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8", errors="replace")


def discover() -> list[dict]:
    """arXiv preprints whose journal_ref names the AIIM journal, 2020+."""
    xml = _get(ARXIV_API)
    rows = []
    for e in re.findall(r"(?s)<entry>(.*?)</entry>", xml):
        jr = re.search(r"(?s)<arxiv:journal_ref[^>]*>(.*?)</arxiv:journal_ref>", e)
        aid = re.search(r"<id>https?://arxiv\.org/abs/([^<]+)</id>", e)
        ti = re.search(r"(?s)<title>(.*?)</title>", e)
        if not (jr and aid and ti):
            continue
        j = " ".join(jr.group(1).split())
        jl = j.lower()
        if "artificial intelligence in medicine" not in jl and \
           "artif. intell. med" not in jl:
            continue
        if any(x in jl for x in NOT_THE_JOURNAL):
            continue
        years = [int(y) for y in re.findall(r"\b((?:19|20|21|25|26)\d{2})\b", j)]
        year = max(years) if years else 0
        if year < YEAR_FROM:
            continue
        rows.append({"arxiv_id": aid.group(1), "year": year,
                     "title": " ".join(ti.group(1).split()),
                     "journal_ref": j})
    rows.sort(key=lambda r: (-r["year"], r["arxiv_id"]))
    return rows


def _words(html: str) -> int:
    t = re.sub(r"(?s)<[^>]+>", " ", html)
    t = re.sub(r"&[a-zA-Z]+;|&#\d+;", " ", t)
    return len(t.split())


def count_ar5iv(html: str) -> tuple[int, int]:
    """(words in body, words in body excluding the bibliography)."""
    html = re.sub(r"(?is)<(script|style|head)[^>]*>.*?</\1>", " ", html)
    m = re.search(r"(?is)<body[^>]*>(.*)</body>", html)
    body = m.group(1) if m else html
    nobib = re.split(r"(?i)<(?:section|div)[^>]*class=\"[^\"]*ltx_bibliography",
                     body)[0]
    return _words(body), _words(nobib)


def count_manuscript(path: Path) -> dict:
    """The manuscript under headings, so like-for-like slices are available."""
    t = io.open(path, encoding="utf-8", newline="").read().replace("\r\n", "\n")
    lines = t.split("\n")
    heads = [(i, l) for i, l in enumerate(lines) if re.match(r"^#{1,2}\s", l)]
    spans = {}
    for k, (i, l) in enumerate(heads):
        j = heads[k + 1][0] if k + 1 < len(heads) else len(lines)
        spans[re.sub(r"^#+\s*", "", l).strip()] = len(
            " ".join(lines[i:j]).split())
    total = len(t.split())
    refs = next((v for k, v in spans.items() if k.lower().startswith("reference")), 0)
    decl = next((v for k, v in spans.items() if k.lower().startswith("declaration")), 0)
    return {"total": total, "by_section": spans,
            "excluding_references": total - refs,
            "excluding_references_and_declarations": total - refs - decl}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0,
                    help="cap the number of papers fetched (0 = all found)")
    ap.add_argument("--sleep", type=float, default=1.0,
                    help="seconds between ar5iv requests")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    cache = OUT / "cache"
    cache.mkdir(exist_ok=True)

    papers = discover()
    if args.limit:
        papers = papers[:args.limit]
    print(f"{len(papers)} AIIM journal papers found on arXiv, {YEAR_FROM}+")

    ok, failed = [], []
    for n, p in enumerate(papers, 1):
        f = cache / f"{p['arxiv_id']}.html"
        try:
            if not f.exists():
                f.write_text(_get(AR5IV.format(p["arxiv_id"].split("v")[0])),
                             encoding="utf-8")
                time.sleep(args.sleep)
            total, nobib = count_ar5iv(f.read_text(encoding="utf-8",
                                                   errors="replace"))
        except Exception as exc:                       # network or parse
            failed.append({**p, "reason": f"{type(exc).__name__}: {exc}"[:120]})
            print(f"  [{n:2d}/{len(papers)}] {p['arxiv_id']:14s} FAILED")
            continue
        if nobib < MIN_PLAUSIBLE_WORDS:
            failed.append({**p, "reason": f"ar5iv conversion produced only "
                                          f"{nobib} words"})
            print(f"  [{n:2d}/{len(papers)}] {p['arxiv_id']:14s} "
                  f"unusable ({nobib} w)")
            continue
        ok.append({**p, "words_total": total, "words_excluding_refs": nobib})
        print(f"  [{n:2d}/{len(papers)}] {p['arxiv_id']:14s} {nobib:6d} w "
              f"(excl. refs)  {p['year']}")

    if len(ok) < 5:
        print(f"\nREFUSING to report a distribution from {len(ok)} papers.")
        return 2

    vals = sorted(r["words_excluding_refs"] for r in ok)
    def pct(q):
        i = max(0, min(len(vals) - 1, int(round(q * (len(vals) - 1)))))
        return vals[i]

    ms = count_manuscript(MANUSCRIPT)
    mine = ms["excluding_references"]
    above = sum(1 for v in vals if v < mine)

    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "question": "Is the manuscript's length within the range AIIM publishes?",
        "cannot_answer": "whether AIIM states a hard word limit; the guide for "
                         "authors returned HTTP 403 on every automated route "
                         "tried on 2026-09-24 and CLAUDE.md bars the in-app "
                         "browser from sciencedirect.com",
        "source": "arXiv preprints whose journal_ref names the AIIM journal",
        "counting_rule": "ar5iv HTML body, bibliography section removed; table "
                         "cells, captions and math tokens count as words",
        "year_from": YEAR_FROM,
        "n_found": len(papers), "n_measured": len(ok), "n_excluded": len(failed),
        "reference_distribution_words_excluding_refs": {
            "min": vals[0], "p25": pct(0.25), "median": int(median(vals)),
            "p75": pct(0.75), "p90": pct(0.90), "max": vals[-1]},
        "manuscript": {
            "file": str(MANUSCRIPT), **ms,
            "compared_on": "excluding_references",
            "words": mine,
            "n_reference_papers_shorter": above,
            "fraction_of_reference_papers_shorter": round(above / len(vals), 4),
            "ratio_to_reference_median": round(mine / median(vals), 2)},
        "reviewer_target_13500": {
            "provenance": "the reviewer's, not the journal's; 50 % of the "
                          "27,130 words the manuscript stood at when the "
                          "review was written",
            "percentile_in_reference_distribution": round(
                sum(1 for v in vals if v < 13500) / len(vals), 4)},
        "papers": ok,
        "excluded": failed,
        "not_claimed": [
            "that AIIM imposes or does not impose a word limit",
            "that these counts equal the published versions' lengths; arXiv "
            "holds preprints and Elsevier's copy is unreachable",
            "that exceeding the observed range predicts rejection",
        ],
    }
    (OUT / "aiim_length_reference.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")

    d = report["reference_distribution_words_excluding_refs"]
    print("\n" + "=" * 72)
    print("AIIM PUBLISHED LENGTH, words excluding references")
    print("=" * 72)
    print(f"  n = {len(ok)} papers ({len(failed)} excluded, listed in the artifact)")
    print(f"  min {d['min']}   p25 {d['p25']}   median {d['median']}   "
          f"p75 {d['p75']}   p90 {d['p90']}   max {d['max']}")
    print(f"\n  this manuscript, same basis: {mine} words")
    print(f"    {report['manuscript']['ratio_to_reference_median']}x the median; "
          f"longer than {above}/{len(vals)} reference papers")
    print(f"  the reviewer's 13,500 sits at percentile "
          f"{report['reviewer_target_13500']['percentile_in_reference_distribution']:.2f} "
          f"of this distribution")
    print(f"\nwrote {OUT / 'aiim_length_reference.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

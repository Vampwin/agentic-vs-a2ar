"""verify_references.py -- does every reference resolve to the paper it claims?

Fabricated or mis-attributed citations are the best-documented failure mode of
language-model-written text, and this manuscript's 61 references have never been
checked. The project audits every *number* against an artifact (G8, 207 checks)
and has never audited a single *citation* against its source.

For each reference carrying a DOI, this resolves it against the Crossref API and
compares the recorded title, first-author surname, year and journal with what
the manuscript states. References without a DOI are reported as
UNVERIFIABLE_NO_DOI with their full text, for a person to check by hand.

**What is sent out:** the DOI only. DOIs are public identifiers whose purpose is
resolution, and no manuscript text, author detail or contact address leaves this
machine. Read-only; nothing is written outside the run directory.

**Judgement is not automated away.** A title comparison is a similarity score,
not a verdict. Anything below `TITLE_OK` is reported as CHECK with both strings
printed, so the difference is visible rather than adjudicated by a threshold.

Usage:
    python verify_references.py --out <a new directory>
"""
from __future__ import annotations

import argparse
import json
import re
import time
from difflib import SequenceMatcher
from pathlib import Path

import requests

TITLE_OK = 0.85            # below this a human looks at both strings
TITLE_CLEARLY_WRONG = 0.45  # below this it is very unlikely to be the same paper
REQUEST_PAUSE_S = 0.35
TIMEOUT_S = 25
USER_AGENT = "AIIM-manuscript-reference-check/1.0 (citation verification)"
CROSSREF = "https://api.crossref.org/works/"
DATACITE = "https://api.datacite.org/dois/"

# References a person has checked by hand against the registry record, with the
# reason the automatic comparison flags them. Same idiom as `evidence_not_used.py`:
# a hand-maintained dict with a required reason, so the judgement is reviewable
# rather than invisible. **A reason here does not change what was retrieved** —
# the registry metadata is still recorded in the output beside the verdict.
ADJUDICATED = {
    12: "The Crossref record for ACM Comput Surv 2018;51(1):4 carries the "
        "truncated title 'Metamorphic Testing'. Author, year and journal all "
        "agree; the registry entry is incomplete, not the citation.",
    18: "10.5281/zenodo.591637 is RDKit's concept DOI, which by design resolves "
        "to whichever release is newest (2026_03_6 today), so its title can "
        "never match a version-independent citation. The two versions this "
        "study used are cited by their own release DOIs alongside it, and both "
        "were confirmed to resolve to the releases they name.",
    63: "A parser artifact, not a discrepancy. The year regex takes the first "
        "four-digit number in the entry, and this title contains one: "
        "'CASF-2016' names the benchmark, not the year of publication. "
        "Crossref's issued date is 2018-11-27, the online date; the article is "
        "in volume 59 issue 2, which is February 2019, and 2019 is what the "
        "publisher and this citation state. Author, journal, volume and pages "
        "all agree.",
}


def _root() -> Path:
    for c in (Path(r"D:\Agentic VS Workflow"), Path("/mnt/d/Agentic VS Workflow")):
        if (c / "00_project_control/08_DECISION_LOG.md").is_file():
            return c
    raise SystemExit("project root not found")


ROOT = _root()
MANUSCRIPT = ROOT / "30_manuscript/CURRENT_v5/manuscript_AIIM_v6_CLEAN.md"


def norm(s: str) -> str:
    s = re.sub(r"[^a-z0-9 ]+", " ", (s or "").lower())
    return re.sub(r"\s+", " ", s).strip()


def similar(a: str, b: str) -> float:
    return SequenceMatcher(None, norm(a), norm(b)).ratio()


def parse_references(text: str) -> list[dict]:
    body = text[text.index("## References"):]
    out = []
    for num, raw in re.findall(r"^\[(\d+)\]\s+(.+?)$", body, flags=re.M):
        doi = None
        m = re.search(r"doi:\s*(10\.\S+?)\s*$", raw, flags=re.I) or \
            re.search(r"(10\.\d{4,9}/[^\s]+)", raw)
        if m:
            doi = m.group(1).rstrip(" .,;")
        year = None
        y = re.findall(r"\b(19|20)\d{2}\b", raw)
        if y:
            year = int(re.findall(r"\b((?:19|20)\d{2})\b", raw)[0])
        # "Surname AB, Surname CD, ... . Title. Journal. Year;..."
        first_author = raw.split(",")[0].split(".")[0].strip()
        surname = first_author.split()[0] if first_author.split() else ""
        parts = [p.strip() for p in raw.split(". ") if p.strip()]
        title = parts[1] if len(parts) > 1 else ""
        out.append({"n": int(num), "raw": raw, "doi": doi, "year_stated": year,
                    "first_author_surname": surname, "title_stated": title})
    return out


def _crossref(doi: str) -> tuple[dict | None, str]:
    try:
        r = requests.get(CROSSREF + requests.utils.quote(doi, safe=""),
                         headers={"User-Agent": USER_AGENT}, timeout=TIMEOUT_S)
    except Exception as exc:
        return None, f"request failed: {type(exc).__name__}"
    if r.status_code == 404:
        return None, "not in Crossref"
    if r.status_code != 200:
        return None, f"Crossref HTTP {r.status_code}"
    try:
        m = r.json()["message"]
    except Exception:
        return None, "unparseable Crossref response"
    dp = (m.get("issued", {}).get("date-parts") or [[None]])[0]
    authors = m.get("author") or []
    return {"title": (m.get("title") or [""])[0], "year": dp[0] if dp else None,
            "first_author": (authors[0].get("family") if authors else None),
            "journal": (m.get("container-title") or [""])[0],
            "type": m.get("type"), "registry": "Crossref"}, "ok"


def _datacite(doi: str) -> tuple[dict | None, str]:
    """Zenodo, Dryad and other data DOIs are registered with DataCite, not
    Crossref. The first version of this script queried Crossref only and
    reported two live Zenodo DOIs as unresolved -- a fault in the checker, not
    in the references."""
    try:
        r = requests.get(DATACITE + requests.utils.quote(doi, safe=""),
                         headers={"User-Agent": USER_AGENT}, timeout=TIMEOUT_S)
    except Exception as exc:
        return None, f"request failed: {type(exc).__name__}"
    if r.status_code != 200:
        return None, f"DataCite HTTP {r.status_code}"
    try:
        a = r.json()["data"]["attributes"]
    except Exception:
        return None, "unparseable DataCite response"
    titles = a.get("titles") or [{}]
    creators = a.get("creators") or []
    first = creators[0] if creators else {}
    return {"title": titles[0].get("title", ""),
            "year": a.get("publicationYear"),
            "first_author": (first.get("familyName") or first.get("name")),
            "journal": (a.get("publisher") if isinstance(a.get("publisher"), str)
                        else (a.get("publisher") or {}).get("name", "")),
            "type": (a.get("types") or {}).get("resourceTypeGeneral"),
            "registry": "DataCite"}, "ok"


def lookup(doi: str) -> tuple[dict | None, str]:
    meta, note = _crossref(doi)
    if meta is not None:
        return meta, note
    time.sleep(REQUEST_PAUSE_S)
    meta2, note2 = _datacite(doi)
    if meta2 is not None:
        return meta2, note2
    return None, f"{note}; {note2}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out = Path(args.out).resolve()
    if out.exists() and any(out.iterdir()):
        raise SystemExit(f"REFUSED: {out} exists and is not empty "
                         "(ADR-AIIM-092: a run directory is written once).")
    out.mkdir(parents=True, exist_ok=True)

    refs = parse_references(MANUSCRIPT.read_text(encoding="utf-8"))
    print(f"{len(refs)} references parsed; "
          f"{sum(1 for r in refs if r['doi'])} carry a DOI\n")

    results, counts = [], {}
    for ref in refs:
        if not ref["doi"]:
            v = {"status": "UNVERIFIABLE_NO_DOI"}
        else:
            time.sleep(REQUEST_PAUSE_S)
            meta, note = lookup(ref["doi"])
            if meta is None:
                v = {"status": "UNRESOLVED", "note": note}
            else:
                # The manuscript style drops subtitles ("Merck molecular force
                # field" for "...field. I. Basis, form, scope..."), so a stated
                # title that is a prefix of the registered one is a truncation,
                # not a mismatch. Scored both ways; the better is used.
                ts = similar(ref["title_stated"], meta["title"])
                nm, nc = norm(ref["title_stated"]), norm(meta["title"])
                is_prefix = bool(nm) and nc.startswith(nm)
                if is_prefix:
                    ts = max(ts, 1.0)
                author_ok = (meta["first_author"] or "").lower().startswith(
                    ref["first_author_surname"].lower()[:4] or "\0")
                year_ok = (ref["year_stated"] is None or meta["year"] is None
                           or abs(ref["year_stated"] - meta["year"]) <= 1)
                if ts >= TITLE_OK and author_ok and year_ok:
                    status = "PASS"
                elif ts < TITLE_CLEARLY_WRONG:
                    status = "MISMATCH"
                else:
                    status = "CHECK"
                if status in ("CHECK", "MISMATCH") and ref["n"] in ADJUDICATED:
                    status = "PASS_ADJUDICATED"
                v = {"status": status, "title_similarity": round(ts, 3),
                     "adjudication": ADJUDICATED.get(ref["n"]),
                     "stated_title_is_prefix_of_registered": is_prefix,
                     "first_author_matches": bool(author_ok),
                     "year_matches": bool(year_ok), "crossref": meta}
        results.append({**ref, **v})
        counts[v["status"]] = counts.get(v["status"], 0) + 1
        mark = {"PASS": "ok", "PASS_ADJUDICATED": "ok (checked)",
                "CHECK": "CHECK", "MISMATCH": "MISMATCH",
                "UNRESOLVED": "UNRESOLVED",
                "UNVERIFIABLE_NO_DOI": "no DOI"}[v["status"]]
        print(f"  [{ref['n']:>2}] {mark:<12} {ref['title_stated'][:62]}")
        if v["status"] == "PASS_ADJUDICATED":
            print(f"        adjudicated: {ADJUDICATED[ref['n']][:96]}")
        if v["status"] in ("CHECK", "MISMATCH"):
            print(f"        manuscript: {ref['title_stated'][:88]}")
            print(f"        crossref  : {v['crossref']['title'][:88]}")
            print(f"        author {ref['first_author_surname']} vs "
                  f"{v['crossref']['first_author']}; year "
                  f"{ref['year_stated']} vs {v['crossref']['year']}")

    report = {
        "question": "Does every reference resolve to the paper it claims?",
        "manuscript": MANUSCRIPT.relative_to(ROOT).as_posix(),
        "sent_externally": "the DOI only; no manuscript text, author detail or "
                           "contact address left this machine",
        "thresholds": {"title_similarity_pass": TITLE_OK,
                       "title_similarity_clearly_wrong": TITLE_CLEARLY_WRONG,
                       "note": "CHECK means a person should compare the two "
                               "strings, which are both recorded here"},
        "n_references": len(refs), "counts": counts, "references": results,
    }
    (out / "reference_check.json").write_text(json.dumps(report, indent=2),
                                              encoding="utf-8")
    print(f"\n{counts}")
    print(f"written -> {out / 'reference_check.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

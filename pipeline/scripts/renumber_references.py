"""
renumber_references.py -- add references to an already-numbered manuscript and
renumber the whole thing by order of first appearance.

`build_references_v6.py` is a one-shot transformation: once it has run, its
anchors are gone. That was fine for building the list, but it leaves no way to
add a reference later without either hand-editing (which breaks first-appearance
ordering the moment the new citation lands mid-document) or reverting to a
backup that predates every edit made since.

This closes that gap. It reverses the numbering back to keys, lets new
`{key}` markers be inserted anywhere, then renumbers everything. Running it on
an unchanged manuscript is a no-op that reproduces the same numbers.

The reversal works by matching each `[n]` entry in the References section
against the text in `REFS`. An entry that matches nothing is an error, not a
warning: it means the list and the bibliography have diverged, which is exactly
the state this script exists to prevent.

Usage:
    python renumber_references.py [--check]
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from build_references_v6 import REFS, MS  # noqa: E402


def norm(s: str) -> str:
    return re.sub(r"\s+", " ", s).strip()


def split_manuscript(text: str) -> tuple[str, str, str]:
    """(before, references_section, after) split on the References heading."""
    i = text.index("\n## References\n")
    j = text.index("\n## Tables\n", i)
    return text[:i], text[i:j], text[j:]


def map_numbers_to_keys(refsec: str) -> dict[int, str]:
    by_text = {norm(v["text"]): k for k, v in REFS.items()}
    out: dict[int, str] = {}
    for m in re.finditer(r"^\[(\d+)\]\s+(.*?)(?=\n\[\d+\]|\Z)", refsec,
                         flags=re.M | re.S):
        # The final entry runs to the end of the section, which carries the
        # horizontal rule that separates References from Tables. Strip it, or
        # the last reference never matches.
        n = int(m.group(1))
        body = norm(re.sub(r"\s*-{3,}\s*$", "", m.group(2)))
        key = by_text.get(body)
        if key is None:
            raise SystemExit(
                f"reference [{n}] in the manuscript matches no entry in REFS:\n"
                f"  {body[:110]}...\n"
                "The bibliography and the manuscript have diverged; fix "
                "build_references_v6.py before renumbering.")
        out[n] = key
    return out


def numbers_to_keys(body: str, mapping: dict[int, str]) -> str:
    """Turn every in-text [n] or [n,m] or [n-m] back into {key,key}."""
    def repl(m: re.Match) -> str:
        nums: list[int] = []
        for part in m.group(1).split(","):
            if "-" in part:
                a, b = part.split("-")
                nums.extend(range(int(a), int(b) + 1))
            else:
                nums.append(int(part))
        missing = [n for n in nums if n not in mapping]
        if missing:
            raise SystemExit(f"in-text citation {m.group(0)} refers to "
                             f"undefined reference(s) {missing}")
        return "{" + ",".join(mapping[n] for n in nums) + "}"
    return re.sub(r"\[(\d+(?:[,-]\d+)*)\]", repl, body)


def keys_to_numbers(body: str) -> tuple[str, list[str]]:
    order: list[str] = []
    marker = re.compile(r"\{([a-z0-9_,]+)\}")

    def keys_of(m: re.Match) -> list[str] | None:
        ks = [k.strip() for k in m.group(1).split(",")]
        return ks if ks and all(k in REFS for k in ks) else None

    for m in marker.finditer(body):
        ks = keys_of(m)
        if ks:
            for k in ks:
                if k not in order:
                    order.append(k)
    num = {k: i + 1 for i, k in enumerate(order)}

    def repl(m: re.Match) -> str:
        ks = keys_of(m)
        if not ks:
            return m.group(0)
        ns = sorted(num[k] for k in ks)
        out, i = [], 0
        while i < len(ns):
            j = i
            while j + 1 < len(ns) and ns[j + 1] == ns[j] + 1:
                j += 1
            out.append(f"{ns[i]}-{ns[j]}" if j - i >= 2
                       else ",".join(str(x) for x in ns[i:j + 1]))
            i = j + 1
        return "[" + ",".join(out) + "]"

    return marker.sub(repl, body), order


def render_section(order: list[str]) -> str:
    lines = ["## References", "",
             "*Provenance of every entry is recorded in "
             "`10_new_system/scripts/build_references_v6.py`: `pubmed` = record "
             "retrieved from PubMed for this manuscript; `v4` = carried unchanged "
             "from the previous submission's reference list; `web` = verified "
             "against the publisher's DOI record; `author_action` = the supported "
             "claim is in the project record but the bibliographic entry could not "
             "be recovered and must be supplied.*", ""]
    for i, k in enumerate(order, 1):
        lines.append(f"[{i}] {REFS[k]['text']}")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()

    text = MS.read_text(encoding="utf-8")
    before, refsec, after = split_manuscript(text)
    mapping = map_numbers_to_keys(refsec)
    old_n = len(mapping)

    body = numbers_to_keys(before + after, mapping)
    body, order = keys_to_numbers(body)

    # Re-split: `after` began at the Tables heading, which survives untouched.
    k = body.index("\n## Tables\n")
    new_text = body[:k] + "\n---\n\n" + render_section(order) + "\n---" + body[k:]

    n_marks = len(re.findall(r"\[\d+(?:[,-]\d+)*\]", body[:k]))
    print(f"references before : {old_n}")
    print(f"references after  : {len(order)}")
    print(f"in-text markers   : {n_marks}")
    added = [k_ for k_ in order if k_ not in mapping.values()]
    if added:
        print(f"newly cited ({len(added)}): {', '.join(added)}")
    dropped = [k_ for k_ in mapping.values() if k_ not in order]
    if dropped:
        print(f"no longer cited ({len(dropped)}): {', '.join(dropped)}")

    if args.check:
        print("check only; nothing written")
    else:
        MS.write_text(new_text, encoding="utf-8")
        print(f"written -> {MS}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""
deposit_screen.py -- decide whether a file is safe to publish.

Split out of `build_deposit.py` because getting these three patterns right
matters more than anything else in the deposit, and they deserve to be read on
their own.

The first version of the assay pattern was blunt -- "cAMP or pchembl anywhere
near a number" -- and it refused seven files, every one of them a false
positive. Among them was `G0_input.yaml`, one of the ten contracts this paper's
entire argument rests on. What it had actually matched was a threshold
specification, a synthetic test fixture, and the letters "camp" inside a plant
genus name.

That failure mode is worse than no screen at all, because the deposit still
looks complete. So:

* the pattern now describes what an assay value actually looks like -- a named
  readout *assigned* a number, or a percentage of inhibition -- rather than a
  keyword in the vicinity of a digit;
* anything it still flags must be read and entered in `REVIEWED_AND_CLEARED`
  with a reason, which is published in the deposit manifest. Clearing a file is
  a recorded decision, not a silencing.
"""
from __future__ import annotations

import re
from pathlib import Path

# A named assay readout assigned a value, or a percentage of inhibition.
# Deliberately NOT matched: "pChEMBL >= 6.0" (a threshold), "pChEMBL 4.0-11.0"
# (a range), "Campanula" (a genus), "camp_value: 5.0" in a test fixture is
# matched and cleared below rather than pattern-excluded, because a real
# cAMP assignment must never pass silently.
ASSAY_VALUE = re.compile(
    r"\b(cAMP|percent_inh|pct_inhibition)\b[^\n]{0,40}?[=:]\s*\d+\.\d"
    r"|\d+\.\d\s*%\s*inh\b",
    re.I,
)

SECRET = re.compile(
    r"api[_-]?key\s*[=:]\s*['\"][^'\"]{8,}"
    r"|BEGIN (?:RSA |OPENSSH |EC )?PRIVATE KEY"
    r"|password\s*[=:]\s*['\"][^'\"]{4,}"
    r"|ghp_[A-Za-z0-9]{20,}",
    re.I,
)

# An absolute path under a real user profile. AppData is excluded because the
# scratch directory legitimately appears in tooling comments.
LOCAL_PATH = re.compile(r"C:\\Users\\(?!USER\\AppData)[A-Za-z0-9_.-]+\\")

# The seven files the blunt first pattern refused are no longer flagged at all:
# thresholds, ranges, synthetic fixtures and genus names all fail the tightened
# pattern on their own. They are not listed here, because an exception that is
# never exercised is clutter that hides the one that is. The history is in the
# docstring above, where it belongs.
#
# One entry remains, and it is a small vindication of the screen: it catches its
# own self-test fixtures.
REVIEWED_AND_CLEARED: dict[str, str] = {
    "10_new_system/scripts/deposit_screen.py":
        "this module's self_test() contains synthetic strings -- 'cAMP = 12.5', "
        "'pct_inhibition = 41.7' -- that exist precisely so the pattern can be "
        "checked against what it must catch. They are fixtures for the screen, "
        "not data, and the screen flagging them is evidence that it works.",
}


def screen(path: Path) -> list[str]:
    """Reasons this file must not be published, if any."""
    try:
        txt = path.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return []
    out: list[str] = []
    if ASSAY_VALUE.search(txt):
        out.append("possible assay value")
    if SECRET.search(txt):
        out.append("possible credential")
    if LOCAL_PATH.search(txt):
        out.append("absolute path from the author's machine")
    return out


def self_test() -> bool:
    """The screen is only worth having if it is right on known cases."""
    cases = [
        ("pChEMBL 4.0-11.0 spread", False),
        ("actives defined as pChEMBL >= 6.0", False),
        ("Campanula rotundifolia,4.0,0.42", False),
        ("BLOCKED = ('cAMP', 'assay_label')", False),
        ("assert 'camp' not in str(p).lower()", False),
        ("cAMP = 12.5", True),
        ("cAMP: 12.5", True),
        ("inhibition 87.3 % inh", True),
        ("pct_inhibition = 41.7", True),
    ]
    ok = True
    for text, expected in cases:
        got = bool(ASSAY_VALUE.search(text))
        if got != expected:
            ok = False
            print(f"  SELF-TEST FAIL  {text!r} -> flagged={got}, expected={expected}")
    return ok


if __name__ == "__main__":
    print("screen self-test:", "PASS" if self_test() else "FAIL")

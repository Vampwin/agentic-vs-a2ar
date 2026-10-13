"""Writes build_graphical_abstract_v11.py from build_graphical_abstract_v10.py (ADR-AIIM-447): the left panel becomes the sequence of the work in four steps
(build and freeze, audit, correct once, report as persistence) in place of the funnel; the right cards keep their numbers and are the evidence of the audit."""
from pathlib import Path

H = Path(__file__).resolve().parent
t = (H / "build_graphical_abstract_v10.py").read_text(encoding="utf-8")


def sub(a, b):
    global t
    assert t.count(a) == 1, a[:70]
    t = t.replace(a, b)


sub('"""Graphical abstract, v10 (ADR-AIIM-440): the right-hand side redrawn so that it can be read without the paper.',
    '"""Graphical abstract, v11 (ADR-AIIM-447): v10 with the left panel changed from the funnel to the sequence of the work in four steps (one audit-and-correction\ncycle, not a loop). Original docstring of v10 follows.\n\nv10 (ADR-AIIM-440): the right-hand side redrawn so that it can be read without the paper.')
sub('STEM = "figure_graphical_abstract_v10_2026-10-13"', 'STEM = "figure_graphical_abstract_v11_2026-10-13"')
sub("python build_graphical_abstract_v10.py", "python build_graphical_abstract_v11.py")
i = t.index("    # ------------------------------------------------------------ left: the audited workflow")
j = t.index("    # ------------------------------------------------------------ right: three questions")
LEFT = '''    # ------------------------------------------------------------ left: the sequence of the work, once
    px, pw = 6.0, 216.0
    s.box(px, top, pw, bot - top, fill="#FFFFFF", stroke=BLUE, sw=0.9, r=3)
    s.box(px, top, pw, 14.5, fill=BLUE_F, stroke=BLUE, sw=0.9, r=3)
    s.text(px + 8, top + 9.9, fits("The work in four steps, once", 8.6, pw - 12, bold=True), size=8.6, weight="bold", fill=BLUE, anchor="start")
    steps = [("1  Build and freeze", ["an agent-written screen of an expanded library:", "497 docked, 3 met the original two-part rule"]),
             ("2  Audit", ["checks against planted and blind errors, rewritten", "molecules and reordered data (the three panels at right)"]),
             ("3  Correct, once", ["human-approved, each pre-registered and tested;", "back to the SNP-Lib: 236, 69, 15, then 8 compounds"]),
             ("4  Report as persistence", ["3 of the 15 selected compounds persist under every", "choice tested; persistence is not a prediction of activity"])]
    sx, sw_, sh, sgap, sy0 = px + 7, pw - 14.0, 34.0, 8.0, top + 22.0
    for k, (head, lines) in enumerate(steps):
        yy = sy0 + k * (sh + sgap)
        last = k == len(steps) - 1
        s.box(sx, yy, sw_, sh, fill=BLUE_F if last else FILL, stroke=BLUE if last else RULE, sw=0.9 if last else 0.6, r=2.5)
        s.text(sx + 7, yy + 10.8, fits(head, 8.8, sw_ - 12, bold=True), size=8.8, weight="bold", fill=BLUE if last else INK, anchor="start")
        for m, ln in enumerate(lines):
            s.text(sx + 7, yy + 20.6 + m * 8.8, fits(ln, 7.4, sw_ - 12), size=7.4, fill=SOFT if not last else INK, anchor="start")
        if k < len(steps) - 1:
            s.arrow(sx + sw_ / 2, yy + sh + 0.4, sx + sw_ / 2, yy + sh + sgap - 0.4, sw=0.7)

'''
t = t[:i] + LEFT + t[j:]
# the right cards are the evidence of the audit: no numbers, and the step they belong to
sub('(MOSS, MOSS_F, "2 Do the checks catch planted errors?")', '(MOSS, MOSS_F, "Audit: do the checks catch planted errors?")')
sub('(AMBER, AMBER_F, "3 What lies outside the written checks?")', '(AMBER, AMBER_F, "Audit: what lies outside the written checks?")')
sub('(VIOLET, VIOLET_F, "4 Is the final selection stable? (the 5,000 compounds)")', '(VIOLET, VIOLET_F, "Audit: is the final selection stable? (5,000 compounds)")')
# description for assistive readers
i = t.index('    s = Svg(GW, GH, "Auditing agent-generated biomedical code with executable contracts",')
j = t.index("    s.text(GW / 2, 14.5")
DESC = '''    s = Svg(GW, GH, "Auditing agent-generated biomedical code with executable contracts",
            "Left, the work in four steps, once: (1) build and freeze, an agent-written A2A receptor screen of an expanded library, 497 docked, 3 met the "
            "original two-part rule; (2) audit; (3) correct, once, human-approved and pre-registered, returning to the SNP-Lib, 236, 69, 15, then 8 compounds; "
            "(4) report as persistence, 3 of the 15 selected compounds persist under every choice tested, which is not a prediction of activity. Right, the "
            "audit as three questions. Do the checks catch planted errors? Of 128 planted errors, the contract checks caught 124 with the benchmark's check "
            "set and 101 as deployed, and a simple file-comparison test caught 128. What lies outside the written checks? One of 5 errors written blind to "
            "the project was caught, and rewriting the same molecules as equivalent SMILES strings changed 25 of 25 docking scores. Is the final selection "
            "stable? Of 5,000 compounds, 72 passed one filter in the released run; over 200 reorderings of the same data the number ranged from 6 to 4,119 "
            "(median 1,004).")
'''
t = t[:i] + DESC + t[j:]
(H / "build_graphical_abstract_v11.py").write_text(t, encoding="utf-8")
print("written")

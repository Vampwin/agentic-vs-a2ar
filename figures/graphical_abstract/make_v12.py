"""Writes build_graphical_abstract_v12.py from build_graphical_abstract_v11.py (ADR-AIIM-454): card 3 shows the two rounds of errors written blind to the project
(1 of 5 caught in the first round, 3 of 6 evaluable cases in the second), as 5 and 6 dots, in place of the first round alone."""
from pathlib import Path

H = Path(__file__).resolve().parent
t = (H / "build_graphical_abstract_v11.py").read_text(encoding="utf-8")


def sub(a, b):
    global t
    assert t.count(a) == 1, a[:70]
    t = t.replace(a, b)


sub('"""Graphical abstract, v11 (ADR-AIIM-447):',
    '"""Graphical abstract, v12 (ADR-AIIM-454): v11 with card 3 showing both rounds of errors written blind to the project (1 of 5, then 3 of 6 evaluable).\nOriginal docstring of v11 follows.\n\nv11 (ADR-AIIM-447):')
sub('STEM = "figure_graphical_abstract_v11_2026-10-13"', 'STEM = "figure_graphical_abstract_v12_2026-10-14"')
sub("python build_graphical_abstract_v11.py", "python build_graphical_abstract_v12.py")
sub('What lies outside the written checks? One of 5 errors written blind to "\n            "the project was caught, and',
    'What lies outside the written checks? Errors written blind to "\n            "the project were caught in 1 of 5 cases in a first round and in 3 of 6 evaluable cases in a second round, and')
old = '''    for i in range(5):
        circle(s, lx + 10 + i * 10.5, r0 - 9.6, 3.5, stroke=col, sw=0.9, fill=col if i == 0 else "#FFFFFF")
    s.text(tx, r0 - 6.6, fits("1 of 5 errors written blind to the project caught", 7.8, lw - 74, bold=True),
           size=7.8, weight="bold", anchor="start")'''
new = '''    for i in range(11):                       # first round: 5 cases, 1 caught; second round: 6 evaluable cases, 3 caught
        gx = lx + 9 + (i + (0.7 if i >= 5 else 0)) * 4.8
        caught = i == 0 or (i >= 5 and i - 5 < 3)
        circle(s, gx, r0 - 9.6, 2.0, stroke=col, sw=0.8, fill=col if caught else "#FFFFFF")
    s.text(tx, r0 - 6.6, fits("errors written blind caught: 1 of 5, then 3 of 6", 7.8, lw - 74, bold=True),
           size=7.8, weight="bold", anchor="start")'''
sub(old, new)
(H / "build_graphical_abstract_v12.py").write_text(t, encoding="utf-8")
print("written")

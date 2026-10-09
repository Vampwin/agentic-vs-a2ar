"""A second, weaker audit for the numbers written into the v8 documents after the freeze.

G8 (10_new_system/scripts/g8_number_traceability.py, unmodified) anchors a fixed list of v6-era numbers to
literal strings in the v6 markdown, so it cannot see a number added since. This script takes every numeric token
inside a TRACKED INSERTION of the two v8 working documents and asks whether some number in the artifacts of the
post-freeze runs, rounded to the token's own number of decimals, equals it. Tier A: the JSON, CSV and MD
artifacts of the run directories. Tier B: only found in the text of decision-log entries ADR-AIIM-311 onward.
Tier C: found in neither, and listed for a person to check.

What it can and cannot show. A match says the number exists somewhere in the artifacts; it does not say the number is
attached to the right quantity, so a transposed pair of correct numbers would pass. A miss says nothing is
recorded under that value. It is a typo and drift detector, not a derivation, and it is not G8.
Small integers (below 100, no decimal point) match almost anything and are counted, not tested.
"""
import csv, json, re, sys, zipfile
from pathlib import Path
ROOT = Path(r"D:\Agentic VS Workflow")
RUNS = ROOT / "20_new_runs"
OUT = RUNS / "g8_extension_2026-10-09"
DIRS = ["diffdock_snplib_2026-10-07", "diffdock_dude_2026-10-05", "md_snplib_duplicate_2026-10-08",
        "md_control_replicates_window50_2026-10-07", "md_window50_2026-10-07", "md_control_spread_2026-10-08",
        "em_reproducibility_2026-10-06", "em_survivors_2026-10-06", "snplib_dock_2026-10-06",
        "snplib_prioritisation_2026-10-05", "snplib_leadlike_2026-10-05", "pki_regression_2026-10-05",
        "chembl251_full_curation_2026-10-05", "figure8_window50_2026-10-07", "md_control_anchor_2026-10-07",
        "md_control_replicates_2026-10-06", "snp_screen_2026-08-04"]
NUM = re.compile(r"(?<![\w.])-?\d+(?:\.\d+)?")


def tokens(text):
    text = re.sub(r"(?<=\d),(?=\d{3}\b)", "", text)
    return NUM.findall(text)


def registry(files):
    reg = {}
    for f in files:
        try:
            if f.suffix == ".json":
                txt = f.read_text("utf-8", "ignore")
            elif f.suffix == ".csv":
                txt = f.read_text("utf-8", "ignore")
            else:
                txt = f.read_text("utf-8", "ignore")
        except Exception:
            continue
        for t in set(tokens(txt)):
            try:
                v = abs(float(t))
            except ValueError:
                continue
            for k in range(0, 7):
                reg.setdefault(k, set()).add(round(v, k))
                if v <= 1.0:
                    reg[k].add(round(v * 100, k))       # a fraction is quoted as a percentage in the text
    return reg


def inserted_text(path):
    x = zipfile.ZipFile(path).read("word/document.xml").decode("utf-8")
    x = re.sub(r"<w:del [^>]*/>", "", x)
    x = re.sub(r"<w:del [^>]*>.*?</w:del>", "", x, flags=re.S)
    out = []
    for blk in re.findall(r"<w:ins [^>/]*>.*?</w:ins>", x, flags=re.S):
        out.append("".join(re.findall(r"<w:t(?: [^>]*)?>([^<]*)</w:t>", blk)))
    # inserted paragraphs and table cells carry their own w:ins wrapper, so the loop above covers them
    return " ¶ ".join(out)


SKIP_BEFORE = re.compile(r"(ADR-AIIM-|Section |Sections |Table S|Tables S|Figure |Figures |S\d|ADR-|ref|\[|n = |e-|20\d\d-)\s*$", re.I)


def main():
    A = registry([f for d in DIRS for f in (RUNS / d).rglob("*") if f.is_file()
                  and f.suffix in (".json", ".csv", ".md") and f.stat().st_size < 30_000_000
                  and "results" not in f.parts[len(RUNS.parts):len(RUNS.parts) + 3][-1:]])
    log = (ROOT / "00_project_control/08_DECISION_LOG.md").read_text("utf-8", "ignore")
    i = log.index("## ADR-AIIM-311 ")
    B = registry_text(log[i:])
    P = ROOT / "30_manuscript/AIIM_submission_package_v8_working_2026-10-03"
    res = {}
    for label, f in (("article", P / "KDSS_A2AR_AIIM_v8_WORKING.docx"), ("supplement", P / "Supplementary_Information_v8_WORKING.docx")):
        txt = inserted_text(f)
        n_small = n_a = n_b = 0; miss = []; logonly = []; logonly = []
        for m in NUM.finditer(re.sub(r"(?<=\d),(?=\d{3}\b)", "", txt)):
            t = m.group(0)
            ctx = txt[max(0, m.start() - 40):m.end() + 25]
            if SKIP_BEFORE.search(txt[max(0, m.start() - 14):m.start()]):
                continue
            if "." not in t and abs(float(t)) < 100:
                n_small += 1; continue
            k = len(t.split(".")[1]) if "." in t else 0
            v = round(abs(float(t)), k)
            if v in A.get(k, ()): n_a += 1
            elif v in B.get(k, ()):
                n_b += 1; logonly.append((t, ctx.replace(chr(10), " ")))
            else: miss.append((t, ctx.replace("\n", " ")))
        res[label] = {"tier_A_artifact": n_a, "tier_B_log_only": n_b, "tier_C_unmatched": len(miss),
                      "small_integers_not_tested": n_small, "unmatched": miss, "log_only": logonly}
        print(f"{label}: A {n_a}  B {n_b}  C {len(miss)}  small integers {n_small}")
    (OUT / "new_number_audit.json").write_text(json.dumps(res, indent=1, ensure_ascii=False), "utf-8")
    for label in res:
        print(f"\n== {label}: unmatched tokens ({len(res[label]['unmatched'])})")
        for t, c in res[label]["unmatched"]:
            print(f"  {t:>10}  ...{c}...")


def registry_text(txt):
    reg = {}
    for t in set(tokens(txt)):
        v = abs(float(t))
        for k in range(0, 7):
            reg.setdefault(k, set()).add(round(v, k))
    return reg


if __name__ == "__main__":
    main()

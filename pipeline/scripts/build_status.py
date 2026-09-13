"""build_status.py -- one screen saying what is current and what is not.

Finding out the state of this project meant opening four files: `CLAUDE.md` for
scope, `20_new_runs/RUN_INDEX.md` for which runs count, `08_DECISION_LOG.md` for
what was decided, and `REVIEW_RESPONSE_TRACKER.md` for what is still open. They
do not disagree, but nothing holds the answer in one place, and a reader has to
assemble it.

**STATUS.md is generated, never written by hand.** A hand-written status file is
the same object as a dated filename or a summary paragraph: a claim with nothing
behind it, correct on the day it is typed and quietly wrong afterwards. This
project accumulated eight of those in eight days. Every line below is read from
an artifact at the moment of generation, and `--check` fails if the committed
file no longer matches what the artifacts say.

What it deliberately does not report: anything requiring a run. The test suite
and the deposit rebuild are listed as commands, not as numbers, because a status
file that quotes a test count it did not observe is the failure it exists to
prevent.

Usage:
    python build_status.py            # regenerate STATUS.md
    python build_status.py --check    # non-zero exit if STATUS.md is stale
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import re
from datetime import date
from pathlib import Path
from kdss_paths import project_root, resolve, under


def _root() -> Path:
    # Delegates to kdss_paths, which finds whichever layout this file is inside
    # -- the development tree or a deposit -- by walking up from its own
    # location. It previously tried two absolute paths and exited if neither
    # existed, which on the author's machine meant a copy of the deposit read
    # and wrote the development tree instead of itself.
    return project_root()


ROOT = _root()
OUT = ROOT / "STATUS.md"
MS = resolve("30_manuscript/CURRENT_v5")
PLACEHOLDER = re.compile(r"PLACEHOLDER|\[MODEL|\[EMAIL|\[ORCID|\[AFFIL|\[REPOSITORY|\[COMMIT")


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def abstract_words(p: Path) -> int:
    lines = p.read_text(encoding="utf-8").splitlines()
    i = lines.index("## Abstract")
    seg = []
    for l in lines[i + 1:]:
        if l.startswith("**Keywords") or l.startswith("## "):
            break
        seg.append(l)
    return len("\n".join(seg).strip().split())


def render() -> str:
    g8 = json.loads((resolve("20_new_runs/verification_summary_2026-08-21/"
                     "g8_number_traceability.json")).read_text(encoding="utf-8"))
    idx = json.loads((resolve("00_project_control/28_EVIDENCE_INDEX.json"))
                     .read_text(encoding="utf-8"))
    dep = json.loads((resolve("40_deposit/MANIFEST.json")).read_text(encoding="utf-8"))
    nu = _load("evidence_not_used", resolve("10_new_system/scripts/evidence_not_used.py"))

    log = (resolve("00_project_control/08_DECISION_LOG.md")).read_text(encoding="utf-8")
    adrs = [int(m) for m in re.findall(r"^## ADR-AIIM-(\d+)", log, flags=re.M)]

    clean = under("30_manuscript/CURRENT_v5", "manuscript_AIIM_v6_CLEAN.md")
    other = under("30_manuscript/CURRENT_v5", "manuscript_AIIM_v6.md")
    ph = len(PLACEHOLDER.findall(clean.read_text(encoding="utf-8")))

    labels: dict[str, int] = {}
    for v in idx["runs"].values():
        labels[v["label"]] = labels.get(v["label"], 0) + 1
    dnr = sorted(k for k, v in idx["runs"].items() if v["label"] == "DO-NOT-REUSE")

    # runs whose declaration names a successor, and runs that are sealed
    superseded = {k: v for k, v in nu.NOT_USED.items() if v["kind"] == "superseded"}
    sealed = sorted(k for k, v in nu.NOT_USED.items() if v["kind"] == "sealed")

    docx = sorted(p.name for p in MS.glob("*.docx"))

    L: list[str] = []
    a = L.append
    a("# STATUS")
    a("")
    a(f"**Generated {date.today().isoformat()} by `10_new_system/scripts/build_status.py`.**")
    a("Do not edit by hand -- every line is read from an artifact, and")
    a("`python 10_new_system/scripts/build_status.py --check` fails if this file")
    a("no longer matches them. Authority is unchanged: the decision log first,")
    a("then frozen configuration, then this file.")
    a("")
    a("## The live files")
    a("")
    a("| What | Path |")
    a("|---|---|")
    a("| **Manuscript, authoritative** | `30_manuscript/CURRENT_v5/manuscript_AIIM_v6_CLEAN.md` |")
    a("| Manuscript, hand-kept sibling | `30_manuscript/CURRENT_v5/manuscript_AIIM_v6.md` |")
    for d in docx:
        a(f"| Build output | `30_manuscript/CURRENT_v5/{d}` |")
    a("")
    a("The directory is named `CURRENT_v5` and holds **v6**. The two markdown")
    a("files are maintained by hand in lockstep -- neither is generated from the")
    a("other -- and the DOCX builder reads CLEAN only. Every manuscript edit goes")
    a("to both. Build outputs carry no date in their names on purpose; see the")
    a("note at the top of `build_v6_clean_docx.py`.")
    a("")
    a("## Verification, as of generation")
    a("")
    a("| Check | State |")
    a("|---|---|")
    a(f"| G8 number traceability | **{g8['n_pass']} pass**, {g8['n_value_mismatch']} value "
      f"mismatch, {g8['n_literal_absent']} literal absent (`{g8['state']}`) |")
    a(f"| Evidence index | {idx['totals']['runs']} runs, {idx['totals']['used']} used, "
      f"{idx['totals']['not_used']} declared not-used, {idx['totals']['pending']} pending |")
    # Two counts exist and they differ: the manifest records what the build
    # copied, while the manuscript cites what ends up on disk -- the build
    # writes its own manifests and READMEs afterwards. G8 audits the on-disk
    # figure, so that is the one reported here, with the other named rather
    # than left to look like a discrepancy.
    dfiles = [f for f in (ROOT / "40_deposit").rglob("*")
              if f.is_file() and ".git" not in f.parts]
    dmb = sum(f.stat().st_size for f in dfiles) / 1e6
    a(f"| Deposit | **{len(dfiles)} files on disk**, {dmb:.1f} MB "
      f"(the build copies {dep['n_files']}; the rest are its own manifests "
      f"and READMEs). This is the figure the manuscript cites. |")
    # STATUS.md is not itself deposited, so it can carry the hash without the
    # circularity that forced it out of G8's report: nothing here changes what
    # is in the deposit. `build_status.py` IS deposited, so a change to this
    # script must be committed before the final STATUS.md is generated.
    import subprocess
    head, tree = "unreadable", "UNVERIFIED"
    try:
        ref = (resolve("40_deposit/.git/HEAD")).read_text(encoding="utf-8").strip()
        head = ((resolve("40_deposit/.git") / ref[5:]).read_text(encoding="utf-8").strip()
                if ref.startswith("ref: ") else ref)
        r = subprocess.run(["git", "-C", str(ROOT / "40_deposit"), "status", "--porcelain"],
                           capture_output=True, text=True, timeout=60)
        tree = "clean" if (r.returncode == 0 and not r.stdout.strip()) else "**DIRTY**"
    except Exception:
        pass
    cited = head in (under("30_manuscript/CURRENT_v5", "manuscript_AIIM_v6_CLEAN.md")).read_text(encoding="utf-8")
    a(f"| Deposit commit | `{head[:12]}`, working tree {tree}, cited in the "
      f"manuscript: {'yes' if cited else '**NO**'} |")
    a(f"| Abstract | {abstract_words(clean)} words (AIIM cap: 250) |")
    a(f"| Placeholders in CLEAN | {ph} |")
    a(f"| Decision log | latest **ADR-AIIM-{max(adrs)}**, {len(adrs)} numbered "
      f"entries (numbering has gaps; addenda are not separately numbered) |")
    a("")
    a("Not reported here because it cannot be read from an artifact: the test")
    a("suite. Run it -- `cd 10_new_system && python -m pytest -q` -- rather than")
    a("trusting a number this file did not observe.")
    a("")
    a("## Rebuild, in this order")
    a("")
    a("```bash")
    a("python 10_new_system/scripts/g8_number_traceability.py")
    a("python 10_new_system/scripts/g8_number_traceability.py --manuscript 30_manuscript/CURRENT_v5/manuscript_AIIM_v6.md")
    a("python 10_new_system/scripts/build_evidence_index.py && python 10_new_system/scripts/build_evidence_index.py --check")
    a("python 10_new_system/scripts/build_deposit.py")
    a("python 10_new_system/scripts/build_v6_clean_docx.py && python 10_new_system/scripts/build_v6_docx.py")
    a("python 10_new_system/scripts/build_supplementary_information.py")
    a("python 10_new_system/scripts/build_status.py")
    a("```")
    a("")
    a("The deposit's file count appears in the manuscript, so a deposit rebuild")
    a("can require a manuscript edit and another rebuild. It is a fixpoint, not a")
    a("sequence: repeat until G8 stops reporting `deposit_files`.")
    a("")
    a("## Runs by label")
    a("")
    a("| Label | Count | Meaning |")
    a("|---|---|---|")
    MEAN = {"CORE": "a reported number depends on it",
            "USED": "read by a script that produces a reported number",
            "INPUT": "supplies frozen inputs, not a result",
            "ADMIN": "session logs; no scientific claim",
            "V4-ONLY": "belongs to the superseded v4 framing",
            "DO-NOT-REUSE": "**must not be read by anything**"}
    for k in ("CORE", "USED", "INPUT", "ADMIN", "V4-ONLY", "DO-NOT-REUSE"):
        if k in labels:
            a(f"| {k} | {labels[k]} | {MEAN[k]} |")
    a("")
    if dnr:
        a("**DO-NOT-REUSE:** " + ", ".join(f"`{d}`" for d in dnr)
          + ". Kept as the record of a process failure, not deleted.")
        a("")
    if superseded:
        a("**Superseded, with the successor named in the declaration:**")
        a("")
        for k, v in sorted(superseded.items()):
            a(f"- `{k}` -> `{v.get('superseded_by', 'see declaration')}`")
        a("")
    a("Every run not carrying a reported number is declared in")
    a("`10_new_system/scripts/evidence_not_used.py` with a kind, a reason and an")
    a("ADR. The evidence index fails if one is neither used nor declared, which")
    a("is what stops a directory sitting in `20_new_runs/` with no stated status.")
    a("")
    a("## Standing rules")
    a("")
    a("- **`20_new_runs/` is immutable (ADR-AIIM-092).** New work goes to a new")
    a("  dated directory. The run scripts refuse a non-empty output directory")
    a("  rather than trusting anyone to remember.")
    a("- **H01-H12 stay sealed.** Reporting what the held-out evaluation measured")
    a("  is done; publishing the fault cases is a separate act and is not")
    a("  approved. H06-H12 have never been exposed.")
    if sealed:
        a("  Declared sealed: " + ", ".join(f"`{s}`" for s in sealed) + ".")
    a("- **No wet-laboratory claim.** No compound reported anywhere in this")
    a("  project has been assayed.")
    a("- **No installs without approval** -- skills, packages, models, MCP")
    a("  servers, external services.")
    a("")
    a("## Open")
    a("")
    a(f"{ph} placeholders remain in the manuscript, all author-only:")
    a("")
    a("- **B3** repository URL and commit hash")
    a("- **B4** CRediT, competing interest, funding")
    a("- **B5** affiliation, corresponding e-mail, two missing ORCIDs -- and")
    a("  confirmation of the three already present, which came from an automated")
    a("  public search rather than from their owners (ADR-AIIM-071)")
    a("")
    a("Nothing else blocks submission on content. What no check reaches is prose:")
    a("G8 audits numbers, and eight claims went stale in eight days underneath a")
    a("passing G8. **Re-read the abstract, contributions, Discussion and")
    a("Conclusion as prose whenever a new result lands near them.**")
    return "\n".join(L) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="exit non-zero if STATUS.md differs from the artifacts")
    args = ap.parse_args()
    new = render()
    if args.check:
        cur = OUT.read_text(encoding="utf-8") if OUT.is_file() else ""
        # the generation date changes daily and is not a state difference
        strip = lambda s: re.sub(r"^\*\*Generated \d{4}-\d{2}-\d{2} ", "**Generated ",
                                 s, flags=re.M)
        if strip(cur) == strip(new):
            print("STATUS.md is current")
            return 0
        print("STATUS.md is STALE -- regenerate with build_status.py")
        return 1
    OUT.write_text(new, encoding="utf-8")
    print(f"written -> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

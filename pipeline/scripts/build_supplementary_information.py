"""
build_supplementary_information.py -- assemble the single readable
Supplementary Information document that accompanies the manuscript.

Why this exists. The supplement bundle (S1-S10) is 178 files of contracts,
evaluators, tests, JSON and CSV -- the right form for reuse, the wrong form
for reading. AIIM requires a concise descriptive caption for every
supplementary file and states that supplementary material is published
exactly as received, unchecked and untypeset. A reviewer therefore needs one
document that says what each section contains and carries the few tables that
are referenced from the main text but too long to sit in it.

What this document is NOT: a second copy of the data. Everything here either
describes a bundled section or reproduces a table the manuscript cites. The
artifacts themselves stay in the bundle and in the repository.

Content that is read from its source rather than retyped -- the incident
taxonomy, the held-out per-case outcome, the environment table -- is pulled in
at build time, so this file cannot drift from the artifacts it summarises.

Usage:
    python build_supplementary_information.py [--out PATH]
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches, Pt


def _root() -> Path:
    marker = "00_project_control/08_DECISION_LOG.md"
    for cand in (Path(r"D:\Agentic VS Workflow"), Path("/mnt/d/Agentic VS Workflow")):
        if (cand / marker).is_file():
            return cand
    raise SystemExit("project root not found")


ROOT = _root()
SUP = ROOT / "30_manuscript/supplement"
S10 = SUP / "S10_verification_2026-08-21"
DEFAULT_OUT = ROOT / "30_manuscript/CURRENT_v5/Supplementary_Information_2026-08-24.docx"

# (label, folder, caption). The caption is what AIIM asks for: concise,
# descriptive, one per supplementary file.
SECTIONS = [
    ("S1", "S1_contracts_and_evaluators",
     "The ten frozen G0-G8 audit contracts as YAML, the fifteen deterministic "
     "evaluator modules that enforce them, and the full unit-test suite "
     "(675 tests across 30 modules). This is the governance layer itself, in "
     "the form required to re-run it."),
    ("S2", "S2_fault_definitions",
     "The eighteen development fault families (F01-F18) with full injector "
     "source, one to three per gate. The held-out fault mechanism is "
     "deliberately not included; see the note below."),
    ("S3", "S3_per_case_benchmark_results",
     "Per-case outcomes behind the aggregate detection figures: the "
     "development full-factorial run, the real-data condition B/D rerun, and "
     "four clean/faulted pairs in which defects were injected into real "
     "project artifacts rather than synthetic fixtures."),
    ("S4", "S4_prompts_and_versions",
     "The authoring brief given to the project-blind session that wrote the "
     "held-out cases, and the resolved version of every external tool used -- "
     "RDKit, AutoDock Vina, DiffDock, PyTorch/CUDA, and the molecular-dynamics "
     "toolchain (GROMACS, AmberTools, ff19SB, Lipid21, OPC, GAFF2, "
     "PACKMOL-Memgen), each read from the run's own logs."),
    ("S5", "S5_human_intervention_log",
     "The complete project decision log: every architecture decision record, "
     "its date, its approver, and the reasoning at the time. This is the "
     "authoritative human-intervention record and includes the decisions that "
     "were later reversed."),
    ("S6", "S6_public_benchmark_curation",
     "The ChEMBL A2AR validation set, the raw training pull it derives from, "
     "and the provenance manifest recording how the benchmark was curated."),
    ("S7", "S7_snplib_ledger",
     "The 236-compound SNP-Lib accountability ledger and its metadata: one row "
     "per source compound, with what happened to it at every stage."),
    ("S8", "S8_evidence_chain_refused_ranking",
     "The full structure-based and ligand-based evidence chain behind the "
     "original decision to refuse a ranking on the 236-compound library, "
     "including the water-dependence and cross-docking work and the MD/MM-GBSA "
     "audit result."),
    ("S9", "S9_environment_and_artifact_manifest",
     "The four per-run artifact registries that exist. Four of the runs "
     "reported here used the registry mechanism; the rest recorded results as "
     "direct artifacts, which is stated rather than glossed."),
    ("S10", "S10_verification_2026-08-21",
     "The verification evidence Sections 2 and 3 rest on: seeded-defect "
     "results, the four metamorphic invariants with their inputs and outputs, "
     "the G4 and G6 contract outcomes, the eight-space applicability sweep "
     "behind Figure 3, the funnel-exclusion control, y-randomisation, the "
     "pose-rescue results behind Table 1, the released prioritisation, the "
     "incident taxonomy reproduced below, the pre-registered orthogonal-engine "
     "check, the held-out per-case outcome, and the control simulation's RMSD "
     "and bilayer traces. Every file carries its SHA-256 in MANIFEST.json."),
]


def add_runs(p, text: str) -> None:
    """Inline markdown -> runs: **bold**, *italic*, `code`."""
    for m in re.finditer(r"(\*\*.+?\*\*|\*[^*\n]+?\*|`[^`\n]+?`|[^*`]+)", text, re.S):
        s = m.group(0)
        if s.startswith("**"):
            p.add_run(s[2:-2]).bold = True
        elif s.startswith("`"):
            r = p.add_run(s[1:-1])
            r.font.name = "Courier New"
            r.font.size = Pt(9)
        elif s.startswith("*") and len(s) > 2:
            p.add_run(s[1:-1]).italic = True
        else:
            p.add_run(s)


def md_table(doc: Document, lines: list[str], font_pt: float = 8.5) -> None:
    rows = []
    for l in lines:
        raw = l.strip().strip("|")
        if re.match(r"^[\s:|-]+$", raw):
            continue
        rows.append([c.strip() for c in raw.split("|")])
    if not rows:
        return
    ncol = max(len(r) for r in rows)
    t = doc.add_table(rows=len(rows), cols=ncol)
    t.style = "Table Grid"
    for ri, row in enumerate(rows):
        for ci in range(ncol):
            cell = t.cell(ri, ci)
            para = cell.paragraphs[0]
            add_runs(para, row[ci] if ci < len(row) else "")
            for run in para.runs:
                run.font.size = Pt(font_pt)
                if ri == 0:
                    run.bold = True


def extract_table(md: str) -> list[str]:
    """The first pipe table in a markdown document."""
    out, seen = [], False
    for line in md.split("\n"):
        if line.lstrip().startswith("|"):
            out.append(line)
            seen = True
        elif seen:
            break
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    args = ap.parse_args()

    doc = Document()
    n = doc.styles["Normal"]
    n.font.name = "Times New Roman"
    n.font.size = Pt(11)
    sec = doc.sections[0]
    sec.page_width, sec.page_height = Inches(8.5), Inches(11)
    for m in ("top", "bottom", "left", "right"):
        setattr(sec, f"{m}_margin", Inches(1))

    h = doc.add_paragraph(style="Heading 1")
    h.add_run("Supplementary Information")
    p = doc.add_paragraph()
    add_runs(p, "*Agentic virtual screening against the adenosine A2A receptor: "
                "verification of the agent-generated code as a reportable result*")
    p.alignment = WD_ALIGN_PARAGRAPH.LEFT

    doc.add_paragraph(
        "This document describes the ten supplementary sections accompanying "
        "the article and reproduces three tables that the main text cites but "
        "is too short to carry. The sections themselves are supplied as a "
        "separate archive; the complete code and verification data are also in "
        "the public repository named in the Data availability statement.")

    # ---- section captions -------------------------------------------------
    doc.add_paragraph(style="Heading 1").add_run("1. Contents of the supplementary sections")
    rows = ["| Section | Contents |", "|---|---|"]
    for label, folder, caption in SECTIONS:
        nfiles = sum(1 for f in (SUP / folder).rglob("*") if f.is_file())
        rows.append(f"| **{label}** ({nfiles} files) | {caption} |")
    md_table(doc, rows, font_pt=8.5)

    doc.add_paragraph()
    p = doc.add_paragraph()
    add_runs(p, "**A note on what is withheld.** The held-out fault "
                "*mechanism* — the `apply_fault()` definitions for cases "
                "H01–H05 — is not supplied. What the held-out evaluation "
                "measured is reported in Section 3.1 of the article and its "
                "per-case outcome is reproduced in Section 3 below. The "
                "definitions are withheld so that the five cases retain "
                "whatever residual value they have as a blind test of a system "
                "that could not have read them. Publishing an outcome and "
                "publishing the test that produced it are separate decisions.")

    # ---- held-out per-case ------------------------------------------------
    doc.add_page_break()
    doc.add_paragraph(style="Heading 1").add_run(
        "2. Held-out fault cases: per-case outcome")
    ho = json.loads((S10 / "heldout_outcome/heldout_evaluation_summary.json")
                    .read_text(encoding="utf-8"))
    doc.add_paragraph(
        f"Five fault cases authored by a session blind to this project, "
        f"evaluated once against the unmodified evaluators. "
        f"{ho['n_detected']} of {ho['n_cases']} were detected. Section 3.1 of "
        f"the article gives the root cause of each miss, read from evaluator "
        f"source. The fifth case is reported as inconclusive rather than "
        f"folded into either count: the reporting schema its fault acts on "
        f"does not yet exist in production.")
    rows = ["| Case | Gate | Check | Expected | Actual | Detected |", "|---|---|---|---|---|---|"]
    for r in ho["results"]:
        exp = r["expected"].split("(")[0].strip()
        rows.append(f"| {r['case']} | {r['gate']} | `{r['check']}` | {exp} | "
                    f"{r['actual_state']} | {'Yes' if r['detected'] else 'No'} |")
    md_table(doc, rows, font_pt=8.5)

    fx = json.loads((S10 / "heldout_outcome/fix_verification_summary.json")
                    .read_text(encoding="utf-8"))
    doc.add_paragraph()
    p = doc.add_paragraph()
    add_runs(p, "**The post-repair result is not a held-out number.** The three "
                "implementation gaps were subsequently repaired and unit-tested, "
                "and a regression check confirms the patches fire: "
                + ", ".join(f"{r['case']} → {r['actual_state']}" for r in fx["results"])
                + ". Re-testing cases after seeing them is regression testing of "
                  "a bug fix, not a second measurement of generalisation, and it "
                  "is not reported as one.")

    # ---- incident taxonomy ------------------------------------------------
    doc.add_page_break()
    doc.add_paragraph(style="Heading 1").add_run(
        "3. Incidents detected and corrected during operation")
    tax = (S10 / "incident_taxonomy/incident_taxonomy.md").read_text(encoding="utf-8")
    doc.add_paragraph(
        "Referenced from Section 3.6. Fifteen incidents recorded in the "
        "decision log during construction and execution, classified by the "
        "layer each occurred in, what detected it, and how it resolved: three "
        "scientific, six infrastructure, six agent-reasoning.")
    md_table(doc, extract_table(tax), font_pt=8.0)

    doc.add_paragraph()
    p = doc.add_paragraph()
    add_runs(p, "**Two patterns hold across the set.** Every agent-reasoning "
                "error was caught by checking a primary source the agent had "
                "not consulted — an event log, a tool's own output, a journal's "
                "back catalogue, the project's frozen labels — and none by the "
                "agent re-reading its own reasoning. And two of the scientific "
                "incidents were caught only by checks that exist because the "
                "workflow requires artifacts to be verified before use; both "
                "would otherwise have silently degraded a downstream result "
                "rather than failing visibly.")

    # ---- environment ------------------------------------------------------
    doc.add_page_break()
    doc.add_paragraph(style="Heading 1").add_run("4. Software environment")
    env = (SUP / "S4_prompts_and_versions/environment_versions.md").read_text(encoding="utf-8")
    doc.add_paragraph(
        "Reproduced from S4. Versions are read from each run's own logs or from "
        "the installed binary, not from memory; the force-field entries are the "
        "lines tleap echoed while building the system that was actually "
        "simulated. The full table, including the per-component notes, is in S4.")
    # Two of these tables are both "Component | Version" and read as duplicates
    # side by side, though one is the Windows analysis stack and the other the
    # WSL docking stack. Label each from its position in the source document.
    LABELS = [
        "RDKit, by environment. Two installs are in play and they differ; the "
        "G3 conformer-generation step shells out to WSL rather than using the "
        "Windows install.",
        "Windows-side analysis stack: evaluators, standardisation, QSAR "
        "training and the test suite.",
        "WSL docking stack.",
        "Molecular-dynamics toolchain (Section 2.9).",
    ]
    tables = re.findall(r"(?:^\|.*\n)+", env, flags=re.M)
    for i, tb in enumerate(tables):
        if i < len(LABELS):
            lp = doc.add_paragraph()
            add_runs(lp, f"**{LABELS[i]}**")
        md_table(doc, tb.strip().split("\n"), font_pt=8.0)
        doc.add_paragraph()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    doc.save(str(out))
    print(f"sections described : {len(SECTIONS)}")
    print(f"tables reproduced  : {len(doc.tables)}")
    print(f"written -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

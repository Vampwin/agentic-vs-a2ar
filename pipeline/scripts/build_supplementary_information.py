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
import csv
import json
import re
import statistics
from pathlib import Path
from kdss_paths import project_root, resolve, under

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches, Pt


def _root() -> Path:
    # Delegates to kdss_paths, which finds whichever layout this file is inside
    # -- the development tree or a deposit -- by walking up from its own
    # location. It previously tried two absolute paths and exited if neither
    # existed, which on the author's machine meant a copy of the deposit read
    # and wrote the development tree instead of itself.
    return project_root()


ROOT = _root()
SUP = resolve("30_manuscript/supplement")
S10 = under("30_manuscript/supplement", "S10_verification_2026-08-21")
# No date in the name, deliberately. A date in a filename is a claim that
# nothing verifies: this file was called `..._2026-08-24.docx` while its
# content was rebuilt daily for twelve days, and the name was wrong for
# eleven of them. The build date is in the file's own mtime, the version
# is in the decision log, and both of those are checkable. An undated name
# also means the artifact has one path, so a stale copy cannot sit beside
# the current one wearing a different date.
DEFAULT_OUT = resolve("30_manuscript/CURRENT_v5/Supplementary_Information.docx")

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
     "direct artifacts."),
    ("S10", "S10_verification_2026-08-21",
     "The verification evidence Sections 2 and 3 rest on: seeded-defect "
     "results, the four metamorphic invariants with their inputs and outputs, "
     "the G4 and G6 contract outcomes, the eight-space applicability sweep "
     "behind Figure 4, the funnel-exclusion control, y-randomisation, the "
     "pose-rescue results behind Table 5, the released prioritisation, the "
     "incident taxonomy reproduced below, the pre-registered orthogonal-engine "
     "check, the held-out per-case outcome, and the control simulation's RMSD "
     "and bilayer traces. 99 of its 101 files carry a SHA-256 in "
     "MANIFEST.json; the manifest hashes neither itself nor this "
     "section's README."),
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
        "the article and reproduces the tables the main text cites or rests on "
        "but is too short to carry. The sections themselves are supplied as a "
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
        "scientific, six infrastructure, six agent-reasoning. A "
        "sixteenth, found during revision rather than construction, is listed "
        "separately: a re-export that dropped a field and made the reported "
        "detection rate too low for three weeks (ADR-AIIM-131).")
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
    env = (under("30_manuscript/supplement", "S4_prompts_and_versions/environment_versions.md")).read_text(encoding="utf-8")
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

    # ---- 5. per-family detection ------------------------------------------
    # The article gives detection per gate. Per family is where the two
    # incomplete gates resolve into the two families that cause them, and a
    # reviewer cannot see it anywhere else.
    doc.add_page_break()
    doc.add_paragraph(style="Heading 1").add_run(
        "5. Seeded-defect detection, per fault family")
    # The recounted file, not its predecessor. Both sit in the S10 bundle and
    # they disagree: the older one counts five vacuous F13 injections as
    # escapes, giving 133 cases at 93.2 %, which is what this document reported
    # until 2026-09-12 while the article reported 128 at 96.9 %. ADR-AIIM-131
    # is the authority; `as_previously_reported` inside the file preserves the
    # superseded numbers.
    mut = json.loads((S10 / "seeded_defects/mutation_score_by_gate_recounted.json")
                     .read_text(encoding="utf-8"))
    ov, cl = mut["overall"], mut["coverage_limitation"]
    vac = mut["vacuous_excluded"]
    doc.add_paragraph(
        f"Reproduced from S10. {ov['n_fault_cases_counted']} fault cases across "
        f"the {len(cl['stub_families']) and 11} implemented families, "
        f"{ov['n_escaped']} escaped, detection {ov['detection_rate']:.1%} at a "
        f"false-escalation rate of {ov['false_escalation_rate_clean_runs']} on "
        f"clean runs.")
    p_vac = doc.add_paragraph()
    add_runs(p_vac,
             f"**The denominator, and only the denominator.** "
             f"{ov['n_fault_cases_injected']} cases were injected and "
             f"{ov['n_fault_cases_counted']} are counted. In "
             f"{vac['n']} {', '.join(vac['all_in_family'])} instances the "
             f"injector matched nothing and the workspace was left "
             f"byte-unchanged, so there was no fault to escape; an export that "
             f"had dropped the vacuity flag had recorded them as escapes "
             f"({mut['authority']}). 124 cases are detected before and after. "
             f"The superseded figures are preserved inside the artifact under "
             f"`as_previously_reported`. Condition C, the learned monitor "
             f"alone, is given for the whole benchmark in the article's "
             f"ablation table; the recount did not re-derive it per family, so "
             f"it is not tabulated here.")
    rows = ["| Family | Gate | Fault mechanism | Cases | Escaped | Detection | "
            "Vacuous, excluded |", "|---|---|---|---|---|---|---|"]
    for fam, v in sorted(mut["by_family"].items()):
        rows.append(f"| {fam} | {v['gate']} | {v['description']} | "
                    f"{v['n_fault_cases']} | {v['n_escaped_condition_A']} | "
                    f"{v['detection_rate']:.1%} | "
                    f"{v['n_vacuous_excluded']} |")
    md_table(doc, rows, font_pt=8.0)
    doc.add_paragraph()
    p = doc.add_paragraph()
    add_runs(p, "**Coverage limitation.** " + cl["note"] + " Gates "
             + ", ".join(cl["gates_with_no_implemented_fault_family"])
             + " carry no implemented family; the stubs awaiting them are "
             + ", ".join(cl["stub_families"]) + ".")

    # ---- 6. the metamorphic invariants ------------------------------------
    doc.add_page_break()
    doc.add_paragraph(style="Heading 1").add_run(
        "6. Metamorphic invariants: inputs, outputs and verdicts")
    meta = json.loads((S10 / "metamorphic/metamorphic_results.json")
                      .read_text(encoding="utf-8"))
    doc.add_paragraph(
        f"Each invariant states a relation that must hold between two runs whose "
        f"inputs differ in a way that cannot change the answer, so no oracle is "
        f"needed. Both variants below are compared against the same 25-compound "
        f"baseline, already shown identical to production and to the 20-thread "
        f"rerun, so run-to-run noise is zero by measurement rather than by "
        f"assumption. Tolerance {meta['tolerance_kcal_mol']} kcal/mol. The "
        f"ledgers and score files are in S10/metamorphic with their SHA-256.")
    rows = ["| Variant | Invariant | Compared | Identical | max |Δ| | mean |Δ| | Verdict |",
            "|---|---|---|---|---|---|---|"]
    for k, v in meta["variants"].items():
        rows.append(f"| {k} | {v['invariant']} | {v['n_compared']} | "
                    f"{v['n_identical']}/{v['n_compared']} | "
                    f"{v['max_abs_deviation']:.4f} | {v['mean_abs_deviation']:.4f} | "
                    f"**{v['verdict']}** |")
    md_table(doc, rows, font_pt=8.0)

    # ---- 7. applicability domain, eight spaces ----------------------------
    doc.add_page_break()
    doc.add_paragraph(style="Heading 1").add_run(
        "7. Applicability domain across eight descriptor spaces")
    # The expanded selection's table, named by its run: depositing the
    # 236-compound library's own sweep beside it made both basenames collide,
    # and S10 prefixes every member of a colliding set rather than only the
    # later one, precisely so a reader cannot mistake which is which.
    cov = list(csv.DictReader(
        (S10 / "applicability_domain"
              / "ad_sweep_expanded_2026-08-21_ad_coverage_by_descriptor.csv")
        .read_text(encoding="utf-8").splitlines()))
    pct = list(csv.DictReader(
        (S10 / "applicability_domain/ad_percentile_position.csv")
        .read_text(encoding="utf-8").splitlines()))
    doc.add_paragraph(
        "One rule in eight representations: a compound is in-domain if its mean "
        "similarity to its five nearest ChEMBL A2AR training neighbours reaches "
        "the tenth percentile of that training set's own self-similarity. The "
        "last column is the same quantity as a rank, which is what makes the "
        "spaces comparable; the admission rate is that rank thresholded at 10.")
    rows = ["| Space | Threshold | Library mean similarity | In-domain | "
            "Median percentile |", "|---|---|---|---|---|"]
    for r in cov:
        s = r["descriptor_space"]
        col = f"{s}_pctile"
        vals = sorted(float(x[col]) for x in pct) if pct and col in pct[0] else []
        med = f"{statistics.median(vals):.2f}" if vals else "--"
        rows.append(f"| {s} | {float(r['threshold']):.4f} | "
                    f"{float(r['snplib_mean_knn_sim']):.4f} | "
                    f"{r['n_inside']}/{r['n_total']} "
                    f"({100 * float(r['coverage_fraction']):.1f} %) | {med} |")
    md_table(doc, rows, font_pt=8.0)

    # ---- 8. conformal diagnostics -----------------------------------------
    doc.add_page_break()
    doc.add_paragraph(style="Heading 1").add_run(
        "8. Conformal prediction: behaviour on and off domain")
    cd = json.loads((S10 / "released_result/conformal_diagnostics.json")
                    .read_text(encoding="utf-8"))
    cs = cd["calibration_set"]
    doc.add_paragraph(
        f"The empty set is how a Mondrian conformal predictor says an object is "
        f"atypical of both classes. The same mechanism was run on the "
        f"calibration data it was built for and on the two natural-product query "
        f"libraries. The per-class guarantee rests on {cs['n_inactive']} inactive "
        f"and {cs['n_active']} active calibration compounds, which is the "
        f"asymmetry the article's Section 3.6 traces the behaviour to.")
    rows = ["| Panel | α | empty {} | {inactive} | {active} | {both} | "
            "Mean set size |", "|---|---|---|---|---|---|---|"]
    for name, panel in cd["panels"].items():
        for sig in ("0.05", "0.1", "0.2"):
            b = panel["by_significance"][sig]
            rows.append(f"| {name.replace('_', ' ')} | {sig} | "
                        f"{100 * b['empty_{}']:.1f} % | "
                        f"{100 * b['singleton_inactive_{0}']:.1f} % | "
                        f"{100 * b['singleton_active_{1}']:.1f} % | "
                        f"{100 * b['both_{0,1}']:.1f} % | {b['mean_set_size']:.3f} |")
    md_table(doc, rows, font_pt=8.0)
    ta = cd["tier_a_criterion"]
    doc.add_paragraph()
    p = doc.add_paragraph()
    add_runs(p, f"**Where the Tier-A criterion sits.** The rule admits a compound "
                f"at conformal *p*(active) above {ta['alpha']}. No compound of the "
                f"library reaches {ta['library_p_active_max']}; the median is "
                f"{ta['library_p_active_median']}; the granularity is "
                f"{ta['p_active_granularity']}. The {ta['n_flagged']} admitted "
                f"compounds occupy "
                f"{ta['flagged_p_active_range'][0]}–{ta['flagged_p_active_range'][1]}.")
    rows = ["| Compound | Calibrated P(active) | Conformal *p*(active) | "
            "Rank of 5,000 |", "|---|---|---|---|"]
    for c in ta["tier_a"]:
        rows.append(f"| {c['compound_id']} | {c['p_active_calibrated']} | "
                    f"{c['conformal_p_active']} | {c['rank_by_p_active_of_5000']} |")
    md_table(doc, rows, font_pt=8.5)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    doc.save(str(out))
    print(f"sections described : {len(SECTIONS)}")
    print(f"tables reproduced  : {len(doc.tables)}")
    print(f"written -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

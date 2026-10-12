"""
build_deposit.py -- assemble the public code-and-data deposit for the manuscript
as a git-ready tree.

Three rules govern what goes in.

**Nothing that is already archived elsewhere is re-deposited.** LOTUS v11 is
107 MB of gzipped CSV that already carries its own Zenodo DOI; ChEMBL and DUD-E
are public releases. Copying them here would add nothing but weight and a
second, divergeable copy. They are cited with their identifiers and the
SHA-256 this project verified them against, so a reader can confirm they have
the same bytes we did.

**Nothing sealed goes in.** The held-out fault families are excluded. The
manuscript reports the RESULT of evaluating H01-H05 once (one of five detected,
with each miss's root cause), and that summary is deposited; the faulted
workdirs, and anything about H06-H12, which have never been exposed, are not.
CLAUDE.md requires a separate human approval to expose them. The v4 supplement
carried some of that material; this deposit does not inherit that decision.

**The manuscript is not deposited.** Neither is any file that holds its text (the
extension's G8 inputs and the rendered paragraphs are skipped by name).

**Nothing is regenerated.** Every artifact is copied as the analysis read it,
with its SHA-256 recorded, because a deposit that recomputes its own numbers
demonstrates nothing about the run that produced them.

Every file is scanned before copying: assay-derived values, credentials, and
absolute paths from the author's machine all abort the build rather than
producing a deposit that has to be retracted.

Usage:
    python build_deposit.py [--out DIR]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import sys
from pathlib import Path
from kdss_paths import project_root, resolve, under, dev_dir


def _root() -> Path:
    # Delegates to kdss_paths, which finds whichever layout this file is inside
    # -- the development tree or a deposit -- by walking up from its own
    # location. It previously tried two absolute paths and exited if neither
    # existed, which on the author's machine meant a copy of the deposit read
    # and wrote the development tree instead of itself.
    return project_root()


ROOT = _root()
DEFAULT_OUT = dev_dir("40_deposit")

# (source, destination, description)
TREE: list[tuple[str, str, str]] = [
    ("10_new_system/contracts", "pipeline/contracts",
     "the ten frozen G0-G8 audit contracts"),
    ("10_new_system/evaluators", "pipeline/evaluators",
     "the deterministic evaluator modules the agent cannot override"),
    ("10_new_system/fault_injection", "pipeline/fault_injection",
     "the eighteen fault-family definitions (eleven implemented, seven declared stubs)"),
    ("10_new_system/tests", "pipeline/tests",
     "the test suite: 682 collected cases from 567 test functions across "
     "30 modules"),
    ("10_new_system/scripts", "pipeline/scripts",
     "every analysis, figure and audit script referenced in the manuscript"),
    ("10_new_system/kdss_audit", "pipeline/kdss_audit",
     "the artifact registry implementation"),
    # These three were missed in the first build: deposited scripts and tests
    # import them, so without them a clean checkout cannot run its own test
    # suite. `benchmark/heldout_interface` stays excluded by SKIP_DIRS.
    ("10_new_system/risk_monitor", "pipeline/risk_monitor",
     "the learned risk monitor, the component the A/C ablation conditions test"),
    ("10_new_system/telemetry", "pipeline/telemetry",
     "the out-of-domain detector and its telemetry features"),
    ("10_new_system/benchmark", "pipeline/benchmark",
     "SEER scoring and the development-fault benchmark harness"),
    ("30_manuscript/supplement/S10_verification_2026-08-21", "verification",
     "the 2026-08-21 verification evidence behind Sections 2 and 3"),
    ("30_manuscript/figures", "figures",
     "figure sources and their provenance records"),
    # Added 2026-10-10 (ADR-AIIM-377, approved by the author). This folder was kept
    # out of the deposit on purpose until then, because the manuscript cites the
    # deposit's file count; the count is now re-cited at each rebuild. Backups
    # (underscore-prefixed) are skipped as scratch.
    ("30_manuscript/v7_graphical_abstract_2026-10-03", "figures/graphical_abstract",
     "the graphical abstract's builders (v7 and the funnel layout of 2026-10-10) and their outputs"),
    ("20_new_runs/verification_summary_2026-08-21", "verification/summary",
     "per-gate mutation scores, structure and protocol checks, the G8 audit"),
    # NOTE: verification/summary and verification/structure_and_protocol share
    # two files; the duplicate is resolved after copying, below.
    ("20_new_runs/metamorphic_2026-08-21", "verification/metamorphic",
     "the four metamorphic invariants, their inputs and their outputs"),
    # DELIBERATELY NOT DEPOSITED: 20_new_runs/heldout_evaluation_2026-08-08.
    # The manuscript reports the held-out RESULT (one of five, with each miss's
    # root cause, §3.1). The faulted workdirs are a different thing: the
    # corrupted artifacts themselves let H01-H05 be reconstructed, which is
    # unsealing, and CLAUDE.md requires separate explicit human approval for
    # that. The Data Availability statement already offers them on request.
    # Depositing them was attempted on 2026-08-23 and reverted (ADR-AIIM-113).
    ("20_new_runs/diffdock_orthogonal_2026-08-23", "verification/orthogonal_engine",
     "the pre-registered DiffDock cross-check: inputs, all predicted poses, "
     "and the machine-readable result the manuscript's numbers derive from"),
]

# ---------------------------------------------------------------------------
# The post-freeze extension (ADR-AIIM-285 onward): the runs behind Section 3.8
# of the article and Section 17 of the supplement. Added 2026-10-09 (ADR-AIIM-346).
#
# These are deposited under `extension/<run>/` and, unlike the verification
# evidence above, through an ALLOW-LIST of suffixes, because the run directories
# hold simulation inputs and trajectories that are hundreds of megabytes and are
# not what a reader needs: the series, tables, scripts, force-field parameter
# files and the machine-readable results are. NOT deposited: trajectories
# (.xtc), run inputs and checkpoints (.tpr .cpt .edr), structures and topologies
# (.gro .pdb .top), logs, and the 3,685 per-complex DiffDock pose files of the
# DUD-E subsample (the scored table that the analysis read is deposited). The
# raw ChEMBL pull is refused by the screen and is a public release anyway.
EXTENSION_SUFFIX = {".json", ".csv", ".md", ".py", ".sh", ".xvg", ".ndx", ".mdp",
                    ".txt", ".yaml", ".yml", ".png", ".svg", ".mol2", ".frcmod",
                    ".sdf", ".in"}
EXTENSION_MAX_BYTES = 8 * 1024 * 1024
# Regular expressions (on the path inside the run directory) of files skipped for one run. The DiffDock runs
# write up to a hundred poses per complex; the rank-1 pose of each complex is what the analysis and the figures
# read, so the other ranks (about 3,500 files) are not deposited. Rank 1 is `rank1.sdf` and `rank1_confidence-*.sdf`.
EXTENSION_SKIP_PATTERNS: dict[str, list[str]] = {
    "extension/docking_matched_search_space_2026-10-09": [r"(^|/)rank(?!1[._])\d+[^/]*\.sdf$"],
    "extension/docking_matched_filter_pocket_2026-10-09": [r"(^|/)rank(?!1[._])\d+[^/]*\.sdf$"],
}
# (source, destination, description, directory names allowed despite SKIP_DIRS,
#  extra directory names skipped for this entry)
EXTENSION_TREE: list[tuple[str, str, str, set, set]] = [
    (f"20_new_runs/{d}", f"extension/{d}", why, allow, skip)
    for d, why, allow, skip in [
        ("pki_regression_2026-10-05", "random-forest regression on pKi", set(), set()),
        ("pki_censored_em_2026-10-05", "censored-record EM refit", set(), set()),
        ("rank_prioritisation_2026-10-05", "rank prioritisation", set(), set()),
        ("orthogonal_intersection_2026-10-05", "orthogonal intersection", set(), set()),
        ("majority_rule_validation_2026-10-05", "majority-rule validation", set(), set()),
        ("snplib_dedup_2026-10-05", "SNP-Lib de-duplication", set(), set()),
        ("snplib_leadlike_2026-10-05", "SNP-Lib lead-likeness gate", set(), set()),
        ("snplib_prioritisation_2026-10-05", "SNP-Lib prioritisation", set(), set()),
        ("pool_comparison_2026-10-06", "pooled against separate libraries", set(), set()),
        ("em_reproducibility_2026-10-06", "the thread-dependent defect of the censored fit", set(), set()),
        ("em_survivors_2026-10-06", "the fifteen survivors and the eight shortlisted", set(), set()),
        ("split_convention_factorial_2026-10-05", "eight split conventions", set(), set()),
        ("label_balance_domain_2026-10-05", "label balance and domain", set(), set()),
        ("ad_distance_curve_2026-10-05", "applicability-domain distance curve", set(), set()),
        ("ad_methods_distance_performance_2026-10-05", "applicability-domain methods", set(), set()),
        ("arm_k_auc_2026-10-05", "k-nearest-neighbour AUC arm", set(), set()),
        ("snplib_dock_2026-10-06", "Vina docking of the 69 in-domain compounds, with poses", {"poses"}, set()),
        ("diffdock_snplib_2026-10-07", "DiffDock-L on the fifteen survivors, with poses", set(), set()),
        ("diffdock_dude_2026-10-05", "DiffDock-L stage 1 on DUD-E, 335 complexes", set(), {"results"}),
        ("md_control_replicates_2026-10-06", "control replicates 2 and 3, series and scripts", set(), set()),
        ("md_control_anchor_2026-10-07", "the original control, re-analysed", set(), set()),
        ("md_window50_2026-10-07", "published MD figures over 50 ns", set(), set()),
        ("md_control_replicates_window50_2026-10-07", "three-run control over 50 ns", set(), set()),
        ("md_control_spread_2026-10-08", "spread of the three control runs", set(), set()),
        ("md_snplib_duplicate_2026-10-08", "SNP00157 and SNP00163 in duplicate, analysis", set(), {"md_text.json", "make_md_text.py", "apply_md_text.py", "build_section17_helpers.py"}),   # the manuscript paragraphs, rendered and as templates
        ("figure8_window50_2026-10-07", "Figure 7 (formerly 8) at 50 ns", set(), set()),
        ("figures_docking_md_2026-10-09", "Figures 6 and 7 with the survivors and the SNP-Lib runs, scripts and checks", set(), {"Figure_6_docking_snplib.svg"}),
        ("g8_extension_2026-10-09", "G8 on the v8 text and the new-number audit", set(), {"inputs"}),   # inputs/ holds the manuscript text itself
        # Added 2026-10-10 (ADR-AIIM-377): the matched-search-space program (ADR-AIIM-353 to 357) and the figure style
        ("docking_matched_search_space_2026-10-09", "Vina over the whole receptor and DiffDock-L filtered to the site box: pre-registration, rank-1 poses, tables and result", set(), set()),
        ("docking_matched_filter_pocket_2026-10-09", "DiffDock-L 100 samples filtered by the box, and DiffDock-Pocket (main, rigid, C-alpha centre): pre-registration, rank-1 poses, tables and result", set(), set()),
        ("figures_style_harmonised_2026-10-10", "the one style of Figures 3 to 7 and the scripts that redraw Figures 3 to 5", set(), set()),
        # Added 2026-10-11 (ADR-AIIM-396): the third replicates of SNP00157 and SNP00163, the post hoc analyses and the redrawn Supplementary Figure S2
        ("md_snplib_replicate3_2026-10-10", "third replicates of SNP00157 and SNP00163: pre-registration, analysis, post hoc pose families, RMSF, hydrogen bonds, the image-jump check and Supplementary Figure S2", set(), set()),
        ("md_candidates_2026-08-22/SNP00157", "SNP00157 build parameters and series", set(), set()),
        ("md_candidates_2026-08-22/SNP00163", "SNP00163 build parameters and series", set(), set()),
        # Added 2026-10-11 (ADR-AIIM-401, approved by the author): the seeded benchmark through the deployed contracts, whose
        # figures the article now states (Section 3.1; Supplementary Table S6); the six sensitivity analyses of ADR-AIIM-399;
        # Figure 6 and Supplementary Figure S1 with the legend wording of the article; the package lists of the environments
        ("deployed_path_benchmark_2026-09-23", "the seeded benchmark under the harness contract and under the contracts as deployed (five arms): per-case outcomes, check coverage and summary", set(), set()),
        ("v9_sensitivity_2026-10-11", "six sensitivity analyses of 2026-10-11: pre-registration, scripts and results (vote per distinct split, applicability threshold, stratified interval, box size, check map, k and percentile of the applicability rule)", set(), set()),
        ("figures_v9_2026-10-11", "Figure 6 and Supplementary Figure S1 redrawn with the legend wording of the article, scripts and provenance", set(), {"Figure_6_docking_snplib.svg"}),
        ("environment_freeze_2026-10-11", "package lists of the analysis, docking and simulation environments, exported on 2026-10-11", set(), set()),
        # Added 2026-10-12 (ADR-AIIM-414, approved by the author): the censored-target fit under four checks (ADR-AIIM-412, 413)
        ("v9_censored_sensitivity_2026-10-12", "the censored-target fit under four checks of 2026-10-12: pre-registration, scripts and results (simulation, convergence, initialisation, a censored-likelihood model, determinism); the arrays (.npz) are not deposited and are rebuilt by analysis/prep.py", set(), set()),
        # Added 2026-10-12 (ADR-AIIM-418, approved by the author): seeded faults for G2 and G7 (ADR-AIIM-416, 417)
        ("v9_g2_g7_extension_2026-10-12", "seeded faults for gates G2 and G7 of 2026-10-12: pre-registration, injectors, scripts and per-case results under the direct and the deployed path", set(), set()),
        # Added 2026-10-12 (ADR-AIIM-422, approved by the author): seeded faults for G3 and G5 (ADR-AIIM-420, 421)
        ("v9_g3_g5_extension_2026-10-12", "seeded faults for gates G3 and G5 of 2026-10-12: pre-registration, injectors, scripts, the redocked poses of the clean and the shifted-box cases and the per-case results under the direct and the deployed path", set(), set()),
    ]
]

# Directories never copied, wherever they appear.
SKIP_DIRS = {"__pycache__", ".pytest_cache", ".ipynb_checkpoints", "held_out",
             "heldout_interface", "poses", "fastpass_poses", "_superseded",
             "_tools"}
SKIP_SUFFIX = {".pyc", ".pyo", ".log", ".tpr", ".trr", ".xtc", ".cpt", ".edr"}

# Files never copied, by name. The two held-out scripts matter most: excluding
# the sealed DIRECTORY was not enough, because these live in scripts/ alongside
# everything else and would have shipped. Sealing a folder does not seal a
# topic.
SKIP_NAMES = {
    "run_heldout_evaluation_2026-08-08.py",
    "verify_h01_h02_h03_fixes_2026-08-08.py",
}

# Working files created while debugging, prefixed with an underscore by
# convention in this project. They are not deliverables.
#
# Evaluator helper modules are the exception, and getting this wrong was not
# cosmetic. They follow the same underscore convention but ARE deliverables:
# excluding them removed _domain_gate_bias_probe.py -- which the manuscript
# offers as a contract "directly transferable to any screening campaign" -- and
# _assay_blindness.py, which implements a stated safety boundary, from every
# deposit ever built. A reader following the Data availability statement could
# not obtain either one. Found by the internal review of 2026-08-25
# (ADR-AIIM-127).
_KEEP_UNDERSCORED_IN = {"evaluators"}


def _is_scratch(f: Path) -> bool:
    if f.name == "__init__.py":
        return False
    if f.parent.name in _KEEP_UNDERSCORED_IN and f.suffix == ".py":
        return False
    return f.name.startswith("_")

# Screening lives in its own module: getting these patterns right matters
# more than anything else here, and it carries a self-test.
sys.path.insert(0, str(Path(__file__).parent))
from deposit_screen import (  # noqa: E402
    REVIEWED_AND_CLEARED, screen, self_test as screen_self_test,
)

# ---------------------------------------------------------------------------
# The pre-registration entries (added 2026-10-10, ADR-AIIM-377, approved by the author).
# The article and the supplement say that the post-freeze analyses were pre-registered in an append-only
# decision log, with the SHA-256 of the predictions recorded before each ran. A reader cannot check that
# without the entries, and the log as a whole is not deposited, so the entries the text relies on are copied
# verbatim into one file, with the SHA-256 of the whole log at the time of copying. Which entries is decided
# here, in the open: the extension (285-336), the matched search-space program (353-356) and the earlier
# pre-registered verification and simulation entries that the text cites.
PREREGISTRATION_ADRS = sorted({131, 150, 156} | set(range(285, 337)) | {353, 354, 355, 356} | {390, 391} | {399, 400} | {412, 413} | {416, 417} | {420, 421} | {214, 215, 229, 246} | {240, 256, 260})   # 412, 413: the censored-target fit under four checks, pre-registration and result (ADR-AIIM-414, approved by the author); ADR-AIIM-111, 112, 161, 171, 193, 194, 236, 237 and 238 were removed on 2026-10-12 (ADR-AIIM-410, approved by the author): the text calls them prespecified, since no hash of the entry, the predictions or the inputs is recorded for them; 214, 215, 229, 246 and 240, 256, 260: the DUD-E enrichment and the spike-in recovery, pre-registration, amendment and results (ADR-AIIM-406, approved by the author); 390, 391: the third replicates of SNP00157 and SNP00163, pre-registration and result (ADR-AIIM-396, approved by the author); 236 to 238: the seeded benchmark through the deployed contracts; 399, 400: the six sensitivity analyses of 2026-10-11 (ADR-AIIM-401, approved by the author)
PREREGISTRATION_DST = "preregistration/PRE_REGISTRATION_ENTRIES.md"
SEALED_MENTION = re.compile(r"\bH(0[6-9]|1[0-2])\b")   # the held-out families that were never exposed


def write_preregistration(out: Path) -> dict:
    log = ROOT / "00_project_control" / "08_DECISION_LOG.md"
    if not log.exists():
        raise SystemExit("missing source: 00_project_control/08_DECISION_LOG.md")
    text = log.read_bytes().decode("utf-8").replace("\r\n", "\n")
    parts = re.split(r"\n(?=## ADR-AIIM-)", text)
    chosen, found = [], set()
    for p in parts:
        m = re.match(r"## ADR-AIIM-(\d+)", p)
        if m and int(m.group(1)) in PREREGISTRATION_ADRS:
            chosen.append(p.rstrip("\n"))
            found.add(int(m.group(1)))
    missing = [n for n in PREREGISTRATION_ADRS if n not in found]
    if missing:
        raise SystemExit(f"pre-registration entries not found in the decision log: {missing}")
    body = "\n\n".join(chosen)
    if SEALED_MENTION.search(body):
        raise SystemExit("a selected pre-registration entry mentions a sealed held-out family; not depositing it")
    head = (
        "# Pre-registration entries\n\n"
        "Verbatim copies of the decision-log entries that the article and the supplement describe as pre-registered, "
        "and of the entries that record their results. The log is append-only; each entry was written before the run it "
        "concerns (or, for a result entry, after it). Each entry fixed the predictions or the decision rule of its "
        "analysis before the run, and recorded the SHA-256 of the predictions or of the inputs frozen for it where the "
        "run had them. Copied as written: first-person wording and the author's quotations "
        "are unchanged. The full log is not deposited.\n\n"
        f"- source: `00_project_control/08_DECISION_LOG.md`, SHA-256 of the whole log (line ends as LF) "
        f"`{hashlib.sha256(text.encode('utf-8')).hexdigest()}`\n"
        f"- entries: {len(chosen)} (ADR-AIIM-131, 150, 156, 214, 215, 229, 240, 246, 256, 260, 285 to 336, 353 to 356, 390, 391, 399, 400, 412, 413, 416, 417, 420, 421)\n\n---\n\n")
    dst = out / PREREGISTRATION_DST
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(head + body + "\n", encoding="utf-8", newline="\n")
    return {"path": PREREGISTRATION_DST, "source": "00_project_control/08_DECISION_LOG.md",
            "bytes": len(as_stored(dst)), "sha256": sha256(dst)}


EXTERNAL = [
    {"name": "LOTUS natural-product knowledge base",
     "version": "v11 frozen snapshot, 13 April 2026",
     "identifier": "doi:10.5281/zenodo.19360665", "licence": "CC BY 4.0",
     "why_not_deposited": "107 MB, already archived with its own DOI",
     "files": ["lotus_260413_frozen.csv.gz", "lotus_260413_frozen_metadata.csv.gz"]},
    {"name": "ChEMBL A2AR bioactivity records",
     "version": "see pipeline/contracts/G0_input.yaml for the frozen extract",
     "identifier": "doi:10.1093/nar/gky1075", "licence": "CC BY-SA 3.0",
     "why_not_deposited": "public database release; the frozen extract's SHA-256 "
                          "is recorded in the run metadata",
     "files": ["chembl_a2a_training_v3.csv", "chembl_a2ar_validation_set.csv"]},
    {"name": "DUD-E A2AR decoy set",
     "version": "as distributed", "identifier": "doi:10.1021/jm300687e",
     "licence": "see source", "why_not_deposited": "public benchmark release",
     "files": ["dude_a2ar_docking_set.csv"]},
    {"name": "A2AR crystal structures",
     "version": "as deposited", "identifier": "PDB 4EIY, 5OLH",
     "licence": "CC0", "why_not_deposited": "public structural archive",
     "files": ["4EIY.pdb", "5OLH.pdb"]},
]


# Added 2026-10-10 (ADR-AIIM-379, approved by the author). MANIFEST.json used to record the hash of the
# development copy, which has Windows (CRLF) line ends, while the deposit's .gitattributes
# (`* text=auto eol=lf`) makes git store -- and every clone check out -- LF. On a fresh clone 342 of 1,021
# recorded hashes did not reproduce. The manifest now records what a clone contains: for a text file the
# SHA-256 and byte count of its content with CRLF turned into LF, for a binary file those of the bytes as
# they are. "Text" follows git's own test for `text=auto`: no NUL byte in the first 8,000 bytes and no
# carriage return that is not part of CRLF, and not a suffix .gitattributes marks binary.
BINARY_SUFFIX = {".png", ".pdf", ".gz", ".docx"}


def as_stored(p: Path) -> bytes:
    """The bytes a clone of the deposit will contain for this file."""
    raw = p.read_bytes()
    if p.suffix.lower() in BINARY_SUFFIX or b"\0" in raw[:8000]:
        return raw
    if re.search(rb"\r(?!\n)", raw):
        return raw      # a lone CR makes git treat the file as binary, and leave it alone
    return raw.replace(b"\r\n", b"\n")


def sha256(p: Path) -> str:
    return hashlib.sha256(as_stored(p)).hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    args = ap.parse_args()
    if not screen_self_test():
        raise SystemExit('deposit screen failed its own self-test; not building a deposit on a screen that does not work')
    out = Path(args.out)
    # Clear the tree but PRESERVE version-control metadata. The first version
    # of this called shutil.rmtree(out), which deletes .git along with
    # everything else: rebuilding the deposit would have destroyed the commit
    # history that makes it citable. It never fired only because Windows marks
    # git object files read-only and the delete raised instead of succeeding.
    # On a POSIX host it would have succeeded silently.
    # .git carries the citable history; the other four are hand-authored and
    # are NOT regenerated by this script, so a rebuild must not delete them.
    KEEP = {".git", ".gitignore", ".gitattributes", "README.md",
            "CITATION_AND_LICENCE.md", "LICENSE", "CITATION.cff"}
    if out.exists():
        for child in out.iterdir():
            if child.name in KEEP:
                continue
            shutil.rmtree(child) if child.is_dir() else child.unlink()
    out.mkdir(parents=True, exist_ok=True)

    manifest: list[dict] = []
    refused: list[tuple[str, str]] = []
    cleared: list[tuple[str, str, str]] = []
    skipped = 0
    total = 0

    for src_rel, dst_rel, desc, *_ext in [(*t, None, None) for t in TREE] + list(EXTENSION_TREE):
        allow, extra_skip = (_ext[0], _ext[1]) if _ext[0] is not None else (set(), set())
        is_ext = _ext[0] is not None
        src = ROOT / src_rel
        if not src.exists():
            raise SystemExit(f"missing source: {src_rel}")
        for f in sorted(src.rglob("*")):
            if not f.is_file():
                continue
            if any(part in SKIP_DIRS and part not in allow for part in f.parts):
                skipped += 1
                continue
            if f.suffix.lower() in SKIP_SUFFIX:
                skipped += 1
                continue
            if is_ext and (f.suffix.lower() not in EXTENSION_SUFFIX
                           or f.stat().st_size > EXTENSION_MAX_BYTES
                           or any(part in extra_skip for part in f.relative_to(src).parts)
                           or any(re.search(pat, f.relative_to(src).as_posix())
                                  for pat in EXTENSION_SKIP_PATTERNS.get(dst_rel, ()))):
                skipped += 1
                continue
            # underscore-prefixed working files are skipped in the verification tree, but in
            # the extension they are the launch and build records (_md_07_snp.py, _launch_build.sh)
            if f.name in SKIP_NAMES or (_is_scratch(f) and not is_ext):
                skipped += 1
                continue
            rel_src = str(f.relative_to(ROOT).as_posix())
            reasons = screen(f)
            if reasons and rel_src in REVIEWED_AND_CLEARED:
                cleared.append((rel_src, "; ".join(reasons),
                                REVIEWED_AND_CLEARED[rel_src]))
                reasons = []
            if reasons:
                refused.append((rel_src, "; ".join(reasons)))
                continue
            rel = f.relative_to(src)
            dst = out / dst_rel / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(f, dst)
            stored = as_stored(dst)
            size = len(stored)
            total += size
            manifest.append({
                "path": str((Path(dst_rel) / rel).as_posix()),
                "source": str(f.relative_to(ROOT).as_posix()),
                "bytes": size, "sha256": hashlib.sha256(stored).hexdigest()})

    pre = write_preregistration(out)
    why = screen(out / PREREGISTRATION_DST)
    if why:
        raise SystemExit(f"the pre-registration file failed the deposit screen: {why}")
    manifest.append(pre)
    total += pre["bytes"]

    if cleared:
        print(f"flagged and cleared on review ({len(cleared)}):")
        for path, why, because in cleared:
            print(f"   {path}")
            print(f"      flagged: {why}")
            print(f"      cleared: {because[:88]}...")
        print()
    # One artifact, one path. Several files arrive by two routes -- the S10
    # supplement bundle and the run directory both carry them -- and two paths
    # for one artifact is exactly the divergence this deposit exists to rule
    # out. Deduplicate by content hash rather than by name, so a pair that
    # differs is kept and reported instead of silently collapsed.
    by_hash: dict[str, list[dict]] = {}
    for m in manifest:
        by_hash.setdefault(m["sha256"], []).append(m)
    # Prefer the topically-named directory over the catch-all summary one.
    def rank(entry: dict) -> tuple[int, str]:
        return (1 if "/summary/" in entry["path"] else 0, entry["path"])
    dupes = 0
    for entries in by_hash.values():
        if len(entries) < 2:
            continue
        entries.sort(key=rank)
        keep = entries[0]
        for drop in entries[1:]:
            # Two sources can map to the SAME destination path -- the S10
            # bundle and the run directory now both supply the DiffDock
            # artifacts under orthogonal_engine/. The second copy simply
            # overwrote the first, so there is one physical file and two
            # manifest entries. Unlinking "the duplicate" then deleted the
            # very file the surviving entry points at, and the deposit lost
            # orthogonal_check_result.json -- the artifact the manuscript
            # says its Section 3.5 numbers derive from -- while still listing
            # it. Only unlink when the paths actually differ.
            # Added 2026-10-10 (ADR-AIIM-381, approved by the author): a file that the S10 bundle
            # lists in its own `verification/MANIFEST.json` is never unlinked as a duplicate of a
            # copy elsewhere. That manifest promises its entries at their paths; two of them
            # (the figure 4 and figure 5 provenance records) were removed because a byte-identical
            # copy sits under figures/, and the nested manifest then listed files that were not there.
            if drop["path"] != keep["path"] and (drop.get("source") or "").replace("\\", "/").startswith(
                    "30_manuscript/supplement/S10_verification_2026-08-21/"):
                continue
            if drop["path"] != keep["path"]:
                (out / drop["path"]).unlink(missing_ok=True)
            manifest.remove(drop)
            dupes += 1

    # A manifest entry pointing at a file that is not there is the failure this
    # deposit exists to prevent, and the dedup step above silently produced
    # exactly that once. Fail loudly rather than write a manifest that lies.
    ghosts = [m["path"] for m in manifest if not (out / m["path"]).is_file()]
    if ghosts:
        raise SystemExit(
            "manifest lists files that are not on disk -- refusing to write it:\n  "
            + "\n  ".join(ghosts))

    if refused:
        print("REFUSED -- these files were not deposited:")
        for path, why in refused:
            print(f"   {path}\n      {why}")

    (out / "MANIFEST.json").write_text(json.dumps({
        "assembled": __import__("datetime").date.today().isoformat(),
        "n_files": len(manifest),
        "total_bytes": total,
        "n_files_excluded_by_rule": skipped,
        "n_duplicates_removed": dupes,
        "n_files_refused_by_screening": len(refused),
        "refused": [{"path": p, "reason": w} for p, w in refused],
        "flagged_and_cleared_on_review": [
            {"path": p, "flagged_because": w, "cleared_because": b}
            for p, w, b in cleared],
        "external_inputs_not_deposited": EXTERNAL,
        "files": manifest,
    }, indent=2), encoding="utf-8")

    # --- PATH_MAP.json ------------------------------------------------------
    # Most of the evidence reaches this deposit through the S10 supplement
    # bundle, so MANIFEST.json records its `source` as a path under
    # `30_manuscript/supplement/`, not as the run directory the analysis
    # scripts actually name. Without the missing link a reader who has only
    # the deposit cannot get from `20_new_runs/<run>/<file>` to the file that
    # answers to it. build_supplement_s10.origin_map() supplies it, derived
    # from that script's own ITEMS table rather than restated here, so the two
    # cannot drift apart.
    from build_supplement_s10 import origin_map                     # noqa: E402
    s10_prefix = "30_manuscript/supplement/S10_verification_2026-08-21/"
    origins = origin_map()
    path_map: dict[str, str] = {}
    for m in manifest:
        src = (m.get("source") or "").replace("\\", "/")
        if src.startswith(s10_prefix):
            dev = origins.get(src[len(s10_prefix):])
            if dev:
                path_map[dev] = m["path"]
        elif src:
            path_map.setdefault(src, m["path"])
    (out / "PATH_MAP.json").write_text(json.dumps({
        "what": "development-tree path -> path inside this deposit",
        "why": "the scripts name their inputs by the layout they were written "
               "against; this deposit reorganises the same files. "
               "kdss_paths.resolve() reads this so a clone resolves its own "
               "layout instead of the author's machine.",
        "derived_from": ["MANIFEST.json source fields",
                         "build_supplement_s10.origin_map()"],
        "n_entries": len(path_map),
        "map": dict(sorted(path_map.items())),
    }, indent=2), encoding="utf-8")

    print(f"\ndeposited {len(manifest)} files, {total / 1024 / 1024:.1f} MB -> {out}")
    print(f"path map: {len(path_map)} development-tree paths resolve here")
    print(f"excluded by rule: {skipped}   refused by screening: {len(refused)}")
    for grp in sorted({m['path'].split('/')[0] for m in manifest}):
        n = sum(1 for m in manifest if m["path"].startswith(grp + "/"))
        b = sum(m["bytes"] for m in manifest if m["path"].startswith(grp + "/"))
        print(f"   {grp:14s} {n:4d} files  {b/1024:8.0f} KB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

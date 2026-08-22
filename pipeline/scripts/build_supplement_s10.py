"""
build_supplement_s10.py -- assemble the supplementary bundle for the
verification work the v6 manuscript rests on.

The existing supplement (S1-S9) was assembled on 2026-08-08 and predates every
run this manuscript's verification claims come from. Rather than edit that
bundle in place -- it is the record of what was submitted as v4 -- this adds a
tenth section carrying only the 2026-08-21 evidence, with a manifest that
records the SHA-256 of every file so a reader can check that the copy in the
supplement is the file the analysis actually read.

Copied, not regenerated: a supplement that re-derives its own numbers proves
nothing about the run that produced them.

Usage:
    python build_supplement_s10.py
"""
from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path


def _root() -> Path:
    marker = "00_project_control/08_DECISION_LOG.md"
    for cand in (Path(r"D:\Agentic VS Workflow"), Path("/mnt/d/Agentic VS Workflow")):
        if (cand / marker).is_file():
            return cand
    raise SystemExit("project root not found")


ROOT = _root()
RUNS = ROOT / "20_new_runs"
DEST = ROOT / "30_manuscript/supplement/S10_verification_2026-08-21"

# (source relative to 20_new_runs, subfolder in S10, what it supports)
ITEMS: list[tuple[str, str, str]] = [
    ("verification_summary_2026-08-21/mutation_score_by_gate.json",
     "seeded_defects", "Table 4 and §3.1: per-gate and per-family detection, the 2x2 ablation"),
    ("WP4_retrain_postfix_2026-08-05/wp4_conditions_ac_results_postfix_2026-08-05.json",
     "seeded_defects", "the underlying run the per-gate table is derived from"),
    ("naive_baseline_2026-08-10/naive_baseline_summary.json",
     "seeded_defects", "§3.1: the ungoverned comparator, SEER 0.477"),
    ("real_artifact_replay_2026-08-05/real_artifact_replay_results.json",
     "seeded_defects", "Table 4: 4/4 real-artifact defect replay, and the 0.379 A redock RMSD"),

    ("metamorphic_2026-08-21/metamorphic_results.json",
     "metamorphic", "§3.1: the four invariants and the SMILES-rewriting violation"),
    ("metamorphic_2026-08-21/metamorphic_inputs_metadata.json",
     "metamorphic", "how the two input variants were constructed, with InChIKey verification"),
    ("metamorphic_2026-08-21/ledger_M1_rowshuffle.csv", "metamorphic", "M1 input"),
    ("metamorphic_2026-08-21/ledger_M2_smiles_rewrite.csv", "metamorphic", "M2 input"),
    ("metamorphic_2026-08-21/scores_M1.csv", "metamorphic", "M1 output"),
    ("metamorphic_2026-08-21/scores_M2.csv", "metamorphic", "M2 output"),

    ("g4_verification_2026-08-21/G4_contract_outcome_passA_cpu6.json",
     "g4_contract", "§3.1: G4 PASS 7/7, the first gate the natural-product arm passed"),
    ("g4_verification_2026-08-21/run_metadata.json",
     "g4_contract", "the seeded 5% subset rule and the two thread settings"),
    ("g4_verification_2026-08-21/rerun_passA_cpu6.csv", "g4_contract", "6-thread rerun"),
    ("g4_verification_2026-08-21/rerun_passB_cpu20.csv", "g4_contract", "20-thread rerun"),

    ("ad_sweep_expanded_2026-08-21/ad_coverage_by_descriptor.csv",
     "applicability_domain", "Figure 3A and §3.4: coverage in all eight spaces"),
    ("ad_sweep_expanded_2026-08-21/ad_per_compound.csv",
     "applicability_domain", "Figure 3B/3C: per-compound similarity and in-domain flags"),
    ("ad_sweep_expanded_2026-08-21/comparison_vs_snplib.json",
     "applicability_domain", "§3.4: coverage before vs after similarity selection"),
    ("ad_sweep_expanded_2026-08-21/run_metadata.json",
     "applicability_domain", "the single rule applied identically in every space"),
    ("ad_sweep_expanded_2026-08-21/figure_ad_multispace_provenance.json",
     "applicability_domain", "Figure 3: input and output digests"),

    ("funnel_exclusion_control_2026-08-21/exclusion_control_summary.json",
     "funnel_control", "§3.4: 7.0% of discarded compounds reach the selection cutoff vs 10.0% kept"),
    ("funnel_exclusion_control_2026-08-21/control_metadata.json",
     "funnel_control", "how the 500 discarded compounds were sampled"),
    ("funnel_exclusion_control_2026-08-21/exclusion_scores.csv",
     "funnel_control", "the control docking scores"),

    ("g6_selection_2026-08-21/G6_contract_outcome_tierA.json",
     "g6_contract", "§3.4: the corrected gate applied to the released selection"),

    ("qsar_lbdd_2026-08-07/qsar_y_randomization.json",
     "y_randomisation", "§2.4: 100 permutations, z = 21.9, p = 0.0099"),

    ("chembl245_ifp_rescoring_2026-08-15/ifp_rescoring_summary.json",
     "pose_rescue", "Table 1: the five pre-registered rescue methods"),
    ("chembl245_ifp_rescoring_2026-08-15/multi_reference_ifp_summary.json",
     "pose_rescue", "Table 1: the second-reference variants"),
    ("ensemble_shell_pilot_2026-08-15/pilot_result.json",
     "pose_rescue", "Table 1: the water-shell ensemble pilot and its stop rule"),
    ("ensemble_shell_pilot_2026-08-15/PREREGISTRATION.json",
     "pose_rescue", "the decision rule, written down before the pilot ran"),

    ("np_library_funnel_2026-08-15/e8_np_library_prioritization.csv",
     "released_result", "Tables 2 and 3: the prioritisation this study releases"),
    ("np_library_funnel_2026-08-15/tierA_convergence_analysis.json",
     "released_result", "§3.5: the count is unremarkable, the identity is not"),
    ("np_library_funnel_2026-08-15/tierA_literature_crosscheck.md",
     "released_result", "§3.5: the post-hoc PubMed cross-check, hits and non-hits"),
    ("np_library_funnel_2026-08-15/conformal_qsar_summary_np_library_dock_ledger.json",
     "released_result", "§3.5: conformal validity and the exchangeability caveat"),

    ("verification_summary_2026-08-21/structure_and_protocol_checks.json",
     "structure_and_protocol", "§2.3/§2.5: receptor constructs, redock RMSD, stereocentre counts"),
    ("verification_summary_2026-08-21/box_and_cutoff_checks.json",
     "structure_and_protocol", "§2.5: ligand extent vs box, centroid-cutoff sensitivity"),
    ("verification_summary_2026-08-21/ligand_extents.csv",
     "structure_and_protocol", "§2.5: per-compound principal extents"),
    ("verification_summary_2026-08-21/g8_number_traceability.json",
     "structure_and_protocol", "the G8 audit of every number in the manuscript"),
]


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    # incident_taxonomy/ is hand-authored here rather than copied from a run
    # directory (Section 3.6 summarises it; the full table lives only in the
    # supplement). Nothing below writes it, so it simply has to survive.
    DEST.mkdir(parents=True, exist_ok=True)
    manifest: list[dict] = []
    missing: list[str] = []
    total = 0

    for rel, sub, purpose in ITEMS:
        src = RUNS / rel
        if not src.is_file():
            missing.append(rel)
            continue
        outdir = DEST / sub
        outdir.mkdir(parents=True, exist_ok=True)
        dst = outdir / src.name
        shutil.copy2(src, dst)
        size = dst.stat().st_size
        total += size
        manifest.append({
            "file": f"{sub}/{src.name}",
            "source": f"20_new_runs/{rel}",
            "supports": purpose,
            "bytes": size,
            "sha256": sha256(dst),
        })

    readme = [
        "# S10 — Verification evidence, 2026-08-21",
        "",
        "Sections S1–S9 were assembled on 2026-08-08 and predate every run the",
        "verification claims in this manuscript rest on. They are left unmodified as",
        "the record of what accompanied the previous version. This section carries the",
        "2026-08-21 evidence and nothing else.",
        "",
        "Every file is a **copy of the artifact the analysis read**, not a regenerated",
        "one, and the manifest records its SHA-256 so a reader can confirm the copy",
        "matches. A supplement that re-derives its own numbers demonstrates nothing",
        "about the run that produced them.",
        "",
        f"**{len(manifest)} files, {total / 1024:.0f} KB.**",
        "",
        "| File | Supports |",
        "|---|---|",
    ]
    for m in manifest:
        readme.append(f"| `{m['file']}` | {m['supports']} |")
    readme += [
        "",
        "## What is deliberately not here",
        "",
        "- **Held-out fault families H01–H12.** Sealed; unsealing requires a separate",
        "  human approval that has not been given.",
        "- **Molecular-dynamics artifacts.** No MD result is reported in this",
        "  manuscript (§4, deferred work), so none is deposited.",
        "- **Full pose archives.** The per-compound pose directories run to several GB.",
        "  They are part of the repository named in the data-availability statement,",
        "  not of this bundle.",
        "",
        "## Scripts that produced these artifacts",
        "",
        "All under `10_new_system/scripts/`:",
        "`run_g4_verification.sh`, `build_metamorphic_ledgers.py`,",
        "`run_metamorphic_docking.sh`, `analyse_metamorphic.py`,",
        "`ad_descriptor_sweep.py`, `funnel_exclusion_control.py`,",
        "`run_g6_on_np_library.py`, `qsar_y_randomization.py`,",
        "`figure_ad_multispace.py`, `g8_number_traceability.py`.",
        "",
    ]
    (DEST / "README.md").write_text("\n".join(readme), encoding="utf-8")
    (DEST / "MANIFEST.json").write_text(
        json.dumps({"assembled": "2026-08-21", "n_files": len(manifest),
                    "total_bytes": total, "missing_sources": missing,
                    "files": manifest}, indent=2), encoding="utf-8")

    print(f"copied {len(manifest)} files, {total / 1024:.0f} KB -> {DEST}")
    if missing:
        print(f"MISSING ({len(missing)}):")
        for m in missing:
            print("   ", m)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

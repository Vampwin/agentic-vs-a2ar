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
import re
import shutil
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
RUNS = dev_dir("20_new_runs")
DEST = resolve("30_manuscript/supplement/S10_verification_2026-08-21")

# (source relative to 20_new_runs, subfolder in S10, what it supports)
def not_here_block() -> list[str]:
    """The exclusions, rendered from the one place they are decided.

    This section was three hand-written bullets until 2026-08-29. The same
    decisions were also recorded in comments in build_deposit.py and in prose in
    two governance files, and the four copies had begun to disagree -- the
    bullets said the candidate simulations "were still running at submission",
    which stopped being true, and said nothing at all about the runs this
    project performed and then set aside. `evidence_not_used.py` is now the one
    place a decision to leave something out is written down, and this is one of
    its two readers (ADR-AIIM-142).

    Only `artifacts_only` entries name something a reader might expect to find
    in this bundle. `whole_run` entries are runs the manuscript does not rest on
    at all; they are listed too, compactly, because a supplement that shows only
    the work that worked is not a complete account of the work.
    """
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "evidence_not_used", Path(__file__).with_name("evidence_not_used.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.self_test()

    withheld = {k: v for k, v in mod.NOT_USED.items()
                if v.get("scope") == "artifacts_only"}
    unused = {k: v for k, v in mod.NOT_USED.items()
              if v.get("scope", "whole_run") == "whole_run"}

    L = ["## What is deliberately not here", "",
         "Generated from `10_new_system/scripts/evidence_not_used.py`, which is "
         "where each of these decisions is recorded with its reason and the "
         "decision record that made it. Nothing here is summarised.", "",
         "### Reported, but the material is withheld", ""]
    for run, e in sorted(withheld.items()):
        L.append(f"**`{run}`** — {e['kind'].replace('_', ' ')}, excluded from "
                 f"{' and '.join(e.get('excluded_from', []))} "
                 f"({', '.join(e.get('adr', []))}).")
        L.append("")
        L.append(e["reason"])
        if e.get("published_anyway"):
            L.append("")
            L.append("Published anyway: "
                     + ", ".join(f"`{q}`" for q in e["published_anyway"]) + ".")
        L.append("")
    L += ["### Runs this manuscript does not rest on", "",
          "Produced during this study and not used. Listed so that the work "
          "which did not survive is visible beside the work that did.", "",
          "| Run | Why |", "|---|---|"]
    for run, e in sorted(unused.items()):
        # Not truncated to a first sentence. The heading three lines above
        # promises nothing here is summarised, and for several of these the
        # reason the run was set aside is in the second sentence, not the first.
        why = " ".join(e["reason"].split())
        if e.get("superseded_by"):
            why += f" Superseded by `{e['superseded_by']}`."
        L.append(f"| `{run}` | {why} |")
    L.append("")
    return L


ITEMS: list[tuple[str, str, str]] = [
    ("seeded_defect_recount_2026-08-29/mutation_score_by_gate_recounted.json",
     "seeded_defects", "Table 1 and §3.1: per-gate and per-family detection and the "
     "2x2 ablation, over the 128 cases in which a fault was actually injected. "
     "Carries the superseded figures beside the corrected ones (ADR-AIIM-131)"),
    ("verification_summary_2026-08-21/mutation_score_by_gate.json",
     "seeded_defects", "SUPERSEDED 2026-08-29. The same analysis divided by 133 "
     "rather than 128, because the export it read had dropped the vacuous flag. "
     "Published so the correction can be checked against what it replaced"),
    ("WP4_retrain_postfix_2026-08-05/wp4_conditions_ac_results_postfix_2026-08-05.json",
     "seeded_defects", "the underlying run the per-gate table is derived from"),
    ("naive_baseline_2026-08-29/naive_baseline_summary.json",
     "seeded_defects", "§3.1: the naive-analyst comparator, SEER 0.477 over the "
     "same 128 cases as the governed arm. Re-run by a driver script "
     "(ADR-AIIM-131); reproduces the 2026-08-10 result exactly, family by family"),
    ("naive_baseline_2026-08-29/naive_baseline_outcomes.csv",
     "seeded_defects", "per-case outcomes of that comparator"),
    ("naive_baseline_2026-08-10/naive_baseline_summary.json",
     "seeded_defects", "the original 2026-08-10 run, which had no driver script. "
     "Kept so the reproduction can be checked against what it reproduces"),
    ("real_artifact_replay_2026-08-05/real_artifact_replay_results.json",
     "seeded_defects", "Table 4: 4/4 real-artifact defect replay, and the 0.379 A redock RMSD"),
    ("control_pose_fragments_2026-08-28/control_pose_fragments.json",
     "structure_checks", "§3.6 and Figure 5: the 1.92 Å control-pose RMSD "
     "decomposed into core 1.03, furanyl 2.88 and tail 2.00 Å"),
    ("pose_disagreement_2026-08-29/pose_disagreement_geometry.json",
     "structure_checks", "§3.6 and Figure 5B: centroid separation beside "
     "in-place RMSD — the 8.80 Å disagreement is reorientation, not displacement"),
    ("WP2.5_WP3A_readiness_2026-07-30/fit_selected_model_summary.json",
     "learned_monitor", "§4: the fitted risk monitor, 188 rows and 24 features"),
    ("WP2.5_WP3A_readiness_2026-07-30/calibration_comparison.json",
     "learned_monitor", "§4: isotonic against sigmoid; isotonic on 188 rows "
     "collapses the probability to four values"),
    ("WP2.5_WP3A_readiness_2026-07-30/compare_models_dev_data.json",
     "learned_monitor", "§4: the family-grouped nested-CV model comparison and "
     "the prespecification override of gradient boosting"),
    ("WP2.5_WP3A_readiness_2026-07-30/policy_config.json",
     "learned_monitor", "§4: the two policy thresholds, which are two of the "
     "four probabilities isotonic can emit"),
    ("smiles_redock_500_2026-08-26/r2_results.json",
     "metamorphic", "§3.2 and §3.6: all 500 re-docked from rewritten SMILES. "
     "Tier-A identity unchanged; 3 of 497 scores identical, ρ = 0.767; the three "
     "reported convergence failures all completed (ADR-AIIM-157)"),
    ("funnel_recall_fullprotocol_2026-09-01/r4_results.json",
     "applicability_domain", "§3.5: the discarded sample at the production "
     "protocol. 10.76 % weighted against the fast pass's 7.41 %, and the "
     "protocol-shift decomposition showing the fast pass under-scored what it "
     "rejected three times as much as what it kept (ADR-AIIM-158)"),
    ("md_pose_persistence_2026-08-31/pocket_occupancy.json",
     "structure_checks", "§3.6: the PBC-proof measurement of whether each "
     "ligand stayed in the site it was docked into (ADR-AIIM-156)"),
    ("trivial_baseline_2026-09-03/trivial_baseline.json",
     "screening_performance", "§3.4, §4, §5: what Table 5's baseline is worth. "
     "Vina 0.601 recomputed as a self-check, heavy-atom count 0.584 on the same "
     "245 rows, paired difference +0.018 (p = 0.722). Not pre-registered and "
     "labelled so (ADR-AIIM-164)"),
    ("pose_validity_2026-09-03/pose_validity.json",
     "structure_checks", "§3.5, §3.6: PoseBusters over all 4,916 "
     "full-protocol poses, tool defaults unmodified. 97.88 % valid; the DiffDock "
     "control that passed the 2.0 A gate overlaps the receptor (ADR-AIIM-166)"),
    ("clinical_reference_2026-09-03/clinical_reference.json",
     "screening_performance", "§3.5, §4: ten A2A ligands of known affinity "
     "docked under the production protocol. None of the seven clinical and "
     "reference antagonists is placed in the Asn253 hydrogen bond; the apparent "
     "score-affinity correlation is carried by two agonists (ADR-AIIM-163)"),
    ("crossdock_zma_2026-09-03/crossdock_zma.json",
     "structure_checks", "§3.6: does redocking accuracy survive a change of "
     "receptor conformation? Control 0.397 A; one clean cross-dock at 0.318 A; "
     "the frozen structure set cannot answer the question (ADR-AIIM-162)"),
    ("ood_sensitivity_2026-09-05/ood_sensitivity.json",
     "seeded_defects", "§2.1, §4: sensitivity of the frozen "
     "applicability-domain detector, which had only ever been measured "
     "in-distribution. Flag rate against displacement factor k; saturates at "
     "5.32 % and does not move at k=1000. Negative control at k=1 reproduces "
     "the recorded training rate (ADR-AIIM-171)"),
    ("redock_9xt_2026-09-05/redock_ligand_transfer.json",
     "structure_checks", "§3.6, §4: does redocking accuracy survive a change of "
     "*ligand*? Every structure in the cross-docking study was solved with "
     "ZM241385. 9XT redocks into its own crystal at 0.216 A and into 4EIY at "
     "9.606 A, with no native-like pose among ten. Criteria and the reading of "
     "each outcome fixed before the run (ADR-AIIM-168)"),
    ("redock_9xt_2026-09-05/prep.json",
     "structure_checks", "receptor preparation for that run: shared-CA "
     "superposition, the water shell selected around each structure's own "
     "co-crystal ligand, and the pdbfixer atom and water counts"),
    ("redock_9xt_2026-09-05/manifest.json",
     "structure_checks", "every input to that run, hashed before anything was "
     "computed"),
    ("validation_matrix_2026-09-05/validation_matrix.json",
     "screening_performance", "Table 7: every technique on five capability "
     "axes, following CASF's decomposition with physical validity added. Each "
     "cell names the artifact it was read from, or says *not measured*; the "
     "file computes nothing itself (ADR-AIIM-169)"),
    ("md_recognition_motif_2026-09-02/recognition_motif.json",
     "structure_checks", "§3.6, §4: Asn253 hydrogen-bond persistence for all "
     "nine candidate replicates against the ZM241385 control simulated "
     "identically. Criteria pre-registered in ADR-AIIM-161 before any "
     "trajectory was read; PARTIAL for all three. Carries both "
     "pre-registered secondaries, neither of which could discriminate"),
    ("g6_selection_2026-09-02/G6_contract_outcome_tierA.json",
     "g4_contract", "§2.6: the G6 selection contract at v1.3.0, the version "
     "that carries the false-refusal check. Supersedes the 2026-08-21 outcome, "
     "which ran v1.2.0 and is kept beside it (ADR-AIIM-159)"),
    ("funnel_recall_reweighted_2026-08-31/funnel_recall_reweighted.json",
     "applicability_domain", "§3.5: the funnel's exclusions weighted by stratum "
     "— 7.41 % against the pooled 7.0 %, and about 1,251 discarded compounds "
     "reaching the cutoff against 500 selected (ADR-AIIM-150)"),
    ("md_pose_persistence_2026-08-31/md_pose_persistence.json",
     "structure_checks", "§3.6, §4 and §5: pose persistence over nine candidate "
     "replicates. NPLIB01185 leaves the pocket in all three; the other two "
     "remain at 4.73 and 5.22 Å from their docked poses against the control's "
     "2.59 Å. No verdict attached — no threshold was pre-registered"),
    ("domain_gate_false_refusal_2026-08-29/domain_gate_false_refusal.json",
     "applicability_domain", "§2.6 and §4: what the applicability-domain gate "
     "does to 1,096 known A2A actives whose scaffolds define no part of the "
     "domain. 5.8-15.0 % refused, against 0 % of the natural-product library "
     "admitted in every structural space (ADR-AIIM-147)"),
    ("conformer_budget_sweep_2026-08-27/r3_results.json",
     "metamorphic", "§3.1: how much of the SMILES dependence is the conformer "
     "budget. Six arms at 8/32/128 conformers, criteria fixed before the run "
     "(ADR-AIIM-133); outcome PARTIAL, ratio 0.491"),
    ("tierA_conditioned_null_2026-08-29/tierA_conditioned_null.json",
     "seeded_defects", "§3.6 and §4: the chemical-coherence claim tested two ways. "
     "The docking-conditioned permutation (unchanged at p = 0.0023) and the "
     "multiplicity-aware version (p = 0.126 at genus, 0.182 at superclass). "
     "Reconstructs the 23-compound reference set the frozen artifact never "
     "documented (ADR-AIIM-145)"),
    ("heldout_case_analysis_2026-08-29/heldout_case_analysis.json",
     "seeded_defects", "§3.1: why the five blind cases are reported as a case "
     "study and not a rate — the Wilson interval, the two cases on gates with no "
     "development fault family, and the two missed at gates scoring 100 % "
     "(ADR-AIIM-135). Per-case gate and observed state only; no fault mechanism"),
    ("gate_reproducibility_2026-08-29/gate_reproducibility.json",
     "seeded_defects", "S3.1 and Table 1: the two-execution comparison. Four "
     "contract outcomes identical, nine of thirteen artifacts byte-identical, "
     "and the receptor non-determinism that leaves the decisions unchanged"),

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

    ("../30_manuscript/figures/figure4_screening_performance_provenance.json",
     "pose_rescue", "Figure 4: input and output digests, and the AUC value under both tie conventions"),

    ("../30_manuscript/figures/figure5_pose_overlay_provenance.json",
     "orthogonal_engine", "Figure 5: input and output digests, RMSD values and their verification "
     "method, and the chain-discontinuity note behind the pocket-local cartoon view"),

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
    ("chembl245_ifp_rescoring_2026-08-15/ifp_rescoring_per_compound.csv",
     "pose_rescue", "Figure 4A and 4C: per-compound scores and labels, n=245"),
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

    # Added 2026-08-24: artifacts the v6 manuscript describes that this bundle
    # did not carry. Each is copied from the run that produced it, like the rest.
    ("diffdock_orthogonal_2026-08-23/orthogonal_check_result.json",
     "orthogonal_engine", "§3.5: the pre-registered DiffDock check -- control gate 1.92 A, "
     "candidate disagreement 6.43/8.80/9.50 A, and the pre-registered prediction scored as failed"),
    ("diffdock_orthogonal_2026-08-23/diffdock_input.csv",
     "orthogonal_engine", "§2.5: the four ligands submitted, as InChIKey-verified SMILES"),
    ("diffdock_orthogonal_2026-08-23/4EIY_protein_only.pdb",
     "orthogonal_engine", "§2.5: the receptor given to DiffDock, waters removed as declared"),

    ("heldout_evaluation_2026-08-08/heldout_evaluation_summary.json",
     "heldout_outcome", "§3.1: which of the five project-blind fault cases each gate detected"),
    ("heldout_evaluation_2026-08-08/fix_verification_summary.json",
     "heldout_outcome", "§3.1: the regression check confirming the H01-H03 patches -- "
     "explicitly NOT a second held-out measurement"),

    ("md_candidates_2026-08-22/NPLIB01185/build_metadata.json",
     "md_candidate_builds", "§2.9: identity check against the frozen ledger and parmchk2 coverage"),
    ("md_candidates_2026-08-22/NPLIB01643/build_metadata.json",
     "md_candidate_builds", "§2.9: identity check against the frozen ledger and parmchk2 coverage"),
    ("md_candidates_2026-08-22/NPLIB03056/build_metadata.json",
     "md_candidate_builds", "§2.9: identity check against the frozen ledger and parmchk2 coverage"),

    ("md_production_2026-08-21/04_run/analysis_prod/ligand_rmsd.xvg",
     "md_control", "§2.9: ZM241385 heavy-atom RMSD over the 100 ns control trajectory"),
    ("md_production_2026-08-21/04_run/analysis_prod/box.xvg",
     "md_control", "§2.9: box dimensions, the bilayer-equilibration check"),

    # --- added by the 2026-08-25/26 remediation. Each is load-bearing for a
    # number in the manuscript and was, until now, carried by no manifest.
    ("ad_sweep_selfsim_2026-08-25/ad_percentile_position.csv",
     "applicability_domain", "§3.5 and Figure 4B: where each compound sits inside its "
     "own space's training self-similarity distribution -- the rank view the "
     "admission rates are a threshold on (ADR-AIIM-126)"),
    ("ad_sweep_selfsim_2026-08-25/comparison_vs_snplib.json",
     "applicability_domain", "Figure 4A: the two libraries side by side, now with a "
     "producer script rather than hand-assembled (ADR-AIIM-136)"),
    ("ad_sweep_selfsim_2026-08-25/run_provenance.json",
     "applicability_domain", "§2.1: tooling versions and command line for the sweep, "
     "which the 2026-08-21 run did not record"),

    ("conformal_diagnostics_2026-08-25/conformal_diagnostics.json",
     "released_result", "§3.6: empty-set rate, prediction-set size and per-class "
     "p-value distributions on and off domain, and where the Tier-A criterion sits "
     "in the distribution it thresholds (ADR-AIIM-137)"),
    ("conformal_diagnostics_2026-08-25/calibration_pvalues.csv",
     "released_result", "§3.6: per-compound cross-conformal p-values on the held-out "
     "set, the in-domain comparator for the query libraries"),

    ("fastpass_correlation_2026-08-25/fastpass_correlation.json",
     "funnel_control", "§3.5: agreement between the fast-pass and full-protocol "
     "scores, and the residual spread against the selection boundary (ADR-AIIM-137)"),
    ("fastpass_correlation_2026-08-25/fastpass_vs_fullprotocol.csv",
     "funnel_control", "§3.5: the per-compound pairing behind that comparison"),

    ("qsar_calibration_2026-08-26/qsar_calibration.json",
     "qsar_calibration", "§2.4: Brier score, reliability curve and the reconciliation "
     "of the three AUC figures the manuscript quotes (ADR-AIIM-140)"),
    ("qsar_calibration_2026-08-26/reliability_curve.csv",
     "qsar_calibration", "§2.4: the ten-bin reliability curve, including the two bins "
     "the Tier-A probabilities fall in"),
    # --- added 2026-09-13: read by g8_number_traceability.py but not shipped -
    ("ad_descriptor_sweep_2026-08-15/ad_coverage_by_descriptor.csv",
     "applicability_domain", "§3.3 and the graphical abstract: per-space coverage "
     "for the 236-compound SNP-Lib itself -- 0 % in four structural spaces and "
     "152 of 236 (64.4 %) physicochemically. Distinct from the expanded "
     "selection's table beside it, and the file the graphical abstract was "
     "corrected against (ADR-AIIM-175)"),
    ("redock_9xt_2026-09-05/scores_control.csv",
     "structure_checks", "§3.6: the ZM241385 control arm that gates the "
     "ligand-transfer redocking"),
    ("redock_9xt_2026-09-05/scores_T1.csv",
     "structure_checks", "§3.6: 9XT redocked into 5OLH, the structure solved "
     "with it -- 0.216 A"),
    ("redock_9xt_2026-09-05/scores_T2.csv",
     "structure_checks", "§3.6: the same ligand into 4EIY, a structure solved "
     "with a different one -- 9.606 A"),
    ("redock_9xt_2026-09-05/poses_T2/9XT.sdf",
     "structure_checks", "§3.6: the failing pose itself, so the 9.606 A is a "
     "geometry a reader can inspect rather than a number to accept"),
    ("WP2.5_WP3A_readiness_2026-07-30/ood_fit_summary.json",
     "learned_monitor", "§3.2: the isolation-forest applicability-domain "
     "detector as fitted, including the contamination setting its saturation "
     "behaviour follows from"),
    ("WP2.5_WP3A_readiness_2026-07-30/wp4_conditions_ac_results.json",
     "seeded_defects", "the pre-recount conditions A and C, published so the "
     "denominator correction of ADR-AIIM-131 can be checked against its input"),
    ("np_library_funnel_2026-08-15/funnel_metadata.json",
     "released_result", "§2.6: the funnel's own record of what entered each "
     "stage and what the cutoffs were"),
    ("md_production_2026-08-21/03_build/system_gmx.gro",
     "structure_checks", "§2.9: the assembled membrane system the simulation "
     "reads. The atom, lipid, water and restraint counts the Methods quote are "
     "counted from this file rather than reported from the build log"),
    ("snp_screen_2026-08-04/G1_standardization/snp_source_ledger.csv",
     "released_result", "the 236-compound SNP-Lib accountability ledger as G1 "
     "produced it. Deposited 2026-09-13 because 42 tests in the fault-injection "
     "suite seed their workspaces from it: without it the suite cannot run from "
     "this deposit, which is the only place a reader can run it"),
    ("external_repro_audit_2026-09-13/coverage_core.json",
     "external_audit", "§4 and the SI's required-analyses section: branch-aware "
     "coverage of the core packages from the 2026-09-13 clean-checkout audit "
     "(ADR-AIIM-189). 75 % overall, but `contract_engine` at 29 % and "
     "`orchestrator` at 25 %, and the lowest evaluator module at 39 %"),
    ("external_repro_audit_2026-09-13/manifest_verification.json",
     "external_audit", "§4: byte-level verification of the deposited manifests. "
     "180 of 369 root entries mismatch, 179 of them explained by CRLF line "
     "endings rather than by changed content, and the nested manifest names "
     "two absent files (ADR-AIIM-189)"),
    ("external_repro_audit_2026-09-13/reference_check.json",
     "external_audit", "the reference audit behind the citation check: 66 "
     "references resolved by DOI, 57 PASS, 3 adjudicated and 6 unverifiable "
     "for want of a DOI. Only the DOI was sent externally (ADR-AIIM-189)"),
]# Hand-authored, not copied from a run directory, but still manifest-listed so
# that the manuscript's claim -- every file in S10 carries a SHA-256 in the
# section manifest -- is true of the whole section rather than most of it.
HAND_AUTHORED: list[tuple[str, str]] = [
    ("incident_taxonomy/incident_taxonomy.md",
     "§3.6: the fifteen incidents classified by layer, what detected each, how it resolved"),
]


# Changed 2026-10-10 (ADR-AIIM-380, approved by the author), with the same rule as build_deposit.py
# (ADR-AIIM-379). This manifest used to hash the development copy, which has Windows (CRLF) line
# ends, while the deposit's .gitattributes (`* text=auto eol=lf`) makes git store, and every clone
# check out, LF: 84 of the 100 entries a clone holds did not reproduce. It now records what a clone
# contains: for a text file the SHA-256 and byte count of its content with CRLF turned into LF, for a
# binary file (a suffix .gitattributes marks binary, a NUL byte in the first 8,000 bytes, or a carriage
# return that is not part of CRLF) those of the bytes as they are.
BINARY_SUFFIX = {".png", ".pdf", ".gz", ".docx"}


def as_stored(p: Path) -> bytes:
    """The bytes a clone of the deposit will contain for this file."""
    raw = p.read_bytes()
    if p.suffix.lower() in BINARY_SUFFIX or b"\0" in raw[:8000]:
        return raw
    if re.search(rb"\r(?!\n)", raw):
        return raw
    return raw.replace(b"\r\n", b"\n")


def sha256(p: Path) -> str:
    return hashlib.sha256(as_stored(p)).hexdigest()


def main() -> int:
    # incident_taxonomy/ is hand-authored here rather than copied from a run
    # directory (Section 3.6 summarises it; the full table lives only in the
    # supplement). Nothing below writes it, so it simply has to survive.
    DEST.mkdir(parents=True, exist_ok=True)
    manifest: list[dict] = []
    missing: list[str] = []
    total = 0

    from collections import Counter
    COLLIDING = Counter((sub, Path(rel).name) for rel, sub, _ in ITEMS)

    for rel, sub, purpose in ITEMS:
        src = RUNS / rel
        if not src.is_file():
            missing.append(rel)
            continue
        outdir = DEST / sub
        outdir.mkdir(parents=True, exist_ok=True)
        # Disambiguate colliding basenames by their source directory. The three
        # candidate builds are all `build_metadata.json`; naming the copies by
        # basename alone silently overwrote two of them. Every member of a
        # colliding set is prefixed, not just the later ones -- otherwise one
        # file keeps the bare name and a reader cannot tell which it belongs to.
        dst = outdir / (f"{src.parent.name}_{src.name}"
                        if COLLIDING[(sub, src.name)] > 1 else src.name)
        shutil.copy2(src, dst)
        size = len(as_stored(dst))
        total += size
        manifest.append({
            "file": f"{sub}/{dst.name}",
            "source": f"20_new_runs/{rel}",
            "supports": purpose,
            "bytes": size,
            "sha256": sha256(dst),
        })

    for rel, purpose in HAND_AUTHORED:
        f = DEST / rel
        if not f.is_file():
            missing.append(f"(hand-authored) {rel}")
            continue
        size = len(as_stored(f))
        total += size
        manifest.append({
            "file": rel,
            "source": "hand-authored for this supplement",
            "supports": purpose,
            "bytes": size,
            "sha256": sha256(f),
        })

    readme = [
        "# S10 — Verification evidence",
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
        *not_here_block(),
        "## Scripts that produced these artifacts",
        "",
        "All under `10_new_system/scripts/`:",
        "`run_g4_verification.sh`, `build_metamorphic_ledgers.py`,",
        "`run_metamorphic_docking.sh`, `analyse_metamorphic.py`,",
        "`ad_descriptor_sweep.py`, `funnel_exclusion_control.py`,",
        "`run_g6_on_np_library.py`, `qsar_y_randomization.py`,",
        "`figure_ad_multispace.py`, `g8_number_traceability.py`,",
        "`diffdock_orthogonal_check.py`, `md_11_analyse_pose_persistence.sh`.",
        "",
    ]
    ghosts = [m["file"] for m in manifest if not (DEST / m["file"]).is_file()]
    if ghosts:
        raise SystemExit("manifest lists files that are not on disk -- "
                         "refusing to write it:\n  " + "\n  ".join(ghosts))
    # Orphans: a file under DEST that the manifest does not list. The ghost
    # check above catches a manifest that promises a file it does not have;
    # this catches the opposite, which is worse, because build_deposit.py
    # copies this whole tree recursively -- an unlisted file gets published
    # with no provenance record and no SHA-256. It happens whenever an ITEMS
    # entry is renamed, removed, or picks up a collision prefix: the previous
    # copy is simply left behind. Found 2026-08-28, when adding a second
    # comparison_vs_snplib.json triggered the collision rule and stranded the
    # unprefixed original.
    listed = {m["file"] for m in manifest}
    orphans = sorted(
        str(p.relative_to(DEST).as_posix())
        for p in DEST.rglob("*")
        if p.is_file() and p.name not in ("MANIFEST.json", "README.md")
        and str(p.relative_to(DEST).as_posix()) not in listed)
    for o in orphans:
        (DEST / o).unlink()
        print(f"removed stale file not in the manifest: {o}")

    dupes = [f for f in {m["file"] for m in manifest}
             if sum(1 for m in manifest if m["file"] == f) > 1]
    if dupes:
        raise SystemExit("manifest lists the same path more than once:\n  "
                         + "\n  ".join(dupes))

    (under("30_manuscript/supplement/S10_verification_2026-08-21", "README.md")).write_text("\n".join(readme), encoding="utf-8")
    (under("30_manuscript/supplement/S10_verification_2026-08-21", "MANIFEST.json")).write_text(
        json.dumps({"assembled": "2026-08-21, extended 2026-08-24",
                    "n_files": len(manifest),
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


def origin_map() -> dict[str, str]:
    """S10-relative path -> the development-tree path it was copied from.

    The deposit reaches most of its evidence through this bundle rather than
    from `20_new_runs/` directly, so a reader who has only the deposit cannot
    get from a path the scripts name to the file that answers to it. This is
    the missing half of that chain; `build_deposit.py` writes it out as
    `PATH_MAP.json` and `kdss_paths.resolve()` reads it.

    The destination-naming rule is the one `main()` applies, including the
    prefixing of colliding basenames, and is derived from `ITEMS` rather than
    restated, so the two cannot drift apart.
    """
    from collections import Counter
    colliding = Counter((sub, Path(rel).name) for rel, sub, _ in ITEMS)
    out: dict[str, str] = {}
    for rel, sub, _ in ITEMS:
        name = Path(rel).name
        dst = f"{Path(rel).parent.name}_{name}" if colliding[(sub, name)] > 1 else name
        out[f"{sub}/{dst}"] = f"20_new_runs/{rel}"
    return out

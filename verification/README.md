# S10 — Verification evidence

Sections S1–S9 were assembled on 2026-08-08 and predate every run the
verification claims in this manuscript rest on. They are left unmodified as
the record of what accompanied the previous version. This section carries the
2026-08-21 evidence and nothing else.

Every file is a **copy of the artifact the analysis read**, not a regenerated
one, and the manifest records its SHA-256 so a reader can confirm the copy
matches. A supplement that re-derives its own numbers demonstrates nothing
about the run that produced them.

**99 files, 5806 KB.**

| File | Supports |
|---|---|
| `seeded_defects/mutation_score_by_gate_recounted.json` | Table 1 and §3.1: per-gate and per-family detection and the 2x2 ablation, over the 128 cases in which a fault was actually injected. Carries the superseded figures beside the corrected ones (ADR-AIIM-131) |
| `seeded_defects/mutation_score_by_gate.json` | SUPERSEDED 2026-08-29. The same analysis divided by 133 rather than 128, because the export it read had dropped the vacuous flag. Published so the correction can be checked against what it replaced |
| `seeded_defects/wp4_conditions_ac_results_postfix_2026-08-05.json` | the underlying run the per-gate table is derived from |
| `seeded_defects/naive_baseline_2026-08-29_naive_baseline_summary.json` | §3.1: the naive-analyst comparator, SEER 0.477 over the same 128 cases as the governed arm. Re-run by a driver script (ADR-AIIM-131); reproduces the 2026-08-10 result exactly, family by family |
| `seeded_defects/naive_baseline_outcomes.csv` | per-case outcomes of that comparator |
| `seeded_defects/naive_baseline_2026-08-10_naive_baseline_summary.json` | the original 2026-08-10 run, which had no driver script. Kept so the reproduction can be checked against what it reproduces |
| `seeded_defects/real_artifact_replay_results.json` | Table 4: 4/4 real-artifact defect replay, and the 0.379 A redock RMSD |
| `structure_checks/control_pose_fragments.json` | §3.6 and Figure 5: the 1.92 Å control-pose RMSD decomposed into core 1.03, furanyl 2.88 and tail 2.00 Å |
| `structure_checks/pose_disagreement_geometry.json` | §3.6 and Figure 5B: centroid separation beside in-place RMSD — the 8.80 Å disagreement is reorientation, not displacement |
| `learned_monitor/fit_selected_model_summary.json` | §4: the fitted risk monitor, 188 rows and 24 features |
| `learned_monitor/calibration_comparison.json` | §4: isotonic against sigmoid; isotonic on 188 rows collapses the probability to four values |
| `learned_monitor/compare_models_dev_data.json` | §4: the family-grouped nested-CV model comparison and the prespecification override of gradient boosting |
| `learned_monitor/policy_config.json` | §4: the two policy thresholds, which are two of the four probabilities isotonic can emit |
| `metamorphic/r2_results.json` | §3.2 and §3.6: all 500 re-docked from rewritten SMILES. Tier-A identity unchanged; 3 of 497 scores identical, ρ = 0.767; the three reported convergence failures all completed (ADR-AIIM-157) |
| `applicability_domain/r4_results.json` | §3.5: the discarded sample at the production protocol. 10.76 % weighted against the fast pass's 7.41 %, and the protocol-shift decomposition showing the fast pass under-scored what it rejected three times as much as what it kept (ADR-AIIM-158) |
| `structure_checks/pocket_occupancy.json` | §3.6: the PBC-proof measurement of whether each ligand stayed in the site it was docked into (ADR-AIIM-156) |
| `screening_performance/trivial_baseline.json` | §3.4, §4, §5: what Table 5's baseline is worth. Vina 0.601 recomputed as a self-check, heavy-atom count 0.584 on the same 245 rows, paired difference +0.018 (p = 0.722). Not pre-registered and labelled so (ADR-AIIM-164) |
| `structure_checks/pose_validity.json` | §3.5, §3.6: PoseBusters over all 4,916 full-protocol poses, tool defaults unmodified. 97.88 % valid; the DiffDock control that passed the 2.0 A gate overlaps the receptor (ADR-AIIM-166) |
| `screening_performance/clinical_reference.json` | §3.5, §4: ten A2A ligands of known affinity docked under the production protocol. None of the seven clinical and reference antagonists is placed in the Asn253 hydrogen bond; the apparent score-affinity correlation is carried by two agonists (ADR-AIIM-163) |
| `structure_checks/crossdock_zma.json` | §3.6: does redocking accuracy survive a change of receptor conformation? Control 0.397 A; one clean cross-dock at 0.318 A; the frozen structure set cannot answer the question (ADR-AIIM-162) |
| `seeded_defects/ood_sensitivity.json` | §2.1, §4: sensitivity of the frozen applicability-domain detector, which had only ever been measured in-distribution. Flag rate against displacement factor k; saturates at 5.32 % and does not move at k=1000. Negative control at k=1 reproduces the recorded training rate (ADR-AIIM-171) |
| `structure_checks/redock_ligand_transfer.json` | §3.6, §4: does redocking accuracy survive a change of *ligand*? Every structure in the cross-docking study was solved with ZM241385. 9XT redocks into its own crystal at 0.216 A and into 4EIY at 9.606 A, with no native-like pose among ten. Criteria and the reading of each outcome fixed before the run (ADR-AIIM-168) |
| `structure_checks/prep.json` | receptor preparation for that run: shared-CA superposition, the water shell selected around each structure's own co-crystal ligand, and the pdbfixer atom and water counts |
| `structure_checks/manifest.json` | every input to that run, hashed before anything was computed |
| `screening_performance/validation_matrix.json` | Table 7: every technique on five capability axes, following CASF's decomposition with physical validity added. Each cell names the artifact it was read from, or says *not measured*; the file computes nothing itself (ADR-AIIM-169) |
| `structure_checks/recognition_motif.json` | §3.6, §4: Asn253 hydrogen-bond persistence for all nine candidate replicates against the ZM241385 control simulated identically. Criteria pre-registered in ADR-AIIM-161 before any trajectory was read; PARTIAL for all three. Carries both pre-registered secondaries, neither of which could discriminate |
| `g4_contract/G6_contract_outcome_tierA.json` | §2.6: the G6 selection contract at v1.3.0, the version that carries the false-refusal check. Supersedes the 2026-08-21 outcome, which ran v1.2.0 and is kept beside it (ADR-AIIM-159) |
| `applicability_domain/funnel_recall_reweighted.json` | §3.5: the funnel's exclusions weighted by stratum — 7.41 % against the pooled 7.0 %, and about 1,251 discarded compounds reaching the cutoff against 500 selected (ADR-AIIM-150) |
| `structure_checks/md_pose_persistence.json` | §3.6, §4 and §5: pose persistence over nine candidate replicates. NPLIB01185 leaves the pocket in all three; the other two remain at 4.73 and 5.22 Å from their docked poses against the control's 2.59 Å. No verdict attached — no threshold was pre-registered |
| `applicability_domain/domain_gate_false_refusal.json` | §2.6 and §4: what the applicability-domain gate does to 1,096 known A2A actives whose scaffolds define no part of the domain. 5.8-15.0 % refused, against 0 % of the natural-product library admitted in every structural space (ADR-AIIM-147) |
| `metamorphic/r3_results.json` | §3.1: how much of the SMILES dependence is the conformer budget. Six arms at 8/32/128 conformers, criteria fixed before the run (ADR-AIIM-133); outcome PARTIAL, ratio 0.491 |
| `seeded_defects/tierA_conditioned_null.json` | §3.6 and §4: the chemical-coherence claim tested two ways. The docking-conditioned permutation (unchanged at p = 0.0023) and the multiplicity-aware version (p = 0.126 at genus, 0.182 at superclass). Reconstructs the 23-compound reference set the frozen artifact never documented (ADR-AIIM-145) |
| `seeded_defects/heldout_case_analysis.json` | §3.1: why the five blind cases are reported as a case study and not a rate — the Wilson interval, the two cases on gates with no development fault family, and the two missed at gates scoring 100 % (ADR-AIIM-135). Per-case gate and observed state only; no fault mechanism |
| `seeded_defects/gate_reproducibility.json` | S3.1 and Table 1: the two-execution comparison. Four contract outcomes identical, nine of thirteen artifacts byte-identical, and the receptor non-determinism that leaves the decisions unchanged |
| `metamorphic/metamorphic_results.json` | §3.1: the four invariants and the SMILES-rewriting violation |
| `metamorphic/metamorphic_inputs_metadata.json` | how the two input variants were constructed, with InChIKey verification |
| `metamorphic/ledger_M1_rowshuffle.csv` | M1 input |
| `metamorphic/ledger_M2_smiles_rewrite.csv` | M2 input |
| `metamorphic/scores_M1.csv` | M1 output |
| `metamorphic/scores_M2.csv` | M2 output |
| `g4_contract/G4_contract_outcome_passA_cpu6.json` | §3.1: G4 PASS 7/7, the first gate the natural-product arm passed |
| `g4_contract/run_metadata.json` | the seeded 5% subset rule and the two thread settings |
| `g4_contract/rerun_passA_cpu6.csv` | 6-thread rerun |
| `g4_contract/rerun_passB_cpu20.csv` | 20-thread rerun |
| `applicability_domain/ad_sweep_expanded_2026-08-21_ad_coverage_by_descriptor.csv` | Figure 3A and §3.4: coverage in all eight spaces |
| `applicability_domain/ad_per_compound.csv` | Figure 3B/3C: per-compound similarity and in-domain flags |
| `applicability_domain/ad_sweep_expanded_2026-08-21_comparison_vs_snplib.json` | §3.4: coverage before vs after similarity selection |
| `applicability_domain/run_metadata.json` | the single rule applied identically in every space |
| `applicability_domain/figure_ad_multispace_provenance.json` | Figure 3: input and output digests |
| `pose_rescue/figure4_screening_performance_provenance.json` | Figure 4: input and output digests, and the AUC value under both tie conventions |
| `orthogonal_engine/figure5_pose_overlay_provenance.json` | Figure 5: input and output digests, RMSD values and their verification method, and the chain-discontinuity note behind the pocket-local cartoon view |
| `funnel_control/exclusion_control_summary.json` | §3.4: 7.0% of discarded compounds reach the selection cutoff vs 10.0% kept |
| `funnel_control/control_metadata.json` | how the 500 discarded compounds were sampled |
| `funnel_control/exclusion_scores.csv` | the control docking scores |
| `g6_contract/G6_contract_outcome_tierA.json` | §3.4: the corrected gate applied to the released selection |
| `y_randomisation/qsar_y_randomization.json` | §2.4: 100 permutations, z = 21.9, p = 0.0099 |
| `pose_rescue/ifp_rescoring_summary.json` | Table 1: the five pre-registered rescue methods |
| `pose_rescue/ifp_rescoring_per_compound.csv` | Figure 4A and 4C: per-compound scores and labels, n=245 |
| `pose_rescue/multi_reference_ifp_summary.json` | Table 1: the second-reference variants |
| `pose_rescue/pilot_result.json` | Table 1: the water-shell ensemble pilot and its stop rule |
| `pose_rescue/PREREGISTRATION.json` | the decision rule, written down before the pilot ran |
| `released_result/e8_np_library_prioritization.csv` | Tables 2 and 3: the prioritisation this study releases |
| `released_result/tierA_convergence_analysis.json` | §3.5: the count is unremarkable, the identity is not |
| `released_result/tierA_literature_crosscheck.md` | §3.5: the post-hoc PubMed cross-check, hits and non-hits |
| `released_result/conformal_qsar_summary_np_library_dock_ledger.json` | §3.5: conformal validity and the exchangeability caveat |
| `structure_and_protocol/structure_and_protocol_checks.json` | §2.3/§2.5: receptor constructs, redock RMSD, stereocentre counts |
| `structure_and_protocol/box_and_cutoff_checks.json` | §2.5: ligand extent vs box, centroid-cutoff sensitivity |
| `structure_and_protocol/ligand_extents.csv` | §2.5: per-compound principal extents |
| `structure_and_protocol/g8_number_traceability.json` | the G8 audit of every number in the manuscript |
| `orthogonal_engine/orthogonal_check_result.json` | §3.5: the pre-registered DiffDock check -- control gate 1.92 A, candidate disagreement 6.43/8.80/9.50 A, and the pre-registered prediction scored as failed |
| `orthogonal_engine/diffdock_input.csv` | §2.5: the four ligands submitted, as InChIKey-verified SMILES |
| `orthogonal_engine/4EIY_protein_only.pdb` | §2.5: the receptor given to DiffDock, waters removed as declared |
| `heldout_outcome/heldout_evaluation_summary.json` | §3.1: which of the five project-blind fault cases each gate detected |
| `heldout_outcome/fix_verification_summary.json` | §3.1: the regression check confirming the H01-H03 patches -- explicitly NOT a second held-out measurement |
| `md_candidate_builds/NPLIB01185_build_metadata.json` | §2.9: identity check against the frozen ledger and parmchk2 coverage |
| `md_candidate_builds/NPLIB01643_build_metadata.json` | §2.9: identity check against the frozen ledger and parmchk2 coverage |
| `md_candidate_builds/NPLIB03056_build_metadata.json` | §2.9: identity check against the frozen ledger and parmchk2 coverage |
| `md_control/ligand_rmsd.xvg` | §2.9: ZM241385 heavy-atom RMSD over the 100 ns control trajectory |
| `md_control/box.xvg` | §2.9: box dimensions, the bilayer-equilibration check |
| `applicability_domain/ad_percentile_position.csv` | §3.5 and Figure 4B: where each compound sits inside its own space's training self-similarity distribution -- the rank view the admission rates are a threshold on (ADR-AIIM-126) |
| `applicability_domain/ad_sweep_selfsim_2026-08-25_comparison_vs_snplib.json` | Figure 4A: the two libraries side by side, now with a producer script rather than hand-assembled (ADR-AIIM-136) |
| `applicability_domain/run_provenance.json` | §2.1: tooling versions and command line for the sweep, which the 2026-08-21 run did not record |
| `released_result/conformal_diagnostics.json` | §3.6: empty-set rate, prediction-set size and per-class p-value distributions on and off domain, and where the Tier-A criterion sits in the distribution it thresholds (ADR-AIIM-137) |
| `released_result/calibration_pvalues.csv` | §3.6: per-compound cross-conformal p-values on the held-out set, the in-domain comparator for the query libraries |
| `funnel_control/fastpass_correlation.json` | §3.5: agreement between the fast-pass and full-protocol scores, and the residual spread against the selection boundary (ADR-AIIM-137) |
| `funnel_control/fastpass_vs_fullprotocol.csv` | §3.5: the per-compound pairing behind that comparison |
| `qsar_calibration/qsar_calibration.json` | §2.4: Brier score, reliability curve and the reconciliation of the three AUC figures the manuscript quotes (ADR-AIIM-140) |
| `qsar_calibration/reliability_curve.csv` | §2.4: the ten-bin reliability curve, including the two bins the Tier-A probabilities fall in |
| `applicability_domain/ad_descriptor_sweep_2026-08-15_ad_coverage_by_descriptor.csv` | §3.3 and the graphical abstract: per-space coverage for the 236-compound SNP-Lib itself -- 0 % in four structural spaces and 152 of 236 (64.4 %) physicochemically. Distinct from the expanded selection's table beside it, and the file the graphical abstract was corrected against (ADR-AIIM-175) |
| `structure_checks/scores_control.csv` | §3.6: the ZM241385 control arm that gates the ligand-transfer redocking |
| `structure_checks/scores_T1.csv` | §3.6: 9XT redocked into 5OLH, the structure solved with it -- 0.216 A |
| `structure_checks/scores_T2.csv` | §3.6: the same ligand into 4EIY, a structure solved with a different one -- 9.606 A |
| `structure_checks/9XT.sdf` | §3.6: the failing pose itself, so the 9.606 A is a geometry a reader can inspect rather than a number to accept |
| `learned_monitor/ood_fit_summary.json` | §3.2: the isolation-forest applicability-domain detector as fitted, including the contamination setting its saturation behaviour follows from |
| `seeded_defects/wp4_conditions_ac_results.json` | the pre-recount conditions A and C, published so the denominator correction of ADR-AIIM-131 can be checked against its input |
| `released_result/funnel_metadata.json` | §2.6: the funnel's own record of what entered each stage and what the cutoffs were |
| `structure_checks/system_gmx.gro` | §2.9: the assembled membrane system the simulation reads. The atom, lipid, water and restraint counts the Methods quote are counted from this file rather than reported from the build log |
| `released_result/snp_source_ledger.csv` | the 236-compound SNP-Lib accountability ledger as G1 produced it. Deposited 2026-09-13 because 42 tests in the fault-injection suite seed their workspaces from it: without it the suite cannot run from this deposit, which is the only place a reader can run it |
| `incident_taxonomy/incident_taxonomy.md` | §3.6: the fifteen incidents classified by layer, what detected each, how it resolved |

## What is deliberately not here

Generated from `10_new_system/scripts/evidence_not_used.py`, which is where each of these decisions is recorded with its reason and the decision record that made it. Nothing here is summarised.

### Reported, but the material is withheld

**`heldout_evaluation_2026-08-08`** — sealed, excluded from deposit and s10 (ADR-AIIM-113, ADR-AIIM-142).

What the held-out evaluation measured is reported in Section 3.1 and is required for the paper to be honest about its own test. The faulted workdirs are a different thing: they let H01-H05 be reconstructed, and H06-H12 have never been exposed. Publishing them would destroy a blind test that still has seven unused cases in it. An attempt to deposit them on 2026-08-23 was caught and reverted; the exclusion is now a standing comment in build_deposit.py. Unsealing requires separate explicit human approval and is not delegated to any script.

Published anyway: `20_new_runs/heldout_evaluation_2026-08-08/heldout_summary.json`.

**`md_candidates_2026-08-22`** — too large to deposit, excluded from deposit (ADR-AIIM-142).

System build metadata and setup parameters are reported and deposited so the runs can be rebuilt. The trajectories themselves are gigabytes of .xtc and .cpt, above the deposit's size cap, and are covered by SKIP_SUFFIX. This is a size decision, not a disclosure decision: the inputs and the protocol are both published.

**`np_library_funnel_2026-08-15`** — too large to deposit, excluded from s10 and deposit (ADR-AIIM-142).

The most-used run in the project, and its per-compound pose directories run to several gigabytes. Every number the manuscript reports from the funnel is published, along with the ledgers those numbers are computed from; the pose archives themselves are in the repository named in the data-availability statement, not in this bundle. A size decision, not a disclosure one.

### Runs this manuscript does not rest on

Produced during this study and not used. Listed so that the work which did not survive is visible beside the work that did.

| Run | Why |
|---|---|
| `WP4_retrain_lofo_2026-08-05` | Leave-one-family-out retraining of the risk monitor. Run to test whether the fitted model survived removing each fault family in turn; the adopted model came from the postfix run instead. No number in the manuscript comes from here. Superseded by `WP4_retrain_postfix_2026-08-05`. |
| `WP4_retrain_logistic_override_2026-08-05` | Retraining with the logistic model forced in place of the selection procedure's choice. A sensitivity check, not a candidate model; not adopted and not reported. Superseded by `WP4_retrain_postfix_2026-08-05`. |
| `confidence_tier_positive_control_2026-08-08` | Positive control for the v4 confidence-tier machinery. v6 reports conformal prediction sets and the learned monitor's three states instead; the tier vocabulary it validates does not appear in the manuscript. |
| `control_representativeness_validation_2026-08-04` | Leave-one-out validation of the control_representativeness check. The manuscript's position is that this check is mis-specified, so validating it is beside the point v6 makes; no machinery references the directory. |
| `exploratory_sbdd_ranking_2026-08-09` | Source of the 10- and 18-compound candidate lists that the v4 manuscript itself withdrew, in its own words, as something that must not be presented as a prioritisation. Superseded entirely by the E8 output. Retained on disk only so that the withdrawal is auditable rather than a claim about a deleted directory. |
| `manuscript_remediation_2026-08-25` | Session logs for the review-remediation work of August 2026. Holds no scientific artifact and feeds no manuscript number. |
| `md_benchmark_2026-08-21` | GROMACS throughput benchmark. The production run settings were chosen from its ns/day figures, but nothing in production reads a file from it -- the choice is recorded in an ADR, not on disk. It measures this machine, not the system under study. |
| `md_mmgbsa_audit_2026-08-07` | Audit of the MM-GBSA rescoring arm. That arm was dropped from v6 and its Methods and Results sections were cut, so nothing here is cited. Kept because the decision to drop it should remain inspectable. |
| `mmgbsa_pilot_2026-08-10` | The MM-GBSA pilot itself, 501 files. Same decision as the audit above. The largest not-used directory in the project. |
| `reference_check_2026-09-03` | A correctness check on the reference list, not a source of any reported number. It found and fixed a citation that pointed at the wrong paper (ADR-AIIM-165); the fix lives in the manuscript, and nothing in §3 depends on this run. |
| `reference_check_2026-09-05` | The reference audit re-run after four citations were added for the capability decomposition. Like its 2026-09-03 predecessor it is a correctness check on the reference list, not a source of any reported number (ADR-AIIM-169). |
| `snplib_source_audit_2026-09-12` | Asks whether reference [13]'s promised list of 38 primary isolation reports can be assembled from anything this project holds. It cannot: the frozen SNP-Lib spreadsheet carries journal, year and source species but no author, title, volume, page or DOI. Reports no number the manuscript states -- its 36 and 39 groupings are the two reconstructions that bracket the stated 38 and match neither, which is the finding rather than a result (ADR-AIIM-180). |
| `validation_matrix_2026-09-05` | A derived index, not a measurement. Table 7 is built from it, but every cell restates a number that G8 already audits against the artifact that produced it, so nothing in the manuscript depends on this run's own output. Kept because the decomposition it applies -- what each technique was validated FOR -- is the argument of section 3.6, and the mapping from claim to artifact should be inspectable (ADR-AIIM-169). |

## Scripts that produced these artifacts

All under `10_new_system/scripts/`:
`run_g4_verification.sh`, `build_metamorphic_ledgers.py`,
`run_metamorphic_docking.sh`, `analyse_metamorphic.py`,
`ad_descriptor_sweep.py`, `funnel_exclusion_control.py`,
`run_g6_on_np_library.py`, `qsar_y_randomization.py`,
`figure_ad_multispace.py`, `g8_number_traceability.py`,
`diffdock_orthogonal_check.py`, `md_11_analyse_pose_persistence.sh`.

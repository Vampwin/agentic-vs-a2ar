# S10 — Verification evidence

Sections S1–S9 were assembled on 2026-08-08 and predate every run the
verification claims in this manuscript rest on. They are left unmodified as
the record of what accompanied the previous version. This section carries the
2026-08-21 evidence and nothing else.

Every file is a **copy of the artifact the analysis read**, not a regenerated
one, and the manifest records its SHA-256 so a reader can confirm the copy
matches. A supplement that re-derives its own numbers demonstrates nothing
about the run that produced them.

**47 files, 1092 KB.**

| File | Supports |
|---|---|
| `seeded_defects/mutation_score_by_gate.json` | Table 4 and §3.1: per-gate and per-family detection, the 2x2 ablation |
| `seeded_defects/wp4_conditions_ac_results_postfix_2026-08-05.json` | the underlying run the per-gate table is derived from |
| `seeded_defects/naive_baseline_summary.json` | §3.1: the ungoverned comparator, SEER 0.477 |
| `seeded_defects/real_artifact_replay_results.json` | Table 4: 4/4 real-artifact defect replay, and the 0.379 A redock RMSD |
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
| `applicability_domain/ad_coverage_by_descriptor.csv` | Figure 3A and §3.4: coverage in all eight spaces |
| `applicability_domain/ad_per_compound.csv` | Figure 3B/3C: per-compound similarity and in-domain flags |
| `applicability_domain/comparison_vs_snplib.json` | §3.4: coverage before vs after similarity selection |
| `applicability_domain/run_metadata.json` | the single rule applied identically in every space |
| `applicability_domain/figure_ad_multispace_provenance.json` | Figure 3: input and output digests |
| `funnel_control/exclusion_control_summary.json` | §3.4: 7.0% of discarded compounds reach the selection cutoff vs 10.0% kept |
| `funnel_control/control_metadata.json` | how the 500 discarded compounds were sampled |
| `funnel_control/exclusion_scores.csv` | the control docking scores |
| `g6_contract/G6_contract_outcome_tierA.json` | §3.4: the corrected gate applied to the released selection |
| `y_randomisation/qsar_y_randomization.json` | §2.4: 100 permutations, z = 21.9, p = 0.0099 |
| `pose_rescue/ifp_rescoring_summary.json` | Table 1: the five pre-registered rescue methods |
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
| `incident_taxonomy/incident_taxonomy.md` | §3.6: the fifteen incidents classified by layer, what detected each, how it resolved |

## What is deliberately not here

- **The held-out fault mechanism (H01–H12).** The `apply_fault()` code that
  defines each case is sealed; unsealing it requires a separate human
  approval that has not been given. What the evaluation *measured* is
  reported: `heldout_outcome/` carries which case each gate detected, and
  §3.1 states the result. Reporting an outcome and publishing the test that
  produced it are different acts (ADR-AIIM-113/115).
- **Molecular-dynamics trajectories.** The control run's analysis output is
  here under `md_control/`; the trajectories themselves run to gigabytes and
  are available from the corresponding author. The three candidate
  simulations were still running at submission (§4).
- **Full pose archives.** The per-compound pose directories run to several GB.
  They are part of the repository named in the data-availability statement,
  not of this bundle.

## Scripts that produced these artifacts

All under `10_new_system/scripts/`:
`run_g4_verification.sh`, `build_metamorphic_ledgers.py`,
`run_metamorphic_docking.sh`, `analyse_metamorphic.py`,
`ad_descriptor_sweep.py`, `funnel_exclusion_control.py`,
`run_g6_on_np_library.py`, `qsar_y_randomization.py`,
`figure_ad_multispace.py`, `g8_number_traceability.py`,
`diffdock_orthogonal_check.py`, `md_11_analyse_pose_persistence.sh`.

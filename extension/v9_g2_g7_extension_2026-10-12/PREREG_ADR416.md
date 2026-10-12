# Pre-registration (ADR-AIIM-416): seeded faults for gates G2 and G7, which carried no implemented fault family

Written 2026-10-12 after the clean workspaces were built and checked, and before any faulted workspace was evaluated. Rules, anchor,
predictions and grading are fixed here; results are recorded afterwards in ADR-AIIM-417, including predictions that are not met.
Nothing in `10_new_system/` is modified: the frozen evaluators and contracts are imported and read, the engine of the frozen tree is
called as in ADR-AIIM-236, and the injectors are written in the run folder, implementing the parameter grids of the stubs F06, F15 and
F16. The seeded benchmark of the article (124 of 128 under the harness contract, 101 of 128 as deployed) is not changed or pooled
with this run: this is an extension that is reported beside it, for two gates it did not cover. The families F07, F08, F11 and F12
(G3 and G5) are not part of this run.

What the author of this file had seen: the source of the G2 and G7 evaluators and contracts (so the predictions below test a reading of
that code, not an unknown behaviour of it); the prepared receptors, the reference ligand and the box of the project; the first lines of
the real ligand-RMSD series and the run parameter files of the molecular-dynamics runs; the outcome of the clean workspaces under both
paths (`anchor_clean.json`: 21 of 21 passed, the anchor below). Not seen: the evaluator state of any faulted workspace.

## Design
Unit: one workspace (a gate's inputs for one system). A faulted workspace is a clean one with one injected fault; the cases of a family
share their base system, so a rate is read per family, not pooled, and "0 of n" is read with the rule of three (upper 95% bound 3/n),
not as a precise zero. Randomisation: the grids of the stubs are fixed; the random choice of F06 uses seed 7 for every instance.
Controls: the clean workspaces below, evaluated by the same two paths.

Paths. B: the frozen evaluator called directly with the frozen YAML contract (arm B of ADR-AIIM-236). C: the deployed path of the frozen
tree, orchestrator, contract engine and report lock (arm C of ADR-AIIM-236). A case is detected when the state is not PASS; under C it
is contained when the workflow halted or the artifact is not reportable. There is no harness contract (arm A) for these two gates.

Clean workspaces (the anchor). G2: ten distinct prepared 4EIY receptors of this project (the registered run and its two reproductions,
the 3, 4 and 5 A water shells, and four copies made for other analyses; one has no water), each with the real reference ligand and box,
the water-retention decision that the file bears, the cofactor decisions and consistent embedding counts (223 of 223). The registry of the
production G2 of the project lists four passed checks, while the contract declares nine; the clean workspaces supply the inputs of
all nine, as a compliant run would. G7: eleven systems built from the real ligand-RMSD series and run parameter files: the five simulated
compounds and the control with all replicates (6), and the five compounds with their first two replicates (5); planned production
50 ns (the common window of the analysis); the claim boundary; the real shortlisted panel; the seed read from each run's parameter
file (the control's first run from the velocity-generating step it continued).
Anchor: all 21 clean workspaces pass under both paths. If not, the faulted cases are not read.

## Faults (nine, ten and ten instances; the grids of the stubs)
- F06 WRONG_PROTONATION (G2), 9 instances (`n_his_flipped` 1 to 3 times `residue_selection` first, random, all): the ring hydrogen of
  the selected histidines of the prepared receptor is moved from ND1 to NE2 (HID to HIE) and the two nitrogens swap AutoDock type;
  line count, serial numbers and residue keys do not change. Every histidine of these receptors (eight) is HID. The three "all"
  instances are the same operation (all eight flipped), so the family holds seven distinct operations. The recorded water, cofactor and
  embedding decisions stay as the clean run wrote them.
- F15 TRUNCATED_TRAJ (G7), 10 instances (`fraction_kept` 0.05, 0.10, 0.15, 0.20, 0.25, 0.33, 0.40, 0.50, 0.67, 0.80 of the planned 50 ns):
  the per-frame series read by the equilibration check is cut to that fraction before analysis; the completion log stays as the run
  wrote it (silent truncation, as in the other families). F15L is the same truncation with the completion log recording the shorter
  production length (10 further instances); it is a control on the evaluator's completion check and is reported apart.
- F16 RESIDUE_INDEX_ERROR (G7), 10 instances (`residue_offset` -50, -30, -20, -10, -5, 5, 10, 20, 30, 50): every residue index that
  defines the binding site (the 52 receptor C-alpha atoms within 10 A of the crystal ligand) in the results table is shifted by the offset.
Vacuity: an injection that changes no input is excluded from the rate (none is expected).

## Predictions
P1. F06: 0 of 9 detected under B and 0 of 9 under C. (The G2 evaluator checks residue counts, chain breaks, the reference ligand, the box
fields, the retained-water count, the cofactor decisions and the embedding counts; nothing in it reads a protonation state.)
P2. F15 (silent): 0 of 10 detected under both paths. (The completion check reads the log; the equilibration check only requires the
series to exist.)
P3. F15L (log updated): 10 of 10 detected under both paths, by `trajectory_completion`.
P4. F16: 0 of 10 detected under both paths. (The evaluator never reads the results table.)
P5. For every case, clean or faulted, the state under C equals the state under B.
P6. Under C, every detected case is contained and every undetected case is not.
P7. Anchor: the 21 clean workspaces pass under both paths (observed before this file was written).

## Reporting
All predictions are graded and reported, met or not, with the checks that fire. The extension enters the supplement beside Section 15 and
the article's Table 2 note is amended: G2 and G7 are no longer without a fault family in the extension, and the detection of the
implemented extension families is stated apart from the 124 of 128 and 101 of 128. A family with no detection is reported as a gap of
the contract (what the check list does not look at), not as a failure of the benchmark. No frozen number changes; the held-out
cases are not used; no compound is called active.

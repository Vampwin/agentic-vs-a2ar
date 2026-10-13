# Pre-registration (ADR-AIIM-427): four added checks for the faults the evaluators did not detect (G2, G7, G5)

Written 2026-10-13 after the clean workspaces were built and evaluated, and before any faulted workspace was built or evaluated. Rules,
anchor, predictions and grading are fixed here; results are recorded afterwards in ADR-AIIM-428, including predictions that are not met.
Nothing in `10_new_system/` is modified (a snapshot of the 242 source, contract and configuration files there is compared before and after
the run). The frozen evaluators are imported and called unchanged by wrapper evaluators in this folder (`analysis/v9ac/`), which append four
checks to their output; the frozen contracts are copied to `contracts_plus/` with the four checks declared and the wrapper named as the
evaluator, so the deployed path (orchestrator, contract engine, report lock) can load them. Integrating the checks into the frozen tree is
a separate act that needs the author's approval and is not part of this run. The seeded benchmark of the article (124 of 128 under the
harness contract, 101 of 128 as deployed) and the extensions of ADR-AIIM-416 to 422 are not changed or pooled with this run.

What the author of this file had seen: the results of ADR-AIIM-417 and 421, so the families that the evaluators missed (F06, F15 silent, F16,
F11 silent) and why; the source of the G2, G7 and G5a evaluators; the clean workspaces of this run (`anchor_clean.json`). Not seen: the state
of any faulted workspace under the added checks. The checks were designed after seeing the faults of ADR-AIIM-416 and 420, so those faults
(set G) are reported as seen and are not evidence of detection of unseen faults. Sets A to F are new instances: other base receptors and other
base systems, other seeds, other severities, other benchmark rankings and other leaked IDs. They are held out from the design of the
checks in that sense only: the fault families are the same, so they test whether the checks generalise across instances, not whether they
detect fault families nobody had thought of.

## The four checks (all fail closed; severity critical, so a failure halts)
- G2 `protonation_matches_record`: the receptor carries, for every histidine, the ring hydrogen (HD1, HE2 or both) that the preparation step
  recorded in `protonation_decision`; a histidine absent from the record, a record naming a histidine absent from the file, and a missing
  record all fail. It compares the file with the record; it does not judge whether the recorded state is the right one chemically.
- G7 `trajectory_covers_plan`: every per-frame series (compound and replicate) reaches at least 0.95 of `planned_production_ns` (47.5 of 50 ns).
  The tolerance of 0.95 is fixed here, before any truncated series is seen; the clean series reach at least 49.9 ns.
- G7 `binding_site_matches_receptor`: the residue indices that define the binding site in the results table equal, for every row, the set of
  receptor C-alpha residues within 10 A of the reference ligand, recomputed from the prepared receptor and the reference ligand files.
- G5a `benchmark_ids_not_exposed`: no file in the access log, read before the gate by a step other than G5a and G5b and not listed in
  `registered_ranking_inputs`, contains a compound ID of the benchmark (any ID, active or inactive) or cannot be opened. A file that is not in
  the access log is not seen by this check (a limit, tested as set F).

## Design
Unit: one workspace. Rates are read per family and set; "0 of n" with the rule of three. Vacuity: an injection that changes no input is excluded.
Paths: F, the frozen evaluator and frozen YAML contract called directly; P, the same with the added checks, called directly; D, the deployed path
on the copied contracts. A case is detected when the gate state is not PASS (for G5, the state of G5a); it is attributed to an added check when
that check is in the failed list. The clean cases are the anchor: 37 workspaces (ten G2 receptors, eleven G7 systems, and sixteen G5 rankings of
the benchmark, five of them with a context file of 40 IDs that are not benchmark compounds read by the ranking step), all PASS in F, P and D
(observed before this file was written; the frozen tree unchanged, 242 files).

Observation made while building the anchor, outside the gate checks: the histidine atom types of two prepared receptors of other analyses of
this project are inconsistent with the hydrogens they carry (the receptor of the cross-docking set holds HIS 333 with HD1 and ND1 typed
acceptor; the cross-docking validation receptor holds eight histidines without any ring hydrogen). This is reported as an observation; it is
not an added gate check and its effect on any result was not evaluated.

Sets (seed 11 selects histidines, seed 13 selects leaked IDs, seed 7 is that of ADR-AIIM-416 and 420):
- A, F06 held out: the nine receptors of ADR-AIIM-416 other than the registered run's, each with 1 and with 3 histidines flipped at random, as
  in ADR-AIIM-416 (hydrogen moved ND1 to NE2 or back, atom types swapped); 18 instances. The cross-docking validation receptor carries no ring
  hydrogens, so its two instances are vacuous.
- B, F15 held out (silent truncation, completion log unchanged): SNP00163, NPLIB01185 and ZM241385, all replicates, with 0.02, 0.30, 0.55,
  0.75, 0.90, 0.94, 0.96 and 0.99 of the 50 ns kept; 24 instances.
- C, F16 held out: the same three systems with the residue offset -100, -3, -2, -1, 1, 2, 3 and 100; 24 instances.
- D, F11 held out (logged leak file): five rankings of the benchmark (`max_ifp_tanimoto`, `asn253_filtered`, `combined_zscore`, `ifp_zm_only`,
  `ifp_max_of_two`) with 1, 3, 10 and 30 public active IDs given to the ranking step and ranked first; 20 instances.
- E, decoy IDs: the same five rankings with 3 and 30 inactive IDs read by the ranking step and ranked first; 10 instances (the check is strict
  about IDs, not about labels; this set measures that strictness).
- F, unlogged leak: the same five rankings with 10 active IDs leaked through a file that the access log does not name; 5 instances.
- G, seen: the faults of ADR-AIIM-416 (F06 9, F15 10, F15L 10, F16 10, on the original bases) and of ADR-AIIM-420 (F11 10 and F11L 10, Vina
  ranking), re-evaluated under F, P and D.

## Predictions
P1. Anchor: all 37 clean cases PASS in F, P and D.
P2. Set A: 16 non-vacuous instances (2 vacuous); F detects 0, P and D detect 16 of 16, all by `protonation_matches_record` alone.
P3. Set B: for the fractions 0.02 to 0.94 (18 instances) F detects 0 and P and D detect 18 of 18, by `trajectory_covers_plan` alone; for 0.96 and 0.99
  (6 instances) none is detected in F, P or D (the truncation is within the registered tolerance).
P4. Set C: F detects 0 of 24; P and D detect 24 of 24, by `binding_site_matches_receptor` alone.
P5. Set D: F detects 0 of 20; P and D detect 20 of 20, by `benchmark_ids_not_exposed` alone.
P6. Set E: F detects 0 of 10; P and D detect 10 of 10.
P7. Set F: none of the 5 is detected in F, P or D.
P8. Set G: every one of F06 (9), F15 (10), F15L (10), F16 (10) and F11 (10) is attributed to its added check in P and D (the fractions of F15L are
  all below 0.95); F11L (10) is detected in P and D and attributed to `benchmark_ids_not_exposed` as well as `label_isolation`.
P9. For every case, the state in P equals the state in D.
P10. Every case detected in D by an added check is halted.
P11. The 242 files of the frozen tree are byte-identical before and after the run.

## Reporting
All predictions are graded and reported, met or not, with the checks that fire. A held-out rate is not read as an unseen-family detection rate;
the families are the ones the checks were written for. The false-escalation count on the 37 clean cases is reported. The result is an
extension that is reported beside the seeded benchmark: the 124 of 128 and 101 of 128 are unchanged, and no contract or evaluator of the frozen
tree is changed. The sealed held-out cases are not used; no compound is called active.

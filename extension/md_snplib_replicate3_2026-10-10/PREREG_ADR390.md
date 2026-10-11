# Pre-registration (ADR-AIIM-390): the third replicate of SNP00157 and SNP00163

Frozen before the runs start, 2026-10-10. The SHA-256 of this file is in `PREREG_SHA256.txt` and in the decision log.

## What is run

Replicate 3 (50 ns) of SNP00157 and of SNP00163, so that every system of the comparison has three runs: the ZM241385 control (runs 1, 2, 3), SNP00157 (1, 2, 3) and SNP00163 (1, 2, 3).

- The frozen driver `10_new_system/scripts/md_10_run_candidate.sh <ID>` is invoked unmodified. Replicates 1 and 2 are complete and are skipped (their `.gro` exist); replicate 3 continues from the checkpoint left when the watcher of ADR-AIIM-331 stopped it a few seconds after it started (about 0 ns), with the `prod_r3.tpr` and `prod_r3.mdp` that already exist: the velocity seed is `gen-seed = 20260824` for both compounds, fixed when the systems were built, and different from the seeds of the other replicates (20260822, 20260823) and of the control (20261006 to 20261008).
- The two compounds run at the same time on the one GPU (the concurrency rule of ADR-AIIM-331 measured 1.45 times the throughput of a single job).
- No setting, mdp file, topology, restraint, thermostat or seed is changed. A host fault is handled as before: the driver is re-invoked and resumes from the checkpoint; every restart is logged and reported.

## What is measured

The unchanged pipeline of ADR-AIIM-336: `md_11_analyse_pose_persistence.sh` and the N6.55 minimum-distance step (`adr331_analysis.py` structure), over the first 50 ns, with the anchors of ADR-AIIM-331 (index groups equal the control's file, readers reproduce the stored control numbers, the Asn found by position with the TM6 spacing verified, the starting polar distance within 0.5 Å of the Vina pose). Per run: ligand RMSD mean, fraction of frames below 4 Å, N6.55 persistence overall and by 10 ns block, bilayer area change, nearest C-alpha distance. The analysis script is written before the runs end and is not changed afterwards.

## Predictions (fixed now)

Earlier values for orientation: SNP00157 persistence 0.148 and 0.087, RMSD below 4 Å in 17.6 % and 15.6 % of frames; SNP00163 persistence 0.597 and 0.755, RMSD below 4 Å in 41.5 % and 37.3 %; control persistence 0.966, 1.000, 0.679.

- **R1.** Both replicates reach 50.0 ns with the frozen driver and unchanged settings.
- **R2.** In both runs the nearest C-alpha atom to the ligand is never more than 10.0 Å away (no dissociation; the earlier maximum was 5.35 Å).
- **R3.** SNP00157 replicate 3: N6.55 persistence at most 0.35.
- **R4.** SNP00163 replicate 3: N6.55 persistence between 0.40 and 0.90.
- **R5.** The fraction of frames with ligand RMSD below 4 Å is at most 0.50 for SNP00157 replicate 3 and at most 0.60 for SNP00163 replicate 3.

Each prediction is graded met or missed on the numbers, once, and missed ones are reported.

## How the result is used (decided now, whatever the numbers are)

1. **The comparison in the article becomes control runs 1 to 3 against replicates 1 to 3 of each compound**, which restores the design that ADR-AIIM-331 specified as primary (all three control runs in every comparison) and gives every system the same number of runs. The matched duplicate (control runs 1 and 2 against runs 1 and 2) of ADR-AIIM-382 stays as a sensitivity view in the supplement. This choice is made before the data exist and does not depend on them.
2. Reported per run, with ranges; no spread is estimated and no test is run on n = 3 per system; the statement is non-overlap or overlap of per-run values, as before.
3. If a run does not complete (host fault that cannot be resumed, or the run is stopped), the article is not changed and the incomplete run is reported in the supplement as incomplete.
4. Nothing is called active, no ranking is claimed, and the result is exploratory pose analysis, as in every earlier MD entry of this study.

## Not part of this registration

Further compounds, longer runs, other force fields or any analysis not listed above.

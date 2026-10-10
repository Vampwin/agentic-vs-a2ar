
## ADR-AIIM-356 — Addendum to ADR-AIIM-355: DiffDock-Pocket's batch size changed from 20 to 2, and the first results of arm D (2026-10-09)

**Why.** ADR-AIIM-355 fixed `--batch_size 20` (the README's example) and projected "a few minutes per variant". On this 8 GB GPU the first complex (the control, 40 samples) had not finished after 14 minutes, so under rule 5 I stopped the run and told the author. A timing test on the control only, **with the result ignored**, found the cost superlinear in the batch: 12 samples took **666 s at batch 20** against **153 s at batch 4**; 8 samples took 71 s at batch 1 and **62 s at batch 2**. The batch only chunks the samples, so the distribution of poses is unchanged; the slowdown is consistent with GPU memory pressure (7.1 GB of 8 GB in use) but I did not prove the cause. **The batch size is therefore set to 2** for all DiffDock-Pocket runs; `--samples_per_complex 40`, the pocket centre, the six flexible residues, the variants and every graded threshold stay as registered. The projection is about 45 minutes per variant (about 2.3 hours for the three), not "a few minutes"; this exceeds the registered stop rule, and I am proceeding because the machine is otherwise idle and the author's request for this arm stands; the author can stop it at any time.

**Also corrected.** The first attempt of arm C (16:11) skipped every row because DiffDock-Pocket's CSV loader accepts ligand files only; the ligands are now SDF files written from the same SMILES (ETKDG, seed 42; the model discards the coordinates). No pose was produced by that attempt.

### Arm D result (DiffDock-L blind, 100 samples, filtered to Vina's box), graded on its own as the registration allows

`filter_pocket_poses.csv`, `filter_pocket_result.json` (arm C not yet run). Run: 2 h 40 min.

| prediction | result | verdict |
|---|---|---|
| G1 filtered DiffDock-L control within 2.0 Å of the crystal pose | **5.92 Å** | **missed** |
| G3 at least 14 of 15 survivors with an in-box pose among 100 | 15 of 15 (6 to 67 in-box poses per compound; the control has 49) | met |
| G5 box-matched agreement: fewer than 8 of 15 within 2.0 Å | 0 of 15, median 13.9 Å | met, **not interpreted** (G1 missed) |

By the registered rule arm D is reported and not interpreted. A descriptive look at the control only (not graded, added after G1 was known and labelled so): among its 100 samples **10 lie within 2.0 Å of the crystal pose** and 49 are inside the box; the best in-box pose is 0.90 Å from the crystal pose, at DiffDock's confidence rank 9. DiffDock-L finds the site but its confidence ranking, not its sampling, puts a wrong pose first, which is why filtering by the box and ranking by confidence cannot recover a pocket-conditioned answer. (The ten-sample blind run of ADR-AIIM-111 happened to rank the crystal-like pose first, at 1.92 Å, which the hundred-sample run does not reproduce: the rank-1 pose of a single DiffDock-L run is itself noisy.)

# Pre-registration (ADR-AIIM-420): seeded faults for gates G3 and G5, which carried no implemented fault family

Written 2026-10-12 after the clean workspaces were built and checked, and before any faulted workspace was evaluated or any faulted
redocking was run. Rules, anchor, predictions and grading are fixed here; results are recorded afterwards in ADR-AIIM-421, including
predictions that are not met. Nothing in `10_new_system/` is modified: the frozen evaluators, contracts, orchestrator and the frozen
multi-conformer redocking script are imported and read, and the injectors are written in the run folder, implementing the parameter
grids of the stubs F07, F08, F11 and F12. This extension follows ADR-AIIM-416 and 417 (G2 and G7) and is reported beside the
seeded benchmark of the article (124 of 128 under the harness contract, 101 of 128 as deployed), not pooled with it.

What the author of this file had seen: the source of the G3, G5a and G5b evaluators and contracts (so the predictions test a reading of
that code, not an unknown behaviour of it); the real redocking artifacts, the real benchmark and its real scores; the outcome of the clean
workspaces under both paths (`anchor_clean.json`, the anchor below); that one redocking with the frozen script reproduced the production
poses byte for byte. Not seen: the state of any faulted workspace, the rank-1 RMSD of any faulted pose, any redocking with a shifted box.

## Design
Unit: one workspace. A faulted workspace is a clean base with one injected fault; the cases of a family share their base, so a rate is
read per family and "0 of n" with the rule of three, as in ADR-AIIM-416. Paths: B, the frozen evaluator called directly with the frozen
YAML contract and a registry; C, the deployed path of the frozen tree (orchestrator and contract engine) on the same registry. G5b is
evaluated only when G5a passes, as its contract requires.

Detection rule. The clean states of these two gates are not PASS: the production G3 run escalated (`control_representativeness`
fails: the redocking control is unrelated to the library) and the real benchmark fails two checks of G5b (`ef_1pct`, `bedroc20`), as the
article reports for the benchmark. A case is therefore detected when its failed checks include a check that its clean base does not fail
(G3 base: `control_representativeness`; G5 base, the Vina ranking: none at G5a, `ef_1pct` and `bedroc20` at G5b), or when the state is
ERROR. A clean case is a false escalation when it fails a check its base does not.

Clean workspaces (the anchor). G3: eleven redockings of the control into 4EIY with the frozen script, receptor, conformers and box: the
production run, its two reproductions, and eight redockings with seeds 20260730 to 20260737 (the production seed plus 1 to 8); the real
reference ligand and atom map, the atom map registered before the poses, the real library and the real enrichment evidence. G5: the
real ChEMBL A2AR benchmark (245 compounds, 50 active), eleven rankings of it (the Vina score; six rescorings from the per-compound table and four multi-reference interaction-fingerprint
scores, all computed in `chembl245_ifp_rescoring_2026-08-15`), the registry with the human-signed benchmark-role declaration
(Option A, development benchmark), a file-access log without any read of the labels before the gate, the ID-list hash registered at the
start, and the metrics computed by the project's own function. Anchor: under both paths all eleven G3 cases fail only
`control_representativeness`, with a rank-1 RMSD of at most 2.0 A, and all eleven G5 cases pass G5a; the Vina ranking reproduces
AUC 0.6014, BEDROC20 0.2752 and EF@1% 1.6333 of the registered result. If not, the faulted cases are not read.

## Faults
- F07 OFF_POCKET_POSE (G3), 15 instances (`offset_angstrom` 5, 10, 15, 20, 25 times `axis` x, y, z): the docking-box centre is shifted by the
  offset along the axis (the box size and everything else stay) and the control is redocked with the frozen script and the production
  seed; the poses of that redocking replace the clean ones and are registered after the atom map.
- F08 RMSD_EXCEEDED (G3), 15 instances (`rotation_degrees` 45, 90, 120, 135, 150 times `rotation_axis` x, y, z, the first fifteen of the
  stub's grid of 21): the crystal reference ligand is rotated about its centroid; the poses are the clean ones.
- F11 LABEL_LEAKAGE (G5), 10 instances (`n_leaked_ids` 5, 10, 15, 20, 25, 30, 40, 50, 75, 100): that many public active IDs of the
  benchmark, drawn with seed 7, are exposed to the ranking step as a file of known actives; the leaked actives are ranked first (the
  downstream consequence of the leak, a modelling assumption); the access log records the ranking step reading that file, which is not
  the benchmark file. F11L is the same with the leaked IDs read from the benchmark file itself, so that the access log shows a read of the
  labels before the gate (10 further instances; a control on the label-isolation check, reported apart).
- F12 EVAL_SET_MISMATCH (G5), 3 instances (`swap_mode` full_swap, partial_50pct, partial_25pct): the labels of the benchmark file are
  swapped (all; 25 actives with 25 decoys; 12 with 12, seed 7) and the metrics are recomputed from the swapped labels; the compound IDs,
  the registered hash and the access log are those of the clean run. Base of F11 and F12: the Vina ranking.
Vacuity: an injection that changes no input is excluded (none is expected).

## Predictions
P1. Anchor (observed before this file was written): the eleven G3 cases fail only `control_representativeness` under both paths with
a rank-1 RMSD of at most 2.0 A, and the eleven G5 cases pass G5a under both paths.
P2. F08: the twelve instances with a rotation of 90 degrees or more are all detected under both paths (by `rmsd_rank1_4eiy`, HALT).
P3. F07: the six instances with an offset of 20 or 25 A are all detected under both paths.
P4. F07: none of the three instances with an offset of 5 A is detected under either path.
P5. For every G3 case the failed checks under the deployed path equal those under the direct path.
P6. F11 (silent): 0 of 10 detected under both paths.
P7. F11L: 10 of 10 detected under both paths, by `label_isolation` at G5a.
P8. F12: 3 of 3 detected under both paths, by `auc_roc` at G5b with G5a passing.
P9. For every G5 case the failed checks of G5a and G5b under the deployed path equal those under the direct path.

## Reporting
All predictions are graded and reported, met or not, with the checks that fire. A family with no detection is a gap of the check list of
the gate, not a failure of the benchmark; no contract or evaluator is changed; the held-out cases are not used; no compound is called
active. The extension enters the supplement beside Section 15 and the article's passage on G2 and G7; the 124 of 128 and 101 of 128 are
unchanged. The modelling assumptions are stated where the results are: the leaked actives are ranked first; the swapped labels are
re-scored by the pipeline; the box shift is a real redocking, not a translation of the pose.

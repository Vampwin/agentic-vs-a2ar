# Pre-registration (ADR-AIIM-431): the four added checks placed in the frozen tree as new files, and an equivalence check

Written 2026-10-13 after the new files were written and their unit tests passed, and before the equivalence run was made. The author approved,
in chat, adding new files to `10_new_system/` without editing any existing file ("ข้อ 1 แบบเพิ่มไฟล์ใหม่"). The files added are: `evaluators/added_checks.py`
(the four checks and the merge of their result into an evaluator's output, the functions of `v9_added_checks_2026-10-13/analysis/v9ac/checks.py`);
`evaluators/G2_evaluator_added_checks.py`, `G7_evaluator_added_checks.py` and `G5a_evaluator_added_checks.py` (each calls the frozen evaluator
unchanged and appends the check or checks); `contracts/G2_preparation_added_checks.yaml`, `G7_md_added_checks.yaml` and
`G5a_validity_added_checks.yaml` (the frozen contract text with a new contract identifier, the new evaluator, the check entries and a version
suffix); and `tests/test_added_checks.py` (eight unit tests). No existing file is opened for writing. The frozen contracts and evaluators
remain those that the seeded benchmark of the article measured; the new contracts are not part of that benchmark.

What the author of this file had seen: the results of ADR-AIIM-428, which this run reproduces from the files that now stand in the tree. Not
seen: the output of the in-tree contracts or evaluators on any workspace.

## Design
The run `v9_added_checks_2026-10-13` is repeated unchanged in everything but the contract directory, the contract file names and the evaluator
modules: the 37 clean and the 160 faulted workspaces are built by the same code, with the same seeds, and evaluated (a) by the frozen evaluators
directly, (b) by the in-tree `*_added_checks` evaluators called directly with the in-tree `*_added_checks` contracts, and (c) deployed (orchestrator,
contract engine, report lock) on the in-tree contracts directory. The grading of ADR-AIIM-427 is applied again. Nothing is written into the tree by
the run; a hash of every file of the tree is compared before and after the files were added.

## Predictions
P1. Of the 271 files that stood in the tree before the new files were written (hash recorded in `tree_before.json`), none has changed; the files
added are the seven above and the test file.
P2. The test suite of the tree passes in full: the 715 existing tests and the eight new ones.
P3. All 37 clean cases PASS in (a), (b) and (c).
P4. For each of the 160 faulted cases the state and the list of failed checks in (a), (b) and (c) equal those recorded in ADR-AIIM-428.
P5. The eleven predictions of ADR-AIIM-427 are met again.

## Reporting
Results are recorded in ADR-AIIM-432, met or not. A difference between this run and the earlier one is reported as such and is not smoothed over.
The 124 of 128 and 101 of 128 of the article are not recomputed by this run: they are the benchmark of the frozen contracts, which are unchanged.
No compound is called active; the sealed held-out cases are not used.

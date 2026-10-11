# Pre-registration (ADR-AIIM-399): six short sensitivity analyses asked for by the internal review of v9

Written 2026-10-11 before any of the analyses below was run. The rules, anchors and predictions are fixed here; the results are
recorded afterwards in ADR-AIIM-400, including predictions that are not met. Nothing in `10_new_system/` is modified; frozen
scripts are imported or invoked unchanged; earlier run folders are read, not written. No compound is called active and no
ranking is changed by any outcome: the released selection, the fifteen and the eight stay what they are, and these analyses
are reported beside them as sensitivity views.

What the author of this file had already seen when writing it (so that the predictions can be weighed): the header and the
first data row of `snplib_dock_2026-10-06/pool69_predictions.csv`; `em_survivors_2026-10-06/em_survivors.json` (the fifteen,
the eight, the five distinct split groups and the number of distinct splits behind each of the fifteen); the summary of
`snplib_dedup_2026-10-05/snplib_dedup.json` (the top-20 lists of an EARLIER run of the eight conventions, superseded by the
one-thread run); `verification_summary_2026-08-21/ligand_extents.csv` rows with L1 > 22 and their stored scores; the counts of
Supplementary Section 11. Not seen: the applicability scores of the other 68 pool compounds, any score of the expansion
library against the 6,284-compound reference other than the count 2,645 at 0.25, any docking result in a larger box, the
contract-to-evaluator map, any applicability sweep at another k or percentile.

## A1. One vote per distinct split (review point P0-09)
Input: `20_new_runs/snplib_dock_2026-10-06/pool69_predictions.csv` (predicted pKi of the 69 pool compounds under each of the
eight conventions, one-thread fits) and the five distinct split groups of `em_survivors.json`.
Rule: unchanged from the registered selection except for the vote. Per split, the top 20 of the pool by predicted pKi; a
compound is selected when it is in the top 20 of at least 3 of the 5 distinct splits (a duplicated convention counts once);
order by (number of distinct splits, mean predicted pKi over the eight conventions), descending; scaffold uniqueness and the
cap of 10 as before.
Anchor (must hold before anything else is read): with the registered vote (at least 5 of 8 conventions) the file reproduces
the fifteen survivors and the eight shortlisted compounds in the registered order.
Predictions. P1: the anchor holds. P2: every one of the fifteen is selected by the 3-of-5 vote and at most three further
compounds are added. P3: at least six of the eight shortlisted compounds are in the shortlist of the 3-of-5 vote.
Reported: the selected set, additions and removals, Jaccard index against the fifteen and against the eight.

## A2. Sensitivity of the focus-library selection to the applicability threshold (P0-08)
Input: as A1, plus the modelling set (6,284 compounds) rebuilt by the registered code path of `pool_comparison.py`
(`rebuild_table`, InChIKey order) and the 5,000-compound expansion library (`np_library_funnel_2026-08-15/np_library_dock_ledger.csv`).
Score: mean ECFP4 (radius 2, 2,048 bits) Tanimoto similarity to the five nearest reference compounds, the function of `pool_comparison.gate`.
Thresholds: 0.25 (the registered gate), 0.30, 0.35 and t*, the quantile-matched value: the score reached by exactly 1,560 of
the 5,000 expansion compounds against the 6,284-compound reference, 1,560 (31.2%) being the number that reached 0.25 against
the 4,566-compound training side on which the competence band was measured. t* is a matching of admission rates, not a
recalibration of the band: no held-out performance is measured against the larger reference.
For each threshold: the pool (those of the 69 at or above it), the selection by the registered vote and by the vote of A1,
and how many of the fifteen and of the eight are still in the pool and still selected. Predictions do not depend on the pool,
so nothing is refitted.
Anchors: 2,645 of 5,000 expansion compounds score at least 0.25; the gated SNP-Lib pool has 69 compounds whose scores equal
those of the stored file to 1e-6.
Predictions. P4: both anchors hold. P5: t* lies between 0.28 and 0.33. P6: at t* the pool holds between 15 and 40 of the 69.
P7: at most 7 of the fifteen are in the pool at t*.

## A3. Interval for the estimate of 1,251 discarded compounds reaching the fast-pass cut-off (P1-19)
Arithmetic on reported counts, no prediction: strata of 13,882 and 3,000 discarded compounds, samples of 249 and 250, the
observed hits read from the artifact of the funnel-exclusion control. Stratified estimator of the total with the
finite-population correction, variance sum over strata of N_h^2 (1 - n_h/N_h) p_h (1 - p_h) / (n_h - 1), normal 95% interval.
Anchor: the point estimate reproduces 7.41% and 1,251.

## A4. The box-constrained compounds in a larger box (P1-17)
The five expansion compounds whose longest principal extent exceeds the 22 Å box edge are NPLIB01997, NPLIB03133,
NPLIB04561, NPLIB03343 and NPLIB04805. The first four completed the full protocol; NPLIB04805 did not and is not docked here.
Run: the frozen `10_new_system/scripts/g4_dock_library.py`, unmodified, at its defaults (eight conformers, exhaustiveness 32,
ten poses, seed 20260729), the receptor and box centre of the production runs, once with the production box (22 Å) and once
with a cubic box of 30 Å (the longest extent, 24.3 Å, plus a margin). Nothing else differs.
Anchor: the 22 Å run reproduces the stored full-protocol scores of the four compounds to 0.001 kcal/mol.
Predictions. P8: the anchor holds. P9: in the 30 Å box at least two of the four score better (more negative) by more than
0.3 kcal/mol. P10: all four keep a rank-1 pose within 8.0 Å of the ZM241385 centroid.
Stop rule: if the docking environment cannot be started, or the two runs together have not finished after 30 minutes, the run
is stopped and reported as not completed; a partial result is not interpreted.

## A5. A static map from declared checks to evaluator source (P0-05)
For every check id declared in the ten contract files of `10_new_system/contracts/`, whether the id occurs as a string literal
in the source of the evaluator module the contract names (and, separately, anywhere under `10_new_system/evaluators/` or
`kdss_audit/`), and whether the seeded benchmark exercised it (`deployed_path_benchmark_2026-09-23/check_coverage.json`,
which covers the five benchmarked gates). This is a string-level inventory of the present tree, which includes the repairs
made after the held-out evaluation; it is not a test that a check is implemented correctly, and an id built at run time
would be missed.
Predictions. P11: the contracts declare 94 checks in total. P12: at most three declared ids are absent from the source of
their own evaluator.

## A6. The applicability rule at other k and percentiles (P2-03)
The frozen `10_new_system/scripts/ad_descriptor_sweep.py` imported unchanged; its module constants K and PCTL are set by the
driver to k = 3, 5, 10 and percentile = 5, 10, 20 (nine settings), SNP-Lib (236 entries) against the same training set, the
same eight spaces.
Anchor: k = 5 and percentile 10 reproduce the stored admissions (0 of 236 in ECFP4, 152 of 236 in the physicochemical space).
Predictions. P13: the anchor holds. P14: in every one of the nine settings ECFP4, FCFP4, topological-torsion and Gobbi
pharmacophore spaces each admit at most 5 of 236. P15: the physicochemical space admits between 50% and 80% of the 236 in
every setting.
Stop rule: if the nine settings have not finished after 40 minutes the run is stopped and the settings completed are reported
as such.

## Use of the results, decided now
Whatever the numbers, the article keeps its registered selection and reports these as sensitivity analyses: A1 and A2 in
Section 3.8.1 and Supplementary Section 17, A3 in Section 3.5, A4 in Section 3.2, A5 in Supplementary Section 15, A6 in
Section 3.3 or Supplementary Section 7. A changed list of compounds is reported as a changed list. No compound is put forward.

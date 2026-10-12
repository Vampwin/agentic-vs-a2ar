# Pre-registration (ADR-AIIM-412): the censored-target fit under four checks (internal review point P1-20)

Written 2026-10-12 before any of the analyses below was run. Rules, anchors, predictions and grading are fixed here; results are
recorded afterwards in ADR-AIIM-413, including predictions that are not met. Nothing in `10_new_system/` is modified; frozen
scripts and earlier run folders are imported or read, never written. No compound is called active and no ranking is changed by
any outcome: the released selection (the fifteen and the eight) stays what it is, and these analyses are reported beside it as
sensitivity views. A censored-likelihood model is added as a comparison, not as a replacement.

The review point (P1-20): "iterative projection of censored pKi values lacks methodological justification and amplifies numerical
noise; the ten-pass feedback loop is not standard censored regression; membership changes with floating-point differences. Remedy:
use a method for censored outcomes that has evidence, or explain the algorithm with simulation validation, a convergence criterion,
one-thread determinism and sensitivity to initialisation."

What the author of this file had already seen: `pki_censored_em_2026-10-05/pki_censored_em.json` (every one of the eight conventions
ran all 10 passes; final maximum change 0.026 to 0.052 pKi against the tolerance 0.01, so no convention stopped by tolerance);
`em_survivors_2026-10-06/em_survivors.json` (the fifteen, the eight, five distinct split groups, wall time 2,218 s);
`em_reproducibility_2026-10-06/em_reproducibility.json` (one-thread fits bit-identical, three-thread fits not); the first rows of
`snplib_dock_2026-10-06/pool69_predictions.csv`. Not seen: any pass beyond the tenth, any other initialisation, any censored-likelihood
fit, any simulation of the algorithm.

## Common set-up
Modelling set, features, scaffold splits, eight conventions, the 69-compound pool and the vote rule are those of
`em_survivors.py` (`rebuild_table`, `qsar.compute_raw_features`, `em.split`, `pc.gate`, `pc.TOP_ABS` = 20, `pc.MAJORITY` = 5,
`pc.SHORTLIST` = 10, scaffold uniqueness). Forest: `VarianceThreshold(0.01)` + `RandomForestRegressor(random_state=42, n_jobs=1)`
inside one `Pipeline`, refitted on the training side of each split; no preprocessing is fitted outside the training rows.
Censored-EM: `em.censored_em` unchanged (initialise at the bound, fit, project into the interval, at most 10 passes, tolerance 0.01),
with the forest of one thread.

## Anchor (must hold before anything else is read)
C-A. The one-thread censored-EM of the registered rule, run on the five distinct splits, reproduces the 69 stored pool predictions
(`pool69_predictions.csv`, stored to four decimals: every difference under 5e-5), and the vote reproduces the fifteen survivors and the eight shortlisted compounds in order.
If it fails, the other results are not read.

## C1. Simulation validation of the projection (truth known)
Features: the real 6,284-row feature matrix. Five replicates (seeds 1 to 5) of two synthetic truths.
- T_lin: a latent value linear in 100 randomly drawn features (weights N(0,1)), rescaled to mean 6.3 and standard deviation 1.3.
- T_rf: the prediction of a one-thread forest fitted on the 5,370 real point records, rescaled to mean 6.3 and standard deviation 1.3.
Observed value = latent value + N(0, 0.35^2). Censoring mirrors the real composition (9.8% upper, 4.7% lower): the observed values
below their 9.8th percentile are reported as upper-censored at that percentile, those above their 95.3rd percentile as lower-censored
at that percentile. Split: the anchor convention (D1 sorted, D2 train first, D3 divert) on the real scaffolds.
Methods: (a) drop the censored records; (b) censored records entered at their bound; (c) censored-EM, the registered rule;
(d) ridge-Tobit (below). Score, on every held-out record against the noise-free latent value: RMSE, Spearman, and the lowest
prediction among held-out records whose latent value is below the lower bound (does the model reach below the bound?).
Predictions. P1: censored-EM has a lower RMSE than (b) in at least 4 of 5 replicates for each truth. P2: censored-EM has a lower RMSE
than (a) in at least 4 of 5 replicates for each truth. P3: ridge-Tobit has a lower RMSE than censored-EM under T_lin in at least
4 of 5 replicates, and censored-EM a lower RMSE than ridge-Tobit under T_rf in at least 4 of 5. P4: for both truths the
censored-EM Spearman is within 0.03 of the best method in at least 4 of 5 replicates.

## C2. Convergence of the real fit
For each of the five distinct splits, censored-EM run for 30 passes with no tolerance stop. Recorded per pass: the maximum and the
mean absolute change of the imputed targets, the number of censored targets that moved by more than 0.01, and the predictions on the
69-compound pool. Reported: the pass at which the maximum change first falls below 0.01 (or "not within 30"); the Spearman
correlation and the top-20 overlap between the pool predictions at pass 10 (the registered stop) and at pass 30.
Predictions. P5: the maximum change falls below 0.01 within 30 passes in at least 4 of 5 splits. P6: Spearman between pass 10 and
pass 30 pool predictions is at least 0.99 and the top-20 overlap is at least 18 in each split. P7: the mean absolute change falls
monotonically (non-increasing from pass 2) in each split.

## C3. Sensitivity to initialisation
Same five splits, 10 passes and tolerance as registered, but the censored targets start (i) at the bound (registered), (ii) 1.0 pKi
inside the interval away from the bound (bound minus 1.0 for upper-censored, plus 1.0 for lower), (iii) at the bound minus a draw from
U(0, 1.5) for upper-censored and plus for lower (seed 42, fixed), (iv) at the mean of the training point records clipped into the
interval. For each: the pool predictions, their Spearman correlation with (i), top-20 overlap with (i) per split, the selection by the
registered vote (survivors and shortlist) and the Jaccard index against the fifteen and against the eight.
Predictions. P8: the Spearman correlation with (i) is at least 0.98 for every initialisation and every split. P9: the Jaccard index of
the survivors against the fifteen is at least 0.8 for every initialisation. P10: at least 6 of the 8 shortlisted compounds are
shortlisted under every initialisation.

## C4. A censored-likelihood model beside the forest
Ridge-Tobit: a linear model, Gaussian censored likelihood (point: normal density; upper-censored: probability that the latent value
is at or below the bound; lower-censored: at or above), L2 penalty on the weights (not on the intercept), parameters (w, b, log sigma),
L-BFGS-B with analytic gradients, 500 iterations at most. Features: the same raw features, `VarianceThreshold(0.01)` and a
`StandardScaler` fitted on the training rows. Penalty chosen from {1, 10, 100, 1000, 10000} by 3-fold `GroupKFold` on the training
scaffolds, maximising the held-out censored log-likelihood, inside the training side only. Run on the eight conventions and on the
69-compound pool; the same vote rule gives its own survivors and shortlist.
Score: Spearman, R2 and RMSE on held-out point records; the share of held-out upper-censored records predicted at or below the bound;
the pool predictions' Spearman correlation with the forest's (mean of the eight conventions); the overlap of its shortlist with the eight.
Predictions. P11: ridge-Tobit has a lower held-out Spearman on point records than the forest in all eight conventions (the forest:
0.728 to 0.773). P12: the Spearman correlation between the Tobit and forest mean pool predictions is at least 0.5. P13: at most 4 of
the Tobit shortlist's compounds are among the eight. P14: ridge-Tobit's share of held-out upper-censored records predicted at or
below the bound is higher than the forest's in all eight conventions.

## C5. Determinism
The anchor run is executed twice with one thread and the pool predictions of the two executions compared byte for byte; the anchor
convention is also run once with three threads and the maximum absolute difference of the pool predictions from the one-thread
result is recorded.
Predictions. P15: the two one-thread executions are byte-identical. P16: the three-thread run differs from the one-thread result in at
least one pool prediction.

## Reporting
All sixteen predictions are graded and reported, met or not. Where the censored-EM does not converge by the registered tolerance the
article says so. The analyses enter the supplement as a section beside the sensitivity analyses of ADR-AIIM-399/400 and the Methods
sentence on the censored fit is amended to match what was found. The fifteen, the eight and every number of the frozen study stay.

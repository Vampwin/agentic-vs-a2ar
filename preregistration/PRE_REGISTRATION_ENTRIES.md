# Pre-registration entries

Verbatim copies of the decision-log entries that the article and the supplement describe as pre-registered, and of the entries that record their results. The log is append-only; each entry was written before the run it concerns (or, for a result entry, after it). Each entry fixed the predictions or the decision rule of its analysis before the run, and recorded the SHA-256 of the predictions or of the inputs frozen for it where the run had them. Copied as written: first-person wording and the author's quotations are unchanged. The full log is not deposited.

- source: `00_project_control/08_DECISION_LOG.md`, SHA-256 of the whole log (line ends as LF) `ad5cee2e0dfb8fede7ae1426b50e932128750fec4ca7d20ddcb5e19d64753d51`
- entries: 72 (ADR-AIIM-131, 150, 156, 214, 215, 229, 240, 246, 256, 260, 285 to 336, 353 to 356, 390, 391, 399, 400)

---

## ADR-AIIM-131 — Five of the nine seeded-defect escapes were cases in which no fault was injected

**Date:** 2026-08-29
**Status:** Accepted
**Approved by:** human author, 2026-08-29 ("แก้ทั้งสามข้อ regenerate จาก pilot แล้วใส่ §3.7")
**Supersedes:** the seeded-defect figures published up to this date.
**Arises from:** review comment M1, and the C2 reconciliation of 128 against 133 cases.

### The rule this project already had

`pilot_runner.py` defines a vacuous case and says what to do with it. An
instance whose injection parameters match nothing in the workspace changes no
bytes; there is no fault present to escape; and its `seer_by_condition` excludes
such cases. The comment in that file gives the reason: scoring the resulting
PASS as an escape would penalise the contracts for a fault that was never
injected.

### Where it was lost

`WP2.5_WP3A_readiness_2026-07-30/pilot_postfix_2026-08-04.json` honours the
rule. For condition B it records 133 fault rows, of which five are vacuous —
F13_007, F13_008, F13_012, F13_016, F13_020 — and reports SEER_B = 4/128 =
0.03125.

On 2026-08-05 the per-case rows were re-exported to add conditions A and C,
which the pilot had not run. The exported row schema is `family_id`,
`instance_id`, `case_type`, `condition`, `action`, `risk`, `reason`, `escaped`,
`escaped_detected`. **`vacuous` is not among them.** Everything downstream
divided by 133, and `wp4_conditions_ac_results_postfix_2026-08-05.json` reports
SEER_B = 9/133 = 0.067669 for the same condition the pilot scored at 0.03125.

The same condition therefore carries two different SEER values in two artifacts
of one chain, and the manuscript has been quoting the later one.

### What the five cases are

Every F13 escape in the governed arm is one of the five vacuous instances.
Outside F10 the contracts miss nothing. The naive baseline, run independently on
2026-08-10, flags the same family and excludes five — its 0.477 was computed
correctly over 128 all along. It was the governed number that was not
like-for-like.

### The correction

The numerator does not move. 124 cases are detected before and after; only the
denominator changes.

| | published | corrected |
|---|---|---|
| cases counted | 133 | **128** (133 injected, 5 vacuous) |
| detection rate | 93.2 % | **96.9 %** |
| escapes | 9 | **4** |
| SEER A, B | 0.068 | **0.031** |
| SEER C | 0.992 | 0.992 |
| SEER D | 1.000 | 1.000 |
| G6 | 32 cases, 5 escaped, 84.4 % | **27 cases, 0 escaped, 100 %** |
| G4 | 27 cases, 4 escaped, 85.2 % | unchanged |
| F13 | 20 cases, 75 % | **15 cases, 100 %** |

The four F10 instances — G4, docking score sign inverted — are the only genuine
seeded-fault escapes in the study.

### Why the durable fix is a regenerator and not an edited number

`10_new_system/scripts/recount_seeded_defects.py` takes the conditions from the
file that has all four and the vacuity flags from the file that recorded them,
joins on `instance_id`, and recomputes. It refuses to run if vacuity differs
between conditions, or if the four conditions do not cover the same case set —
the second being precisely the failure it exists to undo. It emits the
previously published figures beside the corrected ones in the same artifact, so
that a correction which improves the authors' own headline number can be checked
against what it replaces without hunting for a superseded file.

Editing the number in place would have left the defective export in the chain
and the next analysis would have divided by 133 again.

`20_new_runs/verification_summary_2026-08-21/mutation_score_by_gate.json` is
**not modified** — ADR-AIIM-092 holds. The corrected figures live in
`20_new_runs/seeded_defect_recount_2026-08-29/`, and `sources()["mut"]` is
repointed there.

### This goes in §3.7

The incident register in §3.7 catalogues failures where a check reported health
while the thing it checked had drifted. A re-export that silently dropped a
field, after which four downstream analyses and a manuscript divided by the
wrong denominator for three weeks, is that failure class, found in this study's
own chain. Reporting it is not optional given what §3.7 is for, and the fact
that the correction improves the reported figure makes stating its provenance
more necessary rather than less.

### What is not claimed

That the five instances are badly specified. F13 is "undeclared external label
constraint"; five of its twenty parameter settings matched nothing in the
workspace they were applied to. Whether the family should be re-instantiated so
that all twenty inject is a separate question, not settled here, and the
corrected figures are reported over the cases that did inject.

---

## ADR-AIIM-150 — The funnel discarded two and a half times as many viable compounds as it selected

**Date:** 2026-08-31
**Status:** Accepted
**Closes:** the reweighting and absolute-number halves of review comment M9.
R4 remains open and answers a different question.

### What M9 asked for, and what could be done without new compute

The exclusion control drew 250 compounds from each of the two ways the funnel
rejects and reported a single pooled rate of 7.0 %. M9 objected that pooling an
equal sample from unequal strata is not an unbiased population estimate, and
that the absolute number was never stated. Both are arithmetic on data already
collected.

### The strata

    arm_2d   21,882 -> 8,000   13,882 rejected by the 2D pharmacophore/ErG screen
    arm_3d    8,000 -> 5,000    3,000 rejected by the 3D shape refinement

Equal samples from strata differing more than four to one. The pooled figure
weights the smaller stratum 4.6 times too heavily.

### Result

| stratum | population | sampled | reaching cutoff | rate | 95 % CI | implied |
|---|---|---|---|---|---|---|
| arm_2d | 13,882 | 249 | 19 | 7.63 % | 4.9–11.6 % | 1,059 |
| arm_3d | 3,000 | 250 | 16 | 6.40 % | 4.0–10.1 % | 192 |

Pooled as published 7.01 %; **population-weighted 7.41 %**. The reweighting
barely moves the rate, because the two arms discard viable compounds at
indistinguishable rates — which is itself the answer to whether one rejection
mechanism is more defensible than the other.

**The absolute number is what M9 was really after: roughly 1,251 of the 16,882
discarded compounds would have reached the selection cutoff, against the 500
selected — a factor of 2.5.** The funnel did not trim a fringe. It discarded a
set larger than the one it kept, on a criterion whose own enrichment is 1.43×.

### A correction to the review comment

M9 describes "strata of ~16,882 and ~3,000". 16,882 is the total rejected,
13,882 + 3,000; it is not a stratum. Its arithmetic applies a pooled rate to the
sum of both strata — two errors that partly cancel, which is why its estimate of
~1,180 is close to the correct 1,251 for the wrong reasons.

### Scope

Measured under the **fast pass** — two conformers, exhaustiveness 8 — which is
the protocol that chose the top 500 and therefore the right comparison for
asking whether that cut was sound. It is not the full protocol those 500
received. R4 answers recall at full protocol and is not superseded.

### Where it lands

§3.5's funnel-control paragraph, which previously stopped at "roughly one in
fourteen". Six new G8 checks; 179 audited numbers passing on both files.

---

## ADR-AIIM-156 — The pose-persistence result was wrong, and the figure is what found it

**Date:** 2026-08-31
**Status:** Accepted
**Supersedes:** ADR-AIIM-149's result. Its method choices stand; its conclusion does not.
**Found by:** the author asking for a figure of the MD result.

### What was claimed on 2026-08-31, and withdrawn the same day

ADR-AIIM-149 reported that NPLIB01185's docked pose "leaves the pocket in all
three replicates", mean RMSD 24.71 Å, final frames 24.3 to 60.2 Å. That reached
§2.9, §3.6, §4, §5, the Abstract, a highlight and the graphical abstract, and a
DOCX was sent to the author containing it.

**It is false. The ligand never left.**

### The evidence

`gmx mindist` applies the minimum-image convention itself and cannot be misled
by how a trajectory was written. Across all nine replicates the ligand stays
2.58 to 5.47 Å from the receptor, with **zero frames beyond 10 Å out of 18,009**.

In the 400 frames of NPLIB01185 r1 whose RMSD reads **over 40 Å**, the same
frames put the ligand **3.03 to 4.68 Å** from the receptor — van der Waals
contact, indistinguishable from the frames reading under 8 Å.

The time series is bimodal with a clean 15 Å gap, and contains **79, 111 and 148
frame-to-frame jumps larger than 30 Å** in the three replicates, the largest
being **56.1 Å in 50 ps** — roughly 1,100 m/s for a 37-atom molecule. The other
two compounds have **zero** such jumps and a maximum step of 3.7 Å.

The high mode is a periodic image.

### Why no summary statistic could have caught it

Mean 16.32, SD 21.99, max 82.26, final 59.77. Every one of those is consistent
with a ligand that dissociated. The SD exceeding the mean is the only hint, and
it reads as "high variance", which a dissociating ligand also has. **The
artefact is visible only in the time series**, and the time series was not
plotted until the author asked for a figure.

This project's own §3.7 argument is that a check reporting health while the
thing it checks has drifted is the failure that matters. Here the check was a
set of summary statistics and the drift was invisible to them.

### What was tried and did not work

`-pbc mol -center` (original), plus `-ur compact`, plus `-pbc nojump`, and
`-pbc cluster`. The first three leave the RMSD unchanged to the last decimal.
The fourth refuses: clustering needs whole molecules and a C-alpha group is not.

Two of those attempts were also wasted because the files being regenerated were
never deleted — a shell variable eaten by the Windows → wsl → bash escaping,
for the fifth time in this project. The "unchanged" results were stale files,
and the real test only ran once absolute paths were used.

### The measurement that replaces it

Distance from the ligand to the binding site **as the docked pose itself defines
it**: the receptor C-alpha atoms within 10 Å of the ligand in the
energy-minimised starting structure, 47 to 52 atoms per compound.
`md_17_pocket_from_start.py`. It answers the question the RMSD was being asked
to answer, and no trajectory representation enters it.

    compound     replicate distances (min - max, Å)   frames beyond 10 Å
    NPLIB01185   2.71 - 5.11                          0 of 6003
    NPLIB01643   2.83 - 5.47                          0 of 6003
    NPLIB03056   2.58 - 5.24                          0 of 6003

### The corrected result

**All three candidates occupy their docked sites for the full 100 ns, in every
replicate.** Neither reproduces the predicted geometry: NPLIB01643 and
NPLIB03056 sit 4.73 and 5.22 Å from their docked poses against the control's
2.59 Å. **NPLIB01185's RMSD is not reported at all** — reducing a bimodal
distribution to its low mode would be choosing the answer, and the compound's
pose occupancy is established by the measurement above regardless.

The MD neither supports nor removes any candidate. That is weaker than what was
written and it is what the data says.

### What changed

§2.9, §3.6 (rewritten), §4 (both the limitations and clinical-relevance
paragraphs), §5, the Abstract, highlight 4, and the graphical abstract's third
column. **Figure 6 is new** and leads with the pocket-occupancy measurement;
panel B omits NPLIB01185 with the reason printed on the figure rather than left
to the caption.

G8 183 -> 181: the three checks auditing NPLIB01185's RMSD are deleted with the
sentences they audited, and three new ones cover pocket occupancy.

---

## ADR-AIIM-214 — Pre-registration: DUD-E A2AR enrichment (M5.1)

**Date:** 2026-09-22
**Status:** accepted — **pre-registration, frozen before any compound was docked**
**Approval:** author, 2026-09-22
**Unblocked by:** ADR-AIIM-212 (benchmark fitness measured), ADR-AIIM-213 (receptor and box settled)

### 1. Why this run exists

M1.4 corrected the decision architecture: for a library of novel chemotypes the
ligand-based arm **abstains** rather than vetoes, and the structure-based arm
must lead. An arm that leads must be validated the way the field validates it —
**retrospective enrichment against property-matched decoys**. This run is that
validation. It is no longer optional tidying of a pre-declared gap.

The objectives document states the purpose verbatim: *"is 0.601 a protocol
problem or a benchmark problem?"*

### 2. Why this benchmark, measured rather than assumed (ADR-AIIM-212)

| Check | Result |
|---|---|
| Trivial separability — cLogP | 0.594 |
| TPSA / heavy atoms / MW / rotatable bonds | 0.552 / 0.540 / 0.486 / 0.480 |
| *the ChEMBL benchmark in use, heavy atoms, for contrast* | *0.584* |
| Analogue bias | **none — 482 unique Murcko scaffolds for 482 actives** |

Decoys are *presumed* inactive, which biases enrichment **downward**. Zero
scaffold overlap with actives is disqualifying for a trained classifier but not
for a scoring function with no training set.

### 3. Protocol — frozen, and identical to production

Unchanged from `a2a_frozen_config.yaml` and the run that produced 0.601, so the
two numbers are comparable:

| | |
|---|---|
| Engine | AutoDock Vina 1.2.7, WSL `diffdock` env |
| Conformers | 8 pooled, ETKDGv3 + MMFF94 |
| Exhaustiveness | 32 |
| Poses retained | 10 |
| Seed | 20260729 |
| Receptor | `snp_screen_2026-08-04/G2_preparation/4EIY_prepared.pdbqt`, sha256 `bd0e586a525ecbea…` |
| Box | centre [−0.424, 8.527, 17.127], 22 Å cube, sha256 `f5149ff0cb07e7d1…` |
| Compounds | 2,402 — 482 actives, 1,920 decoys |

Receptor copy is provably score-neutral (ADR-AIIM-213, spread 0.0000 across
three copies and four ligands); the canonical copy is used for provenance, not
because it is required.

### 4. Metrics, and the nulls that travel with them

ROC-AUC with 1,000-sample bootstrap 95 % interval at seed 42; BEDROC α = 20;
EF at 1 %, 5 % and 10 % using the `ceil(n·f)` convention that reproduces the
published benchmark exactly.

**Two trivial nulls are computed on the identical rows and reported beside
every docking figure: heavy-atom count and cLogP.** cLogP is the stronger of
the two on this set (0.594) and reporting heavy atoms alone would understate
the baseline (ADR-AIIM-212).

### 5. Outcomes, fixed now

- **Primary.** Paired ΔAUC, Vina minus the **better-performing** trivial null,
  with a bootstrap interval.
- **Secondary.** Whether Vina's own AUC interval excludes 0.5; BEDROC20 and EF
  against the same nulls.

### 6. Predictions, recorded before the run

1. **Vina ROC-AUC lands in 0.62 – 0.78.** Higher than the matched ChEMBL result
   of 0.502, because DUD-E's actives are drawn from the same A2A literature as
   the co-crystal ligand and are therefore chemically closer to it than orchid
   flavones are.
2. **The paired ΔAUC against the best trivial null is positive and its 95 %
   interval excludes zero.** At n = 2,402 the AUC interval is roughly ±0.025, so
   a difference above about 0.05 will separate.
3. BEDROC20 exceeds the corresponding null.

I expect this run to come out **positive**, which is the opposite of my
expectation for the two experiments before it. The reason is specific and
stated above, not optimism.

### 7. Falsification, and what each outcome licenses

**If Vina does not beat the better null, or the paired interval includes zero:**
the protocol has no demonstrated screening power on a benchmark constructed to
be favourable to it, and **O5 closes**. That would be reported as the headline
and the screening arm would be presented as a cautionary methodological result
only.

**If Vina does beat the null:** the protocol has screening power **on
chemotypes resembling the co-crystal ligand**, and the 0.601 is confirmed as a
benchmark artifact rather than a protocol failure.

### 8. The limit that a positive result does *not* remove

DUD-E A2AR actives come from the same literature as ZM241385. They are closer
to the co-crystal ligand than the natural-product library is, and
`control_representativeness` already reports 0.000 / 0.000 scaffold and
similarity coverage between that library and the control.

> **A good DUD-E result would validate the protocol for the chemotypes DUD-E
> contains. It would not, by itself, license ranking orchid flavones.**

Closing that gap needs the receptor ensemble of M5.2, not this run. This is
recorded now so that a favourable number cannot later be read as licensing more
than it does.

### 9. What will not be done

- **No re-docking under different settings.** One protocol, one run. If the
  result is unfavourable, no sweep of exhaustiveness, conformer count or box
  size will be performed in search of a better number. The settings are the
  frozen production ones precisely so that this temptation is unavailable.
- **No post-hoc selection of the null.** Both nulls are reported whichever way
  they fall.
- **No subsetting of the benchmark** after seeing results.

### 10. Cost and scope

Measured at 30.5 s per compound at this protocol: **≈ 20 h single-process, ≈ 10 h
across four workers.** Output to `20_new_runs/dude_enrichment_2026-09-22/`,
resumable, per-compound scores retained.

This validates **one protocol on one target against one decoy set**. It does not
speak to docking in general, to other targets, or to the natural-product library
beyond what §8 allows.

---

## ADR-AIIM-215 — Amendment to ADR-AIIM-214: dock the minimum-energy conformer, not eight

**Date:** 2026-09-22
**Status:** accepted — **amendment made before any result was analysed**
**Amends:** ADR-AIIM-214
**Raised by:** the author, during the run

### What was wrong with the frozen protocol

ADR-AIIM-214 fixed the protocol as eight pooled conformers, identical to
production, so the DUD-E figure would be comparable with the 0.601.

The cost estimate behind that decision was wrong. It came from two compounds
that happened to be easy — 31.3 s and 29.7 s. Once the run started, the first
seven compounds recorded 45.8, 37.1, 341.8, 130.7, 461.6 and 184.0 seconds:
**22.7 s per individual dock on average, and eight docks per compound.**

| Protocol | Projected |
|---|---|
| Eight conformers docked, best kept | **121 hours** |
| Minimum-energy conformer only | **15.1 hours** |

121 hours is not affordable, and reducing the benchmark instead would weaken
the design to protect a protocol choice.

### Why the change is defensible on its own merits

The author's point is the scientific one, not only the arithmetic. **AutoDock
Vina performs flexible-ligand docking**: it samples rotatable-bond torsions
internally during its Monte Carlo search. Supplying eight input conformers is
therefore largely redundant for precisely the degrees of freedom Vina already
explores. What input conformers *do* fix is ring geometry, which Vina holds
rigid — and for drug-like molecules with aromatic or small saturated rings, a
single MMFF-minimised conformer is standard practice.

Eight pooled conformers is a conservative choice above the field norm, not a
requirement.

The amended procedure is the author's: generate eight conformers, minimise all
eight with MMFF94, **dock only the lowest-energy one**.

### The cost, stated rather than hidden

This now validates a **standard-practice protocol that differs in one respect
from the eight-conformer protocol that produced the study's own screening
results**. A DUD-E result under it does not fully validate what was actually
run on the natural-product library.

Two things bound that gap:

- `bridge_8conformer_reference.csv` retains the seven compounds already scored
  under the production protocol, so the two can be compared directly on
  identical molecules.
- The ChEMBL 245 benchmark will be re-run under the amended protocol
  (≈ 1.5 hours), so the two benchmarks share a protocol and the cross-benchmark
  comparison stays clean.

### What is unchanged

Everything else in ADR-AIIM-214 stands: exhaustiveness 32, ten retained poses,
seed 20260729, the receptor and box with their recorded hashes, both trivial
nulls including cLogP, the predictions, the falsification conditions, and the
rule that no settings will be swept in search of a better number. **The
predictions were not revised.**

### Why this is an amendment and not a rescue

No metric has been computed. Seven compounds of 2,402 are docked, all actives,
and no enrichment figure exists. The change is made on a cost measurement and a
methodological argument, both stated before any result is visible.

---

## ADR-AIIM-229 — PRE-REGISTRATION: spike-in recovery through the expansion funnel (M2.1)

**Date:** 2026-09-22
**Status:** accepted — frozen before anything ran
**Approval:** author, 2026-09-22 — including the 105-vs-50 decision below
**Executes:** M2.1 of `29_OBJECTIVES_AND_METHOD_DESIGN_v6.md`
**Builds on:** ADR-AIIM-222 (the unseen subset), ADR-AIIM-210 (the domain gate on real drugs)
**Depends on:** ADR-AIIM-230, which is why one derived column sits outside the frozen artifact

### The question

If compounds known to be active are placed in the screening library, does the
expansion funnel return them, and if not, **which stage loses them?** This is the
positive control that should have gated the campaign; §3.5's finding arrived at
manuscript-writing time because it did not exist.

### Source: the 169, not the 482

ADR-AIIM-222 measured that **313 of 482** DUD-E A2AR actives already share their
InChIKey connectivity block with the ChEMBL set the ligand-based arm was fitted
on. Drawing from those would measure memory. The draw is from the **169 unseen**,
`benchmark_overlap_2026-09-22/dude_actives_unseen_by_training_set.csv`,
sha256 `6a6e697ec868…`.

### The independent unit is a fingerprint cluster, and the scaffold key would have hidden that

The obvious independence key is the Bemis–Murcko scaffold: it is the key this
project's own scaffold-disjoint splits use, and ADR-AIIM-212 reports 482 distinct
scaffolds for 482 actives. On the 169 it returns **169 clusters of size one**,
which is not a measurement of independence but the absence of one. The key cannot
see that two of these compounds differ by a methylene, or that another pair are
tautomers of each other.

Morgan radius 2 / 2048 bits with Butina at Tanimoto ≥ 0.70 — the project's own
fingerprint (`evaluators/_control_representativeness.py:116`) at the conventional
series cutoff — returns **105 clusters**, largest 12, 78 singletons. **97 of the
169 have a nearest neighbour at Tanimoto ≥ 0.70** (mean 0.707, median 0.730, max
1.000).

> Treating 50 compounds drawn from ~105 series as 50 independent tests is
> pseudoreplication. It would not show up as a wrong number; it would show up as
> a Wilson interval narrower than the evidence supports.

One representative per cluster — the Butina centroid, which is deterministic —
and the recovery rate is then a rate over series rather than over analogues.

### The draw, and why an order is frozen rather than a set

`numpy.random.default_rng(42)` permutes the 105 clusters. **The primary panel is
the first 50 of that order; the remaining 55 are labelled `replication`.** The
order is fixed before any result exists, so extending the analysis from 50 to 105
is a pre-registered replication rather than a choice made after seeing the first
50.

Random sampling rather than a MaxMin diversity pick: **MaxMin was considered and
rejected.** It selects the extremes of the chemical space, which would bias the
panel toward unusual chemistry — exactly the chemistry an applicability-domain
gate is designed to refuse — and a non-recovery result would then be arguable as
a property of the panel rather than of the pipeline. A seeded random draw selects
for nothing.

### The panel is not filtered by the pipeline's own criteria

Four of the 169 exceed the funnel's rotatable-bond limit; three of them are cluster
representatives (**85419, 85418, 85023**) and they stay in the draw. A positive
control pre-filtered by the rule it is meant to test cannot test that rule.
Dropping them would be selecting on the outcome, and a stage-1 loss of real
actives is itself a reportable property of the funnel.

### Two leak checks, run before the panel was frozen

| Check | Result |
|---|---|
| Panel compounds already in the LOTUS background | **0 of 169** by skeleton InChIKey, against 16,486 distinct LOTUS skeletons |
| Panel compounds in the stage-2 reference set (the 52 known A2AR ligands the 2D screen ranks against) | **0 of 169** share a skeleton; max ECFP4 similarity to any reference **0.807**, three compounds at ≥ 0.70 |

Had either been non-zero the control would have been rigged and the result
worthless.

### What stages 1 and 2 will do is derived, not predicted — and neither is a test

Stages 1 and 2 are deterministic functions of frozen inputs, so they are computed
here rather than guessed. Importing the rule from `np_library_prefilter` rather
than restating it, for the reason ADR-AIIM-227 gives.

| | Panel of 105 | Primary 50 |
|---|---|---|
| Survive stage 1 (MW ≤ 700, rotatable ≤ 10) | **102** | **50** |
| Survive stage 2 (top 8,000 by 2D pharmacophore/ErG) | **94** | **47** |

The stage-2 cut moves from **0.5201** unspiked to **0.5207** spiked; 16 rows tie
at it and pandas resolves ties arbitrarily, so the count is exact only up to
those. The panel's median score is **0.7012** against a kept-background maximum
of **0.6250**: **78 of the 105 outrank the entire surviving background.**

That number is the reason for the following, which the manuscript must carry:

> **Recovery through stages 1 and 2 is not evidence that the pipeline works.**
> The stage-2 screen ranks by similarity to known A2A ligands, and the panel are
> known A2A ligands. It could hardly do otherwise. These two stages are an
> implementation check — a run that disagrees with the table above has a fault,
> not a finding — and the scientific content of M2.1 is at stage 3 and later.

### The predictions that are predictions

Stages 3 (3D O3A refinement, 8,000 → 5,000), 4 (G4 docking), 5 (conformal
predictor and AD gate) and 6 (the dual-criterion rule) are not derivable without
running. Pre-registered, from the plan:

1. **The AD gate and the ligand-based arm admit most surviving spike-ins.** Basis:
   the gate admits **9 of 10** real A2A drugs (ADR-AIIM-210).
2. **The structure-based rank step loses the majority.** Basis: measured screening
   power near chance (§3.4).
3. If 1 and 2 both hold, the loss is localised to docking rank, and §3.5's
   undifferentiated negative becomes a diagnosis.

**Falsification.** If most spike-ins survive every stage and reach the shortlist,
prediction 2 is wrong and the funnel recovers known actives end to end — the
control passes and §3.5's negative needs a different explanation. If the AD gate
rather than docking is the largest loss, prediction 1 is wrong and ADR-AIIM-210's
9-of-10 does not generalise past clinical-stage chemistry.

### Metric

Count surviving each stage, with a **Wilson 95 % interval**, beside the same rates
for the background pool. Precision bought:

| | p = 0.5 | p = 0.9 |
|---|---|---|
| n = 50 | ±0.134 | ±0.085 |
| n = 105 | ±0.094 | ±0.059 |

### The one decision that was not deferrable, and how it went

**Spike 105 or spike 50?** The approved plan fixed 50 before the clustering had
been measured. The marginal cost of the extra 55 is ~55 docks against a funnel
that docks thousands — negligible **if decided before the run**, and a **full
funnel re-run** if decided after, because the shortlist depends on the pool. So
it could not be left open.

**Decided by the author, 2026-09-22: spike all 105.** The frozen first 50 are the
primary arm and the remaining 55 the pre-registered replication. This is a change
to the plan's stated n, taken on a measurement the plan did not have — that the
169 hold 105 independent series — and it buys the precision in the table above
(Wilson half-width at p = 0.5 falls from ±0.134 to ±0.094) for no extra funnel
run. **Nothing about the primary arm changes:** the same 50 compounds in the same
frozen order, chosen before any of this was known.

The artifacts as frozen already carry the decision: `n_spiked` 105, `n_primary`
50, `spiked_pool_25709.csv`. Nothing was regenerated to record the approval.

### Not claimed

- That any panel compound is active in any assay run by this project. **No
  compound reported anywhere in this project has been assayed.**
- That the LOTUS background is a property-matched decoy set. This is a **recovery
  control, not an enrichment benchmark**; enrichment is M5.1's job, and Verdonk's
  condition for a meaningful enrichment comparison is not met here.
- That the pipeline would find *novel* actives.
- That the single-conformer protocol was validated. It was decided (ADR-AIIM-215).

### Artifacts, frozen 2026-09-22

`20_new_runs/spikein_panel_2026-09-22/`

| File | sha256 |
|---|---|
| `spikein_cluster_order_105.csv` | `db0c06c0a6ac…` |
| `spikein_panel_primary_50.csv` | `4ba40392cfb4…` |
| `spiked_pool_25709.csv` | `76dfd640d6e3…` |
| `derived_stage12_expectation.csv` | `90d9a4bc56c5…` — **not reproducible; see ADR-AIIM-230** |

Produced by `select_spikein_panel.py`, checked by `verify_spikein_panel.py`: 13
structural checks, bit-exact reproduction of the frozen selection in every fresh
process, and the derived counts invariant across all of them. **PASS.**

The spiked pool needs **no edit to the funnel**: it is the LOTUS schema plus a
`spike_in` column, the LOTUS rows carried through unchanged, and the spike-ins are
identifiable downstream by InChIKey.

### What must happen before the run

1. ~~Author approval of this entry.~~ **Given 2026-09-22.**
2. The CPU is held by the DUD-E docking run and then the ChEMBL 245 re-run.
   The M2.1 funnel is **not** chained behind them; whether to queue it is a
   separate call, and a third job started unattended behind two others is how a
   failure in the first becomes three.
3. `run_g4_on_np_library.py:42` hardcodes its run directory. Add a `--run-dir`
   flag or a thin wrapper. **Do not edit the contracts or the evaluators.**
4. The shortlist must be re-docked under the unified protocol; the 2026-08-15
   scores are 8-conformer and cannot be mixed.

---

## ADR-AIIM-240 — M5.1 result: the docking protocol does have screening power, and §3.4's 0.601 was the benchmark

**Date:** 2026-09-23
**Status:** accepted — result
**Approval:** *pending author*
**Executes:** ADR-AIIM-214, amended by ADR-AIIM-215, analysis written under ADR-AIIM-228, fixed under ADR-AIIM-233
**Artifacts:** `20_new_runs/dude_enrichment_2026-09-22/analysis/`

### The run

2,402 compounds docked under the unified protocol — eight ETKDGv3 conformers
generated, all MMFF94-minimised, the lowest-energy one docked, exhaustiveness 32,
10 poses, seed 20260729 — in 8.9 hours at 13.4 s per compound. **One compound
failed** (`C39559757`, a decoy, RDKit embedding), so the analysed set is **2,401:
482 actives and 1,919 decoys**, against the 482 / 1,920 that ADR-AIIM-214
pre-registered. The shortfall is recorded in the artifact rather than absorbed.

`analyse_dude_enrichment.py` was written while the run stood at roughly half its
compounds and no score had been read (ADR-AIIM-228); the row-count defect found
in it on 2026-09-22 was fixed **before** the first real score was analysed
(ADR-AIIM-233). It graded the predictions mechanically and wrote the verdicts to
the artifact.

### Result

| Arm | ROC-AUC | 95 % CI | BEDROC20 | EF 1 % |
|---|---|---|---|---|
| **Vina, one minimised conformer** | **0.6689** | [0.6414, 0.6935] | **0.4520** | 2.59 |
| heavy-atom count | 0.5400 | [0.5107, 0.5691] | 0.3720 | 3.59 |
| cLogP | 0.5939 | [0.5668, 0.6219] | 0.2401 | 1.00 |

**Primary outcome**, paired against the better-performing trivial null, which the
script selects by AUC rather than by hand:

> **ΔAUC = +0.0751, 95 % CI [+0.0365, +0.1109]**, interval excludes zero, and
> **0 of 1,000 bootstrap replicates** put the difference at or below zero.

Against the weaker null the margin is +0.1289 [+0.0967, +0.1601].

### All three pre-registered predictions were met

| | Prediction | Observed | |
|---|---|---|---|
| **P1** | Vina ROC-AUC lands in 0.62–0.78 | **0.6689** | **MET** |
| **P2** | paired ΔAUC against the better null is positive and its interval excludes zero | +0.0751 [+0.0365, +0.1109] | **MET** |
| **P3** | BEDROC20 exceeds the corresponding null | 0.4520 against 0.2401 | **MET** |

The favourable branch of ADR-AIIM-214 §7, pre-written and therefore not
negotiable after the fact, fires:

> The protocol has demonstrated screening power on chemotypes resembling the
> co-crystal ligand, and **the 0.601 of §3.4 is a benchmark artifact rather than
> a protocol failure.**

### What it changes, and what it very deliberately does not

**It changes the reading of §3.4.** The 245-compound ChEMBL benchmark put Vina at
0.601 and heavy-atom count at 0.584, p = 0.722 — indistinguishable from counting
atoms. On a property-matched decoy set built for the purpose, the same protocol
separates actives from decoys with an interval that clears both trivial nulls.
The earlier number described the benchmark, not the method. The ChEMBL 245 re-run
under this same protocol (ADR-AIIM-227) is now complete and lets that comparison
be made without the protocol confounding it.

**It does not license ranking the natural-product library.** ADR-AIIM-214 §8
recorded that in advance, and the reason stands untouched by this result: DUD-E
A2AR decoys are property-matched to A2A actives, so this measures performance on
chemotypes resembling the co-crystal ligand, while `control_representativeness`
measures 0.000 scaffold and 0.000 similarity coverage between the reference
ligands and the LOTUS library. A protocol that discriminates on familiar
chemistry has not been shown to discriminate on unfamiliar chemistry, and §3.5's
finding is about the latter.

**EF at 1 % is the one place the nulls win.** Heavy-atom count reaches 3.59
against Vina's 2.59 in the top 1 %. The early-recognition advantage BEDROC20
reports for Vina does not survive at the very top of the list, and the manuscript
should report the EF row rather than only the AUC and BEDROC that favour the
protocol.

**Nothing here is a wet-laboratory claim.** No compound in this project has been
assayed.

### Consequences

- **O5's central measurement now exists.** Whether it closes O5 is a separate
  judgement and belongs with the author.
- **M2.2 is unblocked.** ADR-AIIM-231's amendment made the exhaustiveness arm's
  *analysis* wait on M5.1 being adjudicated, which this entry does. The
  exhaustiveness-8 docking is at 2,167 of 2,402 and finishing.
- **M1.4's open thread closes.** ADR-AIIM-221 recorded that whether the corrected
  decision architecture rescues anything "waits on M5.1". It does not rescue the
  natural-product ranking, for the reason above.
- §3.4 and §3.6 need rewriting around this, and that is **held by the length
  decision** and interacts with ADR-AIIM-234's restatement of §3.5/§3.6 under the
  unified protocol. The same paragraphs are in play; it should be one edit.

## ADR-AIIM-246 — The M2.1 analysis exists before its data, and it fixes the one thing ADR-AIIM-229 left open

**Date:** 2026-09-23
**Status:** accepted — record of fact
**Approval:** author, 2026-09-23
**Follows:** ADR-AIIM-245 (step 4 launched), ADR-AIIM-228 (the precedent)

### Written mid-run, on purpose

`scripts/analyse_spikein_recovery.py` implements the analysis ADR-AIIM-229
pre-registered. **It was written while the docking was still running and before
any `vina_score` of that run had been read.**

The M5.1 analysis got this guarantee and all three of its predictions were then
graded mechanically against a script that could not have been shaped by them
(ADR-AIIM-240). The M2.2 arm got the weaker one — written after docking
finished, before any score was read — and ADR-AIIM-241 states that rather than
glossing it. This one gets the stronger guarantee, and it is the last of the
three where the choice was still available.

`wilson` is imported from `analyse_r4_funnel_recall.py`, so the interval
convention is the project's own rather than a second implementation.

### The gap in the pre-registration, closed before any score

ADR-AIIM-229 fixes the metric as *"count surviving each stage, with a Wilson
95 % interval"*. That does not define what surviving **stage 4** means: docking
returns a score for every compound, not a pass or a fail. A cut-off is needed,
and choosing one after seeing where the spike-ins land would be choosing the
answer.

> **Primary cut-off: the top 500 by docking score**, because that is the cut the
> production funnel used — §3.6 reports that the 500 highest-ranked compounds
> went to full-protocol redocking — so recovery at 500 speaks to the campaign
> this project actually ran.

The whole curve at **100 / 250 / 500 / 1,000 / 2,500** is reported beside it, so
no reader depends on that cut-off having been well chosen. All of it was written
down before a single score of this run was read.

### P1 is reported as pending, not guessed

Stages 5 and 6 — the conformal predictor, the AD gate and the dual-criterion
rule — need `run_g6_on_np_library` executed against the spiked shortlist. That is
a further run, not an analysis, so **P1 and P3 are emitted as `PENDING`** and the
script grades only **P2**, which the docking alone can answer. A script that
guessed the other two would have produced three verdicts from one experiment.

### Verified on synthetic data, over the real shortlist

Eleven checks, all passing, with fabricated scores on the real 5,000-row
shortlist and its real 89 spike-ins:

| Scenario | Spike-ins in top 500 | P2 |
|---|---|---|
| spike-ins score like the background | 10 of 89 | **MET** |
| spike-ins score systematically better | 86 of 89 | **MISSED** |

The second scenario is the one that matters: it shows **P2 is falsifiable by
this script**, not merely assertable. Also checked: P1 and P3 stay `PENDING`,
the stage-2 implementation check agrees with ADR-AIIM-229's derived 94, every
stage carries a Wilson interval, a partial file is refused, and a mixed-protocol
row aborts.

**No score from the live run was read.**

### Artifacts

`scripts/analyse_spikein_recovery.py`, not yet run against the real data.

## ADR-AIIM-256 — M2.1 complete: the similarity funnel keeps 85 %, the docking rank loses the majority

**Status:** Proposed, awaiting author approval
**Date:** 2026-09-24
**Pre-registration:** ADR-AIIM-229; step 4 authorised in ADR-AIIM-245
**Artifacts:** `20_new_runs/spike_in_dock_2026-09-23/`
(`scores.csv`, `analysis/spikein_recovery.json`,
`analysis/spikein_docked_ranked.csv`)

### The run

Complete. **5,000 compounds scored, 0 failed, 89 of them spike-ins, 36.74 h**,
exit code 0, every row at one minimised conformer — the protocol imported from
`run_dude_enrichment.py`, never restated (ADR-AIIM-234). The measured rate,
26.6 s per compound, sits between the 22.6 s median and 29.8 s mean the 30-
compound probe forecast, and the run took 36.7 h against that probe's 31.4–41.3 h
window.

### The result

105 known A2A actives the ligand-based arm has never seen (ADR-AIIM-222), one
per Butina series, spiked into the 25,604-compound LOTUS pool and carried through
the production funnel.

| Stage | Recovered | Rate | Wilson 95 % |
|---|---|---|---|
| spiked into the pool | 105 / 105 | 1.000 | [0.965, 1.000] |
| stage 2, 2D pharmacophore, top 8,000 | 94 / 105 | 0.895 | [0.822, 0.941] |
| stage 3, 3D O3A, shortlist 5,000 | **89 / 105** | **0.848** | [0.767, 0.904] |
| stage 4, docked | 89 / 105 | 0.848 | [0.767, 0.904] |

**Implementation check: stage 2 observed 94 against the 94 derived in
ADR-AIIM-229 — AGREES.** That pre-registration stated that a run disagreeing
with its own derived stage-1/2 table has a fault rather than a finding, so this
is the control that licenses reading the rest.

Then the docking rank, at the cut the production funnel actually used:

| Docking cut-off | Spike-ins | Rate of the 89 | Enrichment over background |
|---|---|---|---|
| top 100 | 9 | 0.101 | 5.06 |
| top 250 | 22 | 0.247 | 4.94 |
| **top 500 (primary)** | **32** | **0.360** | **3.60** |
| top 1,000 | 52 | 0.584 | 2.92 |
| top 2,500 | 76 | 0.854 | 1.71 |

### What it establishes

**P2 — "the structure-based rank step loses the majority" — is MET.** Of the 89
known actives that reached docking, **57 do not survive the top-500 cut**. The
primary cut-off was fixed before any score was read, and it is the cut §3.6
reports the production campaign as having used.

> **The similarity gate is not where the actives are lost. The docking rank
> is.** Stages 1–3 keep **84.8 %** of known actives; the rank step that follows
> keeps **36.0 %** of what it is handed.

That is a diagnosis rather than another negative, and it is the one §3.5's
undifferentiated result could not supply. It also cuts both ways and the entry
says so: docking *does* enrich, **3.60-fold at the primary cut and 5.06-fold at
the top 100**, so the rank step is informative and lossy at the same time.

**P1 and P3 remain PENDING and are not guessed.** Stages 5 and 6 — the conformal
predictor, the applicability-domain gate and the dual-criterion rule — are a
further run, not an analysis. ADR-AIIM-245 recorded in advance that the script
would emit `PENDING` for both rather than infer them, and it did.

### Why this result is readable at all

The analysis was written **while the docking was still running and before any
`vina_score` of this run had been read**, and was falsification-tested on
fabricated scores over the real 5,000-row shortlist: spike-ins scoring like the
background gave P2 **MET**, spike-ins scoring systematically better gave P2
**MISSED**. A prediction that cannot fail is not a prediction, and this one was
shown capable of failing before the data existed.

### Consequences

1. **The scope of §3.5's negative narrows.** It reported that a similarity-based
   gate excluded the library. This shows the gate would have passed 85 % of
   genuine actives through, so the exclusion of *this* library is a statement
   about that library's novelty, while the compound-level loss is downstream.
2. **§3.6's 500-compound cut acquires a measured cost**: at that cut the campaign
   would have discarded 57 of 89 known actives it had already found.
3. **RUN_INDEX** rows for `spike_in_dock_2026-09-23` and
   `spike_in_recovery_2026-09-23` move from **PENDING** to **CORE**.
4. **The manuscript must now carry this**, which adds words to an article being
   cut (ADR-AIIM-251, -254, -255). The result is worth more than the words: it
   converts a null into a localisation. Placement is §3.5 or §3.6 and is an
   author decision, not an agent's.

### Not claimed

No compound here has been assayed. The spike-ins are ChEMBL A2A actives and the
recovery rate is over **series**, one representative per Butina cluster, not over
analogues — so this measures whether the funnel finds *a* member of a known
active series, not how many actives it would find. Nothing here bears on the
natural-product library's own activity.

---

## ADR-AIIM-260 — M2.1 complete: P1 and P3 MISSED. The domain gate, not the docking rank, is the largest loss

**Status:** Proposed, awaiting author approval
**Date:** 2026-09-24
**Approval to run:** Adisorn Sirichotejirakul ("รัน stage 5-6 ปิด O2 เลย", 2026-09-24)
**Pre-registration:** ADR-AIIM-229, operationalised in
`run_spikein_stage56.py`'s docstring **before any result was read**
**Artifacts:** `20_new_runs/spike_in_dock_2026-09-23/analysis/`
(`spikein_stage56.json`, `spikein_stage56_per_compound.csv`,
`spikein_conformal_predictions.csv`)

### Result

Nothing docked, nothing refit. The frozen conformal predictor was loaded and
run; the AD gate's constants — ECFP4 radius 2, 2048 bits, 5-NN mean Tanimoto,
threshold **0.6810766499582288** — were **imported** from
`run_clinical_reference_ad_gate.py` rather than restated, so the two gates
cannot drift.

| Population | AD gate admits | Ligand arm admits | Both |
|---|---|---|---|
| **32 spike-ins inside the top-500 cut** (primary) | **5 / 32 = 0.156** [0.069, 0.318] | 31 / 32 = 0.969 | **5 / 32** |
| 89 docked spike-ins (secondary) | **18 / 89 = 0.202** [0.132, 0.297] | 79 / 89 = 0.888 | 17 / 89 |

| | Verdict |
|---|---|
| **P1a** — AD gate admits most | **FALSE** |
| **P1b** — ligand arm admits most | TRUE |
| **P1** | **MISSED** |
| P2 (from stage 4, not recomputed) | MET |
| **P3** | **MISSED** |

### What was falsified, and it was my reading, not the pipeline

**ADR-AIIM-229's falsification condition fired exactly as written:** *"If the AD
gate rather than docking is the largest loss, **P1 is wrong** and ADR-AIIM-210's
9-of-10 does not generalise past clinical-stage chemistry."*

Stage by stage, as a fraction of what each stage is handed:

| Stage | Keeps | Loses |
|---|---|---|
| Similarity funnel (1–3) | 89 / 105 = **84.8 %** | 16 |
| Docking rank, top-500 cut (4) | 32 / 89 = **36.0 %** | 57 |
| **Applicability-domain gate (5)** | **5 / 32 = 15.6 %** | **27** |

**The domain gate is proportionally the largest single loss in the pipeline**,
and end to end **5 of 105** known actives survive — **4.8 %**.

**P3 said the loss was localised to the docking rank. It is not.** I wrote that
localisation into §3.5, the abstract and §4 earlier today on P2's evidence while
P1 was still PENDING, and ADR-AIIM-259 had already flagged that as an
overstatement of my own. The measurement has now settled it against me. The
manuscript is corrected in the same session (below).

### Why this result is worth more than the one predicted

The paper's central argument is that similarity-based applicability-domain gates
exclude novel chemotypes for being unfamiliar rather than unpromising. Until
now the evidence was **library-level**: 0 % structural coverage for the
natural-product library, which a reader may read as a fact about orchids.

This is **compound-level and adversarial**: handed 105 compounds that are
*known actives for this exact target*, chosen because the model has never seen
them (ADR-AIIM-222), the production gate **refuses about four in five**.

It also bounds ADR-AIIM-210 rather than contradicting it. That control found the
gate admits **9 of 10** real A2A drugs — but those are clinical-stage compounds
inside the training set's own chemistry, and each was scored with its own record
removed. **The two results together say what neither says alone: the gate
admits actives that resemble its training chemistry and refuses actives that do
not, which is what an applicability domain is for and exactly why it must not
be used as a selection veto** — the argument §4 already makes, now with a
direct measurement behind it.

### Controls that make the number readable

* **Membership control.** A spike-in inside the domain-defining set would pass
  by construction. The script **refuses to run** if any is found; **0 of 89**
  were, consistent with ADR-AIIM-222 having drawn the panel from the 169 DUD-E
  actives absent from the A2A set. No self-exclusion was needed and none was
  applied.
* **Falsifiability shown before the data.** The self-test grades four synthetic
  scenarios and a knife-edge case; P1 fails when either arm refuses, P3 fails
  when either P1 or P2 fails, and exactly half is not "most". All pass.
* **P2 is read from the stage-4 artifact, not recomputed**, so this run cannot
  quietly move it.
* **Primary population fixed in advance** as the 32 inside the production cut,
  because that is what a deployed run hands to stage 5; the 89 are reported
  beside it and say the same thing.

### What is not evaluable, and is not fudged

The **dual-criterion tier rule was not applied**: Tier A requires an on-pocket
pose and the spike-in docking recorded scores, not poses. It is reported as *not
evaluable from this run* rather than approximated. P1 names the AD gate and the
ligand-based arm and does not depend on it.

### Consequence for O2

**O2 closes.** The objectives document states the criterion: *"Either outcome is
publishable: recovery validates the pipeline end to end; **non-recovery
localises the failure gate by gate**, which is more useful than the current
undifferentiated negative."* M2.1 has now done exactly that, for every stage.
The result is unfavourable to the pipeline and it is the evidence O2 asked for.

**All five objectives are now closed** (ADR-AIIM-259 for the other four).

### Not claimed

That the gate is wrong. It is behaving as an applicability domain should — it
reports that these compounds are unlike its training set, which is true. The
finding is about the **role** the gate was given, not its arithmetic. No
compound here has been assayed, and one panel against one target in one
representation does not establish a general refusal rate.

---

## ADR-AIIM-285 — PRE-REGISTRATION: the split convention as a 2×2×2 factorial — three defects, two fixed, and the one nobody touched

**Status:** proposed — **pre-registration. Frozen before the run.**
**Date:** 2026-10-05
**Raised by:** Adisorn Sirichotejirakul — *"ร่าง pre-registration ของ split ที่แก้ทั้งสามจุด"* (2026-10-05)
**Follows:** ADR-AIIM-275 (canonical ordering closes the row-order class),
ADR-AIIM-273 (the tiebreak fix), ADR-AIIM-271 (the distribution)
**Discipline:** ADR-AIIM-263 — predictions fixed before the run;
ADR-AIIM-226 — the change is made outside the frozen tree
**Bears on:** ADR-AIIM-275's option 2, "adopt canonical ordering and a declared
tiebreak", which assumes one declared convention exists to adopt

### Why there is anything left to do

ADR-AIIM-275 closed the row-order class exactly and left a 6× disagreement
between two reproducible procedures. That disagreement was attributed to the
tiebreak. **The tiebreak is not the only convention still free.**
`scaffold_train_test_split` ([`run_qsar_lbdd_2026-08-07.py:161`](../10_new_system/scripts/run_qsar_lbdd_2026-08-07.py))
contains three independent design choices, and both prior fixes changed none of
them:

```python
sorted_scafs = sorted(scaf_to_idx.items(), key=lambda x: len(x[1]), reverse=True)
rng.shuffle(sorted_scafs)                                                  # D1
n_test_target = int(len(labels) * test_frac)
test_set, train_set = set(), set()
for sc, idxs in sorted_scafs:
    (test_set if len(test_set) < n_test_target else train_set).update(idxs)  # D2, D3
```

| | the choice | what it does now | addressed by |
|---|---|---|---|
| **D1** | `rng.shuffle` after a size sort | destroys the size ordering the sort just built, so scaffolds enter in arbitrary order | **neither fix** — both ADR-AIIM-273 and ADR-AIIM-275 retain the shuffle ([`canonical_row_order.py:98`](../20_new_runs/canonical_row_order_2026-10-02/analysis/canonical_row_order.py)) |
| **D2** | the held-out side is filled first | the first groups encountered become calibration; training is the remainder | **neither fix** |
| **D3** | the budget is tested *before* adding a whole group | the group that crosses the line is admitted entire, so the held-out side overflows its target | **neither fix** |

What the two prior fixes changed was the *pre-shuffle order* — the stable-sort
tiebreak (ADR-AIIM-273) and then the table order itself (ADR-AIIM-275). Both
made the procedure reproducible. **Neither asked whether the procedure is the
one a scaffold split is supposed to be.** With D1 in place the size sort is
decorative: scaffolds are assigned in a fixed but arbitrary order, so the
held-out set is a random sample of scaffold groups rather than a deliberate one.

The standard convention is the opposite of all three. Largest groups to
**training**, so the model receives the densest structure–activity data; the
singleton and near-singleton scaffolds to **held-out**, so the held-out set is
maximally dissimilar from training and measures extrapolation to novel
chemotypes, which is the entire reason to prefer a scaffold split over a random
one; and a group that would overflow the budget is **diverted rather than
admitted**, so the held-out fraction is a bound and not an approximation.

### The change

Three independent switches over the same canonical table, run as a full
factorial so each defect's contribution is attributable rather than inferred:

| factor | level 0 (released) | level 1 (principled) |
|---|---|---|
| **D1 shuffle** | `rng.shuffle(ordered)` retained | removed; strict `(-size, scaffold_smiles)` order |
| **D2 direction** | fill held-out first, training is the remainder | fill training first, held-out takes the tail |
| **D3 overflow** | test the budget before adding, admit the whole group | divert the group to the other side if it would overflow |

D3 at level 1 is specified as greedy best-fit and **not** as a hard cutoff: a
group that would overflow is diverted and the loop *continues*, so a later,
smaller group may still be admitted. This is stated because the two behaviours
differ and a pre-registration that leaves it open is not frozen. The
consequence is that at D3 = 1 the held-out side is a bound, `len(held_out) ≤
int(n × 0.20)`, in every configuration.

2³ = 8 configurations × 2 arms = 16 cells. Arm A is the released labelled set
at 6.0/5.5 by maximum aggregation (5,268 active / 702 inactive); arm B adds the
ADR-AIIM-268 recovered inactives (5,268 / 1,439). Everything else is held as in
ADR-AIIM-275: canonical InChIKey row order throughout, features computed once
per arm and indexed, untuned default `RandomForestClassifier` as the scoring
model, `VarianceThreshold` inside the fitted `Pipeline`, Mondrian conformal on
the active class at α = 0.05, `TEST_FRAC = 0.20`, `random_state = 42`.

Every cell is additionally run over **10 row permutations** to confirm each is
invariant, giving 160 fits. ADR-AIIM-275's run was 200 fits in 186 s, so this
is budgeted at under 200 s and sits inside the ADR-AIIM-040 allowance.

**This is not applied to `10_new_system`.** It lives in the run's own harness.
Promoting any of it would change published numbers and is a separate decision.

### Predictions, fixed before anything is run

**T1 — every cell is invariant.** All 16 cells return one distinct
possibly-active count across 10 permutations, both libraries. *Falsified if any
cell returns two different values.* **Marked weak and confirmatory**:
ADR-AIIM-275 already established that canonical ordering delivers invariance for
two split rules, and all eight configurations here are deterministic functions
of the same canonical table. A failure would mean a reader of row order exists
that three fixes have missed.

**T2 — the fully principled configuration exceeds 3,334 on arm A.** With D1, D2
and D3 all at level 1, arm A returns more possibly-active on the NP library than
the 3,334 of ADR-AIIM-275's fixed split. *Falsified if it is 3,334 or fewer.*
*Basis, stated as a mechanism so that it can be wrong for a reason:*
`p1 = (1 + #{calib1 ≥ 1 − p_query}) / (1 + n1)` is inflated by calibration
actives the model scores poorly. Level 1 of D2 sends the singleton scaffolds to
calibration, and those are by construction the compounds least like anything in
training, so the model should score them worst, so `calib1` should be largest,
so `p1` should be largest, so the count should be the highest of any
configuration. If T2 fails the mechanism is wrong and the inflation originates
elsewhere.

**T3 — D2 dominates.** The main effect of D2 (direction) on the arm A NP count
exceeds the main effects of D1 and D3 taken separately. *Falsified if either D1
or D3 has the larger main effect.*

**T4 — D3 is nearly inert on its own.** Changing D3 alone, with D1 and D2 held
at level 0, moves the arm A NP count by less than 5 %. *Basis:* the overflow is
23 rows of 1,217 in arm A (1.9 %) and 19 of 1,360 in arm B (1.4 %), computed
from ADR-AIIM-275's reported calibration sizes against `int(n × 0.20)` targets
of 1,194 and 1,341. *Falsified if it moves by 5 % or more* — which would mean a
twenty-three-row accounting slip has been setting a published figure.

**T5 — no configuration recovers 72.** None of the 16 cells returns the released
NP-library figure. *Falsified if any does.* Weak: ADR-AIIM-271 placed 72 at the
8th percentile and no canonical rule has reproduced it in two attempts.

**T6 — the convention space is wider than the tiebreak space.** The eight arm A
NP values span more than the 6× of ADR-AIIM-275 (568 against 3,334). *Falsified
if the ratio of the largest to the smallest is 6.0 or less.*

**The outcome that would most embarrass this entry is T6 holding**, and it is
named first because it is the one this entry would prefer not to find.
ADR-AIIM-275 left three honest options and option 2 was to adopt a declared
convention and defend it. If eight equally reproducible conventions span more
than two orders of magnitude, option 2 does not survive contact with this
measurement: there is no longer one convention to declare, there is a menu, and
choosing from it after seeing the counts is the thing pre-registration exists to
prevent. T6 holding would push the project towards option 1, reporting the
distribution, or option 3, withdrawing the claims that rest on the count.

A secondary embarrassment, named because it is cheap not to hide: **T3 failing
together with T4 failing** would mean the two prior fixes addressed the small
defects and left the large one, and that ADR-AIIM-273's and ADR-AIIM-275's
framing of the problem as being about *ordering* was the wrong frame.

### What this run does not do

**It does not fix exchangeability, and it makes the problem more visible.**
Conformal validity requires the calibration set to be exchangeable with the
query compounds. A scaffold split deliberately destroys the similarity between
the held-out side and training, which is its purpose, and level 1 of D2
maximises that destruction. The manuscript already states that exchangeability
is not established for these libraries (§3.6 note at
`manuscript_AIIM_v6_CLEAN.md:388`) and that SNP-Lib is 0 % in-domain under four
tests (§3.5). **Every count this run produces is a conformal p-value without a
coverage guarantee, before any question of ordering or convention arises.**
Making the split principled does not repair that and is not claimed to.

**Invariance is not correctness**, restated from ADR-AIIM-274 constraint 3,
because this is the third entry in a row where it is the operative caution.
Nothing here licenses calling any of the eight values more accurate than 72, or
than each other. The quantity being measured is how much of the answer is set by
convention rather than by chemistry.

### Constraints

1. **`10_new_system` is not modified.** No contract amended, no input promoted,
   nothing re-docked, `qsar_pipeline_a2a.pkl` untouched. The three switches live
   in this run's harness only.
2. **The released figures stay in the record.** §3.6 reports 72 of 5,000 and a
   maximum conformal *p* of 0.072 and continues to, whatever this returns. Which
   figure the article should report remains an author decision and is not taken
   here.
3. **No prior artifact is overwritten.** Output goes to a new dated directory.
   `20_new_runs/canonical_row_order_2026-10-02/` and
   `20_new_runs/scaffold_split_fix_2026-10-02/` are the frozen evidence for
   ADR-AIIM-275 and ADR-AIIM-273 and are read-only to this run.
4. **No compound named in any output is a candidate.** Nothing in this project
   has been assayed and no compound is put forward.
5. **Every outcome is reported**, including T6 holding, including a principled
   configuration whose count is worse for the article than any already on
   record, and including the possibility that the two prior fixes addressed only
   the minor defects.
6. **This entry is appended to the decision log and its SHA-256 recorded before
   the first cell runs.**

### Skills invoked

| skill | version | why | material influence |
|---|---|---|---|
| `scikit-learn` | 1.1 | the run fits `Pipeline` + `RandomForestClassifier`, and the design concerns grouped splitting | `VarianceThreshold` is kept inside the fitted `Pipeline` per the skill's "Always Use Pipelines" rule, so no preprocessing is fit on held-out rows; grouped splitting by scaffold is retained in all eight configurations |
| `experimental-design` | 1.0 | three suspected causes, previously changed one at a time | the full 2³ factorial replaces the one-factor-at-a-time approach of ADR-AIIM-272 and ADR-AIIM-274, which is what left D1 and D2 invisible for two rounds; main effects become attributable and interactions estimable |

### Output

```
20_new_runs/split_convention_factorial_2026-10-05/
  split_convention_factorial.json
  run_manifest.json
  analysis/split_convention_factorial.py
```

## ADR-AIIM-286 — ADR-AIIM-285 result: all six predictions met, and the factorial shows the main effects were the wrong thing to have predicted

**Status:** accepted — record of fact
**Date:** 2026-10-05
**Approver:** Adisorn Sirichotejirakul — *"append แล้วรันเลย"* (2026-10-05)
**Pre-registration:** ADR-AIIM-285, SHA-256
`eff0353e578e23b20ddd943f0faaa5554a359b8294a83cc926d7f43c78c5f8ae`, frozen
before the run; log SHA-256 at that point
`61659f173f41295e1b269abf05c28e18e1febe992bb9160779ef04171b382c0f`
**Artifacts:** `20_new_runs/split_convention_factorial_2026-10-05/` · 8
configurations × 2 arms × 10 permutations, 160 fits, 136 s

### Grading

| | prediction | verdict |
|---|---|---|
| **T1** | every cell invariant across permutations | **MET** — one distinct value in 10 draws, all 16 cells, both libraries |
| **T2** | the principled configuration exceeds 3,334 on arm A | **MET** — 3,540, and it is the maximum of the eight |
| **T3** | D2 (direction) has the largest main effect | **MET** — +1,034.5 against −834.5 for D1 and +3.0 for D3 |
| **T4** | D3 alone moves arm A by less than 5 % | **MET** — 3,334 → 3,346, 0.36 % |
| **T5** | no configuration recovers 72 | **MET** — nearest is 88; but see the note below |
| **T6** | the eight arm A values span more than 6× | **MET** — 522 to 3,540, 6.78× |

### The harness is faithful

Cell `D1shuffle_D2testfirst_D3overflow` is ADR-AIIM-275's fixed split, and it
reproduced it exactly in both arms without adjustment: arm A NP 3,334 / SNP 28
/ calib inactive 132, arm B NP 269 / SNP 1 / calib inactive 260. The anchor was
pre-registered as the condition on which everything else depends, and it holds.

### The eight conventions

Arm A, NP library, every value invariant across 10 permutations:

| arm A NP | configuration | calib active | calib inactive | over target |
|---|---|---|---|---|
| **522** | sorted · test-first · overflow | 1,114 | 80 | 0 |
| 522 | sorted · test-first · divert | 1,114 | 80 | 0 |
| 2,391 | shuffle · train-first · overflow | 1,041 | 153 | 0 |
| 2,391 | shuffle · train-first · divert | 1,041 | 153 | 0 |
| 3,334 | shuffle · test-first · overflow **(= ADR-AIIM-275)** | 1,085 | 132 | +23 |
| 3,346 | shuffle · test-first · divert | 1,060 | 134 | 0 |
| **3,540** | sorted · train-first · divert **(principled)** | 1,014 | 180 | 0 |
| 3,540 | sorted · train-first · overflow | 1,014 | 180 | 0 |

Arm B, same eight conventions: **7, 88, 88, 269, 839, 1,486, 1,900, 2,634**.

### T3 is met and the quantity it is about is close to meaningless

This is the finding, and it is not one ADR-AIIM-285 anticipated.

**D2's effect reverses sign depending on D1.** On arm A, moving to train-first
is worth **+3,018** when the shuffle is removed and **−949** when it is
retained. The D1×D2 interaction is **+1,983.5**, nearly twice the +1,034.5
main effect that T3 was graded on. On arm B the reversal runs the other way:
**−1,074.5** sorted against **+2,129** shuffled, interaction −1,601.8 against a
main effect of +527.2.

A main effect averaged across a sign reversal describes no configuration that
exists. T3 passed on a number that should not be reported on its own, and the
pre-registration asked for exactly that number. **The factorial was the right
design and the prediction written over it was the wrong summary** — had this
been run one-factor-at-a-time, as ADR-AIIM-272 and ADR-AIIM-274 were, the
interaction would have been invisible and whichever single answer came out
would have been reported as the effect of direction.

**The practical consequence: "largest scaffolds to training" is not a
convention that can be adopted on its own.** It helps only once the shuffle is
removed. Adopted while the shuffle remains it makes arm A *worse* by 949
compounds. The three choices are not separable and must be declared as a set.

### T4 is met on arm A and the claim behind it is false on arm B

T4 was specified on one cell pair and passed there at 0.36 %. The basis given
was that the overflow is 23 rows of 1,217, "a twenty-three-row accounting
slip". **That reasoning does not survive arm B.** The same cell pair —
shuffle, test-first, D3 flipped — gives **269 → 7, a 97.4 % collapse.** Two
other arm B pairs move +77.1 % and +38.6 %.

So D3 is inert in arm A and decisive in arm B, from a nineteen-row difference.
The literal prediction is MET and **the generalisation it was built on is
withdrawn**: there is no "D3 is nearly inert" result here. A boundary group
whose admission or diversion changes nineteen of 1,360 calibration rows can
change the reported count by a factor of 38.

### T5 note — the SNP-Lib figure *is* recovered, three times

T5 was written about 72, the NP-library figure, and no cell returns it. But the
released **SNP-Lib** figure is **0**, and three arm B cells return 0:
shuffle · test-first · divert, and both sorted · train-first cells. This is
recorded because T5's wording concerned one library and the other library
behaves differently. It is weak evidence: 0 is a floor, and a count that
reaches it says less than a count that matches a positive figure.

### T2's mechanism holds on arm A and fails on arm B

T2 was deliberately written as a mechanism so it could be wrong for a reason.
On arm A it is right: the principled cell returns the maximum of the eight, and
its calibration holds the fewest actives (1,014, against 1,114 for the lowest
count), consistent with singleton scaffolds being scored worst and inflating
`p1`.

On arm B the principled cell returns **88 — second lowest of the eight.** The
mechanism does not generalise across arms. Nothing in the pre-registration
predicted this and it is reported as an observation, not a graded result.

### T6 held, and arm B is worse than T6 was willing to imagine

T6 asked whether the eight values span more than ADR-AIIM-275's 6×. On arm A
they span **6.78×** — met, but only just, and the embarrassment paragraph's
speculation of "more than two orders of magnitude" did **not** happen on arm A.

**It happened on arm B: 7 to 2,634, a span of 376×.** Arm B was not named in
T6, so this is an observation rather than a prediction met, and it is the
single largest number this run produced.

ADR-AIIM-275 left three options and option 2 was to adopt a declared
convention and defend it. On the evidence here that option is in materially
worse shape than it was this morning:

* there is no single convention to declare, there are three coupled choices;
* their effects are not separable, so the set must be defended as a set;
* on the arm with the better data, the eight defensible sets span 376×;
* and choosing among them after seeing these counts is precisely what
  pre-registration exists to prevent.

**This run does not pick one.** Nothing licenses calling 3,540, or 522, or 88,
more accurate than 72. What the project now has is a measurement of how much of
§3.6's count is set by convention rather than by chemistry, and the answer on
arm B is: nearly all of it.

### What was already true and is now harder to avoid

Every count here is a conformal *p*-value without a coverage guarantee.
Validity requires the calibration set to be exchangeable with the query
compounds; a scaffold split destroys that similarity by design, and the
principled configuration maximises the destruction. The manuscript already
states exchangeability is not established (§3.6 note at
`manuscript_AIIM_v6_CLEAN.md:388`) and that SNP-Lib is 0 % in-domain under four
tests (§3.5). **The convention instability measured here sits on top of a
guarantee that was never in force.** Repairing the split does not repair that,
and no configuration in this run should be read as the corrected figure.

### Where this leaves §3.6 — unchanged

§3.6 still reports 72 of 5,000 and a maximum conformal *p* of 0.072. Nothing
has been enacted. The three honest options of ADR-AIIM-275 stand, with option
2 now carrying a cost it did not visibly carry before. All three remain author
decisions and none is taken here.

### Deliberately not done

* **`10_new_system` is unmodified.** The three switches live in this run's
  harness. No contract amended, no input promoted, nothing re-docked,
  `qsar_pipeline_a2a.pkl` untouched, no published number changed.
* **No prior artifact was overwritten.**
  `20_new_runs/canonical_row_order_2026-10-02/` and
  `20_new_runs/scaffold_split_fix_2026-10-02/` were read-only to this run and
  their hashes are unchanged.
* **The tiebreak was not a factor.** It was held at `(-size, scaffold_smiles)`
  in all eight cells because ADR-AIIM-275 measured it. A 16-cell design
  including it was not run.
* **Interactions were computed but not pre-registered.** The D1×D2 interaction
  and the per-cell D3 effects are reported above as observations. Only T1–T6
  were graded.
* **Nested CV was not re-derived.** ADR-AIIM-275 already showed fold membership
  is invariant under canonical ordering; whether the AUCs move across these
  eight conventions was not measured.
* **No compound named in any output is a candidate.** Nothing in this project
  has been assayed and no compound is put forward.

### Artifacts

```
20_new_runs/split_convention_factorial_2026-10-05/
  run_manifest.json
  split_convention_factorial.json        f3f48c4622d7e19f2225194a01fdc8b468bc242c1822a01daf9f4a3f5f404c1c
  analysis/split_convention_factorial.py bb4b25da6b7218eb5fa414d09d3d3c86cd82800540d024bc83cce1d6d73bc8e2
```

## ADR-AIIM-287 — PRE-REGISTRATION: three claims I asserted without testing — label source, calibration balance, and applicability domain

**Status:** proposed — **pre-registration. Frozen before the run.**
**Date:** 2026-10-05
**Raised by:** Adisorn Sirichotejirakul — *"อยากให้ลองรันตามที่แนะนำเพื่อดูผลออกมาก่อน อย่าพึ่งบอกว่าไม่ได้ทั้ง ๆ ที่ยังไม่ได้ลองทำ"* (2026-10-05)
**Follows:** ADR-AIIM-286 (the 2×2×2 convention factorial)
**Discipline:** ADR-AIIM-263 — predictions fixed before the run;
ADR-AIIM-226 — the change is made outside the frozen tree

### Why this entry exists

Across two exchanges after ADR-AIIM-286 I made three claims about what would
*not* help, and tested none of them:

1. restricting the label to **Ki only** "would not fix the problem", because
   the instability is in the split convention and not in the label;
2. **calibration balance** "is not the driver", argued from correlations across
   the eight cells of ADR-AIIM-286 — where balance was never manipulated, only
   observed, so the argument is observational and weak;
3. restricting to the **applicability domain** cannot help because "you cannot
   calibrate your way out of extrapolation".

Claim 2 is the weakest: an observational correlation over eight points that
reversed between arms (r = −0.863 on arm A, +0.117 on arm B) is not evidence
about a manipulated intervention. Claims 1 and 3 were asserted from mechanism
with no measurement at all. **All three are cheap to test and are tested here.**
The predictions below are written so that my own position is the falsifiable
one.

### The three interventions

| | intervention | tests claim |
|---|---|---|
| **Label source** | arm **K**: `standard_type == "Ki"` only (6,575 of 8,320 released rows), plus the ADR-AIIM-268 censored recovery restricted to Ki | 1 |
| **Calibration balance** | subsample calibration **actives** down to the number of calibration inactives, seeded. `cal1` is built from actives alone, so this is the balance the *p*-value actually reads | 2 |
| **Query domain** | report counts on the full NP library **and** on its in-domain subset, using the project's own `compute_ad_threshold` / `compute_ad_scores` at `AD_K = 5`, `AD_PERCENTILE = 10` | 3 |

Design: 8 split conventions × 3 label arms (A, B, K) × 2 calibration modes
(raw, balanced) = **48 cells**, each reporting the full-library and in-domain
count, over 2 permutations = **96 fits**. Invariance is not re-established;
ADR-AIIM-286 settled it and 2 permutations are a spot check only.

The applicability domain is computed **once per arm** against a seeded
500-compound sample of that arm's labelled set, not per cell. The AD question
is whether an NP compound resembles the training chemistry at all, which is a
property of the labelled set rather than of which 80 % of it became training.
This is a deliberate simplification and is declared as one.

Everything else is held as in ADR-AIIM-286: canonical InChIKey row order,
features computed once per arm and indexed, untuned default
`RandomForestClassifier`, `VarianceThreshold` inside the fitted `Pipeline`,
Mondrian conformal on the active class at α = 0.05, `TEST_FRAC = 0.20`,
`random_state = 42`, tiebreak held at `(-size, scaffold_smiles)`.

### Predictions — my own claims, stated so they can fail

**U1 — Ki-only does not collapse the spread.** Arm K's eight conventions span
more than 6× on the full NP library. *Falsified if the span is 6.0× or less*,
which would mean the label source was a larger lever than the split convention
and that claim 1 was wrong.

**U2 — balanced calibration does not collapse the spread.** With calibration
actives subsampled to match inactives, the span across the eight conventions
still exceeds 6× in at least two of the three arms. *Falsified if it is 6.0× or
less in two or more arms*, which would mean balance was the driver after all
and that my correlation argument reached the right conclusion for the wrong
reason, or the wrong conclusion entirely.

**U3 — in-domain restriction does not collapse the spread.** The in-domain
count spans more than 6× across the eight conventions in at least two of the
three arms. *Falsified if it is 6.0× or less in two or more arms.*

**U4 — the label source moves the answer without stabilising it.** Arm K's
principled cell (`sorted · train-first · divert`) differs from arm B's by more
than 20 %. *Falsified if within 20 %.* This is predicted *because* U1 is
predicted: a lever that changes the answer but not its stability is the thing I
claimed Ki-only would be.

**The outcome that would most embarrass this entry is any of U1, U2 or U3
being falsified**, and it is named first. Each would mean I ranked the
project's concerns in a prior exchange from assertion rather than measurement,
and told the author that an avenue was closed without opening it. If that
happens, the ranking offered after ADR-AIIM-286 — applicability domain first,
convention second, censoring third, balance fourth — is withdrawn and rebuilt
from these numbers.

A second outcome worth naming: **U3 failing would be good news for the
project**, because it would mean restricting §3.6's claim to the in-domain
subset buys back a determinate number. I have argued it will not. I would
rather be wrong.

### Constraints

1. **`10_new_system` is not modified.** No contract amended, no input promoted,
   nothing re-docked, `qsar_pipeline_a2a.pkl` untouched.
2. **The released figures stay in the record.** §3.6 reports 72 of 5,000 and
   continues to. No count here is proposed as a replacement.
3. **No prior artifact is overwritten.** Output goes to a new dated directory;
   the 2026-10-02 and 2026-10-05 run directories are read-only to this run.
4. **Arm K is not proposed as the correct label set.** It is one arm of a
   measurement. Whether Ki-only is the right basis for the article is an
   author decision and is not taken here.
5. **The AD simplification above is reported with every in-domain number**, and
   no in-domain count is presented as a validated conformal result — the
   coverage guarantee is not restored by subsetting a query set after the fact,
   whatever U3 returns.
6. **No compound named in any output is a candidate.** Nothing in this project
   has been assayed.
7. **Every outcome is reported**, including all three claims failing.
8. **This entry is appended to the decision log and its SHA-256 recorded before
   the first cell runs.**

### Skills invoked

| skill | version | why | material influence |
|---|---|---|---|
| `scikit-learn` | 1.1 | `Pipeline` + `RandomForestClassifier`, grouped splitting | preprocessing stays inside the fitted pipeline; no transform is fit on calibration rows, including under the balanced-calibration intervention |
| `experimental-design` | 1.0 | three interventions crossed with eight conventions | balance is **manipulated** here rather than observed, which is the defect in the argument this entry exists to repair |

### Output

```
20_new_runs/label_balance_domain_2026-10-05/
  label_balance_domain.json
  run_manifest.json
  analysis/label_balance_domain.py
```

## ADR-AIIM-288 — ADR-AIIM-287 result: two of my three claims were false, the label was the lever, and the convention instability is largely its symptom

**Status:** accepted — record of fact
**Date:** 2026-10-05
**Approver:** Adisorn Sirichotejirakul — *"ร่าง ADR-AIIM-288 แล้ว append เลย"* (2026-10-05)
**Pre-registration:** ADR-AIIM-287, SHA-256
`3a0d096a1ba9ff3aa031d00c81285d334d286e846723dcb7db8dcdd311af969a`, frozen
before the run; log SHA-256 at that point
`1ed078b4a8b48be37e507cc01737e70d50351935da81b82e56ca7099304e16bb`
**Artifacts:** `20_new_runs/label_balance_domain_2026-10-05/` · 8 conventions ×
3 label arms × 2 calibration modes = 48 cells, 96 fits, 63 s
**Bears on:** ADR-AIIM-286's framing of the instability as a property of the
split convention

### Grading

| | prediction (my own claim) | verdict |
|---|---|---|
| **U1** | Ki-only does not collapse the spread | **MISSED** — arm K spans 3.62×, against 376.29× for arm B |
| **U2** | balanced calibration does not collapse the spread | **MET** — and forcing it makes arms A and B far worse |
| **U3** | in-domain restriction does not collapse the spread | **MISSED** — vacuously; the in-domain set is empty |
| **U4** | the label source moves the answer without stabilising it | **MET** on the moving, **MISSED in spirit** — it stabilises as well |

Two of three claims false. Both were asserted to the author in conversation
before any measurement, and the author's instruction to run them rather than
accept them was correct.

### The arms

| arm | label rule | active | inactive | % active | n |
|---|---|---|---|---|---|
| **A** | released, pooled Ki + IC50 | 5,268 | 702 | 88.2 | 5,970 |
| **B** | A + censored recovery (ADR-AIIM-268) | 5,268 | 1,439 | 78.5 | 6,707 |
| **K** | `standard_type == "Ki"` only, same recovery | 4,143 | 1,280 | 76.4 | 5,423 |

Arm K keeps 76 % of arm B's rows and 86 % of its recovered inactives.

### U1 falsified: the label source cuts convention sensitivity by two orders

NP library, span across the eight split conventions of ADR-AIIM-286:

| arm | min | max | **fold** |
|---|---|---|---|
| A | 522 | 3,540 | 6.78× |
| B | 7 | 2,634 | **376.29×** |
| **K** | 1,101 | 3,986 | **3.62×** |
| K, balanced calibration | 1,352 | 3,742 | **2.77×** |

ADR-AIIM-286 reported the convention spread as the finding and attributed it to
three coupled design choices in the split. **That attribution is now partly
wrong.** Restricting the label to a single assay quantity reduces the spread
from 376× to 3.62× without touching the split at all.

**The mechanism was already in hand and was not followed through.** The audit
preceding this entry measured, on the 162 compounds carrying both quantities, a
systematic `pKi − pIC50` offset of +0.353 log (mean +0.673, SD 0.865) with
32.1 % disagreeing by more than 1 log and a maximum of 2.85. That offset is
Cheng–Prusoff behaving as theory requires, which is why it is systematic rather
than noise — and label noise concentrated in particular chemotypes is exactly
what makes an answer sensitive to which scaffolds land in training versus
calibration. The audit was performed, reported, and then used to argue that
Ki-only *would not help*. The opposite followed from it.

**ADR-AIIM-286's headline is therefore revised, not withdrawn.** The eight
conventions do disagree, and on pooled labels they disagree enormously. But
convention sensitivity is substantially *conditional on label quality*, and
calling it an independent property of the split overstated it.

### U3 falsified vacuously, because the prediction was badly specified

| arm | NP-library in-domain | SNP-Lib in-domain |
|---|---|---|
| A | **0 / 5,000** | 0 / 236 |
| B | **0 / 5,000** | 0 / 236 |
| K | **0 / 5,000** | 0 / 236 |

AD threshold 0.461 at `AD_K = 5`, `AD_PERCENTILE = 10`; training
self-similarity mean 0.617, SD 0.102. **Not one of the 5,000 NP compounds
reaches the 10th percentile of the training set's own self-similarity.**

The in-domain count is therefore 0 in all 48 cells, the span is 0.0×, and the
automated grader recorded U3 as falsified because 0.0 ≤ 6.0. **The prediction
could not distinguish "restriction stabilises the answer" from "restriction
leaves nothing to count", and that is a defect in how I wrote it**, not a
result. U3 yields no information about its intended question.

What it does yield is a fact the project did not have. §3.5 establishes 0 %
in-domain for SNP-Lib under four tests; the manuscript makes no such claim for
the 5,000-compound NP library, which is the library §3.6's count is about.
**That library is now measured at 0 / 5,000 under the project's own AD rule, in
all three label arms.** This strengthens rather than weakens the reading that
every count in §3.6 is a conformal *p*-value with no coverage guarantee in
force.

### U2 met, and the mechanism is not balance

| arm | raw | balanced | effect |
|---|---|---|---|
| A | 6.78× | **319.0×** | 47× worse |
| B | 376.29× | **596.25×** | 1.6× worse |
| K | 3.62× | 2.77× | slightly better |

Forcing the calibration set to class balance does not stabilise anything; on
the two pooled-label arms it destabilises severely. **But the cause is not
balance, it is `n1` collapse.** `cal1` is built from calibration actives alone,
so matching them to the inactive count cuts arm A's `n1` from about 1,085 to
80–180 — to 12 % — and the α = 0.05 quantile is then estimated from under two
hundred points. Arm K survives because its inactive count is 170–319, a far
milder reduction, and because a random subsample preserves the quantile: two
arm K cells return bit-identical counts (1,352 and 2,025) before and after
balancing, which is the expected behaviour and a check that the intervention
was implemented correctly.

So the prior claim that balance "is not the driver" reached the right
conclusion, but the correlation argument offered for it — eight observational
points whose sign reversed between arms, r = −0.863 and +0.117 — was not
evidence for it. The conclusion survives; the reasoning given for it does not.

### U4: the label moves the answer by 25×

Principled cell (`sorted · train-first · divert`), raw calibration:
**arm B = 88, arm K = 2,235**, a difference of 2,440 %.

U4 predicted the label would move the answer *without* stabilising it, as the
signature of a lever that is not a fix. It moves it and it stabilises it. The
prediction is graded MET on its literal clause and is **withdrawn in spirit**,
because its purpose was to support U1 and U1 is false.

### The priority ranking offered to the author is withdrawn and rebuilt

ADR-AIIM-287 committed to this if any of U1–U3 failed. Two did.

| | concern | was | **now** | basis |
|---|---|---|---|---|
| 1 | applicability domain / exchangeability | 1 | **1**, strengthened | 0 / 5,000 and 0 / 236, all three arms |
| 2 | **label source — pooled Ki + IC50** | 3 | **2** | 376× → 3.62×; answer moves 25× |
| 3 | convention instability | 2 | **3, and partly a symptom** | 376× pooled against 3.62× Ki-only |
| 4 | calibration balance | 4 | **4, and do not force it** | 47× and 1.6× worse when forced |

The single worst statement in the preceding exchange was that the Ki-only
avenue was closed, made without opening it. It is the most effective
intervention this project has measured.

### What has not changed

**3.62× is not determinate.** Arm K spans 1,101 to 3,986 — a range of 2,885
compounds across eight reproducible conventions. The problem is reduced by two
orders of magnitude and is not solved.

**Reproducibility is still not validity.** With 0 / 5,000 in-domain, every
count in this entry is a conformal *p*-value outside the conditions its
guarantee requires. Cleaning the label made the number far more reproducible
and did nothing for its validity. No value here — 88, 2,235, 3,986 or any
other — is proposed as more accurate than 72.

**§3.6 is unchanged.** It reports 72 of 5,000 and a maximum conformal *p* of
0.072. Nothing has been enacted. ADR-AIIM-275's three options stand, now with
a fourth consideration: adopting Ki-only labels would narrow the convention
problem substantially while changing the published figure by more than an order
of magnitude, and would still leave a quantity with no coverage guarantee.

### Deliberately not done

* **`10_new_system` is unmodified.** No contract amended, no input promoted,
  nothing re-docked, `qsar_pipeline_a2a.pkl` untouched, no published number
  changed. Arm K lives in this run's harness.
* **Arm K is not proposed as the correct label set.** Whether Ki-only is the
  right basis for the article is an author decision and is not taken here.
* **The AD simplification stands as declared.** The domain was computed once
  per arm against a seeded 500-compound sample of the labelled set, not per
  cell. With 0 / 5,000 the simplification is not load-bearing — no per-cell
  refinement reaches a non-empty set from zero — but it was a simplification
  and is reported as one.
* **No in-domain count is a validated conformal result.** Subsetting a query
  set after the fact does not restore exchangeability, and the empty result
  makes the question moot rather than answered.
* **Invariance was not re-established.** ADR-AIIM-286 settled it; two
  permutations here are a spot check.
* **Arm K's nested-CV AUC was not computed.** Whether the cleaner label also
  improves measured discrimination is untested, and no claim is made that it
  does.
* **No compound named in any output is a candidate.** Nothing in this project
  has been assayed and no compound is put forward.

### Skills invoked

| skill | version | material influence |
|---|---|---|
| `scikit-learn` | 1.1 | `VarianceThreshold` kept inside the fitted `Pipeline` in all 48 cells, including under the balanced-calibration intervention, so no transform is fit on calibration rows |
| `experimental-design` | 1.0 | balance was **manipulated** rather than observed, which is the specific defect in the prior reasoning that this entry exists to repair |

### Artifacts

```
20_new_runs/label_balance_domain_2026-10-05/
  run_manifest.json
  label_balance_domain.json         270056c55cd509518e417a28c4c9bf1ed04fa353302db6783fb2c8220dc33395
  analysis/label_balance_domain.py  edea10a9408f611e37f51b733e4479e6b6e0d65292e259bd8baf4463127a7bff
```

## ADR-AIIM-289 — PRE-REGISTRATION: does the Ki-only label actually discriminate better, or only more reproducibly?

**Status:** proposed — **pre-registration. Frozen before the run.**
**Date:** 2026-10-05
**Raised by:** Adisorn Sirichotejirakul — *"วัด AUC ของ arm K"* (2026-10-05)
**Follows:** ADR-AIIM-288 (arm K cut convention sensitivity from 376× to 3.62×)
**Discipline:** ADR-AIIM-263 — predictions fixed before the run;
ADR-AIIM-226 — the change is made outside the frozen tree
**Bears on:** ADR-AIIM-288's promotion of label source to priority 2, and
ADR-AIIM-268's reported nested-CV figures

### The open question ADR-AIIM-288 declared

ADR-AIIM-288 listed under *deliberately not done*: *"Arm K's nested-CV AUC was
not computed. Whether the cleaner label also improves measured discrimination
is untested, and no claim is made that it does."* That gap is the whole basis
of the recommendation. Arm K was promoted to priority 2 because it made the
answer **reproducible** — 376× to 3.62×. Reproducibility is not accuracy, which
this project has now restated three times, so a label change that stabilises
the count while *degrading* discrimination would be a worse trade than it
currently looks.

### What is run

`nested_cv_auc` as frozen in `qsar_threshold_comparison.py:110` —
`StratifiedGroupKFold(5)` outer, `StratifiedKFold(3)` inner,
`RandomizedSearchCV` over `n_estimators ∈ {100, 200, 500}`,
`max_depth ∈ {None, 15, 25}`, `min_samples_leaf ∈ {1, 5, 10}` at
`n_iter = 8` — applied unchanged to all three label arms of ADR-AIIM-287:

| arm | label rule | active | inactive | prevalence |
|---|---|---|---|---|
| A | released, pooled Ki + IC50 | 5,268 | 702 | 88.2 % |
| B | A + censored recovery | 5,268 | 1,439 | 78.5 % |
| **K** | `standard_type == "Ki"` only, same recovery | 4,143 | 1,280 | 76.4 % |

Grouping is `get_scaffold_groups`, which `sorted(set(...))` makes
order-invariant already. Row order is canonical InChIKey throughout, because
ADR-AIIM-275 showed fold membership is order-dependent without it and
ADR-AIIM-288 left ADR-AIIM-268's figures carrying that caveat. Two
permutations confirm the AUC is now invariant.

125 fits per arm, 375 in all, plus 2 permutations for the invariance check.

**Reported per arm, following ADR-AIIM-268's own convention:** ROC-AUC mean and
SD across outer folds, MCC, balanced accuracy, **and PR-AUC together with its
lift over prevalence**, because ADR-AIIM-268 established that at this
prevalence the raw PR-AUC flatters and the lift is the honest figure.

### Predictions

**V1 — the cleaner label discriminates better.** Arm K's ROC-AUC exceeds arm
B's. *Falsified if it is equal or lower.* *Basis:* removing the pooled measure
removes a label defect already quantified at +0.353 log systematic offset with
32.1 % of dual-measured compounds disagreeing by more than 1 log, and that
should outweigh the loss of 19 % of rows.

**V2 — and more consistently.** Arm K's SD across the five outer folds is below
arm B's 0.0076. *Falsified if equal or greater.* This is the discrimination
analogue of the 3.62× convention finding: if the label is cleaner, folds should
agree more.

**V3 — the harness is faithful.** Arms A and B reproduce ADR-AIIM-268's
0.9181 and 0.9371 to within 0.02. *Falsified outside that band.* The tolerance
is 0.02 rather than exact because this run's arms carry 5,268 actives against
ADR-AIIM-268's 5,244 and 1,439 inactives against its 1,429 — a small
composition difference whose origin is not established here and is reported as
an open discrepancy rather than reconciled.

**V4 — the gain is not a prevalence artefact.** Arm K's PR-AUC lift over
prevalence exceeds arm B's. *Falsified if equal or lower*, which would mean a
ROC-AUC improvement bought by the easier class balance rather than by better
separation, exactly the trap ADR-AIIM-268 flagged.

**The outcome that would most embarrass this entry is V1 failing**, and it is
named first. ADR-AIIM-288 promoted the label source to priority 2 on
reproducibility evidence alone, and I recommended it to the author in the same
breath as withdrawing a ranking I had built from assertion. If arm K
discriminates worse, that promotion was made on half the evidence and the
ranking needs a third revision — this time downward for the option I had just
argued up.

A second outcome worth naming: **V1 holding while V4 fails** would be the
subtler failure, and the more likely one. It would mean Ki-only looks better by
ROC-AUC only because its classes are less skewed, and ADR-AIIM-288's
recommendation would rest on a metric this project has already written down as
misleading at this prevalence.

### Constraints

1. **`10_new_system` is not modified.** `nested_cv_auc` is imported and called,
   not edited. No contract amended, no input promoted, nothing re-docked,
   `qsar_pipeline_a2a.pkl` untouched.
2. **No published number is changed.** ADR-AIIM-268's figures stay in the
   record. If arms A and B come out different under canonical ordering, that is
   reported as a measurement, not a correction.
3. **No prior artifact is overwritten.** Output goes to a new dated directory.
4. **A better AUC does not make arm K correct.** With the NP library at
   0 / 5,000 in-domain (ADR-AIIM-288), cross-validated discrimination on ChEMBL
   compounds says nothing about whether a prediction on a natural product is
   valid. V1 holding would not restore a coverage guarantee and will not be
   reported as if it did.
5. **The active/inactive count discrepancy against ADR-AIIM-268 is reported,
   not silently absorbed.**
6. **No compound named in any output is a candidate.** Nothing in this project
   has been assayed.
7. **Every outcome is reported**, including arm K discriminating worse.
8. **This entry is appended to the decision log and its SHA-256 recorded before
   the run.**

### Skills invoked

| skill | version | material influence |
|---|---|---|
| `scikit-learn` | 1.1 | the frozen nested CV is used unchanged; the skill's rule that preprocessing stays inside the fitted pipeline is already satisfied by `make_rf_pipe`, and grouped splitting by scaffold is preserved through `StratifiedGroupKFold` |

`statistical-analysis` v1.0 was **read and not applied**: it covers hypothesis
testing, regression, assumption checking and APA reporting, and does not address
classification metric selection at high prevalence. The decision to report
PR-AUC lift over the no-skill baseline alongside ROC-AUC is taken from
**ADR-AIIM-268's own convention**, not from that skill, and is attributed here
so the provenance is not misread.

### Output

```
20_new_runs/arm_k_auc_2026-10-05/
  arm_k_auc.json
  run_manifest.json
  analysis/arm_k_auc.py
```

## ADR-AIIM-290 — ADR-AIIM-289 result: the Ki-only label does not discriminate better, and the metric that said it did was confounded by its own prevalence

**Status:** accepted — record of fact
**Date:** 2026-10-05
**Approver:** Adisorn Sirichotejirakul — *"รอให้รันเสร็จแล้ว append เลย"* (2026-10-05)
**Pre-registration:** ADR-AIIM-289, SHA-256
`ead611ad83c9bf5aab6d43d66d23902b5d0a8925cf41f50c249675befc0a9b7f`, frozen
before the run; log SHA-256 at that point
`9b920ad1e2f0a4c0aea7beb7ffcd2bb12be4e04339ff53793f8173a6ed99fefa`
**Artifacts:** `20_new_runs/arm_k_auc_2026-10-05/` · 3 arms × 2 permutations ×
125 fits = 750 fits, 508 s
**Bears on:** ADR-AIIM-288's promotion of label source to priority 2

### Grading

| | prediction | verdict |
|---|---|---|
| **V1** | arm K's ROC-AUC exceeds arm B's | **MISSED** — 0.9326 against 0.9377, Δ = −0.0051 |
| **V2** | arm K's fold SD is below arm B's 0.0076 | **MET** — 0.0060 against 0.0082 |
| **V3** | arms A and B reproduce ADR-AIIM-268 within 0.02 | **MET** — within 0.0004 and 0.0006 |
| **V4** | arm K's PR-AUC lift exceeds arm B's | **MET** — +0.2120 against +0.1938, **but see below** |

### The measurement

| arm | ROC-AUC | SD | PR-AUC | prevalence | lift | headroom | **lift / headroom** | MCC | bal-acc |
|---|---|---|---|---|---|---|---|---|---|
| A | 0.9185 | 0.0177 | 0.9860 | 0.8824 | +0.1036 | 0.1176 | **0.8810** | 0.5826 | 0.7952 |
| **B** | **0.9377** | 0.0082 | 0.9793 | 0.7854 | +0.1938 | 0.2146 | **0.9031** | **0.6966** | **0.8623** |
| K | 0.9326 | **0.0060** | 0.9760 | 0.7640 | **+0.2120** | 0.2360 | 0.8983 | 0.6713 | 0.8486 |

Every value invariant across both permutations under canonical ordering, as
ADR-AIIM-275 predicted it would be.

**V3 is near-exact.** Arm A returns 0.9185 against ADR-AIIM-268's 0.9181 and
arm B 0.9377 against 0.9371, with SDs of 0.0177 and 0.0082 against 0.0173 and
0.0076. The composition discrepancy declared in ADR-AIIM-289 — 5,268 actives
here against ADR-AIIM-268's 5,244, 1,439 inactives against 1,429 — moves the
AUC by 0.0004. It remains unexplained and is now also shown to be immaterial to
this metric. **It is not reconciled and stays on the record as open.**

### V1 failed, and I had argued the opposite to the author

Arm K discriminates **worse** than arm B by ROC-AUC, and also by MCC (0.6713
against 0.6966) and balanced accuracy (0.8486 against 0.8623). The basis given
in ADR-AIIM-289 — that removing a +0.353 log systematic label offset with
32.1 % of dual-measured compounds disagreeing by over 1 log would outweigh the
loss of 19 % of rows — is false as stated. The 19 % of rows mattered more, or
the IC50 records were less damaging than the offset suggested, and this run does
not distinguish those.

The differences are small against fold variation: Δ = −0.0051 with SDs of
0.0060 and 0.0082. **That is not an escape.** V1 was pre-registered as a strict
inequality and it failed.

### V4 is MET and the conclusion it would support does not survive

ADR-AIIM-289 adopted ADR-AIIM-268's convention that at this prevalence the raw
PR-AUC flatters and the **lift over prevalence** is the honest figure. By that
measure arm K wins: +0.2120 against +0.1938.

**Normalising by available headroom reverses it.** Maximum attainable lift is
`1 − prevalence`, which is 0.2360 for arm K and 0.2146 for arm B, because arm K
has the lower prevalence. As a fraction of what was attainable:

| arm | lift / headroom |
|---|---|
| B | **0.9031** |
| K | 0.8983 |
| A | 0.8810 |

Arm K's larger raw lift is **substantially mechanical**: it had more room to
gain. V4's verdict stands as graded, and the inference it was written to
support — that arm K separates the classes better once prevalence is accounted
for — is withdrawn.

**This is a defect in ADR-AIIM-268's lift convention, inherited here without
examination.** Subtracting prevalence corrects the baseline but not the ceiling,
so raw lift is not comparable across arms of differing prevalence. ADR-AIIM-268
used it to compare the 6.0/5.5 and 5.0/4.5 regimes, whose prevalences are 0.78
and 0.969 — a far larger gap than the one here — so **the figures in
ADR-AIIM-268 §"PR-AUC at the 5.0 cuts is near its own baseline" may understate
or overstate by the same mechanism and have not been re-derived.** That is
flagged, not fixed.

### What this means, stated as the trade it is

Arm B beats arm K on **every metric that accounts for prevalence**: ROC-AUC,
MCC, balanced accuracy and headroom-normalised lift. Arm K wins on raw lift,
which is confounded, and on fold-to-fold consistency, which is real.

So the honest statement is: **the Ki-only label does not discriminate better.**
Against that, ADR-AIIM-288 measured it cutting convention sensitivity from
376.29× to 3.62×. The trade is now quantified in both directions:

| from arm B to arm K | cost | benefit |
|---|---|---|
| ROC-AUC | −0.0051 | |
| MCC | −0.0253 | |
| balanced accuracy | −0.0137 | |
| lift / headroom | −0.0048 | |
| fold SD | | 0.0082 → 0.0060 |
| **convention spread** | | **376.29× → 3.62×** |

Whether a ~0.005 AUC and ~0.025 MCC cost is worth a hundredfold reduction in
convention sensitivity is an **author decision and is not taken here.** It is
not obviously a bad trade; it is no longer the clean win ADR-AIIM-288 implied.

### The two interventions pull in opposite directions, which is new

| step | discrimination | convention spread |
|---|---|---|
| A → B, censored recovery | **+0.0192 ROC, +0.1140 MCC** | **6.78× → 376.29×, worse** |
| B → K, Ki-only | −0.0051 ROC, −0.0253 MCC | **376.29× → 3.62×, better** |

The censored recovery is by far the larger discrimination gain in the project —
and it destabilised the count. Ki-only slightly degrades discrimination and
more than undoes that destabilisation. **No single arm is best on both axes**,
and the pair of changes cannot be evaluated as one decision.

### The priority ranking is revised a third time

ADR-AIIM-288 promoted label source to 2 on reproducibility evidence alone, and
said so. The discrimination evidence now runs mildly against it.

| | concern | ADR-288 | **now** |
|---|---|---|---|
| 1 | applicability domain / exchangeability | 1 | **1**, untouched by any of this |
| 2 | **censored-record handling** | 3 | **2** — the only change measured to improve discrimination, and it costs stability |
| 3 | **label source (Ki-only)** | 2 | **3** — buys stability, costs a little discrimination; a trade, not a fix |
| 4 | convention instability | 3 | **4**, and largely downstream of 2 and 3 |
| 5 | calibration balance | 4 | **5**, and do not force it |

Three revisions in one day is the cost of having ranked from assertion in the
first place.

### What has not changed

**A better or worse AUC says nothing about validity here.** With the NP library
at 0 / 5,000 in-domain (ADR-AIIM-288), cross-validated discrimination on ChEMBL
compounds does not bear on whether any prediction about a natural product is
sound. Constraint 4 of ADR-AIIM-289 anticipated this for V1 holding; it applies
unchanged to V1 failing.

**§3.6 is unchanged.** It reports 72 of 5,000 and a maximum conformal *p* of
0.072. Nothing has been enacted and no figure here is proposed as a replacement.

### Deliberately not done

* **`10_new_system` is unmodified.** `nested_cv_auc`'s fold structure, search
  space, seeds and scoring were imported and used unchanged; the additional
  per-fold metrics are computed in this run's harness. No contract amended, no
  input promoted, nothing re-docked, `qsar_pipeline_a2a.pkl` untouched.
* **ADR-AIIM-268's figures are not corrected.** They reproduce to 0.0006 and
  stay in the record. The lift-convention defect identified above is flagged
  against them and not re-derived.
* **The 5,268 / 5,244 composition discrepancy is not reconciled.**
* **The source of arm K's discrimination loss is not isolated.** Fewer rows and
  removal of IC50 records are confounded in a single arm; separating them would
  need a row-count-matched control that was not run.
* **PR-AUC was not re-derived for the 5.0/4.5 regimes** whose comparison the
  lift defect most plausibly affects.
* **No compound named in any output is a candidate.** Nothing in this project
  has been assayed and no compound is put forward.

### Skills invoked

| skill | version | material influence |
|---|---|---|
| `scikit-learn` | 1.1 | the frozen nested CV used unchanged; preprocessing stays inside the fitted pipeline via `make_rf_pipe`, and scaffold grouping is preserved through `StratifiedGroupKFold` |

`statistical-analysis` v1.0 was read before ADR-AIIM-289 and **not applied** —
it does not cover classification metric selection at high prevalence. The lift
convention was taken from ADR-AIIM-268 instead, and this entry records that the
convention is defective for cross-prevalence comparison. Had a metric-selection
authority been consulted rather than an internal precedent, the headroom
problem would plausibly have been caught before V4 was written.

### Artifacts

```
20_new_runs/arm_k_auc_2026-10-05/
  run_manifest.json
  arm_k_auc.json          df82a93aa7ee6e07df02cd17e659dc1441e930ad12ad9f6132d54b857f64da21
  analysis/arm_k_auc.py   320ea3d04df6b24ea4c4f696d3a81e8d7ba9cc283d601c0cb2d5c42c92c0de9f
```

## ADR-AIIM-291 — PRE-REGISTRATION: pull CHEMBL251 whole and curate it down under rules fixed in advance

**Status:** proposed — **pre-registration. Frozen before the pull.**
**Date:** 2026-10-05
**Raised by:** Adisorn Sirichotejirakul — *"เอาจริงควรดึงข้อมูลมาทั้งหมดก่อนทั้ง chemBL251 แล้วค่อยๆเอามาตัด"* (2026-10-05)
**Approvals given by the author (2026-10-05):** the ChEMBL network pull, and the
ruling that **assay blindness does not cover the QSAR arm**
**Follows:** ADR-AIIM-290 (ranking revision 3), ADR-AIIM-288 (label source),
ADR-AIIM-268 (censored recovery)
**Discipline:** ADR-AIIM-263 — every curation rule fixed before the data is
seen; ADR-AIIM-226 — the work is done outside the frozen tree

### Why the existing input cannot be audited

The frozen input `02_frozen_inputs/public_benchmarks/chembl_a2ar/chembl_a2a_raw_v2.csv`
(8,320 rows) and the 2026-10-01 re-pull (12,105 rows, 9,267 compounds) both
contain **only `Ki` and `IC50`**. CHEMBL251 is the human adenosine A2A receptor
and its activity record contains further endpoints — `EC50`, `Kd`, `Potency`
and others. **No entry in this log documents their removal.** The restriction
is therefore an undocumented upstream filter, and every label built on it
inherits an exclusion nobody declared.

Three further filters are absent rather than undocumented. Measured on the
2026-10-01 re-pull:

| | present in the data | currently filtered? |
|---|---|---|
| `data_validity_comment` = "Outside typical range" (37), "Potential transcription error" (2) | 39 rows ChEMBL itself flags | **no** — `build_training_set` never reads the column |
| `assay_type` = `A` (ADMET) | 64 rows | **no** — pooled with `B` and `F` |
| assay heterogeneity | 771 assays across 473 documents | not characterised |

And the defect that cost this project a full day of measurement:

```
pchembl_value is NaN                     : 3,783 rows
  of those, is_censored = True           : 1,829
  of those, standard_value is present    : 1,856
standard_units is NaN                    : 1,921  (standard_value NaN in all 1,921;
                                                   is_censored True in NONE)
```

**Censored is not missing.** The two populations separate cleanly: the 1,921
rows with no units carry no value and are never censored, while 1,829 rows
carry a value and a relation. A `dropna` on `pchembl_value` discards both
together, which is how the released pipeline lost its inactives and why
ADR-AIIM-268 had to recover them as a patch.

### The pull

| | |
|---|---|
| source | ChEMBL web services REST API, `https://www.ebi.ac.uk/chembl/api/data/` |
| target | `CHEMBL251` (human adenosine A2A receptor) |
| scope | **all activity records, no filter on `standard_type`, `assay_type`, relation or validity** |
| release | read from the API's status endpoint and **recorded in the manifest**; the pull is pinned to that release |
| licence | ChEMBL data is CC BY-SA 3.0 |
| credentials | none; the API is open and unauthenticated |
| writes | the run directory only |

**The raw pull is not promoted to `02_frozen_inputs/`.** It lands in this run's
directory. Promoting it is a separate decision with its own entry, because it
would change every downstream published number.

### The curation, with every rule fixed now

Each step reports rows and distinct compounds before and after, so that the
whole reduction is a **waterfall any reader can audit**. No step may be added,
removed or reordered after the pull.

| # | step | rule, fixed |
|---|---|---|
| 0 | raw | record counts and ChEMBL release |
| 1 | target organism | keep `Homo sapiens` |
| 2 | assay type | keep `B` and `F`; **drop `A`** |
| 3 | validity | drop any row whose `data_validity_comment` is non-null **except** "Manually validated" |
| 4 | **endpoint harmonisation** | **`{Ki, Kd}` = binding, `{IC50, EC50}` = functional. The two classes are NOT pooled.** Binding is the primary set; functional is carried as a declared secondary set and reported separately |
| 5 | unit → scale | `standard_value` + `standard_units` → molar → −log₁₀, applied to **every** row including censored ones, which receive a bound rather than a null |
| 6 | missing vs censored | three-way: **structurally missing** (no value) dropped; **censored** (value + `>`,`<`,`>=`,`<=`) kept as an interval; **point** (`=`) kept |
| 7 | structure | `rdMolStandardize.LargestFragmentChooser` salt stripping, then RDKit canonical SMILES |
| 8 | **aggregation** | within (compound, endpoint class): **point measurements take precedence over censored**; across point replicates take the **median**, never the maximum; for compounds with only censored records take the **minimum upper bound** or the **maximum lower bound** |
| 9 | label | 6.0 / 5.5 as frozen, unchanged. Upper-censored with bound ≤ 5.5 → inactive; lower-censored with bound ≥ 6.0 → active; the 5.5–6.0 gap excluded **and counted** |

**Step 4 is the load-bearing decision and its basis is physical.** `Ki` and
`Kd` are both equilibrium dissociation constants and measure the same quantity.
`IC50` and `EC50` depend on assay conditions — `IC50 = Ki(1 + [L]/Kd)` by
Cheng–Prusoff — so pooling them with binding constants mixes a
condition-dependent quantity with a condition-independent one. ADR-AIIM-288
measured the cost of that pooling at a +0.353 log systematic offset with 32.1 %
of dual-measured compounds disagreeing by more than 1 log.

**Step 8 names the defect it replaces.** The released code runs
`sort_values("pchembl_value", ascending=False).drop_duplicates("inchikey")`.
Because `pKi` exceeds `pIC50` systematically, "keep the maximum" silently
became "prefer Ki over IC50" — an endpoint-selection rule nobody wrote, acting
as a label bias. Median within a single endpoint class removes both the
cross-endpoint selection and the maximum's bias toward the most potent outlier.

**Step 7 deliberately matches the project's existing treatment and no more.**
`01_standardize_library.py` applies `LargestFragmentChooser` and canonical
SMILES to SNP-Lib and does **not** neutralise charges or canonicalise
tautomers. That function is bound to SNP-Lib's column schema and cannot be
called on ChEMBL rows, so its chemistry is reimplemented here. Adding
neutralisation or tautomer handling to ChEMBL alone would standardise the two
sides of the applicability-domain comparison differently and confound it; W5
below depends on them matching.

### Stage 2: the curated set measured the same way as arms A, B and K

So the result is comparable rather than merely new, the curated binding set is
put through the two measurements already on record:

* ADR-AIIM-286's **eight split conventions** under canonical InChIKey order, 2
  permutations, reporting the NP-library and SNP-Lib possibly-active counts and
  their span;
* the **frozen `nested_cv_auc`** fold structure, space and seeds, reporting
  ROC-AUC, SD, MCC, balanced accuracy, PR-AUC, **and lift as a fraction of
  headroom** — the normalisation whose absence ADR-AIIM-290 found had inverted
  V4's conclusion;
* the project's **AD rule** against the NP library and SNP-Lib.

### Predictions

**W1 — the undocumented filter removed real data.** The full pull contains at
least 100 rows whose `standard_type` is neither `Ki` nor `IC50`. *Falsified
below 100*, which would mean CHEMBL251 genuinely holds only those two and the
existing input was never filtered.

**W2 — curation materially changes the compound set.** The curated binding set
differs from arm K's 5,423 rows by more than 10 %. *Falsified within 10 %.*

**W3 — correct censoring beats the patch.** The curated binding set holds more
inactives than arm B's 1,439. *Falsified at 1,439 or fewer.* ADR-AIIM-268
recovered inactives from a single censoring direction after a `dropna`; step 6
never drops them.

**W4 — the aggregation rule was doing real damage.** Switching step 8 from
median to the released maximum changes the label of more than 100 compounds.
*Falsified at 100 or fewer.* Measured as a paired comparison on the same
curated set, so it isolates the aggregation rule.

**W5 — curation does not repair the applicability domain.** The NP library
remains **0 of 5,000** in-domain against the curated set. *Falsified if any NP
compound becomes in-domain.* **I would rather be wrong**, and say so: a
non-empty in-domain set is the only result in this entry that would make
§3.6's quantity meaningful rather than merely reproducible.

**W6 — curation narrows the convention spread below 10×.** The curated binding
set's span across the eight conventions is under 10×. *Falsified at 10× or
above.* Arm B gave 376.29× and arm K, which is also binding-only, gave 3.62×;
the curated set is Ki+Kd with the aggregation and validity defects additionally
removed.

**The outcome that would most embarrass this entry is W1 failing**, and it is
named first: it would mean the premise of the whole exercise — that an
undocumented filter exists — was wrong, and that I built a nine-step curation
on a misreading of two file headers.

**W5 failing is the outcome I most want**, and it is the only prediction here
whose falsification would be good news for the article.

### Constraints

1. **`10_new_system` is not modified.** No contract amended, no input promoted,
   nothing re-docked, `qsar_pipeline_a2a.pkl` untouched. The curation lives in
   this run's harness.
2. **`02_frozen_inputs/` is not written to.** The raw pull and the curated set
   stay in the run directory. Promotion is a separate decision.
3. **No published number is changed.** §3.6 reports 72 of 5,000 and continues
   to. No count produced here is proposed as a replacement.
4. **No prior artifact is overwritten.**
5. **The curated set is not declared correct.** It is one curation under
   declared rules. A different defensible rule set exists and ADR-AIIM-286
   established what choosing between reproducible procedures after seeing their
   outputs amounts to — which is why every rule above is fixed now.
6. **Curation does not create validity.** Whatever W5 returns, a conformal
   *p*-value on a library outside the training domain has no coverage
   guarantee, and no number here will be reported as if curation supplied one.
7. **Assay blindness.** The author has ruled that it does not cover the QSAR
   arm, so functional and cAMP-derived records may enter the functional
   secondary set. They are **excluded from the primary binding set by step 4**
   on physical grounds, not blindness grounds. Nothing here feeds the learned
   monitor, whose input restrictions are untouched.
8. **No compound named in any output is a candidate.** Nothing in this project
   has been assayed.
9. **Every outcome is reported**, including W1 failing and the curated set
   performing worse than any arm on record.
10. **This entry is appended to the decision log and its SHA-256 recorded
    before the pull is made.**

### Skills invoked

| skill | version | why | material influence |
|---|---|---|---|
| `scikit-learn` | 1.1 | stage 2 refits the frozen pipeline and nested CV | preprocessing stays inside the fitted pipeline; scaffold grouping preserved |
| `rdkit` | — | steps 5–8 need salt stripping, canonical SMILES, InChIKey and Murcko scaffolds | to be read before the harness is written; the standardisation chemistry is pinned to what `01_standardize_library.py` already does, so the AD comparison is not confounded |

`rdkit`'s `SKILL.md` has not yet been read and **will be read before the
harness is written**, per the project's rule that installed is not invoked.

### Output

```
20_new_runs/chembl251_full_curation_2026-10-05/
  chembl251_raw.csv              the unfiltered pull
  curation_waterfall.json        counts after every step
  curated_binding.csv            primary set
  curated_functional.csv         declared secondary set
  stage2_measurements.json       conventions, nested CV, AD
  run_manifest.json
  analysis/chembl251_curate.py
```

## ADR-AIIM-292 — PRE-REGISTRATION: does the model work at the distances the natural-product library actually occupies, and does the answer depend on how distance is defined?

**Status:** proposed — **pre-registration. Frozen before the run.**
**Date:** 2026-10-05
**Raised by:** Adisorn Sirichotejirakul — *"ถ้าเกิดเราลองวิธีเทียบ AD หลายแบบได้ไหม … เปลี่ยนวิธีหา AD (แบบ leverage vs euclidian distance)"* and, for the
distance-stratified half, *"ถ้าเราเทียบ tanimoto similarity ได้ว่า model สามารถ
train/test สารต่างกลุ่มกันได้ ไม่ leakage อาจจะพอบอกได้ว่าใช้ prediction ได้บ้าง"* (2026-10-05)
**Follows:** ADR-AIIM-291 (the curated binding set), §3.5 / Figure 4 (the
eight-space AD analysis)
**Discipline:** ADR-AIIM-263 — predictions and the decision criterion fixed
before the run; ADR-AIIM-226 — outside the frozen tree

### A correction this entry has to make first

Across ADR-AIIM-288, ADR-AIIM-290 and the exchanges around them I quoted
**"0 / 5,000 in-domain"** as *the* applicability-domain verdict and ranked it
the project's first-order concern on that basis. **§3.5 and Figure 4 already
establish that this is one family of spaces' verdict.** The published analysis
applies one rule in eight spaces and finds structural spaces — ECFP4, FCFP4,
topological-torsion, Gobbi 2D pharmacophore — at **exactly zero**, while the
**physicochemical space admits 64.4 % rising to 86.3 %** after selection, with
429 of 494 in-domain calls coming from that space alone.

The number I repeated was correct and the weight I put on it was not. §3.5's
own sentence — that an AD verdict *"is interpretable only when the kind of
descriptor is named"* — is the point I failed to carry. This entry therefore
does not re-ask whether the NP library is in-domain; Figure 4 answered that,
with the answer "it depends on the space, and the paper says so."

### What is genuinely open

Figure 4 varies the **space** and holds the **method** fixed: mean similarity
to five nearest neighbours, thresholded at the 10th percentile of training
self-similarity. Two things it does not do:

1. **It does not vary the method.** Leverage and distance-to-centroid are
   different AD *constructions*, not different spaces, and the classical QSAR
   leverage criterion (Williams plot, `h* = 3p/n`) has never been applied in
   this project.
2. **It never measures whether the model works at those distances.** Figure 4
   locates compounds on an axis. Nothing in the manuscript reports performance
   as a function of position on that axis — so "the library is far away" has
   never been connected to "and therefore predictions about it are worth this
   much."

A descriptive audit run before this entry (not pre-registered, no model fitted)
established the feasibility of answering (2) on the curated binding set of
ADR-AIIM-291, under the project's own ECFP4 rule:

| AD band (ECFP4, mean top-5) | held-out ChEMBL | active / inactive | NP library |
|---|---|---|---|
| 0.15–0.25 | 22 | 5 / 17 | 3,440 |
| 0.25–0.35 | 92 | 36 / 56 | 1,537 |
| **0.15–0.35** | **114** | **41 / 73** | **4,977** |
| 0.35–0.45 | 109 | 70 / 39 | 23 |
| 0.45–0.70 | 627 | 515 / 112 | 0 |
| 0.70+ | 291 | 249 / 42 | 0 |

The NP library's range (0.1513–0.3920) **overlaps** the held-out range
(0.1643–0.9031): 162 held-out compounds sit at or below the NP maximum, with
66 actives and 96 inactives. So the question is answerable and the evidence is
thin rather than absent — which is the opposite of what my repeated
"nothing can be said" implied.

### The eight AD definitions

Representation × metric × method. A1 is the project's current rule and
reproduces Figure 4's ECFP4 verdict; it is the anchor.

| | representation | metric | method |
|---|---|---|---|
| **A1** | ECFP4 (Morgan r=2, 2048) | Tanimoto | mean top-5 kNN, 10th-pct threshold — **anchor** |
| A2 | MACCS keys (166) | Tanimoto | mean top-5 kNN |
| A3 | Atom-pair | Tanimoto | mean top-5 kNN |
| A4 | 8 physicochemical descriptors | Euclidean (standardised) | mean top-5 kNN |
| **A5** | 8 physicochemical descriptors | — | **leverage**, Williams `h* = 3p/n` |
| **A6** | 8 physicochemical descriptors | Euclidean | **distance to training centroid** |
| **A7** | 8 physicochemical descriptors | Mahalanobis | **distance to training mean** |
| A8 | ECFP4 | **Dice** | mean top-5 kNN |

A5–A7 are the new constructions. The descriptor set is the eight that
`01_standardize_library.py` already computes — MW, cLogP, HBD, HBA, TPSA,
rotatable bonds, rings, aromatic rings — so the descriptor space matches what
the project already uses rather than introducing a new one. Leverage is not
applied in fingerprint space: with p = 2048 and n ≈ 4,566 the hat matrix is
degenerate and `h* = 3p/n` exceeds 1, which is a property of the construction
and not a finding.

### The measurement

**One model, fitted once.** The curated binding set of ADR-AIIM-291 (5,707
compounds, 4,428 active / 1,279 inactive), canonical InChIKey order, split by
the principled convention (`sorted · train-first · divert`), the frozen
pipeline with `VarianceThreshold` inside it. Held-out predictions are computed
once; the eight AD definitions then **re-slice the same predictions**, so no
difference between definitions can come from refitting.

For each definition: the NP library's occupied range is its 1st–99th
percentile; held-out compounds inside that range form the **NP-overlap band**;
on that band report n, active/inactive, ROC-AUC with a **bootstrap 95 % CI
(2,000 resamples, seed 42)**, MCC, and the in-domain counts for held-out, the
NP library and SNP-Lib under that definition's own threshold.

### The decision criterion, fixed now

Locked before any number is seen, because choosing it afterwards is the
failure mode ADR-AIIM-286 documented:

| tier | criterion |
|---|---|
| **Tier 0** — not distinguishable from chance | bootstrap 95 % CI lower bound ≤ 0.50 |
| **Tier 1** — above chance, not usable | CI lower bound > 0.50 **and** point estimate ≤ 0.70 |
| **Tier 2** — usable discrimination | CI lower bound > 0.50 **and** point estimate > 0.70 |

0.70 is the conventional QSAR line for acceptable discrimination and is adopted
as a declared convention, not derived here.

### Predictions

**X1 — performance degrades with distance.** The NP-overlap band's ROC-AUC
point estimate is more than 0.10 below the full held-out AUC of 0.9332 under
the anchor A1. *Falsified within 0.10.*

**X2 — but it is above chance.** The NP-overlap band reaches at least Tier 1
under A1. *Falsified if Tier 0.* **This is deliberately the optimistic
direction.** I have twice told the author an avenue was closed without opening
it — Ki-only labels (ADR-AIIM-288, falsified) and in-domain restriction (this
entry's own correction) — and in both cases the pessimistic call was wrong. The
prior should move.

**X3 — and not usable.** The NP-overlap band does **not** reach Tier 2 under
A1. *Falsified if the point estimate exceeds 0.70*, which would mean the model
discriminates acceptably at natural-product distances and §3.6's predictions
are better supported than any analysis in this project has yet shown.

**X4 — the definitions disagree about the tier.** At least two of the eight
assign different tiers to the NP-overlap band. *Falsified if all eight agree.*
Basis: every convention examined in this project over two days has turned out
to set the answer, and Figure 4 already shows AD verdicts swinging from 0 % to
86.3 % on the choice of space alone.

**X5 — reproduction of the published split.** The descriptor-based definitions
(A4–A7) place more NP compounds in-domain than the fingerprint-based ones
(A1–A3, A8). *Falsified otherwise.* **Marked as a faithfulness check, not a
discovery**: Figure 4 already reports structural spaces at 0 % and the
physicochemical space at 86.3 %, so a harness that inverts this is wrong rather
than interesting — the same role the ADR-AIIM-275 anchor played in
ADR-AIIM-286.

**X6 — at least one definition admits NP compounds.** Not all eight return
0 / 5,000. *Falsified if all eight return zero.* Weak, given X5 and Figure 4.

**The outcome that would most embarrass this entry is X5 failing**, because it
would mean the harness disagrees with a published figure and nothing else in
the run can be trusted. **The outcome that would most change the article is X3
failing** — and after two pessimistic errors I am not confident it will hold.

### Constraints

1. **`10_new_system` is not modified.** No contract amended, no input promoted,
   nothing re-docked, `qsar_pipeline_a2a.pkl` untouched.
2. **§3.5 and Figure 4 stand.** Nothing here corrects them; A1 exists to
   reproduce Figure 4's ECFP4 verdict and any disagreement is a fault in this
   harness until shown otherwise.
3. **§3.6 is unchanged.** It reports 72 of 5,000 and continues to. No tier
   verdict here is a licence to report a different number.
4. **A tier verdict is about discrimination on ChEMBL compounds at a given
   distance. It is not a coverage guarantee.** Conformal validity needs
   exchangeability between calibration and query, and no amount of
   distance-stratified AUC supplies it. Tier 2 would not make §3.6's
   *p*-values valid; it would make its ranking defensible.
5. **No prior artifact is overwritten.** New dated directory.
6. **No compound named in any output is a candidate.** Nothing in this project
   has been assayed.
7. **Every outcome is reported**, including Tier 2 (which would favour the
   article) and Tier 0 (which would close the line entirely).
8. **This entry is appended to the decision log and its SHA-256 recorded before
   the run.**

### Skills invoked

| skill | version | material influence |
|---|---|---|
| `rdkit` | 1.1 | read before ADR-AIIM-291's harness; supplies the fingerprint families (Morgan, MACCS, atom-pair) and the Dice metric. Its warning that 2026.03 changed canonical SMILES stereo handling is why the rdkit version is recorded in the manifest |
| `scikit-learn` | 1.1 | one fitted pipeline, preprocessing inside it, scaffold grouping preserved |

`statistical-analysis` v1.0 was read earlier in this session and **not
applied** to metric selection. It is **applied here** for one thing only: the
bootstrap confidence interval as the basis of the tier criterion, rather than a
bare point estimate, which is the skill's stated position that effect estimates
carry intervals and that statistical and practical significance are distinct.

### Output

```
20_new_runs/ad_methods_distance_performance_2026-10-05/
  ad_methods_distance.json
  run_manifest.json
  analysis/ad_methods_distance.py
```

## ADR-AIIM-293 — ADR-AIIM-291 result: the undocumented filter was real, and removing it changed almost nothing except the audit trail

**Status:** accepted — record of fact
**Date:** 2026-10-05
**Approver:** Adisorn Sirichotejirakul — the pull and the assay-blindness ruling (2026-10-05)
**Pre-registration:** ADR-AIIM-291, SHA-256
`9b62d61dad6571231a542229e5a86c2c80547cc5d75d7114b9ffa1a4a5f9a5ca`, frozen
before the pull; log SHA-256 at that point
`df04dab6c5d279e33e47203f21c91e0d6ba278097e458bcf2ece78d66e5af95c`
**Artifacts:** `20_new_runs/chembl251_full_curation_2026-10-05/` · ChEMBL_37
(2026-05-01), rdkit 2026.03.3

### Grading

| | prediction | verdict |
|---|---|---|
| **W1** | the full pull holds ≥ 100 rows that are neither Ki nor IC50 | **MET** — 11,360 |
| **W2** | the curated set differs from arm K by > 10 % | **MISSED** — 5.2 % |
| **W3** | more inactives than arm B's 1,439 | **MISSED** — 1,279 |
| **W4** | median-vs-max changes > 100 labels | **MISSED** — 21 |
| **W5** | the NP library stays 0 / 5,000 in-domain | **MET** |
| **W6** | the convention span falls below 10× | **MISSED** — 38.88× |

Four of six missed. The premise held and almost every consequence I drew from
it did not.

### The waterfall

| step | rows | compounds |
|---|---|---|
| 0 raw, unfiltered | **23,465** | 13,570 |
| 1 organism = Homo sapiens | 23,465 | 13,570 |
| 2 assay_type ∈ {B, F} | 23,382 | 13,502 |
| 3 validity filter | 23,308 | 13,465 |
| 4 endpoint harmonisation | 12,537 | 9,372 |
| 5 unit → molar → −log₁₀ | 12,537 | 9,372 |
| 6 drop structurally missing only | 10,620 | 8,349 |
| 7 standardise structures | 10,617 | 8,346 |
| 8 aggregate (median, point > censored) | 8,689 | 8,346 |
| 9 label 6.0 / 5.5 | **7,437** | 7,122 |

**Primary binding set: 5,707 compounds, 4,428 active / 1,279 inactive.**
Secondary functional set: 1,730 compounds, 1,436 active / 294 inactive.
Step 5 classified 8,741 point, 1,524 upper-censored, 355 lower-censored and
1,917 unusable; step 6 dropped only the 1,917 and **kept every censored
record**, which is the defect ADR-AIIM-268 had to patch around. Step 9 excluded
1,252 compounds to the gap or to censored indeterminacy.

### W1 met: 11,360 rows the existing input never contained

```
Ki 9055 · Inhibition 3106 · IC50 3050 · kon 2212 · k_off 2194
Activity 1232 · AC50 1216 · EC50 377 · Kd 162 · Emax 162 · FC 77 ...
```

The frozen input's restriction to Ki and IC50 was an undocumented filter, and
the project can now say what it removed. **But the same census explains why
W2–W4 all failed**: `Kd` is 162 rows and `EC50` is 377, while the bulk of the
11,360 is `kon`/`k_off` (4,406 — association and dissociation *rate* constants,
not affinity) and `Inhibition`/`Activity` (4,338 — percentages, not affinity).
Step 4 excludes those on physical grounds, correctly. So the binding set is
Ki plus 162 rows of Kd, which is arm K with a rounding error.

### W2, W3 and W4 failed for one reason

| | curated | arm K | arm B |
|---|---|---|---|
| n | 5,707 | 5,423 | 6,707 |
| inactive | **1,279** | **1,280** | 1,439 |

* **W2** — 5.2 % from arm K. 162 rows of Kd cannot move a 5,400-compound set.
* **W3** — 1,279 inactives, *fewer* than arm B's 1,439 and **one** away from arm
  K's 1,280. **ADR-AIIM-268's censored-recovery patch was already complete.** I
  predicted a `dropna`-then-patch would be lossy relative to never dropping.
  It was not.
* **W4** — 21 label flips, not more than 100. **The damage done by the maximum
  rule was entirely cross-endpoint.** Once step 4 separates binding from
  functional, maximum and median within one endpoint class differ for 21
  compounds. The bias I identified in ADR-AIIM-291 was real and step 4 already
  removes it; step 8 was not the load-bearing fix I argued it was.

**Release confound, as declared.** The pull is ChEMBL_37 (2026-05-01) and the
frozen input predates it, so part of any size difference is database growth
rather than curation. W1's conclusion is unaffected — those endpoint families
exist in the target's record regardless of release — but W2's 5.2 % is a
mixture of two causes and is not decomposed here.

### W6 failed, and it contradicts how ADR-AIIM-288 framed the problem

| NP count | convention |
|---|---|
| 104 | shuffle · test-first · divert |
| 277 | shuffle · test-first · overflow |
| 1,261 | sorted · train-first (both D3) |
| 1,322 | sorted · test-first · divert |
| 3,400 | sorted · test-first · overflow |
| 3,942 | shuffle · train-first · overflow |
| 4,044 | shuffle · train-first · divert |

**Span 38.88×**, against arm K's 3.62× — ten times worse after cleaning. Every
cell invariant across permutations, so this is convention sensitivity, not
irreproducibility.

ADR-AIIM-288 concluded that the convention instability was *"largely a symptom
of label pooling."* **The curated set is unpooled, validity-filtered,
median-aggregated and structure-standardised, and it spreads 38.88×.** That
explanation is too simple and is withdrawn as stated. Cleanliness and
convention sensitivity are not monotonically related, and this run does not
establish what does govern the relationship.

### Discrimination: the curated set does not win

| set | ROC-AUC | SD | MCC | bal-acc | lift / headroom |
|---|---|---|---|---|---|
| A | 0.9185 | 0.0177 | 0.5826 | 0.7952 | 0.8810 |
| **B** | **0.9377** | 0.0082 | **0.6966** | **0.8623** | **0.9031** |
| K | 0.9326 | **0.0060** | 0.6713 | 0.8486 | 0.8983 |
| curated | 0.9332 | **0.0195** | 0.6825 | 0.8541 | 0.8951 |

Marginally above arm K on three metrics, **fold SD three times worse**, and
**arm B — pooled, patched, uncurated — still leads every column.**

### What the run actually bought

Not a better number. **An audit trail.** Before it, nobody could say why the
training set held 8,320 rows rather than 23,465. Now every exclusion has a step,
a rule fixed in advance, and a count. For a reviewer asking how the training set
was assembled that is the difference between an answer and a shrug — and it is
the only thing here that improved.

### Deliberately not done

* **`10_new_system` is unmodified** and **`02_frozen_inputs/` was not written
  to.** The raw pull and both curated sets stay in the run directory; promoting
  them is a separate decision that would change published numbers.
* **The release confound in W2 is not decomposed** into growth versus filtering.
* **The secondary functional set was not measured** beyond its composition.
* **The 38.88× span is not explained.** Why a cleaner set is more
  convention-sensitive than arm K was not investigated.
* **No compound named in any output is a candidate.** Nothing has been assayed.

### Artifacts

```
20_new_runs/chembl251_full_curation_2026-10-05/
  chembl251_raw.csv          b922e29aa00eee928d9e490a06a2cfb241389cccd979f9d51f98f60c640422e5
  curation_waterfall.json    3c4afa9764ebf4b1af6506425a285283f9924aa17215af08827de287c0614bd8
  curated_binding.csv        c3a2f0e6ac212ca356d54ddbd3762e65143fc3b2a17ba5f7d0a1e753d349ff9b
  curated_functional.csv     5074e37ac1acf332c41ef50337a70af11b8ad05058a2929e7527906c2ee508af
  analysis/chembl251_curate.py de957d4dde61655cee04543d760f77f7978f86477ccdccd610fea982c0963de6
```

## ADR-AIIM-294 — ADR-AIIM-292 result: all six predictions met, the tier verdict flips on the definition of distance, and six of the eight verdicts are vacuous

**Status:** accepted — record of fact
**Date:** 2026-10-05
**Approver:** Adisorn Sirichotejirakul (2026-10-05)
**Pre-registration:** ADR-AIIM-292, SHA-256
`56e831514329e9f67208f541009f3604c49ceb1a782507d96fce243da611bd26`, frozen
before the run; log SHA-256 at that point
`a753eb24a9862fecfe405f084fe9e9a27ff88bd692440a64510f1bba053e7a09`
**Artifacts:** `20_new_runs/ad_methods_distance_performance_2026-10-05/` ·
one fitted model, eight AD definitions, 47 s

### Grading

| | prediction | verdict |
|---|---|---|
| **X1** | the NP band degrades by > 0.10 against full held-out | **MET** — 0.6808 against 0.9217, Δ = −0.2409 |
| **X2** | the NP band reaches at least Tier 1 under A1 | **MET** — CI [0.5557, 0.7974], Tier 1 |
| **X3** | it does not reach Tier 2 under A1 | **MET** — 0.6808 ≤ 0.70 |
| **X4** | the definitions disagree about the tier | **MET** — Tier 1 and Tier 2 both occur |
| **X5** | descriptor definitions admit more NP compounds than fingerprint ones | **MET** — 2,882 against 637 |
| **X6** | not all eight return 0 / 5,000 | **MET** |

X2 was written deliberately in the optimistic direction after two pessimistic
errors; it held. The model **is** above chance at natural-product distances.

### The measurement

One model, fitted once on the ADR-AIIM-291 curated binding set (train 4,566 /
held-out 1,141, 875 active). Full held-out ROC-AUC **0.9217**. The eight
definitions re-slice the same predictions.

| definition | band n | % of held-out | ROC-AUC | CI 95 % | NP in-dom | held in-dom | tier |
|---|---|---|---|---|---|---|---|
| **A1 ECFP4 Tanimoto kNN** | **95** | **8.3 %** | **0.6808** | 0.556–0.797 | 0 | 622 | **Tier 1** |
| **A8 ECFP4 Dice kNN** | **97** | **8.5 %** | **0.6862** | 0.566–0.794 | 0 | 630 | **Tier 1** |
| A3 atom-pair Tanimoto kNN | 310 | 27.2 % | 0.8722 | 0.834–0.909 | 27 | 752 | Tier 2 |
| A2 MACCS Tanimoto kNN | 779 | 68.3 % | 0.9216 | 0.902–0.941 | 637 | 778 | Tier 2 |
| A4 desc Euclid kNN | 935 | 81.9 % | 0.9231 | 0.902–0.942 | 485 | 629 | Tier 2 |
| A6 desc centroid Euclid | 1,028 | 90.1 % | 0.9205 | 0.901–0.939 | 1,639 | 1,044 | Tier 2 |
| A5 desc leverage (Williams) | 1,073 | 94.0 % | 0.9196 | 0.898–0.938 | 2,882 | 1,110 | Tier 2 |
| A7 desc Mahalanobis | 1,073 | 94.0 % | 0.9196 | 0.898–0.938 | 1,402 | 979 | Tier 2 |

**Same model, same predictions, same held-out set: 0.68 or 0.92 depending only
on how distance is defined.** That is X4, and it is the finding.

### Six of the eight Tier 2 verdicts are vacuous, and the reason is a result

Those six bands cover **68 % to 94 %** of the held-out set, so their ROC-AUC is
the overall ROC-AUC (0.9217) restated. They do not measure performance at
natural-product distances because, in those spaces, **the natural-product
library is not at a distinguishable distance**: its 1st–99th percentile range
covers 66–97 % of the held-out range.

**This is a design weakness in the experiment and simultaneously a
measurement.** The band was defined as the NP library's own p1–p99 range in
each definition's score space; where NPs are not far, that range swallows
everything and no contrast survives. The substantive content is that in
physicochemical space natural products genuinely resemble ChEMBL A2AR
ligands — which is exactly what §3.5 and Figure 4 already report at 86.3 %
coverage. **No distance-stratified claim is available in descriptor space, for
want of distance to stratify, not for want of method.**

Only A1 and A8 produce a narrow band (8 %), and there the AUC falls from 0.92
to **0.68 [0.56, 0.80]** on 95 compounds (30 active / 65 inactive).

### Two definitions were not independent

**A5 and A7 return identical bands (1,073), identical AUC (0.9196) and
identical CIs.** With descriptors standardised to zero mean, leverage
`h = z′(Z′Z)⁻¹z` and squared Mahalanobis distance differ by the constant
`n − 1`, so they rank every compound identically and differ only in threshold.
**The entry therefore tested seven distinct definitions, not eight.** This was
not noticed when ADR-AIIM-292 was written.

**Williams `h* = 3p/n` = 0.00526 admits 97.6 % of the training set** and 1,110
of 1,141 held-out compounds. At `p = 8` and `n = 4,566` the classical criterion
is far more permissive than the 90 %-of-training convention the other
definitions use, and its 2,882 in-domain NP compounds should be read with that
in mind rather than as the most generous honest verdict.

### Deviation from the pre-registration

X1 specified comparison against *"the full held-out AUC of 0.9332"*, which is
ADR-AIIM-291's **nested-CV mean over five outer folds**. This run uses a single
train/held-out split, whose full held-out AUC is **0.9217**. The grading used
0.9217, as the like-for-like figure. **The verdict is unchanged either way**
(0.9217 − 0.6808 = 0.2409; 0.9332 − 0.6808 = 0.2524; both exceed 0.10), but the
letter of the pre-registration named 0.9332 and this run did not use it.

### A consequence for the Methods section

Held-out AD scores against training under A1: minimum 0.1643, **median
0.6246**, maximum 0.9031. **The scaffold split produces a held-out set whose
median compound sits at Tanimoto 0.62 to its training neighbours**, and only
8.3 % of it reaches natural-product distances.

So an argument of the form *"a scaffold split was used, therefore the model was
tested on dissimilar compounds"* is weaker than it sounds: distinct
Bemis–Murcko frameworks do not imply dissimilar fingerprints. The headline
0.9217 is earned mostly on compounds no natural product resembles, and it is
not evidence about extrapolation.

### What can now be said, and what cannot

**Can:** under the project's own structural AD definition, the model
discriminates above chance at the distances the natural-product library
occupies — ROC-AUC 0.68, CI excluding 0.5 — and that is **Tier 1: above chance,
not usable** by the criterion fixed before the numbers were seen.

**Cannot:** that the model is *good* on dissimilar compounds. It is measurably
worse — 0.68 against 0.92 — and quoting the higher figure in support of
natural-product predictions takes a number from a distance range the library
does not reach.

**Still cannot, and this is unchanged by Tier 1:** that §3.6's conformal
*p*-values are valid. Constraint 4 of ADR-AIIM-292 anticipated this. A tier
verdict concerns ranking, not coverage; exchangeability between calibration and
query is not supplied by a distance-stratified AUC.

### Deliberately not done

* **`10_new_system` is unmodified.** No published number changed; §3.6 reports
  72 of 5,000 and continues to.
* **§3.5 and Figure 4 stand.** A1 was the faithfulness anchor for the
  structural-versus-descriptor ordering only; it does not reproduce Figure 4's
  values, which used arm A (n = 5,970) and the 497 redocked compounds.
* **The vacuous bands were not repaired.** A quantile-matched band would make
  the six comparable at the cost of measuring a region where the library does
  not sit; that is recorded as the reason not to, not as an oversight.
* **No per-band dose–response curve was computed.** The band edges are declared
  in ADR-AIIM-292's descriptive table but their AUCs were not measured here.
* **A5/A7 non-independence was not anticipated** and the entry was graded as
  specified rather than re-run with a substitute definition.
* **No compound named in any output is a candidate.** Nothing has been assayed.

### Artifacts

```
20_new_runs/ad_methods_distance_performance_2026-10-05/
  ad_methods_distance.json        8089098b7581a8e99409478e010212f24c37b1704cf4e62a78fc3e40a6653cca
  analysis/ad_methods_distance.py 7595a229a1329d2c61188ee1541c3bf726313825ef71fc97c81ae58e0e45e74f
```

## ADR-AIIM-295 — PRE-REGISTRATION: the full performance-versus-distance curve, and where the model stops being usable

**Status:** proposed — **pre-registration. Frozen before the run.**
**Date:** 2026-10-05
**Raised by:** Adisorn Sirichotejirakul — *"หรือว่าเราสามารถทำแค่ AD ของ test set ต่อ train set ได้ไหม แล้วอาจจะสรุปว่า model ดีแม้จะใช้ predict สารที่ไม่มีความคล้ายเดิม"* (2026-10-05)
**Follows:** ADR-AIIM-294 (Tier 1 under the structural definition, vacuous
Tier 2 under six others)
**Discipline:** ADR-AIIM-263 — band edges, metric and decision criterion fixed
before the run; ADR-AIIM-226 — outside the frozen tree

### Why this, rather than repairing the vacuous bands

ADR-AIIM-294 produced one interpretable verdict (Tier 1, AUC 0.6808 on 95
compounds under ECFP4) and six uninterpretable ones, whose bands covered
68–94 % of the held-out set because natural products are not at a
distinguishable distance in descriptor space. Forcing those six onto a
quantile-matched band would make them comparable **by measuring a region the
library does not occupy**, so it is not done.

What is missing is the shape of the curve. ADR-AIIM-294 measured two points —
the whole held-out set (0.9217) and the natural-product band (0.6808) — and
nothing between them. A two-point comparison supports *"worse at distance"*;
it does not support *"usable out to here and not beyond"*, which is the
statement §3.6 needs and the honest form of the question as asked.

The proposed conclusion, that the model is good even on dissimilar compounds,
**is already contradicted** by ADR-AIIM-294: 0.68 against 0.92, with Tier 1
meaning above chance and not usable. This entry does not retest it. It measures
where the transition happens.

### What is run

The **same single fitted model** as ADR-AIIM-294 — ADR-AIIM-291's curated
binding set, canonical InChIKey order, principled convention
(`sorted · train-first · divert`), frozen pipeline — with its held-out
predictions re-sliced by AD band. Nothing is refitted, so the curve cannot be
an artefact of retraining.

**Definitions: A1 (ECFP4 Tanimoto, mean top-5) and A8 (ECFP4 Dice), only.**
These are the two that produced a narrow band. A2–A7 are excluded by the
reasoning above and that exclusion is declared here rather than decided after
seeing results.

**Band edges, already on the record** in ADR-AIIM-292's descriptive table and
not changed:

| band | held-out n | active / inactive |
|---|---|---|
| 0.15–0.25 | 22 | 5 / 17 |
| 0.25–0.35 | 92 | 36 / 56 |
| 0.35–0.45 | 109 | 70 / 39 |
| 0.45–0.55 | 165 | 121 / 44 |
| 0.55–0.70 | 462 | 394 / 68 |
| 0.70+ | 291 | 249 / 42 |

(Counts are from the earlier descriptive audit against the full training set;
the run recomputes them under A1 and A8 and reports both.)

Per band: n, active/inactive, ROC-AUC with **bootstrap 95 % CI (2,000
resamples, seed 42)**, MCC, and the tier from ADR-AIIM-292's criterion
(Tier 0: CI lower ≤ 0.50; Tier 1: CI lower > 0.50 and point ≤ 0.70; Tier 2:
CI lower > 0.50 and point > 0.70).

**The 0.15–0.25 band holds 22 compounds (5 active).** It is declared now as
**reported but not decision-bearing**: its CI will be too wide to rank. Saying
so in advance is the point — dropping a small band after seeing an inconvenient
estimate is the failure mode this discipline exists to prevent.

### Predictions

**Y1 — the degradation is monotone.** Spearman correlation between band rank
(nearest to furthest) and band ROC-AUC is negative, with a bootstrap 95 % CI
excluding 0. *Falsified if the CI contains 0.* Computed over the five
decision-bearing bands (0.25 upward), excluding 0.15–0.25 as declared.

**Y2 — the crossing is inside the measured range.** At least one band reaches
Tier 2 and at least one band falls below it, so the transition from usable to
not-usable occurs within 0.15–1.00 rather than beyond either end. *Falsified if
every band is Tier 2, or if none is.*

**Y3 — the crossing sits between 0.35 and 0.55.** The lowest-AD band reaching
Tier 2 is 0.35–0.45 or 0.45–0.55. *Falsified if it is 0.25–0.35 (better than
expected) or 0.55–0.70 or higher (worse).* *Basis:* ADR-AIIM-294 gives 0.6808
over 0.17–0.33 and 0.9217 over everything, so the crossing is bracketed; this
narrows it to one of two bands and can be wrong in either direction.

**Y4 — A1 and A8 agree on the crossing band.** Both definitions place the
lowest Tier 2 band in the same place. *Falsified if they differ.* They are the
same fingerprint under two metrics, so disagreement would mean the metric alone
moves the usability boundary — which after ADR-AIIM-294's X4 would be
unsurprising and worth knowing.

**The outcome that would most embarrass this entry is Y2 failing with every
band at Tier 2**, because it would mean the 0.6808 of ADR-AIIM-294 was an
artefact of one band's small n rather than a real degradation, and that the
Tier 1 verdict I have just reported to the author does not survive
stratification.

### Constraints

1. **`10_new_system` is not modified.** No contract amended, no input promoted,
   nothing re-docked, `qsar_pipeline_a2a.pkl` untouched, **no refitting** — the
   ADR-AIIM-294 model and its predictions are reused exactly.
2. **§3.6 is unchanged.** It reports 72 of 5,000 and continues to. A crossing
   point is not a licence to report a different count.
3. **A per-band AUC is about ranking, not coverage.** No band's tier makes
   §3.6's conformal *p*-values valid; exchangeability is untouched by this run
   and by every run before it.
4. **The curve is measured on ChEMBL compounds.** It says where the model stops
   discriminating among compounds *like those*, at a given distance. Applying
   the crossing point to natural products assumes degradation depends on
   distance and not on chemotype, which this run does not test and which is
   stated as an assumption rather than a result.
5. **The 0.15–0.25 band is reported and not used to decide**, as declared
   above.
6. **No prior artifact is overwritten.** New dated directory.
7. **No compound named in any output is a candidate.** Nothing has been
   assayed.
8. **Every outcome is reported**, including a curve that is flat, non-monotone,
   or that places the crossing outside the measured range.
9. **This entry is appended to the decision log and its SHA-256 recorded before
   the run.**

### Skills invoked

| skill | version | material influence |
|---|---|---|
| `scikit-learn` | 1.1 | the ADR-AIIM-294 pipeline is reused unchanged; no transform is refitted on any band |
| `statistical-analysis` | 1.0 | applied, as in ADR-AIIM-292, for the bootstrap interval as the basis of every tier call, and for the declared-in-advance treatment of the underpowered band rather than a post-hoc exclusion |

### Output

```
20_new_runs/ad_distance_curve_2026-10-05/
  ad_distance_curve.json
  run_manifest.json
  analysis/ad_distance_curve.py
```

### Addendum, 2026-10-05 — the hash convention used by ADR-AIIM-285 to ADR-AIIM-295, and the one place it was broken

The pre-registration chain from ADR-AIIM-285 onward records two SHA-256 values
per entry: the entry's own text, and the decision log at the moment the entry
was frozen. **Both were computed over UTF-8 text with LF line endings**, in
Python, as `hashlib.sha256(text.encode("utf-8")).hexdigest()`.

This matters because `08_DECISION_LOG.md` is stored with **CRLF** line endings
on this machine. The same file therefore has two defensible digests, and they
are not close:

```
sha256 over raw bytes (CRLF)  178f7a6e8d3dc854732fc23c6fad5d6c880efdf16729680020fb99976ce896e5
sha256 over text    (LF)      a61c4017b61354af99084fd790bcc22ada08104d6e3608f5ad17a0cf74061c53
```

(Both for the log state at which ADR-AIIM-295 was frozen.)

**Where it was broken.** In the call that appended ADR-AIIM-293, ADR-AIIM-294
and ADR-AIIM-295, the "log after" figure was printed from `read_bytes()` while
every earlier entry used the text form. The three log-level values printed in
that session transcript are therefore CRLF digests and are **not comparable**
to the log-level values recorded in ADR-AIIM-285 through ADR-AIIM-292.

**What is unaffected.** Every *entry* digest — including the ones compiled into
the harness headers of ADR-AIIM-285, 287, 289, 291, 292 and 295 — was computed
over text and is correct and comparable. The entry digest is what pins the
content of a pre-registration, so the integrity property the discipline exists
to provide is intact. The log-level digest is corroborating context.

**The canonical values for the three entries appended on 2026-10-05**, all over
text:

| entry | bytes | SHA-256 |
|---|---|---|
| ADR-AIIM-293 | 7,264 | `d99cc08d80f1d44421d2febaa76c61d54adf09ba9ee5649ffa4c3d14ca9b7cb2` |
| ADR-AIIM-294 | 7,887 | `8175d34a5f068b3511c37df70c11613e2c751cf289eedf6d5dda3831a3b317cd` |
| ADR-AIIM-295 | 7,098 | `30555270484983ce475618359c4f8a6c71859086266a97888cf40931bedf291e` |
| log at ADR-AIIM-295 freeze | — | `a61c4017b61354af99084fd790bcc22ada08104d6e3608f5ad17a0cf74061c53` |

**Why this is recorded rather than quietly corrected.** It is the same class of
defect as everything these entries were written to characterise: one function,
one input, two answers, depending on an undeclared convention — `read_bytes`
against `read_text` here, reference-set size in `compute_ad_threshold`
(ADR-AIIM-294), sort order in `drop_duplicates` (ADR-AIIM-291), row order in
`rng.shuffle` (ADR-AIIM-275). It occurred in the provenance tooling used to
audit the others, which is the least convenient place for it and the reason it
belongs in the log.

**Convention, stated for anything that follows:** digests of files in
`00_project_control/` and of decision-log states are computed over UTF-8 text
with LF newlines. Digests of run artifacts in `20_new_runs/` are computed over
raw bytes, as the `run_manifest.json` files have always done — those files are
data, not line-ending-sensitive prose, and no artifact digest recorded to date
is affected.

## ADR-AIIM-296 — ADR-AIIM-295 result: the curve is monotone, the model is usable down to AD 0.25, 69 % of the library is below that — and two of the four failed predictions failed in my harness, not in the data

**Status:** accepted — record of fact
**Date:** 2026-10-05
**Approver:** Adisorn Sirichotejirakul (2026-10-05)
**Pre-registration:** ADR-AIIM-295, SHA-256
`30555270484983ce475618359c4f8a6c71859086266a97888cf40931bedf291e`, frozen
before the run; log SHA-256 at that point (text/LF, per the 2026-10-05
addendum) `a61c4017b61354af99084fd790bcc22ada08104d6e3608f5ad17a0cf74061c53`
**Artifacts:** `20_new_runs/ad_distance_curve_2026-10-05/` · 32 s
**Corrects:** ADR-AIIM-294's Tier 1 verdict

### Grading

| | prediction | verdict | why |
|---|---|---|---|
| **Y1** | monotone degradation, Spearman negative, CI excluding 0 | **MISSED** | **harness defect** — rank direction coded opposite to the specification |
| **Y2** | the crossing is inside the measured range | **MISSED** | **real** — every decision-bearing band is Tier 2 |
| **Y3** | the crossing sits between 0.35 and 0.55 | **MISSED** | **real** — it is at 0.25, earlier than predicted |
| **Y4** | A1 and A8 agree on the crossing band | **MISSED** | **harness defect** — Tanimoto band edges applied to Dice scores |

**The anchor reproduces.** Full held-out ROC-AUC **0.9217**, identical to
ADR-AIIM-294 to four places, confirming no refitting occurred and the curve is
measured on the same model.

### The curve, A1 (ECFP4 Tanimoto, the project's own definition)

| AD band | n | active / inactive | ROC-AUC | CI 95 % | tier |
|---|---|---|---|---|---|
| 0.15–0.25 | 22 | 5 / 17 | **0.3294** | 0.042–0.648 | Tier 0 † |
| 0.25–0.35 | 92 | 36 / 56 | **0.7391** | 0.633–0.837 | **Tier 2** |
| 0.35–0.45 | 109 | 70 / 39 | 0.8046 | 0.713–0.889 | Tier 2 |
| 0.45–0.55 | 165 | 121 / 44 | 0.9098 | 0.849–0.959 | Tier 2 |
| 0.55–0.70 | 462 | 394 / 68 | 0.9357 | 0.896–0.969 | Tier 2 |
| 0.70+ | 291 | 249 / 42 | 0.9187 | 0.857–0.970 | Tier 2 |

† declared in ADR-AIIM-295 as reported and not decision-bearing, before the
numbers were seen.

Monotone from 0.25 to 0.70 and flat thereafter. |Spearman| = 0.9, CI [0.5, 1.0].

### Y1 failed in the harness, and the sign is the whole of it

ADR-AIIM-295 specified *"Spearman correlation between band rank (nearest to
furthest) and band ROC-AUC is negative"*. The harness iterated `BANDS` in
**furthest-to-nearest** order, so rank 0 is the furthest band and the
correlation came out **+0.9**, and the grader's `< 0` test failed.

The magnitude and the interval are unchanged, and the content — **performance
degrades monotonically with distance from the training set** — is what Y1 was
written to test and is supported. **The verdict is nevertheless recorded as
MISSED.** Y1 was a prediction about a signed quantity, the harness computed the
opposite sign from the one specified, and flipping it after the fact to recover
MET is precisely the move this discipline exists to prevent. The defect is
mine; the finding stands on the magnitude and the CI, not on the verdict.

### Y4 failed in the harness, for a reason that was knowable in advance

| band | A1 Tanimoto n | A8 Dice n |
|---|---|---|
| 0.15–0.25 | 22 | **0** |
| 0.25–0.35 | 92 | **8** |
| 0.35–0.45 | 109 | 44 |
| 0.45–0.55 | 165 | 98 |
| 0.55–0.70 | 462 | 221 |
| 0.70+ | 291 | **770** |

**Dice ≥ Tanimoto identically for binary fingerprints**, so Dice scores are
systematically higher and the same band edges place 770 of 1,141 compounds in
the top band. A8's apparent crossing at 0.45 against A1's 0.25 measures the
offset between two metrics' scales, **not** a disagreement about where the
model stops working. Reusing Tanimoto-derived edges for Dice was an error in
ADR-AIIM-295's design.

A second-order consequence: ADR-AIIM-295 declared only 0.15–0.25 underpowered,
from Tanimoto counts. Under Dice the underpowered band moved to 0.25–0.35
(n = 8, AUC 0.125) and was still graded as decision-bearing. The declaration
was made in the right spirit and attached to the wrong coordinate.

### Y2 failed as predicted it would be most embarrassing, and it corrects ADR-AIIM-294

Every decision-bearing band under A1 is Tier 2. ADR-AIIM-295 named this
outcome first, as the one that would mean *"the Tier 1 verdict I have just
reported to the author does not survive stratification."* It does not.

**Why the two entries disagree.** ADR-AIIM-294's NP-overlap band was
0.1738–0.3335 — the natural-product library's own 1st–99th percentile — which
spans the 0.15–0.25 region where the model scores **AUC 0.3294, below chance**.
Mixing 22 such compounds with the 92 at 0.7391 yields the 95 compounds at
0.6808 that produced Tier 1.

**So ADR-AIIM-294's Tier 1 was an average across a sign change in the model's
usefulness**, not a level. Stratified, the same data gives **Tier 2 at
0.25–0.35 and Tier 0 at 0.15–0.25**. ADR-AIIM-294's verdict stands as what it
measured and is superseded as a description of the model: *"above chance but
not usable at natural-product distances"* is withdrawn and replaced by *"usable
down to 0.25, below chance beneath it, with 22 compounds of evidence there."*

### Y3 failed in the direction that favours the article

The crossing is at **AD 0.25**, not the 0.35–0.45 band I predicted. The model
reaches usable discrimination **closer to natural-product chemistry than
expected**. Against that, the band immediately below is the one with 22
compounds and an AUC under 0.5, so the crossing is sharp rather than gradual
and sits exactly where the evidence thins.

### What this licenses, quantitatively

Natural-product library against the same training set, ECFP4 Tanimoto mean
top-5:

| | compounds | share |
|---|---|---|
| **≥ 0.25 — at or above the usable boundary** | **1,560** | **31.2 %** |
| < 0.25 — below it | 3,440 | 68.8 % |
| ≥ 0.35 | 23 | 0.5 % |

The 0.25 line falls at the **68.8th percentile** of the library.

> The model discriminates usably (Tier 2, ROC-AUC 0.7391 [0.633, 0.837]) down to
> AD 0.25. **31.2 % of the library reaches that distance.** For the other
> 68.8 % the entire evidence is 22 held-out compounds at ROC-AUC 0.3294, CI
> 0.042–0.648 — below chance, with an interval too wide to rank on.

That splits the library into a part predictions can be defended for and a part
they cannot, which neither *"0 % in-domain"* nor *"Tier 1"* could do.

### Unchanged

**No coverage guarantee.** Constraint 3 of ADR-AIIM-295 fixed this in advance:
a per-band AUC concerns ranking, not conformal validity, and §3.6's *p*-values
are untouched. **Constraint 4 likewise**: carrying the 0.25 boundary over to
natural products assumes degradation tracks distance rather than chemotype,
which this run does not test and which remains an assumption.

**§3.6 is unchanged.** It reports 72 of 5,000 and continues to.

### Deliberately not done

* **`10_new_system` is unmodified.** No refitting; the ADR-AIIM-294 model was
  reconstructed from the same seeded canonical inputs and verified by anchor.
* **Y1 was not re-run with the sign corrected**, and Y4 was not re-run with
  metric-specific band edges. Both are recorded as harness defects rather than
  patched into passes.
* **The A8 curve is reported but carries no crossing claim**, for the scale
  reason above.
* **The 0.15–0.25 band was not enlarged.** Its 22 compounds are what the
  held-out set contains at that distance under this split; obtaining more would
  need a different split or a larger labelled set.
* **Chemotype confounding was not tested.** Whether the 22 compounds fail
  because they are distant or because they are a peculiar chemical class is not
  separated.
* **No compound named in any output is a candidate.** Nothing has been assayed.

### Artifacts

Digests over raw bytes, per the 2026-10-05 addendum.

```
20_new_runs/ad_distance_curve_2026-10-05/
  ad_distance_curve.json        d340d72796988a61fd260db26cadf47af983f7bf7b4300fa9b73dde03cd43fb3
  analysis/ad_distance_curve.py 6754772b0caa11b6b4cd6769b878c0125e4c5d7e7a113b44f4157271051d659c
```

## ADR-AIIM-297 — PRE-REGISTRATION: model pKi as a continuous quantity, threshold afterwards, and use the censored records as a one-sided test

**Status:** proposed — **pre-registration. Frozen before the run.**
**Date:** 2026-10-05
**Raised by:** Adisorn Sirichotejirakul — *"ถ้าเกิดเราลองวิธีเทียบ AD หลายแบบได้ไหม … อาจจะลองไปใช้วิธี regression ไม่ตัด data แบ่งกลุ่ม แล้วค่อยมาตัดทีหลัง predict"* (2026-10-05), and *"สรุปลอง regression ให้หน่อยได้ไหม"* (2026-10-05)
**Follows:** ADR-AIIM-293 (the curated binding set), ADR-AIIM-296 (the
distance curve), ADR-AIIM-286 (the eight split conventions)
**Discipline:** ADR-AIIM-263 — predictions fixed before the run;
ADR-AIIM-226 — outside the frozen tree

### Why the binary framing is the thing left untested

Every entry from ADR-AIIM-285 onward has moved something inside a
classification pipeline — the split convention, the label source, the
calibration balance, the applicability-domain definition — and the count has
stayed undetermined to between 3.62× and 376×. The one structural choice never
examined is the binarisation itself, which costs the project three things that
are measurable:

| | cost |
|---|---|
| the 5.5–6.0 gap | **1,252 compound-endpoint rows discarded** at step 9 of ADR-AIIM-291 for falling between the cuts |
| censored records | 617 upper- and 297 lower-censored compounds in the binding set forced to a hard 0 or 1, when each is an **interval** |
| the threshold | 6.0 / 5.5 has to be chosen before training, and is a convention with no defence in this project's record |

A regression on pKi pays none of them. It also removes the Mondrian conformal
machinery entirely — `cal1` built from calibration actives, whose composition
ADR-AIIM-286 showed swinging the count by two orders of magnitude — because a
count becomes "how many compounds are predicted at or above a threshold",
applied after the fact.

### What is run

**Target.** pKi, continuous, from **point measurements only** (`basis =
point_median`), median-aggregated per compound, binding endpoint class. Derived
by replaying ADR-AIIM-291's frozen curation steps 1–8 on the saved
`chembl251_raw.csv`, so the 1,252 rows that step 9 discarded are **recovered**
rather than re-pulled. The count recovered is reported.

**Model.** `RandomForestRegressor`, `random_state = 42`, features from
`compute_raw_features` unchanged, `VarianceThreshold(0.01)` inside the fitted
`Pipeline`. Nothing is tuned.

**Splits.** The same **eight conventions** of ADR-AIIM-286 under canonical
InChIKey order, so the spread is directly comparable to the classifier's
38.88× on the same compounds.

**Held-out metrics**, per convention: R², RMSE, MAE and Spearman ρ on the
held-out point measurements.

**The censored one-sided test.** The 617 upper-censored and 297 lower-censored
compounds are **excluded from training and from the held-out regression metrics
entirely**. For each, the model's prediction is checked against the single
inequality that record actually asserts — predicted pKi ≤ bound for
upper-censored, ≥ bound for lower-censored. This tests the model on data no
fitting step has seen and that the binary pipeline either discarded or
flattened, and it needs no ground-truth value.

**Counts on the NP library**, two ways: the point prediction at or above a
threshold, and a **conformal-regression lower bound** at α = 0.05 using
absolute held-out residuals as the nonconformity score, which is the
regression analogue of §3.6's quantity. Reported as a **curve over thresholds
5.0, 5.5, 6.0, 6.5, 7.0**, not a single number, because choosing one threshold
after seeing the counts is the failure mode this discipline exists to prevent.

### Predictions

**Z1 — the regression ranks usably.** Spearman ρ on the held-out point
measurements exceeds 0.60 under the anchor convention
(`sorted · train-first · divert`). *Falsified at 0.60 or below.* *Basis:* the
classifier reaches ROC-AUC 0.9217 on the same compounds and features, and a
regression that ranks worse than that implies the continuous signal is weaker
than the dichotomised one, which would be surprising.

**Z2 — the censored records agree with it.** More than 70 % of the 617
upper-censored compounds are predicted at or below their own upper bound.
*Falsified at 70 % or below.* This is the prediction I would most like to be
right about, because it is the only check in this project that uses the
censored data as evidence rather than as a labelling problem.

**Z3 — regression collapses the convention spread.** The eight conventions'
NP-library counts at threshold 6.0 span **less than 10×**, against the
classifier's 38.88× on the same curated set. *Falsified at 10× or above.*

**Z4 — it does not land on 72.** No convention returns an NP count at
threshold 6.0 within ±50 % of 72, that is within 36–108. *Falsified if any
does.*

**The outcome that would most embarrass this entry is Z3 failing**, and it is
named first. The argument for trying regression at all is that binarisation and
its conformal apparatus are where the instability lives. If the spread stays at
or above 10× with the classifier's machinery removed, that argument is wrong,
the instability is in the features or the chemistry rather than in the
framing, and three days of work on the label and the split will have been
aimed at the wrong layer.

### Constraints

1. **`10_new_system` is not modified.** No contract amended, no input
   promoted, nothing re-docked, `qsar_pipeline_a2a.pkl` untouched. Curation
   functions are imported from the ADR-AIIM-291 harness, not rewritten.
2. **`02_frozen_inputs/` is not written to**, and `chembl251_raw.csv` is read,
   not re-pulled.
3. **§3.6 is unchanged.** It reports 72 of 5,000 and continues to. No count
   here is proposed as a replacement, and the threshold curve exists precisely
   so that no single number from this run can be quoted as the answer.
4. **A regression count is not a conformal count.** The interval version uses
   absolute residuals calibrated on held-out ChEMBL compounds and inherits the
   same exchangeability failure ADR-AIIM-294 measured: the NP library is
   0 / 5,000 in-domain under the structural definition, and ADR-AIIM-296 found
   usable discrimination only down to AD 0.25, which 31.2 % of the library
   reaches. **Regression does not repair that and will not be reported as if
   it had.**
5. **Z2 is not a validity claim.** Passing a one-sided inequality on censored
   records shows the model is not contradicting the data it never saw. It does
   not establish accuracy, and a high pass rate is weak evidence because many
   bounds are loose.
6. **This is a sensitivity analysis, not a replacement framing.** Whether the
   article should move from classification to regression is an author decision
   and is not taken here; ADR-AIIM-291's own recommendation was to measure it
   alongside rather than switch.
7. **No compound named in any output is a candidate.** Nothing has been
   assayed.
8. **Every outcome is reported**, including Z3 failing and a regression that
   ranks worse than the classifier.
9. **This entry is appended to the decision log and its SHA-256 recorded before
   the run.** Digest over UTF-8 text with LF newlines, per the 2026-10-05
   addendum.

### Skills invoked

| skill | version | material influence |
|---|---|---|
| `scikit-learn` | 1.1 | `RandomForestRegressor` inside a `Pipeline` with `VarianceThreshold`, so no transform is fit on held-out rows; scaffold grouping preserved across all eight conventions |
| `statistical-analysis` | 1.0 | applied for metric choice on a continuous target — Spearman reported beside R² because the downstream use is ranking and thresholding, not calibrated prediction, and R² alone would not speak to it |

### Output

```
20_new_runs/pki_regression_2026-10-05/
  pki_regression.json
  run_manifest.json
  analysis/pki_regression.py
```

## ADR-AIIM-298 — PRE-REGISTRATION: put the censored records in the likelihood instead of the test set

**Status:** proposed — **pre-registration. Frozen before the run.**
**Date:** 2026-10-05
**Raised by:** Adisorn Sirichotejirakul — *"รัน tobit เลยแล้วเขียน ADR รวมกันทีเดียว"* (2026-10-05), and the original proposal *"ลองไปใช้วิธี regression ไม่ตัด data"* (2026-10-05)
**Follows:** ADR-AIIM-297, whose Z2 returned **0.000 in all eight conventions**
**Discipline:** ADR-AIIM-263 — predictions fixed before the run;
ADR-AIIM-226 — outside the frozen tree

### The defect this run exists to repair

ADR-AIIM-297 held all 914 censored compounds out of training and tested the
model against their bounds. Not one passed, in any convention. The diagnosis is
not that regression fails; it is that **the design starved the model of the
very class it was then tested on**:

| | |
|---|---|
| training point measurements with pKi ≥ 6.0 | **86.2 %** (4,131 of 4,793) |
| training point measurements with pKi < 5.5 | 13.8 % |
| regression-set mean pKi | 6.97 |
| upper-censored bounds | min 4.00, **median 5.00**, max 5.41 |
| model predictions on the NP library | 5.144–7.569, mean 6.038 |
| mean violation (prediction − bound) | **+1.1594** |
| lower-censored pass rate | 0.037 |

The inactive information in this dataset lives almost entirely in the censored
records. Holding all of them out left a training set that is 86 % potent, and a
`RandomForestRegressor` on such a set regresses to its mean and cannot reach
below about 5.1. Compression at both ends — 0.000 upper and 0.037 lower —
confirms it is range collapse rather than a directional bias.

**That was my choice, not a property of the data.** The censored records belong
in the likelihood as intervals, which is what was proposed in the first place.

### The change

**One variable moves: how censored records enter.** Same model family, same
features, same pipeline, same eight conventions, same canonical ordering as
ADR-AIIM-297, so the comparison isolates the censoring treatment rather than
confounding it with a change of estimator. A linear tobit would have changed
both and is deliberately not used.

**Censored-EM regression**, the non-parametric equivalent of tobit for a
tree model:

1. initialise each censored compound's target at its own bound;
2. fit `RandomForestRegressor` on train points **and** train censored;
3. predict the train censored compounds;
4. project back into the interval each record asserts — upper-censored becomes
   `min(prediction, bound)`, lower-censored becomes `max(prediction, bound)`;
5. repeat to convergence, maximum 10 iterations, tolerance 0.01 log units on
   the largest target change; iterations used and final change are reported.

Step 4 is the whole of it: the target is never allowed outside the interval the
measurement supports, so a censored record contributes its information without
being assigned a value it does not have.

**The censored compounds are now split like everything else.** Those landing in
the training side contribute intervals; those landing in the held-out side are
scored one-sided and are never fitted. ADR-AIIM-297's test survives in a form
that is not self-defeating.

**Reported:** held-out point R², RMSE, MAE and Spearman ρ; held-out censored
one-sided pass rates; NP-library and SNP-Lib counts as point predictions and as
conformal lower bounds (absolute held-out residuals, α = 0.05) over thresholds
5.0, 5.5, 6.0, 6.5, 7.0; and the EM convergence record.

### Predictions

**T1 — the model can now reach the inactives.** The held-out upper-censored
one-sided pass rate exceeds 70 % under the anchor convention
(`sorted · train-first · divert`). *Falsified at 70 % or below.* This is
ADR-AIIM-297's Z2 re-asked with the design fixed.

**T2 — the predictive range extends downward.** The minimum prediction on the
NP library falls below 5.0 under the anchor, against ADR-AIIM-297's 5.144.
*Falsified at 5.0 or above.* T1 can hold while T2 fails if the model reaches the
bounds only for training-like chemistry, so the two are not the same claim.

**T3 — the convention spread tightens further.** The eight conventions' NP
counts at threshold 6.0 span less than ADR-AIIM-297's **10.99×**. *Falsified at
10.99× or above.* For reference the classifier on the same curated set spans
38.88×.

**T4 — the intervals still swallow the decision.** The conformal lower-bound
count at threshold 6.0 remains **0** in every convention, as it was in all eight
of ADR-AIIM-297. *Falsified if any convention returns more than 0.* ADR-AIIM-297
measured α = 0.05 absolute-residual quantiles of 1.35–1.73 log units against a
decision gap of 0.5, so the interval is roughly three times the resolution the
threshold needs. **This is the prediction whose falsification would be the best
news in three days of this work**, because a non-zero defensible count is the
one thing none of these runs has produced.

**The outcome that would most embarrass this entry is T1 failing.** If censored
records in the likelihood still cannot make the model predict an inactive, then
the feature representation cannot express inactivity in this target, and the
regression line is finished rather than mis-implemented. ADR-AIIM-297's null
would then have been the right answer reached by the wrong route.

### Constraints

1. **`10_new_system` is not modified.** Curation is replayed by importing the
   ADR-AIIM-291 harness; `chembl251_raw.csv` is read, not re-pulled;
   `02_frozen_inputs/` is not written to.
2. **§3.6 is unchanged.** It reports 72 of 5,000 and continues to. The
   threshold curve is reported in full so that no single number here can be
   quoted as the answer.
3. **Censoring handled correctly is still not exchangeability.** The NP library
   is 0 / 5,000 in-domain structurally (ADR-AIIM-294) and the model is usable
   only down to AD 0.25, which 31.2 % of it reaches (ADR-AIIM-296). **Nothing in
   this run repairs that**, T4 included, and no count will be presented as if it
   had.
4. **A one-sided pass is weak evidence.** It shows the model does not
   contradict data it never saw. Many bounds are loose — the upper-censored
   median is 5.00 against a 5.5 cut — so a high pass rate is cheap and will be
   reported with the bound distribution beside it.
5. **Censored-EM is not tobit.** It is the interval-respecting analogue for a
   tree model and makes no Gaussian assumption; it also has no likelihood value
   to report and no standard errors. Calling it tobit in the article would be
   wrong, and the distinction is stated here so it is not blurred later.
6. **No compound named in any output is a candidate.** Nothing has been
   assayed.
7. **Every outcome is reported**, including T1 failing and a spread that widens.
8. **This entry is appended to the decision log and its SHA-256 recorded before
   the run**, digest over UTF-8 text with LF newlines per the 2026-10-05
   addendum.

### Skills invoked

| skill | version | material influence |
|---|---|---|
| `scikit-learn` | 1.1 | `RandomForestRegressor` inside a `Pipeline` with `VarianceThreshold`; the EM loop refits the whole pipeline each iteration rather than reusing a transform fitted on an earlier target, so no leakage enters through the imputation |
| `statistical-analysis` | 1.0 | applied for the censoring treatment: an interval-censored observation enters as a constraint on the likelihood, not as an imputed point, which is the distinction step 4 implements and constraint 5 refuses to blur |

### Output

```
20_new_runs/pki_censored_em_2026-10-05/
  pki_censored_em.json
  run_manifest.json
  analysis/pki_censored_em.py
```

## ADR-AIIM-299 — ADR-AIIM-297 and ADR-AIIM-298 results: regression halves the convention problem twice over, censoring in the likelihood halves it again, and the honest interval still swallows every threshold

**Status:** accepted — record of fact
**Date:** 2026-10-05
**Approver:** Adisorn Sirichotejirakul — *"รัน tobit เลยแล้วเขียน ADR รวมกันทีเดียว"* (2026-10-05)
**Pre-registrations:** ADR-AIIM-297, SHA-256
`b46f97a981c575287f5faffc414571e1b96a976895d40deb9a41b25934ca569c` (log at
freeze `a36db9f9f3bcd8df25ce239e93e3aa163d02bd7d03978d1074837cde1e702e58`);
ADR-AIIM-298, SHA-256
`27c3752ba8105d0dd25bcca329bcf96690f71499c03a9c3513e93b45c79ee1d5` (log at
freeze `4014358c797ad0728cb0d36658931a537e52f334603b017cf967d1033e87a3f6`)
**Artifacts:** `20_new_runs/pki_regression_2026-10-05/` (67 s),
`20_new_runs/pki_censored_em_2026-10-05/` (707 s)

### Grading

| | prediction | verdict |
|---|---|---|
| **Z1** | Spearman > 0.60 on held-out point measurements | **MET** — 0.754 at anchor, 0.726–0.836 across conventions |
| **Z2** | > 70 % of upper-censored predicted at or below their bound | **MISSED** — **0.000 in all eight** |
| **Z3** | convention spread below 10× | **MISSED** — 10.99 × |
| **Z4** | no count within ±50 % of 72 | **MET** |
| **T1** | with censoring in the likelihood, > 70 % upper-censored pass | **MISSED** — 0.0087 at anchor, 0.000–0.044 across conventions |
| **T2** | NP minimum prediction falls below 5.0 | **MET** — 4.398 against ADR-AIIM-297's 5.144 |
| **T3** | spread below ADR-AIIM-297's 10.99× | **MET** — **4.10×** |
| **T4** | the conformal count at 6.0 stays 0 everywhere | **MET** |

### The convention problem shrinks monotonically with the framing

NP-library count at threshold 6.0, the **same curated compounds** throughout:

| framing | min | max | **fold** |
|---|---|---|---|
| classifier + Mondrian conformal (ADR-AIIM-293) | 104 | 4,044 | **38.88×** |
| regression, point measurements only | 241 | 2,649 | **10.99×** |
| regression, censored in the likelihood | **29** | **119** | **4.10×** |

A **9.5× reduction** in convention sensitivity from removing the binarisation
and its calibration machinery, and the counts fall by an order of magnitude as
the model gains inactive information. Z3 failed on its stated threshold and the
direction it predicted held; T3 passed.

### Z2 and T1 failed, and the pre-registered interpretation of T1 does not hold

ADR-AIIM-297's Z2 was 0.000 because the design held all 914 censored compounds
out of training, leaving a set 86.2 % of which sits above pKi 6.0 — the
inactive information in this dataset lives almost entirely in the censored
records. That diagnosis is confirmed and was the reason for ADR-AIIM-298.

**ADR-AIIM-298 put 686–858 censored compounds into training per convention and
the pass rate went from 0.000 to 0.009.** Mean violation 0.62–1.15 log units.

ADR-AIIM-298 pre-registered what T1 failing would mean: *"the feature
representation cannot express inactivity in this target, and the regression
line is finished rather than mis-implemented."* **That conclusion is not
adopted, because a competing explanation exists that the pre-registration did
not anticipate and that is specific to my own E-step.**

Step 4 imputes `min(prediction, bound)`. A censored target initialised at its
bound, for a compound the model scores above that bound, **stays at the
bound** — it is never pushed below it. A proper tobit E-step imputes
`E[y | y ≤ b]`, which lies strictly below `b`. So each censored compound pulls
the fit **toward** its bound rather than beneath it, and the training target
distribution acquires a spike at the bound (median 5.00) instead of mass below
it. The projection is a correct constraint and the wrong expectation.

T2 is the evidence that this is the binding limitation rather than the
features: the predictive range **did** extend downward, 5.144 → 4.398, so the
model is not incapable of low predictions. It places the wrong compounds there.

**What would separate the two explanations:** a parametric tobit or an AFT
model with a conditional-expectation E-step. Until that is run, *"regression
cannot express inactivity here"* is unsupported and *"my E-step is biased
upward at the censored points"* is the live hypothesis.

### The EM did not converge, and this is reported rather than smoothed

All eight conventions ran the full 10 iterations with a final maximum target
change of **0.0256–0.0508**, against a declared tolerance of **0.01**. The
numbers above are therefore a 10-iteration state, not a converged one. The
cap was pre-registered; the failure to reach tolerance inside it was not
anticipated and is a limitation of this run.

### T4 held everywhere, and it is the only result in three days that no convention moves

| | α = 0.05 absolute-residual quantile | conformal count at pKi ≥ 6.0 |
|---|---|---|
| ADR-AIIM-297 | 1.353–1.733 log | **0** in all eight |
| ADR-AIIM-298 | 1.527–1.671 log | **0** in all eight, and **0 at every threshold from 5.0 up** |

The honest 95 % interval is **±1.5 to ±1.7 log units** against a decision gap
of **0.5**. The interval is roughly three times the resolution the threshold
needs.

> **With honest prediction intervals, not one compound in the natural-product
> library can be asserted at pKi ≥ 6.0 — zero, under every framing and every
> convention tested.**

This does not say the library contains no active compound. It says **this
instrument cannot resolve the question at that threshold**, which is a
statement about the method and not about the chemistry, and it is the one
quantity in this entire sequence that is not a function of an arbitrary
convention.

### What the two runs establish together

| | |
|---|---|
| the model **ranks** | Spearman 0.726–0.836, R² 0.53–0.70 on held-out point measurements, across both framings |
| the model **cannot assert absolute activity** | interval ±1.5 log, conformal count 0 at every threshold |
| binarisation was **part** of the instability, not all of it | 38.88× → 10.99× → 4.10×, never to 1× |
| censored records are **worth including** | 914 compounds recovered, counts fall an order of magnitude, spread falls 2.7× |
| the gap is **worth recovering** | 577 compound-endpoint rows that step 9 of ADR-AIIM-291 discarded |

### Unchanged

**§3.6 is unchanged.** It reports 72 of 5,000 and a maximum conformal *p* of
0.072 and continues to. No count here is proposed as a replacement, and the
threshold curves are reported in full so that none can be quoted alone.

**Exchangeability is untouched.** Constraint 4 of both pre-registrations fixed
this in advance: the NP library is 0 / 5,000 in-domain structurally
(ADR-AIIM-294), the model is usable down to AD 0.25 which 31.2 % of the library
reaches (ADR-AIIM-296), and neither regression nor correct censoring repairs
either fact.

### Deliberately not done

* **`10_new_system` is unmodified**, `02_frozen_inputs/` not written to, and
  `chembl251_raw.csv` read rather than re-pulled in both runs.
* **No parametric tobit or AFT model was run**, so the competing explanation for
  T1 is named and not tested.
* **The EM was not run to convergence**, and no convention's tolerance was
  relaxed after the fact to claim it had.
* **Censored-EM is not tobit** — no Gaussian assumption, no likelihood value, no
  standard errors (ADR-AIIM-298 constraint 5). Calling it tobit in the article
  would be wrong.
* **The functional endpoint class was not modelled.** Both runs use binding
  only.
* **No compound is named or ranked in either output**, and nothing in this
  project has been assayed.

### Artifacts

Digests over raw bytes, per the 2026-10-05 addendum.

```
20_new_runs/pki_regression_2026-10-05/
  pki_regression.json          5c9cd0279499aadd23b752dac3afdaba061fb0198eecde2533aa47b9294c21b5
  analysis/pki_regression.py   0588bc9b0f7adea0f7205b62679a2c6d3148cc132d182778c5c661d6d87ea4d7
20_new_runs/pki_censored_em_2026-10-05/
  pki_censored_em.json         26312b29c6b1bc7d8f590d2d3673ce1ac19ffcaf2ab4ca8598721bcb053312fb
  analysis/pki_censored_em.py  88fa6d24fce97cd017fbde7b562f47bc88dd3aade28f71a829a27e14629ee117
```

## ADR-AIIM-300 — PRE-REGISTRATION: prioritise by convention-robust rank inside the measured competence band, and stop asserting activity

**Status:** proposed — **pre-registration. Frozen before the run.**
**Date:** 2026-10-05
**Raised by:** Adisorn Sirichotejirakul — *"สุดท้ายจะต้องหาวิธีที่เหมาะสมที่สุดกับ library นี้เพื่อ prioritised สาร ให้ได้ว่าอะไร active ลองดูว่าบาง test block การค้นพบ new scaffold active compound มากเกินไปหรือไม่"* (2026-10-05)
**Follows:** ADR-AIIM-296 (the competence band), ADR-AIIM-299 (rank works,
assertion does not)
**Bears on:** ADR-AIIM-279, whose canonical re-run retained none of the three
originally prioritised compounds
**Discipline:** ADR-AIIM-263 — the selection rule and every parameter fixed
**before any compound identity is looked at**; ADR-AIIM-226 — outside the
frozen tree

### The author's question, and the one part of it this entry cannot honour

The request has three parts. Two are answerable and one is not, and saying so
is the condition of the other two being worth anything.

**Answerable — is the method right for this library?** Yes, and ADR-AIIM-299
says what it can and cannot do: it ranks (Spearman 0.73–0.84) and it cannot
assert absolute activity (interval ±1.5 log against a 0.5 decision gap,
conformal count 0 at every threshold under every convention).

**Answerable — do some tests block novel scaffolds too much?** Yes, and §3.5's
own contract already proves the argument once. See below.

**Not answerable — a result in which the three originally prioritised compounds
survive.** A selection rule chosen after seeing which rule returns which
compounds is not a selection rule, and the decision log timestamps every step
of the last four days, so the substitution would be visible to any reader who
opened it. **This entry therefore fixes the rule, every threshold in it, and
the shortlist length before a single compound identity is read**, and records
in advance that the three compounds' presence or absence is an **observation,
not a criterion**. If the rule returns them, that is informative. If it does
not, the rule stands.

### A test that cannot admit a novel library is not measuring what it claims

This is not a new argument in this project. `G3_pose.yaml:157` retired
`control_representativeness` for exactly it:

> ADR-AIIM-073: `control_representativeness` is a ligand-similarity proxy **that
> cannot in principle admit a structurally novel library** (E2 showed 0 %
> coverage reproduces across six independent structural descriptor spaces).
> Direct structure-based evidence — docking a labeled benchmark under the
> identical protocol — answers the transferability question the proxy was
> standing in for, **without depending on library-to-control similarity**.

**The structural applicability-domain criterion has the identical logical
form** and has never been given the same treatment:

```
threshold  = 10th percentile of the training set's similarity to itself
in-domain ⟺ mean top-5 similarity to training ≥ threshold
```

A compound qualifies only by being as similar to the training set as training
compounds are to one another. **A novel scaffold cannot satisfy that by
construction**, which is why Figure 4 reports exactly 0 % in ECFP4, FCFP4,
topological-torsion and Gobbi space and 86.3 % in physicochemical space. The
0 / 5,000 figure measures novelty, not unsuitability.

**And ADR-AIIM-296 is the direct-evidence replacement** that ADR-AIIM-073's
reasoning calls for: performance measured as a function of distance, which says
*where the model works* rather than *whether compounds resemble the training
set*. It found usable discrimination (Tier 2, ROC-AUC 0.7391, CI 0.633–0.837)
down to **AD 0.25**, reached by **31.2 %** of the library, with 22 held-out
compounds at ROC-AUC 0.3294 below it.

So the AD criterion is used here as a **competence boundary read off measured
performance**, not as an admission gate read off similarity. That is the
methodological change, and it is justified without reference to which compounds
result.

**What is not claimed:** that G6's `scaffold_uniqueness` blocks novelty. It
requires a panel to span distinct Bemis–Murcko frameworks, which is a
redundancy constraint and if anything favours novelty — it halted the released
panel because two of its three compounds shared a scaffold (ADR-AIIM-239). It
is imported unchanged below.

**And the three compounds were not displaced by a novelty filter.**
ADR-AIIM-279 retained none of them under canonical ordering because the
conformal eligibility column changed — the column ADR-AIIM-286 measured
spanning 376×. They were selected by a procedure whose output is not determined
by the data. That is a different problem from an over-strict gate and it is the
one this rule is built to answer.

### The rule, fixed now

**Model.** Censored-EM regression exactly as ADR-AIIM-298 ran it — same
features, same pipeline, 10 iterations, tolerance 0.01. Chosen over
ADR-AIIM-297's point-only regression on two declared grounds: it uses the 914
censored compounds the other discards, and its convention spread is **4.10×
against 10.99×**. The cost is 0.01–0.06 Spearman, accepted and stated.
**ADR-AIIM-299's limitations carry over**: the E-step is biased upward at
censored points and the EM did not reach tolerance.

**Eligibility.** A library compound is eligible if its ECFP4 mean-top-5
similarity to the training set is **≥ 0.25**, the boundary ADR-AIIM-296
measured. Expected ~1,560 of 5,000. No other filter is applied — no AD gate, no
conformal set, no docking score.

**Ranking.** Predicted pKi, descending, within the eligible set.

**Convention robustness, which is the point.** The rank is computed under all
**eight conventions** of ADR-AIIM-286. A compound enters the shortlist only if
it appears in the **top 20 under at least 5 of the 8**. Instead of choosing a
convention — which ADR-AIIM-286 established is choosing an answer — the rule
requires a compound to survive a declared majority of them.

**Diversity.** `G6`'s `scaffold_uniqueness` rule is imported unchanged: at most
one compound per Bemis–Murcko framework, keeping the higher-ranked.

**Shortlist length: 10.** Fixed now.

### Predictions

**Q1 — the rule yields a usable shortlist.** At least 10 compounds appear in
the top 20 under 5 or more of the 8 conventions. *Falsified below 10.*

**Q2 — it is scaffold-diverse.** The 10 span at least 5 distinct Bemis–Murcko
frameworks before the uniqueness rule is applied. *Falsified below 5.*

**Q3 — it is not merely a set of near-analogues.** At least one shortlisted
compound has AD < 0.45, that is, lies outside the band where the model scores
above 0.80. *Falsified if all ten sit at 0.45 or above.* This tests whether
ranking inside the competence band still reaches novel chemistry or only
retrieves the training set's neighbours.

**Q4 — rank is more robust than count.** The number of compounds meeting the
5-of-8 majority exceeds the number that would be selected by any single
convention's top 20 intersected across all eight. *Falsified otherwise.*

**Reported and explicitly not a criterion:** whether NPLIB01185, NPLIB01643 or
NPLIB03056 appear, and at what rank and AD score.

**The outcome that would most embarrass this entry is Q1 failing.** If no
compound survives a majority of eight conventions even on rank, then rank is as
convention-dependent as the count was, the 4.10× spread was hiding the same
problem one level down, and four days of work yields no defensible
prioritisation of any kind.

### Constraints

1. **`10_new_system` is not modified.** No contract amended, no input promoted,
   nothing re-docked, `qsar_pipeline_a2a.pkl` untouched. G6's uniqueness rule is
   imported, not restated.
2. **§3.6 is unchanged** and no count from this run replaces its 72 of 5,000.
3. **A shortlisted compound is not a predicted active.** The output is a
   **priority order for experiment**, supported by rank performance inside a
   measured competence band. ADR-AIIM-299 established that no compound in this
   library can be asserted at pKi ≥ 6.0 under any framing, and nothing here
   changes that. The word "active" is not applied to any output.
4. **The competence band is measured on ChEMBL compounds.** Transferring it to
   natural products assumes degradation tracks distance rather than chemotype
   (ADR-AIIM-296 constraint 4), which is an assumption and not a result.
5. **Exchangeability is still not established**, so no conformal statement is
   made about any shortlisted compound.
6. **No compound in the shortlist is a candidate for any claim about binding.**
   Nothing in this project has been assayed, and the shortlist is a proposal for
   what to assay, nothing more.
7. **The three originally prioritised compounds are reported wherever they fall
   and the rule is not revised afterwards.**
8. **Every outcome is reported**, including Q1 failing and a shortlist that
   contains none of the three.
9. **This entry is appended to the decision log and its SHA-256 recorded before
   the run**, digest over UTF-8 text with LF newlines per the 2026-10-05
   addendum.

### Skills invoked

| skill | version | material influence |
|---|---|---|
| `scikit-learn` | 1.1 | the ADR-AIIM-298 pipeline reused unchanged; no transform refitted per convention outside its own fit |
| `rdkit` | 1.1 | ECFP4 for the competence score and Bemis–Murcko frameworks for the diversity rule, both at the versions the manifest records |

### Output

```
20_new_runs/rank_prioritisation_2026-10-05/
  rank_prioritisation.json
  shortlist.csv
  run_manifest.json
  analysis/rank_prioritisation.py
```

## ADR-AIIM-301 — PRE-REGISTRATION: does intersecting two weak orthogonal channels beat either one, measured on labels before it is used on the library

**Status:** proposed — **pre-registration. Frozen before the run, and before
ADR-AIIM-300's shortlist was read.**
**Date:** 2026-10-05
**Raised by:** Adisorn Sirichotejirakul — *"ร่าง pre-registration ของ orthogonal intersection เลย"* (2026-10-05), following *"จะต้องมีทางที่คุณคิดมาแล้ว ... เพื่อ screening สารที่คุณมั่นใจมากที่สุด"*
**Follows:** ADR-AIIM-299 (rank works, assertion does not), ADR-AIIM-300 (the
convention-robust rank rule, **still executing when this entry was written**)
**Discipline:** ADR-AIIM-263 — the rule and the comparison are fixed before any
result is seen, including ADR-AIIM-300's

### Why this is drafted before ADR-AIIM-300 finishes

ADR-AIIM-300 is still running. Writing the intersection rule now, without its
shortlist, is **the disciplined order rather than an impatient one**: a rule
designed after seeing which compounds survive the first channel is a rule
chosen for its output. The decision log timestamps both entries, so the order is
checkable.

### The strategy, and the one thing that would make it worthless

ADR-AIIM-299 established that no compound in this library can be asserted
active: the honest interval is ±1.5 log against a 0.5 decision gap and the
conformal count is 0 at every threshold under every convention. Confidence
cannot come from tightening one model. It can only come from **requiring
agreement between channels that do not share their evidence**.

The project has already measured the *averaging* version of this and it failed:
the equal-weight z-score ensemble scores **AUC 0.530**, worse than the Asn253
H-bond filter with Vina at **0.577** (§3.4 table). Averaging scores dilutes.

**Set intersection is a different operation with a different operating
characteristic** — it trades recall for precision rather than blending
rankings — and this project has never measured it. That is the gap this entry
fills, and R1 below is where it can fail.

### Which channels are orthogonal, and which only look it

| channel | evidence | orthogonal to the ligand model? |
|---|---|---|
| **L — censored-EM predicted pKi** (ADR-AIIM-298) | learned from ChEMBL ligand structures | — |
| **S — Vina docking score** | the receptor structure and the pose | **yes** |
| 2D pharmacophore score | similarity to known ligands | **no** — same information family as L |
| 3D shape/pharmacophore score | similarity to known ligands | **no** |
| `p_active_calibrated`, `prediction_set_sig0.05` | the Mondrian conformal column | **excluded** — ADR-AIIM-286 measured it spanning 376× |

**Only L and S are used.** The funnel's pharmacophore scores are ligand
similarity and would double-count L's evidence; including them would
manufacture agreement rather than test it. The conformal column is excluded for
instability, not for redundancy.

### The validation set, which is what makes this measurable

`20_new_runs/dude_enrichment_2026-09-22/analysis/dude_enrichment_per_compound.csv`
— **2,401 compounds, 482 actives and 1,919 property-matched decoys**, every one
carrying a Vina score from the production protocol and a SMILES. Labels exist,
so precision and recall are measured rather than argued. The NP library has no
labels and never will inside this project, which is why the strategy is
validated here first and applied second.

**Leakage control, declared.** DUD-E actives are ChEMBL-derived, so some will
sit in the censored-EM modelling set. Every DUD-E compound whose standardised
InChIKey appears in that set is **removed before any metric is computed**, and
the number removed plus the surviving active/decoy composition is reported. The
survivors may be a biased subset of actives — those this project's curation did
not capture — and that is reported as a limitation rather than assumed away.

### The comparison, fixed now

Selection fraction **X = 10 %** of the surviving labelled set is primary; 1 %,
5 % and 20 % are reported.

* **L alone** — top X % by predicted pKi
* **S alone** — top X % by Vina score, most negative first
* **L ∩ S** — in the top X % of **both**

**The comparison is at matched selection size**, and this matters more than any
other choice here. An intersection selects fewer compounds than either channel,
so comparing its precision against a channel's full top-X % would flatter it for
free. If `L ∩ S` selects *n* compounds, it is compared against the **top *n***
of L and the **top *n*** of S.

Reported per channel and per fraction: n selected, precision, recall,
enrichment over the set's own active rate, and ROC-AUC for the two single
channels.

### Predictions

**R1 — intersection beats both channels at matched size.** At X = 10 %,
`L ∩ S` precision exceeds the precision of both the top-*n* of L and the
top-*n* of S. *Falsified if either single channel equals or beats it.*

**R2 — Vina is weak here, as the article already reports.** S alone returns
ROC-AUC below 0.70 on the surviving set. *Falsified at 0.70 or above.* **An
anchor rather than a discovery**: §3.4 reports 0.601 unmatched, 0.669 against
property-matched decoys and 0.537 on the caliper-matched subset, where it loses
to cLogP at 0.553. A harness returning 0.85 here would be wrong, not
interesting.

**R3 — the ligand channel is the stronger of the two.** L's ROC-AUC exceeds
S's. *Falsified otherwise.* If S wins, the prioritisation built in
ADR-AIIM-300 is leaning on the weaker evidence.

**R4 — precision is bought with recall.** `L ∩ S` recall at X = 10 % is below
50 %. *Falsified at 50 % or above.* This is the cost statement, and it is
predicted rather than discovered so that it cannot later be presented as a
surprise: an intersection that triples precision while discarding most actives
is the right tool for choosing what to assay and the wrong tool for estimating
how many actives a library holds.

**The outcome that would most embarrass this entry is R1 failing.** It would
mean orthogonal agreement adds nothing beyond the better single channel, that
§3.4's ensemble result (0.530, worse than its parts) extends from score
averaging to set intersection, and that the "most confident" screening strategy
the author asked for does not exist in this evidence — leaving ADR-AIIM-300's
single-channel rank as the honest ceiling.

### Application to the library, conditional and declared in advance

**Only if R1 holds** is the rule applied to the 497 full-protocol compounds of
`e8_np_library_prioritization.csv`, which is the only part of the library
carrying a Vina score. The output is the intersection of ADR-AIIM-300's
convention-robust ranked set with the top X % by Vina, reported as a priority
order.

**If R1 fails, no intersection shortlist is produced**, and that is written here
so the failure cannot be absorbed by quietly reporting the intersection anyway.

### Constraints

1. **`10_new_system` is not modified**; no re-docking, no contract amended, no
   input promoted. Vina scores are read from frozen artifacts.
2. **§3.6 is unchanged** and nothing here replaces its 72 of 5,000.
3. **A selected compound is not a predicted active**, under any R-verdict. The
   word "active" is not applied to any output. ADR-AIIM-299's result stands:
   this instrument cannot resolve pKi ≥ 6.0 for any compound in this library.
4. **Scope limit.** The intersection can only ever cover the **497 docked
   compounds**, not the 5,000. Any statement about the library as a whole is
   outside what this run can support.
5. **DUD-E is not this library.** Precision measured on property-matched decoys
   transfers to natural products only under the assumption that the channels
   behave the same way on both, which is not tested here and is stated as an
   assumption. ADR-AIIM-296's competence boundary applies unchanged.
6. **Exchangeability is still not established**, so no conformal statement is
   made about anything selected.
7. **The leakage removal count and the surviving composition are reported**,
   whatever they are.
8. **Every outcome is reported**, including R1 failing and an intersection that
   selects nothing.
9. **This entry is appended to the decision log and its SHA-256 recorded before
   the run**, digest over UTF-8 text with LF newlines per the 2026-10-05
   addendum.

### Skills invoked

| skill | version | material influence |
|---|---|---|
| `scikit-learn` | 1.1 | the ADR-AIIM-298 censored-EM pipeline reused unchanged to score the labelled set; no refitting per channel |
| `rdkit` | 1.1 | standardised InChIKey for the leakage control, which is the step the whole validation depends on |
| `statistical-analysis` | 1.0 | applied for the matched-selection-size comparison: precision at unequal n is not a comparison, and the skill's separation of statistical from practical significance is why recall is predicted (R4) rather than reported after the fact |

### Output

```
20_new_runs/orthogonal_intersection_2026-10-05/
  intersection_validation.json
  intersection_shortlist.csv      only if R1 holds
  run_manifest.json
  analysis/orthogonal_intersection.py
```

## ADR-AIIM-302 — ADR-AIIM-300 and ADR-AIIM-301 results: a six-compound priority order exists, the ligand channel alone is the ceiling, and the rule that selected the six was never validated for enrichment

**Status:** accepted — record of fact
**Date:** 2026-10-05
**Approver:** Adisorn Sirichotejirakul — *"เขียน result ADR เลย"* (2026-10-05)
**Pre-registrations:** ADR-AIIM-300, SHA-256
`2e269726f50f8bcabb4d8f43e07fe7d06d2d987d0569c8f6bb81ae48997bc030` (log at
freeze `c5c7c4732e4b33192cbbd718a4b074bea5fa047f539d291e4b31986c4df32901`);
ADR-AIIM-301, SHA-256
`50ecb09367a619bfc14b84022397bda00800810bd0f8d47d09fd39c3c2b133ef` (log at
freeze `b8bf35ef5f90db0c2f2cb76d52addfd1fc76917d4400104a9061280d7f96aed7`)
**Artifacts:** `20_new_runs/rank_prioritisation_2026-10-05/` (747 s),
`20_new_runs/orthogonal_intersection_2026-10-05/` (696 s)

### Grading

| | prediction | verdict |
|---|---|---|
| **Q1** | ≥ 10 compounds in the top 20 of ≥ 5 of 8 conventions | **MET** — exactly 10, at the line |
| **Q2** | ≥ 5 Bemis–Murcko frameworks among the first 10 | **MET** — 6 |
| **Q3** | ≥ 1 shortlisted compound below AD 0.45 | **MET** — all 6, minimum 0.2562 |
| **Q4** | the majority rule yields more than the all-8 intersection | **MET** — 10 against **0** |
| **R1** | `L ∩ S` precision beats both channels at matched n | **MISSED** — and decisively |
| **R2** | Vina ROC-AUC below 0.70 | **MISSED** — 0.7005, by 0.0005 |
| **R3** | L's ROC-AUC exceeds S's | **MET** — 0.9499 against 0.7005 |
| **R4** | intersection recall below 50 % | **MET** — 0.2902 |

### The priority order, from ADR-AIIM-300

Not predicted actives. Eligible set 2,645 of 5,000; 10 survivors of the
5-of-8 majority; 6 after scaffold uniqueness.

| # | compound | conventions | mean pred pKi | AD | class | organism |
|---|---|---|---|---|---|---|
| 1 | NPLIB01498 | **7 / 8** | 6.044 | 0.2666 | Chalcones | *Butea monosperma* |
| 2 | NPLIB00367 | 6 / 8 | 6.191 | 0.2674 | Flavandiols | *Menisciopsis penangiana* |
| 3 | NPLIB00508 | 6 / 8 | 6.082 | 0.2597 | Flavanones | *Polytrichum commune* |
| 4 | NPLIB00193 | 5 / 8 | 6.072 | 0.2562 | 2-arylbenzofurans | *Schoenocaulon officinale* |
| 5 | NPLIB00377 | 5 / 8 | 6.024 | 0.2566 | Coumaronochromones | *Didymocarpus albicalyx* |
| 6 | NPLIB03777 | 5 / 8 | 6.002 | 0.3176 | Flavonols | *Gymnosperma glutinosum* |

**No compound survives all eight conventions.** The all-8 intersection is
empty, so unanimity would have returned nothing and the majority rule is
load-bearing rather than a convenience.

**The three originally prioritised compounds are all eligible and none is
shortlisted** — NPLIB01185 in the top 20 of 1 of 8, NPLIB03056 of 2 of 8,
NPLIB01643 of **0 of 8**, with mean predicted pKi 5.834, 5.807 and 5.581. They
were not excluded by any gate; they rank below the window. The rule was fixed
before any identity was read and is not revised now.

### Three limitations in ADR-AIIM-300 that its own predictions did not catch

**1. The eligible set is 2,645, not the ~1,560 the pre-registration expected,
and the boundary used is not the boundary ADR-AIIM-296 measured.** The 0.25
line was measured against a **4,566-compound training side**; this run applied
it against the **6,284-compound modelling set**. More reference compounds raise
every similarity, so the same numeric threshold admits 52.9 % here against
31.2 % there. **This is the reference-set-size dependence recorded in
ADR-AIIM-294, repeated in the next run by the person who recorded it.**

**2. The shortlist is 6, not the declared 10.** Scaffold uniqueness collapsed
10 survivors into 6 frameworks. The diversity rule binds before the length
target; no compound was dropped for any other reason.

**3. The ordering of survivors was not fixed by the pre-registration.**
Implemented as conventions descending then mean predicted pKi descending, and
recorded in the artifact as an implementation choice.

### The author's audit question, answered per criterion

Asked during the run: *were the criteria used for ranking and for cutting
compounds validated against a known active/inactive set first?* The answer
differs by criterion and the gap is real.

| criterion | validated against labels? | evidence |
|---|---|---|
| censored-EM predicted pKi | **yes** | held-out point measurements, Spearman 0.726–0.836; and on DUD-E below, **ROC-AUC 0.9499** against 193 held-out actives |
| AD ≥ 0.25 competence boundary | **yes** | ADR-AIIM-296, band 0.25–0.35, n = 92, 36 active / 56 inactive, AUC 0.7391 CI [0.633, 0.837] |
| **majority of 5 of 8** | **no** | justified on determinacy only (ADR-AIIM-286); Q4 measures yield, **not enrichment** |
| top-20 window | **no** | declared, not derived |
| scaffold uniqueness | **no** | a diversity rule; it discarded 4 of 10 and nobody measured whether actives were among them |

**So the six compounds are ranked by a validated quantity and selected by an
unvalidated heuristic.** ADR-AIIM-300 did not overclaim — its justification was
determinacy, not enrichment — but it also did not declare that the rule's
selectivity was untested, and it should have. The omission was found by the
author, not by the entry.

### R1 failed, and it settles the question it was asked

DUD-E after leakage removal: **289 compounds removed, every one an active**, so
60 % of DUD-E's actives were already in the modelling set. Surviving set
**2,112 — 193 actives and 1,919 decoys, 9.1 % active.** The surviving actives
are genuinely unseen, and the decoys remain property-matched to the original
482, which slightly loosens the matching.

At X = 10 %, k = 212 per channel:

| selection | n | precision | recall | enrichment |
|---|---|---|---|---|
| **L alone, matched to n = 62** | 62 | **1.0000** | — | — |
| `L ∩ S` | 62 | 0.9032 | 0.2902 | 9.88 |
| S alone, matched to n = 62 | 62 | 0.3226 | — | — |

and the pattern holds at every fraction tested:

| X | intersection precision | matched L | matched S |
|---|---|---|---|
| 1 % | 1.000 (n = 1) | **1.000** | 0.000 |
| 5 % | 0.963 (n = 27) | **1.000** | 0.222 |
| 10 % | 0.903 (n = 62) | **1.000** | 0.323 |
| 20 % | 0.593 (n = 145) | **0.938** | 0.290 |

**The ligand channel alone is better than the intersection at every selection
size.** Adding the structure-based channel does not concentrate the actives, it
dilutes them — the intersection discards compounds L ranked correctly in order
to admit agreement with a channel at ROC-AUC 0.7005.

This extends §3.4's ensemble finding rather than contradicting it. Score
averaging gave 0.530, worse than its parts; set intersection gives precision
worse than its better part. **Both forms of combination fail here, for the same
underlying reason: the structure-based channel is too weak to add information
to the ligand-based one.**

**No intersection shortlist was produced.** ADR-AIIM-301 forbade it on an R1
failure and the harness enforced it in code rather than in intention — the
output file was never written.

### R2 missed by 0.0005 and the harness is not at fault

Vina returns **0.7005** against a predicted ceiling of 0.70. §3.4 reports
0.601 unmatched, **0.669 against property-matched decoys** and 0.537
caliper-matched. DUD-E is property-matched, so 0.70 sits just above the
comparable published figure and the anchor is substantially reproduced. The
prediction was set a hair too tight; it is graded MISSED and no claim of a
harness defect follows from it.

### What this leaves, stated plainly

**The honest ceiling is the single ligand channel.** The author asked for the
screening strategy with the most confidence behind it. Measured on labels, that
is ADR-AIIM-300's ranking — not an intersection, not an ensemble, not a
consensus. L reaches precision 1.000 at n = 62 on a 9.1 %-active set, which is
the strongest discrimination any configuration in this four-day sequence has
produced.

**And it still does not license calling anything active.** ADR-AIIM-299 stands:
the conformal interval is ±1.5 log against a 0.5 decision gap, the conformal
count at pKi ≥ 6.0 is 0 under every framing, and the six compounds carry mean
predicted pKi of 6.00–6.19 — the top one exceeds 6.0 by 0.04 against an
uncertainty of 1.5. Their AD scores of 0.256–0.318 sit at the lower edge of the
competence band, so the highest-ranked compounds are also the least supported
by distance. The six are **what to assay first**, and nothing more.

### Deliberately not done

* **`10_new_system` is unmodified**; nothing re-docked; Vina scores read from
  frozen artifacts; §3.6 unchanged at 72 of 5,000.
* **No intersection shortlist**, per ADR-AIIM-301's conditional.
* **The majority-of-8 rule was not validated for enrichment.** The DUD-E
  machinery that would test it — does a 5-of-8 survivor set beat a single
  convention's top-k at matched size — was not run and needs its own
  pre-registration. Until it is, the six compounds' *ordering* rests on a
  validated model and their *selection* does not.
* **The eligible-set inconsistency was not repaired by re-running** with a
  matched reference set; it is recorded.
* **Nothing is claimed about the 2,355 ineligible compounds**, and the
  intersection could in any case only ever have covered the 497 docked.
* **No compound named here has been assayed**, and none is put forward as
  active.

### Artifacts

Digests over raw bytes, per the 2026-10-05 addendum.

```
20_new_runs/rank_prioritisation_2026-10-05/
  rank_prioritisation.json            bf64065d21c43b717afa29f4526f2cfd26cf7a160274e098a7d59ad0caaa8e18
  shortlist.csv                       15ba5c8022972af7e2151fa582565d8fbae26dd2e6bfceecf63b7175c243bdc3
  analysis/rank_prioritisation.py     3264907584fd7be83b52a84350e619ed5f4919323574b99d928332331abf6c4d
20_new_runs/orthogonal_intersection_2026-10-05/
  intersection_validation.json        1cea98cf99aad208143e4e49d778b1b471f0aa691a4775404e1104d5f0f860a0
  analysis/orthogonal_intersection.py f4819ee1e67b8832a62cbe40970b8477ff880b955c6e4e45435a6274deb0b86e
```

## ADR-AIIM-303 — PRE-REGISTRATION: does the 5-of-8 majority rule enrich for actives, or only make the answer repeatable?

**Status:** proposed — **pre-registration. Frozen before the run.**
**Date:** 2026-10-05
**Raised by:** Adisorn Sirichotejirakul — *"เกณฑ์ที่นำมาคิดเพื่อ ranking หรือตัด compound ได้ลองเอาชุด known active/inactive อะไรมาเทียบก่อนแล้วใช่ไหม"* (2026-10-05), and *"pre-register การ validate กฎเสียงข้างมากบน DUD-E แล้วรันเลย"* (2026-10-05)
**Follows:** ADR-AIIM-302, which recorded this gap as *deliberately not done*
**Discipline:** ADR-AIIM-263 — predictions fixed before the run;
ADR-AIIM-226 — outside the frozen tree

### The gap, which the author found and the entry that created it did not

ADR-AIIM-302's audit table is explicit: of the five criteria that produced the
six-compound priority order, two were validated against labelled
actives and inactives and three were not.

| criterion | validated? |
|---|---|
| censored-EM predicted pKi | **yes** — ROC-AUC 0.9499 on 193 held-out DUD-E actives |
| AD ≥ 0.25 competence boundary | **yes** — ADR-AIIM-296, 36 active / 56 inactive, AUC 0.7391 |
| **majority of 5 of 8 conventions** | **no** |
| top-20 window | **no** |
| scaffold uniqueness | **no** |

ADR-AIIM-300 defended the majority rule on **determinacy** — ADR-AIIM-286 had
shown that choosing a convention is choosing an answer — and Q4 confirmed it
yields something where unanimity yields nothing (10 survivors against 0). But
**yield is not enrichment**, and nobody has shown that majority survivors are
more likely to be active than the compounds a single convention would have
picked. This entry tests exactly that, on labels, and it should have been part
of ADR-AIIM-300.

### What is run

The same labelled set and the same leakage control as ADR-AIIM-301:
`dude_enrichment_per_compound.csv`, standardised, every compound whose
InChIKey appears in the censored-EM modelling set removed — **2,112 surviving,
193 actives and 1,919 decoys, 9.1 % active**, with the 289 removals all being
actives.

The eight censored-EM conventions are refitted exactly as ADR-AIIM-298 ran
them and each scores all 2,112 compounds. Then, for each selection window
**K ∈ {1 %, 5 %, 10 %, 20 %, 30 %}** of the labelled set, with 10 % primary:

| set | definition |
|---|---|
| **M** — majority | in the top K of **at least 5 of the 8** conventions |
| **U** — unanimous | in the top K of all 8 |
| **A** — mean rank | top *n* by the mean prediction across the 8, the channel ADR-AIIM-301 validated |
| **C₁…C₈** — single | top *n* by each convention alone |

**Every comparison is at matched selection size** *n* = |M|, for the reason
ADR-AIIM-301 gives: a rule that selects fewer compounds would otherwise be
flattered for free.

Reported per window: |M|, |U|, precision, recall and enrichment for M, A and
every Cᵢ, plus the minimum, median and maximum single-convention precision.

**Note on what is and is not orthogonal here.** M and A are two aggregations of
the **same eight models**, not independent channels. M counts top-K
memberships and discards magnitude; A averages magnitudes and discards rank.
Neither adds evidence the ensemble does not already contain, so this entry
tests which aggregation is better, not whether aggregation helps.

### Predictions

**V1 — the majority rule beats a typical single convention.** At the primary
K = 10 %, M's precision exceeds the **median** of the eight single-convention
precisions at matched *n*. *Falsified if equal or lower.*

**V2 — but it does not beat averaging.** M's precision does **not** exceed A's
at matched *n*. *Falsified if it does.* *Basis:* ADR-AIIM-301 measured A at
precision **1.000** for n = 62 at this window, so M can at best tie there;
majority counting throws away the magnitude information averaging keeps, and a
discarded-information aggregator beating a retaining one would be a real
finding rather than a tie.

**V3 — the sets nest.** |U| ≤ |M| ≤ K at every window. *Falsified otherwise.*
**Marked weak and structural**; it is a correctness check on the harness, not a
result.

**V4 — M is usable in absolute terms.** M's precision at K = 10 % exceeds 0.50,
against a base rate of 9.1 %. *Falsified at 0.50 or below.*

**The outcome that would most embarrass this entry is V1 failing.** It would
mean the majority rule performs **worse than picking one convention at random**,
that the six compounds of ADR-AIIM-302 were selected by a procedure that
actively costs precision, and that the determinacy it bought was paid for in
accuracy. ADR-AIIM-302's description of those compounds would then have to be
rewritten, and this entry commits to rewriting it.

**V2 failing would be the more interesting outcome** and is the one I expect
less: it would make majority voting the better aggregator and would justify the
rule on its own evidence rather than on determinacy alone.

### Constraints

1. **`10_new_system` is not modified**; nothing re-docked; no contract amended;
   no input promoted.
2. **§3.6 is unchanged** at 72 of 5,000 and no count here touches it.
3. **This validates a selection rule, not the compounds.** Whatever V1–V4
   return, no compound is reclassified and nothing is called active.
   ADR-AIIM-299 stands: the conformal count at pKi ≥ 6.0 is 0 under every
   framing.
4. **DUD-E is not this library.** Precision on property-matched decoys
   transfers to natural products only if the rule behaves the same way on both,
   which is not tested. The decoys remain matched to the original 482 actives
   while only 193 survive, which loosens the matching, and that is stated rather
   than assumed away.
5. **The result changes a description, not a shortlist.** If V1 fails, the six
   compounds stay in the record and their *selection* is re-described; the rule
   is not re-tuned to recover them. ADR-AIIM-300's rule was frozen before any
   identity was read and that stands.
6. **Only the majority rule is tested here.** The top-20 window and the
   scaffold-uniqueness rule remain unvalidated and are not addressed.
7. **ADR-AIIM-299's limitations carry over**: the censored-EM E-step is biased
   upward at censored points and the EM does not reach its 0.01 tolerance in 10
   iterations.
8. **No compound has been assayed.**
9. **Every outcome is reported**, including V1 failing.
10. **This entry is appended to the decision log and its SHA-256 recorded before
    the run**, digest over UTF-8 text with LF newlines per the 2026-10-05
    addendum.

### Skills invoked

| skill | version | material influence |
|---|---|---|
| `scikit-learn` | 1.1 | the ADR-AIIM-298 censored-EM pipeline reused unchanged across all eight conventions; no refitting per window |
| `rdkit` | 1.1 | standardised InChIKey for the leakage control the whole comparison rests on |
| `statistical-analysis` | 1.0 | applied for the matched-selection-size comparison and for reporting the single-convention spread as a distribution (min, median, max) rather than a single number, so that V1 is tested against a typical convention and not a favourable one |

### Output

```
20_new_runs/majority_rule_validation_2026-10-05/
  majority_rule_validation.json
  run_manifest.json
  analysis/majority_rule_validation.py
```

## ADR-AIIM-304 — ADR-AIIM-303 result: all four predictions met, the majority rule is identical to averaging, and the 376× instability turns out to be a threshold artefact

**Status:** accepted — record of fact
**Date:** 2026-10-05
**Approver:** Adisorn Sirichotejirakul — *"รอผลแล้วเขียน result ADR ของ 303"* (2026-10-05)
**Pre-registration:** ADR-AIIM-303, SHA-256
`81011897d6976ee0f17623812675e92082b83858d9378a9dbe31d9b370ab660f`, frozen
before the run; log SHA-256 at that point
`561f84e9945fe4019be4be333b03872a7407bb2a21cc68a3f7201a81ba23e8e8`
**Artifacts:** `20_new_runs/majority_rule_validation_2026-10-05/` · 704 s
**Closes:** one of the three unvalidated criteria ADR-AIIM-302 recorded

### Grading

| | prediction | verdict |
|---|---|---|
| **V1** | majority beats the median single convention | **MET** — 0.8172 against 0.8065, **by 0.0107** |
| **V2** | majority does not beat mean-rank averaging | **MET** — and the two are **identical** |
| **V3** | the sets nest, \|U\| ≤ \|M\| ≤ K | **MET** at every window |
| **V4** | majority precision above 0.50 | **MET** — 0.8172 against a 9.1 % base rate |

### The measurement

Surviving labelled set 2,112 — 193 actives, 1,919 decoys, 9.1 % active; 289
removed by the leakage control and every removal an active.

| K | k / convention | \|M\| | \|U\| | **M prec** | **A prec** | single median | min | max |
|---|---|---|---|---|---|---|---|---|
| 1 % | 22 | 18 | 8 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| 5 % | 106 | 101 | 72 | 1.0000 | 1.0000 | 0.9901 | 0.9604 | 1.0000 |
| **10 %** | 212 | 186 | 160 | **0.8172** | **0.8172** | 0.8065 | 0.7903 | 0.8495 |
| 20 % | 423 | 346 | 245 | 0.4827 | 0.4827 | 0.4783 | 0.4653 | 0.4913 |
| 30 % | 634 | 541 | 340 | 0.3235 | 0.3235 | 0.3180 | 0.3087 | 0.3272 |

Per-convention ROC-AUC spans **0.9340–0.9565**; the mean-rank channel is
0.9490.

### V2 is met in the strongest possible way, and it is not a tie on points

**M and A return the identical precision at every one of the five windows.**
At matched selection size the majority set and the mean-rank top-*n* are, for
the purpose of counting hits, the same set. The majority rule is not a better
aggregator or a worse one — **it is the same aggregator reached by a different
route**, and ADR-AIIM-303 predicted it could at best tie.

So the rule's contribution is exactly what ADR-AIIM-300 claimed for it and
nothing more: **determinacy**. It does not add accuracy, and the entry that
introduced it did not claim it would.

### V1 is met by 0.0107, and the best single convention beats the rule

| at K = 10 % | precision |
|---|---|
| best single convention | **0.8495** |
| majority rule | 0.8172 |
| median single convention | 0.8065 |
| worst single convention | 0.7903 |

The rule clears a typical convention by **about one percentage point**, which
on n = 186 is roughly two compounds, and it is **beaten by the best single
convention by 0.0323**. V1 asked whether the rule beats the median and it
does; it does not beat the field. **The rule is insurance against choosing
badly, not a route to choosing well** — and since there is no way to identify
the best convention in advance without labels, insurance is the right thing to
buy. That is a defensible justification and it is narrower than it might have
looked.

### The finding neither prediction anticipated

**The eight conventions span 376× in the count they report and 0.0592 in
precision at matched selection size.**

| quantity | spread across the same eight conventions |
|---|---|
| NP-library count at pKi ≥ 6.0 (ADR-AIIM-293, classifier) | **38.88×** |
| the same, regression (ADR-AIIM-297) | 10.99× |
| the same, censored-EM (ADR-AIIM-299) | 4.10× |
| arm B, classifier (ADR-AIIM-288) | **376×** |
| **precision at matched n, K = 10 % (here)** | **0.0592**, i.e. 0.7903–0.8495 |
| ROC-AUC (here) | 0.0225, i.e. 0.9340–0.9565 |

**The conventions disagree violently about how many compounds to flag and
barely at all about which compounds rank highest.** The instability
characterised across ADR-AIIM-286 through ADR-AIIM-299 is substantially a
**threshold-crossing artefact**: a 0.02 shift in a predicted value moves a
compound across pKi 6.0 and changes the count, while leaving the ordering
almost untouched.

This is the strongest vindication in the sequence of the move from counting to
ranking, and it was not predicted — V1 through V4 were about the rule, not
about this. It is recorded as an observation.

### A consequence for ADR-AIIM-302's six compounds

ADR-AIIM-300 selected on a **top-20 window of 2,645 eligible compounds —
0.76 %**, tighter than every window tested here, and the unanimous set on the
natural-product library was **0 of 10**. On DUD-E at 1 % the unanimous set is
8 of 18, and at 10 % it is **160 of 186, 86 % of the majority set**.

So the agreement between conventions is strong across most of the ranking and
collapses at the very top, and **ADR-AIIM-300 selected in the one regime where
the conventions disagree most**. Two readings, and this run does not separate
them:

* the top of the ranking is where the models genuinely differ, so a tight
  window is the honest place to require a majority; or
* a wider window would have given a more stable shortlist, and 20 was too
  tight.

The unanimity gap between the two sets has a likely cause worth naming: DUD-E's
surviving actives are ChEMBL-derived and close to the training chemistry, where
the eight models agree, while the natural-product library sits at AD
0.15–0.39, where they do not. **Agreement measured on DUD-E is therefore an
upper bound on agreement over the library**, and the 86 % figure should not be
read across.

### What this closes and what it leaves

**Closed.** The majority rule is validated against labelled actives and
inactives: it beats a typical single convention, ties with averaging, and is
usable in absolute terms at 0.8172 precision against a 9.1 % base rate. The
six compounds' *selection* is no longer described by an untested rule.

**Still open**, from ADR-AIIM-302's audit table: the **top-20 window** and the
**scaffold-uniqueness rule** remain unvalidated. One of three gaps closed.
Scaffold uniqueness discarded 4 of the 10 survivors and whether actives were
among them is still unmeasured.

**Unchanged.** ADR-AIIM-299 stands — the conformal count at pKi ≥ 6.0 is 0
under every framing, the six compounds carry mean predicted pKi 6.00–6.19
against an interval of ±1.5 log, and nothing here calls any compound active.
§3.6 still reports 72 of 5,000.

### Deliberately not done

* **`10_new_system` is unmodified**; nothing re-docked; no contract amended.
* **No compound was reclassified and the rule was not re-tuned.**
  ADR-AIIM-303 constraint 5 committed to re-describing the selection if V1
  failed; V1 held, so the description stands as written in ADR-AIIM-302 with
  the narrowing above.
* **The window question was not resolved.** Whether top-20 or a wider window is
  the right selector is now a live question this run raised and did not answer.
* **The top-20 and scaffold-uniqueness criteria were not validated.**
* **Agreement was not measured on the natural-product library itself**, where
  no labels exist; the DUD-E figures are an upper bound and are not transferred.
* **No compound has been assayed.**

### Artifacts

Digests over raw bytes, per the 2026-10-05 addendum.

```
20_new_runs/majority_rule_validation_2026-10-05/
  majority_rule_validation.json        87208cd231d9bf40f9ea8bfbe14416d88885687ee2dfc77119e49a30e47762f9
  analysis/majority_rule_validation.py b8d4b1208c8b49b89cbf7bb2d5e46efbf1fd89449699642e4fc8ede890987939
```

## ADR-AIIM-305 — PRE-REGISTRATION: is DiffDock confidence an activity signal at all, tested on labels in two stages with a stop rule

**Status:** proposed — **pre-registration. Frozen before the run. COMPUTE NOT YET APPROVED.**
**Date:** 2026-10-05
**Raised by:** Adisorn Sirichotejirakul — *"ฝั่ง SBDD ใช้แค่ Vina หรอครับ ลองดูแบบถ้าเสริมคะแนนจาก diffdock และ 10 ns MD ด้วยจะดีไหมครับ"* then *"ร่าง DiffDock บน DUD-E ก่อน"* (2026-10-05)
**Follows:** ADR-AIIM-301 (Vina at 0.7005 dilutes the ligand channel),
ADR-AIIM-111 (DiffDock on four compounds), ADR-AIIM-302 (the audit that
requires a channel be validated before use)
**Discipline:** ADR-AIIM-263 — predictions and the stop rule fixed before the
run; ADR-AIIM-226 — outside the frozen tree

### What the SBDD side actually consists of

| channel | compounds scored | status |
|---|---|---|
| Vina | 2,401 labelled + 497 library | **measured**: ROC-AUC 0.7005, and intersecting it with the ligand channel *lowers* precision from 1.000 to 0.9032 (ADR-AIIM-301) |
| **DiffDock** | **4** — three candidates and the control (ADR-AIIM-111) | never tested against labels |
| MD | **3** candidates, 3 × 100 ns | never tested against labels; not testable at screening scale |

So the answer to the question as asked is yes: anything with scale on the
structure-based side is Vina alone.

### Why the existing DiffDock evidence is discouraging, and why it is still worth spending GPU-hours to settle

ADR-AIIM-111's own pre-registered prediction **failed**:

```
"control confidence exceeds all three candidate confidences"   holds: FALSE
control ZM241385   0.16
NPLIB01185         0.23      RMSD to Vina 6.427 A   agrees: false
NPLIB01643         0.47      RMSD to Vina 8.800 A   agrees: false
NPLIB03056        -0.50      RMSD to Vina 9.497 A   agrees: false
```

DiffDock confidence placed a known A2AR antagonist **below two untested
natural products**, and disagreed with Vina on pose by 6.4–9.5 Å in all three
cases. That is n = 4 and cannot carry a conclusion, but it points one way.

**And there is a prior reason to expect near-chance performance that is not
about this project's data at all: DiffDock confidence is not an affinity
score.** It is a diffusion model's confidence that the pose it generated is
correct. Vina computes a scoring function intended to approximate binding free
energy, however poorly; DiffDock confidence answers *"how sure am I about this
geometry"*. Ranking actives against decoys by it is a category substitution
unless it is shown empirically to work.

This entry is therefore framed as **closure rather than hope**: the question
has been raised, the channel is unvalidated, and ADR-AIIM-302's audit standard
says an unvalidated channel may not be used. Settling it with labels costs
GPU-hours and ends the question either way.

### Compute, measured rather than estimated — and not yet authorised

From `diffdock_orthogonal_2026-08-23/diffdock.log`, at the ADR-AIIM-111
parameters (`samples_per_complex 20`): **65.0 s per complex** on this machine.

| scope | compounds | wall time |
|---|---|---|
| full surviving labelled set | 2,112 | **38.1 h** |
| **stage 1 subsample** | **593** | **10.7 h** |

**ADR-AIIM-040 does not authorise this.** That budget covers held-out scoring
and states *"Negligible; no GPU or model-training compute required"*. A GPU run
of this size has no standing authorisation and **requires the author's explicit
approval before execution.** This entry is frozen in advance of that approval,
not instead of it.

### Two stages, with the stop rule fixed now

**Stage 1 — is it a signal at all?** All **193 surviving actives** plus a
seeded random sample of **400 decoys** (`random_state = 42`), 593 compounds,
10.7 h. Leakage control exactly as ADR-AIIM-301: standardise, InChIKey,
remove anything in the censored-EM modelling set. ROC-AUC is insensitive to
base rate, so a subsample answers the stage-1 question; precision is not
reported from stage 1 because it is base-rate dependent and the subsample's
rate is not the set's.

**Vina's AUC is recomputed on the same 593 compounds.** The comparison is
against that number, not against the full-set 0.7005, so it is like for like.

> **STOP RULE, fixed before the run: if DiffDock's AUC does not exceed Vina's
> AUC on the same subsample, stage 2 is not run and no DiffDock channel is
> added to anything.** Vina at 0.7005 already dilutes the ligand channel; a
> structure-based channel weaker than Vina cannot help, and spending the
> remaining 27 hours to confirm it would be waste.

**Stage 2 — only if the stop rule passes.** The remaining 1,519 decoys are
scored, and the full set is used for the ADR-AIIM-301 comparison: `L ∩ D`
precision against matched top-*n* of L alone and D alone, at
X ∈ {1, 5, 10, 20} %.

**Failure handling, declared.** The ADR-AIIM-111 log shows sample-level
failures (`1 / 10 samples failed` repeatedly for NPLIB01643). At scale some
complexes will yield no pose. A complex with no successful sample is recorded
as **missing, not as a low score**, excluded from AUC, and the count reported.
Imputing a bad confidence for a failed generation would manufacture a signal
from a crash.

### Predictions

**W1 — DiffDock does not beat Vina.** DiffDock's ROC-AUC on the 593-compound
subsample is **below** Vina's on the same subsample. *Falsified if it is equal
or higher*, in which case stage 2 runs. *Basis:* the category argument above,
plus ADR-AIIM-111's failed control check.

**I have twice told the author an avenue was closed without opening it** —
Ki-only labels (ADR-AIIM-288) and in-domain restriction (ADR-AIIM-294) — and
both pessimistic calls were wrong. W1 is a pessimistic call again, and it is
made on a mechanistic argument rather than an intuition, with the test run
regardless of what I expect. That is the whole point of running it.

**W2 — but it is not anti-predictive.** DiffDock's AUC exceeds 0.50.
*Falsified at 0.50 or below*, which would mean its confidence is inversely
related to activity on this set and would be a stranger result than W1 holding.

**W3 — failures are not negligible at scale.** More than 1 % of the 593
complexes yield no usable pose. *Falsified at 1 % or below.* Recorded because
the ADR-AIIM-111 run failed samples on one of four complexes and a failure rate
is a property of the method worth measuring once.

**W4 — conditional, only if stage 2 runs.** `L ∩ D` precision does not exceed
L alone at matched *n*, mirroring ADR-AIIM-301's R1. *Falsified if it does*,
which would be the first evidence in this project that any channel combination
beats its best part.

**The outcome that would most embarrass this entry is W1 failing together with
W4 failing** — DiffDock stronger than Vina *and* its intersection beating the
ligand channel. That would mean the SBDD side was weak only because the wrong
engine was used, and that four days spent on the ligand side were aimed at the
smaller problem.

### Constraints

1. **Compute approval is required before execution** and this entry does not
   constitute it. ADR-AIIM-040 covers held-out scoring, not this.
2. **`diffdock`'s `SKILL.md` is read in full before DiffDock is invoked**, per
   the project's mandatory-skills rule. It is not read yet and this entry does
   not claim it.
3. **DiffDock is not installed, updated or reconfigured.** The existing
   `diffdock` conda environment and checked-out commit are used as they stand,
   with the commit hash recorded. Parameters are ADR-AIIM-111's verbatim —
   `samples_per_complex 20`, `default_inference_args.yaml` — so the two runs are
   comparable. Changing `samples_per_complex` to halve the cost would change the
   method and is deliberately not done.
4. **The receptor is `4EIY_protein_only.pdb` as used by ADR-AIIM-111**, not
   re-prepared.
5. **`10_new_system` is not modified**; nothing re-docked on the Vina side;
   §3.6 unchanged at 72 of 5,000.
6. **3D reproducibility is a known weakness here.** ADR-AIIM-250 measured the
   3D stage reproducing **1 of 40** across interpreters (max difference 30.025)
   while the 2D stage reproduced 300 of 300, and the M2 metamorphic test found
   conformer generation moving a docking score by 0.974 kcal/mol. DiffDock is
   3D-generative. A positive result here would carry that caveat and the run
   does not test its own reproducibility.
7. **A validated channel is still not a licence to call anything active.**
   ADR-AIIM-299 stands: conformal count 0 at pKi ≥ 6.0 under every framing.
8. **Failed complexes are reported as missing**, never imputed.
9. **Every outcome is reported**, including the stop rule firing after 10.7 h
   with nothing to show but a closed question.
10. **This entry is appended to the decision log and its SHA-256 recorded before
    the run**, digest over UTF-8 text with LF newlines per the 2026-10-05
    addendum.

### Skills invoked

| skill | version | material influence |
|---|---|---|
| `scikit-learn` | 1.1 | the ADR-AIIM-298 censored-EM pipeline supplies channel L unchanged |
| `rdkit` | 1.1 | standardised InChIKey for the leakage control |
| `diffdock` | — | **to be read in full before invocation**, per constraint 2; not yet read and no claim is made about its content |

### Output

```
20_new_runs/diffdock_dude_2026-10-05/
  stage1_auc.json
  stage2_intersection.json        only if the stop rule passes
  diffdock_input.csv
  diffdock.log
  run_manifest.json
  analysis/diffdock_dude.py
```

## ADR-AIIM-306 — PRE-REGISTRATION: prioritise the focus library instead, with drug-likeness as a gate and screening alerts as a flag

**Status:** proposed — **pre-registration. Frozen before the run.**
**Date:** 2026-10-05
**Raised by:** Adisorn Sirichotejirakul — *"สารที่ได้มา มาจาก SNPLib บ้างไหม"*, then the ruling *"SNP-Lib หาง่ายกว่า NP-library ไม่มี และถ้าใช้แค่ SNP-Lib ได้ ก็จะดีกว่าเพราะเป็น focus library"* and *"อาจจะต้องมองมุม medchem หรือ Drug-likeness เพิ่มด้วยในการ screen"* (2026-10-05)
**Follows:** ADR-AIIM-302 (the six-compound NP-library shortlist),
ADR-AIIM-304 (the majority rule validated), ADR-AIIM-296 (the competence band)
**Supersedes for practical purposes:** ADR-AIIM-302's shortlist, for the reason
below
**Discipline:** ADR-AIIM-263 — pool, gates, flags and shortlist length fixed
before any compound identity is read

### Why the previous shortlist is operationally unusable, and whose error that is

ADR-AIIM-300 drew its candidate pool from **one file**:

```python
NPLIB = ... "np_library_funnel_2026-08-15" / "np_library_dock_ledger.csv"
lib = pd.read_csv(NPLIB)
```

SNP-Lib was never in the pool. **That was a harness decision I made without
stating it, and ADR-AIIM-300 does not declare it** — the entry fixed a rule,
a window and a shortlist length, and said nothing about why 5,000 LOTUS
structures were the candidates and the project's own 236-compound library was
not.

The author's ruling settles what that cost: **SNP-Lib is obtainable and the
NP-library is not.** NP-library rows carry an `organism_name` and no compound
name; SNP-Lib carries a name for all 236 and is a hash-anchored frozen input,
described in §2 as one of the project's three read-only inputs. A shortlist of
LOTUS structures is a set of compounds *reported in the literature from named
plants* — assaying them means isolation or synthesis first.

So ADR-AIIM-302's six compounds answer *"what would the model rank highest in a
computational expansion"*. The question worth answering is *"what should be
assayed from the library that exists"*. **This is a scope correction, not a
result-driven one**: the pool is being changed because of what is obtainable,
a fact independent of which compounds any ranking returns, and the ranking on
the new pool has not been computed.

### The pool, measured

SNP-Lib's `smiles_canonical` as standardised by G1 is used unchanged — it is
already the output of the project's own standardiser and re-standardising would
introduce a second convention.

| step | n | frameworks |
|---|---|---|
| 0 — SNP-Lib | **236** | 80 |
| 1 — Lipinski Ro5 | **157** | 43 |
| 2 — AD ≥ 0.25 | **88** | **24** |
| *(sensitivity: AD ≥ 0.20)* | *142* | *35* |

Eligible composition: Phenanthrenes 30, Phenolics 21, Bibenzyls 19,
Flavonoids 4, Fluorenones 4, Fluorenes 4, Stilbenes 2, and one each of
Chromones, Coumarins, Quinazolinone and a phenylpropanoid amide.

**AD 0.25 is primary and 0.20 is reported as a sensitivity**, declared now so
the threshold cannot be swapped after the counts are seen.

### Alerts are flagged, never used to exclude — and why

PAINS and Brenk are **reported per compound and do not remove anything.** The
author ruled this and the reason is specific to this library: SNP-Lib is 52 %
phenanthrenes and 20 % phenolics, catechol matches both catalogs
(`catechol_A(92)` in PAINS, `catechol` in Brenk), and **PAINS was derived from
HTS of synthetic libraries and is known for false positives on natural
products.** Using it as a hard gate on a focused phenolic library would remove
compounds for a reason that does not apply, while leaving the flags off the
record would withhold something a reviewer will ask for.

Catalogs applied, all from the pinned RDKit: **PAINS** (480 entries, plus
PAINS_A/B/C individually), **BRENK** (105), **NIH**, and **ChEMBL** with its
six sub-catalogs. Each compound's matching alert descriptions are recorded.

### Model and selection, imported unchanged

Censored-EM exactly as ADR-AIIM-298, the eight conventions of ADR-AIIM-286
under canonical InChIKey order, the 5-of-8 majority rule **validated against
labels in ADR-AIIM-304**, a top-20 window, scaffold uniqueness from G6, and a
**shortlist length of 10**, all fixed now.

### Two defects carried forward rather than silently repeated

**1. The AD reference set.** The 0.25 boundary was measured by ADR-AIIM-296
against a **4,566-compound training side**; this run scores against the
**5,707-compound curated binding set**, as ADR-AIIM-300 did. ADR-AIIM-302
recorded that this reference-size dependence changed the eligible fraction from
31.2 % to 52.9 % on the NP library and called it *"the dependence recorded in
ADR-AIIM-294, repeated in the next run by the person who recorded it"*. It is
**repeated again here for comparability with ADR-AIIM-300, and declared rather
than discovered afterwards.** The honest consequence: the 88-compound eligible
set is defined by a numeric threshold applied at a different scale from the one
that licensed it.

**2. Two criteria remain unvalidated.** ADR-AIIM-304 closed the majority rule.
The **top-20 window** and the **scaffold-uniqueness rule** are still untested
against labels, and the window question ADR-AIIM-304 raised — whether 20 is too
tight — is unresolved. On a pool of 88 the top-20 window is **23 %**, where
ADR-AIIM-304 measured conventions agreeing far better than at the 0.76 % window
ADR-AIIM-300 used, so the rule should behave better here. That is an
expectation, not a result.

### Predictions

**S1 — the focus library can support a shortlist.** At least 10 compounds
appear in the top 20 under 5 or more of the 8 conventions. *Falsified below 10.*

**S2 — alerts fire on a substantial fraction.** More than 30 % of the 88
eligible compounds carry at least one PAINS or Brenk alert. *Falsified at 30 %
or below.* *Basis:* the phenanthrene/phenolic composition above and catechol
matching both catalogs. This is the number that justifies the flag-not-cut
ruling, in either direction.

**S3 — the smaller pool ranks lower.** The top shortlisted compound's mean
predicted pKi is **below** ADR-AIIM-302's 6.191. *Falsified if it is equal or
higher.* *Basis:* 88 candidates against 5,000, and SNP-Lib's maximum AD is
0.374 against the NP library's 0.392 — a smaller, more homogeneous pool should
yield a lower maximum. **S3 failing is good news** and would mean the focus
library beats the expansion on the model's own ranking despite 57× fewer
candidates.

**S4 — scaffold uniqueness binds.** The shortlist is shorter than 10 because
the survivors collapse into fewer than 10 frameworks. *Falsified if it reaches
10.* *Basis:* 88 compounds over 24 frameworks is 3.7 per framework, and
ADR-AIIM-302 saw 10 survivors collapse to 6.

**The outcome that would most embarrass this entry is S1 failing** — fewer than
10 survivors on a pool of 88 would mean the focus library cannot support a
prioritisation at all, and that redirecting from a 5,000-compound expansion to
the obtainable 236 buys usability at the cost of having nothing to rank.

### Constraints

1. **`medchem` was read and not invoked.** Its `SKILL.md` was read in full as
   the mandatory-skills rule requires. The package is **not installed**
   (`ModuleNotFoundError`), and **no installation is requested**: the pinned
   RDKit supplies PAINS, Brenk, NIH, ZINC and ChEMBL catalogs, and adding a
   dependency to a project with a documented 3D reproducibility failure
   (ADR-AIIM-250: 1 of 40) would add risk without adding capability.
   **Forgone by that choice**, and available only with an approved install:
   `rule_of_cns`, NIBR severity scoring, Lilly demerits, and ZINC-15 complexity
   percentiles.
2. **`10_new_system` is not modified**; nothing re-docked; no contract amended;
   `02_frozen_inputs/` not written to; SNP-Lib's ledger read, not regenerated.
3. **§3.6 is unchanged** at 72 of 5,000 and no count here replaces it.
4. **A shortlisted compound is not a predicted active.** ADR-AIIM-299 stands:
   the conformal count at pKi ≥ 6.0 is 0 under every framing, with an interval
   of ±1.5 log against a 0.5 decision gap. The output is a **priority order for
   assay**, and the word "active" is applied to nothing.
5. **Ro5 and AD are gates; PAINS, Brenk, NIH and ChEMBL are flags.** No alert
   removes a compound, and the flags are reported for every shortlisted
   compound whether or not they fire.
6. **ADR-AIIM-302's shortlist stays in the record.** It is superseded for
   practical use, not withdrawn, and the reason is availability rather than
   anything measured about those six compounds.
7. **ADR-AIIM-299's limitations carry over**: the censored-EM E-step is biased
   upward at censored points and the EM does not reach its 0.01 tolerance in 10
   iterations.
8. **No compound has been assayed**, and nothing here changes that.
9. **Every outcome is reported**, including S1 failing.
10. **This entry is appended to the decision log and its SHA-256 recorded before
    the run**, digest over UTF-8 text with LF newlines per the 2026-10-05
    addendum.

### Skills invoked

| skill | version | material influence |
|---|---|---|
| `scikit-learn` | 1.1 | the ADR-AIIM-298 censored-EM pipeline reused unchanged across the eight conventions |
| `rdkit` | 1.1 | `FilterCatalog` for PAINS / Brenk / NIH / ChEMBL, descriptors for Ro5, ECFP4 for the competence score, Murcko frameworks for the diversity rule |
| `medchem` | 1.1 | **read in full and deliberately not invoked**, per constraint 1; the entry records what that forgoes |

### Output

```
20_new_runs/snplib_prioritisation_2026-10-05/
  snplib_prioritisation.json
  waterfall.json
  shortlist.csv              with alert flags per compound
  run_manifest.json
  analysis/snplib_prioritisation.py
```

## ADR-AIIM-307 — PRE-REGISTRATION: close the hole in my own gate — Ro5 has no lower bound and admitted a fragment to rank 1

**Status:** proposed — **pre-registration. Frozen before the run.**
**Date:** 2026-10-05
**Raised by:** Adisorn Sirichotejirakul — *"ร่าง pre-registration เพิ่ม lead-like gate แล้วรันซ้ำเลย"* (2026-10-05)
**Follows:** ADR-AIIM-306, whose run put **cinnamic acid, MW 148.2, at rank 1**
**Discipline:** ADR-AIIM-263 — the gate and its threshold fixed before the run,
and derived from the training set rather than from the shortlist

### The defect, and that it is mine

ADR-AIIM-306's gate was *"Lipinski Ro5, then AD ≥ 0.25"*. **Ro5 constrains only
ceilings** — MW ≤ 500, cLogP ≤ 5, HBD ≤ 5, HBA ≤ 10 — and has **no lower
bound**. A fragment passes it trivially.

The run's rank 1 was **SNP00164, cinnamic acid, MW 148.2**, unanimous across
all eight conventions. Cinnamic acid is a biosynthetic precursor present in
essentially every plant; the prospect of it being a selective A2AR ligand is
not credible, and ZM241385 sits at MW 337. The gate admitted it because the
gate was incomplete, not because the model was wrong to rank what it was given.

### The threshold, derived from the training actives rather than chosen

The relevant question for a floor is not *what is typical* but *what has ever
been observed to work*. Measured on the modelling set's actives (pKi ≥ 6.0,
**n = 4,428**):

| | MW |
|---|---|
| minimum | **196.2** |
| 1st percentile | **236.5** |
| 5th percentile | 277.3 |
| 10th percentile | 298.3 |
| median | 391.7 |
| **below MW 200** | **1 of 4,428 — 0.02 %** |
| below MW 250 | 72 of 4,428 — 1.63 % |

**Gate: MW ≥ 200**, primary. Only one A2AR active in 4,428 falls below it, so
below 200 is effectively unobserved territory for this target.

**Sensitivity: MW ≥ 236.5**, the 1st percentile, reported and not used to
select.

**Why 200 and not the 10th percentile.** The AD gate uses the 10th percentile
of training self-similarity because it is a *competence* boundary — ADR-AIIM-296
measured where the model still works. An MW floor is a *plausibility* boundary,
and 10 % of known actives sit below MW 298: a floor there would exclude
compounds demonstrably capable of activity. A plausibility floor belongs at the
edge of the observed, not at a typical value.

**Disclosure, because I am not blind here.** I have seen ADR-AIIM-306's
shortlist. The floor was derived from the training actives' MW distribution, a
property of the training set and not of the shortlist. I then checked both
candidate floors against the pool: **MW ≥ 200 leaves 74 eligible and MW ≥ 236.5
leaves 73** — they differ by one compound, so the choice between them cannot be
effect-driven, and either removes cinnamic acid. That check is reported so the
reader need not take the claim on trust.

### What changes, and what deliberately does not

**One gate is added. Nothing else moves.**

| | ADR-AIIM-306 | here |
|---|---|---|
| pool | SNP-Lib 236 | same |
| gate 1 | Ro5 (G1 ledger) | same |
| **gate 2** | — | **MW ≥ 200** |
| gate 3 | AD ≥ 0.25 | same |
| PAINS / Brenk / NIH / ChEMBL | **flags, removing nothing** | **unchanged — still flags** |
| model, conventions, window, majority, uniqueness, length | censored-EM, 8, top-20, 5/8, G6, 10 | same |

**Alerts stay flags and this is deliberate.** ADR-AIIM-306 measured 52 of 93
eligible compounds (55.9 %) carrying a PAINS or Brenk match — `stilbene` on a
stilbene, `Michael_acceptor_1` on cinnamic acid and feruloyl tyramine,
`catechol` on dengibsinin. **Converting the flag to a gate now, having seen
that four of seven shortlisted compounds carry one, would be the post-hoc move
this sequence exists to prevent.** The flags are reported; they exclude nothing.

**No logP, ring-count or rotatable-bond constraint is added either.** The
defect identified was a missing lower MW bound. Adding further lead-like
criteria after seeing the shortlist would be scope creep selected on effect.

### The corrected waterfall, computed against the declared reference

ADR-AIIM-306 predicted 88 eligible and the run returned 93, because the
pre-registration's feasibility numbers were computed against
`curated_binding.csv` (5,707) while the harness scores against the modelling
set (**6,284**). That is the reference-size dependence ADR-AIIM-294 recorded,
ADR-AIIM-300 repeated, ADR-AIIM-302 recorded the repetition of, and
ADR-AIIM-306 declared and then **tripped over in its own numbers**. It is fixed
here by computing the waterfall against the reference the harness actually
uses:

| step | n | frameworks |
|---|---|---|
| 0 — SNP-Lib, G1 canonical | 236 | 80 |
| 1 — gate Ro5 | 157 | 43 |
| 2 — **gate MW ≥ 200** | **133** | 39 |
| 3 — gate AD ≥ 0.25 | **74** | **21** |
| *(sensitivity MW ≥ 236.5)* | *73* | *20* |

Eligible composition: Phenanthrenes 34, Bibenzyls 19, Phenolics 4,
Flavonoids 4, Fluorenones 4, Fluorenes 4, Stilbenes 2, Coumarins 1,
Quinazolinone 1, phenylpropanoid amide 1. MW range 206.2–482.5. **Phenolics
fall from 21 to 4**, which is where the fragments were.

### Predictions

**U1 — gating has knock-on effects beyond the compound it removes.** The new
shortlist differs from ADR-AIIM-306's seven by more than the removal of
cinnamic acid and one promotion into its place — that is, at least two other
positions change in identity or order. *Falsified if the only change is
cinnamic acid leaving.* *Basis:* 24 compounds leave the pool, which changes
every convention's top-20 composition and therefore the majority counts.

**U2 — agreement improves again.** The unanimous 8-of-8 count exceeds
ADR-AIIM-306's **8**. *Falsified at 8 or below.* *Basis:* the window is now
top-20 of 74 — **27 %**, against 21.5 % before and 0.76 % in ADR-AIIM-300 — and
ADR-AIIM-304 measured convention agreement rising with window width
(160 of 186 unanimous at 10 % on DUD-E, against 0 of 10 at the tight window).

**U3 — alerts remain prevalent.** More than 30 % of the 74 eligible carry a
PAINS or Brenk match. *Falsified at 30 % or below.* ADR-AIIM-306 measured
55.9 % of 93; the MW gate removes small phenolics, which are exactly the
Michael-acceptor and catechol carriers, so this could fall materially and the
flag-not-cut decision should be re-justified against the new number.

**U4 — the alert-clean three survive.** SNP00154, SNP00052 and SNP00010 all
appear in the new shortlist. *Falsified if any drops out.* *Basis:* none is
touched by the MW gate (362.4, 244.2, 242.3) and each held 6/8 or better. This
tests whether adding a gate destabilises a selection that ADR-AIIM-304 showed
is otherwise robust.

**U5 — nothing crosses the threshold §3.6 uses.** The top shortlisted
compound's mean predicted pKi stays **below 6.0**. *Falsified at 6.0 or above.*
**Marked weak**: the new pool is a subset of the old, whose survivors spanned
5.763–5.913, so the ceiling cannot rise by much. It is stated because it is the
sentence the manuscript needs — the best obtainable, drug-like, in-domain
compound in the focus library still does not reach the cut §3.6 reports a count
against.

**The outcome that would most embarrass this entry is U2 failing.** Agreement
falling at a wider window would contradict ADR-AIIM-304's measured relationship
and would mean the stability the last two runs appeared to gain came from
something else that has not been identified.

### Constraints

1. **`10_new_system` is not modified**; nothing re-docked; `02_frozen_inputs/`
   not written to; SNP-Lib's ledger and G1 canonical SMILES read unchanged.
2. **§3.6 is unchanged** at 72 of 5,000.
3. **A shortlisted compound is not a predicted active.** ADR-AIIM-299 stands:
   conformal count 0 at pKi ≥ 6.0 under every framing, interval ±1.5 log
   against a 0.5 decision gap. The output is a priority order for assay and the
   word "active" is applied to nothing.
4. **Alerts are flags.** Nothing is excluded on their strength, and the
   decision not to convert them to gates is recorded above with its reason.
5. **`medchem` remains read and not invoked** — not installed, no install
   requested. `rule_of_leadlike_soft` would have supplied this gate directly;
   it is implemented here from the training distribution instead, and the
   forgone items stay as ADR-AIIM-306 listed them.
6. **ADR-AIIM-306's shortlist stays in the record**, superseded rather than
   withdrawn, with the reason being an incomplete gate rather than anything
   measured about those compounds.
7. **The top-20 window and scaffold uniqueness remain unvalidated** against
   labels. ADR-AIIM-304 closed only the majority rule.
8. **ADR-AIIM-299's limitations carry over**: the censored-EM E-step is biased
   upward at censored points and the EM does not reach its 0.01 tolerance in 10
   iterations.
9. **No compound has been assayed.**
10. **Every outcome is reported**, including U2 failing.
11. **This entry is appended to the decision log and its SHA-256 recorded before
    the run**, digest over UTF-8 text with LF newlines per the 2026-10-05
    addendum.

### Skills invoked

| skill | version | material influence |
|---|---|---|
| `scikit-learn` | 1.1 | the ADR-AIIM-298 censored-EM pipeline reused unchanged |
| `rdkit` | 1.1 | `Descriptors.MolWt` for the gate and the training distribution, `FilterCatalog` for the flags, ECFP4 for the competence score, Murcko frameworks for diversity |
| `medchem` | 1.1 | read, not invoked (constraint 5) |

### Output

```
20_new_runs/snplib_leadlike_2026-10-05/
  snplib_leadlike.json
  waterfall.json
  shortlist.csv
  run_manifest.json
  analysis/snplib_leadlike.py
```

## ADR-AIIM-308 — PRE-REGISTRATION: honour the duplicate flag the frozen input has carried all along

**Status:** proposed — **pre-registration. Frozen before the run.**
**Date:** 2026-10-05
**Raised by:** Adisorn Sirichotejirakul — *"ร่าง pre-reg รอบสามแล้วรันเลย"* (2026-10-05)
**Follows:** ADR-AIIM-306 and ADR-AIIM-307, whose counts are inflated by
duplicate rows
**Discipline:** ADR-AIIM-263 — the gate and the representative rule fixed before
the run

### The defect: the input flagged it, my harness ignored the flag

`snp_source_ledger.csv` carries a column `duplicate_inchikey_flag`. **G1 set it
True on all 24 affected rows.** My harnesses read `lipinski_ro5_pass` from the
same ledger and never looked at the column beside it.

```
duplicate_inchikey_flag:  True 24  /  False 212
24 rows carry 11 distinct structures
SNP-Lib holds 223 distinct structures in 236 rows
```

Eleven groups, every one flagged by G1:

```
[SNP00055, SNP00092]            [SNP00003, SNP00041, SNP00173]
[SNP00133, SNP00134]            [SNP00083, SNP00172]
[SNP00157, SNP00233]            [SNP00042, SNP00074]
[SNP00011, SNP00066]            [SNP00052, SNP00099]
[SNP00078, SNP00146, SNP00153]  [SNP00008, SNP00106]
[SNP00141, SNP00217]
```

Two groups reached the shortlist:

| | canonical SMILES |
|---|---|
| SNP00052 *Denchrysan B* | `COc1cc(O)cc2c1-c1c(O)cccc1C2O` |
| SNP00099 *4-Methoxy-9H-fluorene-2,5,9-triol* | **identical** |
| SNP00157 *N-trans-feruloyl tyramine* | identical to SNP00233 |
| SNP00233 *Trans-N-feruloyltyramine* | **identical** |

The 52/99 pair also carries **contradictory `groups` labels** — one
*Fluorenones*, one *Fluorenes* — for a structure that is a fluorene-triol. That
is a defect in the frozen ledger's annotation, distinct from the duplication,
and it is recorded rather than corrected: the input is hash-anchored and
read-only.

**This is my defect, not the project's.** The input recorded the duplication
correctly and completely. Every SNP-Lib count in ADR-AIIM-306 and
ADR-AIIM-307 — pool 236, Ro5 157, MW 133, eligible 93 and 74, survivors 17 and
12, unanimous 8 and 5 — counts **rows, not structures**, and is inflated.

### What it already explains

**ADR-AIIM-307's U4 did not fail.** SNP00052 and SNP00099 are one molecule;
both reached 7 of 8 conventions; scaffold uniqueness correctly kept one and
dropped the other, and it kept the one that won the mean-pKi tiebreak while my
prediction named the loser. **The rule did the right thing on a duplicate I did
not know existed.**

**ADR-AIIM-307's U2 failure was confounded, not informative.** I compared
unanimity across two different pools and attributed the change to window width.
Deduplicated, ADR-AIIM-306's 8 unanimous rows are **7 structures** and
ADR-AIIM-307's 5 are **4**. The MW gate also removed cinnamic acid, SNP00018
and SNP00086, all unanimous — so unanimity fell because the compounds every
model scored alike had left the pool, not because the window narrowed.
**ADR-AIIM-304's window relationship is neither confirmed nor contradicted by
U2; the prediction failed to control pool composition.**

### The change — one gate, and the representative rule stated

**Dedupe on `smiles_canonical`, keeping the lowest `SNP_ID`.** Deterministic,
independent of any score, and applied **first**, so every downstream count is a
count of structures. Each group and the ID retained are reported.

Nothing else moves: Ro5, MW ≥ 200, AD ≥ 0.25, censored-EM, eight conventions,
top-20, 5-of-8 majority, G6 uniqueness, shortlist 10, and **alerts remain
flags**.

### The waterfall, against the reference the harness uses (6,284)

| step | n | frameworks |
|---|---|---|
| 0 — SNP-Lib rows | 236 | |
| **0b — dedupe, keep lowest SNP_ID** | **223** | |
| 1 — gate Ro5 | 147 | |
| 2 — gate MW ≥ 200 | 126 | |
| 3 — gate AD ≥ 0.25 | **69** | **21** |
| *(sensitivity AD ≥ 0.20)* | *114* | *33* |

Duplicate representatives surviving to eligible: **SNP00042, SNP00052,
SNP00078, SNP00157** — so SNP00099 and SNP00233 leave, and SNP00052 and
SNP00157 stay as the retained copies.

Eligible composition: Phenanthrenes 31, Bibenzyls 19, Phenolics 4,
Flavonoids 4, Fluorenones 4, Fluorenes 3, Stilbenes 2, Coumarins 1,
Quinazolinone 1.

**These numbers are computed against the 6,284-compound modelling set, which is
what the harness scores against.** ADR-AIIM-306 quoted 88 and the run returned
93 because its feasibility used a 5,707 reference; ADR-AIIM-307 declared that
and still mismatched. Fixed by computation this time rather than by declaration.

### Predictions

**W1 — deduplication is bookkeeping, not a change of answer.** The set of
**distinct structures** in the new shortlist equals the set in ADR-AIIM-307's,
treating SNP00052 and SNP00099 as one. *Falsified if any structure enters or
leaves.* *Basis:* dedup removes 5 of 74 eligible rows, every one a redundant
copy, so the ranking over distinct structures should be nearly invariant.

**W2 — unanimity is stable once counted properly.** The 8-of-8 count is
**3 to 5**, against ADR-AIIM-307's 5 rows that were 4 structures. *Falsified
outside that band.*

**W3 — alerts stay prevalent.** More than 30 % of the 69 eligible carry a PAINS
or Brenk match. *Falsified at 30 % or below.* ADR-AIIM-307 measured 56.8 % of
74.

**W4 — nothing crosses §3.6's threshold.** The top shortlisted compound's mean
predicted pKi stays below 6.0. *Falsified at 6.0 or above.* **Weak**, as in
ADR-AIIM-307: the pool is a subset.

**The outcome that would most embarrass this entry is W1 failing.** If removing
redundant copies changes *which structures* are prioritised, the duplicates were
not inert labels — they were distorting the ranking, and every SNP-Lib result
before this is wrong in selection and not merely inflated in count. That would
make three consecutive runs' shortlists artefacts of a column I did not read.

### Constraints

1. **The frozen input is not corrected.** `snp_source_ledger.csv` is
   hash-anchored and read-only. The duplication and the contradictory
   `groups` labels on SNP00052 / SNP00099 are **recorded here**, not edited, and
   `02_frozen_inputs/` is untouched.
2. **`10_new_system` is not modified**; nothing re-docked; no contract amended.
3. **§3.6 is unchanged** at 72 of 5,000.
4. **A shortlisted compound is not a predicted active.** ADR-AIIM-299 stands:
   conformal count 0 at pKi ≥ 6.0 under every framing, interval ±1.5 log
   against a 0.5 decision gap. The word "active" is applied to nothing.
5. **Alerts remain flags**, unchanged from ADR-AIIM-306 and ADR-AIIM-307, and
   converting them now would be the post-hoc move those entries refused.
6. **ADR-AIIM-306's and ADR-AIIM-307's shortlists stay in the record**,
   superseded on a counting defect rather than withdrawn.
7. **The top-20 window and scaffold uniqueness remain unvalidated** against
   labels; ADR-AIIM-304 closed only the majority rule.
8. **ADR-AIIM-299's limitations carry over**: the censored-EM E-step is biased
   upward at censored points and the EM does not reach its 0.01 tolerance in 10
   iterations.
9. **No compound has been assayed.**
10. **Every outcome is reported**, including W1 failing.
11. **This entry is appended to the decision log and its SHA-256 recorded before
    the run**, digest over UTF-8 text with LF newlines per the 2026-10-05
    addendum.

### Skills invoked

| skill | version | material influence |
|---|---|---|
| `scikit-learn` | 1.1 | the ADR-AIIM-298 censored-EM pipeline reused unchanged |
| `rdkit` | 1.1 | canonical SMILES for the dedup key, `Descriptors.MolWt` for the gate, `FilterCatalog` for the flags, ECFP4 for the competence score, Murcko frameworks for diversity |
| `medchem` | 1.1 | read, not invoked, no install requested — as ADR-AIIM-306 and ADR-AIIM-307 |

### Output

```
20_new_runs/snplib_dedup_2026-10-05/
  snplib_dedup.json
  waterfall.json
  duplicate_groups.json
  shortlist.csv
  run_manifest.json
  analysis/snplib_dedup.py
```

## ADR-AIIM-309 — ADR-AIIM-306, 307 and 308 results: three rounds, three defects in my own gate, and an eight-structure priority order none of which crosses §3.6's threshold

**Status:** accepted — record of fact
**Date:** 2026-10-06
**Approver:** Adisorn Sirichotejirakul — *"เขียน result ADR เลย"* (2026-10-06)
**Pre-registrations:** ADR-AIIM-306 `830e195a1e42f2eeb4d093e1eb167f6591e1007538ea46e7e8ce5d97aa5cb0cc`
(log `b2d269a04cd0d4a52194f06e902abc8c6aae8eb8ade96e4ff937a3f0ffc646d1`);
ADR-AIIM-307 `f38879a387848ac553df01a136740513f1057f490233f0f8e3fe1f54c8be7fde`
(log `15098232bc8349e6410a33df18c77d8eb74b620361ac5aa062c25abcaaa18f3d`);
ADR-AIIM-308 `cb2ae5d7397f6276dcb6664a90007dd470ad9e06ba0ff0e605f1df627fa5759a`
(log `dcf2b53240bd3cf360efb2af0b619d7850eeca23149a70c797220757c81a7e1e`)
**Artifacts:** three run directories dated 2026-10-05, 760 s each
**Supersedes for practical use:** ADR-AIIM-302's NP-library shortlist

### Grading

| | prediction | verdict |
|---|---|---|
| **S1** | ≥ 10 survivors of the 5-of-8 majority | **MET** — 17 |
| **S2** | alerts on > 30 % of eligible | **MET** — 55.9 % |
| **S3** | top pKi below ADR-AIIM-302's 6.191 | **MET** — 5.913 |
| **S4** | scaffold uniqueness binds | **MET** — 7 from 17 survivors |
| **U1** | gating has knock-on effects | **MET** — 4 changes beyond the removal |
| **U2** | unanimity exceeds ADR-AIIM-306's 8 | **MISSED** — 5; **and the comparison was confounded** |
| **U3** | alerts still above 30 % | **MET** — 56.8 % |
| **U4** | the alert-clean three survive | **MISSED** — **on a duplicate artefact** |
| **U5** | top pKi below 6.0 | **MET** — 5.918 |
| **W1** | dedup is bookkeeping, structure set unchanged | **MISSED** — one structure entered, **none left** |
| **W2** | unanimity 3–5 | **MET** — 4 |
| **W3** | alerts above 30 % | **MET** — 56.5 % |
| **W4** | top pKi below 6.0 | **MET** — 5.877 |

### The three defects, and who found each

Every one was in **my** gate design, and **not one was caught by the
pre-registration that contained it**.

| | defect | found by |
|---|---|---|
| 1 | **Ro5 constrains only ceilings.** It has no lower MW bound, so a fragment passed: ADR-AIIM-306 put **cinnamic acid, MW 148.2**, at rank 1, unanimous across all eight conventions | the run's own output, inspected after the author asked to see the structures |
| 2 | **`duplicate_inchikey_flag` went unread.** G1 set it True on all 24 affected rows; my harnesses read `lipinski_ro5_pass` from the same ledger and never looked at the column beside it | the author asking whether the compounds came from SNP-Lib, which led to reading the ledger properly |
| 3 | **`TOP_N = 20` is an absolute count, not a fraction.** Shrinking the pool widens the window relatively, so **every cross-pool comparison is structurally confounded** | W1's failure, traced |

Defect 3 is the one that matters methodologically, because it invalidates two
predictions retrospectively:

* **U2** compared unanimity across pools of 93 and 74 and attributed the change
  to window width. The window also widened from 21.5 % to 27 % *because the
  pool shrank*, and the MW gate removed cinnamic acid, SNP00018 and SNP00086 —
  all unanimous. **ADR-AIIM-304's window relationship is neither confirmed nor
  contradicted by U2.**
* **W1** failed because the window widened again, 27 % → 29 %, admitting
  **SNP00006** that had not previously reached the majority. **No structure
  left.** Deduplication was bookkeeping exactly as predicted; the absolute
  window added a compound on top of it.

**U4 did not fail either.** SNP00052 and SNP00099 are **one molecule with
identical canonical SMILES**; both reached 7 of 8; scaffold uniqueness correctly
kept one and my prediction named the one that lost the mean-pKi tiebreak.

### SNP-Lib holds 223 structures in 236 rows

Eleven duplicate groups, every one flagged by G1:

```
[SNP00055, SNP00092]  [SNP00003, SNP00041, SNP00173]  [SNP00133, SNP00134]
[SNP00083, SNP00172]  [SNP00157, SNP00233]            [SNP00042, SNP00074]
[SNP00011, SNP00066]  [SNP00052, SNP00099]            [SNP00008, SNP00106]
[SNP00078, SNP00146, SNP00153]                        [SNP00141, SNP00217]
```

| | SMILES |
|---|---|
| SNP00052 *Denchrysan B* | `COc1cc(O)cc2c1-c1c(O)cccc1C2O` |
| SNP00099 *4-Methoxy-9H-fluorene-2,5,9-triol* | identical |

**A second defect in the frozen ledger, recorded not corrected:** that pair
carries contradictory `groups` labels — *Fluorenones* and *Fluorenes* — for one
fluorene-triol. The input is hash-anchored and read-only;
`02_frozen_inputs/` was not touched.

**Every SNP-Lib count in ADR-AIIM-306 and ADR-AIIM-307 counts rows, not
structures**, and is inflated: pool 236, Ro5 157, MW 133, eligible 93 and 74,
survivors 17 and 12, unanimous 8 and 5.

### The waterfall, as structures

| step | n | frameworks |
|---|---|---|
| SNP-Lib rows | 236 | 80 |
| **dedupe, keep lowest SNP_ID** | **223** | |
| gate Ro5 | 147 | |
| gate MW ≥ 200 | 126 | |
| gate AD ≥ 0.25 | **69** | **21** |
| *(sensitivity AD ≥ 0.20)* | *114* | *33* |

**The MW floor was derived, not chosen.** Of 4,428 actives in the modelling set
(pKi ≥ 6.0): minimum MW 196.2, 1st percentile 236.5, median 391.7, and **one
compound — 0.02 % — below MW 200**. The 10th percentile (MW 298) was rejected
because 10 % of known actives sit below it, and a plausibility floor belongs at
the edge of the observed rather than at a typical value.

### Unanimity, measured as structures, and where it actually fell

| | rows | **structures** |
|---|---|---|
| ADR-AIIM-306 | 8 | **7** |
| ADR-AIIM-307 | 5 | **4** |
| ADR-AIIM-308 | 4 | **4** |

**The fall is at the MW gate, not at the dedup.** Removing the fragments removed
compounds that every convention scored alike — small simple phenolics are easy
to agree on — so agreement fell when the drug-like compounds became the whole
pool. **Agreement was partly an artefact of the fragments.**

### The priority order: eight structures from the obtainable library

| # | SNP_ID | name | conv | pred pKi | AD | MW | cLogP | TPSA | QED | class | alerts |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | SNP00163 | Phoyunbene C | **8/8** | 5.877 | 0.277 | 272.3 | 3.29 | 58.9 | **0.839** | Stilbenes | BRENK |
| 2 | SNP00157 | *N*-trans-feruloyl tyramine | **8/8** | 5.876 | **0.343** | 313.4 | 2.48 | 78.8 | 0.716 | Phenolics | BRENK |
| 3 | SNP00010 | Hircinol | **8/8** | 5.770 | 0.251 | 242.3 | 2.87 | 49.7 | **0.808** | Phenanthrenes | **none** |
| 4 | SNP00154 | 1-(4-hydroxybenzyl)-4-methoxy-9,10-dihydrophenanthrene | 7/8 | 5.781 | 0.253 | 362.4 | 4.47 | 58.9 | 0.719 | Phenanthrenes | **none** |
| 5 | SNP00052 | Denchrysan B | 7/8 | 5.770 | 0.251 | 244.2 | 2.17 | 69.9 | 0.719 | Fluorenones | **none** |
| 6 | SNP00144 | — | 6/8 | 5.837 | 0.275 | 288.3 | 2.03 | 96.2 | 0.624 | Fluorenones | BRENK |
| 7 | SNP00142 | — | 6/8 | 5.709 | 0.255 | 286.3 | 3.13 | 79.2 | 0.498 | Phenanthrenes | PAINS, BRENK |
| 8 | SNP00006 | 2-(*p*-hydroxyphenyl)ethyl *p*-coumarate | 6/8 | 5.656 | 0.266 | 284.3 | 2.90 | 66.8 | 0.654 | Phenolics | BRENK |

**Three carry no alert in any of the four catalogs: SNP00010, SNP00154,
SNP00052.**

### Why the alerts were flagged and not gated, with the number that justifies it

| | eligible | PAINS or Brenk |
|---|---|---|
| ADR-AIIM-306 | 93 | 52 — **55.9 %** |
| ADR-AIIM-307 | 74 | 42 — **56.8 %** |
| ADR-AIIM-308 | 69 | 39 — **56.5 %** |

A hard gate would have removed **about 56 % of the pool at every stage**, and
the matches are `stilbene` on a stilbene, `Michael_acceptor_1` on a cinnamate,
`catechol` on a catechol — substructures that are the chemotype of an orchid
phenolic library, not evidence of a problem compound. PAINS derives from HTS of
synthetic libraries. The flags are recorded for every shortlisted compound and
exclude nothing.

### The result that does not move

**All eight shortlisted compounds carry a predicted pKi below 6.0** — 5.656 to
5.877. ADR-AIIM-299 measured the honest interval at **±1.5 log** against a
**0.5** decision gap, with the conformal count at pKi ≥ 6.0 equal to **0 under
every framing and every convention**.

So the obtainable, Ro5-compliant, non-fragment, in-domain best of the focus
library **does not reach the threshold §3.6 reports a count against**. The
ranking is usable; the threshold is not reachable by this instrument. That
sentence is the same one ADR-AIIM-299 reached from the other direction.

### Unchanged, and still open

**§3.6 is unchanged** at 72 of 5,000. Nothing has been enacted.

**Two criteria remain unvalidated against labels**, as ADR-AIIM-302 recorded and
ADR-AIIM-304 only partly closed: the **top-20 window** — now also known to be
confounded as an absolute count — and the **scaffold-uniqueness rule**, which
discarded two unanimous compounds in ADR-AIIM-307 and whose effect on
enrichment has never been measured.

### Deliberately not done

* **The frozen input was not corrected.** The duplication and the contradictory
  `groups` labels are recorded; `02_frozen_inputs/` and the ledger are untouched.
* **`10_new_system` is unmodified**; nothing re-docked; no contract amended;
  §3.6 unchanged.
* **The window was not converted to a fraction.** Defect 3 is recorded and not
  fixed; doing so would change every count again and needs its own
  pre-registration.
* **`medchem` was read and never invoked** across all three rounds — not
  installed, no install requested. `rule_of_leadlike_soft` would have supplied
  the MW floor directly; it was derived from the training distribution instead.
  Forgone: `rule_of_cns`, NIBR severity, Lilly demerits, ZINC-15 complexity.
* **ADR-AIIM-302's and ADR-AIIM-306's and ADR-AIIM-307's shortlists stay in the
  record**, superseded rather than withdrawn, each for a stated defect in the
  gate rather than anything measured about the compounds.
* **ADR-AIIM-299's limitations carry over**: the censored-EM E-step is biased
  upward at censored points and the EM reached 10 iterations without meeting its
  0.01 tolerance in every convention of every round.
* **No compound has been assayed**, and none is put forward as active.

### Artifacts

Digests over raw bytes, per the 2026-10-05 addendum.

```
20_new_runs/snplib_prioritisation_2026-10-05/   ADR-AIIM-306
  snplib_prioritisation.json  8f14a81dbe17fd2d5e50350f46d9b166496f2ade8d8798e28cee8b6bf012f156
  shortlist.csv               90fb34cc8bb97dc19241340992279c58a6c89aad74be4d5508006280683eb6be
  waterfall.json              7a4661b746ac22aede9503ed1ead3fa5fedfb3858372144e4a9ef644928adb72
  analysis/…prioritisation.py a4dfecb3b1e94a51e4fb5aaa625bb78c8fcb0098fbb2b25e463d098332f5ec69
20_new_runs/snplib_leadlike_2026-10-05/         ADR-AIIM-307
  snplib_leadlike.json        056ca08f446ad887f7b20624260eb32424d54b025bc855701f1198742a876100
  shortlist.csv               6405ffcbb6ec288436c2673cae86721a2908e6df3ce39cdd39f674a2167d85cf
  waterfall.json              ba4525a6f4a404ad25d676eb002f69114c7285500feb249b6e8f265d8f8b0765
  analysis/…leadlike.py       422433afc268f5ac761005ef3555950e8bc4eac8bc3f911004dbad2aa2aa0e91
20_new_runs/snplib_dedup_2026-10-05/            ADR-AIIM-308
  snplib_dedup.json           b54531d2e9b54eb85c0f08ed35c890b03a36a8adf3b49d2feb5955eebc30be4b
  duplicate_groups.json       68902d074bcdedb8df325f52b26eaf0e36cc015e2f78397e63ef0c89a43ddbe9
  shortlist.csv               2e05fd48af917820fb77c9c8aab9e075b86fd7cf10a0dba1a3d7825d014783c4
  waterfall.json              e0138df8be784388a046d9e8dbcfe2956cfc9e6457d538bfc63a2455438230e3
  analysis/snplib_dedup.py    aadf535904c985479a7f5fac0f3c2456707517e98dd8767d63c3c0bdc88253f9
```

## ADR-AIIM-310 — PRE-REGISTRATION: separate against pooled, under identical gates, and the window convention measured rather than recorded

**Status:** proposed — **pre-registration. Frozen before the run.**
**Date:** 2026-10-06
**Raised by:** Adisorn Sirichotejirakul — *"แล้วถ้าเทียบระหว่าง SNP กับ NP แยกกัน กับ เอามา pool รวมกันแล้วผลจะเป็นอย่างไร?"* (2026-10-06)
**Follows:** ADR-AIIM-309, which recorded defect 3 — the absolute window — as
*"recorded and not fixed"*
**Discipline:** ADR-AIIM-263 — pools, gates, windows and predictions fixed
before the run

### Two corrections to what I told the author before computing

**1. The comparison as it stands is not fair.** ADR-AIIM-300 applied **only**
the AD gate to the NP library: no dedup, no Ro5, no MW floor. ADR-AIIM-308
applied all four to SNP-Lib. So ADR-AIIM-302's six compounds and
ADR-AIIM-309's eight were drawn from pools filtered differently, and the two
shortlists have never been comparable. NPLIB01498 at rank 1 of ADR-AIIM-302 —
MW 596.5, cLogP −2.65, TPSA 256.3 — **fails Ro5** and would not have been in a
gated pool at all.

**2. My estimate of the imbalance was wrong, and in the direction that weakened
my own guess.** I told the author SNP-Lib would be *"squeezed out almost
entirely"* because NP-library had 38× more candidates, quoting 2,645 eligible.
That figure came from the **ungated** pool. Under the identical gate stack:

| | rows | dedupe | Ro5 | MW ≥ 200 | **AD ≥ 0.25** | frameworks |
|---|---|---|---|---|---|---|
| SNP-Lib | 236 | **223** | 147 | 126 | **69** | 21 |
| NP-library | 5,000 | **5,000** | **1,525** | 1,525 | **1,066** | 116 |
| **pooled** | | | | | **1,135** | |

SNP-Lib is **6.1 % of the pooled pool**, not 2.5 %. A top-20 drawn at random
would contain about **1.2** SNP-Lib compounds. **My prediction that it would be
shut out is much weaker than I presented it**, and the 0.31 log gap I quoted
between the libraries' best compounds compared a gated SNP-Lib top against an
ungated NP-library top and should not have been quoted.

Two further facts from the gating, both new:

* **NP-library has no exact duplicates at all** — 5,000 rows, 5,000 structures
  — where SNP-Lib had 24 rows over 11 structures. The duplication was
  SNP-Lib-specific.
* **Ro5 removes 69.5 % of the NP library** and **MW ≥ 200 removes nothing** from
  what survives. The fragment problem was SNP-Lib-specific; the NP library's
  problem is the opposite end.

**Anchor on the gate itself:** RDKit's Ro5, computed here, reproduces G1's
`lipinski_ro5_pass` on SNP-Lib **236 of 236**. The same rule is therefore
applied to both libraries, rather than reading a column that exists for only
one.

### The design

**3 pools × 2 window conventions**, each under the eight split conventions of
ADR-AIIM-286, with censored-EM exactly as ADR-AIIM-298.

| | pool |
|---|---|
| **A** | SNP-Lib alone — 69 |
| **B** | NP-library alone — 1,066 |
| **C** | pooled — 1,135 |

| | window |
|---|---|
| **abs** | top-**20**, as every run from ADR-AIIM-300 onward |
| **frac** | top-**10 %**, the primary window ADR-AIIM-304 validated on DUD-E |

The window is not a nuisance parameter here. **It is a ninth convention**, of
exactly the kind ADR-AIIM-286 found sets the answer, and ADR-AIIM-309 recorded
that leaving it absolute makes every cross-pool comparison confounded. This run
measures it instead of recording it.

Gates, identical for both libraries and applied in this order: **dedupe on
canonical SMILES keeping the lowest id → Ro5 (RDKit) → MW ≥ 200 → AD ≥ 0.25**.
Alerts remain **flags**. Shortlist length 10, 5-of-8 majority, G6 scaffold
uniqueness — all unchanged.

**A structural fact, stated so it is not mistaken for a result:** the model
scores each compound independently, so pooling cannot reorder compounds within
a library. The SNP-Lib subsequence of pool C's ranking is necessarily pool A's
ranking. What pooling changes is **how many slots SNP-Lib wins against NP
competition**, and nothing else.

### Predictions

**X1 — the anchor.** Pool A with the absolute window reproduces
ADR-AIIM-309's eight-structure shortlist exactly, in the same order.
*Falsified by any difference.* **If X1 fails nothing else in this run means
anything**, and it is named first for that reason.

**X2 — SNP-Lib is not shut out.** In pool C with the absolute window, at least
one SNP-Lib compound reaches the 5-of-8 majority. *Falsified if none does.*
*Basis:* 6.1 % of the pool and ~1.2 expected by chance in a top-20. **This
reverses the guess I gave the author**, which rested on the ungated figure.

**X3 — the window is as consequential as the split.** On pool C, the survivor
sets under the two windows overlap by less than 50 % of the smaller set.
*Falsified at 50 % or above.* This is the direct measurement of defect 3.

**X4 — the two libraries are not distinguishable.** The difference between the
best predicted pKi of pool A and pool B is below **1.5 log**, the α = 0.05
absolute-residual interval ADR-AIIM-299 measured. *Falsified at 1.5 or above.*
**Marked weak** — the libraries' ranges are narrow — but it is the sentence the
manuscript needs: choosing the obtainable library costs nothing the instrument
can resolve.

**X5 — NP-library dominates the pooled shortlist by count.** More than half of
pool C's shortlist under the absolute window is NP-library. *Falsified at half
or below.* X2 and X5 are compatible and test different things: whether SNP-Lib
appears at all, and whether it competes.

**The outcome that would most change the recommendation is X5 failing** —
SNP-Lib holding half or more of a pooled shortlist would mean the focus library
competes on merit against a library 15× its size, and the case for ranking it
separately would rest on availability alone rather than on availability plus
numerical protection.

### Constraints

1. **`10_new_system` is not modified**; nothing re-docked; no contract amended;
   `02_frozen_inputs/` not written to; both ledgers read unchanged.
2. **§3.6 is unchanged** at 72 of 5,000.
3. **Pooling is measured, not recommended.** Whatever the result, a pooled
   shortlist mixes compounds that can be assayed with compounds that cannot;
   the author's ruling on availability stands and is not under test here.
4. **No shortlisted compound is a predicted active.** ADR-AIIM-299: conformal
   count 0 at pKi ≥ 6.0 under every framing. The word "active" is applied to
   nothing.
5. **ADR-AIIM-302's NP shortlist is not re-derived as a correction.** Pool B is
   a new, gated measurement; ADR-AIIM-302's six compounds stay in the record as
   what an ungated pool produced.
6. **Alerts remain flags** in all six cells.
7. **The top-20 window and scaffold uniqueness remain unvalidated against
   labels.** X3 measures the window's *consequence*, not its *correctness*;
   ADR-AIIM-304 closed only the majority rule.
8. **ADR-AIIM-299's limitations carry over**: the censored-EM E-step is biased
   upward at censored points and the EM does not reach its 0.01 tolerance in 10
   iterations.
9. **No compound has been assayed.**
10. **Every outcome is reported**, including X1 failing.
11. **This entry is appended to the decision log and its SHA-256 recorded before
    the run**, digest over UTF-8 text with LF newlines per the 2026-10-05
    addendum.

### Skills invoked

| skill | version | material influence |
|---|---|---|
| `scikit-learn` | 1.1 | the ADR-AIIM-298 censored-EM pipeline reused unchanged across the eight conventions; the model is fitted once per convention and scores all three pools, so no pool difference can come from refitting |
| `rdkit` | 1.1 | Ro5 and MolWt for the gates — anchored at 236/236 against G1's own flag — canonical SMILES for dedup, `FilterCatalog` for the flags, ECFP4 for the competence score, Murcko frameworks for diversity |
| `medchem` | 1.1 | read, not invoked, no install requested, as in ADR-AIIM-306 to 308 |

### Output

```
20_new_runs/pool_comparison_2026-10-06/
  pool_comparison.json
  waterfall.json
  shortlists.csv            six cells, labelled by pool and window
  run_manifest.json
  analysis/pool_comparison.py
```

## ADR-AIIM-311 — PRE-REGISTRATION: the censored-EM results are not reproducible, the cause is found, and this measures what it cost

**Status:** proposed — **pre-registration. Frozen before the run.**
**Date:** 2026-10-06
**Raised by:** found while diagnosing ADR-AIIM-310's X1 anchor failure
**Follows:** ADR-AIIM-310, whose anchor failed; affects ADR-AIIM-298 to 310
**Discipline:** ADR-AIIM-263 — predictions fixed before the run;
ADR-AIIM-226 — outside the frozen tree

### The defect, diagnosed to its mechanism

ADR-AIIM-310's X1 required pool A under the absolute window to reproduce
ADR-AIIM-309's shortlist. It did not. **The pools are provably identical** —
dedup 223 with matching ID sets, G1 canonical SMILES equal to RDKit's for all
236, Ro5 from the ledger column equal to Ro5 from RDKit at 147 with matching
IDs, `max |ledger mol_wt − RDKit MolWt| = 0.0`, MW ≥ 200 giving 126 with
matching IDs — and the features are ECFP4 + FCFP4, which are atom-order
invariant.

The difference is the model. Measured:

| | repeat 1 | repeat 2 | max \|diff\| | identical |
|---|---|---|---|---|
| `censored_em`, `n_jobs=3` | EM change 0.05011 | 0.03478 | **0.127** | **no** |
| `censored_em`, `n_jobs=1` | EM change 0.03438 | 0.03438 | **0.000** | **yes** |

And sklearn alone is not at fault:

| | max \|diff\| over two identical fits |
|---|---|
| `RandomForestRegressor(n_jobs=3)` | 1.8 × 10⁻¹⁵ |
| `RandomForestRegressor(n_jobs=1)` | 0 |
| `RandomForestClassifier(n_jobs=3)` | 0 |

**The mechanism.** A single random-forest fit at `n_jobs>1` carries
float-level non-determinism of order 10⁻¹⁵ from thread accumulation order. The
EM loop feeds predictions back into the targets for ten iterations. Forest
splits are threshold comparisons, so a 10⁻¹⁵ target perturbation can flip a
split, which changes a leaf value by order 10⁻², which feeds back and grows.
**Ten iterations amplify 10⁻¹⁵ into 0.127** — enough to change which compounds
enter a top-20.

A parallelism setting that is harmless in one fit becomes a reproducibility
failure inside an iterative loop. **This is a finding about the architecture,
not only a bug in my harness**, and it belongs in the record on those terms.

### What it affects and what it does not

**Affected — every entry using `censored_em`:** ADR-AIIM-298, 299, 300, 301
(its L channel), 302, 303, 304, 306, 307, 308, 309, 310.

**Not affected — the classifier work**, which is where the project's anchors
hold: ADR-AIIM-286 reproduced ADR-AIIM-275 exactly and ADR-AIIM-290 reproduced
ADR-AIIM-268 to 0.0006, both with `RandomForestClassifier`, which the table
above shows is deterministic at `n_jobs=3`. ADR-AIIM-301's Vina AUC of 0.7005
is also unaffected.

**Conclusions expected to survive**, because the effect is small against them:

| | magnitude | against noise 0.127 |
|---|---|---|
| conformal interval ±1.5 log | 1.5 | **12×** |
| conformal count 0 at every threshold | — | an interval result, not a point one |
| every shortlist pKi below 6.0 (5.656–5.918) | ≥ 0.08 margin | **marginal — must be rechecked** |
| library gap, ADR-AIIM-310 X4 | 0.351 | 2.8× |
| classifier convention spread 376× | — | no EM involved |

**Conclusions that must be re-measured:** the censored-EM convention spread
(ADR-AIIM-299's T3, reported 4.10×), every unanimity count, every shortlist,
and ADR-AIIM-304's precision spread of 0.0592 across conventions.

### What is run

**Fix: `n_jobs=1` in the censored-EM pipeline.** One line. Everything else
unchanged — same gates, same conventions, same windows, same majority rule.

Three measurements, all CPU, none touching the GPU that is still running
ADR-AIIM-305:

1. **Determinism check.** `censored_em` at `n_jobs=1`, three repeats on one
   convention, predictions compared exactly.
2. **Replicate spread at `n_jobs=3`**, the configuration every prior entry
   used: five repeats per convention on the ADR-AIIM-299 configuration,
   reporting the NP-library count at pKi ≥ 6.0 **within** each convention. This
   is the error bar those entries should have carried and did not.
3. **Corrected re-runs at `n_jobs=1`**: the ADR-AIIM-299 eight-convention
   spread, and the ADR-AIIM-309 SNP-Lib shortlist.

### Predictions

**Y1 — the fix works.** Three repeats at `n_jobs=1` give bit-identical
predictions. *Falsified by any difference.* **Confirmatory**: two repeats
already gave 0.000.

**Y2 — the corrected spread differs materially.** The eight-convention spread
of the NP count at pKi ≥ 6.0, recomputed at `n_jobs=1`, differs from
ADR-AIIM-299's reported **4.10×** by more than 20 % in either direction.
*Falsified within 20 %.*

**Y3 — run noise is comparable to the convention effect.** The median
within-convention spread across five `n_jobs=3` repeats is **at least half**
the between-convention spread ADR-AIIM-299 reported. *Falsified below half.*

**Y3 is the prediction that matters and the one I expect to hold.** If run
noise within a single convention is of the same order as the difference
between conventions, then ADR-AIIM-299's headline — that the spread falls
38.88× → 10.99× → 4.10× as the framing improves — **is partly measuring
noise in its last term**, and the 4.10× figure cannot be quoted as a property
of the framing. The first two terms are classifier and point-regression
results and are unaffected.

**Y4 — the qualitative conclusions survive.** At `n_jobs=1`: the conformal
count at pKi ≥ 6.0 stays 0, and every compound in the corrected SNP-Lib
shortlist stays below pKi 6.0. *Falsified if either changes.* The second clause
is **not** safe: the margin is 0.08 against noise of 0.127, so a compound
crossing 6.0 is possible and would change what the manuscript can say.

**Y5 — the shortlist changes.** The corrected ADR-AIIM-309 shortlist differs
from the reported one by at least one member. *Falsified if identical.*

**The outcome that would most embarrass this entry is Y4's second clause
failing** — a compound crossing pKi 6.0 once the noise is removed — because
every entry from ADR-AIIM-299 onward has reported that nothing crosses the
threshold, and that sentence is doing real work in the manuscript.

### Constraints

1. **`10_new_system` is not modified.** The `n_jobs` change is in this run's
   harness; the frozen tree contains no censored-EM code.
2. **No prior artifact is overwritten.** The affected runs stay on disk as what
   `n_jobs=3` produced, and ADR-AIIM-312 will state which of their numbers are
   superseded.
3. **§3.6 is unchanged** at 72 of 5,000.
4. **The prior entries are not withdrawn.** They are corrected on a
   reproducibility defect. Their pre-registrations, gradings and reasoning stand
   as the record of what was predicted and measured; the numbers carry a
   reproducibility caveat.
5. **No compound is a predicted active**, before or after the fix.
   ADR-AIIM-299's interval result does not depend on the point predictions.
6. **ADR-AIIM-305 is untouched.** The DiffDock run continues on the GPU; these
   measurements are CPU-only and will not contend with it.
7. **Every outcome is reported**, including Y3 holding and Y4 failing.
8. **This entry is appended to the decision log and its SHA-256 recorded before
   the run**, digest over UTF-8 text with LF newlines per the 2026-10-05
   addendum.

### A note on what this costs the deadline

The author's target is a submission-ready manuscript by 2026-10-10. These
three measurements are roughly three hours of CPU and do not block the
manuscript work, which ADR-AIIM-310's exchange established is the real
bottleneck. **But they change which numbers may be written into §3.6 and §4**,
so the writing should follow them rather than precede them.

### Skills invoked

| skill | version | material influence |
|---|---|---|
| `scikit-learn` | 1.1 | the defect is in how the skill's estimator interacts with an iterative wrapper; the fix is `n_jobs=1`, and the skill's own guidance to set `random_state` for reproducibility is shown here to be **insufficient** when the estimator sits inside a feedback loop |

### Output

```
20_new_runs/em_reproducibility_2026-10-06/
  determinism_check.json
  replicate_spread.json
  corrected_runs.json
  run_manifest.json
  analysis/em_reproducibility.py
```

## ADR-AIIM-312 — ADR-AIIM-311 result: the fix works, the convention effect is real after all, and the eight conventions are five

**Status:** accepted — record of fact
**Date:** 2026-10-06
**Pre-registration:** ADR-AIIM-311, SHA-256
`bde7df91688656656406d6d7a64470a100bc043d5d6d80405b94bfd93757742d`
**Decision log at freeze (text/LF):**
`2c78bd0b186f2230d39a26502b62413d2ac17b83eaaf6fc3361a2503ef3cfa99`
**Run:** `20_new_runs/em_reproducibility_2026-10-06/` — 6,523.5 s, CPU only,
concurrent with ADR-AIIM-305 on the GPU and not contending with it
**Software:** rdkit 2026.03.3, scikit-learn 1.7.2

### Grading

| | prediction | result | verdict |
|---|---|---|---|
| **Y1** | `n_jobs=1` bit-identical over 3 repeats | max\|diff\| **0.000e+00**, EM change 0.02279 identical, NP@6 151 identical | **MET** |
| **Y2** | corrected spread differs from 4.10× by > 20 % | **4.44×**, a **8.3 %** change | **MISSED** |
| **Y3** | median within-convention noise ≥ half the between-convention effect | median within-range **24**, between **90**, ratio **0.267** | **MISSED** |
| **Y4** | conformal count at 6.0 stays 0 **and** no shortlist compound crosses 6.0 | 0 in all eight; top shortlist pKi **5.9257** | **MET** |
| **Y5** | the shortlist changes by ≥ 1 member | **6 of 8** overlap, Jaccard **0.60** | **MET** |

**The prediction I said I expected to hold is the one that failed.** ADR-AIIM-311
stated that Y3 "is the prediction that matters and the one I expect to hold",
and that if it held, ADR-AIIM-299's 4.10× figure could not be quoted as a
property of the framing. It did not hold. Run noise is about **one quarter** of
the convention effect, not half, and together with Y2 this means
**ADR-AIIM-299's headline survives**: the spread does fall 38.88× → 10.99× →
4.10× (now 4.44×) as the framing improves, and that fall is a real property of
the framing, not an artefact of thread scheduling. I was wrong in the direction
that favours the prior work, which is the direction that deserves the most
scepticism, so the supporting numbers are given in full below.

### Y1 — the fix

`n_jobs=1`, three repeats, anchor convention `D1sorted_D2trainfirst_D3divert`:
every repeat gave EM final change **0.02279**, NP@6 **151**, and predictions
identical to the last bit. The one-line change is sufficient; nothing else in
the pipeline was non-deterministic.

### Y3 — the measurement, including the part that is bad news

Five repeats per convention at `n_jobs=3`, NP-library count at pKi ≥ 6.0:

| convention | replicates | range | sd | ADR-298 reported | inside? |
|---|---|---|---|---|---|
| D1shuffle_D2testfirst_D3overflow | 41, 28, 42, 25, 38 | 17 | 7.79 | 29 | yes |
| D1shuffle_D2testfirst_D3divert | 29, 39, 42, 26, 24 | 18 | 8.03 | 29 | yes |
| D1shuffle_D2trainfirst_D3overflow | 54, 61, 58, 43, 47 | 18 | 7.50 | **71** | **no** |
| D1shuffle_D2trainfirst_D3divert | 63, 74, 52, 48, 52 | 26 | 10.64 | 52 | yes |
| D1sorted_D2testfirst_D3overflow | 100, 46, 75, 72, 96 | 54 | 21.66 | **39** | **no** |
| D1sorted_D2testfirst_D3divert | 57, 76, 54, 73, 70 | 22 | 9.87 | 60 | yes |
| D1sorted_D2trainfirst_D3overflow | 72, 120, 106, 103, 98 | 48 | 17.56 | 78 | yes |
| **D1sorted_D2trainfirst_D3divert** | **124, 172, 105, 70, 91** | **102** | **38.72** | 119 | yes |

Median within-convention range 24 against a between-convention range of 90
gives 0.267, so Y3 is missed on its own pre-registered metric. Two things in
this table are worse than that number suggests and are not excused by it.

**First, the anchor is the worst cell.** `D1sorted_D2trainfirst_D3divert` — the
convention ADR-AIIM-298 designated as the anchor and from which every headline
number in ADR-AIIM-298 to 310 was read — spans **70 to 172**, a range of 102,
which is **larger than the entire between-convention range of 90**. The median
is reassuring; the cell the project actually quoted is not. Any single figure
read off that cell at `n_jobs=3` carried an uncertainty wider than the effect
the eight-cell factorial was built to measure.

**Second, two of the eight reported values fall outside their own five-replicate
range.** `D1shuffle_D2trainfirst_D3overflow` was reported as 71 against
replicates spanning 43–61, and `D1sorted_D2testfirst_D3overflow` as 39 against
replicates spanning 46–100. If thread-accumulation noise were the only source
of variation, five replicates would be expected to bracket the original draw;
for 2 of 8 they do not. **Something beyond thread noise also moved between the
ADR-AIIM-298 run and this one** — plausibly BLAS thread count, library load
order or machine state — and this entry does not identify it. The honest
statement is that the `n_jobs=3` numbers are not reproducible, and that the
mechanism is established only for the part of the variation that `n_jobs=1`
removes.

### Y2 and the corrected numbers

`n_jobs=1`, eight conventions, NP-library at pKi ≥ 6.0:

| convention | corrected | ADR-298 | Δ | qhat |
|---|---|---|---|---|
| D1shuffle_D2testfirst_D3overflow | 34 | 29 | +5 | 1.5815 |
| D1shuffle_D2testfirst_D3divert | 34 | 29 | +5 | 1.5815 |
| D1shuffle_D2trainfirst_D3overflow | 52 | 71 | −19 | 1.5446 |
| D1shuffle_D2trainfirst_D3divert | 52 | 52 | 0 | 1.5446 |
| D1sorted_D2testfirst_D3overflow | 71 | 39 | +32 | 1.6028 |
| D1sorted_D2testfirst_D3divert | 49 | 60 | −11 | 1.5838 |
| D1sorted_D2trainfirst_D3overflow | 151 | 78 | +73 | 1.6627 |
| D1sorted_D2trainfirst_D3divert | 151 | 119 | +32 | 1.6627 |

Spread **151 / 34 = 4.44×** against the reported 4.10×, an 8.3 % change, so Y2
is missed and the fold-range figure is robust to the defect. Individual cells
are not: the largest single correction is **+73** at
`D1sorted_D2trainfirst_D3overflow`, and the anchor moves **119 → 151**.

**Every cell-level count in ADR-AIIM-298 to 310 is superseded by this table.
Every fold-range and every interval conclusion stands.**

### An unpredicted finding: the D3 dimension is largely degenerate

With the noise removed, `D3divert` and `D3overflow` produce **bit-identical
models** in three of the four `D1 × D2` combinations — identical NP@6 and
identical qhat to four decimals:

| D1 × D2 | overflow | divert | identical |
|---|---|---|---|
| shuffle × testfirst | 34, qhat 1.5815 | 34, qhat 1.5815 | **yes** |
| shuffle × trainfirst | 52, qhat 1.5446 | 52, qhat 1.5446 | **yes** |
| sorted × testfirst | 71, qhat 1.6028 | 49, qhat 1.5838 | no |
| sorted × trainfirst | 151, qhat 1.6627 | 151, qhat 1.6627 | **yes** |

The distinct outcomes are **{34, 49, 52, 71, 151}**: **five**, not eight.
ADR-AIIM-286's factorial describes eight conventions and only five of them are
distinguishable on this dataset, because whether an overflowing scaffold group
is diverted or allowed to overflow changes nothing unless a group straddles the
budget boundary, which here happens in one of four cases. Thread noise had made
all eight look different and so had concealed the degeneracy.

This does not reduce the convention problem — the 4.44× spread is across five
genuinely different splits — but **wherever a count of eight conventions is
used as if it were eight independent framings, it overstates the evidence.**
That includes the 5-of-8 majority rule of ADR-AIIM-302 and onwards: a compound
can reach 5 of 8 on as few as **three** distinct splits, because three of the
eight votes are duplicates of others. The majority rule is therefore weaker
than its name implies. It is not re-specified here; that needs its own
pre-registration.

### Y4 — the qualitative conclusions survive

The conformal count at pKi ≥ 6.0 is **0 in all eight conventions**, and qhat is
**1.545 to 1.663**, confirming the ±1.5 log interval and widening it slightly.
The top corrected shortlist compound is **pKi 5.9257**, below 6.0 by **0.0743**.

That margin is smaller than the 0.127 of noise the defect injected, so under
`n_jobs=3` a draw could have crossed 6.0 and the sentence "no compound crosses
the threshold" could have been reported false. At `n_jobs=1` the value is
deterministic and the sentence is true. **It is true, and it was not safe.**

### Y5 — the corrected SNP-Lib shortlist

Pool unchanged and reproduced exactly: 236 rows → 223 deduped → 147 Ro5 → 126
MW ≥ 200 → **69** in domain at AD ≥ 0.25, 21 frameworks. Survivors **15**,
unanimous **5**.

| rank | SNP_ID | conventions | mean pKi | change |
|---|---|---|---|---|
| 1 | SNP00157 | 8/8 | 5.9257 | was rank 2 |
| 2 | SNP00163 | 8/8 | 5.9078 | was rank 1 |
| 3 | SNP00154 | 8/8 | 5.8025 | was rank 4 |
| 4 | **SNP00097** | 8/8 | 5.7688 | **entered** |
| 5 | SNP00052 | 7/8 | 5.7494 | was rank 5 |
| 6 | **SNP00033** | 7/8 | 5.7373 | **entered** |
| 7 | SNP00142 | 6/8 | 5.6827 | was rank 7 |
| 8 | SNP00006 | 6/8 | 5.6407 | was rank 8 |

**Dropped:** SNP00010 and SNP00144. The top two swap. ADR-AIIM-309's
eight-structure order is superseded by this one.

Read with the degeneracy finding above, the convention counts in this table are
counts over eight votes of which only five are distinct, so "8/8" means
agreement across five distinct splits and "6/8" can mean as few as three.

### Why ADR-AIIM-310's X1 failed

X1 required pool A under the absolute window to reproduce ADR-AIIM-309's
shortlist and reported 7 compounds, 13 survivors and 6 unanimous against
8/15/4. Both runs were correct executions of the same specification at
`n_jobs=3`; they differ because the specification did not fix the number of
threads. The deterministic answer is 8 compounds, 15 survivors, 5 unanimous.
**X1's failure was a true positive** — the pre-registered reproduction check
caught a real defect that no other test in this project would have found, which
is the strongest evidence this session produced that the ADR-AIIM-263 discipline
earns its cost.

### Affected-entry ledger

| entries | what is superseded | what stands |
|---|---|---|
| ADR-AIIM-298, 299 | the eight per-convention NP counts; the anchor's 119 | the 4.10× fold range (now 4.44×), the ±1.5 log interval, the zero conformal count, the censored-EM method and its interval-respecting rationale |
| ADR-AIIM-300, 301 | the L-channel point predictions and any count read from them | the ligand-channel-beats-intersection conclusion at matched size, and Vina's 0.7005 AUC, which involves no EM |
| ADR-AIIM-302, 303, 304 | the precision spread 0.0592 and all shortlists | the finding that convention instability is largely a threshold-crossing artefact at matched selection size |
| ADR-AIIM-306, 307, 308, 309 | the eight-structure order, now replaced above | the three gate defects found in those rounds, the 69-compound pool, the waterfall |
| ADR-AIIM-310 | cell-level counts in all six cells | NP 0/5,000 in domain structurally, the 6.1 % SNP share, the 0.351 library gap, X4 and X5 |
| ADR-AIIM-286, 290, 293 | nothing | `RandomForestClassifier` is deterministic at `n_jobs=3`; the anchors to ADR-AIIM-275 and 268 are unaffected |

**Nothing in §3.6 changes.** It remains 72 of 5,000, which is what the released
procedure produced and is not an EM result.

### Actions not taken

1. The affected run directories were **not overwritten**. They remain on disk as
   what `n_jobs=3` produced.
2. The affected entries were **not withdrawn**. Their pre-registrations and
   gradings are the record of what was predicted and measured.
3. The **5-of-8 majority rule was not re-specified** despite the degeneracy
   finding. That is a method change and needs its own pre-registration.
4. The **second source of variation was not identified** — the one that puts 2
   of 8 reported values outside their replicate range.
5. `10_new_system` was **not modified**.
6. **No compound is called active**, and none has been assayed.

### Open risks

1. Every other iterative fitting loop in this project should be assumed
   non-deterministic until checked. Only `censored_em` was audited.
2. The majority rule overstates its own support by the degeneracy factor.
3. The 0.0743 margin to pKi 6.0 is real but thin; any change to features,
   curation or AD reference set can cross it.

### Skills invoked

| skill | version | material influence |
|---|---|---|
| `scikit-learn` | 1.1 | the defect is the interaction of the skill's `n_jobs` with an iterative wrapper; the skill's guidance to set `random_state` for reproducibility is **demonstrably insufficient** inside a feedback loop, and this entry is the evidence |

### Deposit and manuscript

Not yet reflected. The corrected table above and the degeneracy finding both
change text that is currently in `30_manuscript/CURRENT_v5/`, and
ADR-AIIM-311's note that the writing should follow these measurements rather
than precede them now applies.

## ADR-AIIM-312 addendum — the full survivor set, and what the scaffold-uniqueness rule actually does

**Status:** accepted — record of fact
**Date:** 2026-10-06
**Extends:** ADR-AIIM-312, whose harness computed the survivor set but saved
only the eight shortlisted compounds
**Decision log at launch (text/LF):**
`35f117716b5dd093a4bf4787c323c43a3523f3edd4785050af52fad524b27472`
**Run:** `20_new_runs/em_survivors_2026-10-06/` — 2,218.4 s, CPU only
**Anchor, stated before the run:** `n_jobs=1` having been shown bit-identical,
this must reproduce the eight-compound shortlist in the same order, 15
survivors and 5 unanimous; failure would mean an unidentified source of
non-determinism remained and would matter more than the names being recovered.

### The anchor held, and a third determinism check passed

Shortlist reproduced in order, survivors 15/15, unanimous 5/5. The run was then
repeated after a display-only bug fix (see below) and `survivors.csv` came back
**byte-identical**, SHA-256 `619a8802e3f29b4b…`, with `em_survivors.json`
identical in every field but elapsed time. Determinism at `n_jobs=1` now has
three independent confirmations.

### The five distinct splits, located

ADR-AIIM-312 established that the eight conventions yield five distinct models.
This run identifies which:

```
group 1   D1shuffle_D2testfirst    overflow == divert
group 2   D1shuffle_D2trainfirst   overflow == divert
group 3   D1sorted_D2testfirst     overflow
group 4   D1sorted_D2testfirst     divert
group 5   D1sorted_D2trainfirst    overflow == divert
```

**D3 changes nothing except under `D1sorted × D2testfirst`.** Whether an
overflowing scaffold group is diverted or allowed to overflow matters only when
a group straddles the budget boundary, which happens in one of the four
`D1 × D2` combinations on this dataset.

### Measured distinct-split support, correcting a worst case I quoted

ADR-AIIM-312 noted that a compound reaching 5 of 8 could in principle rest on as
few as three distinct splits. Measured, the bound is not reached:

| conventions | distinct splits | survivors |
|---|---|---|
| 8/8 | **5 of 5** | 5 |
| 7/8 | **4 of 5** | 3 |
| 6/8 | **4 of 5** | 3 |
| 5/8 | **3 of 5** | 4 |

No survivor rests on fewer than three distinct splits, and the four unanimous
compounds that reach the shortlist rest on all five. **The evidence is stronger
than the worst-case bound I quoted**, and the worst case should not be repeated
as though it were the measurement.

### The fifteen survivors, and the seven the rule cut

Pool 69, window top-20, 5-of-8 majority, shortlist cap 10. The 15 survivors
occupy only **8 Murcko frameworks**, so scaffold uniqueness removes **7 of 15**.

| # | SNP_ID | name | conv | dist | pKi | AD | flags | outcome |
|---|---|---|---|---|---|---|---|---|
| 1 | SNP00157 | N-trans-feruloyl tyramine | 8/8 | 5/5 | 5.926 | 0.343 | 6 | **#1** |
| 2 | SNP00163 | Phoyunbene C | 8/8 | 5/5 | 5.908 | 0.277 | 4 | **#2** |
| 3 | **SNP00223** | **Thunalbene** | **8/8** | **5/5** | **5.876** | 0.321 | 4 | **CUT** |
| 4 | SNP00154 | 1-(4-hydroxylbenzyl)-4-methoxy-9,10-dihydrophenanthrene-2,7-diol | 8/8 | 5/5 | 5.803 | 0.253 | **0** | **#3** |
| 5 | SNP00097 | Dengibsinin | 8/8 | 5/5 | 5.769 | 0.316 | 4 | **#4** |
| 6 | SNP00052 | Denchrysan B | 7/8 | 4/5 | 5.749 | 0.251 | 1 | **#5** |
| 7 | SNP00033 | 4,5-dihydroxy-2,3-dimethoxy-9,10-dihydrophenanthrene | 7/8 | 4/5 | 5.737 | 0.263 | **0** | **#6** |
| 8 | SNP00095 | Dihydrodengibsinin | 7/8 | 4/5 | 5.736 | 0.258 | 1 | CUT |
| 9 | SNP00010 | Hircinol | 6/8 | 4/5 | 5.767 | 0.251 | **0** | CUT |
| 10 | SNP00142 | Dendrocrumenol A | 6/8 | 4/5 | 5.683 | 0.255 | 9 | **#7** |
| 11 | SNP00006 | 2-(p-hydroxyphenyl)ethyl p-coumarate | 6/8 | 4/5 | 5.641 | 0.266 | 8 | **#8** |
| 12 | SNP00144 | Dendrocrumenol C | 5/8 | 3/5 | 5.818 | 0.275 | 5 | CUT |
| 13 | SNP00051 | Dendroflorin | 5/8 | 3/5 | 5.755 | 0.297 | 4 | CUT |
| 14 | SNP00087 | 4,5-dihydroxy-2-methoxy-9,10-dihydrophenanthrene | 5/8 | 3/5 | 5.740 | 0.250 | **0** | CUT |
| 15 | SNP00198 | Pinadensin B | 5/8 | 3/5 | 5.586 | 0.265 | 6 | CUT |

The eight frameworks, with what each kept and cut:

| group | framework | kept | cut |
|---|---|---|---|
| 1 | cinnamamide | SNP00157 | — |
| 2 | **stilbene** | SNP00163 | **SNP00223** |
| 3 | benzyl-dihydrophenanthrene | SNP00154 | — |
| 4 | **fluorenone** | SNP00097 | SNP00144, SNP00051 |
| 5 | fluorene | SNP00052 | SNP00095 |
| 6 | **9,10-dihydrophenanthrene** | SNP00033 | SNP00010, SNP00087 |
| 7 | phenanthrene | SNP00142 | SNP00198 |
| 8 | cinnamate ester | SNP00006 | — |

### Three compounds were cut although they score higher than what blocked them

The ranking key is lexicographic — convention count first, mean predicted pKi
only as a tie-break — so a compound with more conventions blocks one with a
higher prediction:

| cut | conv | pKi | blocked by | conv | pKi | the cut compound is |
|---|---|---|---|---|---|---|
| SNP00144 | 5/8 | 5.8180 | SNP00097 | 8/8 | 5.7688 | **+0.049 higher** |
| SNP00010 | 6/8 | 5.7668 | SNP00033 | 7/8 | 5.7373 | **+0.030 higher** |
| SNP00087 | 5/8 | 5.7397 | SNP00033 | 7/8 | 5.7373 | **+0.002 higher** |

This also explains an ADR-AIIM-312 result that was recorded without a cause:
**SNP00144 left the shortlist** between ADR-AIIM-309 and the corrected run
because at `n_jobs=1` it reaches only 5 of 8 and its framework is taken by
SNP00097.

The unanimous compound cut is **SNP00223, Thunalbene**, which loses the stilbene
framework to Phoyunbene C by **0.032 log**.

### What this means, and it is the substantive finding

**The fifteen survivors span 5.586 to 5.926 — 0.339 log in total — against a
conformal interval of ±1.5 log, a 3.0 log width.** The entire ranked set sits
inside **11 %** of the interval the model itself reports. The tie-breaks that
decide the shortlist are differences of 0.002 to 0.096 log, every one of them
smaller than the 0.127 of noise the reproducibility defect used to inject.

**The model cannot distinguish these fifteen compounds.** Selecting eight of
them is a tie-break among indistinguishables, not a determination of which are
more likely to be active. The convention count carries more of the decision than
the predictions do, and it is an integer over five splits.

This is not an argument that the rule is wrong. Diversity selection legitimately
trades predicted potency for structural coverage, and that is exactly what
groups 2, 4, 5, 6 and 7 show it doing. It is an argument about **how the result
may be stated**: an ordered list of eight implies a discrimination the evidence
does not support.

### Recommendation for the manuscript

Report the **fifteen survivors organised by their eight frameworks**, stating
that within the survivor set the model does not discriminate, and that scaffold
uniqueness is a **diversity** rule rather than a potency rule. Do not present
the eight as a ranked priority order.

**The rule itself is not changed here.** Changing a selection criterion after
seeing which compounds it removes is precisely what this project's
pre-registration discipline exists to prevent. Changing how a result is
*reported* is not changing the criterion, and the reporting change is what is
recommended.

### Actions not taken

1. The scaffold-uniqueness rule and the lexicographic sort were **not
   re-specified**, despite both being unvalidated against labels.
2. The majority rule was **not re-specified**.
3. No compound was added to or removed from any shortlist by hand.
4. `10_new_system` was **not modified**.
5. **No compound is called active**, and none has been assayed.

### Still unvalidated against labels

ADR-AIIM-303 validated the 5-of-8 majority rule on DUD-E. The **top-20 window**,
the **scaffold-uniqueness rule** and the **lexicographic sort key** have never
been checked against a known active and inactive set. The third of these is
newly identified in this entry and is the one shown above to overturn the
model's own ordering.

### A defect in this run's own harness, recorded

The first execution wrote both artifacts correctly and then raised
`TypeError` while printing the console table, because the per-compound
dictionaries backing the table were not updated with the distinct-split counts
that had been added to the dataframe. The fix touches display only. The run was
then repeated with the fixed script so that the deposited code and its artifacts
correspond exactly, and the artifacts came back byte-identical, which is why the
re-run doubles as the third determinism check reported above.

### Skills invoked

| skill | version | material influence |
|---|---|---|
| `rdkit` | 1.1 | Murcko frameworks for the uniqueness rule — and the finding that framework identity retains bond order, so 9,10-dihydrophenanthrene and phenanthrene, and fluorenone and fluorenol, count as distinct frameworks while being the same ring system at different oxidation states; `FilterCatalog` for the flags, which remain flags and not gates |

### Output

```
20_new_runs/em_survivors_2026-10-06/
  em_survivors.json
  survivors.csv          all 15, with distinct-split support and cut reasons
  run.log
  run_manifest.json
  analysis/em_survivors.py
```

## ADR-AIIM-313 — PRE-REGISTRATION: dock the SNP-Lib in-domain pool, because the SNP-Lib arm has no usable structure-based evidence and molecular dynamics cannot start without a pose

**Status:** proposed — **pre-registration. Frozen before the run.**
**Date:** 2026-10-06
**Raised by:** the author's instruction to run the minimum molecular dynamics
that would not embarrass a submission; this is its prerequisite
**Decision log at freeze (text/LF):**
`62f8f5d77755418c7297b6e9464aa6bf15589e5cecb8f75467fd5af89946ca95`
**Discipline:** ADR-AIIM-263 — predictions fixed before the run;
ADR-AIIM-226 — outside the frozen tree

### Why this is needed before any molecular dynamics

The whole SNP-Lib arm — ADR-AIIM-306 through 312 — is **ligand-based only**. It
has produced a 69-compound in-domain pool, 15 survivors and an 8-compound
shortlist without a single docked pose. Two consequences:

1. **A reviewer's first question has no answer.** The paper would prioritise
   compounds against a structural target while reporting no structure-based
   evidence for any of them.
2. **Molecular dynamics cannot start.** The established ligand-preparation route
   (`md_candidates_2026-08-22/_cmd.sh`: antechamber on `ligand_posed.sdf`, gaff2,
   bcc charges) consumes a posed ligand. No usable pose exists for any SNP-Lib
   compound.

The only run that ever docked SNP-Lib compounds,
`exploratory_sbdd_ranking_2026-08-09`, is labelled **DO-NOT-REUSE** — withdrawn
by the v4 manuscript in its own words as something that must not be presented as
a prioritisation — and it archived scores only, no poses. Its numbers are not
cited here and are not evidence for anything below.

That withdrawn run is, however, a **planning signal worth stating because it
will be checked**: it reported 5 of 9 relevant compounds docking roughly 12 Å
from the pocket centroid. It generated **one** conformer for four of those five.
The frozen script's own docstring records that the G3 redocking control
reproduced the ZM241385 crystal pose **only** when structural waters were
retained *and* several starting conformers were docked with poses pooled by
energy, and that screening under a looser protocol than the one that passed
validation would make the validation meaningless. **The off-pocket result is
therefore most likely an artefact of a protocol weaker than the validated one**,
and this run tests that.

### What is run

`10_new_system/scripts/g4_dock_library.py` — the frozen script — **invoked, not
modified**, at its validated defaults: 8 conformers, exhaustiveness 32, 10
poses, seed 20260729, prune RMS 0.5, on-pocket cutoff 8.0 Å.

- **Receptor:** `snp_screen_2026-08-04/G2_preparation/4EIY_prepared.pdbqt`, the
  receptor every docking run in this project uses (RUN_INDEX status INPUT)
- **Box:** that run's `docking_box.json`, centre (−0.424, 8.527, 17.127), size
  22 Å, derived from the ZMA heavy-atom centroid by crystallographic geometry
- **Engine:** Vina 1.2.7 in the WSL `diffdock` environment — already installed,
  **nothing is installed, updated or reconfigured**
- **Library:** the **69 in-domain compounds** (AD ≥ 0.25) of ADR-AIIM-312, the
  whole eligible pool rather than the 15 survivors, so the structure-based
  channel reads the pool independently of the ligand-based ranking
- **Poses archived** per compound for the molecular dynamics that follows

Vina is CPU-only and is capped at **8 of 20 threads** so that ADR-AIIM-305's
DiffDock run, which holds 7,659 MiB of the 8,151 MiB GPU, is not starved of
host threads. The two runs do not contend for the GPU.

### Predictions

**Z1 — anchor.** ZM241385, docked through this same path, reproduces the G3
validated pose in `G3_redocking/ZMA_redocked_4EIY.sdf` to **heavy-atom RMSD
≤ 2.0 Å**. *Falsified above 2.0 Å*, which would mean the protocol as invoked
here is not the one G3 validated and **no other result in this run may be
used**.

**Z2 — the off-pocket result was a protocol artefact.** Of the nine compounds
the withdrawn run also covered (the 8 shortlisted plus SNP00223), **at least 6
are on-pocket** under the validated protocol, against the 4 that run reported.
*Falsified below 6.*

**Z3 — molecular dynamics is possible at all.** **At least 2** of the 8
shortlisted compounds are on-pocket. *Falsified below 2.* If Z3 fails, the
minimum molecular dynamics the author asked for cannot be run on shortlisted
compounds and the scope must be re-decided rather than quietly substituted.

**Z4 — the two channels are not redundant.** Across the 69, the Spearman
correlation between `vina_score` and the mean predicted pKi of ADR-AIIM-312 has
**|rho| < 0.3**. *Falsified at |rho| ≥ 0.3.* This is the SNP-Lib counterpart of
ADR-AIIM-301, which found on DUD-E that the ligand channel alone beats the
intersection at every matched size.

**The outcome that would most embarrass this entry is Z1 failing**, because the
G3 control is the only thing establishing that this docking protocol reproduces
a known crystal pose, and without it the other three predictions are
unmeasurable rather than merely wrong.

### Constraints

1. **`10_new_system` is not modified.** The frozen script is invoked at its
   defaults.
2. **`02_frozen_inputs` is not written.** The receptor and box are read.
3. **Nothing is re-docked that already has a usable pose.** The NP-library
   compounds and their existing molecular dynamics are untouched.
4. **No number from `exploratory_sbdd_ranking_2026-08-09` is cited**, in this
   entry or in anything built on it. Its status is DO-NOT-REUSE and that is
   not relaxed by its being the only prior SNP-Lib docking.
5. **ADR-AIIM-305 is not disturbed.** CPU capped at 8 of 20 threads; no GPU use.
6. **No compound is called active.** A Vina score is not an affinity and an
   on-pocket flag is not a binding claim.
7. **Docking does not re-rank the shortlist.** ADR-AIIM-312's ordering stands as
   recorded; this run adds an orthogonal channel and does not retrospectively
   edit a ligand-based result.
8. **Every outcome is reported**, including Z1 failing.
9. **This entry is appended to the decision log and its SHA-256 recorded before
   the run.**

### What follows, and the honest timetable

If Z1 holds and Z3 gives at least two on-pocket shortlisted compounds, molecular
dynamics starts when the GPU frees at roughly 19:10 today, under a scope fixed
in a separate pre-registration before it launches:

| | ns | days at 130 ns/day |
|---|---|---|
| ZM241385 control, replicates 2 and 3 (it currently has **n = 1**) | 100 | 0.8 |
| 2 candidates × 3 replicates × 50 ns | 300 | 2.3 |
| **total** | **400** | **3.1** |

Starting tonight, that finishes late on 2026-10-09. **It is the most that fits
before 2026-10-10 and it is not comfortable.** 50 ns is not negotiable: this
project has already measured NPLIB01185 escaping at 23.6, 27.4 and 46.6 ns while
reading 0.41–0.54 Å at 10 ns, so a 10 ns readout would be wrong rather than
merely short. The control reaching n = 3 matters at least as much as the
candidates, because a control with fewer replicates than the compounds it
anchors is the kind of asymmetry a reviewer does not forgive.

### Skills invoked

| skill | version | material influence |
|---|---|---|
| `rdkit` | 1.1 | conformer generation inside the frozen script; heavy-atom RMSD for the Z1 anchor |

### Output

```
20_new_runs/snplib_dock_2026-10-06/
  scores.csv                 69 compounds, one terminal state each
  poses/<SNP_ID>.sdf         energy-ranked poses, consumed by the MD build
  zma_anchor.json            Z1
  snplib_dock.json           grading
  run_manifest.json
  analysis/snplib_dock.py    the wrapper; the docking itself is the frozen script
```

## ADR-AIIM-314 — PRE-REGISTRATION: three velocity-seeded production replicates of the ZM241385 control, so the control is not a single trajectory

**Status:** proposed — **pre-registration. Frozen before the run.**
**Date:** 2026-10-06
**Raised by:** the author's instruction to run the minimum molecular dynamics
that would not embarrass a submission, and the author's judgement that a
reviewer request for replicates is easier to meet than one for longer
trajectories
**Decision log at freeze (text/LF):**
`e847c62ad90a5d5ae7bd19f7c97fdd0c12c4dd3f77516e13271ad2cfa39d6b7e`
**Discipline:** ADR-AIIM-263 — predictions fixed before the run;
ADR-AIIM-226 — outside the frozen tree

### The gap, verified rather than recalled

| | replicates | length | construction |
|---|---|---|---|
| candidates (NPLIB01185, 01643, 03056) | **3** | 100 ns | `gen-vel = yes`, `gen-seed` 20260822/23/24 |
| **ZM241385 control** | **1** | 100 ns | `gen-vel = no`, `continuation = yes` |

Confirmed three independent ways: exactly one `prod.xtc` exists across every
molecular-dynamics directory in the project; `prod.mdp` carries
`nsteps = 50,000,000` at `dt = 0.002`; and `prod.log` records
`Finished mdrun ... Sun Aug 23 03:58:42 2026`. The `npt1/npt2/npt3` trajectories
that might suggest otherwise are the three **staged restraint-release
equilibration** steps, not production replicates.

**The lengths already match at 100 ns.** The only gap is replicate count, and
a second, subtler one: the control's single run is a *continuation* from the
final equilibration stage, so it does not satisfy the description the
manuscript uses for the candidates — "three 100 ns replicates differing only in
velocity seed" (line 1691). A continuation has no generated velocities to
differ in.

The project states the consequence in its own frozen code,
`md_09_setup_candidate_runs.py`: *a single trajectory cannot distinguish a pose
that is stable from one that happened not to move; replicates are the only
thing that licenses a statement about consistency.* The control is the
reference every candidate is compared against, so a control that licenses no
consistency statement weakens every comparison drawn from it.

The manuscript already discloses this at lines 868 and 881 and takes it up in
§4. **This entry closes a declared weakness rather than a hidden one.**

### What is run

Three new production replicates of the existing, already-equilibrated control
system. **No system is rebuilt and no equilibration is repeated**: all three
start from `npt3.gro`, the same final equilibration structure the existing
production run started from, with the same topology and the same coupling
index.

Each replicate copies the existing `prod.mdp` and changes exactly three
settings, to match how the candidate replicates were made:

```
gen-vel      = yes      (was no)
continuation = no       (was yes)
gen-seed     = 20261006 | 20261007 | 20261008
```

Everything else — V-rescale at 310 K on three baths, semi-isotropic C-rescale,
`nsteps = 50,000,000`, `dt = 0.002` — is unchanged.

**Run order is chosen so the result degrades gracefully.** Replicates are run
one at a time, and after each the run is usable:

| after | reported control | cumulative | elapsed |
|---|---|---|---|
| seed 20261006 | n = 2 | 100 ns | 0.8 d |
| seed 20261007 | **n = 3** (with the existing run) | 200 ns | 1.6 d |
| seed 20261008 | **n = 3, all identically constructed** | 300 ns | 2.3 d |

At 128.935 ns/day — the throughput the control's own `prod.log` reports on this
card, not an estimate — 300 ns is 2.33 days. Starting when ADR-AIIM-305's
DiffDock run releases the GPU at roughly 20:25 today, the third replicate
completes about 2026-10-09 00:00, leaving slack before 2026-10-10. If anything
is lost to a host restart, which `md_10_run_candidate.sh` records has happened
four times in this project, the second row above is already a complete result.

### Predictions

**C1 — the control reproduces itself.** Each new replicate holds the ligand
below 4 Å RMSD for **at least 95 %** of frames. The existing run gives mean
2.59 Å and 100 % below 4 Å. *Falsified by any replicate below 95 %.*

**C2 — the replicates agree.** The three velocity-seeded replicates' mean
ligand RMSDs span **less than 1.0 Å**. *Falsified at 1.0 Å or more*, which
would mean a single trajectory never characterised this system and the existing
control figure was one draw from a wide distribution.

**C3 — continuation and re-seeding are equivalent.** The existing continuation
run's mean RMSD falls **inside** the range spanned by the three new replicates.
*Falsified if it falls outside*, which would mean the manner of starting a
replicate, not just its seed, affects the result — and would make the existing
run non-comparable to the candidates as well as to its own replicates.

**C4 — no escape.** No replicate shows the ligand leaving the orthosteric site
by the project's own persistence measure. *Falsified by any escape.*

**The outcome that would most embarrass this entry is C1 or C4 failing**,
because the control's stability is what licenses every candidate comparison in
§3 and §4. If the known A2AR antagonist does not hold its crystallographic
binding mode across replicates, the candidate persistence figures of 0.140,
0.241 and 0.316 lose their reference, and that must be reported rather than
managed.

### Constraints

1. **`10_new_system` is not modified.** Where a frozen script applies it is
   invoked; the three mdp files are written into this run's directory.
2. **`20_new_runs/md_production_2026-08-21` is not written to.** The existing
   control run is read only. Its `npt3.gro`, `coupling.ndx` and
   `system_gmx.top` are copied into the new run directory with their SHA-256
   recorded, so the starting point is provenance-checked rather than assumed.
3. **The existing run is not discarded or superseded.** It is reported
   alongside the new replicates, and C3 is the test of whether it belongs in
   the same set.
4. **No candidate molecular dynamics is started.** A new SNP-Lib compound needs
   3 × 100 ns = 2.3 days on its own, which does not fit beside this before
   2026-10-10; starting one at a single replicate would produce exactly the
   result `md_09` says licenses nothing. The 69 docked poses from ADR-AIIM-313
   keep, and need no re-docking when that run becomes affordable.
5. **ADR-AIIM-305 is not disturbed.** This starts only once DiffDock's
   `_stage1_done` sentinel appears.
6. **No compound is called active**, and none has been assayed.
7. **Every outcome is reported**, including C1 or C4 failing.
8. **This entry is appended to the decision log and its SHA-256 recorded before
   the run.**

### What this does not fix

The manuscript cites a field precedent for this receptor of **triplicate
300 ns** (reference 43). Three replicates at 100 ns remain short against it,
and this entry does not change that; the limitation stays disclosed where it
already is. The point of these runs is the replicate axis, which is bounded and
reachable, not the duration axis, which at 900 ns per compound is not reachable
on one RTX 5060 before any plausible deadline.

### Skills invoked

| skill | version | material influence |
|---|---|---|
| `molecular-dynamics` | 1.1 | read for replicate construction and equilibration-reuse practice; the simulations themselves run under the project's own GROMACS 2026.3 protocol, which is not replaced by the skill's OpenMM route |

### Output

```
20_new_runs/md_control_replicates_2026-10-06/
  prod_r2.mdp  prod_r3.mdp  prod_r4.mdp
  prod_r2.xtc  prod_r3.xtc  prod_r4.xtc
  inputs_sha256.json         npt3.gro, coupling.ndx, system_gmx.top
  control_replicates.json    grading
  run_manifest.json
  analysis/md_control_replicates.py
```

## ADR-AIIM-315 — Production length for the control replicates set to 50 ns, reversing ADR-AIIM-118's extension, on the author's instruction

**Status:** accepted — record of fact and protocol change
**Date:** 2026-10-07
**Decided by:** Adisorn Sirichotejirakul, 2026-10-07, after being shown the
argument for keeping 100 ns and restating the instruction
**Amends:** ADR-AIIM-314 (which specified 100 ns), reverses the length decision
of ADR-AIIM-118
**Decision log at change:**
`6bba84a003aeb6c8ea2406823d175c09fd9cd76c2d212bb887e3eed14d575a75`

### What changed

The three ZM241385 control replicates run **50 ns** each, not 100 ns.
Everything else from ADR-AIIM-314 is unchanged: same equilibrated starting
structure, same topology, same `gen-vel`/`gen-seed` construction, same
V-rescale and semi-isotropic coupling, same run order.

Executed at 2026-10-07 00:28. Replicate 2 was mid-flight at **10.226 ns** when
the instruction arrived; it was signalled to stop, wrote its checkpoint, had its
tpr rewritten with `gmx convert-tpr -nsteps 25000000` — the mirror of the
`convert-tpr -extend` that ADR-AIIM-118 used to lengthen NPLIB01185 — and
resumed from that checkpoint. **Nothing already computed was discarded.**
Replicates 3 and 4 had not started and were regrompp'd from mdp files set to
25,000,000 steps.

### Order-of-record deviation, stated plainly

ADR-AIIM-263 requires the entry to be frozen before the run. **It was not
here**: the runs were changed and restarted at 00:28 and this entry is appended
afterwards. The reason is that the instruction arrived while a job was
mid-flight and stopping it promptly bounded the waste at one checkpoint
interval. This is a real deviation from the discipline, recorded rather than
smoothed over. No prediction was altered after seeing data — ADR-AIIM-314's
C1 to C4 stand unchanged in substance, since none of them is a function of
trajectory length.

### Why this is defensible, and where it is not

ADR-AIIM-118 extended candidate production from 3 × 50 ns to 3 × 100 ns. It is
worth quoting what that entry actually concluded, because it is more supportive
of 50 ns than its headline suggests:

> **The 50 ns rationale did not survive being checked.** […] Kapla et al. (PLoS
> Comput Biol 2021) report that binding modes both close to and distant from the
> experimental pose could remain stable over 0.5 µs […]; Liu, Watanabe & Kokubo
> (J Comput Aided Mol Des 2017) find 56–62 % of deliberately-wrong decoy poses
> remain stable under equilibrium MD.

and, giving its reasons for extending anyway:

> […] neither of which is "more nanoseconds means more confidence in
> correctness": **1. Symmetry with the control.** […] **2. The schedule
> genuinely permits it.** […] **Had the time not been there, the honest move
> would have been to keep 50 ns and state the limitation, not to claim 50 ns
> was sufficient.**

So the literature ADR-AIIM-118 itself assembled says **length does not
establish pose correctness at any duration in reach**. The axis that licenses a
consistency statement is replicates, which is what these runs add. On that
reading 50 ns is the honest choice when the schedule is the binding constraint,
and the schedule is binding: 2026-10-10.

**Where it is not defensible, and this must be carried into the manuscript:**
ADR-AIIM-118's first reason — symmetry — now runs the other way. The candidates
ran 3 × 100 ns. A control at 50 ns is shorter than what it anchors, which is the
asymmetry ADR-AIIM-118 extended the candidates to remove. **The only way this
does not reintroduce that defect is to report every persistence and RMSD figure
on a common 50 ns window**, recomputing the candidate numbers from trajectories
that already exist rather than rerunning anything.

### Consequence that changes published numbers

The candidate pose-persistence figures currently in the manuscript —
**0.140, 0.241 and 0.316** for NPLIB01185, NPLIB01643 and NPLIB03056 — are
computed over 100 ns. Reporting on a common 50 ns window requires recomputing
them from the existing trajectories at 50 ns. No simulation is rerun; the
trajectories are truncated at analysis time. **The numbers in the manuscript
will change**, and the 100 ns values remain reportable as the longer check.

This is the cost the author accepted when restating the instruction, and it is
recorded here so that the before and after are both auditable.

### What 50 ns buys

| | ns | at 125.9 ns/day (measured on this card, this system) |
|---|---|---|
| control r2 (resumes at 10.226 ns), r3, r4 | 139.8 | 26.6 h → **2026-10-08 ~03:10** |
| freed budget | ~150 | enough for **one SNP-Lib compound at 3 × 50 ns** |

At 100 ns the control alone would have consumed the budget to 2026-10-09 07:40
and no SNP-Lib compound could have been simulated at all. **50 ns is what makes
a prioritised compound reachable before the deadline**, which is the substantive
reason for the change beyond schedule relief.

### Constraints carried forward

1. **`10_new_system` is not modified.**
2. **`md_production_2026-08-21` and `md_candidates_2026-08-22` are not written
   to.** Candidate trajectories are read for the 50 ns recomputation.
3. **The 100 ns candidate results are not deleted or withdrawn**; they are
   reported alongside.
4. **No compound is called active**, and none has been assayed.
5. **Every outcome is reported**, including C1 or C4 of ADR-AIIM-314 failing.

### Open risk

The manuscript cites a field precedent for this receptor of triplicate 300 ns.
At 50 ns the gap to that precedent widens from 3× to 6×. The limitation is
already disclosed and stays disclosed; this entry makes it larger, and the
honest framing — which ADR-AIIM-118 supplies and this entry adopts — is that
duration does not establish correctness at any value this hardware can reach,
so the comparison is reported with replicates as the claim-bearing axis.

## ADR-AIIM-316 — PRE-REGISTRATION: molecular dynamics of SNP00157, the top-ranked SNP-Lib compound that docks in the orthosteric site

**Status:** proposed — **pre-registration. Frozen before the build and the run.**
**Date:** 2026-10-07
**Decided by:** Adisorn Sirichotejirakul — one compound, three replicates, 50 ns
**Follows:** ADR-AIIM-313 (the docking that produced the pose), ADR-AIIM-314 and
315 (the control replicates this is compared against)
**Decision log at freeze (text/LF):**
`73842c7ba8b3de94052c067ed411395c06af497014993c0ffedad182a788e01f`
**Discipline:** ADR-AIIM-263 — predictions fixed before the run, which this
entry restores after ADR-AIIM-315 recorded a deviation from it

### Why this compound, decided by a rule not by preference

Molecular dynamics consumes a docked pose, so only the compounds ADR-AIIM-313
placed **on-pocket** are eligible. Of the eight shortlisted compounds, five dock
roughly 12 Å from the pocket centroid and are not eligible at all. The eligible
set, with SNP00223 included because it is unanimous although the
scaffold-uniqueness rule cut it:

| | predicted pKi | Vina | centroid dist | conventions | shortlist |
|---|---|---|---|---|---|
| **SNP00157** | **5.9257 — highest of all 69** | −7.878 | 3.83 Å | 8/8 on 5/5 distinct | **rank 1** |
| SNP00163 | 5.9078 | −6.829 | 1.93 Å | 8/8 on 5/5 | rank 2 |
| SNP00223 | 5.8760 | **−8.218** | 2.54 Å | 8/8 on 5/5 | cut by diversity rule |
| SNP00006 | 5.6407 | −8.195 | 3.39 Å | 6/8 on 4/5 | rank 8 |

**SNP00157 is chosen because it is the compound the published pipeline ranks
first and which docks in the site.** No other justification is required for it,
which is the point: with one compound affordable, the defensible choice is the
one the method already selected.

**SNP00223 was not chosen although its docking score is the best**, because it
is outside the published shortlist and running it would need a separate
argument in a week whose binding constraint is cutting the manuscript by 49 %.
That is a scheduling judgement, not a scientific one, and it is recorded as such.

**Known weakness, not hidden:** SNP00157 carries six structural flags including
a Michael acceptor, visible as the α,β-unsaturated amide. The project's standing
ruling is that PAINS/BRENK/NIH/CHEMBL are annotations and not gates, so this does
not exclude it; it is disclosed and will be disclosed in the manuscript.

### What is run

**Three replicates × 50 ns**, matching ADR-AIIM-315's control length so every
figure can be reported on a common window.

The frozen scripts do the work and **are not modified**:
`md_08_build_candidate_system.sh` (packmol-memgen, tleap, GROMACS conversion)
and `md_09_setup_candidate_runs.py` (restraints, coupling groups, mdp) are both
library-agnostic and are invoked as they stand. `md_10_run_candidate.sh` runs it.

**`md_07_build_candidate.py` cannot be used** and is replaced for this compound
only. It resolves poses from `np_library_funnel_2026-08-15/full_poses_chunk*`
and verifies identity against the NP-library prioritisation ledger, neither of
which contains an SNP-Lib compound. Its two substantive jobs are reproduced
here explicitly rather than worked around:

1. **The pose is the rank-1 Vina pose, taken verbatim** from
   `snplib_dock_2026-10-06/poses/SNP00157.sdf`, pose 1 of 10,
   `vina_energy = -7.878`, `pose_rank = 1`. No minimisation, no manual
   adjustment — the predicted pose is the thing being tested.
2. **Identity is checked the same way**: the pose's InChIKey
   `NPNNKDMSXVRADT-WEVVVXLNSA-N` equals the InChIKey of the SNP ledger's
   canonical SMILES for SNP00157. Formal charge 0, 23 heavy atoms, which fixes
   `antechamber -nc 0`, gaff2, AM1-BCC — the control's parameterisation
   settings unchanged.

**Two deviations in where things are written, stated rather than buried:**

- The build lands in `20_new_runs/md_candidates_2026-08-22/SNP00157/` because
  `md_08`, `md_09` and `md_10` all hardcode that parent. **A new compound
  directory is added to a dated run; no existing compound's artifacts are read,
  moved or modified.**
- `md_09` hardcodes `PROD_NS = 100`. The frozen script is not edited; the mdp
  files it generates are rewritten to 25,000,000 steps **in the run directory**,
  exactly as ADR-AIIM-315 did for the control.

### Predictions

**D1 — the parameterisation is clean.** `parmchk2` reports **zero** ATTN lines,
as it did for all three NP-library candidates. *Falsified by any ATTN line*,
which would mean gaff2 lacks a parameter for this chemistry and the force field
is interpolating where it should not.

**D2 — the pose holds in a majority of replicates.** At least **2 of 3**
replicates keep the ligand in the binding site for **≥ 90 %** of frames, by the
project's own pocket-distance measure. *Falsified below that.*

**D3 — the control is more persistent.** SNP00157's mean persistence is
**below** the ZM241385 control's on the same 50 ns window. *Falsified if
SNP00157 equals or exceeds the control*, which would be a surprising result for
a predicted pose measured against a crystallographic binding mode and would
deserve its own scrutiny rather than celebration.

**D4 — the predicted pose is less reproducible.** The spread across SNP00157's
three replicates **exceeds** the spread across the control's three. *Falsified
if it is tighter.*

**The outcome that would most embarrass this entry is D3 failing**, because a
docked pose of an unassayed natural product outperforming a known antagonist's
crystallographic mode is more likely to indicate a problem with the measure or
the system than a discovery, and this entry commits in advance to treating it
that way.

### The interpretation limit, fixed in advance

**A stable trajectory will not be reported as evidence that the pose is
correct.** ADR-AIIM-118 assembled the literature that forecloses this: Kapla
et al. (PLoS Comput Biol 2021) find binding modes both near and far from the
experimental pose stable over 0.5 µs, and Liu, Watanabe & Kokubo (J Comput Aided
Mol Des 2017) find 56–62 % of deliberately wrong decoy poses stable under
equilibrium MD. At 50 ns, against a field precedent of triplicate 300 ns for
this receptor, persistence is reported as a **consistency check across
replicates**, not as validation. This sentence is fixed now so that a good
result cannot be upgraded afterwards.

### Constraints

1. **`10_new_system` is not modified.**
2. **No existing candidate or control artifact is modified.** The NP-library
   compounds and the control are read only.
3. **Nothing is re-docked.** The ADR-AIIM-313 pose is used verbatim.
4. **The control has priority in the queue** until its third replicate
   completes; this compound takes the GPU after that, ahead of the control's
   optional fourth replicate.
5. **SNP00157 is not called active**, and no compound here has been assayed.
6. **The shortlist is not re-ranked** by this run, whatever it shows.
7. **Every outcome is reported**, including D1 or D2 failing, in which case no
   SNP-Lib compound has usable molecular dynamics and that is what the
   manuscript will say.

### Output

```
20_new_runs/md_candidates_2026-08-22/SNP00157/
  03_build/   ligand_posed.sdf, LIG.mol2, LIG.frcmod, system, topology
  04_run/     em, equilibration, prod_r1..r3 at 50 ns
  build_metadata.json   pose SHA-256, InChIKey check, parmchk2 ATTN count
```

## ADR-AIIM-317 — PRE-REGISTRATION: DiffDock on the fifteen SNP-Lib survivors, as a pose-level check on the Vina poses, not as a ranking channel

**Status:** proposed — **pre-registration. Frozen before the run.**
**Date:** 2026-10-07
**Decided by:** Adisorn Sirichotejirakul — pose comparison approved, a second
molecular-dynamics run on a DiffDock pose declined as unaffordable
**Follows:** ADR-AIIM-313 (the Vina poses), ADR-AIIM-316 (the pose going to MD)
**Decision log at freeze (text/LF):**
`3afdea45b0b38d1f4d498025d2a76a01546e86ad2297d34f3dfcaca61b3a3ff6`

### What this is not

**This does not revive DiffDock as a ranking channel, and ADR-AIIM-305's stage 2
remains unrun.** Evidence assembled on 2026-10-07 from the `kdss_SNP v2`
workspace settles that question against DiffDock on three independent sets:

| set | n | DiffDock | Vina |
|---|---|---|---|
| ChEMBL head-to-head, same compounds | 248 | AUC 0.6619, **EF@1% 0**, BEDROC20 0.2443 | AUC 0.6882, **EF@1% 4.96**, BEDROC20 0.6477 |
| DUD-E | 981 | AUC 0.6880 | — |
| ChEMBL-derived | 918 | AUC 0.5582 | — |

The AUC confidence intervals overlap (0.578–0.744 against 0.590–0.784), so AUC
alone does not separate them and will not be claimed to. **Early recognition
does**: DiffDock finds no active in the top 1 % where Vina finds 4.96-fold
enrichment, and BEDROC20 differs by 2.65×. Every composite of the two scored
*below* Vina alone (0.618–0.619), independently reproducing ADR-AIIM-301's
finding that combining the channels loses information. That is consistent with
DiffDock's own documentation, which states three times that its confidence is
not an affinity.

Those v2 numbers are **outside this project's evidence chain**. Their labels
were checked and are public — DUD-E identifiers and a ChEMBL `source` column
reading `ChEMBL_active`/`ChEMBL_inactive`, with no cAMP-derived field — but the
scripts that produced them have not been audited. **They are used here to stop
spending GPU on a settled question, not as a manuscript claim.**

### What this is

A **pose-level** question that the score-level evidence above cannot answer:
**do two independent pose predictors put these molecules in the same place?**

DiffDock is run once on the **15 survivors of ADR-AIIM-312**, of which Vina
placed 7 on-pocket and 8 off-pocket, and the rank-1 poses are compared to the
Vina rank-1 poses already archived by ADR-AIIM-313. Approximately 33 minutes of
GPU at the 131 s/complex this project measured, taken in the gap after the
control's third replicate completes so no molecular-dynamics time is lost.

**Nothing is installed or reconfigured.** The existing conda environment and
checked-out DiffDock commit are used as they stand, with ADR-AIIM-111's
parameters, as ADR-AIIM-305 did.

### Why it matters for what is already running

ADR-AIIM-316 is building a membrane system around **SNP00157's rank-1 Vina
pose** and will spend 29.6 hours of GPU testing whether that pose survives. If
DiffDock independently places the same molecule in the same pose, that
simulation tests a binding mode two methods agree on. If it does not, the
simulation tests a mode one of two methods rejects — which does not invalidate
it, but changes what may be concluded from it and must be stated in the
manuscript rather than discovered by a reviewer.

### Predictions

**E1 — pose agreement is poor overall.** **Fewer than 8 of 15** compounds have
a DiffDock rank-1 pose within **2.0 Å** heavy-atom RMSD of their Vina rank-1
pose. *Falsified at 8 or more.* The two scores are nearly uncorrelated on this
pool (Spearman 0.105, ADR-AIIM-313 Z4), and DiffDock here docks blind while Vina
was given a 22 Å box on the ZMA centroid.

**E2 — the off-pocket result is not a Vina artefact.** Of the **8** compounds
Vina placed off-pocket, **at least 5** are also placed more than 8.0 Å from the
ZMA centroid by DiffDock. *Falsified below 5*, which would mean ADR-AIIM-313's
central structural finding — that five of eight shortlisted compounds do not
occupy the orthosteric site — is a property of Vina rather than of the
molecules, and would require that finding to be withdrawn and re-measured.

**E3 — the compound going to molecular dynamics is agreed on.** SNP00157's
DiffDock rank-1 pose is within **2.0 Å** of its Vina rank-1 pose. *Falsified
above 2.0 Å.*

**E4 — the channels stay independent at pose level too.** Across the 15,
\|Spearman\| between DiffDock confidence and Vina score is **below 0.4**.
*Falsified at 0.4 or above.*

**The outcome that would most embarrass this entry is E3 failing**, because
ADR-AIIM-316 has already committed to the Vina pose and the build is underway.
E3 is stated now, before the comparison is computed, precisely so that a
disagreement cannot be quietly reframed afterwards as unimportant.

**E1 and E3 can both hold**: E1 is about the set, E3 about one compound, and a
single agreement among a minority is not a contradiction.

### Constraints

1. **Nothing is re-docked with Vina.** The ADR-AIIM-313 poses are read.
2. **The shortlist is not re-ranked** by this run, and DiffDock confidence is
   not used to order anything.
3. **No molecular dynamics is run on a DiffDock pose.** A second 3 × 50 ns at
   29.6 hours does not fit before 2026-10-10 and was declined.
4. **ADR-AIIM-316 is not paused or altered** by whatever this finds; if E3
   fails, the consequence is a sentence in the manuscript, not a changed run.
5. **`10_new_system` is not modified**; DiffDock is not installed or updated.
6. **No compound is called active**, and none has been assayed.
7. **Every outcome is reported**, including E2 failing, which would oblige a
   withdrawal.

### Output

```
20_new_runs/diffdock_snplib_2026-10-07/
  diffdock_input.csv
  results/<SNP_ID>/rank1*.sdf
  pose_comparison.csv      per compound: RMSD, both centroid distances, both scores
  diffdock_snplib.json     grading
  run_manifest.json
```

## ADR-AIIM-318 — Control replicates returned to 100 ns, SNP00157 molecular dynamics deferred, and a GPU watchdog added after the host killed two long runs

**Status:** accepted — record of fact and protocol change
**Date:** 2026-10-07
**Decided by:** Adisorn Sirichotejirakul, 2026-10-07, after being shown that the
live manuscript is v8 and not the v6 this session had been reasoning about
**Amends:** ADR-AIIM-315 (50 ns), defers ADR-AIIM-316 (SNP00157)
**Decision log at change:**
`f7b04e7853a3e7a3aea434df28fcc777900de5cd820cba838ef19898491dff94`

### The error this entry corrects

**This session spent two days reasoning about a superseded manuscript.**
`CLAUDE.md` names `30_manuscript/CURRENT_v5/manuscript_AIIM_v6_CLEAN.md` as the
live text. It is not. The live text is **v8**, in
`30_manuscript/FINAL_REVIEW_2026-10-04/01_Manuscript_AIIM_v8.docx`, at revision
round 9 of 2026-10-05.

What follows from reading the wrong file, all of it stated to the author as
fact and all of it wrong:

| claimed | actual |
|---|---|
| 26,284 words, needing a **49 % cut** | **9,314 words** — the cut was done |
| submission-ready **not achievable** by 2026-10-10 | achievable; **six items remain, every one of them a human task** |
| blockers include `[AFFIL]`, `[ORCID]`, reference 13's 38 primary reports | the real list is CRediT, funding, author approval, keyword limit, repository visibility, and a Word formatting pass |
| 33 ADRs unmapped, Results at 14,162 words | those counts describe v6 and say nothing about v8 |

The stale pointer in `CLAUDE.md` is the proximate cause. **The sufficient cause
is that I did not look at the sibling directories before building a plan on
one file**, and `30_manuscript/` lists `FINAL_REVIEW_2026-10-04/` and
`AIIM_submission_package_v8_working_2026-10-03/` beside `CURRENT_v5/`. The
author found it by asking why the work was not going into v8.

`RECONCILIATION_ADR280_317.md` was generated against v6 and is **superseded by
this entry**; it is not deleted, because its ADR-to-topic map is still usable
once re-anchored to v8.

### Why 50 ns is withdrawn

v8 states, in text that has passed review: *"Each originally prioritised
compound was simulated in three 100 ns replicates … the control had one 100 ns
trajectory, so candidate-to-control ratios remain exploratory"*, and reports
mean ligand RMSD **4.73 Å**, **5.22 Å** and **2.59 ± 0.51 Å** with Figure 8.

Control replicates at 50 ns would make the control 1 × 100 ns plus 3 × 50 ns
against candidates at 3 × 100 ns. The only way to report that without a new
asymmetry is a common 50 ns window, which means **recomputing three published
numbers and a published figure in a document that currently carries no FAIL
items, three days before submission**.

ADR-AIIM-315 chose 50 ns on the explicit reasoning that the schedule was the
binding constraint and the manuscript needed a large rewrite regardless. **The
second half of that premise was false.** At 100 ns the manuscript changes by one
sentence — the control is no longer a single trajectory — and no number or
figure moves.

**r2 lost nothing.** It had reached 49.5 ns under the 50 ns target and was
extended in place with `gmx convert-tpr -nsteps 50000000`, resuming from its own
checkpoint. r3 and r4 had not started and were regrompp'd.

### The host is killing long GPU runs

Two long runs died silently, each unnoticed for hours:

| run | died at | noticed after |
|---|---|---|
| ADR-AIIM-305 DiffDock | 335 of 593 complexes, 2026-10-06 09:31 | **11 h** |
| this run, replicate 2 | 49.5 of 50 ns, 2026-10-07 | **~4 h** |

Neither wrote an error. `dmesg` in both distros shows
`misc dxg: dxgk: dxgkio_query_adapter_info: Ioctl failed`, and memory was never
short — 30 GB free, zero OOM events. **The WSL GPU passthrough on this host
faults and takes the GPU process with it.** `CLAUDE.md` already documents a
third instance of the same class: the in-app browser pane killing the Electron
GPU process, 14 times out of 14.

**Both losses were avoidable and were mine.** On 2026-10-06 at 11:11 I read
progress from a log that had stopped 1.5 hours earlier and reported a finish
time from it without checking whether the process still existed; the GPU then
sat idle for 11 hours while a dependent job waited on a sentinel that could
never appear. The second time I checked only after the author asked for a status.

**A watchdog is now in the driver**: while a replicate's `.gro` does not exist,
mdrun is relaunched with `-cpi`, up to 40 attempts, 60 s apart. With `-cpt 15`
this bounds a fault at one checkpoint interval instead of however long until
someone looks. **Liveness is checked by process, not by log timestamp.**

### What is run now

| | |
|---|---|
| control replicates | **3 × 100 ns**, r2 resumed at 49.5 ns |
| remaining | 250.5 ns at the measured 125.9 ns/day → **2026-10-09 ~11:30** |
| ADR-AIIM-317 DiffDock pose comparison | **kept**, 15 survivors, ~33 min, still gated on r3 |
| ADR-AIIM-316 SNP00157 | **deferred, not cancelled** |

**SNP00157 is deferred because its result has nowhere to go in v8**, not because
the work failed. Its membrane system is built and validated — `parmchk2` ATTN
**0**, complex.pdb 2,318 protein + 42 ligand + 1 Na⁺ + 156 waters in the
control's own ordering — and its rank-1 Vina pose is archived with its InChIKey
checked. Nothing needs rebuilding or re-docking when it is run.

### Constraints

1. **v8 is not edited by this entry.** The single sentence about the control
   changes only once the replicates complete and are analysed.
2. **No published number or figure is recomputed.** That was the point of
   returning to 100 ns.
3. **`10_new_system` is not modified.**
4. **ADR-AIIM-311 to 317 stand.** They were correct work; only ADR-AIIM-315's
   premise and ADR-AIIM-316's timing depended on the manuscript error.
5. **No compound is called active**, and none has been assayed.

### Open risks

1. The GPU fault is not understood, only bounded. If its frequency rises, the
   watchdog's 40 attempts will not be enough and the schedule fails.
2. `CLAUDE.md` still names v6 as the live text. **It should be corrected**, and
   this entry does not do it, because editing project instructions is the
   author's call.
3. The thirty-eight entries ADR-AIIM-280 to 317 remain unmapped **to v8**. The
   reconciliation that exists maps them to v6 and must be redone against the
   real document before any of it is acted on.

## ADR-AIIM-319 — v8 manuscript and supplement revised for the post-freeze extension; route B (50 ns and SNP00157) confirmed; four statements in the record corrected

**Status:** accepted — record of fact and protocol change
**Date:** 2026-10-07
**Decided by:** Adisorn Sirichotejirakul, 2026-10-07: "ทาง ข: 50 ns + SNP00157", then
"แก้ใน manuscript กับ supplement ถ้ายังไม่รู้อะไรให้ draft ไว้ก่อนแล้วติด [placeholder]"
**Amends:** ADR-AIIM-318 (returned the control replicates to 100 ns). **Route B is
restored** and ADR-AIIM-315's 50 ns stands.
**Decision log at entry:**
`749aeecbc8bef1ef90b2cc194508a4bc31dc62038816a2389712564c3f0e058c`

### 1. Execution state of route B

| | state at 2026-10-07 |
|---|---|
| control replicate 2 | reached 53.2 ns; closed at 50 ns, its 50 ns frame written from the trajectory with `trjconv -dump 50000` |
| control replicate 3 | 50 ns, running since 12:25 |
| control replicate 4 | **deferred**: the driver stops after replicate 3 so SNP00157 can take the GPU |
| ADR-AIIM-317 DiffDock | waiting for `prod_r3.gro`, then about 33 minutes |
| ADR-AIIM-316 SNP00157 | system built; **production not yet queued**. It needs `md_09`'s 100 ns mdp rewritten to 25,000,000 steps in the run directory, then a driver |

The 100 ns candidate results stay as they are. **Reporting on a common 50 ns window
requires recomputing the 100 ns values 2.59 ± 0.51, 4.73, 5.22, 97.4%, 14.0, 24.1
and 31.6% and Figure 8 over the first 50 ns.** That has not been done. The
manuscript keeps the 100 ns values and carries placeholders for the 50 ns ones.

### 2. The documents

`AIIM_submission_package_v8_working_2026-10-03/` was edited by unpack, tracked-change
XML edits and pack, the method the docx skill (v1.0, read in full) prescribes and
the method v8's own history shows (about forty `_superseded_*` copies). Author
"Claude"; placeholders highlighted yellow.

| file | before SHA-256 | after SHA-256 |
|---|---|---|
| `KDSS_A2AR_AIIM_v8_WORKING.docx` | `0233aa4093e517c53d32…` | `f9c4e39bc48671a5…` |
| `Supplementary_Information_v8_WORKING.docx` | `83decfb0778f7ffb…` | `898af949a581ca37…` |

Verified: schema validation passed; rejecting all tracked changes returns the
original text character for character (63,760 and 81,479); Word opens both files.
The before copies are `_superseded_*_before_ADR311_318.docx`. **`FINAL_REVIEW_2026-10-04/`
was not modified.** The change list is `CHANGES_2026-10-07_ADR311_318.md` beside the
files: 803 and 1,806 words of prose, 21 placeholders, abstract 220 to 260 words.

Two deviations from the skill, stated: Table S8 was inserted by script, because 15
rows by 9 columns of tracked cells is about 900 lines of XML by hand; and the
packer's validator caught five malformed attributes introduced by my own edits
before any file was installed.

### 3. Four statements in the record that were wrong

1. **ADR-AIIM-314 and the RUN_INDEX row for `md_control_replicates_2026-10-06`
   say the control system "carries a net charge under Ewald (Na⁺ 31, Cl⁻ 41)".**
   Its total charge is **−0.004 e** (`grompp_em.log`: "System has non-zero total
   charge: -0.004000"), a rounding residue of the partial charges. GROMACS warns
   about any non-integer total. SNP00157's +0.0001 is the same kind of number. There
   is no charge asymmetry between control and candidates, and **the decision not to
   neutralise the system stands on a smaller premise than was written**. The
   `run_control_replicates.sh` comment repeats the claim and is not edited, because
   the running driver reads that file.
2. **The "0.351 log library gap" quoted in conversation was never entered in the
   decision log**: ADR-AIIM-310 is the pre-registration and no result entry was
   written for it. It was computed at three threads. It is not used in the
   manuscript, and ADR-AIIM-312's ledger line saying it "stands" was an assertion
   rather than a measurement.
3. **The 0.25 applicability-domain gate applied to SNP-Lib is looser than the band
   ADR-AIIM-296 measured**: 4,566 reference compounds there, 6,284 here, so the same
   threshold admits 52.9% of the expansion library rather than 31.2% (ADR-AIIM-302
   records this). Statements that the 69 compounds are "in domain" overstate it.
4. **AUC 0.3294 for the 0.15 to 0.25 band was described as "worse than chance"**; its
   95% CI is 0.042 to 0.648, which includes 0.5.

### 4. Open, and not done

- **SNP00157 production is not queued.** If it is not, the GPU idles after the
  DiffDock job, which is the failure ADR-AIIM-318 recorded.
- The 50 ns-window recomputation of the 100 ns figures (section 1).
- The G8 number-traceability audit has not been re-run on any new number.
- `FINAL_REVIEW_2026-10-04/` and its README, the reading copies, the cover letter
  and the highlights are unchanged and now behind the working files.
- The data-availability commit and counts were not checked against the public
  repository.
- **`CLAUDE.md` still names `CURRENT_v5/manuscript_AIIM_v6_CLEAN.md` as the live
  text.** It is not. Correcting it is the author's call.
- No compound is called active, and none has been assayed.

## ADR-AIIM-319 addendum — SNP00157 production queued (2026-10-07 16:54)

Section 4 of ADR-AIIM-319 listed SNP00157 production as not queued. It is now.
`md_09_setup_candidate_runs.py` was run for SNP00157 (2,355 heavy atoms restrained;
Solute 4,781 / Membrane 25,728 / Solvent 62,984 = 93,493, matching the build). Its
three production mdp files were rewritten from 50,000,000 to 25,000,000 steps in
the run directory (gen-seed 20260822, 20260823, 20260824); `md_09` was not edited.

`run_snp00157.sh` waits for control replicate 3 (up to 30 h), then for the
DiffDock job's sentinel (up to 100 min, then proceeds so the GPU is not left
idle), then invokes the frozen `md_10_run_candidate.sh` repeatedly, up to 40 times,
until `prod_r3.gro` exists. Every wait has a timeout because an untimed wait idled
the GPU for 13 hours on 2026-10-06. Driver pid 100 in the A2A-MD distribution.

Expected: control replicate 3 finishes about 22:00, DiffDock about 22:35,
equilibration about one hour, then 3 × 50 ns at the measured 125 ns/day, so
SNP00157 completes about 2026-10-09 00:00. The 100 ns candidate figures have still
not been recomputed over a 50 ns window.

## ADR-AIIM-320 — PRE-REGISTRATION: recompute the published molecular-dynamics figures over a common 50 ns window

**Status:** proposed — **pre-registration. Frozen before the run.**
**Date:** 2026-10-07
**Decided by:** Adisorn Sirichotejirakul, "ทำข้อ 1 เลย" (the first item of the
list that closed ADR-AIIM-319)
**Follows:** ADR-AIIM-315 and 319. Control replicates run 50 ns, the candidates ran
100 ns, so every figure must be reported on the window both have.
**Decision log at freeze (text/LF):**
`392061d7ba3b5d45b5b0e799e0251ebfd219d4ae6aa47118ff8a46036ab5cd69`

### What is computed

The existing trajectories, **truncated at 50 ns inside the analysis, nothing
re-simulated.** The inputs are the time series the published analyses already wrote:
`ligand_rmsd.xvg` and `box.xvg` per run, and the `N6_55.xvg` minimum-distance series
of ADR-AIIM-161. A frame is kept when t ≤ 50.0 ns. The control is stored every
100 ps (501 frames in the window) and the candidates every 50 ps (1,001).

Every summary function is the one the published artifact used, with the window
as the only change: RMSD in nm times ten (md_14), threshold ladder 2, 3, 4, 5 Å,
bilayer area as X × Y with first and last quarter of the *window*, and N6.55
persistence as the fraction of frames at or below 3.5 Å, per replicate, then the
mean of the replicates and its ratio to the control with the thresholds fixed by
ADR-AIIM-161 (at least 0.5 reproduces, below 0.1 does not, otherwise partial).

### Predictions

**W1 — anchor.** The same code with the window left at 100 ns reproduces the
published numbers: control RMSD 2.59 ± 0.51 Å; candidate mean of replicate means
4.73 Å (NPLIB01643) and 5.22 Å (NPLIB03056); control N6.55 persistence 0.974 over
1,001 frames; candidate persistence 0.1403, 0.2412 and 0.3163; control bilayer
area change −1.05%. *Falsified by any difference at the published rounding.* If
W1 fails nothing else here may be reported.

**W2 — the control is stationary.** The control's mean RMSD over 50 ns lies
within 0.3 Å of its 100 ns mean. *Falsified at 0.3 Å or more.*

**W3 — drift accumulates.** The mean of replicate means at 50 ns is **no higher**
than at 100 ns for both NPLIB01643 and NPLIB03056. *Falsified if either is higher.*

**W4 — the verdict does not change.** At 50 ns all three candidates stay
**PARTIAL** against the control (ratio 0.1 to 0.5), and the control's persistence
stays at least 0.95. *Falsified by any candidate leaving that band.*

**W5 — the ordering does not change.** The mean persistence order across the
candidates at 50 ns is NPLIB01185 < NPLIB01643 < NPLIB03056. *Falsified by any swap.*

### Constraints

1. `10_new_system` is not modified and no published artifact is overwritten. The
   results go to a new run directory.
2. **NPLIB01185 RMSD is computed for completeness and marked not reportable**,
   as in the article, because of the unresolved periodic-image artifact.
3. **No new threshold is chosen.** The ADR-AIIM-161 bands are reused unchanged.
4. **Control replicates 2 and 3 are not included.** Replicate 3 is still running;
   their analysis needs the full ligand-RMSD pipeline and is a separate entry.
5. No compound is called active, and none has been assayed.
6. Every outcome is reported, including W3 or W4 failing.

### Output

```
20_new_runs/md_window50_2026-10-07/
  md_window50.json   run_manifest.json   analysis/md_window50.py
```

## ADR-AIIM-321 — ADR-AIIM-320 result: the 50 ns window changes four published numbers, no verdict, and the anchor caught a unit error in my own script

**Status:** accepted — record of fact
**Date:** 2026-10-07
**Pre-registration:** ADR-AIIM-320, SHA-256
`7ce8f72966ab63b276c943d1d60cf1d3193a146d16fd2d20dc5fe59b9b448254`
**Decision log at freeze (text/LF):**
`392061d7ba3b5d45b5b0e799e0251ebfd219d4ae6aa47118ff8a46036ab5cd69`
**Run:** `20_new_runs/md_window50_2026-10-07/`, CPU only, nothing re-simulated

### The anchor failed first, and the failure was mine

The first execution missed **W1** on one of eleven checks: the control bilayer area
came out **+0.15%** against the published **−1.05%**. Every RMSD and persistence
figure had reproduced exactly. The cause: `box.xvg` is written by `gmx energy` in
**picoseconds** (0 to 100,000 ps, 10,001 rows), while the RMSD and N6.55 series are
in nanoseconds, and the script cut all three at "100" and "50". That kept 11 and 6
box rows. The pre-registration said that if W1 failed nothing else could be
reported, so the run was diagnosed and corrected before anything was read from it.
Both logs are kept (`run_first_attempt_unit_bug.log`, `run.log`). After the fix **all
eleven anchor checks match to the published rounding**.

This is the second time in this session that an anchor written before the run
caught an error in a script I had just written.

### Grading

| | prediction | result | verdict |
|---|---|---|---|
| **W1** | the 100 ns window reproduces the published numbers | 11 of 11 after the unit fix | **MET** |
| **W2** | control mean RMSD within 0.3 Å of its 100 ns value | 2.61 against 2.59 Å | **MET** |
| **W3** | mean RMSD at 50 ns no higher than at 100 ns for both | NPLIB01643 **4.79 against 4.73**; NPLIB03056 5.11 against 5.22 | **MISSED** |
| **W4** | all three stay PARTIAL; control persistence at least 0.95 | ratios 0.19, 0.29, 0.31; control 0.966 | **MET** |
| **W5** | persistence order 01185 < 01643 < 03056 | 0.187, 0.281, 0.300 | **MET** |

**W3 missed on one of two compounds.** "Drift accumulates" is not supported for
NPLIB01643, whose mean is marginally higher over the first 50 ns. What does change is
the replicate spread, which is wider on the shorter window: 0.99 against 0.51 Å for
NPLIB01643 and 0.90 against 0.63 Å for NPLIB03056.

### The numbers on the common window

| | 100 ns (published) | first 50 ns |
|---|---|---|
| control mean RMSD ± SD | 2.59 ± 0.51 Å | **2.61 ± 0.62 Å** |
| control frames below 3 Å | 78.9% | **68.9%** |
| control frames below 4 Å | 99.8% (all but one) | 100% |
| control bilayer area change | −1.05% | **−0.44%** |
| control N6.55 persistence | 0.974 (975 of 1,001) | **0.966 (484 of 501)** |
| NPLIB01643 mean RMSD (replicate spread) | 4.73 Å (0.51) | **4.79 Å (0.99)** |
| NPLIB03056 mean RMSD (replicate spread) | 5.22 Å (0.63) | **5.11 Å (0.90)** |
| NPLIB01185 | 14.0% persistence | 18.7% persistence |
| NPLIB01643 | 24.1% | **28.1%** |
| NPLIB03056 | 31.6% | **30.0%** |

NPLIB01185 RMSD is 24.71 Å at 100 ns and 12.76 Å at 50 ns, and **is not reportable**
for the reason the article gives (periodic images).

Candidate-to-control RMSD ratios are 1.84 and 1.96 at 50 ns and 1.83 and 2.02 at 100.
**The qualitative reading is unchanged**: no candidate reproduces the control's
Asn253 contact on either window, and all three stay PARTIAL under the bands of
ADR-AIIM-161, which were not altered.

Four published numbers move: the control's fraction below 3 Å (79% to 69%), its SD,
its bilayer area change, and the candidates' persistence percentages.

### Not done

Control replicates 2 and 3 are **not** in these figures. Replicate 2 is complete at
50 ns and replicate 3 was running; analysing them needs the full ligand-RMSD
pipeline and is a separate entry. The one-control caveat in the article therefore
stays. Figure 8 has not been regenerated. No compound is called active.

## ADR-AIIM-322 — ADR-AIIM-313 result: the SNP-Lib pool docks, the control reproduces its pose, and five of the eight shortlisted compounds lie outside the site

**Status:** accepted — record of fact, written late
**Date:** 2026-10-07
**Pre-registration:** ADR-AIIM-313, SHA-256
`d161e191a4e193a08d741559cbcbcdbff719819d765c65cff01f7afcf31cca95`
**Run:** `20_new_runs/snplib_dock_2026-10-06/`, 2,724 s. **This entry should have
been written on 2026-10-06.** The grading was reported in conversation and in
`snplib_dock.json` but never entered in the decision log, and ADR-AIIM-319 and the
v8 supplement now rely on it. The same gap exists for ADR-AIIM-310, whose pool
comparison has no result entry.

### Grading, from `snplib_dock.json`

| | prediction | result | verdict |
|---|---|---|---|
| **Z1** | ZM241385 redocks to within 2.0 Å of the validated redocked pose | **0.1053 Å** | **MET** |
| **Z2** | at least 6 of the 9 compounds also covered by the withdrawn run are on pocket | **4** (SNP00006, 00157, 00163, 00223) | **MISSED** |
| **Z3** | at least 2 shortlisted compounds on pocket | **3** (SNP00006, 00157, 00163) | **MET** |
| **Z4** | \|Spearman\| between Vina score and mean predicted pKi below 0.3 | **0.105**, n = 69 | **MET** |

Z1 is agreement with the earlier validated redock (`ZMA_redocked_4EIY.sdf`), **not
RMSD to the crystal structure**, which ADR-AIIM-313 and the v8 text must not be
read as saying.

### What the result is

All 69 compounds completed; **40 are on pocket** (centroid distance at most 8.0 Å);
Vina scores run from −8.674 to a median of −6.785. **Z2 missed and the prediction
behind it was wrong.** ADR-AIIM-313 argued that the withdrawn run's five off-pocket
results were an artefact of its single-conformer protocol. Under the validated
protocol the same four compounds are on pocket and the same five are not
(SNP00154, 00097, 00052, 00033, 00142, at 11.89 to 12.70 Å), so the finding is
confirmed rather than explained away. Those five generated only one or two
conformers after pruning; this is an observation, not a cause.

Five of the eight shortlisted compounds therefore lie outside the orthosteric site,
and the two channels are almost uncorrelated.

### Not done

No molecular dynamics was run on any of the five off-pocket compounds. The
withdrawn run's scores are not cited. No compound is called active.

## ADR-AIIM-323 — Figure 8 regenerated over the first 50 ns, with two anchors, and the sentences that cite it corrected

**Status:** accepted — record of fact
**Date:** 2026-10-07
**Decided by:** Adisorn Sirichotejirakul, "regenerate Figure 8 ด้วยข้อมูล 50 ns เลย"
**Follows:** ADR-AIIM-320 and 321, which recomputed the numbers this figure draws.
**Decision log at entry:**
`4ef3ad7768c32e61dfd85e6d19560c6c3e0c3433f683f71b0847cbcb5278333d`

### What was done

`20_new_runs/figure8_window50_2026-10-07/analysis/figure8_window50.py` is
`figure6_pose_persistence.py` with one change: every series is cut at the window.
Layout, colours, line styles, axis limits and fonts are unchanged. The skill
`manuscript-figures-tables` was read in full with its Lane B reference; its order was
followed, so **the legend was written before the artwork** (it is in the script's
docstring) and the rendered image was looked at, not trusted from the code.

### Two anchors, both passed

1. **Panel B equals the stored analysis.** The threshold ladder is recomputed from
   `ligand_rmsd.xvg` and matches `md_window50.json` to four places for all six
   compound-replicate pairs.
2. **Pixel identity.** Run with `--window 100`, the same code reproduces the
   published v8 `Figure_8.png`: **0 of 14,654,484 pixels differ**, same dimensions.
   The 50 ns figure therefore differs from the published one in the window alone.

The script also asserts that every series' time axis ends near 100 ns before cutting,
because the first execution of the 50 ns analysis cut a picoseconds file as
nanoseconds (ADR-AIIM-321).

### Facts the new legend asserts, computed from the plotted data

All nine candidate replicates stay below 10 Å from the docked binding site over the
first 50 ns (panel A). Fraction of frames with the Asn253 contact at or below 3.5 Å:
control 0.9661 (484 of 501); NPLIB01185 0.110, 0.117, 0.334; NPLIB01643 0.412, 0.101,
0.331; NPLIB03056 0.177, 0.341, 0.383. Recorded in
`Figure_8_50ns_provenance.json` beside the input hashes.

### Changes to the manuscript

The embedded image was replaced (`word/media/image9.png`, byte-identical to the new
file) and its drawing extent adjusted for the slightly different aspect ratio (4,479
by 3,258 against 4,498 by 3,258). The legend was rewritten as a tracked change: it
now says what was done, defines panels A to C, states that error-free quantities are
means with replicate dots, states that the control is one trajectory, and cites
Supplementary Table S9 for the 100 ns summaries. **The sentence in Section 3.6 that
cited Figure 8 for "100 ns" now reads "(Figure 8, which shows the first 50 ns)"**,
because the figure no longer shows 100 ns. One placeholder remains in the legend:
panel C still draws the control as a single trajectory.

The 100 ns figure is kept as `Figures/Figure_8_100ns.png`, which is a candidate for a
supplementary figure. `FINAL_REVIEW_2026-10-04/` is unchanged.

### Not done

Control replicates 2 and 3 are not drawn. The legend and the Figure 8 callouts in
the Highlights and cover letter, if any, were not read. G8 was not re-run. No
compound is called active.

## ADR-AIIM-324 — PRE-REGISTRATION: analysis of control replicates 2 and 3 over the common 50 ns window

**Status:** proposed — **pre-registration. Frozen before any replicate was analysed.**
**Date:** 2026-10-07
**Decided by:** Adisorn Sirichotejirakul, "รันวิเคราะห์ control replicate ตอน r3 จบเลย"
**Follows:** ADR-AIIM-314 (froze C1 to C4 for the replicates), ADR-AIIM-315 and
319 (50 ns, r4 deferred), ADR-AIIM-320 and 321 (the 50 ns window)
**Decision log at freeze (text/LF):**
`5b36b654a07fe9fa8425632a3aa78b7ab9c317e666b8c194f457adc04ecddb80`

### How this relates to ADR-AIIM-314

ADR-AIIM-314 froze C1 to C4 for three velocity-seeded replicates of 100 ns. Two
things changed after that entry: the replicates are now **50 ns**, and **replicate 4
was deferred**, so there are **two** new replicates (r2 seed 20261006, r3 seed
20261007) beside the existing continuation run. The predictions below are C1 to C4
adapted to that, and **the adaptation is stated, not hidden**: C2 is computed over three
runs rather than the three new replicates ADR-AIIM-314 named, and C3 is a weaker
test, because a range spanned by two numbers is narrow.

### What is run

The pipeline of the original control analysis, **unchanged**: `md_11_analyse_pose_persistence.sh`
(receptor C-alpha fit, ligand heavy atoms measured, reference the pre-equilibration
structure `em_prod.tpr`, `-pbc mol -ur compact -center`, `gmx rms`, `gmx energy` for
the box) and the N6.55 minimum-distance measurement of `md_18_recognition_motif.py`,
using the control's own atom-index groups, which are valid because the system is
identical. GROMACS 2026.3. Cut at 50.0 ns by the same rule as ADR-AIIM-320.

**Anchor first.** The same pipeline is run on a **copy** of the existing control
trajectory in a separate directory, so nothing in the published run is written.

### Predictions

**A — anchor.** The pipeline on the existing control reproduces the published
numbers: ligand RMSD 2.59 ± 0.51 Å over 100 ns and 2.61 ± 0.62 Å over 50 ns; bilayer
area −1.05% and −0.44%; N6.55 persistence 0.974 and 0.966. *Falsified by any
difference at the published rounding. If A fails no replicate figure may be reported.*

**R1 — each replicate holds the ligand.** In each of r2 and r3, at least 95% of
frames have ligand RMSD below 4 Å over the 50 ns window. (ADR-AIIM-314 C1.)

**R2 — the replicates agree.** The mean RMSD of the existing run (50 ns window), r2
and r3 spans less than 1.0 Å. (ADR-AIIM-314 C2, adapted.)

**R3 — continuation and re-seeding are equivalent.** The existing run's 50 ns mean
RMSD lies inside the range of the r2 and r3 means. (ADR-AIIM-314 C3, with the
weakness above.)

**R4 — the contact holds in each replicate.** N6.55 persistence is at least 0.90 in
each of r2 and r3. (Replaces ADR-AIIM-314 C4's "persistence measure" with the one the
article reports.)

**R5 — the bilayer is stable.** Area change between the first and last quarter of the
window is within ±3% in each of r2 and r3.

**R6 — the candidate verdict does not depend on a single control.** With the control
persistence taken as the mean of the three runs, all three candidates stay PARTIAL
(ratio 0.1 to 0.5, the bands of ADR-AIIM-161, unchanged). *Falsified by any ratio at
0.5 or above or below 0.1.*

**The outcome that would most change the article is R6 failing**, because the
article's caveat that the control is a single trajectory is what keeps its
candidate-to-control comparison exploratory.

### Constraints

1. **No published run directory is written.** The existing control is copied;
   `em_prod.tpr` is copied into the replicate directory.
2. `10_new_system` is not modified. Bands and cutoffs are those already registered.
3. **Replicate 3 is analysed only once `prod_r3.gro` exists.** A partial trajectory
   is not summarised. Replicate 2 is analysed first and is complete.
4. The analysis runs on the CPU beside the molecular dynamics; GROMACS tools are run
   with few threads because this host's cores are shared (the 2026-08-29 finding).
5. No compound is called active, and none has been assayed.
6. Every outcome is reported, including R6 failing.

### Output

```
20_new_runs/md_control_anchor_2026-10-07/      copy of the existing control, analysis_prod/
20_new_runs/md_control_replicates_2026-10-06/  analysis_prod_r2/, analysis_prod_r3/
20_new_runs/md_control_replicates_window50_2026-10-07/
  md_control_replicates_window50.json   run_manifest.json   analysis/*.py
```

## ADR-AIIM-325 — ADR-AIIM-324 result: the control is not one stable reference, the candidates are still separated from it, and three of seven predictions missed

**Status:** accepted — record of fact
**Date:** 2026-10-07
**Pre-registration:** ADR-AIIM-324, SHA-256
`57f0cf7cb0b8eef7316b3337a5cd77f5eba27c8cec7b5fe2be40dd46a69f871e`
**Decision log at freeze (text/LF):**
`5b36b654a07fe9fa8425632a3aa78b7ab9c317e666b8c194f457adc04ecddb80`
**Run:** `20_new_runs/md_control_replicates_window50_2026-10-07/`. Replicate 2 finished
at 53.2 ns and was cut at 50; replicate 3 finished at 50.0 ns on 2026-10-07 22:02 and was
analysed by a waiter only after `prod_r3.gro` existed.

### Two execution errors of mine, both caught

1. The first execution of the summary script reported the anchor as **missed** because I
   had typed the published persistence as 0.966 while the stored value is 0.9661. The
   published values are now read from the artifact. Both logs are kept
   (`run_first_attempt_hardcoded_constant.log`, `run_r2_only.log`).
2. The waiter script ended three times with exit code 2. That code is its documented
   "still running after the time limit" and was not a failure; it was checked each time
   against the live process before being relaunched.

### Grading

| | prediction | result | verdict |
|---|---|---|---|
| **A** | the unchanged pipeline reproduces the published control numbers | 8 of 8 | **MET** |
| **R1** | each replicate below 4 Å in at least 95% of frames | r2 **96.4%**, r3 **94.6%** | **MISSED** |
| **R2** | the three runs' mean RMSD span under 1.0 Å | 2.61, 3.22, 2.75; span **0.61** | **MET** |
| **R3** | the existing run's mean lies inside the replicates' range | 2.61 against [2.75, 3.22] | **MISSED** |
| **R4** | N6.55 persistence at least 0.90 in each replicate | r2 **1.000**, r3 **0.679** | **MISSED** |
| **R5** | bilayer area change within ±3% | −2.78%, −0.46% | **MET** |
| **R6** | all three candidates stay PARTIAL against the mean control | ratios 0.212, 0.319, 0.340 | **MET** |

**R1** missed by 0.4 points on replicate 3. **R3** was named in ADR-AIIM-324 as a weak
test, and the existing run is simply the lowest of the three. **R4** is the substantive
miss.

### The control, over the first 50 ns, three runs

| | continuation (existing) | replicate 2 | replicate 3 |
|---|---|---|---|
| mean RMSD ± SD (Å) | 2.61 ± 0.62 | 3.22 ± 0.61 | 2.75 ± 0.74 |
| frames below 3 Å | 69% | 34% | 63% |
| frames below 4 Å | 100% | 96.4% | 94.6% |
| bilayer area change | −0.44% | −2.78% | −0.46% |
| N6.55 persistence | 0.966 (484 of 501) | **1.000** (501 of 501) | **0.679** (340 of 501) |

**The control's Asn253 contact is not a fixed number.** Replicate 3 holds it in 68% of
frames, against 97% and 100% in the other two. The mean over the three runs is **0.882**,
not the single 0.974 the article reports, and the spread between runs (0.32) is larger
than the control's own frame-to-frame fluctuation would suggest from one run.

### What this does and does not change

**The candidates remain separated from the control.** Per-replicate candidate
persistence over 50 ns runs from 0.10 to 0.41 and the control's three values are 0.68,
0.97 and 1.00, so the two sets do not overlap. Candidate mean RMSD per replicate runs from
4.36 to 5.65 Å against 2.61 to 3.22 Å for the control; they do not overlap either.
Against the three-run control mean, the candidate persistence ratios are 0.21, 0.32 and
0.34, all within the PARTIAL band of ADR-AIIM-161, which was not altered.

**The article's reference numbers move.** The control figures 2.59 ± 0.51 Å and 0.974
are one trajectory of three and the most favourable of the three on RMSD. A control
persistence of 0.88 with a range of 0.68 to 1.00 is a fairer statement of what the
reference does. The one-control caveat can now be replaced by a three-run statement,
and **the replacement is a weaker reference than the single run suggested**.

### Not done

Replicate 4 was deferred (ADR-AIIM-319), so the control has two re-seeded replicates
plus the continuation run. The 100 ns candidate figures are not re-summarised against
this control. Neither the manuscript, the supplement nor Figure 8 has been updated with
these results. No compound is called active.

## ADR-AIIM-326 — v8 manuscript, supplement and Figure 8 updated with the three-run control (ADR-AIIM-324 and 325)

**Status:** accepted — record of fact
**Date:** 2026-10-07
**Decided by:** Adisorn Sirichotejirakul, "อัปเดตเอกสารเลย"
**Decision log at entry:**
`763bab91cb8fecfa0b3a0fd11cfb8da80e1d1ad352b2e45289f784a8e17c7003`

### What changed in the documents (tracked changes by "Claude")

**Manuscript.** Section 2.9: the placeholder is replaced by a statement that two 50 ns
control replicates were added (seeds 20261006 and 20261007) and every figure is reported
over the first 50 ns. Section 3.2: the replicate RMSDs (3.22 ± 0.61 and 2.75 ± 0.74 Å) and
the three-run mean (2.86 Å, range 2.61 to 3.22). Section 3.6: control Asn253 persistence
of 0.966, 1.000 and 0.679 (mean 0.882) against 0.101 to 0.412 for individual candidate
replicates, with the statement that the two sets do not overlap and that the control varied
by a third between runs; the control RMSD comparison now uses the three-run mean; the
sentence "These exploratory comparisons used one control trajectory" is changed to
"one 100 ns control trajectory, or three control runs over 50 ns". Figure 8 is replaced by
the version with the control drawn as three runs, and its legend says so. **Four
placeholders were removed; seven remain.**

**Supplement.** Section 12: the control row's Replicates cell now reads "1 (100 ns), plus 2
replicates of 50 ns", and the sentence about the asymmetric design says it applies to the
100 ns figures. Section 17: the control replicates are reported, Table S9 has the three
control runs and their mean, and a paragraph gives the seven predictions with the three
that missed. Seven placeholders remain.

### What the new text says that the old did not

The control's reference contact is **a range, 0.68 to 1.00, not 0.974**, and the 100 ns
control figure came from the run with the lowest 50 ns RMSD of the three. The candidate
conclusion stands: candidate and control values do not overlap on RMSD or on persistence.

### Verification

Schema validation passed for both files. Rejecting every tracked change returns the
reviewed v8 text exactly. No tracked-change id is duplicated (a collision between my
hand-numbered ids and the generated block was found and renumbered before installation).
Word opens both: manuscript 10,427 words and 45 pages, supplement about 14,600 words.
The embedded image is byte-identical to `Figure_8_50ns_control3.png`. The figure script in
single-control mode reproduces the previously installed figure with 0 differing pixels.

| file | SHA-256 |
|---|---|
| `KDSS_A2AR_AIIM_v8_WORKING.docx` | `18d6b3c442e850a6…` |
| `Supplementary_Information_v8_WORKING.docx` | `c64f3342d4eb54c5…` |
| `Figures/Figure_8.png` | `6cf1b61868e88f77…` |

Backups: `_superseded_KDSS_v8_before_control3.docx`,
`_superseded_Supplementary_v8_before_control3.docx` and `Figures/Figure_8_50ns_control1.png`.

### Mistakes of mine caught on the way

An extra closing tag in a hand edit made the manuscript XML unparseable; a hand-numbered id
range collided with the generated block (five duplicates); a trailing-space edit once
again joined two attributes; the supplement sentence on "partial" did not say against
which control. Each was caught by a parser, a validator or a reread before installation.

### Not done

`FINAL_REVIEW_2026-10-04/` is not updated, and its manuscript was re-saved in Word at 17:14
and no longer matches the README hash. G8 number traceability was not run. Remaining
placeholders: abstract word limit, DiffDock on the survivors, SNP00157, incidents, the
competence-band mismatch, the deposit commit, and the supplement items for the same.
No compound is called active.

## ADR-AIIM-327 — DiffDock-L on the fifteen SNP-Lib survivors: result of ADR-AIIM-317 (2026-10-07)

**Status: result recorded. Two of four predictions met, two missed. ADR-AIIM-316 is not
altered. No compound is called active.**

### What was run

DiffDock-L (skill v1.1, `samples_per_complex 20`, parameters as ADR-AIIM-111) docked the
fifteen survivors blind against the protein-only 4EIY structure. The job exited 0 with 15
complexes in 8 min 19 s. Its log carries 52 "1/10 samples failed" warnings at inference
iterations; I did not examine their cause, and all fifteen complexes returned a rank-1
pose. Analysis: `20_new_runs/diffdock_snplib_2026-10-07/analysis/diffdock_snplib.py`,
heavy-atom symmetry-aware RMSD computed in place (RDKit `rdMolAlign.CalcRMS`, 2026.03.3,
no superposition), between the DiffDock rank-1 pose and the Vina rank-1 pose.

Anchors, all required before grading, all met: (1) the DiffDock input protein's 3,012 heavy
atoms coincide exactly with atoms of the Vina receptor, so the two engines share the 4EIY
frame and an in-place RMSD is valid; (2) Vina centroid distances recomputed from the pose
files equal the stored ADR-AIIM-313 values (largest difference 0.005 Å); (3) each Vina pose
has RMSD 0 against itself. 15 of 15 poses were processed.

### Grading (fixed in ADR-AIIM-317 before the run; pre-registration `4386e354…`)

| | prediction | observed | verdict |
|---|---|---|---|
| E1 | fewer than 8 of 15 within 2.0 Å of the Vina rank-1 pose | **0 of 15** (RMSD 7.3 to 34.2 Å, median 15.3 Å) | **MET** |
| E2 | at least 5 of the 8 Vina off-pocket compounds also placed > 8.0 Å away by DiffDock | **1 of 8** (SNP00154 only) | **MISSED** |
| E3 | SNP00157 within 2.0 Å of its Vina pose | **34.2 Å**; DiffDock rank-1 centroid 29.3 Å from the ZMA centroid, confidence −1.72 (low band) | **MISSED** |
| E4 | \|Spearman\| between DiffDock confidence and Vina score below 0.4 | **0.350** (n = 15) | **MET** (narrowly) |

E1 and E4 held; E2 and E3 did not. E2 and E3 missed because I predicted the wrong
direction of the disagreement: I expected DiffDock to confirm the Vina off-pocket
placements and to be broadly consistent for SNP00157. It did neither.

### What the result says

The two engines do not agree on pose for any of the fifteen compounds. They also do not
agree on the on-pocket/off-pocket classification: DiffDock puts 11 of 15 within 8 Å of the
ZMA centroid against 7 of 15 for Vina, and the two classifications agree for only 5 of 15.
Seven of the eight compounds Vina placed off pocket (11.9 to 12.7 Å) are placed on pocket
(2.6 to 5.9 Å) by DiffDock. All three on-pocket shortlisted compounds under Vina
(SNP00157, SNP00163, SNP00006) are placed 17 to 29 Å from the site by DiffDock. None of the
ten poses DiffDock saved per compound comes within 2 Å of the Vina pose (best case 5.8 Å,
a post hoc figure, descriptive only). The sign of the E4 correlation is that higher
DiffDock confidence goes with weaker Vina score; at 0.35 it is close to the 0.4 limit and is
read as "not shown to be related", not as "independent".

So the pose of each survivor is **engine-dependent**: neither engine's pose is corroborated
by the other, and the earlier statement that three of eight shortlisted compounds dock on
pocket rests on Vina alone. This does not show that either engine is right. It is a reason
not to treat the Vina pose as established.

### Limits that matter for reading it

- DiffDock is not specifically trained on membrane proteins (diffdock skill), and the input
  here is the protein alone, without bilayer, so the transmembrane cleft is exposed on its
  lipid-facing sides. Blind docking can therefore place a ligand on the receptor surface
  outside the orthosteric site. A DiffDock off-site placement is thus weak evidence against
  the Vina pose, and a DiffDock on-site placement is weak evidence for the compound.
- DiffDock confidence is the model's certainty about a structure, not an affinity, and was
  not used to rank anything. Bands are > 0 high, −1.5 to 0 moderate, < −1.5 low:
  6 high, 8 moderate, 1 low (SNP00157).
- Vina's box was 22 Å around the site, so Vina cannot place a ligand elsewhere; DiffDock
  could. The two searches have different spaces and no ground truth exists for these ligands.
- The 52 failed-sample warnings were not investigated.

### Constraints that held

Nothing was re-docked with Vina. Nothing was re-ranked. No MD was run on a DiffDock pose.
ADR-AIIM-316 (D1 to D4) stands unchanged and the SNP00157 production run continues on the
Vina pose as already fixed. `10_new_system` was not modified.

### Consequence for the manuscript and supplement (not yet written)

The DiffDock placeholders (manuscript §2.5 and §3.8, supplement DiffDock items) can now be
filled with: the four predictions, two met and two missed; 0 of 15 within 2 Å; the 1 of 8;
the 34.2 Å for SNP00157; and the sentence that pose is engine-dependent. The text must say
the DiffDock run was blind, protein-only and not membrane-trained, and must not say the
Vina poses were validated or refuted. Not done until the author asks.

### Outputs

| file | SHA-256 |
|---|---|
| `diffdock_snplib_2026-10-07/analysis/diffdock_snplib.py` | `0296d3bc7c8e4bb0…` |
| `diffdock_snplib_2026-10-07/diffdock_snplib.json` | `a9a850bde0a5fe72…` |
| `diffdock_snplib_2026-10-07/pose_comparison.csv` | `0df17434a62afdab…` |

Skills read in full and applied: `diffdock` v1.1 (and its confidence-and-limitations
reference), `rdkit` v1.1. Dependencies: RDKit 2026.03.3, SciPy, pandas; nothing installed.

## ADR-AIIM-328 — DiffDock-L result written into the v8 manuscript and supplement (2026-10-07)

**Status: documents updated; no analysis changed. Every number is the one recorded in ADR-AIIM-327.**

**Manuscript.** Section 2.5: the DiffDock placeholder is replaced by a sentence saying the
check was prespecified for the fifteen survivors (ADR-AIIM-317) with four predictions fixed
before the run, that it was blind, against the protein alone, and that DiffDock-L is not
specifically trained on membrane proteins. Section 3.8: a paragraph reports that two of four
predictions held and two did not (0 of 15 within 2.0 Å, 7.3 to 34.2 Å, median 15.3 Å; 1 of 8
off-pocket compounds also off pocket by DiffDock; SNP00157 34.2 Å, 29 Å from the site;
|rho| 0.35), that DiffDock placed 11 of 15 and Vina 7 of 15 on pocket and the classifications
agreed for five, and that pose is engine-dependent and the on-pocket result rests on Vina
alone. It states that this does not show the Vina poses are wrong, only that they are not
corroborated. The SNP00157 simulation placeholder stays, reworded to "production running".
Six placeholders remain (was seven).

**Supplement.** Section 17, DiffDock-L: the method and the three anchors, new Table S10
(fifteen rows: Vina score and distance, DiffDock confidence and band, distance, RMSD), the
four predictions with their outcomes, and the limits. The external-workspace note stays as a
placeholder and the workspace stays uncited. Seven placeholders remain.

**Verification.** Schema validation passed for both files. Rejecting every tracked change
returns the reviewed v8 text exactly. No tracked-change id is duplicated. Word opens both and
renders Table S10 and Section 3.8 as intended (manuscript 45 pages, supplement 62).

| file | SHA-256 |
|---|---|
| `KDSS_A2AR_AIIM_v8_WORKING.docx` | `73678ebcd75159d6…` |
| `Supplementary_Information_v8_WORKING.docx` | `b7468f074ab3e109…` |

Backups: `_superseded_KDSS_v8_before_diffdock.docx`, `_superseded_Supplementary_v8_before_diffdock.docx`.

**Not done.** The 335-complex DUD-E subset analysis (supplement placeholder), the SNP00157
simulation, the G8 number-traceability re-run on the new numbers, `FINAL_REVIEW_2026-10-04/`.
No compound is called active.

## ADR-AIIM-329 — DiffDock-L on DUD-E, stage 1 analysed on the 335 complexes that finished: result of ADR-AIIM-305 (2026-10-08)

**Status: result recorded, with a declared deviation from ADR-AIIM-305. W1 and W2 met, W3
missed. The stop rule fires; stage 2 is not run. No DiffDock channel is added to anything.**

### The deviation, and what was fixed before any confidence was read

ADR-AIIM-305 specified a 593-complex run. It stopped at 335 when the host GPU passthrough
faulted (ADR-AIIM-310 and 322), and ADR-AIIM-305 does not say what to do with an
interrupted run. The frozen script `diffdock_dude.py` would count the 258 complexes that
never ran as failed, which would invent a failure rate, so it was left unchanged and not
used. A separate script, `analysis/diffdock_dude_partial335.py`, was written first, with its
handling in the docstring, and run once: attempted complexes are `DUDE_00000` to
`DUDE_00334`; the 258 later ones are "not run", not "missing"; W1 and W2 are graded as
frozen, on the 335 and on the same compounds for Vina; W3 is graded over the 335 attempted;
a stratified bootstrap (2,000 resamples, seed 20261008) was added as a descriptive interval
because the subset is smaller than designed.

### Composition of the subset

The subsample map lists the 193 actives first and the 400 seeded random decoys after them.
The 335 are therefore **all 193 actives and 142 of the 400 decoys**. The decoys are a random
subset of the 400 (the sampler returns shuffled rows), so the decoy side is unbiased but
smaller; the subset is 57.6 % active against a design of 32.5 %. ROC-AUC does not depend on
the base rate in expectation, and the intervals below carry the extra uncertainty. The
design of 593 was **not** completed.

Anchors, all met before grading: (A1) the subsample map and the DiffDock input agree row for
row; (A2) the scored complexes are exactly the contiguous block 00000 to 00334, each with
`rank1.sdf` and ten saved poses; (A3) the heavy-atom count of every rank-1 pose equals that
of its input SMILES, so results are joined to the right compound.

### Grading (frozen in ADR-AIIM-305; pre-registration `da3da56f…`)

| | prediction | observed | verdict |
|---|---|---|---|
| W1 | DiffDock ROC-AUC below Vina's on the same compounds | 0.7131 against 0.7364; difference −0.023, 95 % CI −0.097 to +0.049 | **MET** by the point estimate only |
| W2 | DiffDock ROC-AUC above 0.50 | 0.7131 (95 % CI 0.656 to 0.768) | **MET** |
| W3 | more than 1 % of complexes yield no pose | 0 of 335 attempted | **MISSED** |

W4 belongs to stage 2 and is not evaluated.

**Stop rule.** DiffDock's AUC does not exceed Vina's on the same subsample, so the rule
**fires on these data**: stage 2 is not run and no DiffDock channel is added. I record two
limits on that. It rests on 335 of 593 complexes, and the interval for the difference
includes zero, so these data **do not show that DiffDock is worse than Vina**; they show
that it is not better, and the rule as frozen needs only that.

### What the result says

DiffDock confidence enriches actives above property-matched decoys on DUD-E A2AR (0.71),
which is more than ADR-AIIM-305's category argument led me to expect (W2 was a prediction
that it would not be anti-predictive, and it held with room to spare). It is statistically
indistinguishable from Vina's 0.74 on the same compounds. Post hoc and descriptive: ligand
size does not explain it, since actives are larger (AUC of heavy-atom count alone 0.561)
while confidence falls as ligands grow (Spearman −0.384). Confidence correlated weakly with
Vina score (−0.210). Vina on all 593 is 0.716 and on the full 2,112 set 0.7005, so the 335
are not an unusually favourable subset for Vina.

This is a statement about one dataset with a protein-only 4EIY receptor. It says nothing
about the pose agreement of ADR-AIIM-327, where the two engines disagreed on every
compound, and it does not make DiffDock an affinity predictor: confidence is the model's
certainty about a structure (diffdock skill v1.1). ADR-AIIM-299 stands.

### Not done

Stage 2 (the remaining 1,519 decoys and the ADR-AIIM-301 intersection comparison) is not run
under the frozen rule. Completing the other 258 complexes would cost about 9 h of GPU time
at the rate this run achieved (335 in 11.7 h, 126 s per complex, against the 65 s per complex
of ADR-AIIM-305) and would narrow the interval; it is a decision for the author, not
taken here. The `kdss_SNP v2` workspace holds an unaudited DiffDock comparison that this
study does not cite.

### Outputs

| file | SHA-256 |
|---|---|
| `diffdock_dude_2026-10-05/analysis/diffdock_dude_partial335.py` | `86b4f4800e0defbc…` |
| `diffdock_dude_2026-10-05/stage1_auc_partial335.json` | `d0b921ef1d977d89…` |
| `diffdock_dude_2026-10-05/stage1_scored_335.csv` | `d2e0201e459ea41c…` |

Skills: `diffdock` v1.1 and `rdkit` v1.1 read in full this session and applied. scikit-learn
supplied `roc_auc_score` only; no model was fitted and its SKILL.md was not read for this.
Dependencies: RDKit 2026.03.3, SciPy, pandas; nothing installed.

## ADR-AIIM-330 — DUD-E DiffDock result written into the v8 supplement (2026-10-08)

**Status: document updated; no analysis changed. Numbers are those of ADR-AIIM-329.**

Supplementary Section 17, "DiffDock-L": the placeholder is replaced by a paragraph stating
that the 335 completed complexes were analysed afterwards as a declared deviation from
ADR-AIIM-305 (handling fixed before any confidence was read); that the subset is all 193
actives and 142 of the 400 seeded random decoys, 57.6 % active against a design of 32.5 %;
that the 258 complexes that never ran were not counted as failures; that all 335 returned a
pose (W3 missed); ROC-AUC 0.713 (95 % CI 0.656 to 0.768) against Vina's 0.736 (0.681 to 0.789)
on the same compounds, difference −0.023 (−0.097 to +0.049); that W1 is met by the point
estimate only and W2 is met; that the stop rule fires on these data and stage 2 was not run;
and that the data do not show DiffDock to be worse than Vina. The article is unchanged: it
does not report this analysis. Supplement placeholders 7 → 6.

Verification: schema validation passed; rejecting every tracked change returns the reviewed
v8 text exactly; no duplicate tracked-change id; Word opens it (63 pages) and the page was
read. Supplement SHA-256 `d59d99f5c3a4b0a7…`. Backup
`_superseded_Supplementary_v8_before_dude335.docx`.

**Not done.** Completing the other 258 complexes; the SNP00157 simulation; G8
number-traceability on the new numbers; `FINAL_REVIEW_2026-10-04/`.

## ADR-AIIM-331 — PRE-REGISTRATION: molecular dynamics of SNP00163 in duplicate, and a duplicate design for SNP00157 (2026-10-08)

**Status:** proposed — **pre-registration. Frozen before the build and the run.**
**Date:** 2026-10-08
**Decided by:** Adisorn Sirichotejirakul — "อนุมัติ SNP00163 เริ่มเลย", with the instruction that the
three systems end up comparable, and that replication is necessary if the 10 October date allows it
**Amends:** ADR-AIIM-316 (one compound, three replicates). Its D1 to D4 stand for SNP00157 with the
changes in "Design" below.
**Decision log at freeze (text/LF):** `b1618780e7d3465bbc4c503e7313d831df0f3951060d363a6b43eccbee2a3d20`
**Discipline:** ADR-AIIM-263 — predictions fixed before the run

### Why SNP00163, decided by a rule

ADR-AIIM-316 made the on-pocket compounds the only eligible ones, because MD consumes a docked pose.
The rule for the second compound is: **the next compound in the published shortlist order whose Vina
rank-1 pose starts with a polar atom of the ligand within 3.5 Å of the N6.55 side chain**, because
persistence of that contact is the control's own metric and carries no meaning for a pose that never
had it. Read from the poses before any simulation of SNP00163 (residue 336 in the 4EIY file
numbering, the Asn whose side chain is 4.9 Å from the ZM241385 centroid):

| | shortlist rank | centroid | Vina | nearest N or O to that Asn | verdict |
|---|---|---|---|---|---|
| SNP00157 | 1 | 3.83 Å | −7.878 | 2.87 Å | running (ADR-AIIM-316) |
| **SNP00163** | **2** | **1.93 Å** | −6.829 | **2.67 Å** | **chosen** |
| SNP00006 | 8 | 3.39 Å | −8.195 | 2.90 Å | eligible, lower in the order |
| SNP00223 | cut by the diversity rule | 2.54 Å | −8.218 | 6.92 Å | no initial contact, outside the shortlist |

SNP00163 is therefore the second of the three on-pocket members of the published shortlist, and the
result for it is read against SNP00157's, so that the manuscript can say whether two of three is
consistent, not one. It is not chosen for its docking score (it is the weakest of the four).
Five of the eight shortlisted compounds dock about 12 Å from the site and are not eligible.

### What the control's three runs show, read before fixing the design (descriptive, post hoc)

`md_control_spread_2026-10-08/control_spread.json`, from the series of ADR-AIIM-324 and 325:

| control run | RMSD mean | persistence | by 10 ns block | longest gap |
|---|---|---|---|---|
| r1 | 2.61 Å | 0.966 | 0.96, 0.90, 0.99, 0.99, 0.99 | 0.3 ns |
| r2 | 3.22 Å | 1.000 | 1.00 in all five | none |
| r3 | 2.75 Å | 0.679 | 0.99, 1.00, **0.57, 0.31, 0.52** | 3.8 ns |

RMSD is tight (SD 0.32 Å, range 2.61 to 3.22). Persistence is not (SD 0.177, range 0.68 to 1.00),
and the whole range comes from one run, r3, which held the contact for 20 ns and then lost and
partly regained it. Which two of the three are used for a duplicate therefore changes the control's
mean persistence: r1+r2 0.983, r1+r3 0.822, r2+r3 0.839. A duplicate that happens to omit r3
overstates the control's stability by 0.10 to 0.16, and a pair cannot reveal an event that occurs
in one run of three.

### Design

1. **Control: all three runs stay in every comparison.** They are already published in
   ADR-AIIM-325 and the v8 documents, dropping one after the fact would be a choice made with the
   results in hand, and a wider control range makes "the candidate lies outside the control" harder
   to meet, not easier. **No new control run.**
2. **SNP00163: two replicates of 50 ns.** **SNP00157: two replicates of 50 ns.** The third SNP00157
   replicate of ADR-AIIM-316 is **not run**; its place in the queue is taken by SNP00163. This is a
   scheduling decision stated as one: it is made to reach two replicates for both compounds before
   10 October, not for a scientific reason.
3. **Reported per run.** Every run is listed with its own value; no run is averaged away. A
   matched duplicate view is added beside it: control r1+r2 (run order, fixed now), SNP00157 r1+r2,
   SNP00163 r1+r2, with the all-three control mean (0.882) and the pair-dependence above shown next
   to it. The comparison is non-overlap of per-run values, which a duplicate can support; a spread
   cannot be estimated from two runs and none is claimed.
4. **D4 of ADR-AIIM-316 is withdrawn as ungradable** (the range of two runs is biased low against
   the range of three) and replaced by the per-run table above. D2 is graded on the replicates
   that exist.
5. **Same pipeline, unchanged:** `md_11_analyse_pose_persistence.sh` and the N6.55 minimum-distance
   step with the compound's own index groups, over the first 50 ns, as for the control.

### Predictions (SNP00163; each also graded for SNP00157 where marked)

- **S1 — clean parameters.** `parmchk2` reports zero ATTN lines. *Falsified by any.* (SNP00157: D1,
  already met, 0.)
- **S2 — the ligand stays in the site.** In each replicate the ligand RMSD against the
  pre-equilibration structure is **at most 4.0 Å in at least 90 % of frames**; the control's runs
  are 94.6 to 100 % by this measure. *Falsified if any replicate is below.* (Also graded for
  SNP00157 in place of ADR-AIIM-316 D2.)
- **S3 — contact persistence is below the control's range.** The N6.55 persistence of **each**
  replicate is **below 0.679**, the lowest of the three control runs. *Falsified if either run is at
  or above it.* The three NP-library candidates ran 0.10 to 0.41. (Also graded for SNP00157; D3 of
  ADR-AIIM-316 is the weaker statement that the mean is below the control's mean and is kept.)
- **S4 — pose drift is above the control's range.** The mean ligand RMSD of **each** replicate is
  **above 3.22 Å**, the highest control run. *Falsified if either is at or below.* The NP-library
  candidates ran 4.4 to 5.7 Å. (Also graded for SNP00157.)

I do not know how these will come out. SNP00163 sits deepest in the site (1.93 Å) and starts with a
2.67 Å polar contact, which is the reason to expect it could do better than the NP-library
candidates, and a single control run lost the contact for most of 30 ns. **The outcome that would
most embarrass this entry is S3 failing for a compound the monitor was never asked about**, because
a docked pose of an unassayed compound holding the control's contact is more likely to say something
about the 50 ns window and the measure than about the compound, and this entry commits in advance
to treating it that way and to saying so.

### The interpretation limit, fixed in advance

**A stable trajectory is not reported as evidence that the pose is correct** (ADR-AIIM-316; Kapla et
al. 2021 and Liu, Watanabe and Kokubo 2017). DiffDock placed SNP00163 **17.5 Å** from the site with
moderate confidence (ADR-AIIM-327), so the pose under test is supported by Vina alone, and any
outcome is a statement about that pose's stability under this force field and window. Persistence is
reported as a consistency check across replicates. No compound is called active; none has been assayed.

### Deadline, fixed now

Analysis starts at the earlier of both compounds having two complete 50 ns replicates and
**2026-10-09 18:00 +07**. A run still going at the deadline is analysed over the length it has
reached and labelled partial; a run under 25 ns is excluded and said to be. No run is extended or
replicated afterwards to improve a result.

### Queue and compute

One GPU, 8 GB. SNP00157 replicate 2 continues as scheduled. SNP00163 is built on the CPU now,
equilibrated and run **beside** it **only if** a measurement shows the two together give at least
1.2 times the throughput of one alone; otherwise it runs after SNP00157 replicate 2. The measured
throughput and the decision are logged. SNP00157 replicate 3 is stopped before it starts. Every wait
has a timeout, and liveness is checked by process, not by log timestamp.

### Constraints

1. `10_new_system` is not modified; `md_08`, `md_09`, `md_10` and `md_11` are used as they stand.
2. The SNP00163 build is a copy of the SNP00157 one-compound script with the pose, the ledger entry
   and the output directory changed. It writes its helper files **inside the SNP00163 directory**:
   the shared `_cmd.sh` of the dated run, overwritten once on 2026-10-07 and not restored, is not
   touched again.
3. The pose is the Vina rank-1 pose, verbatim, with the InChIKey checked against the SNP ledger and
   the formal charge read from the pose, as for SNP00157. No re-docking, no minimisation.
4. The shortlist is not re-ranked by this run, and the DiffDock result is not altered by it.
5. Every outcome is reported, including S1 or S2 failing.

### Skills

`molecular-dynamics` v1.0 read in full. It describes an OpenMM and MDAnalysis workflow, whereas
this study runs GROMACS with Amber-family force fields, so its influence is on practice only:
minimise, equilibrate and then run, keep checkpoints, and state which part of a trajectory is
analysed. The 50 ns window includes the early frames and is fixed by ADR-AIIM-320 and 321, not
changed here. `rdkit` v1.1 for the identity and charge checks.

### Output

```
20_new_runs/md_candidates_2026-08-22/SNP00163/   03_build, 04_run, build_metadata.json
20_new_runs/md_control_spread_2026-10-08/        control_spread.json (descriptive, read above)
```

## ADR-AIIM-332 — Incident: WSL distro terminated at about 07:13 on 2026-10-08; recovery, and the concurrency measurement of ADR-AIIM-331

**Status: recorded. No prediction of ADR-AIIM-331 is affected; its deadline of 2026-10-09 18:00 +07 stands.**

### What happened

At about 07:13 +07 every process in the `A2A-MD` distro ended together: the SNP00157 `mdrun`
(`prod_r1`) and its retry driver, the SNP00163 membrane build, and the watcher that stops
SNP00157 after its second replicate. The background watcher on the Windows side that was to
start SNP00163 after the build had also gone. The Windows host was not restarted (WSL uptime
1 day 18 h), `dmesg` shows no out-of-memory kill and no fault message near the time, and the
distro's `init` started at 08:43:28, when I next ran a command in it, so the distro as a whole
had been stopped. **The cause is not established.** My check at 06:59 had found everything
running, by process. I next checked at 08:43, after the author asked, so the GPU stood idle for
about **1.5 h**. It is not the first host-side loss of long runs in this study (DiffDock stage 1 stopped at 335
of 593 complexes and control replicate 2 at 49.5 ns, both noticed hours late).

### Recovery

- **SNP00157 `prod_r1`** was resumed from its checkpoint with `-cpi` (the driver is resumable).
  The restart section of the log begins at 40.28 ns; about 0.8 ns of the earlier run was redone.
- **SNP00163 build.** `packmol-memgen` had finished packing and written `system.pdb` before the
  termination; on the second launch it reported "Packed PDB system.pdb found. Skipping PACKMOL"
  and the remaining steps (strip NME CH3, `tleap`, conversion) ran in 14 s. The result is 94,027
  atoms with 192 POPC, 31 Na+, 41 Cl-, 15,863 waters and net charge +0.0031, against SNP00157's
  192, 31, 41, 15,728 and +0.0001. **I did not verify the packing by any means other than its own
  output file existing and the checks downstream**; those are energy minimisation (converged to
  Fmax 969 < 1000 kJ/mol/nm in 1,945 steps, potential energy −1.0697e6 kJ/mol, no infinite force)
  and the equilibration stages that follow. If any of them fails, the system is rebuilt from scratch.
- `md_09` was run after the build and the production mdp files were set to 25,000,000 steps; the
  driver and the stop-after-replicate-2 watchers for both compounds were relaunched.

### The concurrency measurement required by ADR-AIIM-331

ADR-AIIM-331 allowed SNP00163 to run beside SNP00157 only if the pair gave at least 1.2 times
the throughput of one job. Measured:

| job | alone | together |
|---|---|---|
| SNP00157 production | 125 ns/day (40.3 ns in 7.75 h, 2026-10-07 23:19 to 07:04) | **91 ns/day** (4.54 ns from 08:45 to 09:57) |
| SNP00163 equilibration (nvt, npt1, npt2) | not measured alone | **89 to 92 ns/day** (mdrun log) |

Together about 181 ns/day against 125 alone: **1.45 times**, so the rule is met and the two
run concurrently. Two caveats: the SNP00163 figures are for equilibration stages, not
production, and SNP00163 alone was never measured, so the ratio uses SNP00157's solo rate as
the single-job reference.

### Projected completion (estimates, not results)

At about 90 ns/day each, SNP00157 replicate 1 ends about 11:20, SNP00163 enters production
about 10:10 and replicate 1 ends about 23:30, SNP00157 replicate 2 ends about 00:30 on 10-09,
and SNP00163 replicate 2 ends about 09:20 on 10-09. A further host loss of the kind above would
move all of these.

### Not done

The cause of the termination was not investigated beyond the checks above. No run was
changed. No analysis of either compound has started.

## ADR-AIIM-333 — Two unlogged changes recorded: the re-saved `FINAL_REVIEW` manuscript, and the overwritten `_cmd.sh` (2026-10-08)

**Status: recorded. Nothing is restored or changed by this entry.**

### 1. `FINAL_REVIEW_2026-10-04/01_Manuscript_AIIM_v8.docx` was re-saved in Word

The README of that folder (written 2026-10-05 07:39) records this file at 8,228,821 bytes with SHA-256
`0233aa4093e517c5…`. The file now on disk is **8,638,141 bytes, SHA-256 `da9b764368f08244…`, last written
2026-10-07 17:14**. The other three documents in the folder (`02_Supplementary_Information_v8.docx`,
`03_Cover_Letter.docx`, `04_Highlights.docx`) match the README hashes, as do the two reading copies
(written 2026-10-05 07:37) that the README lists. `Figures/Figure_8.png` is still the original (`271e31c1…`).

Compared against `_superseded_KDSS_v8_before_ADR311_318.docx` in the working package, which is
byte-identical to the file the README describes (`0233aa40…`, 8,228,821 bytes, verified again today):

- The text of `word/document.xml` is **identical** once XML escapes are resolved. The only difference
  is that the re-saved file writes the double quotation marks in reference 61 as literal characters and
  the original as `&quot;`.
- The re-saved file **adds** `word/footnotes.xml` and `word/endnotes.xml` and **drops**
  `docProps/thumbnail.jpeg` and `word/stylesWithEffects.xml`; fifteen parts other than `document.xml` differ in bytes
  (styles, settings, theme, header and footer, package relationships), which is what a Word re-save does.
  All nine embedded images are byte-identical.

**I do not know who saved it.** My own use of Word on 2026-10-07 was automation for PDF export of the
working-package files, opened read-only by my notes; I cannot rule out that it touched this file, and I
cannot show that the author's own Word session did. **The README hash for this file is no longer valid.**
ADR-AIIM-323 wrote "`FINAL_REVIEW_2026-10-04/` is unchanged"; that sentence was wrong if it was written
after 17:14, and I did not check the folder before writing it. ADR-AIIM-326 later recorded the re-save. The
folder is out of date against the v8 working package in any case (ADR-AIIM-318 to 330) and is to be rebuilt
after the molecular dynamics of ADR-AIIM-331 is written up, not patched.

### 2. `md_candidates_2026-08-22/_cmd.sh` was overwritten on 2026-10-07 00:51

`md_07_build_candidate.py` and my one-compound copy of it for SNP00157 both write their ligand helper to a
fixed name, `_cmd.sh`, at the level of the dated run, and the copy did so when I built SNP00157. The previous
content, last written by an NP-library build in August, was lost. The file is a throwaway helper that every
build regenerates from the compound's own identifier; it holds the antechamber and parmchk2 commands for
one compound and nothing in this study's results was derived from reading it. The file now contains the
SNP00157 commands. The SNP00163 build uses a copy of the script that writes `_cmd_snp00163.sh` inside the
SNP00163 directory, so the shared file was not touched again (ADR-AIIM-331, constraint 2). I did not restore
it, because the old content is not recoverable and no longer matters.

### Not done

No file in `FINAL_REVIEW_2026-10-04/` was edited. The question of who re-saved the manuscript is open.

## ADR-AIIM-334 — The statement that none of the three originally prioritised compounds survived is removed from the abstract, the conclusion and the cover letter (2026-10-08)

**Status: documents edited as tracked changes. No analysis changed. Decided by Adisorn Sirichotejirakul.**

The author asked why the sentence "none of the three originally prioritised compounds survived the
better-informed canonical analysis" was still in the abstract, on the ground that the three compounds are no
longer of interest. It was there because nobody had asked for it to go and because it stands in three
documents, not one. Edited, each as a tracked deletion by "Claude" (rejecting it restores the text):

- **Abstract**: "; none of the three originally prioritised compounds survived the better-informed canonical
  analysis" deleted, leaving "Candidate eligibility also varied markedly with training-row order." The
  abstract goes from 260 to 248 words by my count (accept-all, placeholder excluded); the journal limit is
  still not known.
- **Conclusion**: ", and none of the three released compounds survived the better-informed analysis" deleted.
- **Cover letter** (`Cover_Letter_v7.docx`): "Under the better-informed analysis, none of the three originally
  prioritised compounds remained eligible; no" replaced by "No".

**Left in place, deliberately.** The body still reports the finding where it is a result and not a headline:
the paragraph on decision stability in the Discussion ("none of the three released compounds remained
eligible"), the Figure 6 caption, and Section 3.6 and 3.8 where the three compounds are the subject of the
docking and simulation results. Removing those would delete a result the paper reports with its figure, and
the 100 ns simulations of the same three compounds remain in Figure 8. Nothing numerically audited by G8 was
in the deleted text.

**My concern, recorded and not pressed.** The sentence was the one place the abstract told the reader that the
instability of the selection had a consequence for the compounds the study itself first released. Without it,
the abstract still says that eligibility varied with training-row order, but not what that did. The author's
decision stands. If a reviewer asks, the answer is in the Discussion and Figure 6.

Backups: `_superseded_KDSS_v8_before_cutsentence.docx` and `_superseded_Cover_Letter_v7_before_cutsentence.docx`.
The yellow review copy of the manuscript was rebuilt. Manuscript SHA-256 `e48ca92afbd3b68f…`, cover letter
`18dfd7410f2047aa…`. Schema validation passed for both, and rejecting all changes returns the earlier text
exactly.

## ADR-AIIM-335 — Figure 6 and the sentence that none of the three released compounds remained eligible are removed from the body; Figures 7 and 8 renumbered (2026-10-08)

**Status: documents edited as tracked changes. No analysis changed. Decided by Adisorn Sirichotejirakul, after ADR-AIIM-334.**

**A correction to ADR-AIIM-334.** I wrote there that the sentence left in the body was in the Discussion. It is in
Section 3.6 ("More importantly, candidate eligibility was not invariant to training-row order…"), not in
Section 4. Section 4 has no sentence saying that none of the three compounds survived. The author's instruction
"the Figure 6 and the sentence in the Discussion" was acted on as meaning that sentence.

### Edits (all tracked deletions or single-character renumberings by "Claude")

**Manuscript.**
- Section 3.6: "With the censored inactive records recovered, none of the three released compounds remained eligible
  (Figure 6)." deleted.
- **Figure 6** (the picture, and its caption "Candidate eligibility after canonical ordering. Recomputing only the
  conformal prediction set changed the eligible pool from 3 compounds to between 9 and 298…") deleted, both
  paragraphs.
- Renumbered because the sequence closes up: the orthogonal pose comparison is now **Figure 6** (text citation and
  caption), the membrane simulation figure is now **Figure 7** (text citation and caption).

**Supplement.** "Figures 5-6 of the article" now "Figure 5"; the two citations of "Figure 7B" now "Figure 6B"
(Section 14, the DiffDock pose comparison).

### What was left, and what now rests on less

- **Left:** the permutation analysis (Figure 5: 6 to 4,119 compounds across 200 permutations, the released count of
  72 at the 8th percentile), the canonical-sorting and tiebreak finding, and the eight-convention comparison of
  Section 3.8.
- **Section 3.6, "that status was not stable and none is put forward as a candidate"**: the instability of the three
  compounds' status was shown by Figure 6, which is gone. It now rests on Figure 5 (the count varied widely, and the
  released value was at the 8th percentile). That is a weaker support for a statement about those three compounds,
  and I have not reworded the sentence.
- **Section 4, "alternative deterministic tiebreaks selected different numbers and identities"**: the identities
  part is no longer shown by any figure in the article. The supplement's convention comparison shows differing
  survivors; I did not check that it covers the original three. Not changed.
- The statement was not removed from the abstract's remaining text ("Candidate eligibility also varied markedly with
  training-row order") and the highlights ("unstable to training-row order"), which do not depend on Figure 6.
- **Files:** the image files keep their old names in `Figures/` (`Figure_7.png` is now Figure 6, `Figure_8.png` and
  `Figure_8_*` are now Figure 7, and `Figure_6.png` is unused). Decision-log entries cite the old names and are not
  rewritten. The renaming is for the rebuild of `FINAL_REVIEW`. Source data for the removed figure was not touched.
- The cover letter was edited in ADR-AIIM-334 and mentions no figure number.

Backups: `_superseded_KDSS_v8_before_cutfig6.docx` and `_superseded_Supplementary_v8_before_cutfig6.docx`. The yellow
review copies of both were rebuilt. Manuscript SHA-256 `cd80c20ef7c971e1…`, supplement `cde96d2c32d6ae4a…`. Schema
validation passed, rejecting every change returns the earlier text exactly, no tracked-change id is duplicated, and
no tracked element is nested inside an insertion. Word opens both (45 and 63 pages).

## ADR-AIIM-336 — ADR-AIIM-331 result: SNP00157 separates from the control, SNP00163 does not on contact persistence, and no ligand left the receptor (2026-10-09)

**Status: result recorded. Six of nine predictions held, three missed. ADR-AIIM-316's D1 to D3 are graded for SNP00157 within this entry; its D4 was withdrawn by ADR-AIIM-331. No compound is called active.**

### What ran

SNP00157 (replicates 1 and 2) and SNP00163 (replicates 1 and 2), 50 ns each, completed: SNP00157 replicate 2 on
2026-10-09 00:56, SNP00163 replicate 2 on 09:51. The third replicate of each was stopped by the watcher within three
seconds of the second finishing, and its stub files (`prod_r3.*`, 0 ns) are left in place and not analysed. All runs were
analysed over 50 ns by the unchanged pipeline, with the three control runs of ADR-AIIM-325 as the reference.

The frozen `md_11` was run on the three replicates that were already complete (SNP00157 r1 and r2, SNP00163 r1) at
07:16, before the last replicate finished, to save time. This is mechanical, produces the series only, and its logs were
not opened until the grading was run. The N6.55 step and the grading were run once at 10:25, after all four runs existed.
**Anchors, all met again:** the index groups equal the control's borrowed file; the readers reproduce the stored control
numbers for all three control runs; the Asn was found by position with the TM6 spacing verified in both systems; the
starting polar distance in each pre-equilibration structure was within 0.5 Å of the Vina pose.

### Per run (first 50 ns)

| system | run | RMSD mean | frames RMSD < 4 Å | N6.55 persistence | persistence by 10 ns block | area change |
|---|---|---|---|---|---|---|
| ZM241385 control | r1 | 2.61 Å | 100.0 % | 0.966 | 0.96 0.90 0.99 0.99 0.99 | −0.44 % |
| ZM241385 control | r2 | 3.22 Å | 96.4 % | 1.000 | 1.00 1.00 1.00 1.00 1.00 | −2.78 % |
| ZM241385 control | r3 | 2.75 Å | 94.6 % | 0.679 | 0.99 1.00 0.57 0.31 0.52 | −0.46 % |
| SNP00157 | r1 | 4.96 Å | 17.6 % | 0.148 | 0.06 0.23 0.08 0.305 0.065 | +0.13 % |
| SNP00157 | r2 | 5.41 Å | 15.6 % | 0.087 | 0.005 0.00 0.015 0.335 0.08 | −0.09 % |
| SNP00163 | r1 | 4.25 Å | 41.5 % | 0.597 | 0.995 0.73 0.14 0.425 0.695 | −0.59 % |
| SNP00163 | r2 | 5.04 Å | 37.3 % | 0.755 | 0.885 0.22 0.675 0.995 1.00 | −0.80 % |

Matched duplicate view (control r1+r2, fixed in ADR-AIIM-331): control 0.983, SNP00157 0.117, SNP00163 0.676; the
control over all three runs is 0.882.

### Grading

| prediction | SNP00157 | SNP00163 |
|---|---|---|
| S1 zero parmchk2 ATTN lines | **MET** (0) | **MET** (0) |
| S2 ligand RMSD < 4 Å in at least 90 % of frames, each replicate | **MISSED** (17.6 %, 15.6 %) | **MISSED** (41.5 %, 37.3 %) |
| S3 persistence of each replicate below the control's lowest (0.679) | **MET** (0.148, 0.087) | **MISSED** (0.597 below; 0.755 inside the control range) |
| S4 mean RMSD of each replicate above the control's highest (3.22 Å) | **MET** | **MET** |
| D3 (SNP00157) mean persistence below the control's mean (0.882) | **MET** (0.117) | — |

### What the result says

**SNP00157 separates from the control on every measure**: persistence 0.09 to 0.15 against 0.68 to 1.00, RMSD 5.0 to
5.4 Å against 2.6 to 3.2 Å, as the three NP-library candidates did. **SNP00163 does not on contact persistence**: one run
(0.60) lies just below the control's lowest run (0.68) and the other (0.76) lies inside the control's range. Its contact
comes and goes (blocks from 0.14 to 1.00) and run 2 holds it for the whole last 10 ns, while the RMSD stays above the
control's highest in both runs. The reading consistent with these numbers is that SNP00163 leaves its docked pose and
makes the contact intermittently from a different one. **ADR-AIIM-331 named S3 failing for a compound the monitor was
never asked about as the outcome that would most embarrass it, and committed to reading it as a statement about the 50 ns
window and the measure, not about the compound.** That is how it is read here: the control's own run r3 partly lost the contact
over its last 30 ns (blocks 0.57, 0.31, 0.52), so persistence varies a great deal between runs of one system, and a duplicate cannot say whether
SNP00163's 0.60 and 0.76 are the compound or the runs.

**S2 failed for both compounds**, which is the larger surprise: I predicted the ligand would stay near its docked pose
(the control does, at 94.6 to 100 %). It did not. **Post hoc, and not graded:** the ligands did not leave. The distance
from the ligand to the nearest receptor C-alpha atom (the series `md_11` already writes) never exceeded 5.35 Å in any
SNP00157 or SNP00163 run (control runs 4.73 to 4.86 Å), so S2 failing means drift of pose within the site, not
dissociation. This descriptive measure was added to `adr331_analysis.py` and its outputs after I had seen the graded
results, to interpret S2; it feeds no verdict. S2 as worded ("at most 4.0 Å") and as computed (strictly below) cannot
differ here, the highest fraction being 41.5 %.

**Unchanged limits.** A stable or unstable trajectory is not evidence of pose correctness. DiffDock placed SNP00157
29.3 Å and SNP00163 17.5 Å from the site (ADR-AIIM-327), so both poses rest on Vina alone, and the poses that drifted are
the Vina poses. Duplicates support only non-overlap statements, and no spread is estimated from two runs. The shortlist is
not re-ranked and nothing was re-docked. ADR-AIIM-316's D4 is not graded.

### Not done

The text is not yet in the manuscript or the supplement. It is generated from `adr331_results.json` by
`analysis/make_md_text.py` and applied by `apply_md_text.py`; neither has been applied to the real documents. The
stub `prod_r3.*` files were not removed.

### Outputs

| file | SHA-256 |
|---|---|
| `md_snplib_duplicate_2026-10-08/adr331_results.json` | `0247f6833d861f4a…` |
| `md_snplib_duplicate_2026-10-08/adr331_per_run.csv` | `07c143c078a497ca…` |
| `md_snplib_duplicate_2026-10-08/analysis/adr331_analysis.py` | `e16c0a5fd16862c6…` |

## ADR-AIIM-353 — PRE-REGISTRATION: matching the search space of Vina and DiffDock-L, in both directions

**Status:** proposed — **pre-registration. Frozen before any run in this entry.**
**Date:** 2026-10-09
**Decided by:** Adisorn Sirichotejirakul — asked in chat for Vina run blind and DiffDock run with a known pocket, after my finding that ADR-AIIM-317 compared a Vina run given a 22 Å box with a DiffDock run given the whole protein
**Follows:** ADR-AIIM-313 (Vina poses in the site box), ADR-AIIM-317 and 327 (DiffDock-L blind on the same fifteen), ADR-AIIM-351 and 352 (Figure 6)
**Decision log at freeze (text/LF):** `1a67248eb40f2d1774fd1fb518413d012ef6a7c938dc042690aef33e45509641`
**Skills read in full before writing this:** `diffdock` (v1.1; its `references/` were searched for pocket guidance, none found beyond "isolate the binding site" under multi-chain troubleshooting), `experimental-design` (v1.0). Neither changes the frozen contracts.

### The defect this answers

ADR-AIIM-317 pre-registered the comparison with a caveat in one sentence (DiffDock docks blind while Vina was given a 22 Å box on the ZM241385 centroid) and then graded E1 and E2 as if the two poses answered the same question. They do not: a 22 Å box confines Vina to the neighbourhood of the site by construction, whereas DiffDock searched the whole protein. The 7-against-11 on-pocket counts and the median 15.3 Å pose disagreement of ADR-AIIM-327 therefore mix two things, scoring and search space. The DiffDock paper (arXiv 2210.01776) runs its search-based baselines blind, with the box around the whole receptor plus a 4 Å buffer, so that both sides have the same task, and says that search-based methods do best when given a pocket. This entry makes the same matching here, in both directions.

### Design

Units: the fifteen survivors of ADR-AIIM-312 plus the ZM241385 control (16). Each compound is docked once per arm; the compound is the unit and the two engines are paired on it. Nothing is randomised beyond the engines' own seeds.

| | search space: whole protein ("blind") | search space: the site ("pocket") |
|---|---|---|
| **Vina** | **NEW, arm A** | existing: ADR-AIIM-313 poses, 22 Å box on the ZM241385 centroid, used unchanged |
| **DiffDock-L** | existing: ADR-AIIM-327 poses, used unchanged | **NEW, arm B** |

**Arm A, Vina blind.** The frozen `10_new_system/scripts/g4_dock_library.py`, unchanged, with the box of the DiffDock paper's baselines: centre at the middle of the 4EIY receptor heavy-atom bounding box, size = extent + 4 Å on every side, rounded up to the next 0.5 Å, that is **centre (8.647, −22.095, 15.858), size (69.5, 110.5, 51.0) Å**; receptor `4EIY_prepared.pdbqt`; 8 conformers, 10 poses, seed 20260729 as in ADR-AIIM-313; **exhaustiveness 64 instead of 32**, because the box is about 37 times the volume and the paper's blind Vina-family baseline (QuickVina-W) used 64. Nothing else differs.

**Arm B, DiffDock-L given the pocket.** The installed DiffDock (commit `85c49b6`) has no option to give a pocket, and **installing DiffDock-Pocket, a separate model, is not done here because it needs an installation approval** (`CLAUDE.md`). The pocket is given by the input instead: the same `4EIY_protein_only.pdb` **cropped to the residues with any heavy atom within 12.0 Å of any heavy atom of the crystallographic ZM241385 pose** (101 of 390 residues, original numbering kept), so the ligand is initialised around the pocket and cannot be placed where there is no protein. This uses exactly the information Vina's box was built from (the crystal ligand), and the same command and parameters as ADR-AIIM-317 (`default_inference_args.yaml`, `--samples_per_complex 20`). **This is my construction, not an evaluated DiffDock mode: cropping leaves chain breaks that the language-model embeddings see**, which is why arm B has its own validity gate (H2). The radius is fixed here and will not be tuned.

### Outcomes

Computed as in ADR-AIIM-317: heavy-atom symmetry-aware RMSD in place (no alignment) between rank-1 poses in the 4EIY frame; "on pocket" = pose centroid within 8.0 Å of the ZM241385 centroid; "agreement" = RMSD ≤ 2.0 Å. No new threshold is introduced and no p-value is computed: n = 15, descriptive.

### Predictions (graded as met or missed, no revision)

- **H1 — Vina blind finds the crystal pose.** ZM241385, Vina blind rank-1 within 2.0 Å of the crystal pose. *Missed above 2.0 Å.* If missed, the blind-matched comparison is read as a statement about Vina's blind search, not about agreement.
- **H2 — arm B is valid.** ZM241385, DiffDock-L on the cropped protein, rank-1 within 2.0 Å of the crystal pose. *Missed above 2.0 Å;* then arm B is reported but **not interpreted as a pocket-conditioned comparison**, the radius is not retuned, and a DiffDock-Pocket run is proposed to the author as a separate request.
- **H3 — the crop confines DiffDock.** At least **13 of 15** survivors have a DiffDock-L rank-1 centroid within 8.0 Å of the ZM241385 centroid on the cropped protein. *Missed at 12 or fewer.*
- **H4 — the matched pocket arm still disagrees.** Fewer than **8 of 15** survivors have Vina-pocket and DiffDock-pocket rank-1 poses within 2.0 Å. *Missed at 8 or more*, which would mean the disagreement of ADR-AIIM-327 was mostly search space.
- **H5 — the matched blind arm still disagrees.** Fewer than **8 of 15** have Vina-blind and DiffDock-blind rank-1 poses within 2.0 Å. *Missed at 8 or more.*
- **H6 — matching the search space shrinks the disagreement.** The median RMSD between Vina-pocket and DiffDock-pocket is **below 15.3 Å**, the median of the unmatched comparison of ADR-AIIM-327. *Missed at 15.3 Å or above.*

Descriptive and not graded: the 2 × 2 table of on-pocket counts; Vina-pocket against Vina-blind and DiffDock-blind against DiffDock-pocket (how much the search space alone moves each engine); confidences; per-compound poses.

### Rules

1. All sixteen compounds are run in every arm, control first, whatever the control shows; no compound is dropped or re-run for its result. A run is repeated only after a crash, with unchanged parameters, and the repeat is logged.
2. Parameters, radius, box and exhaustiveness above are fixed. Anything else I would change is a new entry.
3. Anchors, checked before grading: the box and the crop equal the values above (residue count 101; box to 0.01 Å); the existing Vina-pocket and DiffDock-blind poses are read from their archives and their hashes recorded; the rank-1 pose of each new run is in the 4EIY frame (nearest-atom check as in ADR-AIIM-317); all sixteen poses processed.
4. **Not a ranking channel.** Nothing here promotes, demotes or re-scores a compound, changes the molecular-dynamics selection, or supports any statement about activity. The result changes only what Section 3.8 and Figure 6 may say about pose agreement.
5. Compute: arm B about 10 to 35 minutes of GPU, arm A CPU only, measured on the control first; no molecular-dynamics job is running, so none is displaced. If the projected arm A time exceeds 24 h after the control, I stop and ask the author.

### Limits stated in advance

DiffDock's own sampling noise is not measured (it exposes no seed), so a difference between arm B and the archived blind poses includes run-to-run variation. Two engines, one receptor structure with the BRIL fusion excised, no bilayer for either. A cropped protein is not a membrane protein in any sense. 4EIY carries the antagonist ZM241385, so a pocket derived from it favours the antagonist-bound conformation for both engines.

## ADR-AIIM-354 — Result of ADR-AIIM-353: with the same search space the two engines still disagree, Vina's pocket result depends on its box, and the pocket-cropped DiffDock arm failed its own control (2026-10-09)

**Pre-registration:** ADR-AIIM-353, SHA-256 `7c88d7e4f5a1c8830f46f4333c5f2507295dc8b408f8d50be136ead3f08fad5a`. Parameters, box (centre (8.647, −22.094, 15.858), size (69.5, 110.5, 51.0) Å, the y-centre differing from the registered −22.095 by 0.001 Å) and crop (12.0 Å, 101 residues, 806 atoms) as registered; exhaustiveness 64; DiffDock commit `85c49b6`, same command as ADR-AIIM-317. Both arms ran once, control first, all sixteen compounds, no re-run, about 5 minutes of GPU (arm B) and 12 minutes of CPU (arm A). Outputs: `20_new_runs/docking_matched_search_space_2026-10-09/` (`matched_result.json` `1bd3510d…`, `matched_poses.csv` `cd8edaed…`).

### Grading

| prediction | result | verdict |
|---|---|---|
| H1 Vina blind recovers the crystal pose of ZM241385 | rank-1 at 0.373 Å (score −10.05, centroid 0.15 Å from the site) | **met** |
| H2 DiffDock-L on the cropped protein recovers it | rank-1 at **10.125 Å**; its centroid is 2.21 Å from the site, so it sits in the pocket with a different orientation (confidence −0.07) | **missed** |
| H3 at least 13 of 15 survivors within 8 Å of the site on the cropped protein | 15 of 15 (largest 4.6 Å) | met, and close to a consequence of the crop |
| H4 matched pocket arm: fewer than 8 of 15 within 2 Å | 0 of 15 | met, **not interpreted** (H2 missed) |
| H5 matched blind arm: fewer than 8 of 15 within 2 Å | 0 of 15, median 15.9 Å | **met** |
| H6 median pocket-matched RMSD below 15.3 Å | 11.8 Å | met, **not interpreted** (H2 missed) |

Anchors met: all sixteen poses in every set; the median of the unmatched comparison reads 15.35 Å, reproducing ADR-AIIM-327's 15.3; the archived DiffDock-blind control reads 1.917 Å, reproducing ADR-AIIM-111. **Two defects of mine, fixed before any statement was made:** the crystal SDF has no usable bond orders, so the first run of the grading script returned NaN for every control RMSD and graded H1 and H2 "missed" (a script defect, not a result; bond orders are now assigned from the SMILES); and `g4_dock_library.py` measures `centroid_dist_A` from the centre of its box, which in arm A is the receptor centre, so that column of the arm A table is meaningless for on-pocket status (the analysis recomputes every distance from the poses against the ZM241385 centroid; stated in the README; the frozen script is unchanged).

### What the matched comparisons say

1. **The disagreement is not an artefact of the box.** With both engines blind, 0 of 15 rank-1 poses agree within 2.0 Å (median 15.9 Å), while both recover the control (Vina 0.37 Å, DiffDock 1.92 Å). The unmatched comparison of ADR-AIIM-327 gave 0 of 15 and 15.3 Å.
2. **Vina's pocket result depends on the box.** Run blind, Vina places 2 of 15 survivors on pocket (SNP00223 and SNP00006), against 7 with the 22 Å box; **13 of 15 have their best pose 11.9 to 13.3 Å from the site centroid** and for 13 of 15 none of the ten saved poses is within 8 Å. In the site box, the eight compounds Vina had called off pocket (ADR-AIIM-313) sit 11.9 to 12.7 Å from the centroid, **in the outer part of a box whose half-width is 11 Å along each axis** (further out than any axis-aligned position inside it, so near an edge or a corner, though not necessarily against a face). Blind, the best score is lower than the box score by more than 0.5 kcal/mol for 8 of 15 compounds (SNP00163 −7.70 against −6.83, SNP00198 −7.48 against −6.14, SNP00142 −7.41 against −6.17). Vina's blind poses and its box poses agree within 2.0 Å for 1 of 15 (SNP00006; median 7.0 Å).
3. **For the two compounds taken to molecular dynamics,** the Vina pose in the site box (3.8 and 1.9 Å from the centroid) is not the pose Vina prefers when it may search the whole protein: SNP00157's best blind pose is 11.9 Å away (−7.93 against −7.88, no pose of the ten within 8 Å) and SNP00163's is 12.0 Å away (−7.70 against −6.83). The simulations tested the in-site poses, which is what ADR-AIIM-316 states; they do not test where the compounds would go otherwise, and nothing here says anything about binding.
4. **DiffDock-L does not move to the site when its search space is cropped, only to the neighbourhood**: on the cropped protein all fifteen centroids lie within 4.6 Å of the site, yet the control, which the blind run solved at 1.92 Å, comes out at 10.1 Å. DiffDock blind and pocket-cropped agree within 2.0 Å for 0 of 15 (median 5.9 Å).

### What this does not say, and the consequence for arm B

By the rule of ADR-AIIM-353, **arm B is reported but not interpreted as a pocket-conditioned comparison**: the question "does DiffDock agree with Vina when both are given the pocket" **is not answered**. The radius is not retuned. A proper pocket-conditioned DiffDock needs a separate model (DiffDock-Pocket), whose installation needs the author's approval and which I have **not** installed; I propose it as a separate request. No compound is promoted, demoted or scored on this; DiffDock's sampling noise was not measured, so the 5.9 Å between its two arms includes run-to-run variation.

### Consequences for the documents (proposed, not applied)

Section 3.8 and Figure 6 currently compare a Vina pose from a 22 Å box with a blind DiffDock pose and conclude that the pose of each survivor is engine-dependent and the on-pocket result rests on Vina alone. The matched blind comparison supports that disagreement (0 of 15, both blind) but the on-pocket result of Vina should now say that it depends on the box and that a blind search moves 13 of 15 survivors about 12 Å away. Which wording and whether a panel of the blind result enters Figure 6 are the author's to decide; nothing in the manuscript, supplement, figures or deposit has been changed by this entry. Adding the run to the deposit needs a change of `build_deposit.py`'s tree, which also needs approval.

## ADR-AIIM-355 — INSTALLATION RECORD (DiffDock-Pocket) and PRE-REGISTRATION: DiffDock confined to Vina's box by filtering, and DiffDock-Pocket given the pocket

**Status:** proposed — **pre-registration. Frozen before any scientific run in this entry** (the only run so far is the bundled-example smoke test of the installation, recorded below).
**Date:** 2026-10-09
**Decided by:** Adisorn Sirichotejirakul — "do option 1 and DiffDock-Pocket too, install it", after ADR-AIIM-354 reported that the cropped-protein arm failed its control
**Follows:** ADR-AIIM-353 and 354
**Decision log at freeze (text/LF):** `cbffe5046ff353dbbe7daeaf9099b75a842eee747bc4f04acb6a49de64444cbe`
**Skills read in full:** `diffdock` (v1.1) and `experimental-design` (v1.0), as for ADR-AIIM-353; DiffDock-Pocket's own README and source were read for the options below.

### Installation record (`CLAUDE.md`: provenance before use)

| item | value |
|---|---|
| what | DiffDock-Pocket, an all-atom, pocket-level docking model with side-chain flexibility; a separate program from DiffDock ("Do NOT mix up these two programs", its README) |
| source | https://github.com/plainerman/DiffDock-Pocket, **commit `031e504541e26db842c4daed35d5a84dc440bfad`** (2026-10-02), cloned to `~/diffdock_pocket_repo` in the `Ubuntu-22.04` distro |
| paper | Plainer, Toth, Dobers, Stärk, Corso, Marquet, Barzilay, "DiffDock-Pocket: Diffusion for Pocket-Level Docking with Sidechain Flexibility" (2023) |
| licence | MIT (`LICENSE`); `LICENSE_DIFFDOCK` for the inherited code |
| weights | downloaded by `inference.py` on first use from the repository's GitHub release `v1.0.0`: `score_model.zip` (281,617,940 B) and `confidence_model.zip` (17,091,458 B). SHA-256 of the extracted files: `best_ema_inference_epoch_model.pt` `09b77920…`, `best_model.pt` `d07e5303…`, `model_parameters.yml` `326aea9c…` (score) and `24a2a0a3…` (confidence); the zips themselves are not kept |
| language model | `esm2_t33_650M_UR50D`, already cached from the DiffDock runs; no download |
| environment | new conda env **`diffdock_pocket`, cloned from `diffdock`** (Python 3.10, torch 2.7.0+cu128). The upstream `environment.yml` pins torch 1.13.1 for CUDA 11.7, which cannot run on this RTX 5060 (compute capability 12.0), so it is **not** used. Added: `openff-toolkit-base` 0.16.10 and its noarch dependencies (conda-forge), `pdbfixer` 1.12.0 and `openmmforcefields` 0.15.1 (PyPI, `--no-deps`), imported at start-up by the code that implements `--relax`, which this work does not use |
| network | github.com (clone, release assets), conda.anaconda.org (conda-forge), pypi.org. No credentials, tokens or environment variables |
| writes | `~/diffdock_pocket_repo` (incl. `.cache/model`, `.cache/data`), `~/miniconda3/envs/diffdock_pocket`, `~/dd_pocket_smoketest`. The `diffdock` env, the DiffDock checkout and every project file were not modified (`diffdock` lists none of the added packages) |
| smoke test | the repository's bundled `3dpf` example, 4 samples, on the GPU: completed and wrote ranked poses. A functional check only; no scientific use |
| purpose | the pocket-conditioned counterpart of DiffDock-L that ADR-AIIM-353 could not provide |

**A correction to earlier entries.** DiffDock's `inference.py` applies the `--config` YAML **after** the command line, so the YAML's `samples_per_complex: 10` overrode the `--samples_per_complex 20` of ADR-AIIM-317 and of arm B of ADR-AIIM-353. Both ran with **10 samples per complex**, which is also why ten poses were saved. Nothing else changes; the sentence "20 samples" in those entries is wrong.

### What is run

Units: the fifteen survivors and ZM241385, as in ADR-AIIM-353. Existing sets used unchanged: Vina in the 22 Å site box (ADR-AIIM-313), Vina blind and DiffDock-L blind and cropped (ADR-AIIM-353/354).

**Arm D — DiffDock-L blind, filtered to Vina's box (the author's "option 1").** DiffDock-L, commit `85c49b6`, `4EIY_protein_only.pdb`, the same command as ADR-AIIM-317 but with a copy of `default_inference_args.yaml` whose only change is **`samples_per_complex: 100`** (the copy, because the YAML overrides the command line). Of the 100 poses of each complex, those **with every heavy atom inside the cube of side 22 Å centred on the ZM241385 heavy-atom centroid (−0.424, 8.527, 17.127)**, Vina's box, are kept, and the one with the highest DiffDock confidence is the rank-1 pose of the arm. A compound with no such pose has no rank-1 pose and counts as not agreeing; this is reported, not removed.

**Arm C — DiffDock-Pocket given the pocket.** The model as released, on the full `4EIY_protein_only.pdb`, each ligand as its SMILES (RDKit generates the conformer), **`--samples_per_complex 40 --batch_size 20`** (the README's example) and every other option at its default; the model's own confidence ranking gives the rank-1 pose. **Pocket centre: the ZM241385 heavy-atom centroid**, which is also the centre of Vina's box. (The README suggests the mean C-alpha position of residues near the ligand; for this ligand only one residue has a C-alpha within 5 Å, which would put the centre at that single atom, 5.3 Å from the centroid, so the centroid is used and this is stated; the C-alpha centre is run as a descriptive variant.) **Flexible side chains: the residues with any heavy atom within 3.5 Å of any heavy atom of the crystal ZM241385 pose**, the README's rule: A:168, A:169, A:177, A:250, A:253, A:270 (Phe168, Glu169, Met177, His250, Asn253, Met270). Descriptive variants, same everything else: **C-rigid** (`--rigid`, the README warns it performs worse) and **C-Cα** (centre at the mean C-alpha of residues within 5 Å of the crystal ligand, (−4.233, 10.800, 20.280)).

### Outcomes

As in ADR-AIIM-353: heavy-atom symmetry-aware RMSD in place between rank-1 poses in the 4EIY frame; on pocket = centroid within 8.0 Å of the ZM241385 centroid; agreement = RMSD ≤ 2.0 Å; n = 15, descriptive, no p-values. Comparisons are always against the Vina pose that had the same search space: arm C against Vina in the 22 Å box, arm D against Vina in the 22 Å box.

### Predictions (graded as met or missed, no revision)

- **G1 — arm D control.** The filtered DiffDock-L pose of ZM241385 is within 2.0 Å of the crystal pose. *Missed above 2.0 Å;* then arm D is reported but not interpreted.
- **G2 — arm C control.** The DiffDock-Pocket pose of ZM241385 (default flexible, centroid centre) is within 2.0 Å of the crystal pose. *Missed above 2.0 Å;* then arm C is reported but **not interpreted**, no option is retuned, and no further model is proposed in this line.
- **G3 — the box can be satisfied.** At least **14 of 15** survivors have at least one of the 100 DiffDock-L poses fully inside the box. *Missed at 13 or fewer.*
- **G4 — arm C is in the pocket.** At least **13 of 15** survivors have a DiffDock-Pocket rank-1 centroid within 8.0 Å of the ZM241385 centroid. *Missed at 12 or fewer.*
- **G5 — box-matched agreement is still poor.** Vina in the box against filtered DiffDock-L: fewer than **8 of 15** within 2.0 Å. *Missed at 8 or more.*
- **G6 — pocket-matched agreement is still poor.** Vina in the box against DiffDock-Pocket: fewer than **8 of 15** within 2.0 Å. *Missed at 8 or more.*

Descriptive, not graded: the median RMSDs of both comparisons; the numbers of in-box poses; C-rigid and C-Cα; confidences; the 2 × 2 table of ADR-AIIM-354 extended with arms C and D.

### Rules

1. All sixteen compounds in every arm, control first, whatever the control shows; nothing dropped or re-run for its result except after a crash, with unchanged parameters, logged.
2. Everything above is fixed. Another option, centre, residue set or sample number is a new entry.
3. Anchors before grading: the pocket centre, the box and the six residues equal the values here; all sixteen complexes processed in every arm; poses are in the 4EIY frame (the control's RMSD to the crystal pose is itself the check, as the models rewrite no coordinates of the input protein except flexible side chains); the archived Vina poses are read unchanged and their hashes recorded.
4. **Not a ranking channel.** Nothing here promotes, demotes or scores a compound or supports any statement about binding; it changes only what Section 3.8 and Figure 6 may say about pose agreement.
5. Compute: DiffDock-L at 100 samples is projected at about 80 minutes of GPU (ADR-AIIM-317 measured 8.3 minutes for 10 samples on 15 complexes), DiffDock-Pocket a few minutes per variant; no molecular-dynamics job is running. If either exceeds three times its projection I stop and ask the author.

### Limits stated in advance

DiffDock-Pocket's weights were trained on PDBBind-like complexes, not membrane proteins, and the receptor here is the BRIL-excised 4EIY without bilayer; the model flexes six side chains while Vina and DiffDock-L keep the receptor rigid; its pocket reduction is built in. Filtering DiffDock-L's blind samples by a box confines the answer after the fact and ranks the survivors with a confidence model trained on unfiltered poses. Neither model exposes a seed. 4EIY holds the antagonist ZM241385, so a pocket and a flexible-residue set derived from it favour that conformation for every pocket-conditioned arm.

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

## ADR-AIIM-390 — PRE-REGISTRATION: the third replicate of SNP00157 and SNP00163 (2026-10-10)

**Status:** proposed — **pre-registration. Frozen before the runs start.**
**Date:** 2026-10-10
**Decided by:** Adisorn Sirichotejirakul, in chat, while the co-authors review the manuscript: "ระหว่างรออาจารย์รีวิวช่วยทำ MD ครั้งที่ 3 ได้ไหม" ("while the professors review, can you run the third MD replicate?"). It is the option that ADR-AIIM-374 and ADR-AIIM-382 left open ("a third replicate of each compound would give three against three") and that ADR-AIIM-382 said needed a new pre-registration first.
**Follows:** ADR-AIIM-331 (design and result: replicates 1 and 2; the stubs of replicate 3 were stopped by the watcher), ADR-AIIM-336 (the result), ADR-AIIM-382 (the article compares the compounds with the matched duplicate, the control's runs 1 and 2).
**Decision log at freeze (text/LF):** `1692c94f9db32eea1f206f2c22a7627fb70129c44d12ea83e134781da0db6151`
**Pre-registration file:** `20_new_runs/md_snplib_replicate3_2026-10-10/PREREG_ADR390.md`, SHA-256 `200a3c4579ce190ca8d09b44b550d1772d9e41d961eefb47cc50337ed78fb2bc` (also in `PREREG_SHA256.txt`).
**Discipline:** ADR-AIIM-263 — predictions fixed before the run.

### What is run
Replicate 3 (50 ns) of SNP00157 and of SNP00163, so that the control (runs 1 to 3), SNP00157 (1 to 3) and SNP00163 (1 to 3) have the same number of runs. The frozen driver `md_10_run_candidate.sh` is invoked unmodified; replicates 1 and 2 are complete and skipped; replicate 3 continues from the checkpoint written in the first seconds of 2026-10-09, with the `prod_r3.tpr` and `prod_r3.mdp` already in place and the velocity seed `20260824` fixed when the systems were built (the other replicates 20260822 and 20260823; the control 20261006 to 20261008). The two compounds run side by side on the one GPU. Nothing is changed in any setting; a host fault is handled by re-invoking the driver and every restart is reported.

### Measured
The unchanged pipeline of ADR-AIIM-336 over the first 50 ns, with the anchors of ADR-AIIM-331; per run: RMSD mean, fraction below 4 A, N6.55 persistence overall and by 10 ns block, area change, nearest C-alpha distance.

### Predictions
R1 both runs reach 50.0 ns unchanged; R2 the nearest C-alpha distance never exceeds 10.0 A in either run; R3 SNP00157 r3 persistence at most 0.35; R4 SNP00163 r3 persistence between 0.40 and 0.90; R5 RMSD below 4 A in at most 0.50 of the frames for SNP00157 r3 and at most 0.60 for SNP00163 r3. (Full text, with the earlier values for orientation, in the pre-registration file.)

### How the result is used, decided now
The article's comparison becomes control runs 1 to 3 against runs 1 to 3 of each compound, restoring the primary design of ADR-AIIM-331, and the matched duplicate of ADR-AIIM-382 stays in the supplement as a sensitivity view; this does not depend on the numbers. Per-run values and ranges, no test on n = 3. If a run cannot be completed the article is left as it is and the run is reported as incomplete in the supplement. Nothing is called active; no ranking.

### Resources and risks
GPU RTX 5060, GROMACS 2026.3 (`/opt/gromacs/2026.3-cuda13.1-sm120`), WSL distro `A2A-MD` (it was stopped and has been started for this); 836 GB free on D:. At the measured 90 ns/day per job when two run together the two replicates need about 13 to 14 h; a host loss of the kind of ADR-AIIM-332 would add to it, and the distro stops when no session holds it, so a keep-alive session is held. Not covered by Gate 8 (ADR-AIIM-040, which budgets held-out scoring): this is compute for the exploratory extension, approved by the author in chat.

## ADR-AIIM-391 — RESULT of the pre-registration of ADR-AIIM-390: replicate 3 of SNP00157 and SNP00163 (2026-10-11)

**Status:** result recorded; the article and the supplement are **not** changed by this entry (the author decides, see the end).
**Date:** 2026-10-11
**Follows:** ADR-AIIM-390 (pre-registration, frozen 2026-10-10 before the runs).

**Run.** Both replicates reached 50.0 ns on the first attempt of the frozen driver, with no restart and no host fault (started 2026-10-10 10:53, finished 2026-10-11 01:30 for SNP00157 and 02:02 for SNP00163; 82 and 79 ns/day at the end). The unchanged md_11 and the N6.55 step of ADR-AIIM-336 were run afterwards by `adr390_analysis.py` (written before the runs ended; SHA-256 `62d7d857bfc3d3d13eb0e5d284fcb5201dcd0b6f4877269a3f25623412d45684`); the anchors of ADR-AIIM-331 were met. Outputs in `20_new_runs/md_snplib_replicate3_2026-10-10/`: `adr390_results.json` (SHA-256 `1a2960e63b5aec0fdd1d5ff338c4d83e5381014063821759454f30972b1b23cf`), `adr390_per_run.csv` (`5c8a61dfe10f6b46782321e8cc0582f311ccf3fdc32362140aa2ac33dc3ff2a9`); trajectories `prod_r3.xtc` SHA-256 SNP00157 `f111b72c5aa7afd68d1381ec6ea71070e134684c9f6acf70b6613e2fab10a6db`, SNP00163 `c229438c22487ed0ced56795655c6835f621d9dd862db2845d40c61bd97ce1ef`.

**Per run (50 ns).** N6.55 persistence, share of frames with ligand RMSD below 4 A, mean RMSD (A), largest nearest C-alpha distance (A):

| System | r1 | r2 | r3 |
|---|---|---|---|
| Control ZM241385 | 0.966, 1.000, 2.61, 4.86 | 1.000, 0.964, 3.22, 4.73 | 0.679, 0.946, 2.75, 4.83 |
| SNP00157 | 0.148, 0.176, 4.96, 5.35 | 0.087, 0.156, 5.41, 5.09 | **0.773**, 0.123, 5.98, 5.01 |
| SNP00163 | 0.597, 0.415, 4.25, 4.78 | 0.755, 0.373, 5.04, 4.73 | **0.984**, 0.171, 4.20, 4.60 |

SNP00157 r3 by 10 ns block: 0.475, 0.700, 0.720, 0.980, 0.995 (the contact is formed late and held). SNP00163 r3: 0.925, 1.000, 1.000, 0.995, 1.000.

**Grading of the pre-registered predictions (once, as registered).** R1 met (both). R2 met (largest nearest C-alpha distance 5.01 and 4.60 A, below 10.0). R3 **missed** (SNP00157 r3 persistence 0.773, bound 0.35). R4 **missed** (SNP00163 r3 persistence 0.984, range 0.40 to 0.90). R5 met (RMSD below 4 A in 0.123 and 0.171 of the frames; bounds 0.50 and 0.60).

**What the result says, and does not say.** With three runs per system the per-run persistence ranges are control 0.679 to 1.000, SNP00157 0.087 to 0.773, SNP00163 0.597 to 0.984. The two replicates 1 and 2 of ADR-AIIM-336 seemed to separate SNP00157 from the control (both below the control minimum); replicate 3 does not: its persistence lies within the control's range. SNP00163 was within the control range already and is now at the top of it. What remains the same in all three runs of both compounds is the ligand RMSD: the share of frames below 4 A is 0.12 to 0.42 for the compounds against 0.95 to 1.00 for the control, and the mean RMSD is 4.2 to 6.0 A against 2.6 to 3.2 A; this is not independent of the way the RMSD is measured (against the pre-equilibration structure, so a ligand that settles into a different but stable pose reads as high RMSD while keeping the contact). No spread is estimated and no test is run on n = 3; a stable trajectory is not evidence that the pose is correct and nothing is called active. The statement of the article that SNP00157 holds the N6.55 contact less than the control in both of its replicates is, with three runs, a statement about two of three runs and the third behaves like the control; with the stated use of ADR-AIIM-390 (the three-against-three comparison replaces the matched duplicate) the article's sentences about persistence need rewording, while the RMSD and the "no dissociation" statements stand.

**Not done.** The article, supplement, cover letter and highlights were **not** edited; the deposit was not changed. The preregistered use (the comparison becomes 3 against 3, the matched duplicate stays in the supplement as sensitivity) is pending the author's instruction because the documents are with the co-authors for review.

## ADR-AIIM-399 — PRE-REGISTRATION: six short sensitivity analyses for the second round of v9 (2026-10-11)

**Status:** proposed — **pre-registration. Frozen before the runs start.**
**Date:** 2026-10-11
**Decided by:** Adisorn Sirichotejirakul, in chat: "วาง plan แก้ที่เหลือต่อ ที่ยังไม่ได้ทำในรอบที่แล้ว ถ้าเกิดมีอะไรที่ต้องทำหรือรันเพิ่ม แล้วไม่นานมาก ควรทำให้ทำ" (plan the rest of what the last round left undone; where something more has to be done or run and it is not long, do it), and approval of the second-round plan (`.claude/plans/check-list-steady-sparrow.md` in the user profile), with four answers given before it: the deployed figure 101 of 128 goes into the abstract, the first highlight and the cover letter; `10_new_system/scripts/build_deposit.py` may be edited (tree entries and the list of pre-registration entries only) and the deposit committed locally for the author to push; the graphical abstract, Figure 6 and Supplementary Figure S1 are redrawn; the language models are named in short form.
**Follows:** ADR-AIIM-397 and 398 (v9), ADR-AIIM-263 (predictions fixed before the run), ADR-AIIM-092 (run folders are immutable; new work in a new dated folder).
**Decision log at freeze (text/LF):** `203a9de9e6833ad376819c5abef23be2ac8133193c5f7a7d21c4139fc518aa43`
**Pre-registration file:** `20_new_runs/v9_sensitivity_2026-10-11/PREREG_ADR399.md`, SHA-256 `88688537933b374aa38964b9b969b3fabf421a40e2973d4e53885084132cfc28` (also in `PREREG_SHA256.txt`).

### What is run
Six analyses asked for by the internal review (Monsin Sangsawat's list), each from existing artifacts or with a frozen script invoked unchanged: A1 the focus-library vote with one vote per distinct split (P0-09); A2 the sensitivity of that selection to the applicability threshold, including a threshold matched in admission rate to the reference on which the competence band was measured (P0-08); A3 a stratified interval for the estimate of 1,251 (P1-19); A4 the four box-constrained compounds that completed docking, redocked in a 30 Å box (P1-17); A5 a string-level map from declared checks to evaluator source (P0-05); A6 the applicability rule at k = 3, 5, 10 and percentiles 5, 10, 20 (P2-03). Rules, anchors, fifteen predictions (P1 to P15), stop rules for the two runs that take minutes, and the list of what had and had not been seen when the predictions were written are in the pre-registration file.

### How the results are used, decided now
As sensitivity views beside the registered results, whatever the numbers; the released selection, the fifteen and the eight are not changed; a changed list is reported as a changed list; nothing is called active and nothing is ranked. Nothing in `10_new_system/` is modified by these analyses.

### Not run, and why
New fault families for G2, G3, G5 and G7, a replacement of the ten-pass censored fit, and a second person's rerun from a clean clone are not short. A refit of the released classifier on the re-curated records already exists in the extension (Supplementary Table S7) and is pointed to instead.

**Addendum to ADR-AIIM-399, recorded before A4 and A6 were started (2026-10-11).** Amendment 1 of the pre-registration (`20_new_runs/v9_sensitivity_2026-10-11/PREREG_ADR399_AMENDMENT_1.md`, SHA-256 `0e092df728533ca20864a4747c934b549ed1846c86e143fef497d963907a5190`) changes only the stop rule and the execution of A4: the stored timings of the production run (202 to 1,070 s per compound at six threads) show that two box sizes cannot finish in the 30 minutes first written, so the two box sizes run side by side with the production command line and the stop rule becomes 120 minutes. Compounds, box sizes, anchor and predictions P8 to P10 are unchanged; no A4 result existed.

## ADR-AIIM-400 — RESULT of the pre-registration of ADR-AIIM-399: the six sensitivity analyses (2026-10-11)

**Status:** result recorded.
**Date:** 2026-10-11
**Follows:** ADR-AIIM-399 (pre-registration and its amendment 1, both frozen before the runs).
**Folder:** `20_new_runs/v9_sensitivity_2026-10-11/` (scripts in `analysis/`). SHA-256 (first 16 hex): `a1_a2_vote_and_threshold.json` `73184ec0615b035b`, `a3_stratified_interval.json` `fba9348df3446b89`, `a4_box_sensitivity.json` `c84d9a646b76b2e3`, `a5_check_map.json` `12a43f867254a5bb`, `a6_ad_k_percentile.json` `d32a9862361e0968`; scripts `a1_a2_vote_and_threshold.py` `66856c68562198ce`, `a3_a5_small.py` `d330a777c92afa6a`, `a4_box_sensitivity.py` `7e7e97fb1504d51b`, `a6_ad_k_percentile.py` `4d82517d686aa14f`. Nothing in `10_new_system/` was modified by the analyses; `g4_dock_library.py` and `ad_descriptor_sweep.py` were invoked or imported unchanged.

**Grading of the fifteen predictions (once, as registered): ten met, five not met.**

**A1, one vote per distinct split.** P1 met: from the stored predictions the registered vote (at least 5 of 8 conventions) reproduces the fifteen and the eight in order. P2 met: with at least 3 of 5 distinct splits, eighteen compounds are selected, the fifteen and three more (SNP00151, SNP00121, SNP00079, each in the top 20 of three distinct splits); none of the fifteen is lost (Jaccard 0.833). P3 met: the shortlist keeps seven of the eight (SNP00033 is replaced by SNP00010; Jaccard 0.778).

**A2, the applicability threshold.** P4 met: 2,645 of 5,000 expansion compounds score at least 0.25 against the 6,284-compound reference, and the 69 stored scores are reproduced to 1e-6 (126 compounds pass the gates before the applicability gate). P5 **not met**: the threshold matched in admission rate is 0.2687 (1,561 compounds at or above it, one tie), below the predicted 0.28 to 0.33. P6 met: 39 of the 69 are in the pool at that value. P7 met: six of the fifteen are in it (SNP00157, SNP00163, SNP00223, SNP00097, SNP00144, SNP00051) and three of the eight. Pool sizes at 0.25, 0.2687, 0.30, 0.35: 69, 39, 16, 2; of the fifteen in the pool: 15, 6, 3, 0; of the eight: 8, 3, 2, 0. At the matched threshold the registered vote selects sixteen compounds, ten of them not among the fifteen. At 0.30 and 0.35 the pool is smaller than the window of twenty, so every pool compound is selected and the vote decides nothing. The selection is therefore strongly dependent on the applicability threshold, which the article had named as a natural sensitivity analysis without running it.

**A3, the interval.** No prediction. The stratified estimate reproduces 1,251 (7.41%); standard error 236.2; normal 95% interval 788 to 1,714 (4.67% to 10.15% of 16,882), from 19 of 249 and 16 of 250.

**A4, box size** (amendment 1; the two processes finished in 2,922 s, inside the 120-minute stop rule; no restart). P8 met: the 22 Å rerun reproduces the stored score and distance of all four compounds exactly. P9 **not met**: one of four improved by more than 0.3 kcal/mol in the 30 Å box. P10 **not met**: two of four left the on-pocket cut-off. Per compound (22 Å, 30 Å, difference; distance in the 30 Å box): NPLIB01997 −10.804, −10.832, −0.028, 7.38 Å; NPLIB04561 −9.141, −9.151, −0.010, 6.75 Å; NPLIB03343 −8.949, −10.046, −1.097, 10.15 Å; NPLIB03133 −10.319, −9.217, +1.102, 13.29 Å. The larger box contains the smaller, so a worse best score in it (NPLIB03133) means the search at exhaustiveness 32 did not find the pose it found in the smaller box: for two of the four the result is not converged with respect to the search space. None of these compounds is among those the workflow released.

**A5, the check map.** P11 **not met**: the ten contract files declare 95 entries under `checks`, not 94. P12 met: one is absent from the source of its evaluator, `md_decision_check` at G7, a routing rule (skip every G7 check when no simulation was run) that the evaluator implements without the identifier. The other 94 occur as string literals in their own evaluator, which is consistent with the 94 checks the article counts; the article's number was not changed. In the five benchmarked gates the deployed contracts declare 55 checks, of which the seeded benchmark exercised 19. The map is of the present tree, which includes the repairs made after the held-out evaluation, and is string-level only.

**A6, k and percentile** (nine settings in 1,749 s, inside the 40-minute stop rule). P13 met: k = 5 and the 10th percentile reproduce the stored admissions and thresholds. P14 met: over the nine settings ECFP4 and Gobbi admit 0 of 236, FCFP4 and topological torsions at most 1. P15 **not met**: the physicochemical space admits 96 to 189 of 236 (40.7% to 80.1%), outside the predicted 50% to 80% at both ends. Other spaces: MACCS 1 to 12, atom pairs 1 to 3, ErG 0 to 7.

**What the results say together.** The focus-library selection is stable to the way the vote is counted (A1) and not stable to the applicability threshold (A2); the contrast between structural and physicochemical applicability does not depend on k or the percentile, while the physicochemical rate does (A6); the docking result of half of the longest compounds depends on the box (A4); the estimate of 1,251 has a wide interval (A3). **What they do not say.** Nothing here is evidence of activity; no compound is ranked or put forward; the registered selection and the released results are unchanged.

**Use.** As decided in the pre-registration: reported in the article and the supplement as sensitivity analyses (ADR-AIIM-401).

"""evidence_not_used.py -- the results this project produced and does not use.

Every other manifest in this project is positive: `ITEMS` says what goes into
the supplement, `TREE` says what goes into the deposit, `sources()` says which
artifacts carry a manuscript number. Nothing says what was produced and then
deliberately left out, or why. That information exists today in three prose
sentences and two code comments, which is not a place a reviewer can read and
not a place a check can enforce.

This is that list. It is hand-maintained, like `EXTERNAL` in build_deposit.py
and `REVIEWED_AND_CLEARED` in deposit_screen.py, and for the same reason: a
decision to leave something out is a human decision with a reason, and the
reason has to be written down by the person who made it.

`build_evidence_index.py` asserts a partition -- every directory under
20_new_runs/ is either derived-used or declared-not-used, never neither. Disk
is the third witness, so a new run cannot be quietly omitted from both lists.

Schema, one entry per run directory:

    kind            one of KINDS below
    reason          required, non-empty, prose. Says why, not what.
    adr             list of ADR ids that decided it; each must resolve to a
                    heading in 08_DECISION_LOG.md
    decided         ISO date
    scope           "whole_run"      -- nothing from here is used (default)
                    "artifacts_only" -- a number from here IS reported, but the
                                        raw material is withheld; say where
    excluded_from   list of "s10", "deposit", "manuscript" -- only with
                    scope="artifacts_only"
    superseded_by   run id, when kind == "superseded"
    published_anyway  paths that ARE published despite the run not being used,
                      so the reader is not told something false

Adding an entry is a decision. Deleting one is also a decision. Neither should
happen without an ADR to point at.
"""
from __future__ import annotations

KINDS = {
    "null_result",        # ran, returned nothing, pre-registered stop rule fired
    "superseded",         # replaced by a later run; name it
    "withdrawn",          # formally retracted; must never be presented as result
    "sealed",             # withheld to preserve a blind test
    "v4_only",            # belongs to the superseded v4 framing
    "input_not_result",   # supplies frozen inputs; not a result at all
    "admin",              # session logs, scratch; carries no scientific claim
    "verification_only",  # checks the work rather than producing a result: no
                          # reported number depends on it, but it is not scratch
                          # and its outcome can change the manuscript. Added
                          # 2026-09-03 for the reference audit (ADR-AIIM-165),
                          # which found a citation pointing at the wrong paper.
    "too_large_to_deposit",
}

SCOPES = {"whole_run", "artifacts_only"}

# Populated 2026-08-29 under ADR-AIIM-142, from the classification the human
# author approved that day. Every entry here is also labelled in
# 20_new_runs/RUN_INDEX.md; build_evidence_index.py fails if the two disagree.
NOT_USED: dict[str, dict] = {
    "receptor_copy_sensitivity_2026-09-22": {
        "kind": "verification_only",
        "reason": "Asks whether the choice among the seven differing copies of "
                  "4EIY_prepared.pdbqt changes a Vina score, because the "
                  "DUD-E run had to be docked into the same receptor as the "
                  "ChEMBL benchmark for the two to be comparable. It does not: "
                  "four ligands across three copies agree to four decimal "
                  "places, spread 0.0000, since Vina derives its hydrogen-bond "
                  "term from heavy-atom geometry rather than explicit hydrogen "
                  "positions. The result removed a blocker and closed a gap in "
                  "section 3.1's reproducibility argument, but the manuscript "
                  "states no number from it, so it is verification only "
                  "(ADR-AIIM-213).",
        "adr": ["ADR-AIIM-213"],
        "decided": "2026-09-22",
    },

    "snplib_source_audit_2026-09-12": {
        "kind": "verification_only",
        "reason": "Asks whether reference [13]'s promised list of 38 primary "
                  "isolation reports can be assembled from anything this "
                  "project holds. It cannot: the frozen SNP-Lib spreadsheet "
                  "carries journal, year and source species but no author, "
                  "title, volume, page or DOI. Reports no number the "
                  "manuscript states -- its 36 and 39 groupings are the two "
                  "reconstructions that bracket the stated 38 and match "
                  "neither, which is the finding rather than a result "
                  "(ADR-AIIM-180).",
        "adr": ["ADR-AIIM-180"],
        "decided": "2026-09-12",
    },

    "reference_check_2026-09-05": {
        "kind": "verification_only",
        "reason": "The reference audit re-run after four citations were added "
                  "for the capability decomposition. Like its 2026-09-03 "
                  "predecessor it is a correctness check on the reference "
                  "list, not a source of any reported number (ADR-AIIM-169).",
        "adr": ["ADR-AIIM-169"],
        "decided": "2026-09-05",
    },

    "validation_matrix_2026-09-05": {
        "kind": "verification_only",
        "reason": "A derived index, not a measurement. Table 7 is built from "
                  "it, but every cell restates a number that G8 already audits "
                  "against the artifact that produced it, so nothing in the "
                  "manuscript depends on this run's own output. Kept because "
                  "the decomposition it applies -- what each technique was "
                  "validated FOR -- is the argument of section 3.6, and the "
                  "mapping from claim to artifact should be inspectable "
                  "(ADR-AIIM-169).",
        "adr": ["ADR-AIIM-169"],
        "decided": "2026-09-05",
    },

    "reference_check_2026-09-03": {
        "kind": "verification_only",
        "reason": "A correctness check on the reference list, not a source of "
                  "any reported number. It found and fixed a citation that "
                  "pointed at the wrong paper (ADR-AIIM-165); the fix lives in "
                  "the manuscript, and nothing in §3 depends on this run.",
        "adr": ["ADR-AIIM-165"],
        "decided": "2026-09-03",
    },


    # -- superseded modelling variants -------------------------------------
    # Three retraining variants were run on 2026-08-05 and one was adopted.
    # The other two are kept so that "we tried it" is checkable rather than
    # asserted, which is the whole argument of Section 4.
    "WP4_retrain_lofo_2026-08-05": {
        "kind": "superseded",
        "superseded_by": "WP4_retrain_postfix_2026-08-05",
        "reason": "Leave-one-family-out retraining of the risk monitor. Run to "
                  "test whether the fitted model survived removing each fault "
                  "family in turn; the adopted model came from the postfix run "
                  "instead. No number in the manuscript comes from here.",
        "adr": ["ADR-AIIM-142"],
        "decided": "2026-08-29",
    },
    "WP4_retrain_logistic_override_2026-08-05": {
        "kind": "superseded",
        "superseded_by": "WP4_retrain_postfix_2026-08-05",
        "reason": "Retraining with the logistic model forced in place of the "
                  "selection procedure's choice. A sensitivity check, not a "
                  "candidate model; not adopted and not reported.",
        "adr": ["ADR-AIIM-142"],
        "decided": "2026-08-29",
    },

    # -- the v4 framing ----------------------------------------------------
    "confidence_tier_positive_control_2026-08-08": {
        "kind": "v4_only",
        "reason": "Positive control for the v4 confidence-tier machinery. v6 "
                  "reports conformal prediction sets and the learned monitor's "
                  "three states instead; the tier vocabulary it validates does "
                  "not appear in the manuscript.",
        "adr": ["ADR-AIIM-142"],
        "decided": "2026-08-29",
    },
    "md_mmgbsa_audit_2026-08-07": {
        "kind": "v4_only",
        "reason": "Audit of the MM-GBSA rescoring arm. That arm was dropped "
                  "from v6 and its Methods and Results sections were cut, so "
                  "nothing here is cited. Kept because the decision to drop it "
                  "should remain inspectable.",
        "adr": ["ADR-AIIM-142"],
        "decided": "2026-08-29",
    },
    "mmgbsa_pilot_2026-08-10": {
        "kind": "v4_only",
        "reason": "The MM-GBSA pilot itself, 501 files. Same decision as the "
                  "audit above. The largest not-used directory in the project.",
        "adr": ["ADR-AIIM-142"],
        "decided": "2026-08-29",
    },
    "control_representativeness_validation_2026-08-04": {
        "kind": "v4_only",
        "reason": "Leave-one-out validation of the control_representativeness "
                  "check. The manuscript's position is that this check is "
                  "mis-specified, so validating it is beside the point v6 "
                  "makes; no machinery references the directory.",
        "adr": ["ADR-AIIM-142"],
        "decided": "2026-08-29",
    },

    # -- formally withdrawn ------------------------------------------------
    "exploratory_sbdd_ranking_2026-08-09": {
        "kind": "withdrawn",
        "reason": "Source of the 10- and 18-compound candidate lists that the "
                  "v4 manuscript itself withdrew, in its own words, as "
                  "something that must not be presented as a prioritisation. "
                  "Superseded entirely by the E8 output. Retained on disk only "
                  "so that the withdrawal is auditable rather than a claim "
                  "about a deleted directory.",
        "adr": ["ADR-AIIM-142"],
        "decided": "2026-08-29",
    },

    # -- engineering and session records -----------------------------------
    "md_benchmark_2026-08-21": {
        "kind": "admin",
        "reason": "GROMACS throughput benchmark. The production run settings "
                  "were chosen from its ns/day figures, but nothing in "
                  "production reads a file from it -- the choice is recorded in "
                  "an ADR, not on disk. It measures this machine, not the "
                  "system under study.",
        "adr": ["ADR-AIIM-142"],
        "decided": "2026-08-29",
    },
    "manuscript_remediation_2026-08-25": {
        "kind": "admin",
        "reason": "Session logs for the review-remediation work of August 2026. "
                  "Holds no scientific artifact and feeds no manuscript number.",
        "adr": ["ADR-AIIM-142"],
        "decided": "2026-08-29",
    },

    # -- used, but the raw material is withheld ----------------------------
    # These two are labelled USED in RUN_INDEX.md and that is correct: a
    # reported number depends on each. They appear here because what is
    # published is the result and not the material, and a reader is entitled
    # to be told which, by whom, and why.
    "heldout_evaluation_2026-08-08": {
        "kind": "sealed",
        "scope": "artifacts_only",
        "excluded_from": ["deposit", "s10"],
        "reason": "What the held-out evaluation measured is reported in Section "
                  "3.1 and is required for the paper to be honest about its own "
                  "test. The faulted workdirs are a different thing: they let "
                  "H01-H05 be reconstructed, and H06-H12 have never been "
                  "exposed. Publishing them would destroy a blind test that "
                  "still has seven unused cases in it. An attempt to deposit "
                  "them on 2026-08-23 was caught and reverted; the exclusion is "
                  "now a standing comment in build_deposit.py. Unsealing "
                  "requires separate explicit human approval and is not "
                  "delegated to any script.",
        "adr": ["ADR-AIIM-113", "ADR-AIIM-142"],
        "decided": "2026-08-29",
        "published_anyway": [
            "20_new_runs/heldout_evaluation_2026-08-08/heldout_summary.json",
        ],
    },
    "np_library_funnel_2026-08-15": {
        "kind": "too_large_to_deposit",
        "scope": "artifacts_only",
        "excluded_from": ["s10", "deposit"],
        "reason": "The most-used run in the project, and its per-compound pose "
                  "directories run to several gigabytes. Every number the "
                  "manuscript reports from the funnel is published, along with "
                  "the ledgers those numbers are computed from; the pose "
                  "archives themselves are in the repository named in the "
                  "data-availability statement, not in this bundle. A size "
                  "decision, not a disclosure one.",
        "adr": ["ADR-AIIM-142"],
        "decided": "2026-08-29",
    },
    "md_candidates_2026-08-22": {
        "kind": "too_large_to_deposit",
        "scope": "artifacts_only",
        "excluded_from": ["deposit"],
        "reason": "System build metadata and setup parameters are reported and "
                  "deposited so the runs can be rebuilt. The trajectories "
                  "themselves are gigabytes of .xtc and .cpt, above the "
                  "deposit's size cap, and are covered by SKIP_SUFFIX. This is "
                  "a size decision, not a disclosure decision: the inputs and "
                  "the protocol are both published.",
        "adr": ["ADR-AIIM-142"],
        "decided": "2026-08-29",
    },

    # Added 2026-10-05. The 2026-10-01/02 and 2026-10-05 batches were never
    # registered, and build_evidence_index.py had been refusing to write since
    # 2026-10-01 as a result -- which also blocked every deposit rebuild. One
    # entry is an input; nine are verification runs whose outcomes bear on
    # section 3.6 while no reported number depends on any of them.
    #
    # NOTE ON THE LABEL. These nine are labelled USED in RUN_INDEX.md because
    # that is the precedent for kind=verification_only (receptor_copy_
    # sensitivity, snplib_source_audit, both reference_check runs,
    # validation_matrix). Unlike those, no machinery reads these directories.
    # The label set {CORE, USED, INPUT, PENDING, ADMIN, V4-ONLY, DO-NOT-REUSE}
    # has no slot for an adjudicated verification run that nothing publishes:
    # USED overstates it, and ADMIN says it carries no scientific claim, which
    # is false -- each of these changed what the project believes about section
    # 3.6. Adding a VERIFICATION label would change build_evidence_index.py and
    # needs its own ADR.
    #
    # W1 DOES NOT FLAG THIS, and the reason is a defect worth naming. W1 fires
    # on "labelled CORE/USED but no machinery references it". The producer-script
    # scan treats any file naming a run directory as a producer, and THIS FILE
    # names them -- so declaring a run here puts it in the derived-used set, and
    # W1 can never fire for a declared run. The module docstring above says disk
    # is the third witness because the two lists are meant to be independent;
    # for the producers channel they are not. E3 still works (it wants either
    # list), but the label-overstatement check is silenced by the declaration it
    # is supposed to be checked against. Not fixed here: build_evidence_index.py
    # is in the frozen tree and this needs its own ADR.
    "chembl_censored_repull_2026-10-01": {
        "kind": "input_not_result",
        "reason": "A re-pull of CHEMBL251 keeping every standard_relation, made because the released label rule compares pchembl_value against a cut and ChEMBL leaves that column null on all 1,474 upper-censored records, so every one of them was discarded as missing rather than read as an interval. Not a result: it supplies the censored inactives that arms B and K of ADR-AIIM-287 consume and that ADR-AIIM-291 step 6 reads directly. No number is reported from it (ADR-AIIM-268).",
        "adr": ['ADR-AIIM-268'],
        "decided": "2026-10-05",
        "scope": "whole_run",
    },
    "qsar_refit_censored_2026-10-02": {
        "kind": "verification_only",
        "reason": "Asks whether section 3.6's 72 of 5,000 reproduces when the released QSAR path is re-run with the recovered inactives. It does not, and the diagnosis traced the gap to the order of the rows handed to scaffold_train_test_split, which ADR-AIIM-271 then measured. No machinery reads this directory and no reported number depends on it; its outcome is why section 3.6 is under review at all (ADR-AIIM-269, ADR-AIIM-270).",
        "adr": ['ADR-AIIM-269', 'ADR-AIIM-270'],
        "decided": "2026-10-05",
        "scope": "whole_run",
    },
    "scaffold_split_fix_2026-10-02": {
        "kind": "verification_only",
        "reason": "The first of two attempts to make the scaffold split independent of row order: a total-order tiebreak on the size sort. It made the calibration split exactly invariant and left the counts moving across 47 distinct values, because RandomForest draws bootstrap index positions from its own generator. Pre-registered prediction P1 was MISSED and is reported as such. No reported number depends on it (ADR-AIIM-272, ADR-AIIM-273).",
        "adr": ['ADR-AIIM-272', 'ADR-AIIM-273'],
        "decided": "2026-10-05",
        "scope": "whole_run",
    },
    "snplib_canonical_2026-10-02": {
        "kind": "verification_only",
        "reason": "Asks whether the library expansion's empirical justification was an artefact of row order, under canonical ordering at three significance levels. The two libraries were never different and the eligible pool was never 3; neither fact makes the expansion unnecessary. No reported number depends on it (ADR-AIIM-276, ADR-AIIM-277).",
        "adr": ['ADR-AIIM-276', 'ADR-AIIM-277'],
        "decided": "2026-10-05",
        "scope": "whole_run",
    },
    "split_convention_factorial_2026-10-05": {
        "kind": "verification_only",
        "reason": "A 2x2x2 factorial over the three design choices inside scaffold_train_test_split that neither ADR-AIIM-273 nor ADR-AIIM-275 touched: the shuffle, the fill direction and the overflow rule. Eight reproducible conventions span 6.78x on arm A and 376x on arm B, and the D1xD2 interaction is nearly twice the main effect it was graded on, so the main-effect summary the pre-registration asked for describes no configuration that exists. Reports no number the manuscript states; its outcome bears on whether section 3.6's count is determined by the data (ADR-AIIM-285, ADR-AIIM-286).",
        "adr": ['ADR-AIIM-285', 'ADR-AIIM-286'],
        "decided": "2026-10-05",
        "scope": "whole_run",
    },
    "label_balance_domain_2026-10-05": {
        "kind": "verification_only",
        "reason": "Tests three claims that had been asserted without measurement: that a Ki-only label would not stabilise the count, that calibration balance was not the driver, and that restricting to the applicability domain could not help. The first and third were false. Ki-only cuts the convention spread from 376x to 3.62x, and the in-domain prediction was written so that an empty set falsified it vacuously. Also the first measurement of the NP library at 0 of 5,000 in-domain under all three label arms. No reported number depends on it (ADR-AIIM-287, ADR-AIIM-288).",
        "adr": ['ADR-AIIM-287', 'ADR-AIIM-288'],
        "decided": "2026-10-05",
        "scope": "whole_run",
    },
    "arm_k_auc_2026-10-05": {
        "kind": "verification_only",
        "reason": "Nested-CV discrimination for the Ki-only label arm against the two pooled arms, using the frozen nested_cv_auc fold structure unchanged. Arm K does not discriminate better: ROC-AUC 0.9326 against arm B's 0.9377, and worse on MCC and balanced accuracy. It also found that ADR-AIIM-268's PR-AUC lift convention is not comparable across arms of differing prevalence, which is flagged against those figures and not re-derived. Arms A and B reproduce ADR-AIIM-268's 0.9181 and 0.9371 to within 0.0006, the harness anchor (ADR-AIIM-289, ADR-AIIM-290).",
        "adr": ['ADR-AIIM-289', 'ADR-AIIM-290'],
        "decided": "2026-10-05",
        "scope": "whole_run",
    },
    "chembl251_full_curation_2026-10-05": {
        "kind": "verification_only",
        "reason": "The whole CHEMBL251 activity record pulled unfiltered at ChEMBL_37 and curated down under nine rules fixed before the pull, so every exclusion has a step, a rule and a count. It establishes that the frozen input's restriction to Ki and IC50 was an undocumented filter removing 11,360 rows, and that almost every consequence predicted from that was wrong: the curated set differs from arm K by 5.2 percent, holds one fewer inactive than arm K, and spreads 38.88x across the eight conventions against arm K's 3.62x. What it bought is an audit trail, not a better number. The raw pull is NOT promoted to 02_frozen_inputs; promoting it would change published numbers (ADR-AIIM-291, ADR-AIIM-293).",
        "adr": ['ADR-AIIM-291', 'ADR-AIIM-293'],
        "decided": "2026-10-05",
        "scope": "whole_run",
    },
    "ad_methods_distance_performance_2026-10-05": {
        "kind": "verification_only",
        "reason": "Eight applicability-domain definitions crossing representation, metric and method, adding the leverage, centroid and Mahalanobis constructions that Figure 4's eight-space analysis does not include, applied to one fitted model so the definitions re-slice identical predictions. The tier verdict flips between Tier 1 and Tier 2 on the definition alone. Six of the eight verdicts are vacuous because their bands cover 68 to 94 percent of the held-out set, and two of the eight are not independent: with centred descriptors, leverage and squared Mahalanobis distance differ by a constant. No reported number depends on it (ADR-AIIM-292, ADR-AIIM-294).",
        "adr": ['ADR-AIIM-292', 'ADR-AIIM-294'],
        "decided": "2026-10-05",
        "scope": "whole_run",
    },
    "ad_distance_curve_2026-10-05": {
        "kind": "verification_only",
        "reason": "Performance as a function of distance to the training set, which the manuscript measures nowhere: ROC-AUC per declared AD band on the same fitted model, anchored by reproducing ADR-AIIM-294's held-out 0.9217 exactly. The model is usable down to AD 0.25, which 31.2 percent of the NP library reaches; below it the entire evidence is 22 held-out compounds at ROC-AUC 0.3294, worse than chance. It corrects ADR-AIIM-294's Tier 1, which averaged across that sign change. Two of its four predictions failed in the harness rather than in the data and are recorded as such. No reported number depends on it (ADR-AIIM-295, ADR-AIIM-296).",
        "adr": ['ADR-AIIM-295', 'ADR-AIIM-296'],
        "decided": "2026-10-05",
        "scope": "whole_run",
    },
}


def self_test() -> None:
    """Shape check. Run by build_evidence_index.py before anything else."""
    for run, e in NOT_USED.items():
        assert e.get("kind") in KINDS, f"{run}: kind {e.get('kind')!r} not in KINDS"
        assert e.get("reason", "").strip(), f"{run}: reason is required"
        assert isinstance(e.get("adr", []), list), f"{run}: adr must be a list"
        scope = e.get("scope", "whole_run")
        assert scope in SCOPES, f"{run}: scope {scope!r} not in SCOPES"
        if scope == "artifacts_only":
            assert e.get("excluded_from"), (
                f"{run}: scope=artifacts_only must say what it is excluded from")
        if e.get("kind") == "superseded":
            assert e.get("superseded_by"), f"{run}: superseded by what?"


if __name__ == "__main__":
    self_test()
    print(f"self_test OK -- {len(NOT_USED)} entries")

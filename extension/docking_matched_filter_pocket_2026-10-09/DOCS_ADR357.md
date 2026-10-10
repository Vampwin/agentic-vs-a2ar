
## ADR-AIIM-357 — Results of ADR-AIIM-353 to 356 entered in the v8 working article and supplement, with placeholders for what is not yet known (2026-10-09)

**Decided by:** Adisorn Sirichotejirakul — "put in what can be put in, and a [placeholder] where something is still pending". The Figure 6 option the author was offered was not chosen; I took the lightest one (option C of `RESULTS_SUMMARY_DRAFT.md`): **no figure change**, a new table in the supplement, and a sentence in the legend. The deposit builder was not changed (that needs approval).

### What was entered (tracked changes by "Claude"; edits inside my own earlier insertions made directly)

**Article.** The methods paragraph on orthogonal pose checks: the matched-search-space runs named, with a placeholder for the DiffDock-Pocket citation. Section 3.8: the sentence "the on-pocket result … rests on Vina alone" replaced by the matched comparisons (Vina over the whole receptor 0.37 Å on the control and 2 of 15 survivors within 8.0 Å; none of fifteen poses agreeing within 2.0 Å blind against blind, median 15.9 Å; DiffDock-Pocket 0.35 Å on the control, none of fifteen within 2.0 Å of the Vina pose, median 12.8 Å; the two arms that failed their control named and not interpreted), ending "engine-dependent even when both engines search the same space", plus a placeholder for the two variants still running. The molecular-dynamics paragraph: one sentence that Vina over the whole receptor preferred other places for SNP00157 and SNP00163 (11.9 and 12.0 Å), so the simulations tested the in-site poses. The Figure 6 legend: one sentence that panels D and E compare unmatched search spaces, pointing to Table S12.

**Supplement.** The Section 14 sentence "at 20 samples per complex" corrected to 10 (the YAML overrode the command line; ADR-AIIM-355); the closing sentence of the "do not show which engine is right" paragraph now points to the matched comparison; the molecular-dynamics paragraph of Section 17 gets the same sentence as the article, with the DiffDock-Pocket distances (9.9 and 3.4 Å from the Vina poses); a new block after Table S10's discussion: the methods of the four new runs (including the installation facts of ADR-AIIM-355), **Supplementary Table S12** (16 rows, distances and RMSDs from `matched_poses.csv` and `filter_pocket_poses.csv`), and the results paragraph. Two yellow placeholders: the pending DiffDock-Pocket variants and the citation; and that the new runs are not yet in the deposit.

### Placeholders left (yellow)

Article: the DiffDock-Pocket citation (reference list to be edited by the author), the two descriptive DiffDock-Pocket variants (ADR-AIIM-358), plus the two already open (abstract word limit; the repository still private). Supplement: the same citation and variants, and the deposit of the new runs.

### Checks

Reject-all text equals the reviewed v8 for both documents; tracked ids unique (article 62, supplement 1,206); the article is 49 pages and the supplement 67 in Word, Table S12 rendered and read. G8 on the accepted text: 114 pass, 0 mismatches, 322 literal-not-found (unchanged), `deposit_commit` PASS at `3e282e5` and the deposit working tree clean; the original `g8_number_traceability.json` restored (`661e3bc6…`). The new-number audit, run on a scratch copy with the two new run folders added (the script in the extension folder was not changed), finds every number of the new insertions in an artifact of those runs (0 unmatched in both documents); as before it detects typos, not a wrong quantity. Yellow copies rebuilt. Backups `_superseded_KDSS_v8_before_matchedsearch.docx` and `_superseded_Supplementary_v8_before_matchedsearch.docx`.

### Not done

The two new run folders are not in the deposit; `build_deposit.py` still needs the author's approval to list them, after which the deposit is rebuilt, amended, and the commit hash and counts in both documents updated. Nothing was pushed. The INCIDENTS draft was not extended with this session's new entries (the grading script that read a crystal file without bond orders as NaN; the Vina script's box-centred distance column; the 20-sample statement; the batch size). Article `e8023a99b0df7a3a…`, supplement `35c5b29d43921364…`.

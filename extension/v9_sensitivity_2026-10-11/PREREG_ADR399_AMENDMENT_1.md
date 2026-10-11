# Amendment 1 to the pre-registration of ADR-AIIM-399 (written 2026-10-11, before A4 or A6 was started)

Only the stop rule and the execution details of A4 change. Its compounds, box sizes, anchor and predictions P8 to P10 stand.

Why: the stop rule of 30 minutes was written without reading the stored timings. The production run records, for the four
compounds, 202, 972, 916 and 1,070 s each at `--cpu 6` (`np_library_funnel_2026-08-15/full_scores_chunk*.csv`, column
`elapsed_s`), about 53 minutes for the four in sequence; two box sizes cannot finish in 30 minutes. No docking result of A4
existed when this was written.

Amended A4:
- Two processes run side by side, one per box size (22 Å and 30 Å), each docking the four compounds in sequence with the
  production command line (`--id-col compound_id --smiles-col SMILES --n-conformers 8 --exhaustiveness 32 --n-poses 10 --cpu 6`),
  the production receptor (`snp_screen_2026-08-04/G2_preparation/4EIY_prepared.pdbqt`) and the production box centre; the
  22 Å run uses the production box file itself, the 30 Å run a copy of it in which only `size` is changed to 30, 30, 30.
- The input rows are those of the four compounds in `np_library_funnel_2026-08-15/np_library_dock_ledger.csv`.
- Stop rule: 120 minutes of wall time for the two processes together; a run that has not finished by then is stopped and
  reported as not completed, and a partial result is not interpreted.
- The on-pocket distance of P10 is the one the frozen script reports (centroid distance to the box centre, cut-off 8.0 Å).

A6 is started while A4 runs; its stop rule (40 minutes) is counted from its own start.

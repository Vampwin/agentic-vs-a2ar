# md_snplib_duplicate_2026-10-08 — analysis of SNP00157 and SNP00163 (ADR-AIIM-331)

Scripts only until the replicates finish; no result is recorded here yet.

Run order, once a replicate's `prod_rN.gro` exists (for SNP00157 r1 it already does):

1. `run_candidate_md11.sh <SNP00157|SNP00163> prod_r1|prod_r2`, in WSL (distro A2A-MD). The frozen
   `md_11_analyse_pose_persistence.sh`, unmodified: receptor C-alpha fit, ligand RMSD against `em.tpr`,
   bilayer box. Writes `analysis_<deffnm>/` inside the compound's `04_run`.
2. `python analysis/adr331_analysis.py --anchors` checks the machinery without reading a candidate trajectory.
3. `python analysis/adr331_analysis.py --analyse` builds the compound's own N6.55 index groups, runs the
   mindist, grades S1 to S4 and D3, and writes `adr331_results.json` and `adr331_per_run.csv`.
   It refuses to grade if an anchor fails.

Anchors (all met on 2026-10-08, before any candidate result was read):
A1 the index groups built from the control's structure equal the borrowed control index file;
A2 the readers reproduce the stored control numbers for all three control runs (RMSD mean, fraction
below 4 A, N6.55 persistence, bilayer area change); A3 the Asn is found by sequence position, the TM6
spacing passes and the polar-atom count equals N plus O of the SMILES, for both compounds; A4 the nearest
Anchor-to-Ligand-polar distance in each pre-equilibration structure is within 0.5 A of the Vina-pose value
(3.17 against 2.87 for SNP00157, 2.80 against 2.67 for SNP00163).

Two things to know:
- The series are written every 50 ps for the candidates, against 100 ps for the control, so a candidate
  run has 1,001 frames and a control run 501. Every statistic is a fraction or a mean, so the density
  does not bias it.
- S2 uses "strictly below 4.0 A", the measure behind the control figures quoted in ADR-AIIM-331 (94.6 to
  100 %). That ADR words S2 as "at most 4.0 A"; the two differ only for a frame at exactly 4.000 A.
- While testing the N6.55 plumbing on 2026-10-08 I created `analysis_prod_r1/N6_55.*` for SNP00157 and
  checked only that it had 1,001 frames covering 0 to 50 ns. I did not print its values.

## Text generation for the manuscript and supplement (prepared 2026-10-08, before any result)

- `analysis/make_md_text.py` reads `adr331_results.json` and writes `md_text.json`: the Section 3.8 paragraph, the
  supplement paragraphs, Table S11 rows and caption, and a proposed abstract clause. Every number is read from the
  results; the only words chosen by the script are the ones that depend on which predictions held, in branches for
  all met, none met, mixed, and incomplete runs. `--selftest` prints three synthetic scenarios and reads no real result.
- `analysis/apply_md_text.py --main <unpacked> --supp <unpacked>` writes that text into the unpacked v8 documents as
  tracked changes, replacing the two placeholders about these simulations, and inserts Table S11. The abstract is not edited.
  It was tested on scratch copies with the synthetic text: both files validated, rejecting all changes returned the
  reviewed v8 text exactly, and Word rendered the table. Nothing from that test was installed.
- It refuses to generate if an anchor failed.

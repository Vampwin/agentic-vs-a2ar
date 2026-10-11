# Agentic virtual screening against the adenosine A2A receptor

Code and verification data accompanying the manuscript
*"Auditing LLM-agent-generated computational drug-discovery workflows with executable contracts: fault detection and decision stability in an A2A receptor natural-product virtual screen"*, submitted to *Artificial Intelligence in Medicine*.

**Every candidate compound reported in the manuscript is a computational
hypothesis. No wet-laboratory validation has been performed.**

---

## What this repository is for

The manuscript makes two kinds of claim, and this repository exists so both can
be checked rather than taken on trust:

1. **Standard quality-control gates reject unfamiliar chemistry.** The
   236-compound Thai orchid natural-product library is 0 % in-domain in four
   structural descriptor spaces and 64.4 % in-domain in physicochemical space —
   the same library, the same rule. Under that same rule the 497-compound
   expanded selection reaches 86.3 % physicochemically and stays at 0 %
   structurally. The data behind both is in `verification/`.

2. **The pipeline's own code was verified, and the verification is reported as a
   result.** Contracts, evaluators, the fault-injection suite and the test suite
   are all in `pipeline/`, along with every analysis and audit script the
   manuscript cites.

---

## Layout

| Path | Contents |
|---|---|
| `pipeline/contracts/` | The ten frozen G0–G8 audit contracts (94 checks) |
| `pipeline/evaluators/` | The deterministic evaluator modules the acting agent cannot override |
| `pipeline/fault_injection/` | Eighteen fault-family definitions; eleven implemented, seven declared stubs |
| `pipeline/tests/` | The test suite — 682 tests across 30 modules |
| `pipeline/scripts/` | Every analysis, figure and audit script referenced in the manuscript |
| `pipeline/kdss_audit/` | The append-only artifact registry |
| `verification/` | The evidence behind Sections 2 and 3, copied as the analysis read it |
| `figures/` | Figure sources and their provenance records; `figures/graphical_abstract/` holds the graphical abstract's builders and outputs (added 2026-10-10) |
| `preregistration/` | The decision-log entries the manuscript and its supplement call pre-registered, with the entries that record their results, copied verbatim (added 2026-10-10); the full log is not deposited |
| `extension/` | The post-freeze extension (ADR-AIIM-285 onward, added 2026-10-09): the runs behind Section 3.8 of the manuscript and Section 17 of its supplement, as scripts, series, tables, ligand parameter files and machine-readable results; since 2026-10-10 also the matched search-space program (Vina over the whole receptor, DiffDock-L filtered to the site box, DiffDock-Pocket) and the one style of the figures; since 2026-10-11 also the seeded benchmark under the harness contract and under the contracts as deployed (`extension/deployed_path_benchmark_2026-09-23`, the figures of Section 3.1 and Supplementary Table S6), six pre-registered sensitivity analyses (`extension/v9_sensitivity_2026-10-11`), Figure 6 and Supplementary Figure S1 with the legend wording of the article (`extension/figures_v9_2026-10-11`) and the package lists of the three environments (`extension/environment_freeze_2026-10-11`) |
| `MANIFEST.json` | SHA-256 of every file, plus what was excluded and why. The hash and byte count are those of the content as a clone of this repository holds it: text files with LF line ends (as `.gitattributes` makes git store them), binary files as they are. `sha256sum` on a checkout therefore reproduces them |

### Figure files and the manuscript's figure numbers

The files under `figures/` keep the numbers they were given when they were made,
which are not the manuscript's. Current correspondence:

| Manuscript figure | File |
|---|---|
| 1, architecture | `figures/figure1_architecture_2026-08-21` |
| 2, inside one gate | `figures/figure2_inside_a_gate_2026-08-21` |
| 3, screening performance | `extension/figures_style_harmonised_2026-10-10/out/figure3_screening_performance_harmonised` (script `fig3_screening.py`); the earlier version: `figures/figure4_screening_performance_2026-08-24` |
| 4, applicability domain by representation | `extension/figures_style_harmonised_2026-10-10/out/figure4_ad_multispace_harmonised`; the earlier version: `figures/figure_ad_multispace_2026-08-25` |
| 5, permutation sensitivity | `extension/figures_style_harmonised_2026-10-10/out/figure5_permutation_harmonised`; the earlier version: `figures/figure7_permutation_distribution_2026-10-02` |
| 6, pose comparison with the search space given to both engines (A the control, B and C SNP00157 and SNP00163, PyMOL with one camera; D DiffDock-Pocket against Vina, E Vina in the box against the whole receptor) | `extension/figures_v9_2026-10-11/Figure_6_pocket` (script `analysis/figure6_pocket.py` there, which reads the PyMOL renders made by `extension/figures_docking_md_2026-10-09/analysis/render_pocket_pymol.py`; the version of 2026-10-09, whose legend reads "other survivors": `extension/figures_docking_md_2026-10-09/Figure_6_pocket`; data `extension/docking_matched_filter_pocket_2026-10-09` and `extension/docking_matched_search_space_2026-10-09`) |
| Supplementary Figure S1, DiffDock-L run blind against Vina in the box (the earlier Figure 6) | `extension/figures_v9_2026-10-11/Figure_6_docking_snplib` (the version of 2026-10-09, whose legend reads "other survivors": `extension/figures_docking_md_2026-10-09/Figure_6_docking_snplib`); the earlier two-panel version with NPLIB01643: `figures/figure5_pose_overlay_2026-08-24` |
| Supplementary Figure S2, membrane simulations of all five docked compounds over the first 50 ns, with three runs of SNP00157 and SNP00163 | `extension/md_snplib_replicate3_2026-10-10/figures/Figure_S2_md_all_compounds_3runs` (builder `figure_S2_md_three_runs.py`); the version with two runs of each, drawn by `extension/figures_docking_md_2026-10-09/analysis/figure7_md_snplib.py`: `extension/figures_docking_md_2026-10-09/Figure_S2_md_all_compounds`; the control and SNP00157 and SNP00163 alone, which was Figure 7 until 2026-10-10 and is no longer in the article: `Figure_7_md_snplib` in the same folder; the earlier version: `extension/figure8_window50_2026-10-07` |
| graphical abstract | `figures/graphical_abstract/figure_graphical_abstract_v9_2026-10-11` (builder `build_graphical_abstract_v9.py`; 2.5 to 1, with the seeded-fault result under the harness contract and under the contracts as deployed); the earlier layouts: `figures/graphical_abstract/figure_graphical_abstract_v8_funnel_2026-10-10`, `figures/graphical_abstract/figure_graphical_abstract_v7_2026-10-03` and `figures/figure_graphical_abstract_2026-08-21` |

Two files are **not** in the current manuscript: `figures/figure6_pose_persistence_2026-08-31`
is the earlier 100 ns version of Figure 7, and `figures/figure8_m6_panel_canonical_2026-10-02`
is the source of a figure that was removed (candidate eligibility after canonical ordering).
They are kept because the manuscript's earlier text cited them.

## What can be re-run here, and what cannot

Established by copying this repository to a clean directory and running each
step there, rather than asserted.

**The test suite runs and passes.**

```
cd pipeline
python -m pytest -q          # 682 passed
```

The working directory matters: several tests open contract files by relative
path, so running from the repository root fails at collection.

**The analysis scripts resolve their own layout.** They name their inputs by
the development tree's paths -- `20_new_runs/<run>/<file>` -- while this
repository regroups the same evidence by subject under `verification/`.
`pipeline/scripts/kdss_paths.py` maps between the two using `PATH_MAP.json`,
which `build_deposit.py` writes from the manifest and from the supplement
bundle's own record of where each file came from. Nothing in that chain is
guessed; a path the map does not carry raises, naming the file, rather than
resolving to something that is not it.

Until 2026-09-13 these scripts held two absolute paths on the first author's
machine and exited if neither existed. On that machine a copy of this
repository silently read -- and wrote -- the development tree instead of
itself; anywhere else they stopped. Both are fixed.

**Four things still cannot be re-run here, each for a stated reason.**

| Not runnable | Why |
|---|---|
| `g8_number_traceability.py` end to end | It checks the manuscript's text against the artifacts. The manuscript is not deposited, and one input it reads is a ChEMBL extract that cannot be (below). Its full output is `verification/structure_and_protocol/g8_number_traceability.json` (the same file was also deposited under `verification/summary/` until 2026-10-10; the builder keeps one path for one artifact). |
| Anything reading ChEMBL records | ChEMBL is CC BY-SA 3.0. Depositing records would put a share-alike term on evidence released under CC BY 4.0, so only aggregate statistics are here. |
| Docking and DiffDock | Need Vina, DiffDock, a GPU and the receptor preparations. This repository supplies the resulting poses and scores instead. |
| Molecular dynamics | Trajectories run to gigabytes. The assembled starting structure is deposited, because four numbers in the Methods are counted from it; the trajectories are available from the corresponding author. |

## Environment

Confirmed by interrogating the execution environments directly rather than from
memory. RDKit differs between them, and that is not an error: the conformer
generation step shells out to WSL.

| Component | Version | Where |
|---|---|---|
| Python | 3.12.10 | Windows — evaluators, standardisation, QSAR, tests |
| RDKit | 2026.03.3 | Windows |
| RDKit | 2025.09.5 | WSL `Ubuntu-22.04`, conda env `diffdock` |
| scikit-learn | 1.9.0 | Windows |
| pandas | 3.0.3 | Windows |
| NumPy | 2.4.6 | Windows |
| AutoDock Vina | 1.2.7 | WSL, Forli lab Python bindings |
| DiffDock | commit `85c49b6`, tag `v1.1.3-2-g85c49b6` | WSL |
| PyTorch / CUDA | 2.7.0+cu128 / 12.8 | WSL |

The test suite needs only the Windows-side set. Docking and DiffDock are needed
only to regenerate poses, which this deposit supplies rather than requires.

---

## The post-freeze extension

`extension/<run>/` holds the runs that were done after the verification study was
frozen and that the manuscript reports separately (Section 3.8, Supplementary
Section 17): the censored-fit reproducibility defect, the eight split
conventions, the SNP-Lib funnel, Vina docking of the 69 in-domain compounds,
DiffDock-L on the fifteen survivors and on a DUD-E subsample, the matched search-space
program (Vina over the whole receptor, DiffDock-L filtered to the site box, DiffDock-Pocket),
the three-run molecular-dynamics control, SNP00157 and SNP00163 in duplicate, and their third replicates with the post hoc analyses of them (pose families, receptor flexibility, hydrogen bonds, and the check of the trajectories for periodic-image jumps: `extension/md_snplib_replicate3_2026-10-10`). Each was
pre-registered in the decision log with the SHA-256 of its predictions before it
ran. The full log is not deposited; the entries the manuscript relies on, and the
entries that record their results, are copied verbatim into
`preregistration/PRE_REGISTRATION_ENTRIES.md`, with the SHA-256 of the whole log at
the time of copying.

It is copied from the development tree through an allow-list of suffixes (see
`EXTENSION_SUFFIX` in `pipeline/scripts/build_deposit.py`): scripts, JSON and CSV
results, force-field parameter files, `.xvg` series, `.mdp` and `.ndx` files,
figures. Not deposited: trajectories, run inputs, checkpoints, structures and
topologies, logs, the 3,685 per-complex pose files of the DUD-E DiffDock subsample
(the scored table the analysis read is deposited), and, for the two matched search-space runs,
all poses below rank 1 (about 3,500 files; the rank-1 pose of each complex is what the analysis and the figures read), and the raw ChEMBL pull, which
the screen refuses and which is a public release.

**These scripts do not resolve their own layout.** They were written against the
development tree and name their inputs by absolute paths on the first author's
machine (`D:\Agentic VS Workflow\...`), not through `kdss_paths`. They are
deposited as they were run, so that what produced each number can be read. Running
one elsewhere means changing its root, and most also need inputs that are not here
(the receptor preparations, the trajectories). Their results can be checked
against each other; they cannot be regenerated from this repository alone.

## What is deliberately not here

**Inputs that are already archived elsewhere.** Re-depositing them would add a
second copy that can diverge from the original, and add 107 MB. They are named
with their identifiers in `MANIFEST.json` under
`external_inputs_not_deposited`, so a reader can obtain exactly the same bytes:

| Input | Identifier |
|---|---|
| LOTUS v11 frozen snapshot | doi:10.5281/zenodo.19360665 |
| ChEMBL A2AR bioactivity records | doi:10.1093/nar/gky1075 |
| DUD-E A2AR decoy set | doi:10.1021/jm300687e |
| A2AR crystal structures | PDB 4EIY, 5OLH |

**Held-out fault families.** Sealed. The manuscript reports the *result* of
evaluating H01-H05 once (one of five detected, with the root cause of each miss),
and `verification/heldout_outcome/` holds that summary. The faulted workdirs are
not deposited, because they would let the cases be reconstructed, and H06-H12 have
never been exposed or evaluated. Exposing any of that requires a separate approval
that has not been given.

**Molecular-dynamics trajectories.** The manuscript reports exploratory
membrane simulations of three natural-product candidates and, in the extension, of
two SNP-Lib compounds, each against the ZM241385 control. The trajectories run to
gigabytes and are available from the corresponding author on request; the series
derived from them (ligand RMSD, the N6.55 distance, bilayer box) and the ligand
parameter files are in `extension/`.

**Pose archives.** Several GB of per-compound docking poses from the verification
study, available on request. The extension deposits the Vina poses of the 69
in-domain SNP-Lib compounds and the DiffDock poses of the fifteen survivors.

---

## How the deposit was screened

Nothing here was hand-picked. `pipeline/scripts/build_deposit.py` assembles the
tree, and `pipeline/scripts/deposit_screen.py` reads every file first, looking
for assay-derived values, credentials, and absolute paths from the author's
machine.

The screen's own history is worth knowing, because it is the kind of failure
this project reports rather than hides. Its first version was a blunt
"assay keyword near a number" pattern, and it refused seven files — every one a
false positive, including one of the ten contracts the paper's argument rests
on. What it had matched was a threshold specification, a synthetic test fixture,
and the letters "camp" inside a plant genus name. A screen that silently removes
a contract is worse than no screen, because the deposit still looks complete.

The pattern now describes what an assay value actually looks like, carries a
self-test that the build refuses to proceed without, and anything it still flags
must be read and cleared with a published reason. Exactly one file is currently
cleared: the screen itself, which trips on its own self-test fixtures.

---

## Licence

The pairing is **MIT for the code** under `pipeline/` and **CC BY 4.0 for the data** under `verification/` and `figures/`; `LICENSE` carries both and `CITATION_AND_LICENCE.md` the audit they follow from. Nothing copyleft is linked in: GROMACS, AmberTools, Open Babel and pdb2pqr are run as separate executables, never imported, and every in-process import is BSD, MIT or Apache-2.0. **LOTUS must be attributed** wherever the files under `verification/` are reused.

## Citation

**[AUTHOR ACTION]** Add the manuscript citation once accepted.

If you cite this repository directly, cite a commit rather than a branch — a
branch moves, a commit hash does not:

```
<author list>. Agentic virtual screening against the adenosine A2A receptor: verification of the agent-generated code as a reportable result: code and verification data.
GitHub repository, commit <sha>. <year>.
```

Archiving a release through Zenodo would additionally mint a DOI for that exact
snapshot, which some journal data-availability policies require. That is a
repository setting plus a release; it needs no changes to anything here.

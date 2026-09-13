# Agentic virtual screening against the adenosine A2A receptor

Code and verification data accompanying the manuscript
*"Agentic virtual screening against the adenosine A2A receptor: verification of the agent-generated code as a reportable result"*, submitted to *Artificial Intelligence in Medicine*.

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
| `figures/` | Figure sources and their provenance records |
| `MANIFEST.json` | SHA-256 of every file, plus what was excluded and why |

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
| `g8_number_traceability.py` end to end | It checks the manuscript's text against the artifacts. The manuscript is not deposited, and one input it reads is a ChEMBL extract that cannot be (below). Its full output is `verification/summary/g8_number_traceability.json`. |
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

**Held-out fault families.** Sealed. They are referenced nowhere in this version
of the manuscript, which reports development families only, and exposing them
requires a separate approval that has not been given.

**Molecular-dynamics trajectories.** The manuscript reports the membrane
protocol and its positive control, not a candidate result; trajectory files run
to gigabytes and are available from the corresponding author on request.

**Pose archives.** Several GB of per-compound docking poses. Available on
request.

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

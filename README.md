# Agentic virtual screening against the adenosine A2A receptor

Code and verification data accompanying the manuscript
*"Agentic virtual screening against the adenosine A2A receptor: verification of the agent-generated code as a reportable result"*, submitted to *Artificial Intelligence in Medicine*.

**Every candidate compound reported in the manuscript is a computational
hypothesis. No wet-laboratory validation has been performed.**

---

## What this repository is for

The manuscript makes two kinds of claim, and this repository exists so both can
be checked rather than taken on trust:

1. **Standard quality-control gates reject unfamiliar chemistry.** A 236-compound
   Thai orchid natural-product library is 0 % in-domain in four structural
   descriptor spaces and 86.3 % in-domain in physicochemical space — the same
   library, the same rule. The data behind that is in `verification/`.

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
| `pipeline/tests/` | The test suite — 675 tests across 30 modules |
| `pipeline/scripts/` | Every analysis, figure and audit script referenced in the manuscript |
| `pipeline/kdss_audit/` | The append-only artifact registry |
| `verification/` | The evidence behind Sections 2 and 3, copied as the analysis read it |
| `figures/` | Figure sources and their provenance records |
| `MANIFEST.json` | SHA-256 of every file, plus what was excluded and why |

## Reproducing the audit

`pipeline/scripts/g8_number_traceability.py` is the one to start with. It
recomputes every load-bearing number in the manuscript from the artifact it came
from and checks it against the literal text of the draft — 57 checks, and it
fails rather than warns when a number cannot be re-derived.

```
python pipeline/scripts/g8_number_traceability.py --manuscript <path to the draft>
```

Others worth running directly: `analyse_metamorphic.py` (the four invariants,
including the one that fails), `figure_ad_multispace.py` (Figure 3 from its
source data), `build_references_v6.py --check` (the reference list).

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

**[AUTHOR ACTION — not yet set.]** The intended pairing is MIT for the code
under `pipeline/` and CC BY 4.0 for the data under `verification/` and
`figures/`, which is the usual arrangement for this kind of deposit and is
compatible with the licences of the external inputs listed above. Confirm and
add `LICENSE` before making the repository public.

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

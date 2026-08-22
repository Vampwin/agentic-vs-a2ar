# What must be cited, and under what terms

Audited 22 August 2026 against the tools and data this work actually uses, not
against a generic list. Every entry was traced from an import statement, a
subprocess call, or a file in the deposit.

**Two things are separable and are kept separate here.** A *citation*
obligation is academic — it is what the authors ask for in return for the work.
A *licence* obligation is legal — it constrains what may be redistributed and
under what terms. Most entries below carry the first; only a few carry the
second, and those are the ones that decide this repository's licence.

---

## 1. The licence question, answered first

**Nothing copyleft is linked into this code.** GROMACS (LGPL-2.1), AmberTools
(GPL-3.0-or-later in parts), Open Babel (GPL-2.0) and pdb2pqr are all invoked
as separate executables through `subprocess`, never imported. Running a GPL
program from your own program does not make your program a derivative of it,
so no copyleft term propagates to the code in `pipeline/`.

**The only in-process imports are permissive**, with one to confirm:

| Imported | Licence | Effect on this repository |
|---|---|---|
| RDKit | BSD-3-Clause | none |
| scikit-learn, NumPy, pandas, SciPy, joblib | BSD-3-Clause | none |
| matplotlib | PSF-derived (BSD-compatible) | none |
| PyYAML, pytest | MIT | none |
| AutoDock Vina Python bindings | Apache-2.0 | none; retain the NOTICE if Vina source is ever vendored |
| Meeko | **[CONFIRM]** — believed LGPL-2.1 | if LGPL, importing it is dynamic linking, which LGPL permits for non-LGPL code; confirm before publishing |

**The data question is the one that actually bites.**

| Source | Licence | In this deposit? | Consequence |
|---|---|---|---|
| **LOTUS** v11 | CC BY 4.0 | **Yes** — SMILES, InChIKeys, NPClassifier taxonomy and `organism_name` appear in seven files under `verification/` | Attribution is **mandatory**. CC BY 4.0 is not share-alike, so it does not dictate this repository's licence, but LOTUS must be credited wherever those files are used. |
| **ChEMBL** | CC BY-SA 3.0 | **No** — checked: no `CHEMBL*` identifiers, no `canonical_smiles` column, no bioactivity records. Only aggregate statistics from models trained on it (AUC values, similarity thresholds). | Share-alike is **not triggered** by this deposit. **But if the frozen ChEMBL extract is ever added — for instance because a reviewer asks — that portion must be released under CC BY-SA 3.0, not CC BY 4.0.** |
| DUD-E | see source | No — cited only | citation only |
| PDB 4EIY, 5OLH | CC0 | No — no coordinates deposited | citation only |

**Recommended licensing, and why it is safe:**

- `pipeline/` → **MIT**. No copyleft dependency is linked; every in-process
  import is BSD/MIT/Apache.
- `verification/`, `figures/` → **CC BY 4.0**, with LOTUS attributed in the
  README and in `MANIFEST.json`. Compatible with LOTUS's own CC BY 4.0.
- **Do not** apply CC BY 4.0 to any future addition of ChEMBL records.

---

## 2. Citation obligations, by pipeline stage

### Input curation and standardisation
| Used for | Cite |
|---|---|
| All cheminformatics: parsing, standardisation, descriptors, fingerprints, conformers | **RDKit**, doi:10.5281/zenodo.5548890 |
| Natural-product library source | **LOTUS**, Rutz et al., *eLife* 2022;11:e70780, doi:10.7554/eLife.70780 — and the v11 snapshot, doi:10.5281/zenodo.19360665 |
| Bioactivity records | **ChEMBL**, Mendez et al., *Nucleic Acids Res* 2019;47:D930, doi:10.1093/nar/gky1075 |
| Property and artefact filters | Lipinski doi:10.1016/S0169-409X(96)00423-1 · Veber doi:10.1021/jm020017n · Brenk doi:10.1002/cmdc.200700139 · PAINS doi:10.1021/jm901137j · QED doi:10.1038/nchem.1243 |
| Curation practice | Fourches doi:10.1021/ci100176x and doi:10.1021/acs.jcim.6b00129 |

### Descriptors and similarity
| Used for | Cite |
|---|---|
| ECFP4 / FCFP4 | Rogers & Hahn, doi:10.1021/ci100050t |
| MACCS keys | Durant et al., doi:10.1021/ci010132r |
| ErG | Stiefl et al., doi:10.1021/ci050457y |
| Gobbi 2D pharmacophore | Gobbi & Poppinger, doi:10.1002/(SICI)1097-0290(199824)61:1<47::AID-BIT9>3.0.CO;2-Z |
| Scaffold grouping | Bemis & Murcko, doi:10.1021/jm9602928 |

### Ligand-based modelling
| Used for | Cite |
|---|---|
| All model fitting, CV, calibration | **scikit-learn**, Pedregosa et al., *JMLR* 2011;12:2825 |
| Modelling practice | Tropsha, doi:10.1002/minf.201000061 |
| Y-randomisation | Rücker et al., doi:10.1021/ci700157b |
| Conformal prediction | Norinder doi:10.1021/ci5001168 and doi:10.1016/j.yrtph.2014.12.021 |
| Mondrian conformal | Sun et al., doi:10.1021/acs.jcim.7b00159 |

### Docking
| Used for | Cite |
|---|---|
| Docking engine | **AutoDock Vina** — Trott & Olson doi:10.1002/jcc.21334 **and** Eberhardt et al. doi:10.1021/acs.jcim.1c00203. Apache-2.0. The Vina authors ask that both be cited when 1.2.x is used. |
| Ligand preparation to PDBQT | **Meeko** (Forli lab) — cite the Vina 1.2.0 paper; check the Meeko repository for a preferred citation |
| Conformer generation | ETKDGv3, Riniker & Landrum doi:10.1021/acs.jcim.5b00654 |
| Conformer minimisation | MMFF94, Halgren doi:10.1002/(SICI)1096-987X(199604)17:5/6<490::AID-JCC1>3.0.CO;2-P |
| Receptor protonation | **Open Babel** and **pdb2pqr** — invoked as binaries; both GPL, so cite but no licence propagation |
| Enrichment metrics | BEDROC, Truchon & Bayly doi:10.1021/ci600426e |
| Redocking threshold convention | Warren et al., doi:10.1021/jm050362n |
| Benchmark decoys | DUD-E, Mysinger et al., doi:10.1021/jm300687e |

### Membrane molecular dynamics
| Used for | Cite |
|---|---|
| Simulation engine | **GROMACS**, Abraham et al., *SoftwareX* 2015;1–2:19, doi:10.1016/j.softx.2015.06.001. LGPL-2.1. |
| System construction, parameterisation | **AmberTools** (tleap, antechamber, parmchk2, pdb4amber) — doi:10.1021/acs.jcim.3c01153. GPL-3.0-or-later in parts. **Note: `pmemd` carries a separate, more restrictive licence — this work does not use it.** |
| Bilayer building | **PACKMOL-Memgen**, Schott-Verdugo & Gohlke doi:10.1021/acs.jcim.9b00269, built on PACKMOL |
| Protein force field | ff19SB, Tian et al. doi:10.1021/acs.jctc.9b00591 |
| Lipid force field | Lipid21, Dickson et al. doi:10.1021/acs.jctc.1c01217 |
| Water model | OPC, Izadi et al. doi:10.1021/jz501780a |
| Ligand parameters | GAFF, Wang et al. doi:10.1002/jcc.20035 · AM1-BCC, Jakalian et al. doi:10.1002/jcc.10128 |
| Bilayer choice precedent | Ng, Laughton & Doughty doi:10.1021/ci400463z |
| Receptor structure and its ordered waters | Liu et al., *Science* 2012, doi:10.1126/science.1219218 |

### Verification methodology
| Used for | Cite |
|---|---|
| Seeded-defect testing | Jia & Harman, doi:10.1109/TSE.2010.62 |
| Metamorphic testing | Chen et al., doi:10.1145/3143561 |
| The agent | Anthropic, Claude Code — **[AUTHOR ACTION: state the exact model identifier and version used for the reported runs, per the journal's generative-AI policy]** |

---

## 3. Items still open

1. **Confirm Meeko's licence** before publishing. It is the only in-process
   import whose terms are asserted here rather than verified.
2. **Set the licence files.** MIT in `pipeline/`, CC BY 4.0 for the data
   directories, with LOTUS attribution.
3. **State the model identifier** for the generative-AI declaration.
4. **If ChEMBL records are ever added to this deposit**, that portion becomes
   CC BY-SA 3.0. Do not fold it into a CC BY 4.0 release.

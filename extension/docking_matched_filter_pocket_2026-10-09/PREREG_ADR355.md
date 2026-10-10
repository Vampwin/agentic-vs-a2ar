
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

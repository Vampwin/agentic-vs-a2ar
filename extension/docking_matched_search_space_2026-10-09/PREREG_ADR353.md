
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

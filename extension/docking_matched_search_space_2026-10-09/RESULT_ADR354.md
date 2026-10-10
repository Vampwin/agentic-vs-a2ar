
## ADR-AIIM-354 — Result of ADR-AIIM-353: with the same search space the two engines still disagree, Vina's pocket result depends on its box, and the pocket-cropped DiffDock arm failed its own control (2026-10-09)

**Pre-registration:** ADR-AIIM-353, SHA-256 `7c88d7e4f5a1c8830f46f4333c5f2507295dc8b408f8d50be136ead3f08fad5a`. Parameters, box (centre (8.647, −22.094, 15.858), size (69.5, 110.5, 51.0) Å, the y-centre differing from the registered −22.095 by 0.001 Å) and crop (12.0 Å, 101 residues, 806 atoms) as registered; exhaustiveness 64; DiffDock commit `85c49b6`, same command as ADR-AIIM-317. Both arms ran once, control first, all sixteen compounds, no re-run, about 5 minutes of GPU (arm B) and 12 minutes of CPU (arm A). Outputs: `20_new_runs/docking_matched_search_space_2026-10-09/` (`matched_result.json` `1bd3510d…`, `matched_poses.csv` `cd8edaed…`).

### Grading

| prediction | result | verdict |
|---|---|---|
| H1 Vina blind recovers the crystal pose of ZM241385 | rank-1 at 0.373 Å (score −10.05, centroid 0.15 Å from the site) | **met** |
| H2 DiffDock-L on the cropped protein recovers it | rank-1 at **10.125 Å**; its centroid is 2.21 Å from the site, so it sits in the pocket with a different orientation (confidence −0.07) | **missed** |
| H3 at least 13 of 15 survivors within 8 Å of the site on the cropped protein | 15 of 15 (largest 4.6 Å) | met, and close to a consequence of the crop |
| H4 matched pocket arm: fewer than 8 of 15 within 2 Å | 0 of 15 | met, **not interpreted** (H2 missed) |
| H5 matched blind arm: fewer than 8 of 15 within 2 Å | 0 of 15, median 15.9 Å | **met** |
| H6 median pocket-matched RMSD below 15.3 Å | 11.8 Å | met, **not interpreted** (H2 missed) |

Anchors met: all sixteen poses in every set; the median of the unmatched comparison reads 15.35 Å, reproducing ADR-AIIM-327's 15.3; the archived DiffDock-blind control reads 1.917 Å, reproducing ADR-AIIM-111. **Two defects of mine, fixed before any statement was made:** the crystal SDF has no usable bond orders, so the first run of the grading script returned NaN for every control RMSD and graded H1 and H2 "missed" (a script defect, not a result; bond orders are now assigned from the SMILES); and `g4_dock_library.py` measures `centroid_dist_A` from the centre of its box, which in arm A is the receptor centre, so that column of the arm A table is meaningless for on-pocket status (the analysis recomputes every distance from the poses against the ZM241385 centroid; stated in the README; the frozen script is unchanged).

### What the matched comparisons say

1. **The disagreement is not an artefact of the box.** With both engines blind, 0 of 15 rank-1 poses agree within 2.0 Å (median 15.9 Å), while both recover the control (Vina 0.37 Å, DiffDock 1.92 Å). The unmatched comparison of ADR-AIIM-327 gave 0 of 15 and 15.3 Å.
2. **Vina's pocket result depends on the box.** Run blind, Vina places 2 of 15 survivors on pocket (SNP00223 and SNP00006), against 7 with the 22 Å box; **13 of 15 have their best pose 11.9 to 13.3 Å from the site centroid** and for 13 of 15 none of the ten saved poses is within 8 Å. In the site box, the eight compounds Vina had called off pocket (ADR-AIIM-313) sit 11.9 to 12.7 Å from the centroid, **in the outer part of a box whose half-width is 11 Å along each axis** (further out than any axis-aligned position inside it, so near an edge or a corner, though not necessarily against a face). Blind, the best score is lower than the box score by more than 0.5 kcal/mol for 8 of 15 compounds (SNP00163 −7.70 against −6.83, SNP00198 −7.48 against −6.14, SNP00142 −7.41 against −6.17). Vina's blind poses and its box poses agree within 2.0 Å for 1 of 15 (SNP00006; median 7.0 Å).
3. **For the two compounds taken to molecular dynamics,** the Vina pose in the site box (3.8 and 1.9 Å from the centroid) is not the pose Vina prefers when it may search the whole protein: SNP00157's best blind pose is 11.9 Å away (−7.93 against −7.88, no pose of the ten within 8 Å) and SNP00163's is 12.0 Å away (−7.70 against −6.83). The simulations tested the in-site poses, which is what ADR-AIIM-316 states; they do not test where the compounds would go otherwise, and nothing here says anything about binding.
4. **DiffDock-L does not move to the site when its search space is cropped, only to the neighbourhood**: on the cropped protein all fifteen centroids lie within 4.6 Å of the site, yet the control, which the blind run solved at 1.92 Å, comes out at 10.1 Å. DiffDock blind and pocket-cropped agree within 2.0 Å for 0 of 15 (median 5.9 Å).

### What this does not say, and the consequence for arm B

By the rule of ADR-AIIM-353, **arm B is reported but not interpreted as a pocket-conditioned comparison**: the question "does DiffDock agree with Vina when both are given the pocket" **is not answered**. The radius is not retuned. A proper pocket-conditioned DiffDock needs a separate model (DiffDock-Pocket), whose installation needs the author's approval and which I have **not** installed; I propose it as a separate request. No compound is promoted, demoted or scored on this; DiffDock's sampling noise was not measured, so the 5.9 Å between its two arms includes run-to-run variation.

### Consequences for the documents (proposed, not applied)

Section 3.8 and Figure 6 currently compare a Vina pose from a 22 Å box with a blind DiffDock pose and conclude that the pose of each survivor is engine-dependent and the on-pocket result rests on Vina alone. The matched blind comparison supports that disagreement (0 of 15, both blind) but the on-pocket result of Vina should now say that it depends on the box and that a blind search moves 13 of 15 survivors about 12 Å away. Which wording and whether a panel of the blind result enters Figure 6 are the author's to decide; nothing in the manuscript, supplement, figures or deposit has been changed by this entry. Adding the run to the deposit needs a change of `build_deposit.py`'s tree, which also needs approval.

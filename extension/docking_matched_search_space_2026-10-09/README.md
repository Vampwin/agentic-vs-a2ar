# docking_matched_search_space_2026-10-09 — ADR-AIIM-353 (pre-registration) and ADR-AIIM-354 (result)

Vina and DiffDock-L compared with the same search space, in both directions.

| | whole protein ("blind") | the site ("pocket") |
|---|---|---|
| Vina | **arm A, new**: `vina_blind/` | existing: `snplib_dock_2026-10-06/poses` (22 Å box) |
| DiffDock-L | existing: `diffdock_snplib_2026-10-07/results` | **arm B, new**: `diffdock_pocket/` (protein cropped to 12 Å around the crystal ZM241385 pose) |

- `PREREG_ADR353.md`, `PREREG_SHA256.txt` — the frozen text and its SHA-256 (LF).
- `analysis/build_inputs.py` → `inputs/` (box, ledger, cropped PDB, DiffDock input); `run_vina_blind.sh`; `run_diffdock_pocket.sh`.
- `analysis/matched_analysis.py` → `matched_result.json`, `matched_poses.csv` (grading of H1–H6, anchors, descriptive).
- Vina's frozen `g4_dock_library.py` was used unchanged. **Its `centroid_dist_A` and `on_pocket` columns in `vina_blind/vina_blind_scores.csv` are measured from the centre of the box, which in arm A is the receptor centre, not from the ZM241385 site; ignore them.** The analysis recomputes every distance from the pose files against the ZM241385 centroid.
- The crystal SDF has no usable bond orders; `matched_analysis.py` assigns them from the ligand SMILES before the RMSD to the crystal pose (a first run without this gave NaN, which the script graded as MISSED; that was a defect of the script, fixed, and not a result).
- Absolute `D:\Agentic VS Workflow` paths are in the scripts, as in the other extension runs.
- DiffDock exposes no seed; arm B is a single run.

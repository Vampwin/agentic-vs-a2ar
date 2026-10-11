# Replicate 3 run log (ADR-AIIM-390)

- 2026-10-10 10:53: pre-registration frozen (PREREG_SHA256.txt) and ADR-AIIM-390 appended; both jobs launched with `launch.sh`
  -> `run_r3.sh <ID>` -> frozen `md_10_run_candidate.sh` (unmodified), in WSL distro `A2A-MD`, kept alive by a hidden `sleep infinity`.
- Speed at 44 min: about 2.7 ns per job (about 90 ns/day each), so about 13 h for 50 ns.
- Scripts (SHA-256): adr390_analysis.py 62d7d857bfc3d3d1..., run_r3.sh 0e0c72f85fff2d0..., launch.sh 77ec3ce21b3b7f65...
- Analysis, after both `prod_r3.gro` exist: `run_candidate_md11.sh <ID> prod_r3` per compound (in WSL), then `python adr390_analysis.py --analyse`.
- Host restarts (ADR-AIIM-332 kind) are to be listed here as they occur.
- 2026-10-11 01:30 (SNP00157) and 02:02 (SNP00163): both reached 50 ns on attempt 1, no restart, no host fault (82 and 79 ns/day at the end).
- 2026-10-11 02:02: md_11 (`run_candidate_md11.sh`, via `analyse_md11.sh`) then `adr390_analysis.py`; anchors met. Results SHA-256: adr390_results.json 1a2960e63b5aec0f..., adr390_per_run.csv 5c8a61dfe10f6b46...
- 2026-10-11 post hoc (not pre-registered): `posthoc_pose_clusters.py SNP00157` -> `posthoc_pose_clusters.json`. Pooled-frame clustering of the ligand heavy-atom pose in the receptor frame (cut 3.0 A; 2.0 and 4.0 A as sensitivity), see ADR-AIIM-392.
- 2026-10-11 post hoc (ADR-AIIM-393): posthoc_pose_clusters_all.py for SNP00157, SNP00163 and CONTROL (control topology em_prod.gro, 93,203 atoms); SNP00157 reproduces the ADR-392 clusters exactly.
- 2026-10-11 post hoc (ADR-AIIM-394): run_hbond.sh (gmx hbond, 9 runs, outputs in hbond/), posthoc_rmsf_hbond.py, check_image_jumps.py, plot_rmsf.py; image-jump finding for SNP00163 r2.
- 2026-10-11 post hoc (ADR-AIIM-395): measure_icl3_distance.py -> measure_icl3_distance.json (nearest ligand-to-ICL3-gap distance 33.0 A over all nine runs).

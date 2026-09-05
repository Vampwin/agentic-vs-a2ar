#!/bin/bash
# md_13_analyse_all_candidates.sh -- pose persistence for all nine candidate
# replicates, using the analysis whose three choices were fixed on the control
# before any candidate run existed (md_11_analyse_pose_persistence.sh).
#
# Fit on receptor C-alpha, measure the ligand, reference the pre-equilibration
# structure that holds the docked pose, and report a ladder of thresholds rather
# than one pass/fail. None of those is re-decided here.
#
# A file on disk, not an inline command: four background jobs in this project
# were passed through Windows -> wsl -> bash and had every shell variable eaten
# by the escaping layers (ADR-AIIM-148).
set -u
cd "/mnt/d/Agentic VS Workflow" || exit 1
BASE="20_new_runs/md_candidates_2026-08-22"
ANALYSE="10_new_system/scripts/md_11_analyse_pose_persistence.sh"

for CID in NPLIB01185 NPLIB01643 NPLIB03056; do
  for R in 1 2 3; do
    echo "=== $CID prod_r$R"
    bash "$ANALYSE" "$BASE/$CID/04_run" LIG "prod_r$R" "${CID}_r${R}"
  done
done
echo "### all nine replicates analysed"

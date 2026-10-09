#!/bin/bash
# ADR-AIIM-331 -- mechanical pre-processing of replicates that are already complete, before the deadline.
# Runs the frozen md_11 only; output goes to log files that are NOT read until --analyse is run.
set -u
D="/mnt/d/Agentic VS Workflow/20_new_runs/md_snplib_duplicate_2026-10-08"
for job in "SNP00157 prod_r1" "SNP00157 prod_r2" "SNP00163 prod_r1"; do
  set -- $job
  echo "[$(date '+%F %T')] start $1 $2" >> "$D/precompute_logs/_progress.log"
  bash "$D/analysis/run_candidate_md11.sh" "$1" "$2" > "$D/precompute_logs/$1_$2.log" 2>&1
  echo "[$(date '+%F %T')] done  $1 $2 exit=$?" >> "$D/precompute_logs/_progress.log"
done
echo "[$(date '+%F %T')] ALL DONE" >> "$D/precompute_logs/_progress.log"

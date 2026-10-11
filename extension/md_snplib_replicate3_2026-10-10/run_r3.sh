#!/bin/bash
# ADR-AIIM-390: replicate 3 of one compound. The frozen driver md_10_run_candidate.sh is invoked unmodified; replicates 1 and 2
# are complete and skipped, replicate 3 resumes from its checkpoint. The driver is re-invoked when mdrun dies (the host's GPU
# passthrough has faulted before, ADR-AIIM-332), and every attempt is logged.
set -u
CID="${1:?usage: run_r3.sh <COMPOUND_ID>}"
R="/mnt/d/Agentic VS Workflow"
S="$R/10_new_system/scripts"
B="$R/20_new_runs/md_candidates_2026-08-22/$CID/04_run"
LOG="$R/20_new_runs/md_snplib_replicate3_2026-10-10/run_r3_$CID.log"
log() { echo "[$(date '+%F %T')] $CID $*" | tee -a "$LOG"; }
cd "$B" || exit 1
ATTEMPT=0
while [ ! -s prod_r3.gro ] && [ $ATTEMPT -lt 40 ]; do
  ATTEMPT=$((ATTEMPT + 1))
  log "md_10 attempt $ATTEMPT"
  bash "$S/md_10_run_candidate.sh" "$CID" >> "$LOG" 2>&1
  log "md_10 attempt $ATTEMPT exit=$?"
  [ -s prod_r3.gro ] || sleep 60
done
if [ -s prod_r3.gro ]; then
  log "$CID replicate 3 COMPLETE after $ATTEMPT attempt(s)"
  echo "DONE" > "$R/20_new_runs/md_snplib_replicate3_2026-10-10/_${CID}_r3_done"
else
  log "gave up after $ATTEMPT attempts"; exit 1
fi

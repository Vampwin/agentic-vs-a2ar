#!/bin/bash
# ADR-AIIM-331 -- SNP00163: wait for the build, set up, then equilibrate and run two
# 50 ns replicates through md_10_run_candidate.sh (frozen, unmodified), retried
# through host GPU faults. stop_after_r2.sh ends it before replicate 3.
# Every wait has a timeout.
set -u
B="/mnt/d/Agentic VS Workflow/20_new_runs/md_candidates_2026-08-22/SNP00163"
S="/mnt/d/Agentic VS Workflow/10_new_system/scripts"
LOG="$B/run_snp00163.log"
log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }
set +u; source /etc/profile.d/a2a-cuda.sh 2>/dev/null; source /opt/gromacs/2026.3-cuda13.1-sm120/bin/GMXRC 2>/dev/null; set -u

log "waiting for the build (gate: 'build complete' in build_system.log, up to 5 h)"
for i in $(seq 1 300); do
  grep -q "build complete" "$B/build_system.log" 2>/dev/null && break
  sleep 60
done
grep -q "build complete" "$B/build_system.log" || { log "build did not complete in 5 h"; exit 1; }
log "build complete"

cd "$B/04_run" 2>/dev/null || { log "04_run missing: md_09 has not been run"; exit 1; }
ATTEMPT=0
while [ ! -s prod_r2.gro ] && [ $ATTEMPT -lt 40 ]; do
  ATTEMPT=$((ATTEMPT + 1))
  log "md_10 attempt $ATTEMPT"
  bash "$S/md_10_run_candidate.sh" SNP00163 >> "$LOG" 2>&1
  log "md_10 attempt $ATTEMPT exit=$?"
  [ -s prod_r2.gro ] || sleep 60
done
if [ -s prod_r2.gro ]; then log "SNP00163 two replicates COMPLETE"; echo DONE > "$B/_snp00163_done"
else log "gave up after $ATTEMPT attempts"; exit 1; fi

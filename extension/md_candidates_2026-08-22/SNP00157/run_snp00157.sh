#!/bin/bash
# ADR-AIIM-316 / 319 -- SNP00157 production, queued behind the control and DiffDock.
#
# md_10_run_candidate.sh is the FROZEN driver, invoked unmodified. It is
# resumable (a stage whose .gro exists is skipped, an interrupted mdrun
# continues from its checkpoint) and exits non-zero when mdrun dies. The
# production mdp files were rewritten to 25,000,000 steps (50 ns) in this
# directory; md_09 writes 100 ns and is not edited.
#
# Two lessons from 2026-10-06/07 are built in:
#   1. The host's WSL GPU passthrough faults and kills GPU processes silently,
#      so md_10 is re-invoked until the last replicate's .gro exists.
#   2. Waiting on a sentinel that never appears idled the GPU for 13 hours, so
#      every wait here has a timeout.
set -u
B="/mnt/d/Agentic VS Workflow/20_new_runs/md_candidates_2026-08-22/SNP00157"
S="/mnt/d/Agentic VS Workflow/10_new_system/scripts"
R3="/mnt/d/Agentic VS Workflow/20_new_runs/md_control_replicates_2026-10-06/prod_r3.gro"
DD="/mnt/d/Agentic VS Workflow/20_new_runs/diffdock_snplib_2026-10-07/_diffdock_done"
LOG="$B/run_snp00157.log"
log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }

cd "$B/04_run" || exit 1

# ---- gate 1: the control's third replicate (up to 30 h) --------------------
log "waiting for control replicate 3 (gate: prod_r3.gro)"
for i in $(seq 1 360); do
  [ -s "$R3" ] && break
  sleep 300
done
if [ -s "$R3" ]; then log "control replicate 3 complete"
else log "control replicate 3 still absent after 30 h; proceeding anyway"; fi

# ---- gate 2: DiffDock, which needs the whole card (up to 100 min) ----------
log "waiting for the ADR-AIIM-317 DiffDock job (gate: _diffdock_done)"
for i in $(seq 1 20); do
  [ -f "$DD" ] && break
  sleep 300
done
if [ -f "$DD" ]; then log "DiffDock finished"
else log "DiffDock sentinel absent after 100 min; proceeding so the GPU is not left idle"; fi
sleep 60

# ---- production, retried through GPU faults --------------------------------
ATTEMPT=0
while [ ! -s prod_r3.gro ] && [ $ATTEMPT -lt 40 ]; do
  ATTEMPT=$((ATTEMPT + 1))
  log "md_10 attempt $ATTEMPT"
  bash "$S/md_10_run_candidate.sh" SNP00157 >> "$LOG" 2>&1
  log "md_10 attempt $ATTEMPT exit=$?"
  [ -s prod_r3.gro ] || sleep 60
done

if [ -s prod_r3.gro ]; then
  log "SNP00157 COMPLETE after $ATTEMPT attempt(s)"
  echo "DONE" > "$B/_snp00157_done"
else
  log "gave up after $ATTEMPT attempts; the GPU fault is not transient"
  exit 1
fi

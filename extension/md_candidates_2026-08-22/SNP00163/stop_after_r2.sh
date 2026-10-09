#!/bin/bash
# ADR-AIIM-331 -- stop a candidate after its second 50 ns replicate.
# md_10_run_candidate.sh is FROZEN and loops over three replicates; this
# watcher ends it as soon as prod_r2.gro exists, before the third starts.
# Usage: stop_after_r2.sh <COMPOUND_ID> <driver-script-name-to-kill-first>
set -u
CID="$1"; DRV="$2"
B="/mnt/d/Agentic VS Workflow/20_new_runs/md_candidates_2026-08-22/$CID"
LOG="$B/stop_after_r2.log"
log() { echo "[$(date '+%F %T')] $CID $*" >> "$LOG"; }
cd "$B/04_run" || exit 1
log "watching for prod_r2.gro (timeout 60 h)"
for i in $(seq 1 21600); do
  [ -s prod_r2.gro ] && break
  sleep 10
done
if [ ! -s prod_r2.gro ]; then log "timeout without prod_r2.gro; stopping nothing"; exit 1; fi
log "prod_r2.gro present; stopping the driver chain before replicate 3"
# pkill -f would also match this script, whose command line names the driver
for pat in "$DRV" "md_10_run_candidate.sh $CID"; do
  for p in $(pgrep -f "$pat"); do [ "$p" != "$$" ] && [ "$p" != "$PPID" ] && kill "$p" 2>/dev/null; done
done
sleep 1
for p in $(pgrep -f "deffnm prod_r3"); do
  if [ "$(readlink /proc/$p/cwd 2>/dev/null)" = "$B/04_run" ]; then kill "$p"; log "killed mdrun prod_r3 pid $p"; fi
done
sleep 3
log "remaining processes for $CID: $(pgrep -af "$CID" | grep -v stop_after | wc -l)"
log "files from a started replicate 3, if any: $(ls prod_r3.* 2>/dev/null | tr '\n' ' ')"

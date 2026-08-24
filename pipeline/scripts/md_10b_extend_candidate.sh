#!/bin/bash
# md_10b_extend_candidate.sh <COMPOUND_ID>
#
# Extends a candidate's three production replicates from an already-completed
# 50 ns to 100 ns, continuing from each one's own checkpoint rather than
# rerunning from scratch. Written for NPLIB01185, whose replicates were built
# and (for r1/r2) completed before the target duration was raised from 50 to
# 100 ns (ADR-AIIM-118). NPLIB01643 and NPLIB03056 never had a 50 ns
# checkpoint to extend -- their production mdp files were rewritten to 100 ns
# before their equilibration reached that stage, so they run the ordinary
# driver (md_10_run_candidate.sh) once instead of going through this script.
#
# `gmx convert-tpr -extend` followed by `mdrun -cpi` on the SAME deffnm is the
# documented way to lengthen a finished run: the tpr is rewritten in place so
# its checksum still matches the checkpoint, and mdrun appends to the existing
# trajectory rather than starting a new one. A run still short of its
# original target is left alone here -- it will reach that target on its own,
# and this script is meant to run once the standard driver reports the
# candidate's three 50 ns replicates complete.
set -u
CID="${1:?usage: md_10b_extend_candidate.sh <COMPOUND_ID>}"
EXTEND_PS=50000   # +50 ns, to reach 100 ns from the 50 ns already run

G=/opt/gromacs/2026.3-cuda13.1-sm120/bin/gmx
R="/mnt/d/Agentic VS Workflow"
B="$R/20_new_runs/md_candidates_2026-08-22/$CID"
cd "$B/04_run" || exit 1
mkdir -p "$B/logs"

log() { echo "[$(date '+%F %T')] $CID (100ns extend) $*"; }

for r in 1 2 3; do
  d="prod_r${r}"
  marker="${d}.extended_to_100ns"
  if [ -f "$marker" ]; then
    log "${d}: already extended"; continue
  fi
  if [ ! -f "${d}.gro" ]; then
    log "${d}: FATAL, no ${d}.gro -- the 50 ns run this script extends has not finished"
    exit 1
  fi
  # Idempotent against a re-run that was interrupted mid-extension: if the
  # tpr's own step count already reflects 100 ns, skip straight to mdrun.
  nsteps_now=$("$G" dump -s "${d}.tpr" 2>/dev/null | awk -F= '/^ *nsteps/{gsub(/ /,"",$2); print $2; exit}')
  if [ "${nsteps_now:-0}" -lt 50000000 ]; then
    log "${d}: convert-tpr -extend ${EXTEND_PS} ps"
    "$G" convert-tpr -s "${d}.tpr" -extend "$EXTEND_PS" -o "${d}.tpr" \
      > "$B/logs/convert_tpr_${d}.log" 2>&1 || {
        log "${d}: convert-tpr FAILED"; tail -15 "$B/logs/convert_tpr_${d}.log"; exit 1; }
  else
    log "${d}: tpr already extended, resuming"
  fi
  cpi=""; [ -f "${d}.cpt" ] && cpi="-cpi ${d}.cpt"
  log "${d}: mdrun (resuming to 100 ns)"
  "$G" mdrun -deffnm "${d}" ${cpi} -ntmpi 1 -ntomp 8 \
    -nb gpu -pme gpu -bonded gpu -cpt 15 >> "$B/logs/mdrun_${d}_extend.log" 2>&1
  # After extension the trajectory should cover 100 ns / 2000 frames at the
  # 50 ps/frame the production mdp writes; check the tpr's own final step
  # count rather than trusting the run to have gotten there.
  nsteps_final=$("$G" dump -s "${d}.tpr" 2>/dev/null | awk -F= '/^ *nsteps/{gsub(/ /,"",$2); print $2; exit}')
  if [ "${nsteps_final:-0}" -lt 50000000 ]; then
    log "${d}: mdrun FAILED to reach 100 ns"; tail -20 "$B/logs/mdrun_${d}_extend.log"; exit 1
  fi
  touch "$marker"
  log "${d}: done  $(grep -E '^Performance:' "$B/logs/mdrun_${d}_extend.log" | tail -1)"
done

log "ALL THREE REPLICATES NOW AT 100 ns"

#!/bin/bash
# md_10_run_candidate.sh <COMPOUND_ID>
#
# Minimise, equilibrate and run three production replicates for one candidate.
# Same staged restraint release as the control, then three 50 ns runs that
# differ only in their velocity seed.
#
# Resumable, like the control's driver: a stage whose .gro exists is skipped,
# and a stage interrupted mid-flight continues from its checkpoint. This
# project has lost long runs to four host restarts, and nine production runs
# is a lot to lose.
set -u
CID="${1:?usage: md_10_run_candidate.sh <COMPOUND_ID>}"

G=/opt/gromacs/2026.3-cuda13.1-sm120/bin/gmx
R="/mnt/d/Agentic VS Workflow"
B="$R/20_new_runs/md_candidates_2026-08-22/$CID"
TOP="$B/03_build/system_gmx.top"
NDX="$B/04_run/coupling.ndx"
cd "$B/04_run" || exit 1
mkdir -p "$B/logs"

log() { echo "[$(date '+%F %T')] $CID $*"; }

# ---- minimisation --------------------------------------------------------
if [ ! -f em.gro ]; then
  log "em: grompp"
  "$G" grompp -f em.mdp -c "$B/03_build/system_gmx.gro" \
    -r "$B/03_build/system_gmx.gro" -p "$TOP" -n "$NDX" \
    -o em.tpr -maxwarn 5 > "$B/logs/grompp_em.log" 2>&1 || {
      log "em: grompp FAILED"; tail -15 "$B/logs/grompp_em.log"; exit 1; }
  log "em: mdrun"
  "$G" mdrun -deffnm em -ntmpi 1 -ntomp 8 -nb gpu > "$B/logs/em.log" 2>&1
  fmax=$(grep -A3 "Steepest Descents" "$B/logs/em.log" | grep "Maximum force" | tail -1)
  log "em: done  $fmax"
  case "$fmax" in *inf*) log "em: FATAL, force is infinite -- packing failed"; exit 1;; esac
fi

# ---- staged equilibration ------------------------------------------------
prev="em"
for stage in nvt npt1 npt2 npt3; do
  if [ -f "${stage}.gro" ]; then prev="$stage"; continue; fi
  if [ ! -f "${stage}.tpr" ]; then
    "$G" grompp -f "${stage}.mdp" -c "${prev}.gro" -r "${prev}.gro" \
      -p "$TOP" -n "$NDX" -o "${stage}.tpr" -maxwarn 3 \
      > "$B/logs/grompp_${stage}.log" 2>&1 || {
        log "${stage}: grompp FAILED"; tail -15 "$B/logs/grompp_${stage}.log"; exit 1; }
  fi
  cpi=""; [ -f "${stage}.cpt" ] && cpi="-cpi ${stage}.cpt"
  log "${stage}: mdrun"
  "$G" mdrun -deffnm "${stage}" ${cpi} -ntmpi 1 -ntomp 8 \
    -nb gpu -pme gpu -bonded gpu -cpt 15 >> "$B/logs/mdrun_${stage}.log" 2>&1
  [ -f "${stage}.gro" ] || { log "${stage}: mdrun FAILED"; tail -20 "$B/logs/mdrun_${stage}.log"; exit 1; }
  log "${stage}: done  $(grep -E '^Performance:' "$B/logs/mdrun_${stage}.log" | tail -1)"
  prev="$stage"
done

# ---- three production replicates -----------------------------------------
# All three start from the same equilibrated structure (npt3.gro) with fresh
# velocities. Same geometry, different thermal noise -- which is the only
# thing a replicate is supposed to vary.
for r in 1 2 3; do
  d="prod_r${r}"
  if [ -f "${d}.gro" ]; then log "${d}: already complete"; continue; fi
  if [ ! -f "${d}.tpr" ]; then
    "$G" grompp -f "${d}.mdp" -c npt3.gro -p "$TOP" -n "$NDX" \
      -o "${d}.tpr" -maxwarn 3 > "$B/logs/grompp_${d}.log" 2>&1 || {
        log "${d}: grompp FAILED"; tail -15 "$B/logs/grompp_${d}.log"; exit 1; }
  fi
  cpi=""; [ -f "${d}.cpt" ] && cpi="-cpi ${d}.cpt"
  log "${d}: mdrun (50 ns)"
  "$G" mdrun -deffnm "${d}" ${cpi} -ntmpi 1 -ntomp 8 \
    -nb gpu -pme gpu -bonded gpu -cpt 15 >> "$B/logs/mdrun_${d}.log" 2>&1
  [ -f "${d}.gro" ] || { log "${d}: mdrun FAILED"; tail -20 "$B/logs/mdrun_${d}.log"; exit 1; }
  log "${d}: done  $(grep -E '^Performance:' "$B/logs/mdrun_${d}.log" | tail -1)"
done

log "ALL STAGES COMPLETE (3 replicates x 50 ns)"

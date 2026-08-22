#!/bin/bash
# md_06_run_equilibration.sh -- drive the staged equilibration and production.
#
# Resumable by construction, because this project has already lost long runs to
# four host restarts (ADR-AIIM-089). Each stage is skipped if its .gro exists,
# and a stage interrupted mid-flight is continued from its checkpoint with
# -cpi rather than restarted. Killing this script and running it again loses at
# most the steps since the last checkpoint.
#
# Stages, restraints released gradually rather than at once:
#   nvt    250 ps  NVT   heavy atoms at 1000 kJ/mol/nm^2, velocities generated
#   npt1     1 ns  NPT   still 1000; the bilayer finds its area per lipid
#   npt2     2 ns  NPT   200
#   npt3     2 ns  NPT   50
#   prod   100 ns  NPT   unrestrained
#
# -update gpu is deliberately NOT used: OPC water carries an M-site virtual
# site and GROMACS refuses GPU update for virtual sites. Tested, not assumed.
set -u

B="/mnt/d/Agentic VS Workflow/20_new_runs/md_production_2026-08-21"
G=/opt/gromacs/2026.3-cuda13.1-sm120/bin/gmx
TOP="$B/03_build/system_gmx.top"
NDX="$B/04_run/coupling.ndx"
cd "$B/04_run" || exit 1
mkdir -p "$B/logs"

prev="em_prod"
for stage in nvt npt1 npt2 npt3 prod; do
  if [ -f "${stage}.gro" ]; then
    echo "[$(date '+%F %T')] ${stage}: already complete, skipping"
    prev="$stage"
    continue
  fi

  if [ ! -f "${stage}.tpr" ]; then
    echo "[$(date '+%F %T')] ${stage}: grompp from ${prev}"
    # -r takes the restraint reference coordinates. For every stage after nvt
    # that is the previous stage's output, so the restraints follow the system
    # as it relaxes instead of dragging it back to the minimised geometry.
    "$G" grompp -f "${stage}.mdp" -c "${prev}.gro" -r "${prev}.gro" \
      -p "$TOP" -n "$NDX" -o "${stage}.tpr" -maxwarn 3 \
      > "$B/logs/grompp_${stage}.log" 2>&1
    rc=$?
    if [ $rc -ne 0 ]; then
      echo "[$(date '+%F %T')] ${stage}: grompp FAILED (exit $rc)"
      tail -20 "$B/logs/grompp_${stage}.log"
      exit 1
    fi
  fi

  cpi=""
  [ -f "${stage}.cpt" ] && cpi="-cpi ${stage}.cpt" && \
    echo "[$(date '+%F %T')] ${stage}: continuing from checkpoint"

  echo "[$(date '+%F %T')] ${stage}: mdrun starting"
  "$G" mdrun -deffnm "${stage}" ${cpi} -ntmpi 1 -ntomp 8 -nb gpu -pme gpu \
    -bonded gpu -cpt 15 >> "$B/logs/mdrun_${stage}.log" 2>&1
  rc=$?
  if [ $rc -ne 0 ] || [ ! -f "${stage}.gro" ]; then
    echo "[$(date '+%F %T')] ${stage}: mdrun FAILED (exit $rc)"
    tail -25 "$B/logs/mdrun_${stage}.log"
    exit 1
  fi

  perf=$(grep -E "^Performance:" "$B/logs/mdrun_${stage}.log" | tail -1)
  echo "[$(date '+%F %T')] ${stage}: done   ${perf}"
  prev="$stage"
done

echo "[$(date '+%F %T')] ALL STAGES COMPLETE"

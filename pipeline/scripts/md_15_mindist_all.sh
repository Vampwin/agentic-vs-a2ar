#!/bin/bash
# md_15_mindist_all.sh -- minimum ligand-receptor distance for all nine
# candidate replicates.
#
# The RMSD-to-docked-pose said NPLIB01185 reached 60 A and its mean was 24.7 A.
# The time series showed that value oscillating between about 5 and 60 A
# hundreds of times, which no molecule does, and `gmx mindist` -- which applies
# the minimum-image convention itself instead of trusting how the trajectory was
# written -- puts every frame of that same replicate within 4.88 A of the
# receptor. The ligand never left. The RMSD was reading a periodic image.
#
# This measurement is the one that answers "is it still bound". It is computed
# from the raw .xtc, so no trjconv choice can corrupt it.
set -u
cd "/mnt/d/Agentic VS Workflow" || exit 1
G=/opt/gromacs/2026.3-cuda13.1-sm120/bin/gmx
BASE="20_new_runs/md_candidates_2026-08-22"

for CID in NPLIB01185 NPLIB01643 NPLIB03056; do
  for R in 1 2 3; do
    D="$BASE/$CID/04_run/analysis_prod_r${R}"
    RUN="$BASE/$CID/04_run"
    [ -f "$D/analysis.ndx" ] || { echo "$CID r$R: no index, skipped"; continue; }
    if [ ! -f "$D/mindist.xvg" ]; then
      printf 'Calpha\nLigand\n' | "$G" mindist \
        -s "$RUN/prod_r${R}.tpr" -f "$RUN/prod_r${R}.xtc" \
        -n "$D/analysis.ndx" -od "$D/mindist.xvg" -tu ns \
        > "$D/mindist.log" 2>&1
    fi
    [ -f "$D/mindist.xvg" ] || { echo "$CID r$R: mindist FAILED"; continue; }
    echo "$CID r$R done"
  done
done
echo "### mindist complete for all replicates"

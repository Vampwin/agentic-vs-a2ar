#!/bin/bash
# POST HOC (ADR-AIIM-394): ligand-protein hydrogen bonds with the GROMACS geometric definition, defaults: D-A <= 0.35 nm, H-D-A angle <= 30 degrees.
# gmx hbond is the new implementation (2024+): -r Ligand, -t Protein, per-frame output (-pf) so occupancy and partners can be counted.
set -u
G=/opt/gromacs/2026.3-cuda13.1-sm120/bin/gmx
R="/mnt/d/Agentic VS Workflow/20_new_runs"
O="$R/md_snplib_replicate3_2026-10-10/hbond"
export GMX_MAXBACKUP=-1
one() {  # label tpr analysis_dir
  local L="$1" T="$2" A="$3" RES="$4"
  [ -s "$O/$L.ndx" ] && { echo "$L: exists"; return; }
  "$G" hbond -s "$T" -f "$A/whole.xtc" -r "resname $RES" -t 'group "Protein"'  -pf -o "$O/$L.ndx" -num "$O/$L.num.xvg" > "$O/$L.log" 2>&1 \
    && echo "$L: done" || { echo "$L: FAILED"; tail -5 "$O/$L.log"; }
}
for C in SNP00157 SNP00163; do
  D="$R/md_candidates_2026-08-22/$C/04_run"
  for r in 1 2 3; do one "${C}_r$r" "$D/prod_r$r.tpr" "$D/analysis_prod_r$r" LIG & done
done
one CONTROL_r1 "$R/md_control_anchor_2026-10-07/prod.tpr" "$R/md_control_anchor_2026-10-07/analysis_prod" ZMA &
one CONTROL_r2 "$R/md_control_replicates_2026-10-06/prod_r2.tpr" "$R/md_control_replicates_2026-10-06/analysis_prod_r2" ZMA &
one CONTROL_r3 "$R/md_control_replicates_2026-10-06/prod_r3.tpr" "$R/md_control_replicates_2026-10-06/analysis_prod_r3" ZMA &
wait
ls -la "$O" | head -30

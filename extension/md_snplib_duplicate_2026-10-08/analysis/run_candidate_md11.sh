#!/bin/bash
# ADR-AIIM-331 -- step 1 of the candidate analysis: the frozen md_11 on one replicate.
#   run_candidate_md11.sh <COMPOUND_ID> <DEFFNM>        e.g. SNP00157 prod_r1
# md_11_analyse_pose_persistence.sh is invoked unmodified, as for the control (ADR-AIIM-324):
# receptor C-alpha fit, ligand RMSD against em.tpr (the pre-equilibration structure),
# bilayer box. Few threads: the host's cores are shared with the running simulations.
set -u
CID="${1:?compound}"; DEFFNM="${2:?deffnm}"
R="/mnt/d/Agentic VS Workflow"
RUN="$R/20_new_runs/md_candidates_2026-08-22/$CID/04_run"
export OMP_NUM_THREADS=2 GMX_MAXBACKUP=-1
echo "[$(date '+%F %T')] $CID $DEFFNM: md_11"
bash "$R/10_new_system/scripts/md_11_analyse_pose_persistence.sh" "$RUN" LIG "$DEFFNM" "$CID-$DEFFNM" \
  || { echo "md_11 FAILED"; exit 1; }
O="$RUN/analysis_${DEFFNM}"
for f in ligand_rmsd.xvg box.xvg mindist.xvg; do [ -s "$O/$f" ] || { echo "FATAL: $f missing"; exit 1; }; done
echo "[$(date '+%F %T')] $CID $DEFFNM: md_11 done ($(grep -vcE '^[#@]' "$O/ligand_rmsd.xvg") RMSD frames)"

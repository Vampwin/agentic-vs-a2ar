#!/bin/bash
# ADR-AIIM-324 -- the original control analysis, applied to one trajectory.
#
#   run_control_analysis.sh <RUN_DIR> <DEFFNM> <LABEL>
#
# RUN_DIR holds <DEFFNM>.xtc/.tpr/.edr and em_prod.tpr (the pre-equilibration
# reference). Steps are the two frozen scripts' own, unchanged:
#   1. md_11_analyse_pose_persistence.sh: receptor C-alpha fit, ligand RMSD against
#      the pre-equilibration structure, bilayer box, all with -pbc mol -ur compact
#   2. the N6.55 minimum-distance series, as md_18_recognition_motif.py writes it,
#      using the control's own atom-index groups. They are valid because the system
#      is identical; group sizes are asserted below rather than assumed.
#
# Few threads: this host's cores are shared with the running molecular dynamics.
set -u
RUN="${1:?run dir}"; DEFFNM="${2:?deffnm}"; LABEL="${3:-$DEFFNM}"
G=/opt/gromacs/2026.3-cuda13.1-sm120/bin/gmx
R="/mnt/d/Agentic VS Workflow"
S="$R/10_new_system/scripts"
NDX="$R/20_new_runs/md_recognition_motif_2026-09-02/control_ZM241385/N6_55.ndx"
export OMP_NUM_THREADS=2
export GMX_MAXBACKUP=-1

echo "[$(date '+%F %T')] $LABEL: md_11"
bash "$S/md_11_analyse_pose_persistence.sh" "$RUN" ZMA "$DEFFNM" "$LABEL" || { echo "md_11 FAILED"; exit 1; }

OUT="$RUN/analysis_${DEFFNM}"
# group sizes of the borrowed index file, asserted
na=$(awk '/^\[ Anchor \]/{f=1;next} /^\[/{f=0} f{n+=NF} END{print n+0}' "$NDX")
nl=$(awk '/^\[ Ligand \]/{f=1;next} /^\[/{f=0} f{n+=NF} END{print n+0}' "$NDX")
echo "[$(date '+%F %T')] $LABEL: N6.55 index groups  Anchor $na atoms  Ligand(polar) $nl atoms (control: 2 and 9)"
[ "$na" = "2" ] && [ "$nl" = "9" ] || { echo "FATAL: index group sizes differ from the control's"; exit 1; }

if [ ! -s "$OUT/N6_55.xvg" ]; then
  printf 'Anchor\nLigand\n' | "$G" mindist -s "$RUN/${DEFFNM}.tpr" -f "$RUN/${DEFFNM}.xtc" \
     -n "$NDX" -od "$OUT/N6_55.xvg" -tu ns > "$OUT/N6_55.log" 2>&1
fi
[ -s "$OUT/N6_55.xvg" ] || { echo "FATAL: N6.55 series not produced"; tail -5 "$OUT/N6_55.log"; exit 1; }
echo "[$(date '+%F %T')] $LABEL: done  ($(grep -vcE '^[#@]' "$OUT/N6_55.xvg") N6.55 frames, $(grep -vcE '^[#@]' "$OUT/ligand_rmsd.xvg") RMSD frames)"

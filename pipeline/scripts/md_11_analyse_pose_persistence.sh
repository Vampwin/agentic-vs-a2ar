#!/bin/bash
# md_11_analyse_pose_persistence.sh <RUN_DIR> <LIGAND_RESNAME> <DEFFNM> [LABEL]
#
# Does the ligand stay where it was placed?
#
# Three choices decide whether the answer means anything, and all three were
# fixed on the control's partial trajectory BEFORE any candidate run existed,
# so that no threshold here is chosen after seeing the data it judges.
#
# 1. **Fit on the receptor, measure the ligand.** Without fitting, the number
#    is dominated by the whole system tumbling and diffusing in the box, which
#    says nothing about binding. Fitted, it is how far the ligand moved
#    relative to the pocket.
#
# 2. **Reference the PRE-equilibration structure, not the production start.**
#    `prod.tpr` is built from `npt3.gro`, so using it measures drift since
#    equilibration -- a different and more flattering question than "was the
#    docked/crystallographic pose reproduced". On the control this understated
#    deviation by ~0.3 A (mean 2.31 vs 2.64 A). The reference here is `em*.tpr`,
#    built from the system as packed, which holds the crystallographic pose for
#    the control and the docked pose for a candidate. A `.gro` cannot be used:
#    it carries no masses and gmx rms refuses it.
#
# 3. **Report a ladder of thresholds, not one pass/fail.** The 2.0 A criterion
#    belongs to redocking (Section 2.5), where it compares a predicted pose to
#    a crystal structure solved at cryogenic temperature. Applied to a 310 K
#    trajectory it is far too strict: on the control only 17% of frames sit
#    below 2.0 A, while 71% are below 3.0 A and **100% below 4.0 A** -- a
#    ligand fluctuating in its pocket, not one leaving it. Reducing that to
#    "17% retained" would have been actively misleading. The ladder is printed
#    and the reader draws the line.
#
# No MDAnalysis or mdtraj in this environment, and installing either needs
# approval, so this is gmx plus awk.
set -u
RUN="${1:?usage: md_11_analyse_pose_persistence.sh <RUN_DIR> <LIGRES> <DEFFNM> [LABEL]}"
LIGRES="${2:?ligand residue name, e.g. ZMA or LIG}"
DEFFNM="${3:?trajectory basename, e.g. prod or prod_r1}"
LABEL="${4:-$DEFFNM}"

G=/opt/gromacs/2026.3-cuda13.1-sm120/bin/gmx
# Resolve the script's own directory BEFORE changing into the run directory.
# $0 is a relative path, so `dirname "$0"` after the cd pointed at
# <run>/10_new_system/scripts, which does not exist: the bilayer check silently
# failed on all nine candidate replicates on 2026-08-31 while the RMSD numbers
# beside it were fine. Found because the error was printed, not because
# anything checked that the check had run.
SCRIPTDIR="$(cd "$(dirname "$0")" && pwd)"
AWKDIR="$SCRIPTDIR"
cd "$RUN" || exit 1
OUT="analysis_${DEFFNM}"
mkdir -p "$OUT"

[ -f "${DEFFNM}.xtc" ] || { echo "  no ${DEFFNM}.xtc in $RUN"; exit 1; }
[ -f "${DEFFNM}.tpr" ] || { echo "  no ${DEFFNM}.tpr in $RUN"; exit 1; }

# The pre-equilibration tpr is the reference. Named em_prod.tpr for the
# control, em.tpr for a candidate.
REF=""
for c in em_prod.tpr em.tpr; do [ -f "$c" ] && REF="$c" && break; done
[ -n "$REF" ] || { echo "  FATAL: no em tpr to reference the original pose against"; exit 1; }

# --- index: receptor C-alpha to fit on, ligand to measure ------------------
cat > "$OUT/ndx.in" <<NDXIN
keep 0
r ${LIGRES}
name 1 Ligand
a CA
name 2 Calpha
q
NDXIN
"$G" make_ndx -f "${DEFFNM}.tpr" -o "$OUT/analysis.ndx" \
  < "$OUT/ndx.in" > "$OUT/make_ndx.log" 2>&1

nlig=$(awk '/^\[ Ligand \]/{f=1;next} /^\[/{f=0} f{n+=NF} END{print n+0}' "$OUT/analysis.ndx")
nca=$(awk '/^\[ Calpha \]/{f=1;next} /^\[/{f=0} f{n+=NF} END{print n+0}' "$OUT/analysis.ndx")
[ "$nlig" -gt 0 ] || { echo "  FATAL: ligand group '${LIGRES}' is empty"; exit 1; }
[ "$nca" -gt 0 ] || { echo "  FATAL: no C-alpha atoms to fit on"; exit 1; }

# --- PBC: make molecules whole, compact cell, centre on the receptor -------
# Order and flags both matter, and getting them wrong is not visible in any
# summary statistic.
#
# `-pbc mol -center` alone writes a bound ligand into whichever image the
# default rectangular cell puts it in. On NPLIB01185 that produced an RMSD
# oscillating between about 5 and 60 A hundreds of times across 100 ns. It was
# not the ligand moving: in the 400 frames whose RMSD exceeded 40 A, gmx
# mindist -- which applies the minimum-image convention itself -- put the
# ligand 3.03 to 4.68 A from the receptor, in van der Waals contact, exactly as
# in the frames whose RMSD read under 8 A. A molecule cannot be sixty angstrom
# away and touching at the same time. `-ur compact` is what keeps it adjacent.
#
# `-pbc nojump` was tried and is wrong here: it preserves the unwrapped path, so
# coordinates stay outside the box and `-ur compact` can no longer bring them
# back. The two treatments fight each other. nojump belongs to a genuinely
# dissociating ligand; this is a bound one.
# `-pbc cluster` was tried and refused: it requires the clustering group to
# contain whole molecules, and C-alpha atoms are not molecules.
if [ ! -f "$OUT/whole.xtc" ]; then
  printf 'Calpha\nSystem\n' | "$G" trjconv -s "${DEFFNM}.tpr" -f "${DEFFNM}.xtc" \
    -n "$OUT/analysis.ndx" -o "$OUT/whole.xtc" -pbc mol -ur compact -center \
    > "$OUT/trjconv.log" 2>&1
fi
[ -f "$OUT/whole.xtc" ] || { echo "  FATAL: trjconv produced nothing"; tail -5 "$OUT/trjconv.log"; exit 1; }

# --- minimum ligand-receptor distance, which PBC cannot confuse ------------
# An RMSD says how far the ligand is from where it was placed. This says
# whether it is in contact with the protein at all, and gmx mindist applies
# the minimum-image convention itself rather than trusting how the trajectory
# was written. It is the measurement that decides bound or not; the RMSD
# describes the pose while it is bound.
if [ ! -f "$OUT/mindist.xvg" ]; then
  printf 'Calpha\nLigand\n' | "$G" mindist -s "${DEFFNM}.tpr" -f "${DEFFNM}.xtc" \
    -n "$OUT/analysis.ndx" -od "$OUT/mindist.xvg" -tu ns \
    > "$OUT/mindist.log" 2>&1
fi
[ -f "$OUT/whole.xtc" ] || { echo "  FATAL: trjconv produced nothing"; tail -5 "$OUT/trjconv.log"; exit 1; }

# --- ligand RMSD against the pose as placed --------------------------------
printf 'Calpha\nLigand\n' | "$G" rms -s "$REF" -f "$OUT/whole.xtc" \
  -n "$OUT/analysis.ndx" -o "$OUT/ligand_rmsd.xvg" -tu ns > "$OUT/rms.log" 2>&1
[ -f "$OUT/ligand_rmsd.xvg" ] || { echo "  FATAL: gmx rms produced nothing"; tail -8 "$OUT/rms.log"; exit 1; }

# --- box dimensions, an equilibration check independent of the ligand ------
printf 'Box-X\nBox-Y\n' | "$G" energy -f "${DEFFNM}.edr" -o "$OUT/box.xvg" \
  > "$OUT/energy.log" 2>&1

# --- summarise --------------------------------------------------------------
echo "  ${LABEL}: ligand ${nlig} atoms, fitted on ${nca} C-alpha, ref ${REF}"
awk -v label="$LABEL" '
  /^[@#]/ {next}
  { t=$1; r=$2*10               # gmx reports nm; Angstrom is the convention here
    n++; sum+=r; sq+=r*r
    if (r>max) {max=r; tmax=t}
    if (r<=2.0) b2++; if (r<=3.0) b3++; if (r<=4.0) b4++; if (r<=5.0) b5++
    last=r }
  END {
    if (n==0) {print "    no frames parsed"; exit 1}
    mean=sum/n; sd=sqrt(sq/n - mean*mean)
    printf "    frames %4d over %.1f ns\n", n, t
    printf "    RMSD  mean %.2f +/- %.2f A   max %.2f A (at %.1f ns)   final %.2f A\n",
           mean, sd, max, tmax, last
    printf "    fraction of frames below:  2.0 A %3.0f%%   3.0 A %3.0f%%   4.0 A %3.0f%%   5.0 A %3.0f%%\n",
           100*b2/n, 100*b3/n, 100*b4/n, 100*b5/n
  }' "$OUT/ligand_rmsd.xvg"

# The bilayer check lives in its own awk file: shell heredocs in this project
# have repeatedly mangled backslash escapes, and a silently broken analysis is
# worse than none.
true  # AWKDIR resolved before the cd, see SCRIPTDIR above
if [ -f "$OUT/box.xvg" ]; then
  awk -f "$AWKDIR/md_box_trend.awk" "$OUT/box.xvg"
fi

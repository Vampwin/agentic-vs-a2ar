#!/bin/bash
# md_16_pocket_distance.sh -- did the ligand stay in the orthosteric pocket?
#
# `gmx mindist` to the whole receptor says the ligand is always within about
# 5 A of *some* C-alpha, which a ligand that had migrated across a 50 A protein
# surface would also satisfy. It cannot distinguish "still in the pocket" from
# "somewhere else on the protein", and the RMSD-to-docked-pose cannot be trusted
# to make that call either, because it reported 60 A excursions on a trajectory
# whose ligand never left contact.
#
# This measures the distance from the ligand to Asn253, the orthosteric residue
# the whole study is built around. It is computed from the raw .xtc with the
# minimum-image convention applied by gmx itself, so no trjconv choice can
# corrupt it, and it answers the question the other two cannot.
set -u
cd "/mnt/d/Agentic VS Workflow" || exit 1
G=/opt/gromacs/2026.3-cuda13.1-sm120/bin/gmx
BASE="20_new_runs/md_candidates_2026-08-22"
POCKET_RESID=253

for CID in NPLIB01185 NPLIB01643 NPLIB03056; do
  for R in 1 2 3; do
    RUN="$BASE/$CID/04_run"
    D="$RUN/analysis_prod_r${R}"
    [ -f "$RUN/prod_r${R}.tpr" ] || { echo "$CID r$R: no tpr"; continue; }
    mkdir -p "$D"
    if [ ! -f "$D/pocket.ndx" ]; then
      printf 'keep 0\nr LIG\nname 1 Ligand\nr %s & a CA\nname 2 PocketCA\nq\n' \
        "$POCKET_RESID" \
        | "$G" make_ndx -f "$RUN/prod_r${R}.tpr" -o "$D/pocket.ndx" \
        > "$D/pocket_ndx.log" 2>&1
    fi
    nlig=$(awk '/^\[ Ligand \]/{f=1;next} /^\[/{f=0} f{n+=NF} END{print n+0}' "$D/pocket.ndx")
    npoc=$(awk '/^\[ PocketCA \]/{f=1;next} /^\[/{f=0} f{n+=NF} END{print n+0}' "$D/pocket.ndx")
    if [ "$nlig" -eq 0 ] || [ "$npoc" -eq 0 ]; then
      echo "$CID r$R: empty group (ligand $nlig, pocket $npoc)"; continue
    fi
    if [ ! -f "$D/pocket_dist.xvg" ]; then
      printf 'PocketCA\nLigand\n' | "$G" mindist \
        -s "$RUN/prod_r${R}.tpr" -f "$RUN/prod_r${R}.xtc" \
        -n "$D/pocket.ndx" -od "$D/pocket_dist.xvg" -tu ns \
        > "$D/pocket_dist.log" 2>&1
    fi
    [ -f "$D/pocket_dist.xvg" ] || { echo "$CID r$R: mindist failed"; continue; }
    echo "$CID r$R ok  (ligand $nlig atoms, pocket $npoc atoms)"
  done
done
echo "### pocket distances complete"

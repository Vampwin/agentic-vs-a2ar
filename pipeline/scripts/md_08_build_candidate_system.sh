#!/bin/bash
# md_08_build_candidate_system.sh <COMPOUND_ID>
#
# Membrane build for one Tier-A candidate. Every packmol-memgen and tleap
# setting is copied from the ZM241385 control (md_03b_build_production.sh, and
# the leap input it produced) so the candidate systems differ from the control
# in the ligand and nothing else. That is what makes the control a control.
#
# Three failures from the control build are pre-empted here rather than
# rediscovered, because each cost real time the first time:
#
#   * --nloop 60 --nloop_all 400. At the defaults this composition left 119
#     inter-residue contacts under 0.7 A after parameterisation and drove
#     minimisation to an infinite force (ADR-AIIM-094).
#   * The NME CH3 atoms are stripped before tleap. tleap's NME template builds
#     its own CH3 and rejects a supplied one with "does not have a type".
#   * addionsrand ... 0 neutralises; the retained crystallographic Na+ is not
#     counted by packmol-memgen's own salt calculation.
set -u
CID="${1:?usage: md_08_build_candidate_system.sh <COMPOUND_ID>}"

A=/home/a2amd/A2A_MD/ambertools26_env
export AMBERHOME="$A"
export PATH="$A/bin:$PATH"

R="/mnt/d/Agentic VS Workflow"
B="$R/20_new_runs/md_candidates_2026-08-22/$CID"
cd "$B/03_build" || exit 1
mkdir -p "$B/logs"

echo "[$(date '+%F %T')] $CID: packmol-memgen"
packmol-memgen \
  --pdb complex.pdb \
  --lipids POPC --ratio 1 \
  --dist 15 --dist_wat 17.5 \
  --salt --saltcon 0.15 --salt_c "Na+" --salt_a "Cl-" \
  --ffprot ff19SB --fflip lipid21 --ffwat opc --gaff2 \
  --notprotonate \
  --keepligs \
  --nloop 60 --nloop_all 400 \
  --noprogress --keep \
  -o system.pdb \
  --log "$B/logs/packmol_memgen.log" > "$B/logs/packmol_stdout.log" 2>&1
echo "[$(date '+%F %T')] $CID: packmol-memgen exit=$? (STOP 173 is expected;"
echo "                     it means 'best solution found', not failure)"
[ -f system.pdb ] || { echo "FATAL: no system.pdb"; exit 1; }

# tleap rejects a supplied NME methyl; let its own template build it.
cp -f system.pdb system_withCH3.pdb
awk '{ n=substr($0,13,4); gsub(/ /,"",n); r=substr($0,18,3);
       if (!(n=="CH3" && r=="NME")) print }' system_withCH3.pdb > system.pdb
echo "[$(date '+%F %T')] $CID: stripped $(grep -c 'CH3 NME' system_withCH3.pdb) NME CH3 atoms"

cat > leap_cand.in <<LEAP
source leaprc.lipid21
source leaprc.protein.ff19SB
source leaprc.water.opc
source leaprc.gaff2

LIG = loadmol2 LIG.mol2
loadamberparams LIG.frcmod

SYS = loadpdb system.pdb
# Neutralise. The retained crystallographic Na+ contributes +1 that
# packmol-memgen's ion placement does not account for.
addionsrand SYS Cl- 0
setBox SYS vdw
savepdb SYS system_lipid.pdb
saveAmberParmNetcdf SYS system_lipid.top system_lipid.crd
quit
LEAP

echo "[$(date '+%F %T')] $CID: tleap"
tleap -f leap_cand.in > "$B/logs/leap.log" 2>&1
echo "                     tleap exit=$?"
grep -i "^FATAL" "$B/logs/leap.log" | head -3
[ -s system_lipid.top ] || { echo "FATAL: tleap produced no topology"; exit 1; }

echo "[$(date '+%F %T')] $CID: converting to GROMACS"
"$A/bin/python" - <<PY
import parmed
p = parmed.load_file("system_lipid.top", xyz="system_lipid.crd")
p.save("system_gmx.top", overwrite=True)
p.save("system_gmx.gro", overwrite=True)
q = sum(a.charge for a in p.atoms)
print(f"                     {len(p.atoms)} atoms, net charge {q:+.4f}")
if abs(q) > 0.01:
    raise SystemExit(f"FATAL: net charge {q:+.4f} is not neutral")
PY
echo "[$(date '+%F %T')] $CID: build complete"

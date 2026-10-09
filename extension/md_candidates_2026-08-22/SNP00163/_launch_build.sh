#!/bin/bash
# ADR-AIIM-331 -- detached membrane build for SNP00163.
#
# md_08_build_candidate_system.sh is the FROZEN script, invoked unmodified. It
# is library-agnostic: it reads 03_build/LIG.mol2, LIG.frcmod and
# ligand_posed.sdf and builds the bilayer with the control's packmol-memgen and
# tleap settings. Mostly CPU, so it runs beside the control's GPU replicates.
#
# Detached because packmol-memgen with --nloop 60 --nloop_all 400 takes a while
# and must outlive this caller.
D="/mnt/d/Agentic VS Workflow/20_new_runs/md_candidates_2026-08-22/SNP00163"
S="/mnt/d/Agentic VS Workflow/10_new_system/scripts"

set +u
source /etc/profile.d/a2a-cuda.sh 2>/dev/null
source /opt/gromacs/2026.3-cuda13.1-sm120/bin/GMXRC 2>/dev/null
export AMBERHOME=/home/a2amd/A2A_MD/ambertools26_env
export PATH=$AMBERHOME/bin:$PATH
set -u

cd "$D" || exit 1
setsid nohup bash "$S/md_08_build_candidate_system.sh" SNP00163 \
    > "$D/build_system.log" 2>&1 < /dev/null &
PID=$!
disown "$PID" 2>/dev/null
sleep 5
echo "md_08 pid $PID"
ps -p "$PID" -o pid=,cmd= 2>/dev/null || echo "WARNING: build not running"

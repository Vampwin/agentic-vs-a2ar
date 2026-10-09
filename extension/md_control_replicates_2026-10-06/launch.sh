#!/bin/bash
# ADR-AIIM-314 launcher. Detaches the replicate driver so it outlives this
# caller: the driver waits hours for ADR-AIIM-305 to release the GPU and then
# runs for over two days, which no foreground shell here can hold open.
#
# setsid + nohup, and a sentinel the driver writes when all replicates finish.
D="/mnt/d/Agentic VS Workflow/20_new_runs/md_control_replicates_2026-10-06"

# GMXRC references variables that may be unset, so `set -u` stays off until
# after it is sourced -- with it on, sourcing GMXRC kills this script silently.
set +u
source /etc/profile.d/a2a-cuda.sh 2>/dev/null
source /opt/gromacs/2026.3-cuda13.1-sm120/bin/GMXRC 2>/dev/null
set -u

cd "$D" || exit 1
setsid nohup bash "$D/run_control_replicates.sh" "$@" > "$D/driver.log" 2>&1 < /dev/null &
PID=$!
disown "$PID" 2>/dev/null
sleep 3
echo "driver pid $PID"
ps -p "$PID" -o pid=,cmd= 2>/dev/null || echo "WARNING: driver not running"

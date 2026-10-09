#!/bin/bash
# Stop the driver and any running mdrun cleanly, then confirm.
# GROMACS honours SIGTERM only at a neighbour-search step, so it is given time
# before SIGKILL. The .log lock is released when the process actually exits.
pkill -f run_control_replicates.sh
sleep 2
PIDS=$(pgrep -f 'gmx mdrun')
for p in $PIDS; do
  echo "sending TERM to $p"
  kill -TERM "$p" 2>/dev/null
done
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24; do
  sleep 5
  pgrep -f 'gmx mdrun' > /dev/null || break
done
PIDS=$(pgrep -f 'gmx mdrun')
for p in $PIDS; do
  echo "TERM ignored by $p, sending KILL"
  kill -KILL "$p" 2>/dev/null
done
sleep 5
if pgrep -f 'gmx mdrun|run_control_replicates' > /dev/null; then
  pgrep -af 'gmx mdrun|run_control_replicates'
else
  echo "all stopped"
fi

#!/bin/bash
D="/mnt/d/Agentic VS Workflow/20_new_runs/diffdock_snplib_2026-10-07"
setsid nohup bash "$D/run_diffdock_snplib.sh" > "$D/diffdock.log" 2>&1 < /dev/null &
PID=$!
disown "$PID" 2>/dev/null
sleep 3
echo "diffdock waiter pid $PID"
ps -p "$PID" -o pid=,cmd= 2>/dev/null || echo "WARNING: not running"

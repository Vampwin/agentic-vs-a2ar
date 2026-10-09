#!/bin/bash
D="/mnt/d/Agentic VS Workflow/20_new_runs/md_candidates_2026-08-22/SNP00157"
setsid nohup bash "$D/run_snp00157.sh" > "$D/driver_snp00157.log" 2>&1 < /dev/null &
PID=$!
disown "$PID" 2>/dev/null
sleep 3
echo "SNP00157 driver pid $PID"
ps -p "$PID" -o pid=,cmd= 2>/dev/null || echo "WARNING: not running"

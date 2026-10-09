#!/bin/bash
D="/mnt/d/Agentic VS Workflow/20_new_runs/md_candidates_2026-08-22/SNP00163"
cd "$D" || exit 1
setsid nohup bash run_snp00163.sh > driver_snp00163.log 2>&1 < /dev/null &
sleep 3; pgrep -af run_snp00163 | head -2

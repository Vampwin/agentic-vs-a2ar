#!/bin/bash
D="/mnt/d/Agentic VS Workflow/20_new_runs/md_candidates_2026-08-22/SNP00163"
cd "$D" || exit 1
setsid nohup bash stop_after_r2.sh SNP00157 run_snp00157.sh > /dev/null 2>&1 < /dev/null &
sleep 2
pgrep -af "stop_after_r2.sh" | cut -c1-100

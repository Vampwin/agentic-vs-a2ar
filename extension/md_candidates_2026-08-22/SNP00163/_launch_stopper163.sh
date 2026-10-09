#!/bin/bash
cd "/mnt/d/Agentic VS Workflow/20_new_runs/md_candidates_2026-08-22/SNP00163" || exit 1
setsid nohup bash stop_after_r2.sh SNP00163 run_snp00163.sh > /dev/null 2>&1 < /dev/null &
sleep 2; pgrep -af stop_after_r2 | cut -c1-100

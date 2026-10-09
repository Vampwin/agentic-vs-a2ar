#!/bin/bash
cd "/mnt/d/Agentic VS Workflow/20_new_runs/md_snplib_duplicate_2026-10-08" || exit 1
setsid nohup bash analysis/precompute.sh > /dev/null 2>&1 < /dev/null &
sleep 3; pgrep -af precompute.sh | head -2

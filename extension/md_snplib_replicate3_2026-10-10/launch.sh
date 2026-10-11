#!/bin/bash
D="/mnt/d/Agentic VS Workflow/20_new_runs/md_snplib_replicate3_2026-10-10"
for CID in SNP00157 SNP00163; do
  setsid nohup bash "$D/run_r3.sh" "$CID" > "$D/driver_r3_$CID.log" 2>&1 < /dev/null &
  disown "$!" 2>/dev/null
  echo "$CID driver pid $!"
done
sleep 5
ps -eo pid,cmd | grep -E "[r]un_r3|[g]mx mdrun" | cut -c1-140

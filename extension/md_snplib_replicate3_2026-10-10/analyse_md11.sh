#!/bin/bash
D="/mnt/d/Agentic VS Workflow/20_new_runs/md_snplib_replicate3_2026-10-10"
S="/mnt/d/Agentic VS Workflow/20_new_runs/md_snplib_duplicate_2026-10-08/analysis/run_candidate_md11.sh"
for c in SNP00157 SNP00163; do
  bash "$S" "$c" prod_r3 > "$D/md11_$c.log" 2>&1 &
done
wait
tail -n 3 "$D/md11_SNP00157.log" "$D/md11_SNP00163.log"

#!/bin/bash
# run_fullprotocol_chunk.sh <chunk_number>
# Full-protocol redock (8 conformers, exhaustiveness 32, 10 poses -- matching
# the validated SNP-Lib/benchmark protocol exactly) of the top-500-by-
# fast-pass-score shortlist, split into 3 parallel chunks. Same 30-minute
# restart-loop pattern as run_fastpass_chunk.sh, for the same reason (a
# memory leak in g4_dock_library.py was observed to accumulate ~30MB/compound
# and cause OOM kills after several hundred compounds; see decision log,
# 2026-08-16/17). Poses are archived (--poses-dir) for downstream
# interaction-level analysis and evidence-chain reporting.
set -u
CHUNK="$1"
cd "/mnt/d/Agentic VS Workflow" || exit 1

LEDGER="20_new_runs/np_library_funnel_2026-08-15/full_ledger_chunk${CHUNK}.csv"
OUT="20_new_runs/np_library_funnel_2026-08-15/full_scores_chunk${CHUNK}.csv"
POSES="20_new_runs/np_library_funnel_2026-08-15/full_poses_chunk${CHUNK}"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate diffdock

LTOTAL=$(( $(wc -l < "$LEDGER") - 1 ))
echo "full-protocol chunk ${CHUNK}: ledger has ${LTOTAL} compounds"

while true; do
  if [ -f "$OUT" ]; then
    DONE=$(( $(wc -l < "$OUT") - 1 ))
  else
    DONE=0
  fi
  echo "=== full chunk ${CHUNK} restart cycle: ${DONE}/${LTOTAL} done ==="
  if [ "$DONE" -ge "$LTOTAL" ]; then
    echo "full chunk ${CHUNK} FULLY DONE"
    break
  fi
  timeout 1800 python 10_new_system/scripts/g4_dock_library.py \
    --ledger "$LEDGER" \
    --receptor "20_new_runs/snp_screen_2026-08-04/G2_preparation/4EIY_prepared.pdbqt" \
    --box "20_new_runs/snp_screen_2026-08-04/G2_preparation/docking_box.json" \
    --out "$OUT" \
    --poses-dir "$POSES" \
    --id-col compound_id --smiles-col SMILES \
    --n-conformers 8 --exhaustiveness 32 --n-poses 10 --cpu 6
done

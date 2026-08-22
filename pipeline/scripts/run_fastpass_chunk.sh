#!/bin/bash
# run_fastpass_chunk.sh <chunk_number>
# Restart-loop wrapper around g4_dock_library.py: relaunches every 30 minutes
# to bound the per-process memory footprint (a leak was observed causing OOM
# kills after several hundred compounds; see ADR log 2026-08-16). Relies
# entirely on g4_dock_library.py's own --out-based resume logic, so no work
# is lost across restarts.
set -u
CHUNK="$1"
cd "/mnt/d/Agentic VS Workflow" || exit 1

LEDGER="20_new_runs/np_library_funnel_2026-08-15/fastpass_ledger_chunk${CHUNK}.csv"
OUT="20_new_runs/np_library_funnel_2026-08-15/fastpass_scores_chunk${CHUNK}.csv"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate diffdock

LTOTAL=$(( $(wc -l < "$LEDGER") - 1 ))
echo "chunk ${CHUNK}: ledger has ${LTOTAL} compounds"

while true; do
  if [ -f "$OUT" ]; then
    DONE=$(( $(wc -l < "$OUT") - 1 ))
  else
    DONE=0
  fi
  echo "=== chunk ${CHUNK} restart cycle: ${DONE}/${LTOTAL} done ==="
  if [ "$DONE" -ge "$LTOTAL" ]; then
    echo "chunk ${CHUNK} FULLY DONE"
    break
  fi
  timeout 1800 python 10_new_system/scripts/g4_dock_library.py \
    --ledger "$LEDGER" \
    --receptor "20_new_runs/snp_screen_2026-08-04/G2_preparation/4EIY_prepared.pdbqt" \
    --box "20_new_runs/snp_screen_2026-08-04/G2_preparation/docking_box.json" \
    --out "$OUT" \
    --id-col compound_id --smiles-col SMILES \
    --n-conformers 2 --exhaustiveness 8 --n-poses 3
done

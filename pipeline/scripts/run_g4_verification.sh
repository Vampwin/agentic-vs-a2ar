#!/bin/bash
# run_g4_verification.sh <pass: A|B>
#
# Re-docks the G4 seed_repeatability subset (25 compounds, 5% of the 497
# successfully docked, selected under the contract's own seeded rule) using
# parameters identical to production, at two different thread counts:
#
#   pass A : --cpu 6   (the setting used for roughly the last 180 compounds)
#   pass B : --cpu 20  (all cores; the setting used for roughly the first 320)
#
# Two questions are answered by the same runs:
#   1. G4 `seed_repeatability` — does a rerun at a frozen seed reproduce the
#      original score to within 0.01 kcal/mol?
#   2. Concern 6.1 — does Vina's thread count change the result? The --cpu
#      flag was introduced partway through the production redock, so if
#      thread count matters the 500 reported scores are not homogeneous.
#
# Same 30-minute restart wrapper as production, for the same memory-leak
# reason, and the same --out-based resume logic.
set -u
PASS="$1"
cd "/mnt/d/Agentic VS Workflow" || exit 1

case "$PASS" in
  A) CPU=6  ;;
  B) CPU=20 ;;
  *) echo "usage: $0 A|B"; exit 2 ;;
esac

RUN="20_new_runs/g4_verification_2026-08-21"
LEDGER="$RUN/g4_subset_ledger.csv"
OUT="$RUN/rerun_pass${PASS}_cpu${CPU}.csv"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate diffdock

LTOTAL=$(( $(wc -l < "$LEDGER") - 1 ))
echo "G4 verification pass ${PASS} (cpu=${CPU}): ${LTOTAL} compounds"

while true; do
  if [ -f "$OUT" ]; then DONE=$(( $(wc -l < "$OUT") - 1 )); else DONE=0; fi
  echo "=== pass ${PASS} cycle: ${DONE}/${LTOTAL} ==="
  if [ "$DONE" -ge "$LTOTAL" ]; then echo "pass ${PASS} FULLY DONE"; break; fi
  timeout 1800 python 10_new_system/scripts/g4_dock_library.py \
    --ledger "$LEDGER" \
    --receptor "20_new_runs/snp_screen_2026-08-04/G2_preparation/4EIY_prepared.pdbqt" \
    --box "20_new_runs/snp_screen_2026-08-04/G2_preparation/docking_box.json" \
    --out "$OUT" \
    --id-col compound_id --smiles-col SMILES \
    --n-conformers 8 --exhaustiveness 32 --n-poses 10 --cpu "$CPU"
done

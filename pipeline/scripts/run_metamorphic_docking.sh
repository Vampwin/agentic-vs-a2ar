#!/bin/bash
# run_metamorphic_docking.sh <variant: M1|M2>
#
# Docks a metamorphic input variant with the command used for production and
# for the two seed_repeatability passes -- same receptor, same box, same
# conformer count, same exhaustiveness, same seed. Only the input ledger
# differs, and it differs in a way that cannot change the answer:
#
#   M1  the 25 rows are in a different order
#   M2  each SMILES is a different valid non-canonical string for the same
#       molecule (InChIKey identity enforced when the ledger was built)
#
# Thread count is 8 rather than the 6/20 used for passes A and B. That is
# deliberate: thread invariance was already established at 6 vs 20 with zero
# deviation, so running the metamorphic variants at a third setting costs
# nothing and widens the evidence. It also leaves headroom -- two variants at
# 8 threads plus the MD packing job stays under the 20 available cores, and
# CPU oversubscription has already cost this project one run.
#
# Same 30-minute restart wrapper as production, for the same memory-leak
# reason, and the same --out-based resume logic: killing this script and
# restarting it loses nothing.
set -u
VARIANT="$1"
cd "/mnt/d/Agentic VS Workflow" || exit 1

RUN="20_new_runs/metamorphic_2026-08-21"
case "$VARIANT" in
  M1) LEDGER="$RUN/ledger_M1_rowshuffle.csv"     ;;
  M2) LEDGER="$RUN/ledger_M2_smiles_rewrite.csv" ;;
  *)  echo "usage: $0 M1|M2"; exit 2 ;;
esac
OUT="$RUN/scores_${VARIANT}.csv"
CPU=8

source ~/miniconda3/etc/profile.d/conda.sh
conda activate diffdock

LTOTAL=$(( $(wc -l < "$LEDGER") - 1 ))
echo "metamorphic ${VARIANT} (cpu=${CPU}): ${LTOTAL} compounds"

while true; do
  if [ -f "$OUT" ]; then DONE=$(( $(wc -l < "$OUT") - 1 )); else DONE=0; fi
  echo "=== ${VARIANT} cycle: ${DONE}/${LTOTAL} ==="
  if [ "$DONE" -ge "$LTOTAL" ]; then echo "${VARIANT} FULLY DONE"; break; fi
  timeout 1800 python 10_new_system/scripts/g4_dock_library.py \
    --ledger "$LEDGER" \
    --receptor "20_new_runs/snp_screen_2026-08-04/G2_preparation/4EIY_prepared.pdbqt" \
    --box "20_new_runs/snp_screen_2026-08-04/G2_preparation/docking_box.json" \
    --out "$OUT" \
    --id-col compound_id --smiles-col SMILES \
    --n-conformers 8 --exhaustiveness 32 --n-poses 10 --cpu "$CPU"
done

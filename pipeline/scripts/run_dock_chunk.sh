#!/bin/bash
# run_dock_chunk.sh --ledger F --out F --poses DIR --conformers N --exhaustiveness N [--cpu N] [--seed N]
#
# Generalised replacement for run_fullprotocol_chunk.sh and run_fastpass_chunk.sh,
# which hard-coded their ledger, output and pose paths into
# 20_new_runs/np_library_funnel_2026-08-15/ -- a directory declared immutable by
# ADR-AIIM-092. R2, R3 and R4 all reuse this runner, so it is parameterised once
# here rather than copied three times (ADR-AIIM-133).
#
# The restart loop is not decorative and is kept verbatim: g4_dock_library.py was
# observed to accumulate roughly 30 MB per compound and to be OOM-killed after
# several hundred, so the loop re-enters and resumes from the rows already in
# --out. Resumption is by row count, which is why --out must be unique per run.
set -u

LEDGER=""; OUT=""; POSES=""; NCONF=8; EXH=32; NPOSES=10; CPU=6; SEED=20260729
# Per-restart-window timeout. The original hard-coded 1800 was sized for the
# 8-conformer production protocol at roughly 6 minutes per compound. It does
# not scale, and at 32 conformers (~26 min/compound) a single compound can
# exceed the window: the process is killed before writing its row, the loop
# restarts, retries the same compound, and thrashes forever with every health
# signal looking normal -- processes alive, CPU saturated, no errors. That is
# what happened on 2026-08-27 and cost 18 hours (ADR-AIIM-134). Default now
# scales with the conformer budget; override with --timeout.
TIMEOUT=""
# Consecutive restart cycles that may produce zero new rows before aborting.
# Without this the failure above is silent, which is the whole problem.
MAX_STALL_CYCLES=3
RECEPTOR="20_new_runs/snp_screen_2026-08-04/G2_preparation/4EIY_prepared.pdbqt"
BOX="20_new_runs/snp_screen_2026-08-04/G2_preparation/docking_box.json"

while [ $# -gt 0 ]; do
  case "$1" in
    --ledger) LEDGER="$2"; shift 2 ;;
    --out) OUT="$2"; shift 2 ;;
    --poses) POSES="$2"; shift 2 ;;
    --conformers) NCONF="$2"; shift 2 ;;
    --exhaustiveness) EXH="$2"; shift 2 ;;
    --n-poses) NPOSES="$2"; shift 2 ;;
    --cpu) CPU="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    --timeout) TIMEOUT="$2"; shift 2 ;;
    --max-stall-cycles) MAX_STALL_CYCLES="$2"; shift 2 ;;
    --receptor) RECEPTOR="$2"; shift 2 ;;
    --box) BOX="$2"; shift 2 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done

if [ -z "$LEDGER" ] || [ -z "$OUT" ]; then
  echo "usage: run_dock_chunk.sh --ledger F --out F [--poses DIR] [--conformers N] [--exhaustiveness N] [--cpu N] [--seed N]" >&2
  exit 2
fi

cd "/mnt/d/Agentic VS Workflow" || exit 1
[ -f "$LEDGER" ] || { echo "missing ledger: $LEDGER" >&2; exit 1; }
mkdir -p "$(dirname "$OUT")"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate diffdock

if [ -z "$TIMEOUT" ]; then
  TIMEOUT=$(( 900 * NCONF ))
  [ "$TIMEOUT" -lt 1800 ] && TIMEOUT=1800
fi

LTOTAL=$(( $(wc -l < "$LEDGER") - 1 ))
echo "ledger=$LEDGER  n=$LTOTAL  conformers=$NCONF  exhaustiveness=$EXH  cpu=$CPU  seed=$SEED"
echo "timeout=${TIMEOUT}s per restart cycle  max_stall_cycles=$MAX_STALL_CYCLES"
echo "out=$OUT"

POSE_ARGS=()
[ -n "$POSES" ] && POSE_ARGS=(--poses-dir "$POSES")

PREV=-1
STALL=0
while true; do
  if [ -f "$OUT" ]; then DONE=$(( $(wc -l < "$OUT") - 1 )); else DONE=0; fi
  echo "=== restart cycle: ${DONE}/${LTOTAL} done ==="
  if [ "$DONE" -ge "$LTOTAL" ]; then echo "FULLY DONE: $OUT"; break; fi
  # A cycle that adds nothing means the window is too short for the compound
  # it is stuck on. Retrying cannot fix that, so fail loudly instead of
  # burning CPU indefinitely while every other signal looks healthy.
  if [ "$DONE" -eq "$PREV" ]; then
    STALL=$(( STALL + 1 ))
    echo "!!! no progress in this cycle (${STALL}/${MAX_STALL_CYCLES})"
    if [ "$STALL" -ge "$MAX_STALL_CYCLES" ]; then
      echo "ABORT: ${MAX_STALL_CYCLES} consecutive cycles added no rows at ${DONE}/${LTOTAL}."
      echo "       The ${TIMEOUT}s window is too short for at least one compound"
      echo "       at ${NCONF} conformers. Raise --timeout and resume; the rows"
      echo "       already in $OUT are kept."
      exit 3
    fi
  else
    STALL=0
  fi
  PREV=$DONE
  timeout "$TIMEOUT" python 10_new_system/scripts/g4_dock_library.py \
    --ledger "$LEDGER" \
    --receptor "$RECEPTOR" \
    --box "$BOX" \
    --out "$OUT" \
    "${POSE_ARGS[@]}" \
    --id-col compound_id --smiles-col SMILES \
    --n-conformers "$NCONF" --exhaustiveness "$EXH" --n-poses "$NPOSES" \
    --cpu "$CPU" --seed "$SEED"
done

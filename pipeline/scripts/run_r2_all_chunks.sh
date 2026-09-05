#!/bin/bash
# run_r2_all_chunks.sh -- R2: re-dock all 500 top-ranked compounds under
# rewritten, InChIKey-verified SMILES, at the production protocol exactly
# (8 conformers, exhaustiveness 32, 10 poses, 6 threads per chunk).
#
# Threads per chunk default to 3, not the 6 the original run used, because the
# molecular dynamics on the other WSL distro needs 8 cores of the same 20 and
# GROMACS does not stop needing CPU just because it offloads to the GPU. Running
# three chunks at 6 alongside an mdrun took that run from 115 ns/day to about 21
# (2026-08-29). Override with R2_CPU.
#
# Three chunks in parallel at 6 threads each was the configuration the original
# full-protocol run used: 18 of 20 cores. Running them unbounded caused load
# 38-52 on 20 cores once already (decision log, CPU oversubscription incident),
# which is why --cpu is passed explicitly rather than left to default.
#
# Written as a file rather than an inline command because the chunk variable
# was being eaten by shell-escaping layers when passed through from Windows.
set -u
cd "/mnt/d/Agentic VS Workflow" || exit 1
RUN="20_new_runs/smiles_redock_500_2026-08-26"

for C in 1 2 3; do
  LEDGER="${RUN}/ledger_R2_chunk${C}.csv"
  if [ ! -f "$LEDGER" ]; then
    echo "missing ledger: $LEDGER" >&2
    exit 1
  fi
done

for C in 1 2 3; do
  nohup bash 10_new_system/scripts/run_dock_chunk.sh \
    --ledger "${RUN}/ledger_R2_chunk${C}.csv" \
    --out "${RUN}/scores_R2_chunk${C}.csv" \
    --conformers 8 --exhaustiveness 32 --n-poses 10 --cpu "${R2_CPU:-3}" \
    > "${RUN}/r2_chunk${C}.log" 2>&1 &
  echo "launched chunk ${C}, pid $!"
done
wait
echo "### R2 all chunks complete"

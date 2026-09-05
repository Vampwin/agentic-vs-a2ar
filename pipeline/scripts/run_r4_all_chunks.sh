#!/bin/bash
# run_r4_all_chunks.sh -- R4: the funnel-exclusion control at the FULL protocol.
#
# ADR-AIIM-150 answered M9's reweighting from the 2026-08-21 control, which was
# docked at the fast pass. That is the right comparison for asking whether the
# funnel's own cut was sound, because the fast pass is what made the cut. It is
# not the protocol the 500 selected compounds finally received, so it cannot say
# what fraction of the discarded compounds would have reached the top-500 cutoff
# had they been given the same treatment as the compounds that were kept.
#
# This re-docks the same 500 sampled compounds -- same ledger, same seeded
# sample, no re-draw -- at 8 conformers and exhaustiveness 32, the production
# protocol. The stratum weights are unchanged from ADR-AIIM-150: 13,882 rejected
# by the 2D screen and 3,000 by the 3D refinement, 250 sampled from each.
#
# Three chunks at 6 threads is the validated production configuration, 18 of 20
# cores. Nothing else is running now; if an mdrun starts, drop to 3 threads per
# chunk (R4_CPU) or GROMACS will be starved -- see ADR-AIIM-146.
set -u
cd "/mnt/d/Agentic VS Workflow" || exit 1
SRC="20_new_runs/funnel_exclusion_control_2026-08-21/exclusion_control_ledger.csv"
RUN="20_new_runs/funnel_recall_fullprotocol_2026-09-01"

[ -f "$SRC" ] || { echo "missing ledger: $SRC" >&2; exit 1; }
mkdir -p "$RUN"

# Split into three chunks, preserving the arm column. Deterministic: rows keep
# their ledger order, so a re-run produces the same chunks.
if [ ! -f "$RUN/ledger_R4_chunk1.csv" ]; then
  python3 - <<'PY'
import pandas as pd, pathlib
src = pathlib.Path("20_new_runs/funnel_exclusion_control_2026-08-21/exclusion_control_ledger.csv")
out = pathlib.Path("20_new_runs/funnel_recall_fullprotocol_2026-09-01")
d = pd.read_csv(src)
n = len(d)
bounds = [0, (n + 2) // 3, 2 * (n + 2) // 3, n]
for i in range(3):
    part = d.iloc[bounds[i]:bounds[i + 1]]
    part.to_csv(out / f"ledger_R4_chunk{i+1}.csv", index=False)
    print(f"chunk {i+1}: {len(part)} rows  ({part.arm.value_counts().to_dict()})")
PY
fi

for C in 1 2 3; do
  [ -f "$RUN/ledger_R4_chunk${C}.csv" ] || { echo "missing chunk $C" >&2; exit 1; }
done

for C in 1 2 3; do
  nohup bash 10_new_system/scripts/run_dock_chunk.sh \
    --ledger "${RUN}/ledger_R4_chunk${C}.csv" \
    --out "${RUN}/scores_R4_chunk${C}.csv" \
    --conformers 8 --exhaustiveness 32 --n-poses 10 --cpu "${R4_CPU:-6}" \
    > "${RUN}/r4_chunk${C}.log" 2>&1 &
  echo "launched chunk ${C}, pid $!"
done
wait
echo "### R4 all chunks complete"

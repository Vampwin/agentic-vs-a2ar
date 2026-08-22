#!/bin/bash
# ensure_fullprotocol_running.sh
#
# Idempotent supervisor for the 3 full-protocol redock chunks. Safe to run
# repeatedly (e.g. every status check, or on a schedule): for each chunk it
# checks whether a g4_dock_library.py process for that specific chunk's
# ledger is already alive; only relaunches chunks that are actually down.
# This replaces the previous manual routine (check ps, relaunch 1/2/3 in
# separate sequential commands with visual verification) with one atomic,
# self-verifying, logged operation -- written after two full-WSL-VM-reset
# events in one session both silently killed all 3 wrapper processes.
#
# Also runs a duplicate-compound_id check on each score CSV after (re)launch,
# as a integrity guard against any double-write race (e.g. from an orphaned
# process surviving a kill, as previously happened during a smoke test).
#
# Usage: bash ensure_fullprotocol_running.sh
set -u
cd "/mnt/d/Agentic VS Workflow" || exit 1

RUNDIR="20_new_runs/np_library_funnel_2026-08-15"
SUPLOG="${RUNDIR}/supervisor_events.log"
ts() { date '+%Y-%m-%dT%H:%M:%S%z'; }

for n in 1 2 3; do
  ledger="${RUNDIR}/full_ledger_chunk${n}.csv"
  if pgrep -f "g4_dock_library.py --ledger ${ledger}" >/dev/null 2>&1; then
    echo "$(ts) chunk${n}: already running, skipped" >> "$SUPLOG"
    continue
  fi
  echo "$(ts) chunk${n}: NOT running -- relaunching" >> "$SUPLOG"
  nohup bash "10_new_system/scripts/run_fullprotocol_chunk.sh" "$n" \
    >> "${RUNDIR}/full_chunk${n}_wrapper.log" 2>&1 &
  disown
  sleep 5
  if pgrep -f "run_fullprotocol_chunk.sh $n" >/dev/null 2>&1; then
    echo "$(ts) chunk${n}: relaunch confirmed (wrapper alive)" >> "$SUPLOG"
  else
    echo "$(ts) chunk${n}: WARNING relaunch NOT confirmed after 5s" >> "$SUPLOG"
  fi
done

for n in 1 2 3; do
  f="${RUNDIR}/full_scores_chunk${n}.csv"
  [ -f "$f" ] || continue
  dupes=$(tail -n +2 "$f" | cut -d',' -f1 | sort | uniq -d | wc -l)
  if [ "$dupes" -gt 0 ]; then
    echo "$(ts) chunk${n}: WARNING ${dupes} duplicate compound_id rows in ${f}" >> "$SUPLOG"
  fi
done

echo "$(ts) supervisor check complete" >> "$SUPLOG"
echo "--- summary ---"
tail -n 12 "$SUPLOG"

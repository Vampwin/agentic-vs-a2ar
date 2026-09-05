#!/bin/bash
# watch_r4.sh -- notify when R4's three chunks have docked all 500 compounds.
#
# A file on disk, not an inline string: four background jobs in this project had
# every shell variable consumed by the Windows -> wsl -> bash escaping layers,
# and three of them ran for hours looking healthy while being incapable of
# reporting anything (ADR-AIIM-148).
#
# Checks liveness and stall as well as progress. A docking loop that restarts
# forever without writing a row looks identical to a slow one from outside,
# which is what cost eighteen hours on 2026-08-27 (ADR-AIIM-134).
set -u
cd "/mnt/d/Agentic VS Workflow" || exit 1
RUN="20_new_runs/funnel_recall_fullprotocol_2026-09-01"
EXPECTED=(167 167 166)
STALL_CHECKS=8          # 8 x 15 min = 2 h with no new row anywhere

rows() {
  local f="$RUN/scores_R4_chunk$1.csv"
  if [ -f "$f" ]; then echo $(( $(wc -l < "$f") - 1 )); else echo 0; fi
}

last_total=-1
stalled=0

while true; do
  total=0; done_n=0; line=""
  for c in 1 2 3; do
    n=$(rows "$c"); want=${EXPECTED[$((c - 1))]}
    total=$((total + n)); line="$line c$c=$n/$want"
    [ "$n" -ge "$want" ] && done_n=$((done_n + 1))
  done

  if [ "$done_n" -ge 3 ]; then
    echo "R4 COMPLETE: $total of 500 compounds -$line"; exit 0
  fi
  if ! pgrep -f g4_dock_library.py > /dev/null \
     && ! pgrep -f run_dock_chunk.sh > /dev/null; then
    echo "R4 STOPPED at $total of 500 -$line (no docking process alive)"; exit 1
  fi
  if [ "$total" -eq "$last_total" ]; then
    stalled=$((stalled + 1))
    if [ "$stalled" -ge "$STALL_CHECKS" ]; then
      echo "R4 STALLED at $total of 500 -$line (no new row in 2 h, processes alive)"
      exit 2
    fi
  else
    stalled=0; last_total=$total
  fi
  sleep 900
done

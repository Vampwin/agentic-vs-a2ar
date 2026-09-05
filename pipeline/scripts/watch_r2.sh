#!/bin/bash
# watch_r2.sh -- notify when R2's three chunks have docked all 500 compounds.
#
# The fourth waiter in this project written as a file rather than an inline
# string, after the first three were passed through Windows -> wsl -> bash and
# had every shell variable consumed by the escaping layers. The R2 waiter this
# one replaces spent a day printing
#
#     /bin/bash: line 1: /scores_R2_chunk1.csv: No such file or directory
#     bash: line 9: [: -ge: unary operator expected
#
# on a loop, and would never have reported completion. See ADR-AIIM-148.
#
# Progress alone is not the check. A docking loop that restarts forever without
# writing a row looks identical to a slow one from the outside, which is what
# cost eighteen hours on 2026-08-27 (ADR-AIIM-134). Liveness is checked too, and
# a stall is reported rather than waited through.
set -u
cd "/mnt/d/Agentic VS Workflow" || exit 1
RUN="20_new_runs/smiles_redock_500_2026-08-26"
EXPECTED=(167 167 166)
STALL_CHECKS=8          # 8 x 15 min = 2 h of no new rows anywhere

rows() {
  local f="$RUN/scores_R2_chunk$1.csv"
  if [ -f "$f" ]; then
    echo $(( $(wc -l < "$f") - 1 ))
  else
    echo 0
  fi
}

last_total=-1
stalled=0

while true; do
  total=0
  done_n=0
  line=""
  for c in 1 2 3; do
    n=$(rows "$c")
    want=${EXPECTED[$((c - 1))]}
    total=$((total + n))
    line="$line c$c=$n/$want"
    [ "$n" -ge "$want" ] && done_n=$((done_n + 1))
  done

  if [ "$done_n" -ge 3 ]; then
    echo "R2 COMPLETE: $total of 500 compounds -$line"
    exit 0
  fi

  if ! pgrep -f g4_dock_library.py > /dev/null \
     && ! pgrep -f run_dock_chunk.sh > /dev/null; then
    echo "R2 STOPPED at $total of 500 -$line (no docking process alive)"
    exit 1
  fi

  if [ "$total" -eq "$last_total" ]; then
    stalled=$((stalled + 1))
    if [ "$stalled" -ge "$STALL_CHECKS" ]; then
      echo "R2 STALLED at $total of 500 -$line (no new row in 2 h, processes alive)"
      exit 2
    fi
  else
    stalled=0
    last_total=$total
  fi

  sleep 900
done

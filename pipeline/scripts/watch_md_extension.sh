#!/bin/bash
# watch_md_extension.sh -- notify when NPLIB01185's three replicates reach 100 ns.
#
# Watches the extension queued on 2026-08-24 by
# _launch/wait_then_extend_NPLIB01185.sh, which takes each of the three 50 ns
# replicates to 100 ns with `convert-tpr -extend` and a checkpoint continuation.
#
# Written as a file because the previous waiter was passed inline through
# Windows -> wsl -> bash and every shell variable in it was consumed by the
# escaping layers: its body read `[ -f "/prod_r.gro" ]` and `[ "" -ge 3 ]`, so
# it could never have reported completion. It looked healthy in `ps` the whole
# time, which is the same failure this project has now hit three times.
set -u
cd "/mnt/d/Agentic VS Workflow" || exit 1
RUN="20_new_runs/md_candidates_2026-08-22/NPLIB01185/04_run"
TARGET_NS=99.9          # a shade under 100 to tolerate the last frame's rounding

ns_of() {
  local f="$RUN/prod_r$1.log"
  [ -f "$f" ] || { echo 0; return; }
  grep -A1 "^ *Step *Time" "$f" | grep -E "^ +[0-9]+" | tail -1 \
    | awk '{printf "%.1f", $2/1000}'
}

while true; do
  done_n=0
  line=""
  for r in 1 2 3; do
    ns=$(ns_of "$r")
    line="$line r$r=${ns}ns"
    if awk -v a="$ns" -v b="$TARGET_NS" 'BEGIN{exit !(a>=b)}'; then
      done_n=$((done_n + 1))
    fi
  done

  if [ "$done_n" -ge 3 ]; then
    echo "MD EXTENSION COMPLETE: NPLIB01185 3/3 replicates at 100 ns -$line"
    exit 0
  fi

  # Liveness, not just progress: the driver may have exited between replicates.
  if ! pgrep -f "gmx mdrun" > /dev/null \
     && ! pgrep -f "md_10b_extend_candidate" > /dev/null \
     && ! pgrep -f "wait_then_extend_NPLIB01185" > /dev/null; then
    echo "MD EXTENSION STOPPED at $done_n/3 -$line (no mdrun, no driver alive)"
    exit 1
  fi
  sleep 900
done

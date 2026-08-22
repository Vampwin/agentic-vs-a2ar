#!/bin/bash
# md_04_benchmark.sh — energy-minimise the built system and measure GROMACS
# throughput (ns/day) on this hardware.
#
# Purpose is timing, not science: the number produced here sets the replicate
# length that the production campaign can afford, which the protocol design
# deliberately left open ([BENCHMARK FIRST] in
# 00_project_control/25_MD_PROTOCOL_DESIGN.md §6) rather than guessing.
set -u

G=/opt/gromacs/2026.3-cuda13.1-sm120/bin/gmx
B="/mnt/d/Agentic VS Workflow/20_new_runs/md_benchmark_2026-08-21"
cd "$B/04_run" || exit 1

TOP="$B/03_build/system_gmx.top"
GRO="$B/03_build/system_gmx.gro"

echo "=== index groups (SOLU / MEMB / SOLV for 3-group coupling) ==="
# Amber-converted systems have no default membrane group, so build one.
printf '%s\n' \
  "1 | 13" \
  "name 0 SYSTEM" \
  "q" | "$G" make_ndx -f "$GRO" -o index.ndx > make_ndx.log 2>&1 || true

echo "=== energy minimisation ==="
"$G" grompp -f em.mdp -c "$GRO" -p "$TOP" -o em.tpr -maxwarn 5 > grompp_em.log 2>&1
if [ ! -f em.tpr ]; then echo "grompp(em) FAILED"; tail -25 grompp_em.log; exit 1; fi
"$G" mdrun -deffnm em -ntmpi 1 -ntomp 8 -nb gpu -pin on > mdrun_em.log 2>&1
tail -4 em.log | sed 's/^/   /'

echo "=== benchmark run ==="
"$G" grompp -f bench.mdp -c em.gro -p "$TOP" -o bench.tpr -maxwarn 10 > grompp_bench.log 2>&1
if [ ! -f bench.tpr ]; then echo "grompp(bench) FAILED"; tail -30 grompp_bench.log; exit 1; fi
"$G" mdrun -deffnm bench -ntmpi 1 -ntomp 8 -nb gpu -bonded gpu -pme gpu -pin on > mdrun_bench.log 2>&1

echo "=== PERFORMANCE ==="
grep -A3 "Performance:" bench.log | sed 's/^/   /'
grep -E "^\s+Time:|ns/day|hour/ns" bench.log | tail -4 | sed 's/^/   /'

#!/bin/bash
# ADR-AIIM-314 — three velocity-seeded production replicates of the ZM241385
# control, run one at a time so the result is usable after each.
#
# Pre-registration SHA-256
#   214e646bf96ffe14e001e3957e21d1f18e67f87935642e16383daa4c095a9a18
# Decision log at freeze (text/LF)
#   6bba84a003aeb6c8ea2406823d175c09fd9cd76c2d212bb887e3eed14d575a75
#
# The existing control run in md_production_2026-08-21 is READ ONLY: its
# npt3.gro, coupling.ndx, topology and position-restraint files were copied
# here with their SHA-256 recorded in inputs_sha256.json, and nothing is
# written back to it -- not even grompp's mdout.mdp, which -po redirects here.
#
# Waits for ADR-AIIM-305's DiffDock run to release the GPU before any mdrun.
# Resumable: a replicate whose .gro exists is skipped, and one interrupted
# mid-flight continues from its checkpoint with -cpi. This project has lost
# long runs to four host restarts.
#
#   usage: run_control_replicates.sh [--grompp-only] [--no-wait]
set -u

G=/opt/gromacs/2026.3-cuda13.1-sm120/bin/gmx
D="/mnt/d/Agentic VS Workflow/20_new_runs/md_control_replicates_2026-10-06"
SENTINEL="/mnt/d/Agentic VS Workflow/20_new_runs/diffdock_dude_2026-10-05/_stage1_done"
REPS="prod_r2 prod_r3"   # r4 deferred: SNP00157 takes the GPU after r3 (ADR-AIIM-316 constraint 4)
GROMPP_ONLY=0
WAIT=1

for a in "$@"; do
  case "$a" in
    --grompp-only) GROMPP_ONLY=1 ;;
    --no-wait)     WAIT=0 ;;
  esac
done

cd "$D" || exit 1
log() { echo "[$(date '+%F %T')] $*"; }

# ---- grompp: CPU only, safe to do while the GPU is busy -------------------
for r in $REPS; do
  if [ -s "$r.tpr" ]; then
    log "$r.tpr exists, skipping grompp"
    continue
  fi
  log "grompp $r"
  # -maxwarn 1 admits exactly one known, pre-existing warning: the published
  # control system carries a net charge under Ewald (Na+ 31, Cl- 41, plus the
  # conserved allosteric Na+). That is a property of the system this run must
  # replicate, so it is NOT neutralised here -- adding counter-ions would make
  # these runs something other than replicates of the published trajectory.
  # The original run passed -maxwarn 5 and -maxwarn 10; 1 is used instead so
  # any NEW warning still stops the build rather than being masked.
  "$G" grompp -f "$r.mdp" -c npt3.gro -r npt3.gro -p system_gmx.top \
       -n coupling.ndx -o "$r.tpr" -po "mdout_$r.mdp" -maxwarn 1 \
       > "grompp_$r.log" 2>&1
  rc=$?
  log "grompp $r exit=$rc"
  if [ $rc -ne 0 ] || [ ! -s "$r.tpr" ]; then
    log "grompp $r FAILED -- last lines:"
    tail -20 "grompp_$r.log"
    exit 1
  fi
done

if [ $GROMPP_ONLY -eq 1 ]; then
  log "grompp-only: three tpr files ready, not starting mdrun"
  ls -la ./*.tpr
  exit 0
fi

# ---- wait for the GPU ------------------------------------------------------
if [ $WAIT -eq 1 ]; then
  log "waiting for ADR-AIIM-305 to release the GPU (sentinel: _stage1_done)"
  while [ ! -f "$SENTINEL" ]; do
    sleep 300
  done
  log "sentinel seen; DiffDock stage 1 finished"
  sleep 60
fi

# ---- production, one replicate at a time ----------------------------------
for r in $REPS; do
  if [ -s "$r.gro" ]; then
    log "$r already complete, skipping"
    continue
  fi
  # WATCHDOG. This host's WSL GPU passthrough faults: `dmesg` shows
  # `dxgk: dxgkio_query_adapter_info: Ioctl failed` and the GPU process dies
  # silently, writing no error. It killed ADR-AIIM-305's DiffDock at 335/593
  # and this run's replicate 2 at 49.5 ns, each time unnoticed for hours.
  # GROMACS checkpoints every 15 minutes (-cpt 15), so a retry loop bounds the
  # loss at one checkpoint interval instead of however long until a human looks.
  ATTEMPT=0
  while [ ! -s "$r.gro" ] && [ $ATTEMPT -lt 40 ]; do
    ATTEMPT=$((ATTEMPT + 1))
    [ $ATTEMPT -gt 1 ] && log "$r attempt $ATTEMPT (previous one died; resuming)"
  CPI=""
  [ -f "$r.cpt" ] && CPI="-cpi $r.cpt"
  log "mdrun $r starting ${CPI:+(resuming from checkpoint)}"
  # Flags are the published runs' verbatim -- prod.log for the control and
  # prod_r1.log for the candidates both record exactly this line. In
  # particular `-update gpu` is NOT used: the water model carries virtual
  # sites, and GPU update refuses them. Matching the published invocation is
  # the point, since these are meant to be replicates of that trajectory.
  # shellcheck disable=SC2086
  "$G" mdrun -deffnm "$r" $CPI -ntmpi 1 -ntomp 8 \
       -nb gpu -pme gpu -bonded gpu -cpt 15 >> "$r.mdrun.log" 2>&1
  rc=$?
  log "mdrun $r exit=$rc"
    [ -s "$r.gro" ] || sleep 60      # let the GPU settle before retrying
  done

  if [ -s "$r.gro" ]; then
    grep -E "Performance:" "$r.log" | tail -1
    log "$r COMPLETE after $ATTEMPT attempt(s)"
  else
    log "$r gave up after $ATTEMPT attempts -- the GPU fault is not transient"
    exit 1
  fi
done

log "all replicates complete"
echo "DONE" > "$D/_replicates_done"

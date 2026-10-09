#!/bin/bash
# ADR-AIIM-305 stage 1. Parameters are ADR-AIIM-111's verbatim; changing
# samples_per_complex to halve the cost would change the method and is
# deliberately not done (ADR-AIIM-305 constraint 3).
#
# DiffDock is NOT installed, updated or reconfigured here. The existing conda
# env and checked-out commit are used as they stand and the commit is recorded.
#
# Detached with nohup because the job is ~10.7 h and must outlive the caller.
set -u
SRC="/mnt/d/Agentic VS Workflow/20_new_runs/diffdock_dude_2026-10-05"
W="$HOME/diffdock_dude_stage1"

source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || \
  source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate diffdock

mkdir -p "$W"
cp "$SRC/4EIY_protein_only.pdb" "$W/"
# DiffDock resolves protein_path against its own CWD, not the CSV's location,
# so the path is made absolute here rather than left relative (as ADR-AIIM-111).
sed "s|4EIY_protein_only.pdb|$W/4EIY_protein_only.pdb|" \
    "$SRC/diffdock_input.csv" > "$W/diffdock_input.csv"

cd ~/diffdock || exit 1
COMMIT=$(git rev-parse --short HEAD)
echo "$COMMIT" > "$SRC/diffdock_commit.txt"

{
  echo "[$(date '+%F %T')] ADR-AIIM-305 stage 1 start  commit $COMMIT"
  echo "[$(date '+%F %T')] 593 complexes, samples_per_complex 20"
  nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
  python -m inference \
    --config default_inference_args.yaml \
    --protein_ligand_csv "$W/diffdock_input.csv" \
    --out_dir "$W/results" \
    --samples_per_complex 20
  echo "[$(date '+%F %T')] DiffDock exit=$?"
  mkdir -p "$SRC/results"
  cp -r "$W/results/." "$SRC/results/" 2>/dev/null
  echo "[$(date '+%F %T')] results copied back"
  echo "DONE" > "$SRC/_stage1_done"
} > "$SRC/diffdock.log" 2>&1 &

echo "launched pid $! -- log: $SRC/diffdock.log"

#!/bin/bash
# ADR-AIIM-111 pre-registered orthogonal-engine check.
# Parameters fixed by that ADR, identical for control and candidates.
set -u
SRC="/mnt/d/Agentic VS Workflow/20_new_runs/diffdock_orthogonal_2026-08-23"
W="$HOME/diffdock_a2ar_run"
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || \
  source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate diffdock

rm -rf "$W"; mkdir -p "$W"
cp "$SRC/4EIY_protein_only.pdb" "$W/"
# DiffDock resolves protein_path against its own CWD, not the CSV's location,
# so the path is made absolute here rather than left relative.
sed "s|4EIY_protein_only.pdb|$W/4EIY_protein_only.pdb|" \
    "$SRC/diffdock_input.csv" > "$W/diffdock_input.csv"
head -2 "$W/diffdock_input.csv"
cd ~/diffdock || exit 1

echo "[$(date '+%F %T')] DiffDock start  commit $(git rev-parse --short HEAD)"
python -m inference \
  --config default_inference_args.yaml \
  --protein_ligand_csv "$W/diffdock_input.csv" \
  --out_dir "$W/results" \
  --samples_per_complex 20
echo "[$(date '+%F %T')] DiffDock exit=$?"

mkdir -p "$SRC/results"
cp -r "$W/results/." "$SRC/results/" 2>/dev/null
echo "[$(date '+%F %T')] results copied back"

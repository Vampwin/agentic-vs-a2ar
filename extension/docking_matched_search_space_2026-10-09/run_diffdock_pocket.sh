#!/bin/bash
# ADR-AIIM-353 arm B -- DiffDock-L on the protein cropped to 12 A around the crystal ZM241385 pose.
# Same command and parameters as ADR-AIIM-317 (default_inference_args.yaml, --samples_per_complex 20); only the
# protein file differs. The installed DiffDock is not installed, updated or reconfigured; its commit is recorded.
set -u
ROOT="/mnt/d/Agentic VS Workflow"
RUN="$ROOT/20_new_runs/docking_matched_search_space_2026-10-09"
W="$HOME/diffdock_pocket"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate diffdock
mkdir -p "$W"
cp "$RUN/inputs/4EIY_pocket12A.pdb" "$W/"
sed "s|4EIY_pocket12A.pdb|$W/4EIY_pocket12A.pdb|" "$RUN/inputs/diffdock_pocket_input.csv" > "$W/diffdock_input.csv"
cd ~/diffdock || exit 1
git rev-parse --short HEAD > "$RUN/diffdock_commit.txt"
echo "[$(date '+%F %T')] DiffDock commit $(cat "$RUN/diffdock_commit.txt"), 16 complexes, samples_per_complex 20"
nvidia-smi --query-gpu=name,memory.used,memory.total --format=csv,noheader
python -m inference \
  --config default_inference_args.yaml \
  --protein_ligand_csv "$W/diffdock_input.csv" \
  --out_dir "$W/results" \
  --samples_per_complex 20
echo "[$(date '+%F %T')] DiffDock exit=$?"
mkdir -p "$RUN/diffdock_pocket/results"
cp -r "$W/results/." "$RUN/diffdock_pocket/results/" 2>/dev/null
echo "[$(date '+%F %T')] results copied back: $(ls "$RUN/diffdock_pocket/results" | wc -l) entries"
echo DONE > "$RUN/diffdock_pocket/_done"

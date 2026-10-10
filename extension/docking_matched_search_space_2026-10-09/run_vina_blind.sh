#!/bin/bash
# ADR-AIIM-353 arm A -- Vina blind. The frozen g4_dock_library.py, unchanged, with the whole-receptor box
# of inputs/vina_blind_box.json and exhaustiveness 64; everything else as ADR-AIIM-313.
# Pre-registration SHA-256: see PREREG_SHA256.txt. Nothing is installed or reconfigured.
set -u
ROOT="/mnt/d/Agentic VS Workflow"
RUN="$ROOT/20_new_runs/docking_matched_search_space_2026-10-09"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate diffdock
mkdir -p "$RUN/vina_blind/poses"
echo "[$(date '+%F %T')] start Vina blind, 16 compounds, control first"
python "$ROOT/10_new_system/scripts/g4_dock_library.py" \
  --ledger "$RUN/inputs/vina_blind_ledger.csv" \
  --receptor "$ROOT/20_new_runs/snp_screen_2026-08-04/G2_preparation/4EIY_prepared.pdbqt" \
  --box "$RUN/inputs/vina_blind_box.json" \
  --out "$RUN/vina_blind/vina_blind_scores.csv" \
  --poses-dir "$RUN/vina_blind/poses" \
  --exhaustiveness 64 \
  --cpu 12
echo "[$(date '+%F %T')] Vina blind exit=$?"
echo DONE > "$RUN/vina_blind/_done"

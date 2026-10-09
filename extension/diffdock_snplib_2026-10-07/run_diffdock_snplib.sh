#!/bin/bash
# ADR-AIIM-317 -- DiffDock on the fifteen SNP-Lib survivors, as a pose-level
# check on the Vina poses. NOT a ranking channel: ADR-AIIM-305 stage 2 stays
# unrun and DiffDock confidence orders nothing.
#
# Pre-registration SHA-256
#   4386e354d4ca7bd1e8a08a36fbfb50ab795f24741853e3e25b66f1957469a265
# Decision log at freeze (text/LF)
#   3afdea45b0b38d1f4d498025d2a76a01546e86ad2297d34f3dfcaca61b3a3ff6
#
# Parameters are ADR-AIIM-111's verbatim, as ADR-AIIM-305 used them. DiffDock
# is NOT installed, updated or reconfigured; the existing env and commit are
# used as they stand and the commit is recorded.
#
# WAITS for the control's third replicate so no molecular-dynamics time is
# lost: DiffDock needs ~7.6 GB and GROMACS ~1.8 GB of an 8 GB card, so they
# cannot share it. ~15 complexes at the 131 s/complex measured on 2026-10-06
# is roughly 33 minutes.
set -u
SRC="/mnt/d/Agentic VS Workflow/20_new_runs/diffdock_snplib_2026-10-07"
GATE="/mnt/d/Agentic VS Workflow/20_new_runs/md_control_replicates_2026-10-06/prod_r3.gro"
W="$HOME/diffdock_snplib"

echo "[$(date '+%F %T')] waiting for control replicate 3 to finish (gate: prod_r3.gro)"
while [ ! -s "$GATE" ]; do
  sleep 300
done
echo "[$(date '+%F %T')] gate cleared"
sleep 90   # let GROMACS release the card before claiming it

source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || \
  source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate diffdock

mkdir -p "$W"
cp "$SRC/4EIY_protein_only.pdb" "$W/"
# DiffDock resolves protein_path against its own CWD, so the path is made
# absolute here rather than left relative (as ADR-AIIM-305 did).
sed "s|4EIY_protein_only.pdb|$W/4EIY_protein_only.pdb|" \
    "$SRC/diffdock_input.csv" > "$W/diffdock_input.csv"

cd ~/diffdock || exit 1
COMMIT=$(git rev-parse --short HEAD)
echo "$COMMIT" > "$SRC/diffdock_commit.txt"
echo "[$(date '+%F %T')] DiffDock commit $COMMIT, 15 complexes, samples_per_complex 20"
nvidia-smi --query-gpu=name,memory.used,memory.total --format=csv,noheader

python -m inference \
  --config default_inference_args.yaml \
  --protein_ligand_csv "$W/diffdock_input.csv" \
  --out_dir "$W/results" \
  --samples_per_complex 20
echo "[$(date '+%F %T')] DiffDock exit=$?"

# Copy back immediately. The ADR-AIIM-305 run died before its copy step and its
# results survived only because they were recovered from the distro by hand.
mkdir -p "$SRC/results"
cp -r "$W/results/." "$SRC/results/" 2>/dev/null
echo "[$(date '+%F %T')] results copied back: $(ls "$SRC/results" | wc -l) entries"
echo "DONE" > "$SRC/_diffdock_done"

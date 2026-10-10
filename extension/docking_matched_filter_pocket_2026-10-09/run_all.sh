#!/bin/bash
# ADR-AIIM-355: arm C (DiffDock-Pocket: flexible, rigid, C-alpha centre) then arm D (DiffDock-L blind, 100 samples).
# Sequential on the one GPU; control first within every CSV. Results are copied back after every run.
set -u
ROOT="/mnt/d/Agentic VS Workflow"
RUN="$ROOT/20_new_runs/docking_matched_filter_pocket_2026-10-09"
PROT="$ROOT/20_new_runs/diffdock_snplib_2026-10-07/4EIY_protein_only.pdb"
source ~/miniconda3/etc/profile.d/conda.sh

conda activate diffdock_pocket
cd ~/diffdock_pocket_repo || exit 1
git rev-parse HEAD > "$RUN/diffdock_pocket_commit.txt"
W="$HOME/ddp_run"; mkdir -p "$W"
for V in C Cca; do
  sed "s|4EIY_protein_only.pdb|$PROT|" "$RUN/inputs/pocket_$V.csv" > "$W/pocket_$V.csv"
done
run_pocket () {   # name csv extra-args
  echo "[$(date '+%F %T')] DiffDock-Pocket $1"
  python inference.py --protein_ligand_csv "$2" --samples_per_complex 40 --batch_size 20 --out_dir "$W/$1" ${3:-}
  echo "[$(date '+%F %T')] $1 exit=$?"
  mkdir -p "$RUN/$1" && cp -r "$W/$1/." "$RUN/$1/" 2>/dev/null
}
run_pocket pocket_C "$W/pocket_C.csv" ""
run_pocket pocket_Crigid "$W/pocket_C.csv" "--rigid"
run_pocket pocket_Cca "$W/pocket_Cca.csv" ""
echo DONE > "$RUN/_done_pocket"

conda activate diffdock
mkdir -p "$HOME/ddl_run"; cp "$PROT" "$HOME/ddl_run/"
sed "s|4EIY_protein_only.pdb|$HOME/ddl_run/4EIY_protein_only.pdb|" "$RUN/inputs/diffdock_l_input.csv" > "$HOME/ddl_run/diffdock_l_input.csv"
cp "$RUN/inputs/config_100samples.yaml" "$HOME/diffdock/config_100samples_adr355.yaml"
cd ~/diffdock || exit 1
git rev-parse --short HEAD > "$RUN/diffdock_l_commit.txt"
echo "[$(date '+%F %T')] DiffDock-L blind 100 samples, commit $(cat "$RUN/diffdock_l_commit.txt")"
python -m inference --config config_100samples_adr355.yaml \
  --protein_ligand_csv "$HOME/ddl_run/diffdock_l_input.csv" --out_dir "$HOME/ddl_run/results" --samples_per_complex 100
echo "[$(date '+%F %T')] DiffDock-L exit=$?"
mkdir -p "$RUN/diffdock_l_100/results" && cp -r "$HOME/ddl_run/results/." "$RUN/diffdock_l_100/results/" 2>/dev/null
echo DONE > "$RUN/_done_ddl"

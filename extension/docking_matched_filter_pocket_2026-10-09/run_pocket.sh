#!/bin/bash
# ADR-AIIM-355 arm C rerun: DiffDock-Pocket only, after the first attempt skipped every row (SMILES in the CSV ligand
# column is not accepted by its loader; ligands are now SDF files made from the same SMILES). Parameters unchanged.
set -u
ROOT="/mnt/d/Agentic VS Workflow"
RUN="$ROOT/20_new_runs/docking_matched_filter_pocket_2026-10-09"
PROT="$ROOT/20_new_runs/diffdock_snplib_2026-10-07/4EIY_protein_only.pdb"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate diffdock_pocket
cd ~/diffdock_pocket_repo || exit 1
W="$HOME/ddp_run"; rm -rf "$W"; mkdir -p "$W"
for V in C Cca; do sed "s|4EIY_protein_only.pdb|$PROT|" "$RUN/inputs/pocket_$V.csv" > "$W/pocket_$V.csv"; done
run_pocket () {
  echo "[$(date '+%F %T')] DiffDock-Pocket $1"
  python inference.py --protein_ligand_csv "$2" --samples_per_complex 40 --batch_size 2 --out_dir "$W/$1" ${3:-}
  echo "[$(date '+%F %T')] $1 exit=$?"
  rm -rf "$RUN/$1"; mkdir -p "$RUN/$1" && cp -r "$W/$1/." "$RUN/$1/" 2>/dev/null
}
run_pocket pocket_C "$W/pocket_C.csv" ""
run_pocket pocket_Crigid "$W/pocket_C.csv" "--rigid"
run_pocket pocket_Cca "$W/pocket_Cca.csv" ""
echo DONE > "$RUN/_done_pocket2"

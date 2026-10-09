#!/bin/bash
export AMBERHOME=/home/a2amd/A2A_MD/ambertools26_env
export PATH=$AMBERHOME/bin:$PATH
cd "/mnt/d/Agentic VS Workflow/20_new_runs/md_candidates_2026-08-22/SNP00163/03_build" || exit 1
antechamber -i ligand_posed.sdf -fi sdf -o LIG.mol2 -fo mol2 -c bcc -nc 0 -at gaff2 -rn LIG -s 0 > antechamber.log 2>&1
echo "antechamber exit=$?"
parmchk2 -i LIG.mol2 -f mol2 -o LIG.frcmod -s gaff2 > parmchk2.log 2>&1
echo "parmchk2 exit=$?"
grep -c ATTN LIG.frcmod || true

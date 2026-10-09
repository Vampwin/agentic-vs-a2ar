#!/bin/bash
# ADR-AIIM-316 -- ligand parameterisation for SNP00157.
#
# This is the step md_07_build_candidate.py would have emitted. That script is
# hardcoded to the NP-library funnel and its prioritisation ledger, neither of
# which contains an SNP-Lib compound, and 10_new_system is not modified. The
# settings below are its _cmd.sh verbatim: gaff2, AM1-BCC, -nc 0.
#
# -nc 0 is correct and checked, not assumed: the rank-1 pose parses to formal
# charge 0 with 23 heavy atoms, and its InChIKey NPNNKDMSXVRADT-WEVVVXLNSA-N
# equals the InChIKey of the SNP ledger's canonical SMILES for SNP00157.
set -u
export AMBERHOME=/home/a2amd/A2A_MD/ambertools26_env
export PATH=$AMBERHOME/bin:$PATH

B="/mnt/d/Agentic VS Workflow/20_new_runs/md_candidates_2026-08-22/SNP00157/03_build"
cd "$B" || exit 1

echo "[$(date '+%F %T')] SNP00157: antechamber (AM1-BCC, gaff2)"
antechamber -i ligand_posed.sdf -fi sdf -o LIG.mol2 -fo mol2 -c bcc -nc 0 \
            -at gaff2 -rn LIG -s 0 > antechamber.log 2>&1
echo "[$(date '+%F %T')] antechamber exit=$?"

if [ ! -s LIG.mol2 ]; then
  echo "FAILED: LIG.mol2 not produced"
  tail -20 antechamber.log
  exit 1
fi

echo "[$(date '+%F %T')] SNP00157: parmchk2"
parmchk2 -i LIG.mol2 -f mol2 -o LIG.frcmod -s gaff2 > parmchk2.log 2>&1
echo "[$(date '+%F %T')] parmchk2 exit=$?"

ATTN=$(grep -c ATTN LIG.frcmod || true)
echo "[$(date '+%F %T')] parmchk2 ATTN count = $ATTN   (ADR-AIIM-316 D1 predicts 0)"
grep ATTN LIG.frcmod | head -5
echo "DONE" > "$B/_ligand_done"

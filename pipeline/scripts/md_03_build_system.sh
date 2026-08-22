#!/bin/bash
# md_03_build_system.sh — scripted membrane-protein system build.
#
# No GUI at any step. PACKMOL-Memgen is used because its own publication
# positions it as the command-line alternative to web-interface builders
# (Schott-Verdugo & Gohlke, J Chem Inf Model 2019, 10.1021/acs.jcim.9b00269),
# which matters here because this project's central claim is that an agent
# executed the pipeline reproducibly.
#
# Force fields, all available in the installed AmberTools:
#   protein  ff19SB
#   lipids   Lipid21   (Dickson et al., JCTC 2022, 10.1021/acs.jctc.1c01217)
#   ligand   GAFF2 + AM1-BCC, parameters from md_02_prepare_ligand.py
#   water    OPC
#
# NOTE ON THE LIPID RATIO: POPC:CHL1 is set below, but the cholesterol mole
# fraction is still an OPEN DECISION requiring a citation (see
# 00_project_control/25_MD_PROTOCOL_DESIGN.md §2). The value used here is for
# the timing benchmark only and must not be carried into production without
# a reference.
set -u

A=/home/a2amd/A2A_MD/ambertools26_env
export AMBERHOME="$A"
export PATH="$A/bin:$PATH"

B="/mnt/d/Agentic VS Workflow/20_new_runs/md_benchmark_2026-08-21"
cd "$B/03_build" || exit 1

LIPIDS="POPC:CHL1"
RATIO="9:1"            # <-- benchmark placeholder, needs a citation
SALTCON=0.15

packmol-memgen \
  --pdb complex.pdb \
  --lipids "$LIPIDS" --ratio "$RATIO" \
  --dist 15 --dist_wat 17.5 \
  --salt --saltcon "$SALTCON" --salt_c "Na+" --salt_a "Cl-" \
  --ffprot ff19SB --fflip lipid21 --ffwat opc --gaff2 \
  --ligand_param "$B/02_ligand/ZMA.mol2:$B/02_ligand/ZMA.frcmod" \
  --pdb2pqr_pH 7.4 \
  --keepligs \
  --parametrize \
  --noprogress --keep \
  -o system.pdb \
  --log "$B/logs/packmol_memgen.log"

echo "exit=$?"
ls -la "$B/03_build" | head -20

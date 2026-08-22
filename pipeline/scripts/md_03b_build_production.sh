#!/bin/bash
# md_03b_build_production.sh — production system build.
#
# Differences from the benchmark build, each with a reason:
#   * PURE POPC, no cholesterol. Ng, Laughton & Doughty (J Chem Inf Model
#     2014;54:573-81) simulated A2AR in POPC (300 ns, triplicate) — a direct
#     precedent for this receptor. Cholesterol was considered and rejected for
#     this experiment because Yang & Lyman (Biochemistry 2019;58:4096-4105)
#     show the local lipid environment around A2AR, cholesterol in particular,
#     needs microseconds to equilibrate; a 50-100 ns pose-persistence run
#     cannot deliver that, so including it would add an unequilibrated
#     variable rather than realism.
#   * ICL3 chain break and construct termini capped with ACE/NME.
#   * 156 crystallographic waters within 5 A of the receptor reinstated.
#   * --notprotonate: pdb2pqr rejects the NME caps ("Missing atoms or failed
#     protonation for NME 211/300"), so protonation is delegated to tleap,
#     which handles ACE/NME natively. Consequence to state in Methods:
#     histidine tautomers are assigned by tleap's default (HIE, Ne-H) rather
#     than by a pH-based predictor. This is a declared choice, not an
#     oversight; His250 and His264 sit near the orthosteric pocket, so the
#     assignment should be checked against the built topology before
#     production.
#
# The benchmark's 9:1 POPC:CHL1 placeholder is gone: it had no citation.
#
# --nloop/--nloop_all raised above the defaults (20/100): at the defaults this
# system stopped at GENCAN loop 39 with the objective function still at 1.2e4
# and left 119 inter-residue contacts under 0.7 A (worst 0.117 A), which drove
# the first minimisation to an infinite force. The benchmark system, packed
# from the same defaults, had only 25 such contacts and minimised cleanly --
# so this is a convergence failure specific to this composition, not a
# systematic fault.
set -u

A=/home/a2amd/A2A_MD/ambertools26_env
export AMBERHOME="$A"
export PATH="$A/bin:$PATH"

B="/mnt/d/Agentic VS Workflow/20_new_runs/md_production_2026-08-21"
cd "$B/03_build" || exit 1

packmol-memgen \
  --pdb complex.pdb \
  --lipids POPC --ratio 1 \
  --dist 15 --dist_wat 17.5 \
  --salt --saltcon 0.15 --salt_c "Na+" --salt_a "Cl-" \
  --ffprot ff19SB --fflip lipid21 --ffwat opc --gaff2 \
  --notprotonate   \
  --keepligs \
  --nloop 60 --nloop_all 400 \
  --noprogress --keep \
  -o system.pdb \
  --log "$B/logs/packmol_memgen.log"

echo "packmol-memgen exit=$?"
ls -la system.pdb 2>/dev/null

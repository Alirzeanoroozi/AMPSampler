#!/bin/bash
# Stage 3 - interface energetics with Rosetta InterfaceAnalyzer. Needs Rosetta.
# Emits dG_separated (binding ddG proxy), dSASA_int (buried SASA), sc (shape complementarity).
set -euo pipefail
COMPLEX_DIR="${1:?dir of binder-target complex PDBs}"
OUT="${2:-../../results/stage3_validation/ddg.sc}"
$ROSETTA/bin/InterfaceAnalyzer.default.linuxgccrelease \
  -s "${COMPLEX_DIR}"/*.pdb \
  -interface A_P \
  -compute_packstat true -pack_separated true -tracer_data_print false \
  -out:file:score_only "${OUT}"
echo "Parse ${OUT}: keep favorable dG_separated, high sc, buried-SASA in mini-binder range; key by design_id."

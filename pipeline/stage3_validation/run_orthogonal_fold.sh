#!/bin/bash
# Stage 3 - orthogonal refold of binder-target complexes with an INDEPENDENT model.
# Needs GPU. Use <design_id> as the prediction name so scores join the manifest.
set -euo pipefail
TARGET="${1:?usage: run_orthogonal_fold.sh NDM5|KPC3 <designs.fasta>}"
DESIGNS="${2:?path to designed binder sequences FASTA}"
TGT_SEQ="../../targets/${TARGET}/design_domain.fasta"

# Option A: ColabFold / AF2-multimer (target+binder as a 2-chain complex per design)
#   build per-design A3M/FASTA (target : binder), then:
#   colabfold_batch complexes/ out_${TARGET}/ --num-models 5
#   -> parse iptm + pae_interaction from out_${TARGET}/*_scores.json
#
# Option B: Boltz-2 with affinity head (independent of BoltzGen design head)
#   for each design write a Boltz YAML (target chain + binder chain + 'affinity: binder: B')
#   boltz predict yamls/ --out_dir out_${TARGET} --use_msa_server
#   python parse_boltz2.py --boltz_out out_${TARGET} --out ../../results/stage3_validation/boltz2_${TARGET}.csv
echo "Refold ${DESIGNS} vs ${TGT_SEQ} with AF2-multimer or Boltz-2; keep pae_interaction<10, iptm>0.6."

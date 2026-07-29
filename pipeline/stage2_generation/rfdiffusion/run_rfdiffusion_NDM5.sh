#!/bin/bash
# RFdiffusion binder design - NDM5 (New Delhi metallo-beta-lactamase 5)
# Active-site epitope hotspots (PDB author numbering): A65,A67,A73,A93,A120,A122,A123,A124,A125,A189,A190,A208,A211,A216,A217,A218,A219,A220,A249,A250
set -euo pipefail

# 1) RFdiffusion: scaffold binders against the epitope
python $RFDIFFUSION/scripts/run_inference.py \
  inference.output_prefix=out_NDM5/NDM5 \
  inference.input_pdb=../../../targets/NDM5/structures/NDM5_target.pdb \
  'contigmap.contigs=[A43-270/0 12-45]' \
  'ppi.hotspot_res=[A65,A67,A73,A93,A120,A122,A123,A124,A125,A189,A190,A208,A211,A216,A217,A218,A219,A220,A249,A250]' \
  inference.num_designs=1000 \
  denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0

# 2) ProteinMPNN: assign sequences to the diffused backbones (binder chain only)
python $PROTEINMPNN/protein_mpnn_run.py \
  --pdb_path out_NDM5 --pdb_path_chains P \
  --num_seq_per_target 8 --sampling_temp 0.1 --out_folder out_NDM5/mpnn

# 3) AF2 filter: refold complex, keep low pAE_interaction / high pLDDT (feeds Stage 3)
#    (use e.g. dl_binder_design predict.py or ColabFold; threshold pae_interaction < 10)
echo "RFdiffusion+MPNN done; run AF2 interface filtering, then Stage 3."

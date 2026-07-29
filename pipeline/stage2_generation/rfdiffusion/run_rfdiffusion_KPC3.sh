#!/bin/bash
# RFdiffusion binder design - KPC3 (Klebsiella pneumoniae carbapenemase 3)
# Active-site epitope hotspots (PDB author numbering): A69,A70,A73,A104,A105,A130,A132,A166,A167,A170,A216,A220,A234,A235,A236,A237,A238,A274
set -euo pipefail

# 1) RFdiffusion: scaffold binders against the epitope
python $RFDIFFUSION/scripts/run_inference.py \
  inference.output_prefix=out_KPC3/KPC3 \
  inference.input_pdb=../../../targets/KPC3/structures/KPC3_target.pdb \
  'contigmap.contigs=[A30-295/0 12-45]' \
  'ppi.hotspot_res=[A69,A70,A73,A104,A105,A130,A132,A166,A167,A170,A216,A220,A234,A235,A236,A237,A238,A274]' \
  inference.num_designs=1000 \
  denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0

# 2) ProteinMPNN: assign sequences to the diffused backbones (binder chain only)
python $PROTEINMPNN/protein_mpnn_run.py \
  --pdb_path out_KPC3 --pdb_path_chains P \
  --num_seq_per_target 8 --sampling_temp 0.1 --out_folder out_KPC3/mpnn

# 3) AF2 filter: refold complex, keep low pAE_interaction / high pLDDT (feeds Stage 3)
#    (use e.g. dl_binder_design predict.py or ColabFold; threshold pae_interaction < 10)
echo "RFdiffusion+MPNN done; run AF2 interface filtering, then Stage 3."

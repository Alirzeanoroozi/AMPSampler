#!/usr/bin/env python3
"""
Stage 2 (RFdiffusion + ProteinMPNN + AF2) - emit a ready-to-run command script.

RFdiffusion binder design uses author/PDB numbering for ppi.hotspot_res (e.g. A120).
The contig is the target chain range + the binder length block. After diffusion you
sequence with ProteinMPNN and filter by AF2 (pAE_interaction etc.). Needs a GPU + the
RFdiffusion/ProteinMPNN/AF2 installs (see ../README.md); not run here.
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
BINDER_LEN = (12, 45)
NUM_DESIGNS = 1000


def chain_range(tkey, chain):
    import csv
    rows = list(csv.DictReader(open(os.path.join(REPO, "targets", tkey, "numbering_map.csv"))))
    nums = [int(r["scaffold_resnum"]) for r in rows]
    return min(nums), max(nums)


def build(tkey, tdef):
    epi = json.load(open(os.path.join(REPO, "targets", tkey, "epitope", "boltzgen_hotspots.json")))
    chain = epi["scaffold_chain"]
    hotspots = ",".join(epi["hotspot_residues_scaffold_numbering"])  # already 'A120' style
    lo, hi = chain_range(tkey, chain)
    target_pdb = os.path.relpath(os.path.join(REPO, "targets", tkey, "structures", f"{tkey}_target.pdb"), HERE)
    script = f"""#!/bin/bash
# RFdiffusion binder design - {tkey} ({tdef['name']})
# Active-site epitope hotspots (PDB author numbering): {hotspots}
set -euo pipefail

# 1) RFdiffusion: scaffold binders against the epitope
python $RFDIFFUSION/scripts/run_inference.py \\
  inference.output_prefix=out_{tkey}/{tkey} \\
  inference.input_pdb={target_pdb} \\
  'contigmap.contigs=[{chain}{lo}-{hi}/0 {BINDER_LEN[0]}-{BINDER_LEN[1]}]' \\
  'ppi.hotspot_res=[{hotspots}]' \\
  inference.num_designs={NUM_DESIGNS} \\
  denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0

# 2) ProteinMPNN: assign sequences to the diffused backbones (binder chain only)
python $PROTEINMPNN/protein_mpnn_run.py \\
  --pdb_path out_{tkey} --pdb_path_chains P \\
  --num_seq_per_target 8 --sampling_temp 0.1 --out_folder out_{tkey}/mpnn

# 3) AF2 filter: refold complex, keep low pAE_interaction / high pLDDT (feeds Stage 3)
#    (use e.g. dl_binder_design predict.py or ColabFold; threshold pae_interaction < 10)
echo "RFdiffusion+MPNN done; run AF2 interface filtering, then Stage 3."
"""
    out = os.path.join(HERE, f"run_rfdiffusion_{tkey}.sh")
    open(out, "w").write(script)
    os.chmod(out, 0o755)
    print(f"[{tkey}] contig {chain}{lo}-{hi} + {BINDER_LEN} ; hotspots {hotspots[:40]}... -> {os.path.relpath(out, REPO)}")


def main():
    targets = {k: v for k, v in json.load(open(os.path.join(REPO, "targets", "targets.json"))).items()
               if not k.startswith("_")}
    for tkey, tdef in targets.items():
        build(tkey, tdef)


if __name__ == "__main__":
    main()

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=kutem_gpu
#SBATCH --account=kutem
#SBATCH --gres=gpu:tesla_a100:1
#SBATCH --time=72:00:00
#SBATCH --mem=200G
#SBATCH --output=logs/boltzgen.out
#SBATCH --error=logs/boltzgen.out

# Initialise environment and modules
CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/bin/activate bg_macrel
export LD_LIBRARY_PATH=${CONDA_BASE}/lib

# python src/1_design/boltzgen_design.py --yaml_path boltzgen_inputs/yamls/KPC3.yaml --output_dir boltzgen_outputs/KPC3/
# python src/1_design/boltzgen_design.py --yaml_path boltzgen_inputs/yamls/NDM5.yaml --output_dir boltzgen_outputs/NDM5/

# python src/1_design/boltzgen_inverse_folding.py --design_dir boltzgen_outputs/KPC3 --output_dir boltzgen_outputs/KPC3_inverse_folding/
# python src/1_design/boltzgen_inverse_folding.py --design_dir boltzgen_outputs/NDM5 --output_dir boltzgen_outputs/NDM5_inverse_folding/

python src/1_design/boltzgen_extract_binders.py
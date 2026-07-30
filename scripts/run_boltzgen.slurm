#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=kutem_gpu
#SBATCH --account=kutem
#SBATCH --gres=gpu:tesla_a100:7
#SBATCH --time=72:00:00
#SBATCH --mem=200G
#SBATCH --output=logs/comp.out
#SBATCH --error=logs/comp.out

# Initialise environment and modules
CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/bin/activate bg
export LD_LIBRARY_PATH=${CONDA_BASE}/lib

nvidia-smi

# KPC3 on GPUs 0-2, NDM5 on GPUs 3-5 (same commands, parallel)
CUDA_VISIBLE_DEVICES=0,1,2 boltzgen run boltzgen_inputs/yamls/KPC3.yaml --output boltzgen_generations/KPC3 --protocol peptide-anything --num_designs 5000 --budget 200 &
CUDA_VISIBLE_DEVICES=3,4,5 boltzgen run boltzgen_inputs/yamls/NDM5.yaml --output boltzgen_generations/NDM5 --protocol peptide-anything --num_designs 5000 --budget 200 &
wait
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=kutem_gpu
#SBATCH --account=kutem
#SBATCH --gres=gpu:tesla_a100:1
#SBATCH --time=72:00:00
#SBATCH --mem=200G
#SBATCH --output=logs/3_amp_filter.out
#SBATCH --error=logs/3_amp_filter.out

# Initialise environment and modules
CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/bin/activate bg_macrel
export LD_LIBRARY_PATH=${CONDA_BASE}/lib

python src/3_amp_filter/developability.py --fasta fastas/NDM5.fasta --out results/developability_NDM5.csv
python src/3_amp_filter/developability.py --fasta fastas/KPC3.fasta --out results/developability_KPC3.csv

python src/3_amp_filter/plot_distributions.py

python src/3_amp_filter/create_manifest.py

python src/3_amp_filter/apply_filters.py --manifest results/manifest_NDM5.csv --out results/filtered_NDM5.csv
python src/3_amp_filter/apply_filters.py --manifest results/manifest_KPC3.csv --out results/filtered_KPC3.csv
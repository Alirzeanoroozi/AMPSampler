#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=kutem_gpu
#SBATCH --account=kutem
#SBATCH --gres=gpu:tesla_a100:1
#SBATCH --time=72:00:00
#SBATCH --mem=200G
#SBATCH --output=logs/4_boltz2.out
#SBATCH --error=logs/4_boltz2.out

CONDA_BASE="$(conda info --base)"
source "${CONDA_BASE}/bin/activate" latest_bg
export LD_LIBRARY_PATH="${CONDA_BASE}/lib:${LD_LIBRARY_PATH:-}"

# python src/4_structure_prediction/prepare_configs.py
# boltz predict boltz_inputs --out_dir boltz_results --write_full_pae

PRED=boltz_results/boltz_results_boltz_inputs/predictions

python src/4_structure_prediction/parse_boltz2.py --boltz_out boltz_results --target NDM5
python src/4_structure_prediction/parse_boltz2.py --boltz_out boltz_results --target KPC3

python src/4_structure_prediction/active_site_overlap.py --target NDM5 --complexes "${PRED}"
python src/4_structure_prediction/active_site_overlap.py --target KPC3 --complexes "${PRED}"

python src/4_structure_prediction/run_ipsae.py --boltz-out boltz_results --workers 8

python src/4_structure_prediction/build_manifest.py

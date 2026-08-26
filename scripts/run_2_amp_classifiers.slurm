#!/bin/bash
#SBATCH --job-name=classifiers
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --partition=kutem_gpu
#SBATCH --account=kutem
#SBATCH --output=logs/classifiers.out
#SBATCH --error=logs/classifiers.out
#
# Run AMP Scanner, Macrel, HydrAMP, APEX, and ToxinPred on FASTA files.
#
# Tunables (override at submit time):
#   FASTAS_DIR (default: fastas)          directory of input .fasta files
#   OUTPUT_DIR (default: results/classifiers)
#   ORCHESTRATOR_ENV (default: gmconda_py3923)  conda env for run_classifiers.py
#
# PREREQUISITE done on a LOGIN NODE (compute nodes usually have no internet):
#   cd src/2_amp_classifiers && ./download_models.sh && ./setup_envs.sh

FASTAS_DIR="${FASTAS_DIR:-fastas}"
OUTPUT_DIR="${OUTPUT_DIR:-results/classifiers}"

REPO_DIR="$(pwd)"
mkdir -p "${REPO_DIR}/${OUTPUT_DIR}"

# Initialise environment and modules
CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/bin/activate gmconda_py3923
export LD_LIBRARY_PATH=${CONDA_BASE}/lib

FASTA_FILES=( "${REPO_DIR}/${FASTAS_DIR}"/*.fasta )

cd "${REPO_DIR}/src/2_amp_classifiers"
for fasta in "${FASTA_FILES[@]}"; do
    echo "=== $(basename "${fasta}") ($(date)) ==="
    ./run_classifiers.sh \
        --fasta "${fasta}" \
        --output-dir "${REPO_DIR}/${OUTPUT_DIR}"
done

echo "Merging classifier outputs ..."
python "${REPO_DIR}/src/2_amp_classifiers/merge_classifier_results.py" \
    --results-dir "${REPO_DIR}/${OUTPUT_DIR}" \
    --fastas-dir "${REPO_DIR}/${FASTAS_DIR}" \
    --output "${REPO_DIR}/${OUTPUT_DIR}/merged_classifier_results.csv"

echo "Classifiers finished at $(date). Outputs in ${REPO_DIR}/${OUTPUT_DIR}/"

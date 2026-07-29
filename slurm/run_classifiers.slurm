#!/bin/bash
#SBATCH --job-name=classifiers
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --partition=kutem_gpu
#SBATCH --account=kutem
#SBATCH --output=slurm/logs/classifiers-%j.out
#SBATCH --error=slurm/logs/classifiers-%j.out
#
# Run AMP Scanner, Macrel, HydrAMP, APEX, and ToxinPred on FASTA files.
# SUBMIT FROM THE REPO ROOT:  sbatch slurm/run_classifiers.slurm
#
# Tunables (override at submit time):
#   FASTAS_DIR (default: fastas)          directory of input .fasta files
#   OUTPUT_DIR (default: results/classifiers)
#   ORCHESTRATOR_ENV (default: gmconda_py3923)  conda env for run_classifiers.py
#
# PREREQUISITE done on a LOGIN NODE (compute nodes usually have no internet):
#   cd classifiers && ./download_models.sh && ./setup_envs.sh
set -euo pipefail

FASTAS_DIR="${FASTAS_DIR:-fastas}"
OUTPUT_DIR="${OUTPUT_DIR:-results/classifiers}"
ORCHESTRATOR_ENV="${ORCHESTRATOR_ENV:-gmconda_py3923}"

REPO_DIR="$(pwd)"
mkdir -p "${REPO_DIR}/slurm/logs" "${REPO_DIR}/${OUTPUT_DIR}"

echo "========== Slurm Job =========="
echo "Job ID: ${SLURM_JOB_ID:-N/A} | Node: $(hostname) | Date: $(date)"
echo "Repo: ${REPO_DIR} | fastas=${FASTAS_DIR} | output=${OUTPUT_DIR}"
echo "==============================="

if command -v conda >/dev/null 2>&1; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${ORCHESTRATOR_ENV}"
else
    echo "conda not found" >&2
    exit 1
fi

export PYTHONUNBUFFERED=1

FASTAS_PATH="${REPO_DIR}/${FASTAS_DIR}"
if [[ ! -d "${FASTAS_PATH}" ]]; then
    echo "FASTA directory not found: ${FASTAS_PATH}" >&2
    exit 1
fi

shopt -s nullglob
FASTA_FILES=( "${FASTAS_PATH}"/*.fasta )
if [[ ${#FASTA_FILES[@]} -eq 0 ]]; then
    echo "No .fasta files found in ${FASTAS_PATH}" >&2
    exit 1
fi

cd "${REPO_DIR}/classifiers"
for fasta in "${FASTA_FILES[@]}"; do
    echo "=== $(basename "${fasta}") ($(date)) ==="
    ./run_classifiers.sh \
        --fasta "${fasta}" \
        --output-dir "${REPO_DIR}/${OUTPUT_DIR}"
done

echo "Merging classifier outputs ..."
python "${REPO_DIR}/classifiers/merge_classifier_results.py" \
    --results-dir "${REPO_DIR}/${OUTPUT_DIR}" \
    --fastas-dir "${REPO_DIR}/${FASTAS_DIR}" \
    --output "${REPO_DIR}/${OUTPUT_DIR}/merged_classifier_results.csv"

echo "Classifiers finished at $(date). Outputs in ${REPO_DIR}/${OUTPUT_DIR}/"

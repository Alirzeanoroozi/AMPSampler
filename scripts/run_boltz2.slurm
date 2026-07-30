#!/bin/bash
#SBATCH --job-name=boltzgen_inhib
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --partition=kutem_gpu
#SBATCH --account=kutem
#SBATCH --gres=gpu:1
#SBATCH --output=slurm/logs/boltzgen_inhib-%j.out
#SBATCH --error=slurm/logs/boltzgen_inhib-%j.out
#
# Stage 2 - epitope-conditioned BoltzGen binder design for NDM5 and KPC3.
# SUBMIT FROM THE REPO ROOT:  sbatch slurm/run_boltzgen.slurm
#
# Tunables (override at submit time, e.g. `NUM_DESIGNS=2000 BUDGET=100 sbatch ...`):
#   CONDA_ENV_NAME (default: boltzgen)   NUM_DESIGNS (5000)   BUDGET (200)   TARGETS ("NDM5 KPC3")
#
# PREREQUISITE done on a LOGIN NODE (compute nodes usually have no internet):
#   1) conda env with BoltzGen installed
#   2) BoltzGen model weights pre-downloaded/cached (first `boltzgen run` downloads them)
#   3) targets/<T>/structures/<T>_target.pdb and pipeline/stage2_generation/boltzgen/<T>_design.yaml
#      present (committed with the repo, or regenerated via Stages 0-2 on the login node)
set -euo pipefail

CONDA_ENV_NAME="${CONDA_ENV_NAME:-boltzgen}"
NUM_DESIGNS="${NUM_DESIGNS:-5000}"
BUDGET="${BUDGET:-200}"
TARGETS="${TARGETS:-NDM5 KPC3}"

REPO_DIR="$(pwd)"                       # submission dir = repo root
mkdir -p "${REPO_DIR}/slurm/logs"

echo "========== Slurm Job =========="
echo "Job ID: ${SLURM_JOB_ID:-N/A} | Node: $(hostname) | Date: $(date)"
echo "Repo: ${REPO_DIR} | env: ${CONDA_ENV_NAME} | num_designs=${NUM_DESIGNS} budget=${BUDGET}"
echo "==============================="

if command -v conda >/dev/null 2>&1; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${CONDA_ENV_NAME}"
else
    echo "conda not found" >&2; exit 1
fi
command -v boltzgen >/dev/null 2>&1 || { echo "boltzgen CLI not found in env '${CONDA_ENV_NAME}'" >&2; exit 1; }
command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi || echo "nvidia-smi not found; continuing"
export PYTHONUNBUFFERED=1

cd "${REPO_DIR}/pipeline/stage2_generation/boltzgen"   # yaml paths are relative to this dir
for T in ${TARGETS}; do
    echo "=== BoltzGen: ${T} ($(date)) ==="
    boltzgen run "${T}_design.yaml" \
        --output "out_${T}" \
        --num_designs "${NUM_DESIGNS}" \
        --budget "${BUDGET}"
    echo "=== ${T} done; designs in out_${T}/final_ranked_designs/ ==="
done
echo "Stage 2 finished at $(date). Next: pipeline/run_downstream.sh (CPU)."

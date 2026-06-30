#!/usr/bin/env bash
# Run toxinpred3 on every .fasta in FASTAS_DIR.
# Usage: from repo root or from src/ — paths are relative to this script.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
# Directory containing input .fasta files (use ../genbio/fastas if you keep them in a subfolder).
FASTAS_DIR="${SCRIPT_DIR}/../genbio"
OUTPUT_DIR="${SCRIPT_DIR}/../genbio/toxinpred"
VENV_DIR="${REPO_ROOT}/toxinpred-venv"
mkdir -p "$OUTPUT_DIR"

if [[ ! -d "$FASTAS_DIR" ]]; then
  echo "Input directory not found: ${FASTAS_DIR}" >&2
  exit 1
fi

# Use an isolated env and pin ABI-compatible deps for toxinpred3 model loading.
if [[ ! -d "$VENV_DIR" ]]; then
  python3 -m venv "$VENV_DIR"
fi
source "${VENV_DIR}/bin/activate"

python -m pip install --upgrade pip setuptools wheel >/dev/null
# sklearn==1.0.2 requires NumPy < 2 for binary compatibility.
python -m pip install --upgrade "numpy<2" "scikit-learn==1.0.2" toxinpred3 >/dev/null

shopt -s nullglob
FASTA_FILES=( "${FASTAS_DIR}"/*.fasta )

if [[ ${#FASTA_FILES[@]} -eq 0 ]]; then
  echo "No .fasta files found in ${FASTAS_DIR}" >&2
  exit 1
fi

for FASTA in "${FASTA_FILES[@]}"; do
  BASENAME=$(basename "$FASTA" .fasta)
  echo "Running toxinpred3: ${FASTA} -> ${OUTPUT_DIR}/${BASENAME}_toxinpred.csv"
  "${VENV_DIR}/bin/toxinpred3" -i "$FASTA" -o "${OUTPUT_DIR}/${BASENAME}_toxinpred.csv"
done

echo "Done (${#FASTA_FILES[@]} file(s))."

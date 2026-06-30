#!/usr/bin/env bash
# Run HydrAMP AMP + MIC classifiers on every FASTA under genbio/ (same as predict_sample_args.sh).
# Usage: from the hydramp/ directory:
#   chmod +x run_predict_genbio.sh && ./run_predict_genbio.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

GENBIO_DIR="${SCRIPT_DIR}/genbio"
OUT_DIR="${GENBIO_DIR}/predictions"
mkdir -p "$OUT_DIR"

shopt -s nullglob
fastas=( "${GENBIO_DIR}"/*.fasta )
if [[ ${#fastas[@]} -eq 0 ]]; then
  echo "No .fasta files found in ${GENBIO_DIR}" >&2
  exit 1
fi

for fasta in "${fastas[@]}"; do
  base="$(basename "${fasta}" .fasta)"
  echo "=== ${base} ==="
  python -m amp.inference.scripts.predict_if_amp \
    --model_path models/amp_classifier/ \
    --sequence_path "${fasta}" \
    --format fasta \
    --output_csv "${OUT_DIR}/${base}_amp.csv"

  python -m amp.inference.scripts.predict_if_amp \
    --model_path models/mic_classifier/ \
    --sequence_path "${fasta}" \
    --format fasta \
    --output_csv "${OUT_DIR}/${base}_mic.csv"
done

echo "Done. Outputs under ${OUT_DIR}/"

#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GENBIO_DIR="${SCRIPT_DIR}/genbio"
MODEL_PATH="${SCRIPT_DIR}/trained-models/OriginalPaper_081917_FULL_MODEL.h5"
PREDICT_SCRIPT="${SCRIPT_DIR}/amp_scanner_v2_predict_tf1.py"

if [[ ! -d "$GENBIO_DIR" ]]; then
  echo "Input directory not found: $GENBIO_DIR" >&2
  exit 1
fi

if [[ ! -f "$MODEL_PATH" ]]; then
  echo "Model not found: $MODEL_PATH" >&2
  exit 1
fi

shopt -s nullglob
FASTA_FILES=("${GENBIO_DIR}"/*.fasta)
if [[ ${#FASTA_FILES[@]} -eq 0 ]]; then
  echo "No FASTA files found in $GENBIO_DIR" >&2
  exit 1
fi

for FASTA in "${FASTA_FILES[@]}"; do
  BASENAME="$(basename "$FASTA" .fasta)"
  CANDIDATES_PATH="${GENBIO_DIR}/${BASENAME}_candidates.fasta"
  PREDS_PATH="${GENBIO_DIR}/${BASENAME}.csv"

  echo "Running AMP Scanner for ${BASENAME}"
  python "$PREDICT_SCRIPT" \
    -fasta "$FASTA" \
    -model "$MODEL_PATH" \
    -candidates "$CANDIDATES_PATH" \
    -preds "$PREDS_PATH"
done

echo "Done. Processed ${#FASTA_FILES[@]} FASTA files in $GENBIO_DIR"

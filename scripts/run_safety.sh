#!/bin/bash
# Stage 4 - safety scores (external tools). Emit design_id-keyed CSVs for the manifest.
set -euo pipefail
FASTA="${1:?usage: run_safety.sh designs.fasta TARGET}"
T="${2:?target tag, e.g. NDM5}"
OUT="../../results/stage4_filtering"; mkdir -p "$OUT"

# ToxinPred3 (safety): pip install toxinpred3
toxinpred3 -i "$FASTA" -o "$OUT/toxinpred_${T}.csv"     # keep Prediction == Non-Toxin

# Macrel Hemo head (hemolysis safety) - use classifiers/run_classifiers.sh macrel output
# (hemo_prob column in *_macrel.csv); keep hemo_prob < 0.5.
echo "Wrote $OUT/toxinpred_${T}.csv ; compute hemolysis via macrel Hemo head (hemo_prob)."

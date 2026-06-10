#!/bin/bash
pip3 install toxinpred3
# pip install scikit-learn==1.0.2

FASTAS_DIR="../data/fastas"
OUTPUT_DIR="../results/toxinpred"
mkdir -p "$OUTPUT_DIR"

FASTA_FILES=(
    "boltzgen_X_hallucination.fasta"
    "boltzgen_KPC3_binders.fasta"
    "boltzgen_NDM5_binders.fasta"
    "dplm1.fasta"
    "dplm2.fasta"
    "random_peptides.fasta"
)

for FASTA in "${FASTA_FILES[@]}"; do
    BASENAME=$(basename "$FASTA" .fasta)
    # echo "Running toxinpred3 on $FASTA -> $OUTFILE"
    toxinpred3 -i "$FASTAS_DIR/$FASTA" -o "$OUTPUT_DIR/${BASENAME}_toxinpred.csv"
done
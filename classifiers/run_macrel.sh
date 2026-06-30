#!/bin/bash

FASTAS_DIR="../data/fastas"

# List of fasta files to process (update/add as needed)
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
    OUTDIR="../results/${BASENAME}_macrel"

    rm -rf "$OUTDIR"
    macrel peptides --fasta "$FASTAS_DIR/$FASTA" --output "$OUTDIR" --keep-negatives
    gzip -d "$OUTDIR/macrel.out.prediction.gz"
    tail -n +2 "$OUTDIR/macrel.out.prediction" > "$OUTDIR/macrel.csv"
    rm "$OUTDIR/macrel.out.prediction"

    # Sort macrel.csv by AMP_probability (5th column)
    (head -n 1 "$OUTDIR/macrel.csv" && tail -n +2 "$OUTDIR/macrel.csv" | sort -t$'\t' -k5,5gr) > "$OUTDIR/macrel.sorted.csv"
    mv "$OUTDIR/macrel.sorted.csv" "$OUTDIR/macrel.csv"
done


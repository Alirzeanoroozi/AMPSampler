#!/usr/bin/env python3
"""
Merge HydrAMP AMP/MIC prediction CSVs into genbio/*.csv (root only).

Predictions are aligned to sequences using the same rules as predict_if_amp:
standard amino acids only, length in [MIN_LENGTH, MAX_LENGTH] (amp.config).

Writes ``<name>_with_hydramp.csv`` next to each input CSV (does not overwrite).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from amp.config import MAX_LENGTH, MIN_LENGTH
from amp.data_utils.sequence import check_if_std_aa


def prediction_order_sequences(fasta_path: Path) -> list[str]:
    """Sequences in the order passed to the model (matches predict_if_amp)."""
    out: list[str] = []
    for record in SeqIO.parse(str(fasta_path), "fasta"):
        seq = str(record.seq)
        if not check_if_std_aa(seq):
            continue
        if not (MIN_LENGTH <= len(seq) <= MAX_LENGTH):
            continue
        out.append(seq)
    return out


def load_prediction_column(path: Path) -> pd.Series:
    df = pd.read_csv(path)
    if "Prediction" not in df.columns:
        raise ValueError(f"{path}: expected a 'Prediction' column, got {list(df.columns)}")
    return df["Prediction"]


def merge_one(
    genbio_dir: Path,
    pred_dir: Path,
    csv_path: Path,
    dry_run: bool,
) -> None | str:
    stem = csv_path.stem
    fasta = genbio_dir / f"{stem}.fasta"
    amp_csv = pred_dir / f"{stem}_amp.csv"
    mic_csv = pred_dir / f"{stem}_mic.csv"
    out_path = genbio_dir / f"{stem}_with_hydramp.csv"

    for p in (fasta, amp_csv, mic_csv):
        if not p.exists():
            return f"skip {csv_path.name}: missing {p}"

    seqs = prediction_order_sequences(fasta)
    amp = load_prediction_column(amp_csv)
    mic = load_prediction_column(mic_csv)

    if len(amp) != len(mic):
        return f"error {stem}: amp rows ({len(amp)}) != mic rows ({len(mic)})"
    if len(seqs) != len(amp):
        return (
            f"error {stem}: fasta-filtered sequences ({len(seqs)}) != predictions ({len(amp)}). "
            "Re-run predictions or check FASTA matches the run used for preds."
        )

    pred_map = pd.DataFrame(
        {
            "Sequence": seqs,
            "hydramp_amp": amp.values,
            "hydramp_mic": mic.values,
        }
    ).drop_duplicates(subset=["Sequence"], keep="last")

    main = pd.read_csv(csv_path)
    if "Sequence" not in main.columns:
        return f"error {stem}: no 'Sequence' column in {csv_path.name}"

    merged = main.merge(pred_map, on="Sequence", how="left")
    if dry_run:
        print(f"OK {stem}: rows {len(merged)} -> {out_path.name}")
        return None

    merged.to_csv(out_path, index=False)
    print(f"Wrote {out_path} ({len(merged)} rows)")
    return None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "--genbio_dir",
        type=Path,
        default=Path(__file__).resolve().parent / "genbio",
        help="Directory containing genbio CSVs and FASTAs (default: hydramp/genbio)",
    )
    parser.add_argument(
        "--pred_dir",
        type=Path,
        default=None,
        help="Predictions directory (default: <genbio_dir>/predictions)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Check alignment only; do not write files",
    )
    args = parser.parse_args()
    genbio_dir = args.genbio_dir.resolve()
    pred_dir = (args.pred_dir or (genbio_dir / "predictions")).resolve()

    csvs = sorted(genbio_dir.glob("*.csv"))
    if not csvs:
        print(f"No CSV files in {genbio_dir}", file=sys.stderr)
        sys.exit(1)

    errs: list[str] = []
    for csv_path in csvs:
        # Only merge "source" tables; skip our own outputs and anything in subdirs already excluded
        if csv_path.name.endswith("_with_hydramp.csv"):
            continue
        msg = merge_one(genbio_dir, pred_dir, csv_path, args.dry_run)
        if msg:
            errs.append(msg)
            print(msg)

    if errs and any(e.startswith("error") for e in errs):
        sys.exit(2)


if __name__ == "__main__":
    main()

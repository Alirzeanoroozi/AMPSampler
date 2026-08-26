#!/usr/bin/env python3
"""Merge per-classifier CSV outputs into one table keyed by sequence.

HydrAMP outputs only a single ``Prediction`` column with no sequence IDs.
Predictions are aligned to FASTA order using the same filters as predict_if_amp:
standard amino acids only, length in [MIN_LENGTH, MAX_LENGTH] (see amp.config).
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import pandas as pd

CLASSIFIERS_DIR = Path(__file__).resolve().parent
REPO_ROOT = CLASSIFIERS_DIR.parent

STD_AA = set("ACDEFGHIKLMNPQRSTVWY")
HYDRAMP_MIN_LEN = 0
HYDRAMP_MAX_LEN = 25
APEX_MAX_LEN = 50

FASTA_SPECS = (
    ("KPC3.fasta", "KPC3"),
    ("NDM5.fasta", "NDM5"),
)


def read_fasta_records(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    seq_id = None
    parts: list[str] = []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_id is not None:
                    records.append((seq_id, "".join(parts).upper()))
                seq_id = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.upper())
    if seq_id is not None:
        records.append((seq_id, "".join(parts)))
    return records


def hydramp_eligible(sequence: str) -> bool:
    return all(aa in STD_AA for aa in sequence) and HYDRAMP_MIN_LEN <= len(sequence) <= HYDRAMP_MAX_LEN


def hydramp_eligible_ids(records: list[tuple[str, str]]) -> list[str]:
    return [seq_id for seq_id, seq in records if hydramp_eligible(seq)]


def load_hydramp_predictions(csv_path: Path, eligible_ids: list[str]) -> dict[str, float]:
    if not csv_path.is_file():
        return {}
    preds = pd.read_csv(csv_path)["Prediction"].tolist()
    if len(preds) != len(eligible_ids):
        raise ValueError(
            f"HydrAMP row count mismatch for {csv_path.name}: "
            f"{len(preds)} predictions vs {len(eligible_ids)} eligible sequences"
        )
    return dict(zip(eligible_ids, preds))


def load_target_table(
    fasta_path: Path,
    target: str,
    results_dir: Path,
    stem: str,
) -> pd.DataFrame:
    records = read_fasta_records(fasta_path)
    rows = [
        {
            "seq_id": seq_id,
            "sequence": seq,
            "target": target,
            "length": len(seq),
        }
        for seq_id, seq in records
    ]
    df = pd.DataFrame(rows)

    ampscanner = pd.read_csv(results_dir / f"{stem}_ampscanner.csv")
    df = df.merge(
        ampscanner.rename(
            columns={
                "Prediction_Class": "ampscanner_class",
                "Prediction_Probability": "ampscanner_prob",
            }
        )[["SeqID", "ampscanner_class", "ampscanner_prob"]],
        left_on="seq_id",
        right_on="SeqID",
        how="left",
    ).drop(columns=["SeqID"])

    macrel = pd.read_csv(results_dir / f"{stem}_macrel.csv", comment="#")
    macrel = macrel.rename(
        columns={
            "Access": "seq_id",
            "is_AMP": "macrel_is_amp",
            "AMP_probability": "macrel_amp_prob",
            "Hemolytic": "macrel_hemolytic",
            "Hemolytic_probability": "macrel_hemo_prob",
        }
    )
    df = df.merge(
        macrel[
            [
                "seq_id",
                "macrel_is_amp",
                "macrel_amp_prob",
                "macrel_hemolytic",
                "macrel_hemo_prob",
            ]
        ],
        on="seq_id",
        how="left",
    )

    apex_path = results_dir / f"{stem}_apex.csv"
    if apex_path.is_file():
        apex = pd.read_csv(apex_path).rename(columns={"apex_mic_mean_uM": "apex_mic_mean_uM"})
        df = df.merge(apex[["seq_id", "apex_mic_mean_uM"]], on="seq_id", how="left")
    else:
        df["apex_mic_mean_uM"] = pd.NA

    toxin = pd.read_csv(results_dir / f"{stem}_toxinpred.csv")
    df = df.merge(
        toxin.rename(
            columns={
                "Subject": "seq_id",
                "Prediction": "toxinpred_class",
                "Hybrid Score": "toxinpred_hybrid_score",
                "ML Score": "toxinpred_ml_score",
            }
        )[["seq_id", "toxinpred_class", "toxinpred_hybrid_score", "toxinpred_ml_score"]],
        on="seq_id",
        how="left",
    )

    eligible_ids = hydramp_eligible_ids(records)
    amp_preds = load_hydramp_predictions(results_dir / f"{stem}_hydramp_amp.csv", eligible_ids)
    mic_preds = load_hydramp_predictions(results_dir / f"{stem}_hydramp_mic.csv", eligible_ids)

    df["hydramp_amp_prob"] = df["seq_id"].map(amp_preds)
    df["hydramp_mic_prob"] = df["seq_id"].map(mic_preds)
    df["hydramp_scored"] = df["seq_id"].isin(eligible_ids)
    df["hydramp_skip_reason"] = pd.NA
    too_long = (df["length"] > HYDRAMP_MAX_LEN) & df["hydramp_amp_prob"].isna()
    non_std = ~df["sequence"].map(lambda s: all(aa in STD_AA for aa in s))
    df.loc[too_long, "hydramp_skip_reason"] = f"length>{HYDRAMP_MAX_LEN}"
    df.loc[non_std & df["hydramp_amp_prob"].isna(), "hydramp_skip_reason"] = "non_standard_aa"

    df["apex_scored"] = df["length"] <= APEX_MAX_LEN
    df.loc[df["length"] > APEX_MAX_LEN, "apex_mic_mean_uM"] = pd.NA

    return df


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=REPO_ROOT / "results" / "classifiers",
        help="Directory containing per-classifier CSV outputs",
    )
    parser.add_argument(
        "--fastas-dir",
        type=Path,
        default=REPO_ROOT / "fastas",
        help="Directory containing input FASTA files",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_ROOT / "results" / "classifiers" / "merged_classifier_results.csv",
        help="Merged output CSV path",
    )
    args = parser.parse_args()

    frames = []
    for fasta_name, target in FASTA_SPECS:
        fasta_path = args.fastas_dir / fasta_name
        if not fasta_path.is_file():
            print(f"FASTA not found: {fasta_path}", file=sys.stderr)
            return 1
        stem = fasta_path.stem
        frames.append(load_target_table(fasta_path, target, args.results_dir.resolve(), stem))

    merged = pd.concat(frames, ignore_index=True)
    column_order = [
        "seq_id",
        "target",
        "sequence",
        "length",
        "ampscanner_class",
        "ampscanner_prob",
        "macrel_is_amp",
        "macrel_amp_prob",
        "macrel_hemolytic",
        "macrel_hemo_prob",
        "hydramp_scored",
        "hydramp_skip_reason",
        "hydramp_amp_prob",
        "hydramp_mic_prob",
        "apex_scored",
        "apex_mic_mean_uM",
        "toxinpred_class",
        "toxinpred_hybrid_score",
        "toxinpred_ml_score",
    ]
    merged = merged[column_order]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.output, index=False)

    n = len(merged)
    hydramp_scored = int(merged["hydramp_scored"].sum())
    print(f"Wrote {args.output} ({n} sequences)")
    print(
        f"HydrAMP scored {hydramp_scored}/{n} sequences "
        f"({n - hydramp_scored} skipped: length > {HYDRAMP_MAX_LEN} aa)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

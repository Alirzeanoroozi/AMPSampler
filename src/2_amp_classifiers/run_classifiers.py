#!/usr/bin/env python3
"""Run all AMPSampler classifiers on a single FASTA file.

Each classifier runs in its own conda environment. Outputs are written as CSV
files under the chosen output directory.

Usage:
  python run_classifiers.py --fasta ../fastas/amp_peptide.fasta
  python run_classifiers.py --fasta test_sample.fasta --classifiers macrel apex
"""
import argparse
import csv
import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

CLASSIFIERS_DIR = Path(__file__).resolve().parent
REPO_ROOT = CLASSIFIERS_DIR.parent

ALL_CLASSIFIERS = ("ampscanner", "macrel", "hydramp_amp", "hydramp_mic", "apex", "toxinpred")


def conda_run(env, *cmd, **kwargs):
    cwd = kwargs.get("cwd")
    full = ["conda", "run", "-n", env] + list(cmd)
    subprocess.run(full, cwd=str(cwd) if cwd else None, check=True)


def read_fasta_records(path):
    records = []
    seq_id = None
    seq_parts = []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_id is not None:
                    records.append((seq_id, "".join(seq_parts)))
                seq_id = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.upper())
    if seq_id is not None:
        records.append((seq_id, "".join(seq_parts)))
    return records


def run_ampscanner(fasta, out_dir):
    env = "ampsampler-ascan2"
    model = CLASSIFIERS_DIR / "amp-scanner-v2/trained-models/OriginalPaper_081917_FULL_MODEL.h5"
    script = CLASSIFIERS_DIR / "amp-scanner-v2/amp_scanner_v2_predict_tf1.py"
    if not model.is_file() or model.stat().st_size == 0:
        raise FileNotFoundError(f"AMP Scanner model missing: {model}")

    with tempfile.TemporaryDirectory(prefix="ascan_") as tmp:
        tmp = Path(tmp)
        preds = tmp / "preds.csv"
        candidates = tmp / "candidates.fasta"
        conda_run(
            env,
            "python",
            str(script),
            "-fasta",
            str(fasta.resolve()),
            "-model",
            str(model.resolve()),
            "-candidates",
            str(candidates),
            "-preds",
            str(preds),
            cwd=CLASSIFIERS_DIR / "amp-scanner-v2",
        )
        out_csv = out_dir / f"{fasta.stem}_ampscanner.csv"
        shutil.copy2(preds, out_csv)
    return out_csv


def run_macrel(fasta, out_dir):
    env = "ampsampler-macrel"
    work = out_dir / "{}_{}_macrel".format(fasta.stem, os.getpid())
    if work.exists():
        shutil.rmtree(work)

    conda_run(
        env,
        "macrel",
        "peptides",
        "--fasta",
        str(fasta.resolve()),
        "--output",
        str(work.resolve()),
        "--keep-negatives",
        "--force",
    )

    pred_gz = work / "macrel.out.prediction.gz"
    pred_txt = work / "macrel.out.prediction"
    with gzip.open(pred_gz, "rt") as src, pred_txt.open("w") as dst:
        dst.write(src.read())

    out_csv = out_dir / f"{fasta.stem}_macrel.csv"
    with pred_txt.open() as src, out_csv.open("w", newline="") as dst:
        reader = csv.reader(src, delimiter="\t")
        writer = csv.writer(dst)
        for row in reader:
            writer.writerow(row)
    shutil.rmtree(work)
    return out_csv


def run_hydramp(fasta, out_dir, classifier):
    env = "ampsampler-hydramp"
    hydramp_dir = CLASSIFIERS_DIR / "hydramp"
    model_subdir = "amp_classifier" if classifier == "hydramp_amp" else "mic_classifier"
    model_path = hydramp_dir / "models" / model_subdir
    if not (model_path / "model_config.json").is_file():
        raise FileNotFoundError(f"HydrAMP model missing: {model_path}")

    suffix = "hydramp_amp" if classifier == "hydramp_amp" else "hydramp_mic"
    out_csv = out_dir / f"{fasta.stem}_{suffix}.csv"

    conda_run(
        env,
        "python",
        "-m",
        "amp.inference.scripts.predict_if_amp",
        "--model_path",
        str(model_path.resolve()),
        "--sequence_path",
        str(fasta.resolve()),
        "--format",
        "fasta",
        "--output_csv",
        str(out_csv.resolve()),
        cwd=hydramp_dir,
    )
    return out_csv


def run_apex(fasta, out_dir):
    env = "ampsampler-apex"
    apex_dir = CLASSIFIERS_DIR / "apex"
    models_dir = apex_dir / "APEX_pathogen_models"
    if not models_dir.is_dir() or not any(models_dir.glob("APEX_*")):
        raise FileNotFoundError(f"APEX models missing under {models_dir}")

    runner = out_dir / "_apex_runner.py"
    out_csv = out_dir / f"{fasta.stem}_apex.csv"
    runner.write_text(
        f"""import csv
import glob
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from Bio import SeqIO

sys.path.insert(0, {str(apex_dir)!r})
from APEX_models import AMP_model  # noqa: F401
from utils import make_vocab, onehot_encoding

fasta_path = {str(fasta.resolve())!r}
out_csv = {str(out_csv.resolve())!r}
models_dir = {str(models_dir.resolve())!r}

pathogen_list = [
    'A. baumannii ATCC 19606',
    'E. coli ATCC 11775',
    'E. coli AIC221',
    'E. coli AIC222',
    'K. pneumoniae ATCC 13883',
    'P. aeruginosa PA01',
    'P. aeruginosa PA14',
    'S. aureus ATCC 12600',
    'S. aureus (ATCC BAA-1556) - MRSA',
    'vancomycin-resistant E. faecalis ATCC 700802',
    'vancomycin-resistant E. faecium ATCC 700221',
]

max_len = 52
word2idx, idx2word = make_vocab()

models = []
for model_path in sorted(glob.glob(str(Path(models_dir) / 'APEX_*'))):
    model = torch.load(model_path, map_location='cpu', weights_only=False)
    model.eval()
    models.append(model)

seq_ids = []
seq_list = []
for record in SeqIO.parse(fasta_path, 'fasta'):
    sequence = str(record.seq).upper()
    if len(sequence) > 50:
        continue
    seq_ids.append(record.id)
    seq_list.append(sequence)

if not seq_list:
    raise SystemExit('No sequences <= 50 aa for APEX')

batch_size = 256
seq_arr = np.array(seq_list)
pred_sum = None
for model in models:
    batch_preds = []
    for i in range(int(math.ceil(len(seq_arr) / float(batch_size)))):
        batch = seq_arr[i * batch_size : (i + 1) * batch_size]
        x = torch.LongTensor(onehot_encoding(batch, max_len, word2idx))
        mic = model(x).detach().numpy()
        mic = 10 ** (6 - mic)
        batch_preds.append(mic)
    model_pred = np.vstack(batch_preds)
    pred_sum = model_pred if pred_sum is None else pred_sum + model_pred

avg_pred = pred_sum / float(len(models))
mean_mic = avg_pred.mean(axis=1)

rows = []
for seq_id, seq, mic in zip(seq_ids, seq_list, mean_mic):
    rows.append({{
        'seq_id': seq_id,
        'sequence': seq,
        'apex_mic_mean_uM': float(mic),
    }})

pd.DataFrame(rows).to_csv(out_csv, index=False)
print(f'Wrote {{out_csv}} ({{len(rows)}} sequences)')
"""
    )
    conda_run(env, "python", str(runner))
    if runner.exists():
        runner.unlink()
    return out_csv


def run_toxinpred(fasta, out_dir):
    env = "ampsampler-toxinpred"
    out_csv = out_dir / f"{fasta.stem}_toxinpred.csv"
    conda_run(
        env,
        "toxinpred3",
        "-i",
        str(fasta.resolve()),
        "-o",
        str(out_csv.resolve()),
    )
    return out_csv


RUNNERS = {
    "ampscanner": run_ampscanner,
    "macrel": run_macrel,
    "hydramp_amp": lambda f, o: run_hydramp(f, o, "hydramp_amp"),
    "hydramp_mic": lambda f, o: run_hydramp(f, o, "hydramp_mic"),
    "apex": run_apex,
    "toxinpred": run_toxinpred,
}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fasta", required=True, type=Path, help="Input FASTA file")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "results" / "classifiers",
        help="Directory for output CSV files",
    )
    parser.add_argument(
        "--classifiers",
        nargs="+",
        choices=ALL_CLASSIFIERS,
        default=list(ALL_CLASSIFIERS),
        help="Which classifiers to run (default: all)",
    )
    args = parser.parse_args()

    fasta = args.fasta.resolve()
    if not fasta.is_file():
        print(f"FASTA not found: {fasta}", file=sys.stderr)
        return 1

    out_dir = args.output_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    summary = {}
    failures = {}

    for name in args.classifiers:
        print(f"\n=== Running {name} ===")
        try:
            out_csv = RUNNERS[name](fasta, out_dir)
            summary[name] = str(out_csv)
            print(f"OK: {out_csv}")
        except subprocess.CalledProcessError as exc:
            failures[name] = f"exit code {exc.returncode}"
            print(f"FAILED: {name} ({failures[name]})", file=sys.stderr)
        except Exception as exc:
            failures[name] = str(exc)
            print(f"FAILED: {name} ({exc})", file=sys.stderr)

    manifest = out_dir / f"{fasta.stem}_classifier_manifest.json"
    manifest.write_text(json.dumps({"outputs": summary, "failures": failures}, indent=2) + "\n")
    print(f"\nManifest: {manifest}")

    if failures:
        print(f"\n{len(failures)} classifier(s) failed.", file=sys.stderr)
        return 1
    print(f"\nAll {len(summary)} classifiers completed successfully.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

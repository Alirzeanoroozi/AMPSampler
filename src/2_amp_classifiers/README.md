# Classifier environments and batch runner

One conda environment per classifier, plus a single entry point to score any FASTA.

## Setup (first time)

```bash
cd classifiers

# 1) Download model weights and sync HydrAMP/APEX sources (~500 MB first run)
./download_models.sh

# 2) Create conda envs (10–20 min)
./setup_envs.sh
```

Environments created:

| Env | Classifier |
|-----|------------|
| `ampsampler-ascan2` | AMP Scanner v2 |
| `ampsampler-macrel` | Macrel |
| `ampsampler-hydramp` | HydrAMP AMP + MIC classifiers |
| `ampsampler-apex` | APEX (mean MIC across pathogens) |
| `ampsampler-toxinpred` | ToxinPred3 |

## Run on a FASTA

```bash
./run_classifiers.sh --fasta ../fastas/boltzgen_NDM5_binders.fasta
```

Optional flags:

```bash
./run_classifiers.sh \
  --fasta test_sample.fasta \
  --output-dir ../results/my_run \
  --classifiers macrel apex ampscanner
```

Outputs (one CSV per classifier) go to `results/classifiers/` by default, plus a `{basename}_classifier_manifest.json` summary.

## Run via Slurm

From the repo root (after setup):

```bash
sbatch slurm/run_classifiers.slurm
```

Override input/output directories at submit time:

```bash
FASTAS_DIR=fastas OUTPUT_DIR=results/classifiers sbatch slurm/run_classifiers.slurm
```

This runs all classifiers on every `.fasta` file in `fastas/` (default).

## Test

A 10-sequence sample FASTA is included. After setup:

```bash
./run_classifiers.sh --fasta test_sample.fasta --output-dir ../results/classifiers_test
```

All six outputs should be written with no failures in the manifest.

## Merge results

Combine KPC3 + NDM5 predictions into one table (sequence, target, all classifier columns):

```bash
python merge_classifier_results.py
# -> results/classifiers/merged_classifier_results.csv
```

HydrAMP leaves blank values for sequences longer than 25 aa (see script docstring).

## Binder selection analysis

Merge classifiers with BoltzGen + Boltz2 folding scores, rank binders, and plot:

```bash
conda run -n ampsampler-hydramp python analyze_binder_selection.py
```

Outputs:
- `results/classifiers/merged_binder_selection.csv` — full table (200 rows)
- `results/classifiers/selected_binders_top20.csv` — top 20 per target after safety gates
- `results/classifiers/plots/*.png` — 12 selection plots

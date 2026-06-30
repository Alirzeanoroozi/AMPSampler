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
./run_classifiers.sh --fasta ../fastas/amp_peptide.fasta
```

Optional flags:

```bash
./run_classifiers.sh \
  --fasta test_sample.fasta \
  --output-dir ../results/my_run \
  --classifiers macrel apex ampscanner
```

Outputs (one CSV per classifier) go to `results/classifiers/` by default, plus a `{basename}_classifier_manifest.json` summary.

## Test

A 10-sequence sample FASTA is included. After setup:

```bash
./run_classifiers.sh --fasta test_sample.fasta --output-dir ../results/classifiers_test
```

All six outputs should be written with no failures in the manifest.

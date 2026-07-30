#!/bin/bash
# Chain Stages 3->5 for one target after BoltzGen (Stage 2) has produced designs.
# CPU-only, fast, no internet needed. Run on a login node or a small CPU job.
#
# Usage:
#   PY=python bash pipeline/run_downstream.sh <TARGET> <COMPLEX_DIR> <METRICS_CSV> [N_PANEL]
# where:
#   TARGET      = NDM5 | KPC3
#   COMPLEX_DIR = BoltzGen final designs dir containing the binder-target .cif complexes
#   METRICS_CSV = final_designs_metrics_<budget>.csv (has id + designed_sequence columns)
#   N_PANEL     = number of designs to select for the wet lab (default 32)
set -euo pipefail

T="${1:?usage: run_downstream.sh TARGET COMPLEX_DIR METRICS_CSV [N_PANEL]}"
COMPLEX_DIR="${2:?path to BoltzGen final designs dir (.cif complexes)}"
METRICS_CSV="${3:?path to final_designs_metrics_<budget>.csv}"
N_PANEL="${4:-32}"
PY="${PY:-python}"

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
SCORES="results/scores_${T}"
mkdir -p results/stage2_designs "$SCORES"

echo "== Stage 2->3: designs -> FASTA (design_id == .cif stem) =="
$PY pipeline/stage2_generation/designs_to_fasta.py --csv "$METRICS_CSV" --method boltzgen --target "$T"
FASTA="results/stage2_designs/boltzgen_${T}.fasta"

echo "== Stage 3: BoltzGen own metrics -> joinable scores (iPTM, refold etc.) =="
$PY pipeline/stage3_validation/boltzgen_metrics_to_scores.py \
    --metrics "$METRICS_CSV" --out "$SCORES/boltzgen_metrics.csv"

echo "== Stage 3: active-site overlap (THE inhibitor gate) =="
$PY pipeline/stage3_validation/active_site_overlap.py --target "$T" --complexes "$COMPLEX_DIR" \
    --out "$SCORES/active_site_overlap.csv"

echo "== Stage 4: developability / delivery proxy =="
$PY pipeline/stage4_filtering/developability.py --fasta "$FASTA" --out "$SCORES/developability.csv"

# Optional: drop boltz2_<T>.csv / af_<T>.csv / ddg_<T>.csv / toxinpred_<T>.csv into $SCORES
# before this point and they will be joined automatically.

echo "== Stage 3: build manifest (join all scores on design_id) =="
$PY pipeline/stage3_validation/build_manifest.py --target "$T" --designs "$FASTA" --scores "$SCORES"

echo "== Stage 4: gates + ranking =="
$PY pipeline/stage4_filtering/apply_filters.py --manifest "results/manifest_${T}.csv"

echo "== Stage 5: select diverse panel + alanine controls =="
$PY pipeline/stage5_selection/select_candidates.py --target "$T" --n "$N_PANEL"
$PY pipeline/stage5_selection/alanine_scan.py --target "$T" \
    --selected "results/selected_${T}.fasta" --complexes "$COMPLEX_DIR"

echo "== Quality summary (read this before shipping to the wet lab) =="
$PY pipeline/stage4_filtering/summarize_manifest.py \
    --manifest "results/manifest_${T}.csv" \
    --out "results/quality_report_${T}.md"

echo
echo "DONE [$T]. Deliverables:"
echo "  results/manifest_${T}.csv          (all designs + all scores)"
echo "  results/filtered_${T}.csv          (passed gates, ranked)"
echo "  results/selected_${T}.fasta        (wet-lab panel)"
echo "  results/alanine_controls_${T}.fasta(specificity controls)"
echo "  results/quality_report_${T}.md     (decision-grade summary)"
echo "  -> hand the panel + controls + docs/wetlab_assays.md to the wet lab."

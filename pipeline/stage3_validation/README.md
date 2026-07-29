# Stage 3 — Independent in-silico validation

Stage 2 generators score their own designs and are biased toward their own distribution
(this is why the old repo's Boltz2 iPTMs looked implausibly good — Boltz2 was grading
BoltzGen). Stage 3 adds **independent** evidence and the **active-site requirement**, then
joins everything into one traceable table.

| Step | Script | Runs here? | Question it answers |
|------|--------|-----------|---------------------|
| Active-site overlap | `active_site_overlap.py` | ✅ yes | Does the binder cover the catalytic epitope (not just bind somewhere)? |
| Orthogonal refold | `run_orthogonal_fold.sh` | ❌ GPU | Does an *independent* predictor agree on the bound pose? |
| Boltz-2 affinity | `parse_boltz2.py` | ✅ parses | One affinity/iptm ranking signal (not the only one) |
| Interface ΔΔG / SASA | `run_rosetta_ddg.sh` | ❌ Rosetta | Is the interface energetically real? |
| **Manifest** | `build_manifest.py` | ✅ yes | One row per design, all scores joined on `design_id` |

## 1. Active-site overlap (executable here)
```bash
python active_site_overlap.py --target NDM5 \
  --complexes results/stage2_designs/boltzgen_NDM5/ --cutoff 4.5
```
Outputs `epitope_recall`, `interface_precision`, `n_catalytic_contacts`, `catalytic_ok`.
**A design that does not contact the catalytic core (`catalytic_ok=False`) is not an
inhibitor candidate, no matter how good its iPTM.** Numbering is resolved by alignment, so
it works for BoltzGen (ordinal-numbered) and BindCraft/RFdiffusion (author-numbered) alike.

## 2. Orthogonal refold (needs GPU)
Refold every binder–target complex with a model *different from the generator*:
- **ColabFold / AF2-multimer**: keep `pae_interaction < 10` and `iptm > 0.6`.
- **AF3** (if available) or **Boltz-2** (`boltz predict --use_msa_server`, enable the
  affinity head). Use `<design_id>` as the prediction name so scores join the manifest.
Then `python parse_boltz2.py --boltz_out <dir> --out results/stage3_validation/boltz2_NDM5.csv`.
**Keep designs where the independent model reproduces the designed pose** (low binder RMSD
to the design, high iPTM). This is what kills the self-consistency bias.

## 3. Interface ΔΔG / shape complementarity / buried SASA (needs Rosetta)
`run_rosetta_ddg.sh` wraps Rosetta InterfaceAnalyzer (`dG_separated`, `dSASA_int`, `sc`).
Optionally a short OpenMM/GROMACS MD for complex stability. Keep favorable `dG_separated`
and buried SASA in the typical mini-binder range.

## 4. Build the manifest (executable here)
```bash
python build_manifest.py --target NDM5 \
  --designs results/stage2_designs \
  --scores  results/stage3_validation
```
→ `results/manifest_NDM5.csv`, one row per design with sequence + method + every score,
keyed by `design_id`. Stage 4 filters this; Stage 5 selects from it.

# AMPBinderDesign — peptide inhibitors for NDM-5 and KPC-3

AMPBinderDesign is a sequence-driven pipeline for designing short peptide binders against antimicrobial-resistance proteins.
The first targets are the carbapenemases NDM-5 and KPC-3, which contribute to carbapenem resistance in Klebsiella pneumoniae.

>NDM-5
MELPNIMHPVAKLSTALAAALMLSGCMPGEIRPTIGQQMETGDQRFGDLVFRQLAPNVWQHTSYLDMPGFGAVASNGLIVRDGGRVLLVDTAWTDDQTAQILNWIKQEINLPVALAVVTHAHQDKMGGMDALHAAGIATYANALSNQLAPQEGLVAAQHSLTFAANGWVEPATAPNFGPLKVFYPGPGHTSDNITVGIDGTDIAFGGCLIKDSKAKSLGNLGDADTEHYAASARAFGAAFPKASMIVMSHSAPDSRAAITHTARMADKLR

Close to 4EYL

>KPC-3
MSLYRRLVLLSCLSWPLAGFSATALTNLVAEPFAKLEQDFGGSIGVYAMDTGSGATVSYRAEERFPLCSSFKGFLAAAVLARSQQQAGLLDTPIRYGKNALVPWSPISEKYLTTGMTVAELSAAAVQYSDNAAANLLLKELGGPAGLTAFMRSIGDTTFRLDRWELELNSAIPGDARDTSSPRAVTESLQKLTLGSALAAPQRQQFVDWLKGNTTGNHRIRAAVPADWAVGDKTGTCGVYGTANDYAVVWPTGRAPIVLAVYTRAPNKDDKYSEAVIAAAARLALEGLGVNGQ

Close to 3DW0

## Pipeline

1. run Boltzgen, Pepmlm on each sequence as target and generate binders.
2. run AMPScanner, Macrel, ToxinPred, Apex and Hydramp.
3. Since we have two targets we want to find a peptide that binds to both.
4. filter 50 best.

## Sequence-Based Generation Methods

- Boltzgen
- PepMLM

## Filtering Metrics other that Boltz

1. APEX
2. HydrAMP
3. Macrel
4. AMPScanner

5. Fitness
6. ToxinPred
7. Perpelixity

# Stage 5 — Wet-lab panel selection + specificity controls

| Step | Script | Runs here? |
|------|--------|-----------|
| Diverse panel selection (24–48/target) | `select_candidates.py` | ✅ |

```bash
python select_candidates.py --target NDM5 --n 32 --max-identity 0.8
```

Selection walks designs best-rank-first and keeps one only if <`max-identity` to those
already kept (removes near-duplicate motifs).

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

# Stage 4 — Inhibitor-binder filtering (not AMP filtering)

What changed from the old pipeline:
- **Dropped** "is it an AMP?" (AMPScanner/Macrel-AMP) and **similarity to an AMP database**
  as ranking objectives. The designs are folded enzyme inhibitors; AMP-likeness is the wrong
  target and AMP-DB similarity penalised novelty.
- **Repurposed** the cationic/amphipathic signals as a **periplasmic-delivery proxy**
  (`delivery_proxy` in `developability.py`): NDM-5/KPC-3 are periplasmic, so a binder must
  cross the outer membrane. This informs delivery, not efficacy, and is a soft ranking term.
- **Kept as safety gates**: ToxinPred (toxicity) and Macrel's **Hemo** head (hemolysis) —
  note we now use `hemo_prob`, not `amp_prob`.
- **Added** developability + synthesis liabilities + an aggregation proxy.

| Step | Script | Runs here? |
|------|--------|-----------|
| Developability / delivery proxy | `developability.py` | ✅ |
| Safety (ToxinPred + hemolysis) | `run_safety.sh` | ❌ external installs |
| Gates + ranking | `apply_filters.py` (+ `filters.json`) | ✅ |

```bash
python developability.py --manifest results/manifest_NDM5.csv --out results/stage4_filtering/developability_NDM5.csv
bash   run_safety.sh results/stage2_designs/boltzgen_NDM5.fasta NDM5
# re-build manifest so the new scores join, then:
python apply_filters.py --manifest results/manifest_NDM5.csv
```

`apply_filters.py` applies each gate **only if its column is present**, so it runs at any
stage of completion (it reports which gates were active). Edit thresholds/weights in
`filters.json`. The primary gate is `catalytic_ok` (from Stage 3): a design that does not
contact the catalytic core is not an inhibitor candidate.

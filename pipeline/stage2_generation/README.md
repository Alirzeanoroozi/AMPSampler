# Stage 2 — Epitope-conditioned binder generation

Generate peptide binders **aimed at the Stage-1 active-site epitope** with ≥2 independent
methods, so candidates aren't trusting a single model's biases. All three consume the same
epitope (`targets/<T>/epitope/boltzgen_hotspots.json`) and the cleaned target
(`targets/<T>/structures/<T>_target.pdb`, NDM keeps its catalytic Zn).

> **These tools need a GPU + model weights and are NOT runnable in this environment**
> (no `torch`, no CUDA here). The scripts below are correct and ready to run on a GPU box.
> The config *generators* (`make_*.py`) DO run here and have produced the input files.

| Method | Generator (runs here) | Config produced | Hotspot numbering |
|--------|----------------------|-----------------|-------------------|
| **BoltzGen** | `boltzgen/make_boltzgen_config.py` | `boltzgen/<T>_design.yaml` | **chain ordinal** (1-based) |
| **BindCraft** | `bindcraft/make_bindcraft_settings.py` | `bindcraft/<T>_bindcraft.json` | author/PDB |
| **RFdiffusion** | `rfdiffusion/make_rfdiffusion_cmd.py` | `rfdiffusion/run_rfdiffusion_<T>.sh` | author/PDB (`A120`) |

The numbering differs by tool on purpose — BoltzGen re-indexes each chain from 1, the others
use author/PDB numbers. `boltzgen/<T>_hotspot_map.csv` records the ordinal↔author mapping.

## Run (on a GPU machine)

```bash
# BoltzGen  (pip install boltzgen; downloads weights on first run)
cd boltzgen && bash run_boltzgen.sh

# BindCraft (clone BindCraft; conda env per its README)
cd bindcraft && bindcraft --settings NDM5_bindcraft.json
                bindcraft --settings KPC3_bindcraft.json

# RFdiffusion + ProteinMPNN + AF2 (set $RFDIFFUSION, $PROTEINMPNN)
cd rfdiffusion && bash run_rfdiffusion_NDM5.sh && bash run_rfdiffusion_KPC3.sh
```

Binder length is set to **12–45 aa** in every config (active-site occluder range; edit in the
generators). Suggested scale: ~5k intermediate / ~200 final per target per method.

## Output → Stage 3
Each method yields ranked designs + a metrics CSV. Collect the designed binder sequences into
`results/stage2_designs/<method>_<T>.fasta` and the binder–target complexes (`.cif`/`.pdb`)
into `results/stage2_designs/<method>_<T>/`, then run Stage 3. Keep the **design→sequence→score**
identity intact (the old pipeline lost it — see Stage 3's manifest).

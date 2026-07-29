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

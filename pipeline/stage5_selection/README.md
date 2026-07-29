# Stage 5 — Wet-lab panel selection + specificity controls

| Step | Script | Runs here? |
|------|--------|-----------|
| Diverse panel selection (24–48/target) | `select_candidates.py` | ✅ |
| Alanine-scan specificity controls | `alanine_scan.py` | ✅ |
| Assay request / per-candidate dossier | `../../docs/wetlab_assays.md` | doc |

```bash
python select_candidates.py --target NDM5 --n 32 --max-identity 0.8
# with complexes -> focused interface-Ala controls; without -> full single-Ala scan
python alanine_scan.py --target NDM5 --selected ../../results/selected_NDM5.fasta \
       --complexes ../../results/stage2_designs/boltzgen_NDM5
```

Selection walks designs best-rank-first and keeps one only if <`max-identity` to those
already kept (removes near-duplicate motifs). The alanine controls mirror the BoltzGen
GyrA validation: mutating the designed interface residues to Ala should abolish binding
and inhibition for an on-target binder. Hand `results/selected_<T>.fasta` +
`results/alanine_controls_<T>.fasta` + the dossier to the wet lab; assays are in
`docs/wetlab_assays.md`.

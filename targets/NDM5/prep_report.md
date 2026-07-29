# Stage 0 target prep - NDM5 (New Delhi metallo-beta-lactamase 5)

- Class: Ambler class B1 metallo-beta-lactamase (di-zinc)
- Precursor length: 270 aa
- Predicted signal peptide (literature): residues 1-28. NDM-1 lipoprotein signal peptide ~1-28 (lipidation near Cys26); soluble crystallized constructs typically start ~residue 36-43.
- Design scaffold: 5YPM chain A (seq identity to precursor 0.991)
- Folded/design domain (crystallographically resolved): precursor 43-270 (228 aa)
- Catalytic metals retained in target structure: ['ZN'] (2 ion(s))
- Mutations vs template (already in target sequence; apply to structure with apply_mutations.py if needed): V88L, M154L

## Artifacts
- `design_domain.fasta` - design-domain sequence
- `numbering_map.csv` - precursor <-> scaffold residue map
- `structures/NDM5_target.pdb` - cleaned design target (protein + metals, inhibitor/waters removed)

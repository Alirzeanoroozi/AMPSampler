# Stage 0 target prep - KPC3 (Klebsiella pneumoniae carbapenemase 3)

- Class: Ambler class A serine carbapenemase
- Precursor length: 293 aa
- Predicted signal peptide (literature): residues 1-24. KPC-2 Sec signal peptide ~1-24; mature enzyme 25-293.
- Design scaffold: 3DW0 chain A (seq identity to precursor 0.996)
- Folded/design domain (crystallographically resolved): precursor 30-293 (264 aa)
- Catalytic metals retained in target structure: none (0 ion(s))
- Mutations vs template (already in target sequence; apply to structure with apply_mutations.py if needed): H274Y

## Artifacts
- `design_domain.fasta` - design-domain sequence
- `numbering_map.csv` - precursor <-> scaffold residue map
- `structures/KPC3_target.pdb` - cleaned design target (protein + metals, inhibitor/waters removed)

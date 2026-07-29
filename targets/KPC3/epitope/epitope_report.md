# Stage 1 epitope report - KPC3 (Klebsiella pneumoniae carbapenemase 3)

- Class: Ambler class A serine carbapenemase
- Precursor length: 293 aa
- First crystallographically-ordered residue: precursor pos 30 (biological signal-peptide cleavage annotated in Stage 0)
- Design scaffold: 3DW0 chain A (identity to precursor 0.996)
- Mutations vs template (already in target sequence): H274Y
- Epitope cutoff: 5.0 A around bound inhibitor (+2.8 A around catalytic metal)
- Inhibitor complexes used: 4ZBE, 6D16

## Validation: PASS
active-site residues in epitope: SER=2, LYS=True, GLU=True (expect catalytic Ser + Lys/Glu)

## Source complexes

| PDB | inhibitor | chain | identity | #epitope | title |
|-----|-----------|-------|----------|----------|-------|
| 4ZBE | NXL/A302 | A | 0.996 | 16 | CRYSTAL STRUCTURE OF KPC-2 BETA-LACTAMASE COMPLEXED WITH AVI |
| 6D16 | FUJ/A301 | A | 0.989 | 13 | CRYSTAL STRUCTURE OF KPC-2 COMPLEXED WITH COMPOUND 2 |

## Epitope / hotspot residues (18 total)

Scaffold numbering is what Stage 2 (BoltzGen) consumes.

| precursor | mature | scaffold | residue | min_dist (A) | metal-ligand | sources |
|-----------|--------|----------|---------|--------------|--------------|---------|
| 68 | 39 | 69 | CYS | 3.62 | False | 4ZBE |
| 69 | 40 | 70 | SER | 1.5 | False | 4ZBE|6D16 |
| 72 | 43 | 73 | LYS | 3.59 | False | 4ZBE |
| 103 | 74 | 104 | PRO | 3.8 | False | 6D16 |
| 104 | 75 | 105 | TRP | 3.25 | False | 4ZBE|6D16 |
| 129 | 100 | 130 | SER | 2.62 | False | 4ZBE|6D16 |
| 131 | 102 | 132 | ASN | 3.15 | False | 4ZBE|6D16 |
| 165 | 136 | 166 | GLU | 3.49 | False | 4ZBE |
| 166 | 137 | 167 | LEU | 3.49 | False | 4ZBE|6D16 |
| 169 | 140 | 170 | ASN | 3.42 | False | 4ZBE |
| 215 | 186 | 216 | THR | 3.3 | False | 4ZBE|6D16 |
| 219 | 190 | 220 | ARG | 3.77 | False | 4ZBE|6D16 |
| 233 | 204 | 234 | LYS | 3.31 | False | 4ZBE|6D16 |
| 234 | 205 | 235 | THR | 2.78 | False | 4ZBE|6D16 |
| 235 | 206 | 236 | GLY | 3.47 | False | 4ZBE|6D16 |
| 236 | 207 | 237 | THR | 2.42 | False | 4ZBE|6D16 |
| 237 | 208 | 238 | CYS | 4.14 | False | 4ZBE |
| 272 | 243 | 274 | HIS | 4.59 | False | 6D16 |

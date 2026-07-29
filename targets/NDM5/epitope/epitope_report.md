# Stage 1 epitope report - NDM5 (New Delhi metallo-beta-lactamase 5)

- Class: Ambler class B1 metallo-beta-lactamase (di-zinc)
- Precursor length: 270 aa
- First crystallographically-ordered residue: precursor pos 43 (biological signal-peptide cleavage annotated in Stage 0)
- Design scaffold: 5YPM chain A (identity to precursor 0.991)
- Mutations vs template (already in target sequence): V88L, M154L
- Epitope cutoff: 5.0 A around bound inhibitor (+2.8 A around catalytic metal)
- Inhibitor complexes used: 5YPM, 4EYL, 4EXS, 4RL2

## Validation: PASS
distinct metal-coordinating residues in epitope: 4xHIS, 1xASP, 1xCYS (expect >=3 HIS, >=1 ASP, >=1 CYS for the B1 di-Zn site)

## Source complexes

| PDB | inhibitor | chain | identity | #epitope | title |
|-----|-----------|-------|----------|----------|-------|
| 5YPM | 8YL/A303 | A | 0.991 | 19 | CRYSTAL STRUCTURE OF NDM-1 BOUND TO HYDROLYZED MEROPENEM REP |
| 4EYL | 0RV/A301 | A | 0.991 | 15 | CRYSTAL STRUCTURE OF NDM-1 BOUND TO HYDROLYZED MEROPENEM |
| 4EXS | X8Z/B301 | B | 0.991 | 13 | CRYSTAL STRUCTURE OF NDM-1 BOUND TO L-CAPTOPRIL |
| 4RL2 | 3S3/A303 | A | 0.992 | 17 | STRUCTURAL AND MECHANISTIC INSIGHTS INTO NDM-1 CATALYZED HYD |

## Epitope / hotspot residues (20 total)

Scaffold numbering is what Stage 2 (BoltzGen) consumes.

| precursor | mature | scaffold | residue | min_dist (A) | metal-ligand | sources |
|-----------|--------|----------|---------|--------------|--------------|---------|
| 65 | 23 | 65 | LEU | 4.0 | False | 4RL2 |
| 67 | 25 | 67 | MET | 3.84 | False | 4RL2|5YPM |
| 73 | 31 | 73 | VAL | 3.73 | False | 4EXS|4EYL|4RL2|5YPM |
| 93 | 51 | 93 | TRP | 3.54 | False | 4EXS|4EYL|4RL2|5YPM |
| 120 | 78 | 120 | HIS | 2.05 | True | 4EXS|4EYL|4RL2|5YPM |
| 122 | 80 | 122 | HIS | 1.97 | True | 4EXS|4EYL|4RL2|5YPM |
| 123 | 81 | 123 | GLN | 2.9 | False | 4EYL|4RL2|5YPM |
| 124 | 82 | 124 | ASP | 2.01 | True | 4EXS|4EYL|4RL2|5YPM |
| 125 | 83 | 125 | LYS | 4.73 | False | 4EXS|5YPM |
| 189 | 147 | 189 | HIS | 2.02 | True | 4EXS|4EYL|4RL2|5YPM |
| 190 | 148 | 190 | THR | 4.27 | False | 4EXS|4EYL|4RL2|5YPM |
| 208 | 166 | 208 | CYS | 2.28 | True | 4EXS|4EYL|4RL2|5YPM |
| 211 | 169 | 211 | LYS | 2.94 | False | 4EYL|4RL2|5YPM |
| 216 | 174 | 216 | LYS | 3.98 | False | 5YPM |
| 217 | 175 | 217 | SER | 4.64 | False | 5YPM |
| 218 | 176 | 218 | LEU | 4.0 | False | 4EYL|4RL2|5YPM |
| 219 | 177 | 219 | GLY | 3.3 | False | 4EXS|4EYL|4RL2|5YPM |
| 220 | 178 | 220 | ASN | 2.53 | False | 4EXS|4EYL|4RL2|5YPM |
| 249 | 207 | 249 | SER | 4.18 | False | 4EXS|4EYL|4RL2|5YPM |
| 250 | 208 | 250 | HIS | 1.99 | True | 4EXS|4EYL|4RL2|5YPM |

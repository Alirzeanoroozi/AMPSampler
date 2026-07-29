# AMPBinderDesign — peptide inhibitors for NDM-5 and KPC-3

AMPBinderDesign designs short peptide binders that sit on the **catalytic active site** of two
carbapenem-resistance β-lactamases — **NDM-5** (B1 di-zinc metallo-β-lactamase) and
**KPC-3** (class A serine carbapenemase) — so that meropenem and related carbapenems
work against resistant *K. pneumoniae* again. These are **not** antimicrobial peptides
(membrane lytics). They are folded mini-binders aimed at the enzyme's catalytic
machinery, designed to occlude it.

The conceptual template is the BoltzGen GyrA case study (`image.png`): specify the
functional site as the binding epitope, design peptides to occlude it, validate by
binding + functional inhibition + an alanine specificity control.

## Targets

| | NDM-5 | KPC-3 |
|---|---|---|
| Class | Ambler B1, **metallo**-β-lactamase (di-Zn²⁺) | Ambler A, **serine** carbapenemase |
| Catalytic residues | His120/122/189, Asp124, Cys208, His250 + 2 Zn²⁺ | Ser70, Lys73, Ser130, Glu166 + KTG/oxyanion |
| Inhibitor templates used | 5YPM (meropenem), 4EYL, 4EXS (captopril), 4RL2 | 4ZBE (avibactam), 6D16 |
| Mature domain | precursor 43–270 (228 aa) | precursor 30–293 (264 aa) |

>NDM-5
MELPNIMHPVAKLSTALAAALMLSGCMPGEIRPTIGQQMETGDQRFGDLVFRQLAPNVWQHTSYLDMPGFGAVASNGLIVRDGGRVLLVDTAWTDDQTAQILNWIKQEINLPVALAVVTHAHQDKMGGMDALHAAGIATYANALSNQLAPQEGLVAAQHSLTFAANGWVEPATAPNFGPLKVFYPGPGHTSDNITVGIDGTDIAFGGCLIKDSKAKSLGNLGDADTEHYAASARAFGAAFPKASMIVMSHSAPDSRAAITHTARMADKLR

Close to 4EYL

>KPC-3
MSLYRRLVLLSCLSWPLAGFSATALTNLVAEPFAKLEQDFGGSIGVYAMDTGSGATVSYRAEERFPLCSSFKGFLAAAVLARSQQQAGLLDTPIRYGKNALVPWSPISEKYLTTGMTVAELSAAAVQYSDNAAANLLLKELGGPAGLTAFMRSIGDTTFRLDRWELELNSAIPGDARDTSSPRAVTESLQKLTLGSALAAPQRQQFVDWLKGNTTGNHRIRAAVPADWAVGDKTGTCGVYGTANDYAVVWPTGRAPIVLAVYTRAPNKDDKYSEAVIAAAARLALEGLGVNGQ

Close to 3DW0

## Pipeline

1. run Boltzgen, ProteinHunter on each sequence as target and generate binders.
2. run AMPScanner, Macrel, ToxinPred, Apex and Hydramp.
3. Since we have two targets we want to find a peptide that binds to both.
4. filter 50 best.

## Methods

- Boltzgen
- ESMFold2
- ProteinHunter

## Metrics

1. Amplify
2. APEX
3. HydrAMP
4. Macrel
5. AMPScanner

10. Fitness
12. ToxinPred


AMPBinderDesign panel summary: NDM5  /scratch/anoroozi25/AMPBinderDesign/results/ranked_NDM5.csv  (157 designs)

INHIBITOR GATE  catalytic_ok=True: 95/157  (60.5%)   [STRONG]

STRUCTURE GATE  catalytic_ok and iPTM>=0.5: 42/157  (26.8%)

STRUCTURE GATES  pass: 42/157 (26.8%)

metric                                         n       min       q25    median       q75       max
----------------------------------------------------------------------------------------------------
Boltz-2 iPTM                                 157     0.161      0.36     0.457     0.624     0.855
ipSAE_min (A↔B)                              157         0         0    0.0117    0.0214      0.24
Active-site coverage                         157         0         0      0.25       0.6      0.85
Interface focus on epitope                   157         0         0     0.289      0.52     0.833
# catalytic-core contacts                    157         0         0         1         4         6
pDockQ                                       157    0.0618     0.133      0.18     0.259     0.527
LIS                                          157         0    0.0244    0.0839     0.219     0.575
Boltz-2 complex pLDDT                        157     0.841     0.873     0.887     0.902     0.939
AMPScanner P(AMP)                            157     0.512     0.943     0.999         1         1
Macrel P(AMP)                                157     0.505     0.515     0.525     0.564     0.743
Macrel P(hemolytic)                          157     0.119     0.386     0.515     0.564     0.921
Binder length (aa)                           157        12        26        31        41        49
Net charge at pH 7.4                         157     -2.25     -1.02     -0.99      0.97      5.01
Periplasmic-delivery proxy                   157     0.022     0.124      0.16     0.432     0.918
Aggregation proxy                            157      -0.4      1.77       1.8      2.36      3.06
# synthesis liabilities                      157         0         1         1         1         3
Composite rank score                         157     0.243     0.452      0.53     0.658      0.86

TOP 10 BY rank_score:
  design_id                       rank   iPTM   ipSAE   eRec  #cat    L  charge
  NDM5_1471                      0.860   0.67   0.025   0.75     5   26    +1.0
  NDM5_266_0                     0.857   0.82   0.138   0.75     5   27    -0.0
  NDM5_1416                      0.846   0.84   0.099   0.70     5   22    -1.0
  NDM5_4682                      0.831   0.74   0.075   0.75     5   15    -0.0
  NDM5_790                       0.811   0.77   0.073   0.60     2   29    +1.6
  NDM5_630_0                     0.809   0.75   0.055   0.70     3   18    -0.0
  NDM5_737                       0.806   0.72   0.042   0.85     6   21    -1.4
  NDM5_0637                      0.804   0.76   0.051   0.70     4   26    -1.0
  NDM5_434                       0.803   0.75   0.092   0.75     6   36    -2.2
  NDM5_659                       0.799   0.74   0.042   0.80     6   20    +3.0

WET-LAB READINESS:
  Pool is already AMP-positive (AMPScanner + Macrel). Rank by Boltz-2 iPTM,
  ipSAE_min, and epitope coverage; require catalytic_ok for the shipped panel.
  - catalytic_ok + iPTM>=0.5 yield large enough for n=25: ship that panel.
  - ipSAE_min is typically low here; treat it as a tie-breaker, not a hard cut.
  - Hemolysis is common among AMPs; prefer NonHemo but do not empty the panel.

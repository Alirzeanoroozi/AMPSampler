AMPBinderDesign panel summary: NDM5 selected panel  /scratch/anoroozi25/AMPBinderDesign/results/selected_NDM5.csv  (25 designs)

INHIBITOR GATE  catalytic_ok=True: 25/25  (100.0%)   [STRONG]

STRUCTURE GATE  catalytic_ok and iPTM>=0.5: 25/25  (100.0%)

STRUCTURE GATES  pass: 25/25 (100.0%)

metric                                         n       min       q25    median       q75       max
----------------------------------------------------------------------------------------------------
Boltz-2 iPTM                                  25     0.514     0.632     0.737     0.798     0.855
ipSAE_min (A↔B)                               25    0.0138    0.0329    0.0551    0.0987      0.24
Active-site coverage                          25       0.4       0.6       0.7      0.75      0.85
Interface focus on epitope                    25     0.306       0.5     0.579     0.667     0.833
# catalytic-core contacts                     25         1         4         5         6         6
pDockQ                                        25     0.109     0.148      0.22     0.274     0.527
LIS                                           25     0.205     0.266     0.356     0.499     0.575
Boltz-2 complex pLDDT                         25     0.862     0.888     0.895     0.914     0.939
AMPScanner P(AMP)                             25     0.597     0.912      0.98         1         1
Macrel P(AMP)                                 25     0.505     0.525     0.554     0.584     0.683
Macrel P(hemolytic)                           25     0.188     0.347     0.475     0.624     0.772
Binder length (aa)                            25        12        21        26        32        43
Net charge at pH 7.4                          25     -2.25     -1.02     -0.02      1.01      3.05
Periplasmic-delivery proxy                    25     0.022     0.149     0.221     0.546     0.754
Aggregation proxy                             25      -0.4      1.17      1.61       1.8      3.04
# synthesis liabilities                       25         0         1         1         1         2
Composite rank score                          25     0.715     0.748     0.791     0.806      0.86

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

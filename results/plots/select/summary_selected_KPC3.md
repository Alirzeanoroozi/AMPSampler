AMPBinderDesign panel summary: KPC3 selected panel  /scratch/anoroozi25/AMPBinderDesign/results/selected_KPC3.csv  (25 designs)

INHIBITOR GATE  catalytic_ok=True: 25/25  (100.0%)   [STRONG]

STRUCTURE GATE  catalytic_ok and iPTM>=0.5: 25/25  (100.0%)

STRUCTURE GATES  pass: 25/25 (100.0%)

metric                                         n       min       q25    median       q75       max
----------------------------------------------------------------------------------------------------
Boltz-2 iPTM                                  25     0.615     0.705     0.782     0.805     0.923
ipSAE_min (A↔B)                               25    0.0147    0.0463    0.0596    0.0892     0.302
Active-site coverage                          25     0.611     0.722     0.833     0.889     0.944
Interface focus on epitope                    25     0.341     0.457     0.524     0.577     0.684
# catalytic-core contacts                     25         2         3         3         3         4
pDockQ                                        25     0.176     0.267       0.3     0.346     0.462
LIS                                           25     0.247     0.344     0.395     0.439     0.569
Boltz-2 complex pLDDT                         25     0.846      0.88     0.896     0.911     0.943
AMPScanner P(AMP)                             25      0.52     0.797     0.951     0.984     0.999
Macrel P(AMP)                                 25     0.505     0.525     0.535     0.554     0.594
Macrel P(hemolytic)                           25     0.198     0.317     0.396     0.485     0.812
Binder length (aa)                            25        13        16        22        28        43
Net charge at pH 7.4                          25     -1.02      -0.1     -0.02      0.97      2.01
Periplasmic-delivery proxy                    25     0.018     0.085     0.162     0.258     0.572
Aggregation proxy                             25    -0.343     0.343     0.943      1.29         3
# synthesis liabilities                       25         0         1         1         1         3
Composite rank score                          25     0.713     0.725     0.765     0.792     0.841

TOP 10 BY rank_score:
  design_id                       rank   iPTM   ipSAE   eRec  #cat    L  charge
  KPC3_573_0                     0.841   0.83   0.123   0.83     3   22    +1.0
  KPC3_414_0                     0.826   0.75   0.046   0.89     3   23    -1.0
  KPC3_474                       0.821   0.73   0.052   0.83     3   15    -0.1
  KPC3_661_0                     0.820   0.88   0.186   0.83     3   22    +1.0
  config_KPC3_809                0.814   0.83   0.061   0.72     3   16    -0.0
  KPC3_115                       0.794   0.82   0.156   0.89     4   34    +0.9
  KPC3_862_0                     0.792   0.78   0.054   0.83     3   16    +2.0
  KPC3_697_1                     0.787   0.92   0.302   0.61     3   13    -0.0
  KPC3_662                       0.780   0.68   0.053   0.94     4   22    -1.0
  KPC3_559                       0.778   0.79   0.060   0.89     3   17    +2.0

WET-LAB READINESS:
  Pool is already AMP-positive (AMPScanner + Macrel). Rank by Boltz-2 iPTM,
  ipSAE_min, and epitope coverage; require catalytic_ok for the shipped panel.
  - catalytic_ok + iPTM>=0.5 yield large enough for n=25: ship that panel.
  - ipSAE_min is typically low here; treat it as a tie-breaker, not a hard cut.
  - Hemolysis is common among AMPs; prefer NonHemo but do not empty the panel.

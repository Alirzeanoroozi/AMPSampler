AMPBinderDesign panel summary: KPC3  /scratch/anoroozi25/AMPBinderDesign/results/ranked_KPC3.csv  (210 designs)

INHIBITOR GATE  catalytic_ok=True: 194/210  (92.4%)   [STRONG]

STRUCTURE GATE  catalytic_ok and iPTM>=0.5: 125/210  (59.5%)

STRUCTURE GATES  pass: 122/210 (58.1%)

metric                                         n       min       q25    median       q75       max
----------------------------------------------------------------------------------------------------
Boltz-2 iPTM                                 210     0.256      0.46      0.57     0.676     0.923
ipSAE_min (A↔B)                              210         0    0.0107    0.0134    0.0313     0.302
Active-site coverage                         210         0     0.444     0.667     0.778         1
Interface focus on epitope                   210         0     0.286       0.4     0.485     0.737
# catalytic-core contacts                    210         0         1         2         3         4
pDockQ                                       210    0.0689     0.208     0.262     0.311     0.509
LIS                                          210         0    0.0715     0.174     0.282     0.569
Boltz-2 complex pLDDT                        210     0.831     0.858     0.872      0.89     0.943
AMPScanner P(AMP)                            210     0.502     0.908     0.993     0.999         1
Macrel P(AMP)                                210     0.505     0.525     0.545     0.584     0.752
Macrel P(hemolytic)                          210     0.188     0.356     0.485     0.594     0.822
Binder length (aa)                           210        12        22        29        36        45
Net charge at pH 7.4                         210     -2.14     -1.02     -0.02     -0.02      7.86
Periplasmic-delivery proxy                   210     0.008     0.118     0.163     0.257         1
Aggregation proxy                            210      -0.4     0.729      1.17      1.74         3
# synthesis liabilities                      210         0         1         1         2         5
Composite rank score                         210     0.167     0.443     0.546     0.666     0.841

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

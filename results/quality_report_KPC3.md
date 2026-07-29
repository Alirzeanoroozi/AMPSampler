```
PepGen quality summary: results/manifest_KPC3.csv  (200 designs)

INHIBITOR GATE  catalytic_ok=True: 165/200  (82.5%)   [STRONG]

metric                                         n       min       q25    median       q75       max
----------------------------------------------------------------------------------------------------
Active-site coverage (fraction)              200     0.222     0.556     0.667     0.778     0.944
# catalytic-core residues contacted          200         0         1         2         3         4
Interface focus on active site               200     0.222      0.48     0.542       0.6     0.833
iPTM                                         200     0.231      0.34     0.441     0.575     0.737
Binder length (aa)                           200        12        33        38        41        45
Net charge at pH 7.4                         200     -9.02     -5.02     -3.03     -2.02      1.01
Hydrophobic moment (amphipathicity)          200     0.048     0.524     0.599     0.672     0.948
Periplasmic-delivery proxy [0..1]            200      0.04     0.437       0.5       0.5     0.584
Aggregation proxy (KD window max)            200    -0.114     0.771      1.14      1.71      3.01
# synthesis liabilities                      200         0         0         0         1         2

WET-LAB READINESS HEURISTIC:
  - catalytic_ok yield >= 30%:           strong; ship the top-32 panel + alanine controls.
  - catalytic_ok yield 10..30%:          marginal; ship but expect a lower hit rate in vitro.
  - catalytic_ok yield < 10%:            weak; rerun BoltzGen with a tighter hotspot list,
                                         or filter to catalytic_ok=True before selection.
  - median iPTM >= 0.7 in top-32:        good binder-target interface confidence.
  - median net charge in top-32 +2..+8:  reasonable periplasmic-delivery profile.
```

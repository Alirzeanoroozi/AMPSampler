```
PepGen quality summary: results/manifest_NDM5.csv  (200 designs)

INHIBITOR GATE  catalytic_ok=True: 200/200  (100.0%)   [STRONG]

metric                                         n       min       q25    median       q75       max
----------------------------------------------------------------------------------------------------
Active-site coverage (fraction)              200      0.35      0.65      0.75       0.8       0.9
# catalytic-core residues contacted          200         1         4         5         6         6
Interface focus on active site               200     0.467       0.6      0.65     0.706     0.941
iPTM                                         200     0.315     0.379     0.448     0.518     0.674
Binder length (aa)                           200        12        34        38        42        45
Net charge at pH 7.4                         200     -9.02     -4.02     -2.03     -1.02      2.98
Hydrophobic moment (amphipathicity)          200     0.083     0.497     0.609     0.716      1.04
Periplasmic-delivery proxy [0..1]            200     0.069     0.414       0.5       0.5     0.748
Aggregation proxy (KD window max)            200    -0.657     0.771      1.16      1.54      2.53
# synthesis liabilities                      200         0         0         1         1         3

WET-LAB READINESS HEURISTIC:
  - catalytic_ok yield >= 30%:           strong; ship the top-32 panel + alanine controls.
  - catalytic_ok yield 10..30%:          marginal; ship but expect a lower hit rate in vitro.
  - catalytic_ok yield < 10%:            weak; rerun BoltzGen with a tighter hotspot list,
                                         or filter to catalytic_ok=True before selection.
  - median iPTM >= 0.7 in top-32:        good binder-target interface confidence.
  - median net charge in top-32 +2..+8:  reasonable periplasmic-delivery profile.
```

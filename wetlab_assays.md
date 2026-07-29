# Wet-lab validation plan — NDM-5 / KPC-3 peptide inhibitors

The claim is **"binder restores carbapenem susceptibility by inhibiting the β-lactamase."**
That is three separate measurements — binding, enzyme inhibition, and cellular
re-sensitization — plus a specificity control. This mirrors the BoltzGen GyrA campaign
(`image.png`): designed binders were validated by *binding* **and** *functional inhibition*,
with an *alanine-mutation* control proving the effect was via the designed interface.

Panel: ~24–48 designs/target from `results/selected_<T>.fasta`, each with its alanine
controls from `results/alanine_controls_<T>.fasta`.

## Target biochemistry (drives assay choice)
| | NDM-5 | KPC-3 |
|---|---|---|
| Class | B1 **metallo**-β-lactamase (di-Zn²⁺) | A **serine** carbapenemase |
| Mechanism | Zn-activated hydroxide; no covalent intermediate | Ser70 acyl-enzyme |
| Construct | mature soluble domain **with 2 Zn²⁺** (do not strip metal) | mature soluble domain (Ser70 active) |
| Known inhibitor mode | zinc engagement (captopril), cyclic boronates (taniborbactam) | covalent Ser70 (avibactam), boronates (vaborbactam) |
| Active site we targeted | His120/122/189, Asp124, Cys208, His250 + pocket | Ser70/Lys73/Ser130/Glu166 + KTG/oxyanion |

## 1. Binding (does it bind, and how tightly?)
- **BLI or SPR**: immobilize purified target (or biotinylated peptide); titrate to get K_D, k_on, k_off.
- **ITC** for a label-free K_D + stoichiometry on the strongest hits.
- Expect: selected designs bind in the nM–low-µM range if the in-silico interface is real.

## 2. Enzyme inhibition (does binding block catalysis? — the readout the old pipeline lacked)
- **Chromogenic**: nitrocefin hydrolysis (Δabs 482 nm) ± peptide → IC₅₀, then Kᵢ and mode
  (competitive/uncompetitive) from Michaelis–Menten ± inhibitor.
- **Clinically relevant substrate**: imipenem/meropenem hydrolysis (Δabs ~297–300 nm).
- NDM-5: include physiological Zn²⁺; check inhibition isn't merely nonspecific metal stripping
  (compare ± excess Zn²⁺). KPC-3: a slow/tight covalent vs reversible binder both acceptable.

## 3. Re-sensitization (the actual clinical claim)
- **Checkerboard MIC**: meropenem × peptide on a KPC-3⁺ and an NDM-5⁺ *K. pneumoniae* isolate;
  report MIC shift and **FIC index** (synergy if ≤0.5).
- **Delivery caveat**: both enzymes are periplasmic, so a peptide must cross the outer membrane.
  If enzyme inhibition is strong but cellular rescue is weak, that is a *delivery* problem, not a
  binding failure — test with a permeabilizer (e.g., sub-MIC polymyxin B nonapeptide) or an
  OM-penetrating tag. `delivery_proxy` in the manifest pre-flags likely-poor penetrators.

## 4. Specificity control (proves on-target binding)
- For each active design, test its **alanine-interface mutant** (`*_ifaceAla`) and key single-Ala
  mutants: binding **and** inhibition should drop sharply if the effect is via the designed epitope.

## 5. Counter-screens (safety & selectivity)
- Mammalian cytotoxicity + hemolysis (confirm the in-silico ToxinPred / Macrel-Hemo gates).
- Selectivity: test top hits against a **non-target** β-lactamase (e.g., a class C AmpC or TEM-1)
  — a good active-site binder should prefer its design target.

## Per-candidate dossier (send to wet lab)
From `results/selected_<T>.csv` + `results/manifest_<T>.csv`, per design: sequence, length,
generation method, predicted complex (CIF), **epitope residues contacted** + `catalytic_ok`,
orthogonal-fold iPTM / pAE_int, interface ΔΔG, developability + delivery_proxy, safety flags,
and its paired alanine controls.

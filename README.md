1. We have selected these targets

We could maybe try with NTUH isolate of K. pneumonia (NCBI RefSeq GCF_000009885.1, Assembly ASM988v1)


---
OXA-48, a carbapenemase-resistance gene, are either 6P96 or 4S2P. Antibiotic-bound forms are 6P99, 6P98 and 6P97

https://card.mcmaster.ca/

## Targets
>NDM-5
MELPNIMHPVAKLSTALAAALMLSGCMPGEIRPTIGQQMETGDQRFGDLVFRQLAPNVWQHTSYLDMPGFGAVASNGLIVRDGGRVLLVDTAWTDDQTAQILNWIKQEINLPVALAVVTHAHQDKMGGMDALHAAGIATYANALSNQLAPQEGLVAAQHSLTFAANGWVEPATAPNFGPLKVFYPGPGHTSDNITVGIDGTDIAFGGCLIKDSKAKSLGNLGDADTEHYAASARAFGAAFPKASMIVMSHSAPDSRAAITHTARMADKLR

Close to 4EYL

>KPC-3
MSLYRRLVLLSCLSWPLAGFSATALTNLVAEPFAKLEQDFGGSIGVYAMDTGSGATVSYRAEERFPLCSSFKGFLAAAVLARSQQQAGLLDTPIRYGKNALVPWSPISEKYLTTGMTVAELSAAAVQYSDNAAANLLLKELGGPAGLTAFMRSIGDTTFRLDRWELELNSAIPGDARDTSSPRAVTESLQKLTLGSALAAPQRQQFVDWLKGNTTGNHRIRAAVPADWAVGDKTGTCGVYGTANDYAVVWPTGRAPIVLAVYTRAPNKDDKYSEAVIAAAARLALEGLGVNGQ

Close to 3DW0

2. run boltzgen on each sequence
3. run boltz2 for each peptide and target
4. run AMPScanner, MAcrel, ToxinPred, Similarity test, ModlAMP, ppl
5. filter 50 best

6. Do this for ProteinHunter and DPLM then
7. Add the Hallucitation effect also that 5 models Bg and Ph and Dplm and AMp-Diffusion and my model ASD(AMP Sequecne Diffusion) also generates good AMPS.
put the places for results and add more introduction about tagerts. and also conculsion


# Steps
Look what AMP-diffusion did. Compare the Boltzgen, Proteinhunter, DPLM in each target design and also in hallucination method
Also add the AMP-diffusion method with ASD


1. SoTA protein design approach
You have two targets:
- Use boltzgen, bindcraft, … to design binders with target
- Toxin, Macrel, AMPScanner.
- Calculte some metrics for similarity with AMPs.
- Do the inlab experiment.

2. Without targets we got better results
How these works without target
- Input xxxxxxxx to all three methods and get AMPs.
- Run boltz2, Af2 to get affinity bindings to specific target (NDM).
- Toxin, Macrel, AMPScanner.
- Calculte some metrics for similarity with AMPs.
- Do the inlab experiment.

3. Integrate the conditioning step to guide the generation
- Using Macrel to improve the performance  
- Toxin, Macrel, AMPScanner.
- Calculte some metrics for similarity with AMPs.
- Do the inlab experiment.

# AMPScanenr
used this https://www.dveltri.com/ascan/v2/ascan.html

# Metrics
1. Amplify
2. APEX
3. HydrAMP (MIC)
4. Macrel
5. AMPScanner

6. ppl
7. Diversity
8. Entropy
9. K-Mer Similarity
10. Fitness
11. Uniqness

12. ToxinPred
13. AminoAcid Frequency
14. ModlAMP


python amp_scanner_v2_predict_tf1.py -fasta dit_amp_sequences.fasta -model trained-models/OriginalPaper_081917_FULL_MODEL.h5 -candidates My_AMP_Candidates.fasta -preds dit_amp_sequences.csv
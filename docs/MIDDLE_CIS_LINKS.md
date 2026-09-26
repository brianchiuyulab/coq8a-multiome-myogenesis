# Expanded cis-target exploration

## Design and reproduction

All 34 externally defined middle-opening peaks were retained. Candidate targets were all protein-coding genes with at least one GENCODE v48 transcript TSS within 500 kb of a peak midpoint, irrespective of membership in the initial 221-gene set. There were 593 candidate peak–gene pairs present in each library's RNA matrix. Distance specifies candidates and does not establish regulatory targeting.

The four existing COQ8A/TSS comparisons and their within-library matched nuclei were unchanged. Within each library, binary peak accessibility and log1p(CP10K) RNA were residualized against COQ8A group and log1p RNA/ATAC depths. Correlations required at least 10 open nuclei and 20 RNA-detected nuclei per library. Links with at least three eligible libraries were combined using Fisher-z weights n−5. Two-sided normal p values were adjusted by BH across all testable links within each setting. These are exploratory nucleus-level associations; libraries from the same source line are not independent biological replicates. This custom screen is not Signac LinkPeaks and does not estimate its GC/background-peak null.

```powershell
python scripts/30_middle_cis_links.py --tables results/tables --temporal results/temporal --gtf ../data/public/gencode.v48.annotation.gtf.gz --h5-root ../data/public/GSE208248_processed
python scripts/31_core_atac_summary.py --tables results/tables --out results/temporal
```

## Main comparison: TSS ≥3, COQ8A ≥3 versus exactly 1 UMI

The same 201 matched pairs were used. Of 593 candidates, 38 links were evaluable in at least three libraries. None reached BH q<0.05. Shared links reproduce the previous narrower screen's correlations within 2×10−7.

- LDB3, chr10:86674717–86675570: partial r=0.13081, p=0.01512, q=0.19152, positive in all three eligible libraries. This link statistic is distinct from the peak's high/low accessibility comparison (FC=1.9375, p=0.03154, within-middle-34 q=0.53617).
- COL15A1, chr9:98916421–98917309: COL15A1 RNA was detected in 24 of 402 nuclei, with no library reaching the 20-detected-nucleus eligibility floor. Its RNA link is not evaluable under this rule. Nearby ANKS6 and TGFBR1 did not provide significant alternative links.
- The smallest positive-link p was for POLR1A (r=0.15894, p=0.003075, q=0.07507).

No MYOD1, MYOG, MYF5, MYF6, MEF2A, MEF2C or MEF2D cis candidate occurred within 500 kb of these 34 peaks. MEF2B was a distance candidate. This search does not assess trans regulation or TF binding at distant motifs.

## Sensitivity comparisons

| TSS minimum | COQ8A high versus 1 UMI | Testable links | Positive q<0.05 links |
|---|---|---:|---:|
| 2 | ≥2 | 374 | 6 |
| 2 | ≥3 | 75 | 0 |
| 3 | ≥2 | 329 | 5 |
| 3 | ≥3 | 38 | 0 |

The TSS≥3, COQ8A≥2 comparison identified positive links to MRAS, SRPK3, TNNC1, L1CAM and HRC. These belong to that sensitivity comparison and are not substitutes for the main comparison. Positive peak–RNA correlation also does not imply higher accessibility in the COQ8A-high group.

## Broad core-programme results

`core_programme_ATAC_summary.tsv` and `core_TF_ATAC_summary.tsv` extract the unchanged main comparison. Programme ratios average gene-region accessibility equally before dividing high by low; they are not fragment-count TSS-profile ratios. ATAC Hallmark myogenesis FC=1.00444 (p=0.47047); MRF FC=1.06050 (p=0.15106); MEF2 FC=0.99314 (p=0.83449). None reaches its existing BH threshold. Core gene-region results quantify regions near the named genes, not genome-wide occurrences of their TF motifs.

These broad ATAC findings coexist with the RNA Hallmark programme association (q=0.00990) and the externally defined middle-opening ATAC module (FC=1.22293, ten-module q=0.01785). They measure different, explicitly defined features.

# Results: parallel COQ8A ≥2 and ≥3 versus 1 UMI screens

## Main 201-pair comparison

Both contrasts use within-library depth-matched pairs at ATAC TSS enrichment ≥3 and the same fixed Hallmark/Reactome myogenesis union of 221 measured genes. The broader COQ8A ≥2 versus 1 contrast has 958 pairs and 6,871 of 7,132 candidate peaks passing the pooled prevalence floor. The stronger ≥3 versus 1 contrast has 201 pairs (46/74/58/23 across four libraries) and 3,731 eligible peaks. The biological source-line count is two.

Under ≥3 versus 1, RNA programme scores rise with higher COQ8A: Hallmark Myogenesis +0.03519 log1p(CP10K) per gene (four-programme BH q=0.00990; 4/4 library directions), Reactome Myogenesis +0.02257 (q=0.03158; 3/4), and four MRF loci +0.05899 (q=0.00779; 4/4). The corresponding ATAC programme summaries do not have q<0.05. Among all 221 ATAC gene regions, the minimum q is 0.30246 at ≥2 and 0.14537 at ≥3; among eligible peaks, the minimum q is 0.11899 at ≥2 and 0.58398 at ≥3. Neither complete contrast identifies an FDR-positive ATAC region or peak.

## Candidate ranking and MYOD1 locus

The exploratory ranking combines positive peak–RNA links learned from the more numerous COQ8A ≥2-versus-1 nuclei with ATAC effect direction **separately for each contrast**. MYOD1 ranks first under ≥2 versus 1 (five qualifying links). Under ≥3 versus 1, CKB ranks first (four), and MYOD1 second (three). MYOD1 is followed as a myogenic transcriptional regulator, not because it ranks first under every setting. The rankings use overlapping nuclei and are not independent validation.

For all 19 four-library-common peaks near MYOD1, the ≥2-versus-1 result is **1.0688** (12.801% versus 11.977%; 958 pairs; 221-gene BH q=0.6354). The ≥3-versus-1 result is **1.0934** (13.485% versus 12.333%; 201 pairs; nominal paired p=0.1625; 3/4 positive library directions; 221-gene BH q=0.9210). Both are MYOD1 gene-region results within their respective complete 221-gene screens.

The six MYOD1 peaks with positive RNA links learned at the TSS≥2 link-discovery gate and exploratory MRF PWM score <0.95 give 13.267% versus 10.365% open in the **same main 201 pairs**: +2.902 percentage points, ratio **1.2800**, 4/4 positive directions and nominal paired p=0.02558. The six peaks were chosen using overlapping data, so this p value describes the observed subset and is not a fresh confirmatory test. The 20 region/count/TSS settings are sensitivity checks; no BH q value is calculated across them.

Across the displayed 20-setting grid, **1.2800 is the largest fold ratio**. The smallest nominal p occurs under a different setting: 10 positive-linked peaks at TSS≥2 with COQ8A ≥2 versus 1 (ratio 1.1055, p=0.01248, 1,020 pairs). Thus the setting with the largest estimated effect is not the one with the smallest p value. Neither is an independent confirmation of the peak set selected from the same data.

Using TSS≥3 during link discovery leaves five such peaks. Their open fraction in the same 201 pairs is 13.831% versus 11.144%, ratio **1.2411**, with 3/4 positive library directions. The extra sixth peak is open in 21 high versus 13 low nuclei in those pairs. Under the broader ≥2-versus-1 COQ8A contrast, the six-peak ratio is **1.1097** in 958 pairs, with 3/4 positive directions. Excluding pairs that touch the highest 2.5% RNA Scrublet-score tail lowers the main six-peak ratio to **1.192** in 187 pairs (nominal p≈0.12). The complete settings and pair-level records are in `results/tables/`.

## TSS-aligned fragment context

Figure 1 uses the **same 201 main-contrast pairs** at all 221 gene TSSs. Within ±5 kb, high nuclei have 57,556 Tn5 cuts and low nuclei 56,638, a descriptive ratio of **1.0162**. This nearly overlapping mean profile does not select the MYOD1 peaks. It is a different readout from the six-peak open-nucleus ratio of 1.280.

## Interpretation

The data support an RNA myogenesis association and an exploratory, localized MYOD1 accessibility signal under the stronger COQ8A contrast. **Both** complete ATAC screens are negative at their stated FDR thresholds, and the local effect weakens under the RNA doublet-score challenge. These observations support a target hypothesis for independent testing; they do not establish COQ8A-driven chromatin opening. Figure 2 shows both full 221-gene screens and explicitly contrasts candidate rankings; Figure 3 shows the locus and sensitivity settings.

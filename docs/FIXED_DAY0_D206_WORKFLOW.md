# Fixed Day-0 COQ8A high/undetected analysis

## Question and fixed design

Test whether COQ8A-high nuclei differ in local accessibility near externally defined myoblast differentiation/fusion genes, then assess possible local RNA targets. The comparison is D206 throughout: GSE208248 undifferentiated cultures only, TSS enrichment >=3, COQ8A >=2 versus 0 raw UMI, within-library RNA/ATAC depth matching, caliper 0.30. No marker-positive population restriction or state matching is imposed in the primary contrast. Zero means undetected RNA, not proven absence of expression. These are cultures derived from human muscle stem cells, not newly annotated pure quiescent MuSC nuclei.

The setting was fixed following sensitivity exploration. External gene membership does not depend on current effect estimates; parameter selection was not blinded or preregistered. Current results never substitute a better setting for an individual gene. GSE240061, Day 7, private C2C12 RNA and prior temporal/motif screens do not determine eligibility here.

## 1. External scope and peak construction

Frozen Hallmark/Reactome myogenesis memberships define 221 measured genes. Intersect with the union of GO myoblast differentiation and GO myoblast fusion: 17 + 12 - 4 shared = 25 genes. These are process annotations, not a list exclusively of positive regulators. All members and provenance are in `results/differentiation_fusion25/genes.tsv` and `reference/myogenesis_221_gene_sources.tsv`.

Start from deposited GRCh38 author-called peaks. Anchor on GSM6339597 and match the other three original catalogs one-to-one with >=50% reciprocal overlap. Retain canonical chromosomes, 200–2,000-bp intervals, no ENCODE hg38 v2 blacklist overlap and peak midpoints within 100 kb of an annotated protein-coding transcript TSS. The broad scope has 5,097 common peaks; restricting neighborhood membership to the 25 genes gives 565 distinct peaks. A peak may belong to more than one neighborhood. Catalog construction uses all four original libraries; effect tests below use only Day 0. No FASTQ alignment or new peak calling is claimed.

## 2. Joint QC and pairing

Use GSM6339597 and GSM6339601, two donor-derived cell lines. Both modalities are paired by library and barcode. Apply RNA >=500 UMI; ATAC >=500 detected peaks; exclude the top 5% of either depth within library after depth floors; FRiP >=0.25; nucleosome signal <4; blacklist fragment fraction <0.05; TSS enrichment >=3. TSS scores inherit the archived GENCODE v48 central +/-500-bp/flank-density implementation, including its library-flank fallback. Comprehensive doublet removal was not part of this baseline. No Harmony correction, imputation, cell-type relabeling or trajectory fitting is used for count tests.

Match without replacement within library, requiring absolute differences <=0.30 in both log1p RNA UMI and log1p detected-peak depth, using the archived 100-neighbor greedy matching algorithm. There are 226 pairs from source 1 and 301 from source 2: 527 pairs / 1,054 nuclei. Biological-source N is two. The barcode list is unchanged for all peaks and RNA high/low comparisons.

## 3. Complete ATAC screen and nomination

For all 565 peaks, a positive author matrix entry defines accessibility. Report the ratio of accessible-nucleus fractions, high/low. Exact two-sided paired-binomial tests use high-only versus low-only discordant pairs; BH correction covers all 565 peaks.

For each of all 25 neighborhoods, calculate the maximum standardized paired-difference statistic across constituent peaks. Opening and closing are separate directional families. Within-pair label permutations preserve dependence between peaks. D206 uses 2,000,000 permutations, seed 20261134; p=(exceedances+1)/(B+1), followed by BH25 within direction. Nominate q<0.05 regions; retain 0.05<=q<0.10 as a separately labelled exploratory tier. Follow each nominated region's maximum-statistic peak, deduplicating shared coordinates. No minimum FC, favored gene, motif or RNA change is an eligibility gate. The region q is not the constituent peak q.

| Direction | Neighborhood label | Peak | Peak FC | Region q25 | Peak q565 |
|---|---|---|---:|---:|---:|
| Opening | MYOD1 | chr11:17653349–17654252 | 1.510 | 0.00796 | 0.0303 |
| Opening | CSRP3 | chr11:19218592–19219518 | 1.717 | 0.0135 | 0.0476 |
| Opening | CACNA1H | chr16:1078248–1079160 | 2.474 | 0.0229 | 0.0434 |
| Opening | CAV3 | chr3:8757835–8758737 | 1.316 | 0.0401 | 0.0907 |
| Closing | BOC | chr3:113117234–113118143 | 0.233 | 0.0302 | 0.0434 |
| Closing | NOTCH1 | chr9:136586858–136587780 | 0.459 | 0.0302 | 0.0434 |

All six top peaks agree in direction across sources. The additional tier consists of NOS1, MYH9 and MYF5/MYF6 neighborhoods; MYF5/MYF6 share a peak. Ten neighborhood labels therefore nominate nine distinct peaks. Their complete counts, p/q and directions are retained. Large FC=11 at NOS1 and MYF5/MYF6 reflects only 11 versus 1 accessible nuclei, not broadly accessible sites.

## 4. Local target models

For every nominated peak, evaluate all uniquely mapped deposited RNA features with a GENCODE transcript TSS within +/-500 kb, including noncoding and out-of-scope genes: 170 peak–RNA pairs / 170 RNA genes across nine peaks. Neighborhood naming is not target assignment. For example, the CACNA1H-neighborhood peak lies near SSTR5/SSTR5-AS1 transcript starts; it is not the previously examined chr16:1311478–1312392 peak.

Extract raw RNA and binary ATAC from the same deposited H5 files. Fit source-specific partial correlations between accessibility and log1p CP10k RNA, controlling log1p RNA depth, log1p detected-peak depth and log1p COQ8A CP10k. Fit a second model adding the archived myogenesis-state score (broad gene set excluding the original 25 scope genes). Combine source estimates by Fisher z weighted by residual degrees of freedom minus one. Require >=10 RNA-positive and >=10 peak-positive observations in each source for combined estimation; retain all unavailable tests explicitly. BH covers all 170 local pairs within each population/model, using p=1 for unavailable hypotheses and displaying their p/q as NA.

Fit both all eligible Day-0 nuclei (7,535 + 9,346 = 16,881) and the exact 1,054 matched nuclei. The former improves target mapping coverage; the latter assesses the selected comparison population. These are exploratory nucleus-level associations, not independent-donor causal tests or proof of enhancer action. No RNA significance gate is applied to ATAC nomination.

Separately compare candidate RNA on the exact 527 pairs: FC uses mean linear CP10k; paired t tests use log1p CP10k differences, followed by BH across all 170 local RNA genes. Per-source RNA means and detection counts are also retained.

### Main interpretations

- **MYOD1:** full-population link r=0.157, state-adjusted r=0.118; both sources positive. In exact matched nuclei, r=0.171 (q=3.95e-6), state-adjusted r=0.118 (q=0.0232). This is the strongest positive link among the assessed local pairs. Same-pair RNA FC=1.107, p=0.174, q=0.558: no significant immediate RNA increase. The locus is a plausible MYOD1-associated regulatory candidate, not a proven driver. Primary COQ8A-associated ATAC enrichment itself is state-sensitive in the archived grid; persistence of the peak–RNA link does not change that fact.
- **CAV3:** full-population r=0.0591, state-adjusted r=0.0612, both positive sources. RNA FC=1.253, p=0.000215, q=0.00405, both source means increase. Matched-only link is weaker (r=0.0339, p=0.274; state-adjusted r=0.0349, p=0.260). Thus it has concordant ATAC/RNA group changes and a weak full-population cis association, not a robust link in every population.
- **CSRP3:** RNA FC=3.082, p=0.00316, q=0.0316; full-population link r=0.0289 becomes approximately zero after state adjustment. Retain its accessibility and RNA associations without claiming a direct cis chain.
- **CACNA1H neighborhood:** the named gene has no convincing positive local link. SOX8 and LMF1 have small positive associations; SOX8 r=0.0555 falls to 0.0223 with state adjustment. They remain alternative candidate targets, not confirmed targets.
- **Closing loci:** retain BOC and NOTCH1. The BOC-neighborhood peak has a negative association with BOC RNA and positive association with CCDC80; closing cannot automatically be translated into lower BOC RNA. NOTCH1 has a small positive link but no significant same-pair RNA change.

## 5. Figures and reading guide

- **Figure 1:** A, external scope and sample flow; B, all 565 peaks, four source-by-group columns, one common 0–100% accessibility scale; C, all peak tests with q565<0.05 highlighted. Heatmap rows are peaks, not TSS-centered genomic bins. The volcano uses a 0.5-count offset only to display zero-denominator ratios; test statistics and reported FC remain unmodified.
- **Figure 2:** one panel per nominated distinct peak; lines connect source-specific low/high accessibility fractions. Values are proportions of nuclei, not average read coverage. Titles identify regional q, which may differ from peak q.
- **Figure 3:** local target assessment. Display each neighborhood's named gene(s), plus the two candidates with largest absolute estimable full-population partial r; all 170 tests remain in source tables. A shows effect sizes before/after state control, B the corresponding link q, C RNA high/low FC on the 527 pairs. A colored RNA FC point alone does not signify a significant RNA contrast; full p/q are in the table. The display selection does not redefine correction families.
- **Figure S1:** opening and closing regional q for all 25 candidates.

## 6. Code and reproduction

Existing input acquisition, common peak construction and fragment QC are documented in `TASK1_ANALYSIS_WORKFLOW.md`. The frozen sensitivity grid was generated with scripts 88–89. The current executable workflow consumes those archived outputs and the deposited H5 files; it does not pretend to regenerate upstream alignment/QC.

```powershell
python scripts/92_fixed_day0_workflow.py --work <day0_sensitivity_grid_cache> --h5 <GSE208248_processed>
python scripts/93_plot_fixed_day0.py
```

Script 92 validates 565/25 test-family sizes and BH values, checks the 527 high/low barcodes, re-extracts all selected ATAC/RNA features, reproduces source counts from raw H5, and writes all target-model and RNA-contrast outputs. Regional permutation results are the validated archived D206 run, not newly permuted in script 92. Figure source tables and model tables are under `results/fixed_day0_D206`; PNG and vector PDF outputs are under `figures/fixed_day0_D206`. No private C2C12 data are included.

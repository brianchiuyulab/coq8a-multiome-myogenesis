# Results and decision record

## Primary 958-pair comparison

At TSS≥3 and COQ8A≥2 versus 1 UMI, the two muscle differentiation gene sets rise on RNA in all four libraries: Hallmark Myogenesis +0.01069 log1p(CP10K) per gene (`p=0.02684`, four-programme `q=0.03579`), Reactome Myogenesis +0.01145 (`p=0.02026`, `q=0.03579`), and four MRF loci +0.02764 (`p=0.001845`, `q=0.00738`). These are modest normalized-RNA programme effects and nucleus-level descriptive p values.

Across all nearby common ATAC peaks, Hallmark Myogenesis is +0.0675 percentage points with only 2/4 libraries positive (`p=0.0519`, `q=0.207`). Reactome Myogenesis is −0.0303 percentage points (`p=0.613`). No one of the 221 gene-region summaries has ATAC `q<0.05`, and no individual candidate peak has `q<0.05` under its full test family. The global ATAC analysis therefore does **not** show a programme-wide accessibility increase with higher COQ8A.

## How MYOD1 emerged

The fixed 221-gene search contains 7,699 eligible peak–gene pairs. RNA/peak association detected 15 positive MYOD1 links at same-data `q<0.05` across the link family; ten were evaluable in all four libraries and had a positive pooled link, and five of these ten did not have a strong scored MRF motif. MYOD1 led the exploratory combined gene ranking after link and COQ8A-direction filters. This nomination uses the same nuclei as the final locus test and is therefore a discovery finding. MYF5 was in the original 221 genes and was evaluated; the main all-four-common-peak ATAC summary did not show a signal like MYOD1.

For all 19 four-library common MYOD1 nearby peaks, the primary high/low open fractions are 12.801% and 11.977%, an absolute +0.824 percentage points and ratio **1.0688**, with 4/4 positive library directions (`p=0.0395`, 221-gene `q=0.635`). This is a small, uncorrected locus-level indication. It cannot be presented as an FDR-significant gene-region discovery.

## Sensitivity rather than a single favorable number

The linked five-peak subset under COQ8A≥3 versus 1 UMI and TSS≥3 gives 13.831% versus 11.144% open, absolute +2.687 percentage points, ratio **1.2411**, 201 matched pairs, 3/4 libraries positive, pair-level `p=0.05295`, and 16-setting exploratory `q=0.0847`. Under TSS≥2 the same definition gives ratio 1.2522, 4/4 positive, `p=0.04153`, `q=0.08153`. Under the primary ≥2 contrast and TSS≥3, the five-peak ratio is 1.1278, 3/4 positive, `p=0.03442`, `q=0.08153`.

After removing pairs touching the top 2.5% or 5% RNA Scrublet score tail, the TSS≥3/≥3 five-peak ratios are approximately 1.18 and 1.20 and pair-level p values ~0.16. This subset is therefore sensitive both to region selection and to high-scoring potential multiplets. The detailed 16-row grid and score-tail table are in `results/tables/`.

The earlier exploration reported **1.28**, `p≈0.026`, 4/4 directions for a **different inherited six-peak subset**. The clean 221-gene restart produces a five-peak subset and the **1.241**, `p≈0.053`, 3/4 result above. The previous number is not reused as a result of this pipeline.

## Interpretation for a paper

The data support an RNA association between COQ8A detection and a myogenesis programme. They provide a **local MYOD1 chromatin hypothesis** that can be visualized and followed up, but no stable broad accessibility claim across the knowledge-based myogenesis programme and no FDR-positive individual ATAC gene region or peak. The strongest same-data locus result is exploratory and weaker after a doublet-score challenge. These are two source lines, so neither nucleus-level p values nor four library directions constitute independent-donor replication. No causal inference that COQ8A drives chromatin opening follows from this observational comparison.

These results suggest a defensible presentation: Figure 1 shows the complete search and its broad negative ATAC finding; Figure 2 presents MYOD1 as a nominated, sensitivity-qualified locus. If the manuscript requires a decisive epigenetic mechanism claim, this public data set alone does not supply it.

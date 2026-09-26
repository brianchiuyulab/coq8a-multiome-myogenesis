# Results and decision record

## Primary 958-pair comparison

At TSS≥3 and COQ8A≥2 versus 1 UMI, the two muscle differentiation gene sets rise on RNA in all four libraries: Hallmark Myogenesis +0.01069 log1p(CP10K) per gene (`p=0.02684`, four-programme `q=0.03579`), Reactome Myogenesis +0.01145 (`p=0.02026`, `q=0.03579`), and four MRF loci +0.02764 (`p=0.001845`, `q=0.00738`). These are modest normalized-RNA programme effects and nucleus-level descriptive p values.

Across all nearby common ATAC peaks, Hallmark Myogenesis is +0.0675 percentage points with only 2/4 libraries positive (`p=0.0519`, `q=0.207`). Reactome Myogenesis is −0.0303 percentage points (`p=0.613`). No one of the 221 gene-region summaries has ATAC `q<0.05`, and no individual candidate peak has `q<0.05` under its full test family. The global ATAC analysis therefore does **not** show a programme-wide accessibility increase with higher COQ8A.

## How MYOD1 emerged

The fixed 221-gene search contains 7,699 eligible peak–gene pairs. Under TSS≥3 link discovery, RNA/peak association detected 15 positive MYOD1 links at same-data `q<0.05` across the link family; ten were evaluable in all four libraries and had a positive pooled link, and five of these ten did not have a strong scored MRF motif. Repeating link discovery at TSS≥2 gives **11** four-library positive links and **six** without a strong MRF motif. MYOD1 led the exploratory combined gene ranking after link and COQ8A-direction filters. This nomination uses overlapping nuclei as the final locus test and is therefore a discovery finding. MYF5 was in the original 221 genes and was evaluated; the main all-four-common-peak ATAC summary did not show a signal like MYOD1.

For all 19 four-library common MYOD1 nearby peaks, the primary high/low open fractions are 12.801% and 11.977%, an absolute +0.824 percentage points and ratio **1.0688**, with 4/4 positive library directions (`p=0.0395`, 221-gene `q=0.635`). This is a small, uncorrected locus-level indication. It cannot be presented as an FDR-significant gene-region discovery.

## Sensitivity rather than a single favorable number

**The earlier 1.28 result is reproducible.** The TSS≥2-linked six-peak subset, tested in TSS≥3 nuclei under COQ8A≥3 versus 1 UMI, gives 13.267% versus 10.365% open, absolute +2.902 percentage points, ratio **1.2800**, 201 matched pairs, 4/4 libraries positive, pair-level `p=0.02558`, and 20-setting exploratory `q=0.07551`. The same six peaks under TSS≥2 effect QC give ratio 1.2769, 212 pairs, 4/4 positive, `p=0.02522`. Under the primary COQ8A≥2 versus 1 contrast and TSS≥3, the six-peak ratio is 1.1097, 3/4 positive, `p=0.05304`.

Why the clean first pass reported five peaks: it required TSS≥3 **also during link discovery**. The sixth peak, `chr11:17652327-17652866`, is significant as a MYOD1 link when TSS≥2 matched nuclei are used for link discovery, but at TSS≥3 it is open in only eight matched nuclei in GSM6339603, below the fixed minimum of ten for a library-specific link test (14 at TSS≥2). The [gate-by-gate audit](../results/tables/myod1_link_gate_reconciliation.tsv) lists the count and link status for every MYOD1 candidate. The resulting five-peak subset, tested on the **same 201 TSS≥3 pairs**, gives 13.831% versus 11.144% open, ratio **1.2411**, absolute +2.687 percentage points, 3/4 positive, `p=0.05295`, 20-setting `q=0.07578`. In those same pairs, the extra sixth peak is open in **21 high versus 13 low** nuclei; adding it explains the 1.241→1.280 change without changing the test nuclei or fold formula. The TSS≥2 six-peak choice is a valid visible exploratory sensitivity, but its peak membership is learned from overlapping nuclei rather than an independent discovery cohort.

When pairs touching the top 2.5% or 5% RNA Scrublet-score tail are removed, the six-peak extreme-contrast ratio falls to **1.192** (187 or 173 pairs; p≈0.12 or 0.15). Thus the 1.280 signal is reproducible under its original definition, but sensitive to this QC challenge. The complete 20-row region/count/TSS grid and doublet-score table are in `results/tables/`. The previous exploration's BH q≈0.047 used a **different 32-test family that included COQ8A-zero comparisons**; the current displayed 20-test family gives q≈0.076. Those q values are not interchangeable.

## Interpretation for a paper

The data support an RNA association between COQ8A detection and a myogenesis programme. They provide a **local MYOD1 chromatin hypothesis** that can be visualized and followed up, but no stable broad accessibility claim across the knowledge-based myogenesis programme and no FDR-positive individual ATAC gene region or peak. The strongest same-data locus result is exploratory and weaker after a doublet-score challenge. These are two source lines, so neither nucleus-level p values nor four library directions constitute independent-donor replication. No causal inference that COQ8A drives chromatin opening follows from this observational comparison.

These results suggest a defensible presentation: Figure 1 shows the complete search and its broad negative ATAC finding; Figure 2 presents MYOD1 as a nominated, sensitivity-qualified locus. If the manuscript requires a decisive epigenetic mechanism claim, this public data set alone does not supply it.

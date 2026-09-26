# Figure legends

## Shared target comparison

Unless specified otherwise, GSE208248 COQ8A-high nuclei have ≥3 raw RNA UMI
and low nuclei have exactly 1 UMI. Joint RNA/ATAC QC includes TSS enrichment ≥3.
Within-library matching yields 201 pairs: 46, 74, 58 and 23 in line 1
undifferentiated, line 1 day 7, line 2 undifferentiated and line 2 day 7,
respectively. U denotes undifferentiated and D7 denotes day 7.
Low is grey and high is blue in paired point/line panels.
Tests use nucleus pairs; library points summarize measurements and do not
represent four independent donors. See Methods for all QC and test definitions.

## Figure 1. Externally defined differentiation regions and COQ8A associations

**A**, Biological and genomic selection. The 221 measured myogenesis genes
define candidate neighborhoods; 5,097 regions occur in all four target peak
catalogs. Independent differentiation data define 410 dynamic regions.
External functional annotations nominate a 54-region subset.
**B**, GSE109828 experiment-1 accessibility at 0, 24, 48 and 72 hours for all
410 regions. Each row is one region; accessibility is depth-standardized and
scaled within the row from its minimum (0) to maximum (1). Regions are ordered
by temporal class, external half-rise time and genomic coordinate. Early,
middle, late and closing contain 191, 34, four and 181 regions, respectively.
Class definition also requires concordant 0-to-72-hour direction in experiment 2.
**C**, Mean unscaled accessibility percentages across the member regions at
each external sampling time. These are descriptive mean profiles, not fitted
single-cell trajectories. Retained experiment-1 cell counts are 808, 995,
1,049 and 904; experiment 2 contributes 1,645 and 2,866 at its two endpoints.
**D**, COQ8A-high/low mean fraction of accessible member peaks in each target
library, aligned to the corresponding temporal class. Each connecting segment
joins high and low summaries within one library. FC is the ratio of pooled
group mean module scores. Two-sided paired t-test p and BH q across the ten
temporal programme tests in this fixed setting are shown. The panel is a
target association comparison, separate from the external time-course panels.

## Figure 2. Complete accessibility profiles of functional myogenic regions

**A**, ATAC fragment insertion profiles for all 54 functional dynamic regions,
centered on each peak midpoint (±2 kb, 100-bp bins). Each endpoint of a
deduplicated fragment contributes one insertion; fragment read multiplicity
is not used. Deposited Cell Ranger endpoints are already Tn5 adjusted and
receive no additional shift. Counts are normalized to 100 nuclei per 100 bp.
Gaussian smoothing with a one-bin standard deviation is applied only for
display. Both groups share the same absolute blue scale, capped at their
joint 99.5th percentile. Darker blue means stronger normalized insertion
signal. It is not a row-wise z score or an RNA expression value.
**B**, High-minus-low difference in the fraction of accessible nuclei at each
peak, in percentage points, for the four libraries. Red denotes increased
accessibility in high nuclei; blue denotes decreased accessibility. Values
beyond ±20 percentage points saturate the common color scale.
**C**, Pooled log2 ratio of high/low open-nucleus fractions, displayed on a
−3 to +3 color scale; larger magnitudes and infinite ratios saturate the
scale. A region inaccessible in both groups has an undefined ratio and is grey.
All numerical values are retained in source tables.

The same region order is used in all panels and is determined by external
timing, not COQ8A effect magnitude. Sidebars indicate early (teal, 24 regions),
middle (purple, five) and closing (orange, 25). Region IDs R01–R54 map uniquely
to coordinates in `functional_54_order.tsv`. Repeated gene labels identify
distinct nearby regions, not duplicated tests; a semicolon indicates multiple
candidate neighboring genes. These labels are proximity annotations, not
validated enhancer-to-gene assignments. All 54 regions are shown without
filtering on target p value, fold change or direction.

## Figure 3. Candidate loci and paired RNA measurements

Rows show CSRP3, CAV3, MYOD1 and CACNA1H follow-up loci.
**A**, Actual fragment insertion profiles using the same normalization and
display smoothing as Figure 2. The highlighted interval is the tested peak.
GRCh38 genomic coordinates are shown in kilobases. The gene model uses the
GENCODE v48 transcript whose TSS is nearest the peak midpoint; ties are
resolved by transcript ID. Exons are boxes and the arrow indicates strand.
Models are clipped to the plotted window. The same y scale is used for high
and low within a locus, but differs between loci.
**B**, High and low open-nucleus percentages within each library. FC and
two-sided exact paired-binomial p are calculated from the pooled matched
nuclei; q is BH adjusted across the **54 functional peaks** in the same setting.
**C**, Mean RNA loge(1 + counts per 10,000 RNA UMI) in the same groups and
libraries. Δ is the pooled high-minus-low difference on this scale; p is
from the two-sided paired t test and q is BH adjusted across **221 genes**.

The four rows are exploratory follow-up examples from the complete screen,
selected for biological interpretation and cross-dataset comparison, including
the strongest FDR-supported decrease within the functional family. They are
not the complete set of nominally significant peaks and are not represented
as four FDR-positive peak-to-gene links. Nearby position, differential
accessibility and RNA association are distinct measurements.

## Supplementary Figure S1. Myogenesis gene TSS profiles

ATAC insertion profiles are aligned to one annotated GENCODE gene-feature TSS
for each of the 221 genes, ±5 kb in 100-bp bins; negative-strand genes are
reversed. Both groups use the same absolute color scale and the same row
order, determined by pooled signal. The matched nuclei are the same 201
pairs as in the main figures. Smoothing is for display only. This promoter
context does not determine peak selection. The all-window insertion ratio
is distinct from the ratio of open-nucleus fractions at tested peaks.

## Supplementary Figure S2. Module sensitivity to QC and COQ8A definitions

Temporal classes (**A**) and eight function-by-opening/closing modules (**B**)
are evaluated at TSS ≥2 or ≥3 and COQ8A ≥2 or ≥3 versus exactly 1 UMI.
Cells display the module FC and nominal paired-test p value; color denotes
log2 FC on a shared scale. Each setting uses its own within-library matched
pairs (1,020, 958, 212 and 201 in displayed column order). Complete q values
are available in the source tables: ten temporal programmes or eight
functional modules per fixed setting. Sensitivity settings are not combined
into a single BH family. The two external criteria remain fixed at a 2-kb
marker-promoter window and 50% reciprocal peak overlap.

## Visual-design references

The combination of aligned signal heatmaps, metaprofiles, locus tracks and
paired-condition summaries follows established chromatin-data presentation:
[Martini et al., Nature (2026)](https://www.nature.com/articles/s41586-026-10791-2),
[muscle chromatin study, Nature Communications (2023)](https://www.nature.com/articles/s41467-023-42313-3),
and [Signac, Nature Methods (2021)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9255697/).
The figures here are generated from this repository's data; no published
image is reused. ATAC insertion signals are not labelled as histone-mark signals.

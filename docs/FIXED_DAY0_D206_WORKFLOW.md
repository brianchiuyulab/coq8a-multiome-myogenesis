# Fixed Day-0 COQ8A high/undetected analysis

## Question and fixed design

Test whether COQ8A-high nuclei differ in local accessibility near externally defined myoblast differentiation/fusion genes, then assess possible local RNA targets. The comparison is D206 throughout: GSE208248 undifferentiated cultures only, TSS enrichment >=3, COQ8A >=2 versus 0 raw UMI, within-library RNA/ATAC depth matching, caliper 0.30. No marker-positive population restriction or state matching is imposed in the primary contrast. Zero means undetected RNA, not proven absence of expression. These are cultures derived from human muscle stem cells, not newly annotated pure quiescent MuSC nuclei.

The same setting is applied to every candidate. External gene membership defines the search scope. GSE240061, Day 7, private C2C12 RNA and temporal/motif labels do not determine eligibility here. Parameter comparisons are available in the separate sensitivity analysis report.

## 1. External scope and peak construction

Frozen Hallmark/Reactome myogenesis memberships define 221 measured genes. Intersect with the union of GO myoblast differentiation and GO myoblast fusion: 17 + 12 - 4 shared = 25 genes. These are process annotations, not a list exclusively of positive regulators. All members and provenance are in `results/differentiation_fusion25/genes.tsv` and `reference/myogenesis_221_gene_sources.tsv`.

Start from deposited GRCh38 author-called peaks. Anchor on GSM6339597 and match the other three original catalogs one-to-one with >=50% reciprocal overlap. Apply coordinate QC to the anchor intervals: canonical chromosomes, widths 200–2,000 bp, no ENCODE hg38 v2 blacklist overlap, and peak midpoints within 100 kb of an annotated protein-coding transcript TSS. GTF TSS positions are converted to 0-based coordinates before distance calculation. The broad scope has 5,097 common peaks; restricting neighborhood membership to the 25 genes gives 565 distinct peaks. A peak may belong to more than one neighborhood. Catalog construction uses all four original libraries; effect tests below use only Day 0. No FASTQ alignment or new peak calling is claimed.

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
- **Figure 2:** six q<0.05 neighborhoods, one panel per representative peak. Lines connect source-specific low/high accessibility fractions. Blue is source 1 (226 pairs); red is source 2 (301 pairs). Low is 0 COQ8A UMI; high is >=2 UMI. FC is the pooled accessible-nucleus fraction ratio, and q is the regional directional BH25 value. Neighborhood names indicate proximity, not a confirmed target. Secondary neighborhoods are shown separately in Figure S2.
- **Figure 3:** local target assessment. Display each neighborhood's named gene(s), plus the two candidates with largest absolute estimable full-population partial r; all 170 tests remain in source tables. Row labels mean genomic neighborhood to candidate RNA, not a causal arrow. Panel a uses all 16,881 eligible nuclei; panel b uses the exact 1,054 matched nuclei. Both show partial r before/after state control on the same color scale; an asterisk denotes link q<0.05. NA marks insufficient detection in at least one source. Panel c shows RNA FC on the 527 pairs: filled points denote RNA q<0.05, hollow points q>=0.05. Red/blue indicate increased/decreased RNA. Display selection does not redefine correction families.
- **Figure S1:** opening and closing regional q for all 25 candidates. **Figure S2:** secondary neighborhoods (0.05<=regional q<0.10); MYF5/MYF6 share a peak.
- **Figure 4:** positional profiles based on deduplicated ATAC fragments from the same 527 pairs. Left: the 25 gene-feature TSSs, strand aligned, +/-5 kb in 100-bp bins. Right: all 565 peak midpoints, +/-2 kb in 50-bp bins. Upper panels show average insertion profiles at bin centers, using identical high/low y limits; lower panels use identical row order and color limits for high/low. Rows are sorted by mean signal across both groups. Colors clip at the pooled 99th percentile within each panel family; raw values are retained. Counts are insertions per 100 nuclei per 100 bp, averaged equally over the two sources after source normalization. These are positional profiles, not tests of a genome-wide opening shift. A gene's TSS window need not include its nominated distal peak.
- **Figure 5:** MYOD1, CSRP3 and CAV3 locus tracks from the same fragment extraction, with common high/low axes at each locus. Gold shading marks the nominated peak; the dotted line marks the gene-feature TSS. Bin width is 100 bp, with one-bin Gaussian smoothing for display only. The tracks show source-normalized insertion density, a different measurement from the accessible-nucleus fraction used for FC tests. The panel structure follows the aggregate-profile, positional heatmap and locus-track organization of Martini et al.'s supplied SASP paper, Figure 3o-p; no ChIP signal is substituted for ATAC.

## 6. Code and reproduction

The current entry point is `run_fixed_day0.py`. It applies the fixed design to
all candidates, recomputes the complete peak screen and 2,000,000 regional
permutations, fits local RNA models, and renders the current figure set.

```bash
python run_fixed_day0.py --mode figures
python run_fixed_day0.py --mode analysis --work /work/day0_cache --h5 /data/GSE208248
python run_fixed_day0.py --mode audit --work /work/day0_cache --h5 /data/GSE208248 --gtf /data/gencode.v48.annotation.gtf.gz --permutations
```

[Code availability](CODE_AVAILABILITY.md) documents upstream rebuilding and
fragment extraction. [The computational audit](FIXED_DAY0_AUDIT.md) records
raw-matrix checks and coordinate corrections. Figure source tables and model
outputs are in `results/fixed_day0_D206`; figures are in
`figures/fixed_day0_D206`. The older temporal runner is documented separately.

Fragment plots count each deduplicated row's start and end-1 once, respecting
BED half-open coordinates. PCR multiplicity is not multiplied and no second
Tn5 shift is added. Minus-strand TSS windows are strand aligned using the same
relative-bin convention as plus-strand windows. Values are plotted at bin
centers; a final partial locus bin uses its actual width for normalization.

Figures use sans-serif type at approximately 5.5–7 pt at their saved dimensions,
with larger panel letters, 400-dpi PNG previews and editable PDF/SVG exports.
The layout follows the aggregate profile / positional heatmap / locus track
organization in the supplied Martini et al. SASP Figure 3o-p. Typography and
export choices follow the [Nature figure guide](https://research-figure-guide.nature.com/figures/building-and-exporting-figure-panels/).

Private C2C12 cross-references use script 95 and an enforced external output
path. Same-symbol matches are recorded explicitly; unmatched features remain
unresolved. Private counts and derived figures are not included here.

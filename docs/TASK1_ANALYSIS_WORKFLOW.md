# Task 1: functional-scope screen of COQ8A-associated accessibility

The subsequently requested Day-0-only primary analysis is complete in
[TASK1_DAY0_RESULTS.md](TASK1_DAY0_RESULTS.md). The pooled 958-pair results
documented below retain all four libraries and are not the Day-0 statistics.

## Study question and analysis status

Do COQ8A-high nuclei show increased accessibility at local regulatory candidates
near genes involved in myoblast differentiation or fusion, and do the nominated
genes show concordant RNA expression?

This document specifies the current GSE208248 workflow from processed inputs to
candidate nomination. It reorganizes completed analyses; it does not describe a
new download, a new run, or a prospectively registered discovery. The external
memberships and genomic rules do not use COQ8A effect estimates. Selection of the
reported COQ8A contrast and the regional test occurred during exploration.

GSE240061 is retained as a separate completed exploration and is not a selection
requirement or a successful replication of this workflow.

## Workflow

```text
External biological scope                         Paired target measurements
Hallmark / Reactome myogenesis union              GSE208248: four libraries
221 measured genes                               RNA + ATAC, same nucleus
             |                                                |
Intersect with GO myoblast differentiation        Joint nucleus QC
OR GO myoblast fusion                            Within-library COQ8A high/low
17 + 12 - 4 = 25 genes                            Depth matching
             |                                                |
Author peak catalogs -> genomic QC -> common four-library map |
5,097 peaks near the broad 221-gene scope                       |
             |                                                |
Retain all peaks near the 25 genes: 565 distinct peaks ---------+
             |
Test all 25 gene neighborhoods with the same local-opening test
Within-region max-statistic permutation -> BH over 25 regions
             |
Rank all regions; inspect constituent peaks and library directions
CSRP3 ranks first in TSS >=3, COQ8A >=2 vs 1
             |
For nominated peaks: all local RNA candidates -> peak-RNA association
Separate same-pair RNA response of candidate targets
             |
Nominate CSRP3 as an exploratory candidate
```

There is no requirement for the average of all 565 peaks to increase before a
local association can be examined. The broad inventory and the local test answer
different questions. No observed broad-screen p value is used to construct the
25-gene membership. Early/middle/late labels, MYOD1 occupancy, MRF motifs and
external muscle-selective enhancer tracks are not eligibility filters here.

## 1. Inputs and sample design

GSE208248 contains two human donor-derived muscle stem cell lines, RH-hMuSC-1
and RH-hMuSC-2, each assayed undifferentiated and after seven days of
differentiation. Use the four paired Multiome libraries:

| Source | Condition | RNA / matrix accession | ATAC accession | Matched pairs in the reported setting |
|---|---|---|---|---:|
| Line 1 | Undifferentiated | GSM6339597 | GSM6339598 | 211 |
| Line 1 | Day 7 | GSM6339599 | GSM6339600 | 326 |
| Line 2 | Undifferentiated | GSM6339601 | GSM6339602 | 296 |
| Line 2 | Day 7 | GSM6339603 | GSM6339604 | 125 |
| Total | | | | 958 |

Start with deposited filtered_feature_bc_matrix.h5 files, per-barcode metrics
and indexed ATAC fragments, using the existing local copies. File hashes are
recorded in input_sha256.csv. The reference genome is GRCh38. Author alignment
and peak calling are inherited; this workflow does not rerun FASTQ processing.
Biological source N=2, library N=4. The 958 pairs are 1,916 distinct nuclei,
not 958 independent donors.

RNA and ATAC are linked by library plus barcode. Comparisons are within each
library, so high/low matching never crosses cell line or differentiation stage.
The analysis asks about COQ8A association within these conditions; it is not
simply an undifferentiated-versus-day-7 contrast.

## 2. External functional scope, before target-effect ranking

Use frozen MSigDB memberships rather than manually naming candidate genes:

- HALLMARK_MYOGENESIS and REACTOME_MYOGENESIS form a union of 221 measured genes:
  198 Hallmark genes plus 29 Reactome genes minus six shared members.
  Hallmark lists 200 genes; two are absent from the target RNA feature matrix.
- Intersect that universe with the union of
  GOBP_MYOBLAST_DIFFERENTIATION (GO:0045445) and
  GOBP_MYOBLAST_FUSION (GO:0007520).
- The resulting sets have 17 and 12 genes, sharing ITGB1, MAPK14, MYOD1 and
  MYOG: 25 distinct genes in total.

The 25 genes are ADAM12, ANKRD2, BOC, CACNA1H, CAV3, CDON, CSRP3, IFRD1,
IGF1, IGFBP3, ITGB1, KCNH1, MAPK12, MAPK14, MEF2C, MYF5, MYF6, MYH9,
MYOD1, MYOG, NEO1, NOS1, NOTCH1, RB1 and TGFB1. CSRP3 enters through
differentiation membership; it is not inserted because of its measured effect.

These annotations define processes, not uniformly positive regulators of
differentiation. Increased RNA of every member is therefore not automatically
interpreted as improved differentiation.

Membership records: reference/myogenesis_221_gene_sources.tsv,
reference/muscle_subprogrammes.json, and
results/differentiation_fusion25/genes.tsv. Use the frozen files to reproduce
the counts; a future database release may have different membership.

## 3. Joint nucleus QC and expression grouping

Apply RNA and ATAC QC to the same nucleus before COQ8A matching:

| Parameter | Fixed implementation |
|---|---|
| Total RNA depth | At least 500 UMI |
| ATAC depth | At least 500 detected peaks, not 500 fragments |
| Upper depth tails | Remove above the within-library 95th percentile of either depth, calculated after both lower-depth filters |
| FRiP | At least 0.25; peak-region fragments / total ATAC fragments |
| Nucleosome signal | Below 4; [147,294)-bp fragments / fragments shorter than 147 bp |
| Blacklist fragment fraction | Below 0.05, ENCODE hg38 v2 |
| TSS enrichment | At least 3 for the reported setting; at least 2 retained in the existing sensitivity screen |
| COQ8A high / low | Raw RNA UMI at least 2 / exactly 1; zero excluded from this contrast |

TSS enrichment uses the archived custom implementation: gene-level GENCODE v48
protein-coding TSSs, central +/-500-bp insertion density divided by density in
the two 100-bp flanks at distances 901-1,000 bp. Library-pooled flank density is
used when a nucleus has zero flank insertions. This is not a claim that TSS=3
is interchangeable across software implementations.

Mitochondrial RNA fraction, detected-gene counts and Scrublet tails were separate
sensitivities, not baseline exclusions. Comprehensive doublet removal is not
part of the archived baseline. Harmony, trajectory inference and HMA annotation
transfer are not used in the count comparisons.

The QC values, 100-kb window and COQ8A UMI definitions are analyst-specified
settings, not all exact copies of the source paper. Their rationale is coverage,
technical comparability and candidate regulatory proximity; the exact values
remain explicit methodological choices. A two-versus-one raw UMI contrast does
not imply a twofold normalized-expression difference for every matched pair.

## 4. Within-library matching

Pair high and low nuclei 1:1 without replacement within each library. Require
absolute differences <=0.10 on both log1p(total RNA UMI) and log1p(number of
detected ATAC peaks). Search up to 100 nearest low nuclei in this two-dimensional
space, process high nuclei from largest nearest-low distance downward, and
accept the first unused candidate satisfying both limits. This yields 958 pairs
for TSS>=3 and COQ8A>=2 vs 1. Freeze the barcode list for all corresponding ATAC
and RNA contrasts. Do not rematch separately for favorable genes.

## 5. Peak processing and gene-neighborhood assignment

1. Read author-called peak coordinates from each of the four H5 matrices.
2. Anchor coordinates on GSM6339597. Match each other catalog one-to-one using
   at least 50% reciprocal overlap. These are corresponding author intervals,
   not peaks called again on pooled data.
3. Retain canonical chromosomes, interval widths 200-2,000 bp and no ENCODE
   blacklist overlap.
4. For the broad 221-gene universe, retain peak midpoints within 100 kb of at
   least one annotated protein-coding transcript TSS. Require representation
   in all four catalogs: 5,097 distinct peaks.
5. Retain all peak-gene memberships for the fixed 25-gene scope: 565 distinct
   peaks. Do not remove peaks because of weak effect, unfavorable direction,
   absent MRF motif or lack of a temporal label.

A peak can belong to multiple neighborhoods. Count each coordinate once in the
565-peak screen, while preserving its membership in each applicable regional
test. Neighborhood membership is a candidate association, not proof of a target
gene or a peak-RNA link. The common-catalog rule does not require every nucleus
to have an accessible signal at every peak.

## 6. Regional tests and multiplicity

The regional question is whether a neighborhood contains a local accessibility
increase, rather than whether all its peaks increase on average.

For each peak, let D be the vector of binary-accessibility differences
(high minus low) across matched pairs. A matrix count >0 is accessible. Compute
Z=sum(D)/sqrt(sum(D^2)); use zero when there are no discordant pairs. For each
gene neighborhood, the opening statistic is max(0, max Z over its peaks).

Generate 2,000,000 within-pair high/low swaps with seed 20260928. Apply the same
swap vector across all peaks to preserve their dependence. Recompute each
neighborhood's maximum on every permutation. Use (exceedances+1)/(B+1) for its
p value, then apply Benjamini-Hochberg across all 25 regional p values.

This provides two distinct multiplicity steps: the regional permutation accounts
for searching among its constituent peaks; BH accounts for searching across
25 neighborhoods. It is not BH over 25 selectively retained single-peak p values.
No region is removed before correction for having a weak p value or negative
direction. The exploratory q<0.10 reporting level is not q<0.05.

Rank all regions by the same opening p value. Report all 25, alongside a
two-sided max-absolute-statistic sensitivity and the alternate >=3 vs 1 contrast.
Existing single-peak results retain their separate BH correction over 565 peaks.
These are paired-nucleus association tests, not population-level donor tests.
Within-setting q values do not adjust for the preceding choice among analytical
strategies or COQ8A definitions.

## 7. Candidate nomination and RNA follow-up

### Separate locus nomination from target identification

The 25 scope genes define where to search for ATAC associations. They do not
restrict which RNA a nominated peak can regulate. Once a locus is nominated,
target identification examines all measured local transcripts with an annotated
TSS within 500 kb, including genes outside the 25-gene scope and lncRNAs. The
25-gene high/low RNA screen is a separate supplementary expression comparison;
it does not identify the target of a peak.

For chr11:19218592-19219518, the archived unrestricted map contains 12 RNA
candidates. At the 5% peak/RNA detection criterion, only NAV2 and ZDHHC13 have
combined link estimates and neither has a significant association. At the 1%
sensitivity, CSRP3 and E2F8 additionally become evaluable. With technical and
COQ8A adjustment, CSRP3 r=0.02259, p=0.000346 and E2F8 r=0.02519, p=0.000066.
After additional myogenesis-state adjustment, CSRP3 r=0.00248, p=0.695 and
E2F8 r=0.02298, p=0.000273. These small correlations do not establish a
regulatory target. E2F8 is an alternative target-association candidate, not a
replacement proven mechanism.

The source link q values are retained from the complete archived link families
(61,939 links for each 1% model; 27,720 for each 5% model), not recalculated over
the four displayed genes. Missing link estimates are not evidence of absence.
All 12 annotations and existing estimates are exported in
results/local_signal_decision/CSRP3_neighborhood_focal_all_local_RNA_candidates.tsv
and CSRP3_neighborhood_focal_existing_RNA_links.tsv. This is a provenance
extraction from script 48 outputs, not a new 958-pair link fit.

The interpretation sequence is therefore: nominate an ATAC locus, evaluate its
possible RNA targets, then examine those targets' COQ8A-associated RNA response.
CSRP3-neighborhood nomination alone does not identify CSRP3 as the regulated RNA.

### Regional ranking and supplementary RNA evidence

At TSS>=3, COQ8A>=2 vs 1, the regional opening ranking begins:

| Region | Constituent peaks | Opening p | Regional BH q |
|---|---:|---:|---:|
| CSRP3 | 10 | 0.003825 | 0.095625 |
| TGFB1 | 27 | 0.010905 | 0.136312 |
| MAPK12 | 27 | 0.024477 | 0.203979 |
| CDON | 24 | 0.046559 | 0.290997 |
| CAV3 | 19 | 0.060006 | 0.300030 |

CSRP3 is nominated after this complete regional comparison. Its strongest
constituent peak is chr11:19218592-19219518: 139/958 high versus 92/958 low
nuclei are accessible, or 14.51% versus 9.60%, a 4.91-percentage-point difference
and FC=1.51087. All four library directions are positive. Its individual-peak
two-sided p=0.001190 and q565=0.3361 remain distinct from regional q=0.095625.
The CSRP3 two-sided regional q is 0.19325; in the >=3 vs 1 opening analysis,
regional q is 0.6718.

Compare RNA for all 25 genes in the same 958 pairs, normalize to CP10k, apply a
paired t test on log1p(CP10k), and BH-correct over 25 RNA tests. Report FC as the
ratio of linear mean CP10k values. CSRP3 RNA FC=2.43259, p=0.012362,
q25=0.104054. Its RNA direction is positive in two of four libraries.

Peak-RNA partial correlations are a separate follow-up from the archived
QC-passing-nucleus link analysis; they are not the 958-pair RNA contrast. The
CSRP3 link is evaluable only with the 1% detection sensitivity, not the 5%
criterion: technical/COQ8A-adjusted r=0.0226, and additionally state-adjusted
r=0.00248, p=0.695. Thus the nearby peak, the high/low RNA difference and direct
cis association are not interchangeable evidence. The current nomination is a
CSRP3-neighborhood accessibility/RNA candidate, not an established enhancer-to-
CSRP3 mechanism.

Private C2C12 differentiation-time contrasts are checked after nomination and
remain outside the public repository. They provide a direction-of-expression
comparison, not a COQ8A overexpression test or a human-peak/mouse-peak mapping.

## 8. Presentation order

Recommended panels for this Task 1 result, to be rendered from this setting:

1. External-set and peak-selection flow: 221 -> 25 genes; 5,097 -> 565 peaks.
2. Joint QC and within-library depth balance for the 958 matched pairs.
3. Complete 25-region ranking with regional q values and an accompanying
   all-565-peak effect overview; no significance filter before display.
4. All 10 CSRP3 peaks, followed by the focal locus and four-library open rates.
   A peak-centered fragment heatmap describes peak profiles; a TSS-centered
   heatmap describes promoter profiles. They are not labeled interchangeably.
5. Same-pair CSRP3 RNA with the corresponding RNA statistics and library detail.

These are panel specifications, not a claim that the repository's earlier
temporal figures already display the new regional results. High and low fragment
heatmaps must share row order, normalization and color scale. Figures carry
scientific labels; analysis history belongs in this document and legends.

## 9. Source code and reproduction boundary

| Stage | Source |
|---|---|
| Reference and fragment QC | scripts/00_prepare_reference.py; scripts/00_fragment_qc_reference.py |
| Nucleus extraction | scripts/01_extract_nuclei.py |
| Joint QC and all retained pairs | scripts/02_qc_and_matching.py |
| Common peak map and broad neighborhoods | scripts/03_define_regions.py |
| Complete original peak tests | scripts/05_peak_effects.py |
| Broad local RNA/ATAC extraction and links | scripts/48_unstratified_links.py |
| Frozen functional membership and 565-peak comparisons | results/differentiation_fusion25/genes.tsv; scripts/76_differentiation_fusion_analysis.py |
| Regional permutation and same-pair RNA tests | scripts/84_gene_region_permutation.py |

With the archived extraction directory containing peaks.txt, genes.txt and
per-library *_meta.tsv, *_atac.npz and *_rna.npz, reproduce regional/RNA results:

```bash
python scripts/84_gene_region_permutation.py --work GSE208248_LINK_WORK --permutations 2000000 --seed 20260928
```

Outputs are results/local_signal_decision/GSE208248_gene_region_tests.tsv and
GSE208248_2plus_vs_1_RNA.tsv (plus the >=3 contrast). The command consumes the
committed membership and matched-pair tables; it does not rebuild every upstream
stage. run_paper.py still implements the earlier temporal figure workflow and
must not be described as a one-command runner for this new regional analysis.

The functional membership, neighborhood map and matched-pair files are the
frozen interfaces between stages. Existing inputs and statistics have not been
overwritten to make this workflow appear prospectively selected.

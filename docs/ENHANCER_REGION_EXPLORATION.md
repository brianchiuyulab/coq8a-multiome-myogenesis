# Enhancer selectivity and region-level COQ8A association

## Sample size and the existing MYOD1 peak test

GSE208248 contributes two donor-derived cell lines, each sampled in the stem
and day-7 differentiated conditions. Four libraries are therefore not four
independent donors. Fixed TSS >=3 and COQ8A >=3 versus exactly 1 UMI matching
provides 201 pairs, or 402 distinct nuclei. Matching is within library without
reusing nuclei; neither age nor a temporal class selects these pairs.

| Library | Source | Condition | Pairs | Focal high open | Focal low open | FC |
|---|---|---|---:|---:|---:|---:|
| GSM6339597 | Line 1 | Stem | 46 | 11 | 8 | 1.375 |
| GSM6339599 | Line 1 | Differentiated | 74 | 11 | 4 | 2.750 |
| GSM6339601 | Line 2 | Stem | 58 | 21 | 14 | 1.500 |
| GSM6339603 | Line 2 | Differentiated | 23 | 1 | 1 | 1.000 |

For chr11:17649919-17650798 (GRCh38), 44 high nuclei and 27 low nuclei have
detectable accessibility. There are 33 high-only and 16 low-only pairs.
The exact two-sided paired binomial/McNemar test uses these 49 discordant
pairs and gives p=0.0212941. It is a conditional matched-nucleus association
test, not a test across 201 independent donor samples. The three positive
library directions and one tie should not be described as four independently
significant replications. The earlier three-library RNA eligibility is a
different analysis and does not reduce this ATAC comparison to three libraries.

Among the previously defined 23 core-MRF peaks, this p ranks first and the
BH adjusted value is 0.489765. In this particular table the BH running-minimum
calculation leaves 23 x 0.0212941 as the adjusted value. BH is not generally
equivalent to multiplying every p by the number of tests; its adjustment also
depends on the ranked values of the other tests. The 23 peaks are outcomes,
not the biological sample size.

## A TF-independent external reference filter

The initial 1,777 peaks overlap HSMM strong-enhancer states 4/5. To evaluate
cell-context selectivity without a MYOD occupancy filter, the same published
hg19 ChromHMM states were obtained for all eight other cell types in this
ENCODE/Broad track collection: GM12878, H1-hESC, HepG2, HMEC, HUVEC, K562,
NHEK and NHLF. Source URLs and SHA256 checksums are stored in
`results/enhancer_regions/references.json`.

The strict set consists of HSMM strong-enhancer peaks with no positive
overlap with a strong-enhancer call in any of the eight comparators. A relaxed
sensitivity set allows one comparator. Both rules were specified in code before
their COQ8A associations were calculated. These are exploratory choices, not
retrospectively preregistered analyses. Absence of a strong state in these eight
references is not absence of all enhancer activity in all non-muscle tissues.

| External set | Peaks | Mean-accessibility high/low FC | Exact paired p | BH q across these three sets |
|---|---:|---:|---:|---:|
| All HSMM strong enhancers | 1,777 | 1.0043 | 0.6058 | 0.6058 |
| No strong overlap in the eight comparators | 401 | 1.0289 | 0.1707 | 0.2560 |
| At most one comparator | 750 | 1.0316 | 0.03858 | 0.1157 |

The MYOD1 focal peak survives the strict 401-peak set: it has no strong-enhancer
overlap in any of the eight comparator tracks. Its nominal p and FC do not
change, and its peak-level q across 401 peaks is 1. No peak is selected by its
COQ8A effect, RNA correlation, or MYOD binding in this filter.

## Defining the unit of inference before testing

Two complementary region definitions are reported in full:

1. **External enhancer blocks:** merge overlapping or directly adjacent HSMM
   state-4/state-5 intervals, with no extension across gaps. Assign a target peak
   to its largest-overlap block, with a genomic-start tie-break. This prevents
   double counting a peak that spans more than one reference block. A target
   peak contributes its existing binary detection, not newly clipped fragment
   counts inside the overlap.
2. **Gene neighborhoods:** use all retained peaks associated with each gene in
   the original candidate-neighborhood table. This is a spatial assignment,
   not proof that the enhancer regulates the named gene. All neighborhoods are
   tested, rather than selecting MYOD1 or the four MRF genes in advance.

For each unit, a nucleus's score is the fraction of its member peaks detected.
FC is the ratio of mean scores in the 201 matched high and low nuclei. The
two-sided test exchanges high/low jointly for all peaks within each pair,
preserving dependence among those peaks. Dynamic programming enumerates the
exact weighted sign-sum distribution using integer count differences, so no
Monte Carlo resolution or floating-point tie decision affects p values.
This tests a regional score and does not establish significance of every peak
inside that region. Pair exchangeability is an assumption of this observational
conditional test; it does not account for uncertainty in sampling donor lines.

A Simes combination of constituent exact peak p values is also reported as a
sparse-signal endpoint, in contrast to the mean score. It tests whether at least
one constituent is associated under the usual Simes dependence assumptions.
BH correction is calculated across all units within each region definition
and endpoint. These are distinct exploratory families, not a search for the
smallest q across methods. Neighborhood sets can overlap.

Source-level effect tables summarize the two stages with equal library weights
within each source. They are descriptive checks of direction; no donor-level
p value is claimed with two sources. Existing p values are retained in their
original files and are not renamed as donor-level results.

| Scope | Peaks | External blocks | Gene neighborhoods |
|---|---:|---:|---:|
| All HSMM strong enhancers | 1,777 | 1,479 | 209 |
| Strict comparator-selective set | 401 | 370 | 148 |

The code checks matched-nucleus uniqueness and reproduces the original
per-peak high/low detection totals. Exact score tests for one-peak units agree
with the existing binomial tests. The dynamic-programming calculation was also
checked against complete sign enumeration on small weighted examples.

## Results and priorities

The MYOD1 distal external block is hg19 chr11:17670824-17676224. It includes
three GRCh38 target peaks: chr11:17649919-17650798,
chr11:17652327-17652866 and chr11:17653349-17654252. All three survive the
strict external selectivity rule. Its regional FC is **1.28235**, exact
mean-score **p=0.070328**, and **q=1** across either 1,479 or 370 blocks.
Simes p=0.063882. Both donor-derived sources have positive regional effects.
The preliminary floating-point permutation p (~0.061) was superseded by this
integer-exact result, which handles tied permutation statistics exactly.
This three-peak block is unrelated to the historical six-peak FC of about 1.28.

Testing all 12 MYOD1-neighborhood strong-enhancer peaks gives FC=1.11475 and
p=0.18546. Within the strict selective set, the ten MYOD1-neighborhood peaks
give FC=1.13359 and p=0.15240. Neither regional definition establishes a
significant MYOD1 accessibility programme. The focal peak remains a biologically
annotated exploratory candidate, rather than an FDR-positive MYOD1 discovery.

The strongest strict-set block is **hg19 chr11:1900624-1903824**, containing
GRCh38 chr11:1880073-1880969 and chr11:1882161-1883026. Its regional
FC is **0.388889**, p=**3.52160e-5**, and q=**0.013030** across the 370 blocks.
In the broader 1,479-block analysis, the same p yields q=0.052084. The original
neighborhood annotation lists **LSP1, TNNI2 and TNNT3**; these are candidates,
not experimentally established target genes. Equal-stage source-level ratios
are 0.3368 and 0.3760, both decreasing. This is the only q<0.05 block in the
strict-set regional mean-score analysis, and is an exploratory nucleus-level
closing association, not a population-level donor inference or a demonstrated
mechanism of improved differentiation.

The next analytical priority is to examine enhancer-to-gene evidence and RNA
coupling for both directions, retaining MYOD ChIP as a supporting annotation
after the broad screen. The new external selectivity rule supports a coherent
401-peak search space, but does not guarantee larger effects or rescue MYOD1 FDR.

## Reproduction

```bash
python scripts/54_enhancer_region_tests.py --hmm /references/hsmm_HMM.bed.gz --cache /work/unstratified --out results/enhancer_regions
python scripts/55_external_enhancer_selectivity.py --references /references/nonmuscle_hmm --cache /work/unstratified --out results/enhancer_regions
python scripts/54_enhancer_region_tests.py --hmm /references/hsmm_HMM.bed.gz --cache /work/unstratified --selectivity results/enhancer_regions/peak_selectivity.tsv --out results/enhancer_selective_regions
```

## Methodological references

- [Best practices for differential accessibility analysis in single-cell epigenomics, Nature Communications 2024](https://www.nature.com/articles/s41467-024-53089-5): biological-replicate-aware inference is preferred for population-level claims. The matched-nucleus tests here are a separate exploratory estimand, not a substitute for replicate-level pseudobulk analysis.
- [csaw: correction for multiple testing](https://bioconductor.org/books/release/csawBook/correction-for-multiple-testing.html): distinguishes region-level from window-level claims and describes Simes combination. This analysis uses existing binary peak calls and custom exact paired score tests; it is not a csaw reanalysis of raw ChIP fragments.
- [ENCODE/Broad ChromHMM track](https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&g=wgEncodeBroadHmm): independent chromatin-state definitions and comparator panel.

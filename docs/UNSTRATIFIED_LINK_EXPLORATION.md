# Temporally unstratified peak–RNA exploration

## Design

All **5,097 candidate peaks shared by four libraries** were retained from the
existing 221-gene myogenesis-neighborhood definition. External early/middle/late
or closing labels are not used for inclusion, ranking or testing. The 2,035
additional original candidates present in only three libraries are outside this
fixed four-library universe. This is a complete screen within the specified
myogenesis neighborhoods, not a whole-genome screen.

COQ8A-high remains ≥3 raw RNA UMI and low remains exactly 1 UMI. Existing joint
QC with TSS enrichment ≥3, within-library depth matching and all **201 pairs**
are unchanged. Each peak is tested using the exact paired binomial test on
discordant pairs. Opening and closing directions are both retained. Reported
ATAC FC is the fraction of high nuclei with a detected peak divided by the
corresponding fraction of low nuclei. Infinite ratios are not replaced with
arbitrary finite effect sizes; the complete-screen figure uses the percentage-point
difference so that zero counts do not require pseudocounts.

### Peak–RNA inference

Target linking uses the same **32,977 QC-pass nuclei** and library definitions as
the [middle-region follow-up](MIDDLE_LINK_REASSESSMENT.md). They represent two
donor-derived lines in undifferentiated and day-7 cultures, not four independent
donors. The target search includes every measured unique-symbol gene with a
GENCODE v48 transcript TSS within 500 kb of each peak midpoint: **94,296 pairs**.
Nearby targets may therefore be outside the original 221 genes.

Four partial-correlation models are used: technical depth only, plus COQ8A,
plus the fixed 184-gene myogenesis expression score, and plus both. If a target
belongs to that score, its own expression is removed before calculating the
score for its association. The fixed reference avoids redefining the score
according to favorable results. All 2,080 evaluable middle-region model/detection
results reproduce within 1e-7 under this implementation.

Primary eligibility requires RNA and ATAC detection in >5% of nuclei within
a library, at least 20 RNA-positive and ten peak-positive nuclei, and at least
three eligible libraries. This produces **27,720 pairs per model**. The already
specified 1% detection sensitivity produces **61,939 pairs per model**. Combined
correlations use within-library residualization and Fisher transformation;
library directions and leave-one-library-out estimates are retained.

### Count-based follow-up

The raw-count Poisson models use the same technical and biological covariates
and robust within-library standard errors as the middle follow-up. There are
**284 follow-up pairs** selected by the explicitly recorded rules:

- Nominal ATAC p<0.05, in either direction, with nominal peak–RNA p<0.05 in the
  technical-plus-COQ8A correlation model, at >5% detection.
- Nominal RNA links at previously discussed MYOD1, CSRP3, CAV3 or CACNA1H
  neighborhoods, regardless of the ATAC p.
- The previously discussed genes paired with their original candidate-region
  assignments, regardless of the RNA-link p, when evaluable at >5% detection.

Positive and negative links are retained. Count-model p values are follow-up
statistics on the same data, not independent replication. The primary display
uses nominal p values as an exploratory readout. Complete tables retain BH
values for the full 5,097-peak ATAC family and each complete correlation family.
`q_selected_followup` in the count tables is specific to the selected follow-ups;
it is not full-screen discovery FDR.

The sparse CSRP3 pair at chr11:19201752–19202603 is separately followed up at
the specified 1% detection sensitivity. It is not represented as passing the
5% primary gate. Its four raw-count models and bootstrap results are saved
under `results/unstratified/detection_1pct`.

For the prior focal genes MYOD1 and CAV3, the evaluable original candidate pair
with the smallest ATAC p is subjected to a 1,000-replicate, library-stratified
SCENT-style bootstrap. CSRP3 is separately bootstrapped at 1% detection. These
are selected follow-ups. The pinned SCENT source, basic p-value method, seed
208248 and finite p resolution of 0.002 are as documented in the middle follow-up.
This is a custom adaptation, not the complete adaptive SCENT pipeline.

## Results

Of the 5,097 peaks, **173 have nominal ATAC p<0.05**: 79 increase and 94 decrease.
These are exploratory counts, not claims of 173 FDR discoveries.

### Previously discussed loci

| Nearby/linked gene | GRCh38 interval | COQ8A high/low ATAC FC | Raw ATAC p | High/low RNA direction | Raw RNA p |
|---|---|---:|---:|---|---:|
| MYOD1 | chr11:17649919–17650798 | 1.630 | 0.02129 | Higher | 0.05104 |
| CSRP3 | chr11:19201752–19202603 | 2.429 | 0.03088 | Higher | 0.01210 |
| CAV3 | chr3:8733438–8733968 | 1.765 | 0.05325 | Higher | 0.04316 |

The ATAC fractions are 44/201 versus 27/201 for MYOD1, 17/201 versus 7/201 for
CSRP3, and 30/201 versus 17/201 for CAV3. The CSRP3 interval shown is the previously
discussed functional candidate, not the CSRP3-neighborhood peak with the smallest
ATAC p; that separate interval has FC=1.941 and p=0.01659 and remains in the full
screen. Figure gene labels identify explicit intervals, not averages over every
peak near that gene.

For the MYOD1 interval, the original candidate annotation places the peak midpoint
69,213 bp from the nearest annotated MYOD1 TSS. It is not described as the MYOD1
promoter. Among its evaluable nearby transcripts, MYOD1 has the strongest positive
fully adjusted partial correlation (r=0.0999, nominal p=9.45e-57). NUCB2 has a
smaller positive association, so the target assignment is not asserted to be exclusive.

| Target | RNA count ratio, technical model | RNA count ratio, COQ8A and score included | 95% robust interval | Raw count-model p |
|---|---:|---:|---:|---:|
| MYOD1 | 1.835 | 1.461 | 1.388–1.539 | 6.62e-47 |
| CAV3 | 1.383 | 1.185 | 1.130–1.243 | 3.26e-12 |

These ratios compare RNA counts in peak-detected versus peak-undetected nuclei;
they are not COQ8A high/low ATAC FCs. Both pairs are positive in all three eligible
libraries (25,097 nuclei). Their fully adjusted count ratios remain above one
when any eligible library is omitted. At 1% detection, partial correlations
are positive in all four libraries: MYOD1 r=0.0993 and CAV3 r=0.0451.

The fully adjusted stratified bootstrap gives MYOD1 RNA ratio 1.475
(95% interval 1.394–1.554) and CAV3 1.207 (1.149–1.267). Both bootstrap p values
reach the finite-resolution floor of 0.002; this is not an exact p smaller than 0.002.

CSRP3 has a positive 1%-gate partial correlation (r=0.0189, p=0.00283), but its
raw-count coupling attenuates from 2.584 to 1.003 when the expression score and
COQ8A are included (p=0.987). Thus its apparently favorable ATAC and RNA directions
are compatible with shared expression state; the different association methods
do not uniformly support a residual locus-specific link.
The fully adjusted CSRP3 bootstrap also gives a ratio near one (1.010;
95% interval 0.658–1.502; nominal p=0.954).

### Additional candidates

Five follow-up peak–gene pairs have concordant positive directions and nominal
p<0.05 for all three measured associations: COQ8A high/low ATAC, fully adjusted
peak–RNA counts and COQ8A high/low RNA. Their genes are **CREBRF, DMD, PFN2,
HSPB8 and RTN2**. This is a same-data exploratory intersection, not a new combined
significance test. Exact coordinates, effects and all three p values are in
`nominal_positive_chain.tsv`. Their biological roles must be evaluated separately;
passing this intersection does not establish promotion of differentiation.

MYOD1 and CAV3 remain interpretable focal candidates despite one high/low p value
being near 0.05. The threshold is not used to turn p=0.051 or p=0.053 into a
categorically contradictory biological conclusion.

### Independent MYOD1 occupancy annotation

All three displayed focal intervals overlap author-called MYOD1 ChIP peaks in
both human myoblast and 72-hour myotube GSE50413 samples (GSM1218849 and GSM1218850).
Intersection uses the existing unique hg19 mappings and the author hg19 ChIP
calls. This supports a myogenic regulatory context; occupancy neither establishes
the RNA target by itself nor proves COQ8A-dependent MYOD1 binding.

## Figure legend

**Unstratified candidate-region exploration.** A, all 5,097 peaks, with paired
accessibility differences in percentage points and nominal paired p values.
Orange and blue denote nominal p<0.05 in increasing and decreasing directions;
the dotted line marks 0.05. B, observed peak-detection fractions in 201 matched
COQ8A-high and low nuclei. The displayed FC and nominal p refer to the exact
intervals in the table above. Bars are descriptive fractions, not independent
donor means. C, within-library count-model peak–RNA ratios and robust 95%
intervals, combined across three eligible libraries. Blue adjusts technical
variables; orange also includes COQ8A and the target-excluded expression score.
The ordinate separates MYOD1 and CAV3. Original high/low RNA effects and full
model p values are supplied in the figure source tables. The exploratory figure
uses raw p; correction families and selection rules are specified above.

## Validation and reproduction

All 5,097 high/low peak count vectors and paired p values are recomputed from the
extracted nuclei and match the existing results. Existing overlapping RNA effects
are reproduced within 1e-7. Full validation counts are in `validation.json`.
The standard main figures and the previously reported middle results are unchanged.

After preparing the middle extraction cache and pinned SCENT source as documented:

```bash
python scripts/48_unstratified_links.py --h5 /data/GSE208248 --gtf /data/gencode.v48.annotation.gtf.gz --middle /work/middle_links --out /work/unstratified
python scripts/50_select_unstratified_followup.py --work /work/unstratified
Rscript scripts/49_unstratified_count_models.R /work/unstratified /work/middle_links/SCENTfunctions.R 1000
Rscript scripts/49_unstratified_count_models.R /work/unstratified /work/middle_links/SCENTfunctions.R 1000 /work/unstratified/detection_1pct 0.01
python scripts/51_report_unstratified_links.py --work /work/unstratified --external /work/middle_links/external
```

The four legacy libraries and score are reused for comparability, but no temporal
membership is used in this analysis. Render the figure directly from committed
tables with `python scripts/51_report_unstratified_links.py`.

The [middle follow-up](MIDDLE_LINK_REASSESSMENT.md) provides SCENT and independent
MYOD1 ChIP source links and software requirements. Private C2C12 passage/time
contrasts remain outside this repository and are not COQ8A perturbation contrasts.

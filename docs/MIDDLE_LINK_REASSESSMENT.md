# Middle-region peak–RNA association follow-up

## Question and analysis populations

Do any of the 34 externally defined middle-opening regions have an RNA target
association that persists after accounting for technical depth, COQ8A expression,
and a separate myogenesis expression score?

The region set and original COQ8A comparison are fixed. The original accessibility
analysis uses 201 within-library matched high/low pairs (COQ8A ≥3 versus 1 UMI,
TSS enrichment ≥3). Those effects and their p/q values are unchanged.
Target linking instead uses all 32,977 nuclei passing the existing TSS ≥3 joint
QC, including COQ8A-zero nuclei. This estimates peak–RNA coupling within each
culture library; it does not replace the high/low accessibility comparison.

| Library | Donor-derived line | Culture condition | QC-pass nuclei |
|---|---|---|---:|
| GSM6339597 | 1 | Undifferentiated | 7,535 |
| GSM6339599 | 1 | Day 7 | 8,216 |
| GSM6339601 | 2 | Undifferentiated | 9,346 |
| GSM6339603 | 2 | Day 7 | 7,880 |

The four libraries represent two donor-derived lines. No additional cell-type
annotation or purified subpopulation is imposed. Inference describes associations
among measured nuclei conditional on these libraries, rather than four independent
donor replications.

## Candidate universe and covariates

For every middle peak, all GENCODE v48 transcript TSSs within 500 kb of the peak
midpoint were considered, retaining measured protein-coding and noncoding RNA
features. Duplicate RNA gene symbols were excluded explicitly. The complete
universe contains **647 peak–gene pairs**.

The main detection filter requires both RNA and binary ATAC detection in >5%
of nuclei within a library, at least 20 RNA-positive nuclei and at least three
evaluable libraries for combination. Partial-correlation calculations additionally
require at least ten peak-positive nuclei; count models require at least 20
peak-negative nuclei. A 1% detection sensitivity is supplied for correlations.
These gates select evaluable pairs without using their association p values.
There are **107 evaluable pairs per primary model**.

Four models are reported together:

1. Technical depth covariates.
2. Technical covariates plus log1p COQ8A raw UMI.
3. Technical covariates plus the myogenesis expression score.
4. Technical covariates plus both terms.

The score is the mean log1p CP10k expression of 184 measured myogenesis genes
after excluding COQ8A and every candidate target gene from the 221-gene reference.
It is a covariate for shared expression state, not a newly inferred pseudotime.
Adjusting for it can remove a biological programme shared by peak and RNA; attenuation
does not establish that an unadjusted association was artifactual.

## Estimation and multiplicity

Partial correlations residualize binary peak detection and log1p CP10k RNA against
log1p total RNA UMI and log1p total open peaks, with the additional terms above.
Library estimates are combined by Fisher transformation, weighting by
`n - rank(covariate design) - 2`. Library signs and leave-one-library-out correlations
are retained. The original matched reference is reproduced using its original
binary COQ8A term and `n - 5` weighting: **38 non-missing reference correlations
reproduce within 1e-7**.

The complementary count analysis models raw gene RNA counts using Poisson
regression with binary peak detection. Technical terms are log1p RNA depth,
log1p open-peak depth and mitochondrial RNA percentage. HC0 sandwich standard
errors are calculated within libraries; log ratios are combined using inverse
variance weighting. The reported ratio is **expected RNA count in peak-detected
versus peak-undetected nuclei**, conditional on covariates. It is not the COQ8A
high/low ATAC fold change. Undetected peaks are an operational count category,
not proof of absolute physical closure.

For each fixed primary model, BH correction includes all 107 evaluable pairs.
The complete per-library and combined tables are retained. Detection sensitivities
and different adjustment models are reported separately, not searched for the
smallest p value. The source tables label the exact model and detection filter.

A SCENT-style raw-count bootstrap provides a follow-up for the previously
discussed LDB3 and TNNC1 pairs. It uses library intercepts, resamples nuclei
within each eligible library, and imports `basic_p` unchanged from SCENT commit
`e80b5ba6b445f972c7fe28fb41e24ef4f5b2e373`. Each fit has 1,000 valid bootstrap
replicates; seed 208248. The two-sided p-value resolution floor is 0.002.
Percentile 95% intervals and basic bootstrap p values are reported. This is a
custom fixed-budget, stratified adaptation, **not the full official adaptive
SCENT pipeline**. Follow-up bootstrap p values are nominal and do not substitute
for the 107-pair discovery q values.

## Findings

| Target | Technical-model RNA ratio | Technical-model q | RNA ratio with COQ8A and score | Adjusted-model q | Adjusted partial r |
|---|---:|---:|---:|---:|---:|
| LDB3 | 1.532 | 0.0000266 | 1.042 | 0.819 | 0.0044 |
| TNNC1 | 1.589 | 2.51e-14 | 1.107 | 0.00865 | 0.0287 |

Both pairs pass detection in three libraries (25,097 nuclei); GSM6339603 fails
the RNA detection gate. Thus the increased precision comes from more nuclei,
not additional donors. Adding COQ8A alone has little effect. The myogenesis
score accounts for most LDB3 coupling, whereas TNNC1 retains a small positive
association. For TNNC1, all three eligible libraries retain positive directions.

At the 1% detection sensitivity, TNNC1 is evaluable in all four libraries and
retains a positive adjusted partial correlation (r=0.0233, q=0.000957, four positive
directions). A second, sparse LDB3 region at chr10:86661433–86662314 becomes
evaluable at 1% and has r=0.0192, q=0.0110, with four positive directions.
Its original COQ8A high/low ATAC FC is 1.625 (p=0.383, q across 34 peaks=1.00).
It is retained as a separate sensitivity finding; it does not restore the
coupling of the originally highlighted LDB3 region. The 1% q values use their
own complete evaluable-pair family, supplied in the combined table.
That sensitivity family contains 413 pairs per model.

The pooled, library-stratified bootstrap gives fully adjusted RNA ratios of
1.055 for LDB3 (95% interval 0.928–1.194, p=0.430) and 1.075 for TNNC1
(1.004–1.154, p=0.042). Pooled regression and inverse-variance combination are
different estimators, so their ratios need not be identical.

The original COQ8A high/low accessibility effects remain:

| Region near | ATAC high/low FC | Paired p | BH q across 34 middle peaks |
|---|---:|---:|---:|
| LDB3 | 1.938 | 0.0315 | 0.536 |
| TNNC1 | 1.765 | 0.0660 | 0.748 |

MRAS and SVIL have additional positive peak–RNA links after adjustment. However,
their original COQ8A accessibility tests have q=1.00; SVIL has FC=0.909.
They therefore do not provide a stronger positive COQ8A-opening-to-RNA chain.
All results, including negative associations, remain in the complete tables.

Descriptive RNA comparisons across the author-labelled culture conditions are
provided separately. Mean normalized LDB3 and TNNC1 RNA is lower in day-7 than
undifferentiated libraries in both lines. Consequently, these target libraries
do not independently establish monotonic induction of either RNA over culture
time. Within-library peak–RNA coupling and culture-condition expression are
different comparisons; neither is a differentiation trajectory.

## Independent regulatory annotation

Author-called MYOD1 ChIP peaks from human MB135 myoblasts and 72-hour myotubes
were obtained from GSE50413 (GSM1218849 and GSM1218850). Existing unique hg19
mappings of the 34 regions were intersected with the author hg19 peak calls.
Any positive genomic overlap defines membership; overlap sizes and input hashes
are supplied. No COQ8A result determines this annotation.

Six middle regions overlap myoblast calls and 21 overlap myotube calls. The
LDB3 region overlaps both; the TNNC1 region overlaps neither. This supports a
myogenic regulatory annotation for the LDB3 interval, without identifying its
target solely from occupancy. Different numbers of called peaks do not by
themselves demonstrate quantitative changes in MYOD1 binding.

## Figure legend

**Middle-region peak–RNA associations.** A, LDB3 at chr10:86674717–86675570.
B, TNNC1 at chr3:52446715–52447619 (GRCh38). Points and bars represent combined
RNA count ratios and 95% robust intervals across the three eligible libraries.
Each added term is relative to the technical model; “Both” includes COQ8A and
the 184-gene myogenesis score. The vertical line denotes ratio 1. q values are
BH-adjusted across 107 evaluable peak–gene pairs within each displayed model.
Peak open/closed denotes ATAC detected/undetected. Source values are in
`results/link_followup/figure_source.tsv`.

## Reproduction

The existing processed QC, peak correspondence and temporal membership tables
are prerequisites. Python dependencies are the repository requirements; executed
versions are recorded in `results/link_followup/software_environment.json`. R 4.4.2
and Matrix 1.7-1 were used for the count/stratified-bootstrap script. SCENT source
provenance and its SHA256 are in `results/link_followup/scent_source.json`.

```bash
python scripts/44_expand_middle_links.py --h5 /data/GSE208248 --gtf /data/gencode.v48.annotation.gtf.gz --out /work/middle_links
python -c "import json,urllib.request,pathlib; m=json.load(open('results/link_followup/scent_source.json')); pathlib.Path('/work/middle_links/SCENTfunctions.R').write_bytes(urllib.request.urlopen(m['source_url']).read()); pathlib.Path('/work/middle_links/scent_source.json').write_text(json.dumps(m))"
Rscript scripts/45_middle_count_models.R /work/middle_links /work/middle_links/SCENTfunctions.R 1000
python scripts/46_middle_myod_annotation.py --out /work/middle_links/external
python scripts/47_report_middle_links.py --work /work/middle_links
```

Render the new figure from committed source tables, without large inputs:

```bash
python scripts/47_report_middle_links.py
```

This follow-up is additive; it does not change the main figure pipeline or its
original effects. Private passage/time RNA comparisons remain outside this
repository.

## References

- [SCENT primary study](https://www.nature.com/articles/s41588-024-01682-1)
- [SCENT source](https://github.com/immunogenomics/SCENT)
- [MYOD1 ChIP study GSE50413](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50413)
- [Myoblast MYOD1 ChIP GSM1218849](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1218849)
- [72-hour myotube MYOD1 ChIP GSM1218850](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1218850)

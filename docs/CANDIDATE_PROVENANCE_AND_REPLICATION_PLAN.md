# Candidate provenance and independent multiome assessment

The author object has now been inspected. Actual group counts, barcode checks,
matching feasibility and missing QC metrics are reported in
`GSE240061_FEASIBILITY_RESULTS.md`; this document retains the prior design and
candidate-selection audit.

## What selected CAV3

CAV3 was already a named exploratory candidate before the 1,777-peak and
401-peak enhancer scopes were constructed. The functional-dynamic figure code
(`39_prepare_figure_sources.py`) explicitly lists its focal coordinate.
The later count-model selection (`50_select_unstratified_followup.py`) retains
previously discussed MYOD1, CSRP3, CAV3 and CACNA1H gene pairs in addition to
nominal ATAC/link candidates. Its CAV3 bootstrap target is the eligible original
CAV3 pair with the smallest accessibility p. Scripts 51 and 58 explicitly display
that focal interval. Thus this was not an automatic 401-peak-to-CAV3 selection.
External enhancer and differentiation evidence was subsequently added to an
existing candidate. These two operations must not be represented as a single
prospectively defined discovery funnel.

## Uniform comparison of all candidates

Script 60 annotates every original neighborhood pair with the same evidence:
accessibility FC/p/q, two-source direction, fully adjusted RNA-link r/p/q,
COQ8A high/low RNA contrast, external GO function and independent C2C12 RNA
time-course results. No focal-gene exception is applied to the comparison.
The complete output includes 5,484 peak-gene neighborhood pairs for 5,097 peaks
and 221 genes. Missing RNA links remain missing rather than being replaced by
zero or by a more favorable model.

The focal CAV3 peak ranks 71st by two-sided ATAC p among 1,777 enhancer peaks
and 15th among 401 muscle-selective peaks (minimum ranks for ties; opening and
closing both included). Its accessibility p is 0.05325. It does not pass a nominal
p<0.05 accessibility gate, much less an accessibility FDR gate.

For descriptive comparison, the following cumulative intersections were applied
identically to every candidate; they are not a new combined statistical test:

| Cumulative evidence | 1,777 enhancer scope: distinct peaks | 401 selective scope: distinct peaks |
|---|---:|---:|
| Full scope | 1,777 | 401 |
| Positive ATAC FC and nominal ATAC p<0.05 | 36 | 6 |
| Also positive adjusted peak-RNA link, full link-family q<0.1 | 5 | 1 |
| Also positive COQ8A high/low RNA, nominal RNA p<0.05 | 2 | 0 |
| Also external differentiation RNA up at 48h, q<0.1 | 2 | 0 |

The two remaining 1,777-scope pairs involve HSPB8 and DMD. Their accessibility
tests still do not pass the full ATAC FDR threshold. Crossing several nominal
tests does not change that fact. Dmd's independent RNA baseline is zero in the
three GM replicates, so its very large modeled RNA FC is baseline-sensitive.
The fixed thresholds in this descriptive audit are not a claim that a p of
0.049 is biologically distinct from 0.051. Complete continuous evidence is kept
alongside the threshold intersections.

CAV3 remains biologically motivated by fusion annotations, a TSS-overlapping
interval, external MYOD1 occupancy, weak positive RNA coupling and differentiation
induction. That makes it a hypothesis-led candidate for replication. It does not
make it the statistically strongest discovery in the current dataset.

## Independent GSE240061 design

GSE208248 is the source of the current COQ8A high/low analyses. GSE240061 is an
independent same-nucleus RNA/ATAC study, with six participants and twelve biopsies.
Four participants exercise and two rest; each is sampled before and 3.5h after
the exercise/rest period. GEO assay accessions include technical subdivisions.
They are not additional donors. The author paper describes 37,154 jointly
QC-passing nuclei in fourteen cell types.

The GEO source/tissue fields say PBMC, inconsistent with the paper and series
summary identifying vastus lateralis muscle. The verbatim GEO fields are retained
in the manifest; tissue interpretation follows the paper and requires confirmation
against the actual object annotations. The article's nucleosome-signal cutoff
text is also ambiguous and is not copied as a new QC threshold.

### What transfers

- Biological gene universe and external annotations.
- Same-nucleus raw RNA grouping and raw ATAC measurement.
- Within-sample depth/QC matching, plus continuous-expression sensitivity.
- Full candidate reporting, both directions, and separate RNA-link evaluation.

The old 401 or 1,777 feature counts cannot simply be asserted for the new matrix.
Genome build, coordinates, peak overlap and observability must first be checked.
An overlapping author peak is a region-level correspondence; exact old intervals
require recounting fragments. No overlap is unmeasured, not proof of closure.

### Primary population and statistical units

Preserve the author's cell-type labels and check their marker support. MuSC/
satellite cells address the closest cell context to myoblast differentiation.
Mature myonuclear subtypes are analyzed separately and address maintenance or
maturation context; they are not substituted for myoblasts after inspecting p.
Do not pool these types to create a COQ8A association.

Start with the six baseline biopsies. Within each participant and cell type,
inventory raw COQ8A detection and the existing >=3 versus exactly-1 comparison.
Nuclei with zero observed RNA are not automatically called biological low
expression. Inventory >=2 versus 1 and nonzero-expression quantiles as declared
sensitivities; quantiles that collapse because of discrete UMI ties are reported
as such. Do not force a fixed 25% split by arbitrarily breaking expression ties.

For loci with adequate count information, donor-by-expression-group pseudobulk
ATAC counts support a donor-blocked negative-binomial comparison, e.g.
`~ donor + COQ8A_group`, with library normalization and dispersion estimation
over a suitable broad set of measured peaks. The external biological scope
defines the reported candidate family. Matching handles RNA depth, ATAC depth
and available QC metrics before aggregation. Existing matched binary accessibility
fractions can be retained as an interpretable effect summary; they are a different
estimand from normalized count FC. Report usable donors and group counts for each
cell type. A donor missing a comparison group cannot supply its paired contrast;
no arbitrary minimum of 20 MuSC nuclei is introduced.

Post samples provide a secondary within-person consistency analysis. They are
not six additional independent donors and not an independent validation of the
baseline result. A combined repeated-measure analysis must retain donor dependence
and exercise/time structure; otherwise analyze baseline and post separately and
compare source-specific effects. Pseudobulk or replicate-aware mixed models are
supported by published differential-accessibility benchmarks.

### A function-based secondary family fixed before new outcomes

The external GO myoblast-fusion annotation contains twelve genes in the existing
221-gene universe: ADAM12, CACNA1H, CAV3, CDON, ITGB1, KCNH1, MAPK14, MYH9, MYOD1,
MYOG, NEO1 and NOS1. All are retained, regardless of prior COQ8A p values.
The existing shared peaks with midpoint within 2 kb of an annotated transcript
TSS yield 31 distinct promoter candidates. Coordinates and gene assignments are
frozen in `external_fusion_promoter_candidates.tsv` before inspecting GSE240061
outcomes. This secondary family gives CAV3 a defined biological context without
requiring CAV3 to win. GO membership itself does not specify that increased
expression always promotes fusion.

This is a proposed independent replication design, not retroactive relabeling
of the original GSE208248 analysis. The broad enhancer assessment and this
secondary fusion-promoter family retain their own explicitly declared scopes.
No significant new COQ8A result is assumed in advance.

## Reproduction and sources

```bash
python scripts/60_audit_candidate_provenance.py
python scripts/62_fetch_replication_design.py --out results/replication_GSE240061
Rscript scripts/61_inspect_replication_multiome.R /data/GSE240061_integrated11723.rds /path/to/repository /work/GSE240061_inventory
```

The GEO download has an outer gzip layer around a gzip-compressed RDS. Remove
the outer layer before using the inspection command. The inspection preserves
author labels, verifies RNA/ATAC barcode identity, inventories raw counts and
exports only the candidate matrices for further analysis. Downloaded large
inputs and per-nucleus intermediates stay outside the release repository.

- [GSE240061](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE240061)
- [Rubenstein et al., Genome Research 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12212352/)
- [Best practices for differential accessibility, Nature Communications 2024](https://www.nature.com/articles/s41467-024-53089-5)

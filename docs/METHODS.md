# Methods: COQ8A-associated accessibility at myogenic regions

## 1. Objective and design

We ask whether COQ8A-high nuclei exhibit altered accessibility at regions that
change during human myoblast differentiation, and whether nearby genes show
concordant RNA associations. This is an exploratory secondary analysis of
processed public data. It contains no COQ8A perturbation experiment.

The analysis uses three distinct information sources:

1. MSigDB Hallmark/Reactome and GO annotations define biological scope.
2. GSE109828 defines region dynamics using experimentally recorded times.
3. GSE208248 provides paired RNA/ATAC measurements for COQ8A association tests.

Private C2C12 passage/time RNA results provide a separate biological comparison.
They are not used to define public temporal classes and are not COQ8A OE/KD
contrasts. Their sample-level data are not distributed in this repository.

The analysis is retrospective and exploratory. External memberships are
computed without target COQ8A effect values; choices developed during
exploration are not described as prospective preregistration.

## 2. Biological universe: union versus intersection

The measured union includes 198 Hallmark Myogenesis genes and 29 Reactome
Myogenesis genes, with six shared genes: 198 + 29 - 6 = 221. Hallmark has 200
listed genes; DENND2B and MYL11 are absent from the target RNA feature matrix.
Exact frozen membership is `reference/myogenesis_221_gene_sources.tsv`.

The intersection is MAPK12, MEF2A, MEF2C, MEF2D, MYF6 and MYOG. It has 146 common
candidate peaks, including eight early-opening and six closing dynamic peaks,
with no middle/late peaks under the frozen external definitions. These counts
are an inventory, not a newly fitted intersection association analysis.

CSRP3 and CAV3 are Hallmark-only members; MYOD1 is a Reactome-only member.
Consequently the intersection does not include them as nominated genes.
Union and intersection represent different biological scopes, not different
degrees of statistical validity. Union is appropriate for broad myogenesis;
intersection asks about functions jointly represented in both particular lists.
Changing to the intersection requires distinct analysis outputs and cannot
inherit the union's FC or q values.

## 3. Target multiome: GSE208248

[GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208248) and the
[source paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10123345/) describe
human skeletal-muscle stem cell lines RH-hMuSC-1 and RH-hMuSC-2, each sampled
before and after differentiation. The paper identifies their origins as a
56-year-old male's tibialis anterior and a 62-year-old male's anterior thigh,
respectively. They are not a young-versus-old comparison.

GEO describes expansion for six passages and seven days of differentiation,
with 5% FBS replaced by 2% horse serum in the culture medium. The deposited
sample protocols contain the complete growth-factor/inhibitor medium recipe.
We use the four untreated culture/differentiation multiome libraries, not the
paper's separate cancer-conditioned-medium RNA experiments.

The 10x Multiome assay pairs RNA and ATAC through the same nucleus barcode.
The eight GEO assay accessions represent four paired multiome libraries, not
eight biological samples. Barcodes are scoped by library; identical barcode
strings from different libraries do not identify the same nucleus.

| Line | Condition | H5/RNA accession | ATAC fragment accession | Input filtered barcodes | Main QC-pass nuclei | Main matched pairs |
|---|---|---|---|---:|---:|---:|
| 1 | Undifferentiated | GSM6339597 | GSM6339598 | 8,346 | 7,535 | 46 |
| 1 | Differentiated, day 7 | GSM6339599 | GSM6339600 | 9,221 | 8,216 | 74 |
| 2 | Undifferentiated | GSM6339601 | GSM6339602 | 10,659 | 9,346 | 58 |
| 2 | Differentiated, day 7 | GSM6339603 | GSM6339604 | 13,396 | 7,880 | 23 |
| Total | | | | 41,622 | 32,977 | 201 |

Biological source N is two donor-derived lines, library N is four, and the main
association sample comprises 201 high/low pairs (402 distinct nuclei).

### Starting files and preprocessing boundary

We start from author-deposited `filtered_feature_bc_matrix.h5`, barcode metrics
and indexed ATAC fragment files. We do not rerun alignment or peak calling from
FASTQ. The target reference genome is GRCh38; the deposited metadata identify
the 10x ARC GRCh38-2020-A-2.0.0 reference. GEO uses inconsistent pipeline naming
(`cellranger-atac v2.0` and ARC outputs); we therefore identify the actual
processed inputs/checksums rather than claim a independently verified aligner
execution. `docs/input_sha256.csv` records target source-file hashes.

No Harmony correction, new clustering, trajectory inference or transfer of HMA
annotations enters the current comparisons. Association values derive from
counts and matching within each original library/condition.

## 4. Joint nucleus QC and COQ8A grouping

Apply the following fixed QC to the paired nucleus:

- RNA UMI >=500 and number of detected ATAC peaks >=500.
- Remove nuclei above the within-library 95th percentile of either depth;
  percentiles are calculated after both lower-depth requirements.
- FRiP >=0.25, using 10x peak-region fragments / ATAC fragments.
- Nucleosome signal <4, defined as fragments with lengths [147,294) bp divided
  by fragments shorter than 147 bp.
- ENCODE hg38 v2 blacklist fragment fraction <0.05.
- TSS enrichment >=3 in the main analysis; >=2 in the sensitivity analysis.

The custom TSS implementation uses GENCODE v48 protein-coding gene TSSs.
Central insertion density is measured within +/-500 bp (1,001 bp); background
is the combined 200 bp at distances 901-1,000 bp on either side. A nucleus
with no flank insertions receives the pooled flank-density estimate for its
library. Fragment positions in this QC script are start and end-1. This exact
custom implementation, rather than a named package default, defines the
archived QC scores. See `00_fragment_qc_reference.py`.

Mitochondrial RNA fraction, detected-gene count and Scrublet tail exclusions
are evaluated as additional sensitivities in `25_temporal_rna_qc.py`. They
are not baseline filters; comprehensive doublet removal is not part of this
baseline.

COQ8A-high is raw RNA UMI >=3; low is exactly 1. Zero-count nuclei are outside
this specific comparison. The high >=2 definition is retained as a sensitivity.
These are raw UMI thresholds, not normalized RNA values or expression quantiles.

Match high to low 1:1 without replacement within each library, requiring
absolute differences <=0.10 in both log1p(total RNA UMI) and log1p(open ATAC
peak count). The implementation searches up to 100 nearest low-group nuclei
in Euclidean log-depth space, processes high nuclei from largest nearest-low
distance downward, and accepts the first unused candidate within both calipers.
The paired barcode table is the reproducibility record.

| TSS threshold | High versus low COQ8A UMI | Pairs |
|---|---|---:|
| >=3 | >=3 versus 1 | 201 |
| >=2 | >=3 versus 1 | 212 |
| >=3 | >=2 versus 1 | 958 |
| >=2 | >=2 versus 1 | 1,020 |

All locus, RNA and functional-module comparisons retain their own stated
setting consistently. Matching never pairs an undifferentiated nucleus to a
differentiated nucleus, or one line to the other line.

## 5. Origin of the 5,097 common candidate peaks

These peaks come from GSE208248, not from MSigDB or GSE109828.

1. Extract author-called peak coordinates from each target H5 feature list.
2. Anchor on GSM6339597. Match regions from each other library one-to-one with
   at least 50% reciprocal interval overlap, using `03_define_regions.py`.
3. Retain anchor intervals on canonical chromosomes, 200-2,000 bp wide and
   outside the ENCODE blacklist.
4. Find intervals whose midpoint is within 100 kb of any protein-coding
   transcript TSS for a gene in the fixed 221-gene universe (GENCODE v48).
5. Retain regions represented in all four library peak catalogs, and count
   each anchor interval once even if it is near multiple genes.

The resulting universe contains 5,097 distinct regions. "Common" means a
corresponding interval exists in each library's catalog. It does not mean every
nucleus is accessible there, or that any COQ8A difference is significant.
Genomic proximity nominates potential targets; it does not establish regulatory
links. The broader >=3-library candidate catalog contains 7,132 distinct peaks.
The temporal and functional analyses use the all-four-library subset.

## 6. External timing reference: GSE109828

[Pliner et al., Molecular Cell 2018](https://doi.org/10.1016/j.molcel.2018.06.044)
used human skeletal muscle myoblasts (HSMM), Lonza CC-2580, lot 257130, derived
from a healthy 17-year-old female quadriceps biopsy. Cells were used within
five passages of purchase and grown in SKGM-2. At approximately 80-90%
confluence, differentiation was induced by switching to alpha-MEM with 2%
horse serum. The selected datasets contain no COQ8A perturbation.

We use only the two sci-ATAC experiments listed in the
[author's data inventory](https://cole-trapnell-lab.github.io/cicero-release/data/).
Other GEO entries, including bulk ATAC/CRISPR experiments, are not included.
The selected sci-ATAC data are not same-nucleus RNA/ATAC multiome measurements.

| Experiment | GEO accession | Hours after serum switch | Input cells | Retained cells in our primary QC |
|---|---|---:|---:|---:|
| 1 | GSM2970930 | 0 | 1,674 | 808 |
| 1 | GSM2970930 | 24 | 1,574 | 995 |
| 1 | GSM2970930 | 48 | 1,762 | 1,049 |
| 1 | GSM2970930 | 72 | 1,501 | 904 |
| 2 | GSM2970931 | 0 | 2,971 | 1,645 |
| 2 | GSM2970931 | 72 | 3,885 | 2,866 |

There are 13,367 input cells and 8,267 retained cells, two experiments from the
same source lot, and one donor source. Experiment 2 checks the 0-to-72-hour
direction; it cannot independently estimate the 24/48-hour onset.

### External processing and coordinate alignment

Use the author's sparse counts and barcode-to-time metadata, originally mapped
to hg19. GEO documents Bowtie2 v2.2.3 mapping, per-cell duplicate removal and
MACS2 v2.1.1 peak calling. These upstream operations are not rerun here.

Require >=1,000 accessible sites per cell and accessibility at a promoter of
at least one of MYOG, MYOD1, DMD, TNNT1, MYH1, MYH3 or TPM2. Our promoter
definition uses all annotated transcript TSSs +/-2 kb; +/-1 kb is a sensitivity.
The resulting cell counts are our reconstruction, not the paper's published
pseudotime cell subset.

Lift target hg38 intervals to hg19 with the UCSC chain, requiring unique,
concordant endpoint/midpoint mapping and <=5% length change. 5,095 of 5,097
pass. Select the external interval with greatest minimum reciprocal overlap;
ties use overlap length and coordinate. Require >=50% overlap in both
directions; 25% is a mapping sensitivity. Require external accessibility in
>=1% and >=10 retained cells. This yields 3,324 eligible mapped target peaks.

Within each experiment, divide retained cells into five pooled open-site-depth
quantiles. At each time, calculate binary accessibility rates within depth
strata and standardize to the pooled stratum weights, retaining strata present
at every time. Rates remain unsmoothed. Jeffreys-smoothed binomial rates are
used only to stabilize variance estimates at boundaries.

## 7. Defining early, middle, late and closing regions

In experiment 1, a 3-df Wald test compares accessibility at 24/48/72 hours
against 0 hours, accounting for the shared baseline covariance. BH adjustment
is across 4,580 eligible external intervals among the 4,710 overlapping
intervals prepared for this analysis, before the final 50% target mapping cut.

An opening region requires external q <0.05, increased 72-versus-0-hour
accessibility in both experiments, and sufficient detection in experiment 2.
A closing region requires corresponding decreases in both experiments.

For opening regions, calculate the first time that accessibility crosses half
the observed rise from baseline to the maximum of 24/48/72-hour values, using
linear interpolation between sampled times. This is actual-time interpolation,
not pseudotime and not a claim of exact molecular activation time.

| Temporal class | Definition | Target peaks |
|---|---|---:|
| Early | Half-rise <=24 hours | 191 |
| Middle | Half-rise >24 and <=48 hours | 34 |
| Late | Half-rise >48 hours | 4 |
| Closing | Significant dynamic profile, decreased endpoints in both experiments | 181 |
| All dynamic | All four nonoverlapping peak classes | 410 |

These labels use no target COQ8A effects. Opening/closing describe the external
differentiation experiment; either class can show increased or decreased
accessibility in the target COQ8A-high comparison. Genes can have peaks in more
than one class. Memberships are saved before target-effect testing.

## 8. Functional subset definition

Intersect the 221 genes with the complete external GO sets for myoblast
differentiation, myoblast fusion, positive regulation of muscle cell
differentiation and negative regulation of muscle cell differentiation. Retain
the already-defined 410 dynamic peaks near those genes. This combines two
external criteria: biological function and temporal accessibility dynamics.

The four sets contain 17, 12, 14 and four genes, respectively; their union is
30 genes, of which 22 have retained dynamic regions. There are 54 distinct
peaks: 24 early, five middle, zero late and 25 closing. All four functional
sets are retained; none is selected because of its observed COQ8A P value.
Full memberships, source URLs and hashes are in `results/functional/`.

This step does not require MRF motifs or MYOD binding. Those are different
questions. The CSRP3, CAV3 and MYOD1 loci discussed subsequently are early
regions under the external timing definition.

## 9. Target association tests and adjustment families

For a peak, binary accessibility is 1 if the author count matrix contains a
positive entry and 0 otherwise. Open fraction is open nuclei / evaluated
nuclei. FC is the high-group open fraction divided by the low-group fraction;
absolute differences are also reported in percentage points.

For each matched pair, record whether only high or only low is open. The
two-sided exact binomial test of these discordant counts against probability
0.5 is the paired single-peak test. Unchanged pairs do not enter the discordant
count. All preselected dynamic peaks remain in the 410/54-peak screens; no
nominal-significance or favorable-direction filter is applied before BH.

Module scores are the fraction of distinct member peaks accessible per
nucleus. A paired t test evaluates the mean high-minus-low score, with a 95%
t interval. Module FC is the ratio of group mean scores. The gene-averaged programme score in the broad 221-gene screen instead
weights genes equally; these endpoints are reported separately.

| Endpoint | BH family in a fixed QC/count setting |
|---|---|
| Dynamic individual peaks | All 410 opening/closing peaks |
| Functional individual peaks | All 54 distinct functional dynamic peaks |
| Functional modules | Eight function-by-opening/closing modules |
| Temporal programmes | Ten programme tests per fixed external/QC/count setting |
| Original gene RNA comparisons | All 221 genes |
| Expanded dynamic peak-RNA links | All evaluable links; 1,022 in main setting |

Full temporal programme analyses and phase-specific peak families remain
available with their own explicit labels. Sensitivity settings are not pooled
into one BH family. These are nucleus-level tests; separate library directions
are reported, and four libraries do not constitute four independent donors.

RNA uses log1p(10,000 * raw count / total RNA UMI). Test the paired difference
under the same nucleus IDs, gate and COQ8A definition. A difference on this
log1p scale is not a raw-expression FC and is not labelled as one.

## 10. Peak-to-gene associations and candidate nomination

Test potential cis genes with a protein-coding transcript TSS within 500 kb of
each dynamic peak midpoint. For each library, residualize binary ATAC and
normalized RNA on intercept, COQ8A group, log RNA depth and log ATAC depth.
Correlate residuals. Eligibility requires >=10 accessible nuclei and >=20
RNA-positive nuclei per library; at least three eligible libraries are needed.
The implemented Fisher-z combination weights libraries by n_nuclei - 5.
`30_middle_cis_links.py --region-family opening-closing` evaluates this expanded
screen. This custom partial-correlation method is not Signac LinkPeaks and does
not model GC-matched background peaks.

Nearby position, high/low RNA difference and within-nucleus peak-RNA
correlation are separately reported. None is silently substituted for another.
Complete outputs include all directions and all eligible targets. CSRP3/CAV3/
MYOD1 are exploratory follow-up priorities based on these measurements,
independent passage/time RNA patterns and function; no implemented automatic
filter has been shown to leave exactly those three genes.

### Existing private C2C12 RNA comparison

The supplied expression workbook contains P11, P22 and P33 at differentiation
days 0, 1, 2 and 6, with three sample columns per passage/day: 36 samples in
total. Existing contrast outputs are joined by case-insensitive human/mouse
gene symbols, not by mapping human ATAC peaks onto the mouse genome.
The archived RNA analysis used an expression filter of normalized count >=10
in at least three samples and limma on log2(DESeq-normalized count + 1).
We retain its contrast estimates and original transcriptome-wide FDRs without
refitting or recalculating FDR only for nominated genes. Model FC is 2^log2FC.
P11 D2/D0, P22/P11 at D2 and P33/P11 at D2 address different comparisons and are
not substituted for each other. Sample-column replicates must be described
according to the original experimental records; they are not human donors.
Reproduction of this private comparison requires the owner's existing RNA
results and `35_crosscheck_existing_rna.py`, with output directed outside the
public repository. Raw RNA values and private comparison figures are not
published here.

## 11. Visualization and source data

Figure 1 combines the selection workflow, external temporal heatmap, mean
time profiles and target COQ8A module comparisons. Figure 2 shows all 54
functional regions as peak-centered insertion heatmaps, per-library open
fraction differences and pooled fold changes. Figure 3 shows four candidate
loci with fragment tracks, gene models, library-level accessibility and RNA.

ATAC profiles use deduplicated deposited fragments, counting each start/end
once without read-multiplicity weighting or an additional Tn5 offset. Counts
are normalized per 100 nuclei per 100 bp; low and high share the same scale
and row order. Gaussian smoothing (standard deviation 100 bp) is for display
only. Peak matrices span midpoint +/-2 kb. Locus windows span the tested peak
and nearest annotated transcript TSS, with 5-kb flanks. Transcript selection
uses distance followed by transcript ID; exon models are clipped to the window.

The 221-gene TSS supplement spans +/-5 kb and 100-bp bins, with strand reversal,
a shared absolute scale and row order based on pooled signal. TSS-window
plots do not define candidate peak membership. Exact panel units, color-scale
limits, ordering and correction families are in [Figure legends](FIGURE_LEGENDS.md).
Source data are indexed in [Data availability](DATA_AVAILABILITY.md).

## 12. Reproduction from deposited processed inputs

Use Python 3.13 and the versions in `requirements.txt`. From the repository root:

```bash
python -m pip install -r requirements.txt
python run_paper.py --mode figures
```

For complete numerical analysis from author-deposited processed H5 matrices,
barcode metrics, versioned references and GENCODE v48:

```bash
python run_paper.py --mode analysis --h5-root /data/GSE208248 --gtf /data/gencode.v48.annotation.gtf.gz --external-raw /data/GSE109828 --external-work /work/GSE109828
```

`19_fetch_temporal_reference.py` downloads the external data and liftOver chain.
The pipeline computes both external experiments, temporal definitions, target
associations, functional membership/tests, same-setting links and figures.
`--reuse-external` skips external preparation only when the prepared files
already exist with the intended source data and parameters.

The default analysis uses versioned per-nucleus fragment-QC tables and cached
figure profiles. Barcode/window hashes must agree with current figure sources;
validation fails if the figure profiles belong to different matched nuclei.
To regenerate all figure profiles, install `pysam` and add
`--fragments /data/fragments` to the analysis command. On Windows, indexed
fragment processing can be run with Python/pysam under WSL. Frozen source
files allow figure reproduction without fragment processing.

### Regenerating fragment QC

After downloading indexed fragments, regenerate reference inputs and the three
QC measurements for each ATAC library:

```bash
python -m pip install -r requirements_reference.txt
python scripts/01_extract_nuclei.py --h5-root /data/GSE208248 --out results/tables/nuclei.tsv.gz
python scripts/00_prepare_reference.py --nuclei results/tables/nuclei.tsv.gz --gtf /data/gencode.v48.annotation.gtf.gz --out reference
python scripts/00_fragment_qc_reference.py GSM6339598 tss --reference-root reference --fragment-root /data/fragments --out-root reference/fragment_qc
python scripts/00_fragment_qc_reference.py GSM6339598 nuc --reference-root reference --fragment-root /data/fragments --out-root reference/fragment_qc
python scripts/00_fragment_qc_reference.py GSM6339598 blacklist --reference-root reference --fragment-root /data/fragments --out-root reference/fragment_qc
```

Repeat the last three commands for GSM6339600, GSM6339602 and GSM6339604,
then execute the analysis command. All threshold values are defined in the
QC section; using a different TSS implementation requires a separate analysis.

### Optional RNA-quality sensitivity

```bash
python scripts/scrublet_provenance.py --h5-root /data/GSE208248 --barcodes reference/gse208248_qc_barcodes.tsv --out reference/rna_scrublet_qc.tsv.gz
python scripts/25_temporal_rna_qc.py --h5-root /data/GSE208248 --tables results/tables --temporal results/temporal --scrublet reference/rna_scrublet_qc.tsv.gz
```

The [code map](CODE_AVAILABILITY.md) identifies the responsible script for each
step. Source hashes are recorded in `docs/input_sha256.csv`,
`results/temporal/download_manifest.json` and
`results/functional/analysis_manifest.json`. Membership hashes are checked
before functional testing. Validation recomputes single-peak p/q and temporal
programme statistics from stored paired observations.

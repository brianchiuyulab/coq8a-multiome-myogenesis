# Computational and figure audit

## Scope and outcome

The fixed Day-0 comparison was audited against the deposited GSE208248 H5
matrices and GENCODE v48 annotations. All 33 raw-data audit checks passed.
All 565 peak count contrasts, FCs, exact paired p values and BH565 values were
reproduced. Repeating 2,000,000 permutations with seed 20261134 reproduced
every exceedance count for all 25 opening and 25 closing regional tests.
The six primary neighborhoods and their statistics are unchanged.

The current executable entry point, `run_fixed_day0.py --mode analysis`, was
also run through matrix extraction, matching, complete statistical testing,
local RNA models and figure generation. It now recomputes the fixed regional
tests directly rather than requiring a sensitivity-grid run.

## Checks performed

| Component | Check | Outcome |
|---|---|---|
| Biological scope | Rebuild GO differentiation/fusion intersection with the 221-gene universe | 25 genes |
| Peak scope | Rebuild transcript-TSS neighborhoods with 0-based coordinates | 5,097 broad and 565 restricted peaks; membership unchanged |
| Peak mapping | Canonical anchor coordinates, widths, blacklist flags, one-to-one reciprocal overlap | Consistent with documented rules |
| RNA candidates | Independently parse all GENCODE transcript TSSs within 500 kb and intersect unique RNA features | All 170 pairs reproduced |
| Joint QC | Reapply RNA/ATAC depth floors, within-library upper caps, FRiP, nucleosome, blacklist and TSS gates | 7,535 + 9,346 eligible Day-0 nuclei |
| Raw depths and counts | Independently extract all 565 peaks and all cached RNA features from both H5 files | Exact agreement |
| Grouping and matching | Check raw COQ8A UMI, library, barcode uniqueness and both depth calipers; rerun matcher | Same 226 + 301 pairs |
| Peak statistics | Recalculate all discordant-pair counts, FC, exact p and BH565 | Agreement |
| Regional statistics | Check all maxima and BH25, then repeat all 2 million swaps | All 50 exceedance counts identical |
| RNA models | Recalculate every source/model/population using QR residual projection independently of the original least-squares implementation | All partial correlations agree |
| RNA integration | Check Fisher-z aggregation, degrees of freedom, detection support and four BH170 families | Agreement |
| RNA contrasts | Recalculate all 170 paired RNA contrasts from raw counts | FC, p and q agree |
| QC score formulas | Recalculate archived TSS and nucleosome scores from their insertion/fragment numerators and denominators | Agreement |
| Fragment source tables | Check all 274,412 rows for unique keys and per-source/per-bin normalization | Agreement |
| Figures | Render all seven figures and inspect labels, scales, legends and page boundaries | Revised and checked |

The audit validates existing fragment-QC measurements and their formulas; it
does not claim that this audit realigned FASTQ reads or recounted every
genome-wide fragment for QC. The original QC regeneration code is included.

## Corrections made

1. **Coordinate conversion.** Candidate-neighborhood annotation previously
   compared 1-based GTF transcript starts with 0-based peak coordinates. The
   conversion is now explicit. Every peak–gene membership is unchanged;
   distance annotations change by at most one base.
2. **Fragment display coordinates.** The current plotting extractor used
   `end` for the right fragment endpoint, whereas the QC implementation used
   `end-1`. It now follows the same BED half-open convention. Minus-strand
   TSS windows use the corresponding orientation, and profiles are plotted
   at bin centers. Positional profiles were re-extracted from the original
   indexed fragments. This changes profile-bin allocations, not H5-based
   accessibility counts, high/low membership, FC, p or q. See the
   [10x fragment specification](https://www.10xgenomics.com/support/software/cell-ranger-atac/latest/analysis/outputs/fragments-file).
3. **Analysis entry point.** The homepage and code-availability page previously
   mixed the current comparison with the older temporal runner. The current
   workflow has one entry point, one source-data directory and one figure
   directory. Historical routes remain in a separate index.
4. **Figure interpretation.** The RNA panel now displays eligible and matched
   populations side by side. Filled/hollow points distinguish RNA q<0.05
   from q>=0.05; link significance and unavailable estimates are explicit.
   Primary and secondary accessibility loci are in separate figures.
   Axis spacing, type sizes, color scales and legend placement were revised.

One nonsignificant IGF1 closing neighborhood contains tied maximum-statistic
peaks. Annotation row order can change which tied peak is listed first. Its
regional statistic, p and q are unchanged; none of the selected loci are
affected. The audit accepts a tied representative only if its statistic is
equal to the reported maximum.

## Parameters and interpretation

The audit found no mismatch between the fixed definitions and their execution:
Day 0, TSS>=3, COQ8A>=2 versus 0 raw UMI, within-source matching and caliper
0.30 apply to every candidate. The ±100-kb candidate neighborhood and ±500-kb
RNA search window are explicit analysis choices. They are not stated to be
identical to the original authors' complete analysis pipeline.

Statistical outputs describe paired nuclei within two biological sources.
Opening and closing use separate BH25 families; individual peaks use BH565.
The q values apply within this fixed comparison. Parameter-search results and
state sensitivity are documented in the separate sensitivity report. Local
neighborhood naming is not substituted for RNA target identification.

## Reproduction records

- `scripts/97_audit_fixed_day0.py`: raw H5, annotation and statistical audit.
- `scripts/99_validate_fixed_tables.py`: published-table and normalization checks.
- `results/fixed_day0_D206/audit/audit_report.json`: all 33 checks and runtime versions.
- `results/fixed_day0_D206/audit/regional_permutation_rerun.tsv`: repeated permutation counts.
- `results/fixed_day0_D206/audit/published_table_validation.json`: source-table checks.
- `results/fixed_day0_D206/audit/source_table_sha256.tsv`: figure and model source-table hashes.

Text-table hashes normalize CRLF to LF so that Windows and Linux checkouts
produce the same digest. Gzip tables are hashed as their exact stored bytes;
the manifest records the convention for each file.

Figures are provided as 400-dpi PNG plus editable PDF/SVG. The
[Nature figure guide](https://research-figure-guide.nature.com/figures/building-and-exporting-figure-panels/)
informed the typography and export format. Figure legends specify the signal
units, statistical families and source counts.

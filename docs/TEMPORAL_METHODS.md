# External temporal definition of myogenic accessibility

## Design and data

The question is whether COQ8A-high nuclei have greater accessibility at regions that open during differentiation. The fixed universe remains the existing **221 measured Hallmark/Reactome myogenesis genes**, with 5,097 distinct candidate peaks present in all four target libraries. Regions are assigned temporal classes using an independent dataset before testing COQ8A associations.

[Pliner et al., Molecular Cell 2018](https://doi.org/10.1016/j.molcel.2018.06.044) measured human skeletal muscle myoblast differentiation by sci-ATAC. Its [author data inventory](https://cole-trapnell-lab.github.io/cicero-release/data/) identifies GSE109828. GSM2970930 provides 0, 24, 48 and 72 h; GSM2970931 provides 0 and 72 h. These experiments use the same source cell lot, not independent human donors. The second experiment checks endpoint direction, not onset timing.

The complete processed count matrices contain 59,283,741 and 67,971,331 nonzero region/cell records. All barcodes map to metadata; no duplicate cell/region entries occur among retained overlapping regions. Exact URLs, byte counts and SHA256 hashes are in `results/temporal/download_manifest.json`.

This is a reconstruction of **actual-time accessibility profiles**, not a reproduction of the paper's published cell assignments, Cicero links, pseudotime or changepoints. Early chromatin opening may precede later RNA expression; genes are not excluded merely because their functions concern mature muscle.

## External preprocessing

1. Require >=1,000 accessible sites per cell and accessibility at at least one promoter of MYOG, MYOD1, DMD, TNNT1, MYH1, MYH3 or TPM2, following the study's marker-filtering rationale. Here promoters use **all GENCODE v48 transcript TSSs +/-2 kb**, with +/-1 kb tested in parallel. This explicit annotation/geometry is our implementation.
2. With +/-2 kb, experiment 1 retains 808, 995, 1,049 and 904 cells at 0/24/48/72 h; experiment 2 retains 1,645 and 2,866 cells at 0/72 h. With +/-1 kb the counts are 797, 985, 1,034, 891 and 1,604, 2,816.
3. Convert target hg38 intervals to hg19 using the UCSC hg38ToHg19 chain. Both endpoints and midpoint must map uniquely with concordant chromosome/strand and <=5% length change. **5,095/5,097** intervals pass.
4. Assign each target peak the external interval with greatest minimum reciprocal overlap; break ties by overlap length and coordinate. Require >=50% overlap in both intervals; >=25% is a mapping sensitivity. No COQ8A effect is used. Timing quantiles use unique external intervals when multiple target peaks share a match.
5. Require external detection in >=1% of retained cells and >=10 cells in the relevant experiment. The primary definition yields **3,324** eligible mapped target peaks.

Within each experiment, split pooled cells into five accessible-site-depth quantiles. Calculate each region's binary accessibility rate by time and stratum, then standardize to the pooled stratum distribution, using strata represented at every time. Rates remain unsmoothed; Jeffreys-smoothed stratum rates stabilize binomial variance at zero/one boundaries. All cell counts and weights are exported.

## Temporal classes

A three-degree-of-freedom Wald test compares the 24/48/72 h rates to 0 h, retaining the covariance from their shared baseline. BH correction is across eligible overlapping external regions. Opening requires external q<0.05 and a positive 72-minus-0 h change in **both experiments**; closing uses negative directions.

For opening regions, half-rise time is the first crossing of `p(0) + 0.5 * [max(p(24),p(48),p(72)) - p(0)]`, linearly interpolated between sampling times.

| Class | Definition | Peaks | Genes represented |
|---|---|---:|---:|
| All opening | Dynamic, positive endpoints in both experiments | 229 | 128 |
| Early | Half-rise by 24 h | 191 | 117 |
| Middle | Half-rise after 24 h and by 48 h | 34 | 30 |
| Late | Half-rise after 48 h | 4 | 4 |
| Closing | Dynamic, negative endpoints in both experiments | 181 | 97 |
| Low-change profile | Largest/smallest observed time-point rate <=1.2 | 470 | 192 |

Low-change is a descriptive comparator, not an equivalence-test conclusion; it does not overlap opening/closing in the primary analysis. Earliest-20%, -25% and -33% timing quantiles contain 48, 58 and 76 peaks, retaining ties. These memberships depend only on the external reference.

`21_define_temporal_regions.py` finishes all memberships before target-effect testing. `run_temporal.py` hashes the membership files and checks that effect testing leaves them unchanged.

## Target preprocessing and statistical analysis

Use existing matched GSE208248 nuclei, preserving the established analysis in [Methods](METHODS.md). Four same-nucleus RNA/ATAC libraries represent two source lines, each in an undifferentiated and differentiated state. HMA annotations and HMA/Barthel comparisons are not changed.

Inherited joint QC: RNA UMI>=500, >=500 open ATAC peaks, removal of the upper 5% RNA/ATAC-depth tails within library, FRiP>=0.25, nucleosome signal<4, blacklist fraction<0.05, and TSS enrichment>=2 or >=3. TSS enrichment is the existing insertion-enrichment quality metric; it is distinct from promoter radius and from +/-5 kb plotting windows. Existing peak-coordinate, blacklist, cross-library mapping and 100-kb gene-neighborhood rules remain fixed.

COQ8A-high is raw RNA count **>=2 or >=3 UMI**; low is **exactly 1 UMI**. Zero-count nuclei are excluded from these inherited contrasts. Match high and low nuclei 1:1 without replacement within library on log1p RNA UMI and log1p open-peak count, requiring <=0.10 difference in each coordinate.

| TSS gate | COQ8A contrast | Matched pairs |
|---|---|---:|
| >=2 | >=2 vs 1 | 1,020 |
| >=3 | >=2 vs 1 | 958 |
| >=2 | >=3 vs 1 | 212 |
| >=3 | >=3 vs 1 | 201 |

For each fixed temporal set, each nucleus's score is its fraction of distinct peaks open. Fold is mean high fraction / mean low fraction. Report paired differences, 95% paired intervals, paired-nucleus p values, and library directions. Display >=3 versus 1 at TSS>=3 alongside all four settings.

The ten programme summaries form one BH family **within each QC/contrast/stage/mapping setting**, without correction across sensitivity settings. Gene screens use a separate BH family within each temporal module. Gene-neighborhood assignment is not an experimentally established enhancer-target link.

Direct module comparisons use paired differences scaled to each module's low-group baseline and 2,000 within-library paired bootstrap draws for fold-difference intervals. A source-adjusted HC3 OLS model of paired differences tests the undifferentiated-versus-differentiated interaction. These are nucleus-level associations from two source lines.

RNA gene count, mitochondrial UMI percentage, predicted-doublet exclusion and Scrublet tail exclusions are additional fixed-pair QC challenges. Remove a pair if either member fails. Explicitly retain zero-pair settings as unevaluable. These challenges do not silently alter the baseline.

## Reproduction

After producing the existing target tables with `run_all.py`:

```powershell
python run_temporal.py --raw C:\data\GSE109828 --gtf C:\data\gencode.v48.annotation.gtf.gz --h5-root C:\data\GSE208248_processed --work C:\work\GSE109828
```

Use `--reuse-prepared` only for completed external matrices from identical inputs. The wrapper downloads/checks inputs, prepares external matrices, freezes temporal sets, tests target contrasts, audits RNA QC, generates figures and validates summaries against direct H5 spot checks.

- `results/temporal/`: memberships, QC, complete effect/gene tables, input hashes and validation.
- `figures/temporal/`: PNG and vector PDF figures.
- Raw count files and intermediate sparse matrices remain in supplied raw/work directories, outside the release repository.

The 34-peak middle module is separate from the previous six-peak MYOD1 analysis and imports none of that analysis's selection rules. Results and figure-reading instructions: [Temporal results](TEMPORAL_RESULTS.md).

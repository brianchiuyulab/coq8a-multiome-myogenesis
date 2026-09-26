# COQ8A-associated myogenesis in a public same-nucleus multiome

Reproducible reanalysis of [GSE208248](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208248): four RNA/ATAC libraries from two human muscle source lines, each in a stem and differentiated state. This dataset is separate from the young/old HMA atlas. Its nuclei do not inherit HMA cell-type labels or Barthel Index values.

## Analysis sequence

The [extended exploratory screen](docs/EXTENDED_EXPLORATION.md) evaluates both link-discovery contrasts across the full gene set, promoter/distal and linked-peak definitions, and seven TSS windows. Add `--extended-exploration` to the reproduction command to generate its complete tables and [six-panel figure](figures/exploration/Figure_exploratory_route.png).

The biological search space is the measured union of MSigDB Hallmark and Reactome Myogenesis (**221 genes**). Two COQ8A count contrasts are analysed **in parallel through the complete 221-gene screen**, each matching nuclei within library on RNA and ATAC depth after joint QC at TSS enrichment ≥3: ≥2 versus 1 RNA UMI has 958 pairs and 6,871 eligible peaks; ≥3 versus 1 has 201 pairs and 3,731. The stronger contrast is used for the focal six-peak MYOD1 effect, with its full-screen result displayed beside the broader contrast.

The subsequent peak–RNA linkage and MYOD1 analyses are exploratory. Links are learned in the more numerous ≥2-versus-1 nuclei at TSS≥3 and TSS≥2; their ATAC direction is assessed **separately at both count thresholds** in two labelled full 221-gene rankings. MYOD1 ranks first under ≥2; under ≥3, CKB ranks first and MYOD1 second. The six-peak MYOD1 ratio of **1.280** uses TSS≥2 link discovery and ≥3-versus-1 effect testing at TSS≥3. Region/count thresholds and RNA doublet-score exclusions are sensitivity analyses, with no BH correction across those settings.

## Reproduction

### External differentiation timing

The [middle-region follow-up](docs/MIDDLE_FOLLOWUP.md) presents the [three phase panels](figures/temporal/Figure_three_temporal_phases.png) followed by an [all-34-peak heatmap](figures/temporal/Figure_middle_34_peak_heatmap.png). COL15A1 is the leading middle-region site (2.091-fold, p=0.000936, exploratory middle-only q=0.0318); the analysis record retains the full opening-region correction and explains the distinct families.

The [single-peak follow-up](docs/TEMPORAL_SINGLE_PEAKS.md) screens **all 191 early + 34 middle + 4 late sites together**, with the same membership in all four COQ8A/TSS settings. At >=3 vs 1 UMI/TSS>=3, early GNAO1 and middle COL15A1 neighborhood peaks each show about 2.09-fold accessibility (raw paired p=0.000346 and 0.000936; both BH q=0.0715 across 229 peaks). See the [complete peak-screen figure](figures/temporal/Figure_temporal_single_peaks.png). These local results and their common test family are separate from the module and gene-neighborhood summaries.

The completed [external temporal analysis](docs/TEMPORAL_RESULTS.md) defines region timing using GSE109828 before testing COQ8A. Of the 221-gene candidate space, a **34-peak middle-opening module (24-48 h half-rise)** has a **1.223** high/low accessibility ratio (paired p=0.001785; ten-module q=0.01785) at COQ8A>=3 versus 1 UMI and TSS>=3. The same module is 1.344 in undifferentiated libraries; the earliest-opening subsets are much weaker. All four COQ8A/TSS settings, external-definition sensitivities and RNA/doublet QC challenges are retained.

[Temporal Methods](docs/TEMPORAL_METHODS.md) describes the complete processing and `run_temporal.py` command. See [external timing](figures/temporal/Figure_external_timing.png), [COQ8A/TSS sensitivity](figures/temporal/Figure_COQ8A_temporal_sensitivity.png), and [states and genes](figures/temporal/Figure_temporal_states_and_genes.png). Source data are in `results/temporal/`. This result is separate from the earlier six-peak MYOD1 selection.

### Functional differentiation and fusion subsets

The [complete current analysis design](docs/CURRENT_ANALYSIS_DESIGN.md) details
both GEO cohorts, sample sizes, union/intersection membership, coordinate
mapping, temporal classes, QC, statistical families and reproduction commands.
Machine-readable [sample and scope inventories](results/design/) accompany it.

The [functional subset analysis](docs/FUNCTIONAL_SUBSETS.md) intersects four
external GO annotations with the 221-gene universe and fixed external dynamic
regions. It tests 54 distinct peaks and eight opening/closing modules using the
unchanged matched nuclei. Complete memberships and all four sensitivity
settings are in [results/functional](results/functional/).

### Original target-data preprocessing

Use Python 3.13 and install [`requirements.txt`](requirements.txt). Download the four `filtered_feature_bc_matrix.h5` and four `per_barcode_metrics.csv.gz` files for GSM6339597, GSM6339599, GSM6339601 and GSM6339603 from GSE208248. Download GENCODE v48 hg38 `annotation.gtf.gz`. Exact source-file checksums are in [`docs/input_sha256.csv`](docs/input_sha256.csv). With the versioned reference tables present, the command below rebuilds every analysis table, all three main figures and the QC supplement from those downloaded inputs.

```powershell
python -m pip install -r requirements.txt
python run_all.py --h5-root "C:\path\to\GSE208248_processed" --gtf "C:\path\to\gencode.v48.annotation.gtf.gz"
```

The command writes tables to `results/tables/` and vector PDF plus PNG figures to `figures/`. The versioned `reference/` tables and the committed 221-gene TSS fragment profile permit figure reproduction without downloading the large ATAC fragment files. With `--out-root`, the versioned profile is copied into the new output tree automatically. To rebuild that profile from indexed GSE208248 fragments in the same run, add `--fragments-dir` and install `pysam` from [`requirements_reference.txt`](requirements_reference.txt). The other reference inputs are documented in Methods. [`docs/SELECTION_AUDIT.md`](docs/SELECTION_AUDIT.md) enumerates every gene, peak and exploratory-locus selection. [`docs/FIGURE_LEGENDS.md`](docs/FIGURE_LEGENDS.md) defines each panel; [`docs/RESULTS.md`](docs/RESULTS.md) records the results and their statistical scope.

To regenerate the TSS profile on a Linux or macOS Python environment with `pysam`, after producing `matched_pairs.tsv.gz`:

```bash
python scripts/14_tss_fragment_profiles.py \
  --fragments-dir /path/to/indexed-GSE208248-fragments \
  --gtf /path/to/gencode.v48.annotation.gtf.gz \
  --genes reference/myogenesis_221_gene_sources.tsv \
  --pairs results/tables/matched_pairs.tsv.gz \
  --out results/tables/tss_fragment_profile_221.tsv.gz
python scripts/12_make_figures.py --tables results/tables --figures figures
```

`scripts/13_validate_release.py` checks pair identities and gates, the focal full-test-family summary, both complete 221-gene rankings, the TSS profile shape and group sizes, and the six-region fold against pair-level data. GitHub Actions runs these checks on the released tables.

| Figure | Purpose |
|---|---|
| [Figure 1](figures/main/Figure_1_Multiome_and_TSS.pdf) | Four-library comparison and TSS-aligned ATAC fragment profiles for all 221 genes |
| [Figure 2](figures/main/Figure_2_Global_Screen_and_Ranking.pdf) | Both complete peak and 221-gene region screens, with separately labelled same-data rankings |
| [Figure 3](figures/main/Figure_3_MYOD1_Locus_and_Sensitivity.pdf) | MYOD1 locus, region-set and contrast sensitivity, library effects and doublet-score challenge |
| [Supplementary QC](figures/supplement/Supplementary_Figure_QC.pdf) | TSS enrichment and matched-depth balance |

Both complete contrasts show an RNA myogenesis association, but no ATAC gene region or candidate peak passes its corresponding full-family BH q < 0.05 threshold. This observational dataset contains two independent source lines; paired-nucleus p values describe within-dataset associations and do not establish COQ8A-driven chromatin opening.

### Bidirectional temporal follow-up

[Peak-to-RNA results and reproduction](docs/BIDIRECTIONAL_RESULTS.md) cover all
410 externally defined opening/closing peaks, including CSRP3 and CAV3 as
phenotype-relevant exploratory candidates. Scripts 32, 30 and 33 generate the
public evidence tables. Scripts 35 and 34 optionally cross-check existing private
RNA contrasts and plot the comparison; private RNA values and figures are not
included in this repository.

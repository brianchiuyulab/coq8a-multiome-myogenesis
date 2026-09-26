# COQ8A-associated myogenesis in a public same-nucleus multiome

Reproducible reanalysis of [GSE208248](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208248): four RNA/ATAC libraries from two human muscle source lines, each in a stem and differentiated state. This dataset is separate from the young/old HMA atlas. Its nuclei do not inherit HMA cell-type labels or Barthel Index values.

## Analysis sequence

The biological search space is the measured union of MSigDB Hallmark and Reactome Myogenesis (**221 genes**). The pipeline first applies joint nucleus QC, matches COQ8A ≥2 UMI nuclei to 1 UMI nuclei within each library on RNA and ATAC depth, and tests every eligible gene region and candidate peak. The primary gate is TSS enrichment ≥3. The complete primary test families are saved before any locus is ranked.

The subsequent peak–RNA linkage and MYOD1 analyses are exploratory. Link discovery at TSS≥2 and TSS≥3, stronger COQ8A-count contrasts, and RNA doublet-score exclusions are reported as sensitivity analyses. The six-peak MYOD1 ratio of **1.280** is reproducible under its stated TSS≥2 link-discovery / TSS≥3 effect-test setting; it is not the result of the primary 221-gene ATAC scan.

## Reproduction

Use Python 3.13 and install [`requirements.txt`](requirements.txt). Download the four `filtered_feature_bc_matrix.h5` and four `per_barcode_metrics.csv.gz` files for GSM6339597, GSM6339599, GSM6339601 and GSM6339603 from GSE208248. Download GENCODE v48 hg38 `annotation.gtf.gz`. Exact source-file checksums are in [`docs/input_sha256.csv`](docs/input_sha256.csv).

```powershell
python -m pip install -r requirements.txt
python run_all.py --h5-root "C:\path\to\GSE208248_processed" --gtf "C:\path\to\gencode.v48.annotation.gtf.gz"
```

The command writes tables to `results/tables/` and vector PDF plus PNG figures to `figures/`. The versioned `reference/` tables permit reproduction without reprocessing the large ATAC fragment files. [`docs/METHODS.md`](docs/METHODS.md) describes how to regenerate these inputs from fragments, genome sequence and the public matrices. [`docs/SELECTION_AUDIT.md`](docs/SELECTION_AUDIT.md) enumerates every gene, peak and exploratory-locus selection. [`docs/FIGURE_LEGENDS.md`](docs/FIGURE_LEGENDS.md) defines each panel; [`docs/RESULTS.md`](docs/RESULTS.md) records the results and their statistical scope.

`scripts/13_validate_release.py` checks pair identities and gates, the full primary test-family summary, and the six-region fold against pair-level data. GitHub Actions runs these checks on the released tables.

| Figure | Purpose |
|---|---|
| [Figure 1](figures/main/Figure_1_global_discovery.pdf) | Four-library paired design and directly aligned Hallmark RNA/ATAC effects |
| [Figure 2](figures/main/Figure_2_Global_ATAC_Scan.pdf) | Complete primary peak and gene-region ATAC test families |
| [Supplementary MYOD1](figures/supplement/Supplementary_Figure_MYOD1_Exploratory.pdf) | Locus-level link discovery, threshold sensitivity and library effects |
| [Supplementary QC](figures/supplement/Supplementary_Figure_QC.pdf) | TSS enrichment and matched-depth balance |

The primary analysis shows an RNA myogenesis association, but no ATAC gene region or candidate peak passes its full-family BH q < 0.05 threshold. This observational dataset contains two independent source lines; paired-nucleus p values describe within-dataset associations and do not establish COQ8A-driven chromatin opening.

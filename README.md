# COQ8A and myogenesis in GSE208248 same-nucleus RNA + ATAC

Exploratory, fully scripted reanalysis of four 10x Multiome libraries from two
human muscle source lines, each sampled in a stem and differentiated state.
The data set is **not** the young/old HMA donor atlas and does not inherit its
cell type annotations, BI comparisons, or donor-level estimates.

For the fixed, non-targeted main analysis and its actual decision, start with
[`docs/PRIMARY_ANALYSIS_zh-TW.md`](docs/PRIMARY_ANALYSIS_zh-TW.md). The
six-peak 1.280 MYOD1 result is reproduced separately as an exploratory
sensitivity, not used to define the main search.

## Question and analysis order

1. Start with independently named MSigDB Hallmark Myogenesis and Reactome
   Myogenesis sets (`reference/myogenesis_221_gene_sources.tsv`). Their union
   contains 221 measured genes, including MYOD1, MYOG, MYF5, and MYF6.
2. Apply joint RNA/ATAC nucleus QC and match each COQ8A-high nucleus to one
   COQ8A-low nucleus within its original library on RNA and ATAC depth.
3. Test RNA programmes and nearby ATAC accessibility across **all 221 genes**.
   Define candidate ATAC regions by genome position and peak recurrence before
   examining the COQ8A effect. Test all eligible peaks, not only MRF motifs.
4. Link candidate peaks to RNA of their nearby genes with depth/group-adjusted
   within-library association at TSS≥3 and, separately, TSS≥2; then rank
   exploratory loci. Show the full search and the resulting MYOD1 example.
5. Challenge the locus with region-set, COQ8A count, TSS, and RNA doublet-score
   sensitivity analyses. The peak-to-gene-linked subset is learned in the same
   nuclei, so its apparent effect is **discovery**, not independent validation.

No fold-change cutoff was used to select a result. `fold_open` is the high/low
ratio of the **fraction of nuclei with an open region**. The accompanying
`delta_pp` is the absolute percentage-point difference.

## Reproduce

Python 3.13 was used with versions recorded in `requirements.txt`. Obtain the
four GEO filtered feature-barcode H5 files and matching per-barcode-metrics
files for GSM6339597, GSM6339599, GSM6339601, GSM6339603 from
[GSE208248](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208248).
Provide the directory that contains those eight files. Download GENCODE v48
hg38 `annotation.gtf.gz` and provide its path.
`docs/input_sha256.csv` records the exact nine input files used in this run.
In this workspace they are preserved under `data/public/`, outside the
pending-manual-deletion directory.

```powershell
python -m pip install -r requirements.txt
python run_all.py --h5-root "C:\path\to\GSE208248_processed" --gtf "C:\path\to\gencode.v48.annotation.gtf.gz"
```

The run uses the committed, barcode-indexed reference QC tables and motif
scores in `reference/`. The larger raw H5, fragments, 2bit genome and GTF are
not in Git. `scripts/00_fragment_qc_reference.py` reproduces the fragment
quality tables from indexed ATAC fragments; `scripts/motif_score_provenance.py`
reproduces JASPAR scores with an hg38 2bit genome; and
`scripts/scrublet_provenance.py` regenerates RNA doublet scores. The exact
commands and source links are in `docs/METHODS.md`.

## Results at a glance

Primary comparison: TSS enrichment ≥3, COQ8A 2+ versus 1 UMI, depth matched
within each library, 958 pairs. RNA Hallmark Myogenesis rises in all four
libraries (paired difference +0.0107 log1p(CP10K), exploratory `p=0.0268`,
four-programme `q=0.0358`). RNA MRF loci rise similarly (`p=0.00185`,
`q=0.00738`). Nearby ATAC across the full Hallmark set has no consistent
effect (two of four libraries; `p=0.0519`, `q=0.207`). No individual gene
region or candidate ATAC peak survives its full multiple-testing family at
`q<0.05`.

MYOD1's 19 four-library common nearby peaks show a modest high/low open
fraction ratio of 1.069 (absolute +0.824 percentage points; pair-level
`p=0.0395`, 221-gene `q=0.635`). A six-peak subset nominated with TSS≥2
peak–RNA links reproduces the earlier **1.280** ratio under COQ8A 3+ versus
1 UMI and TSS≥3 (+2.90 percentage points, 201 pairs, 4/4 libraries positive,
pair-level `p=0.02558`, 20-setting exploratory `q=0.0755`). Requiring TSS≥3
already at the link-discovery stage removes one low-support peak, leaving five
peaks and a **1.241** ratio (`p=0.05295`, 3/4 positive). The six-peak ratio
drops to ~1.19 after removing the highest 2.5–5% RNA doublet-score nuclei.
This is a **candidate locus signal**, not robust evidence for a programme-wide
COQ8A-associated accessibility increase.

`docs/RESULTS.md` gives the full interpretation and sensitivity summary.
The superseded exploratory GSE208248 files were moved, without deletion, to
`tmp/COQ8A_unused_sensitivity_pending_manual_delete_20260926/prior_GSE208248_exploratory_files/`;
that directory contains a file manifest. Raw H5 and GTF inputs were moved out
of the pending-deletion area first.

## Layout

| Directory | Contents |
|---|---|
| `scripts/` | Numbered analysis, reconciliation and plotting code plus shared methods |
| `reference/` | Fixed gene list, QC, motif and blacklist inputs |
| `results/tables/` | All nucleus, region, link, effect and sensitivity tables |
| `figures/main/` | Figure 1: discovery funnel and full-gene results; Figure 2: MYOD1 locus and sensitivity |
| `figures/supplement/` | RNA/ATAC QC, matching diagnostics and complete peak-level primary scan |
| `docs/` | Methods, results, figure legends and source provenance |

The key sources for figure and method conventions are the
[GSE208248 article](https://pmc.ncbi.nlm.nih.gov/articles/PMC10123345/),
[Signac QC/visualization](https://stuartlab.org/signac/1.13.0/articles/pbmc_vignette),
[ArchR peak-to-gene linkage](https://www.archrproject.com/bookdown/peak2genelinkage-with-archr.html),
and [single-cell ATAC differential-analysis best practices](https://www.nature.com/articles/s41467-024-53089-5).

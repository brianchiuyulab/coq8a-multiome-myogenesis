# Data availability and source-data index

## Current Day-0 source data

The current high/low tests use only GSM6339597 and GSM6339601 H5 matrices,
with fragment files GSM6339598 and GSM6339602. The four original libraries
contribute to the common peak catalog. No GSE109828 temporal labels or
GSE240061 measurements enter this fixed comparison.

| Location | Current contents |
|---|---|
| `results/fixed_day0_D206/peak_membership.tsv` | 25 gene neighborhoods and 565 distinct peaks |
| `results/fixed_day0_D206/matched_pairs.tsv` | Exact 527 high/low barcode pairs |
| `results/fixed_day0_D206/all_peaks.tsv` | Complete 565-peak screen |
| `results/fixed_day0_D206/all_regions.tsv` | 25 regions in each direction |
| `results/fixed_day0_D206/local_RNA_candidates.tsv` | Complete 170 local RNA candidates |
| `results/fixed_day0_D206/links_by_source.tsv`, `links_combined.tsv` | Source-specific and combined models |
| `results/fixed_day0_D206/RNA_high_low.tsv` | Same-pair RNA contrasts |
| `results/fixed_day0_D206/fragment_profiles.tsv.gz`, `profile_windows.tsv` | Figure insertion counts and window coordinates |
| `results/fixed_day0_D206/audit/` | Audit results, software versions and source-table SHA256 |

## Public source data

| Resource | Data used | Access |
|---|---|---|
| GSE208248 | Four paired RNA/ATAC filtered H5 matrices and barcode metrics; four indexed ATAC fragment files | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208248) |
| GSE109828 | GSM2970930 and GSM2970931 sci-ATAC sparse counts and time metadata | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109828) |
| GENCODE v48 | GRCh38 gene/transcript/exon coordinates | [GENCODE release 48](https://www.gencodegenes.org/human/release_48.html) |
| MSigDB | Hallmark/Reactome myogenesis and the four specified GO sets | [MSigDB human collections](https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp) |
| ENCODE blacklist | hg38 v2 intervals | [Boyle laboratory blacklist](https://github.com/Boyle-Lab/Blacklist) |
| UCSC | hg38-to-hg19 liftOver chain | [UCSC hg38 liftOver](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/) |

Target H5/RNA accessions are GSM6339597, GSM6339599, GSM6339601 and GSM6339603.
For each, retain the original `*_filtered_feature_bc_matrix.h5` and
`*_per_barcode_metrics.csv.gz` names in the `--h5-root` directory.
Fragment accessions are GSM6339598, GSM6339600, GSM6339602 and GSM6339604;
retain original `*_atac_fragments.tsv.gz` and `.tbi` names in `--fragments`.
Download links are the supplementary files in the corresponding GEO records.
`19_fetch_temporal_reference.py` downloads the external count/time files and
chain automatically. External and target file hashes are recorded in
`results/temporal/download_manifest.json` and `docs/input_sha256.csv`.

## Reusable source tables

| Location | Content |
|---|---|
| `reference/myogenesis_221_gene_sources.tsv` | Frozen biological gene universe |
| `reference/fragment_qc/` | Per-nucleus fragment-QC measurements |
| `results/design/` | Samples, cell counts, union/intersection and peak inventories |
| `results/tables/matched_pairs.tsv.gz` | Library-scoped high/low barcode identities for every setting |
| `results/tables/candidate_peak_gene.tsv.gz` | Full candidate peak-to-nearby-gene map |
| `results/temporal/temporal_set_memberships.tsv.gz` | Externally defined temporal region sets |
| `results/temporal/temporal_effect_grid.tsv` | Complete temporal-module effects and statistics |
| `results/temporal/bidirectional_peak_effects.tsv` | All 410 dynamic-region effects in every setting |
| `results/functional/functional_peak_memberships.tsv` | Function-defined 54-region membership |
| `results/functional/functional_peak_effects.tsv` | All 54-region tests and q values |
| `results/temporal/dynamic_cis_links.tsv` | Expanded partial peak–RNA association screen |
| `results/figure_source/external_410_regions.tsv` | Figure 1 external profiles and region ordering |
| `results/figure_source/functional_54_order.tsv` | Figure 2 row IDs, coordinates and effects |
| `results/figure_source/fragment_profiles.tsv.gz` | Per-library insertion counts and normalized profiles |
| `results/figure_source/focal_statistics.tsv` | Figure 3 peak/RNA statistics and link evidence |
| `results/figure_source/focal_transcripts.tsv` and `focal_exons.tsv` | Figure 3 gene models |
| `results/tables/tss_fragment_profile_221.tsv.gz` | Supplementary TSS profiles |

Figure legends specify which correction family each displayed q value belongs
to. Full tables retain all evaluated directions, including nonsignificant
results. Candidate gene names are genomic-proximity annotations unless a
separate statistical link is explicitly reported.

## Private C2C12 data

The owner's P11/P22/P33 differentiation RNA data and derived private figures
are not included in this public repository. They comprise a separate passage/time
comparison and must not be described as the public COQ8A-high/low comparison
or as an OE/KD RNA experiment.

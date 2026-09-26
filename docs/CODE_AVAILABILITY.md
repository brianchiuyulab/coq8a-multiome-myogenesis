# Code availability

Analysis and figure-generation code is available at
[brianchiuyulab/coq8a-multiome-myogenesis](https://github.com/brianchiuyulab/coq8a-multiome-myogenesis).
The repository includes frozen biological memberships, processed source tables,
matched barcode identities, normalization definitions, software requirements,
figure legends and validation code. Record the Git commit used with a manuscript
submission; the repository does not currently provide an archival DOI.

## Entry point

`run_paper.py` is the supported entry point. `--mode figures` regenerates
three main and two supplementary figure sets from committed data. `--mode analysis`
recomputes the analysis from author-deposited processed matrices and frozen
references. Raw-read alignment and peak calling are outside this pipeline.
Use Python 3.13 and `requirements.txt`; fragment regeneration additionally
requires `pysam`, supplied in `requirements_reference.txt`.

## Code map

| Stage | Scripts |
|---|---|
| Input extraction and reference QC | 00_prepare_reference, 00_fragment_qc_reference, 01_extract_nuclei |
| Joint nucleus QC and depth matching | 02_qc_and_matching |
| Peak correspondence and gene neighborhoods | 03_define_regions |
| Complete RNA and candidate peak effects | 04_gene_programme_effects, 05_peak_effects |
| External download, preparation and temporal definition | 19_fetch_temporal_reference, 20_prepare_temporal_reference, 21_define_temporal_regions |
| Temporal-module and single-region association | 22_test_temporal_accessibility, 27_temporal_single_peaks, 32_bidirectional_peak_screen |
| External functional membership and tests | 37_functional_subsets |
| Same-setting cis association and evidence | 07_peak_gene_links, 30_middle_cis_links, 33_dynamic_evidence |
| Sample/selection inventory and figure source preparation | 38_design_inventory, 39_prepare_figure_sources |
| Fragment profiles and rendering | 14_tss_fragment_profiles, 40_extract_fragment_windows, 36_plot_tss_heatmap, 41_make_submission_figures |
| Numerical validation | 26_validate_temporal, 42_validate_submission |
| Optional RNA-quality sensitivity | 25_temporal_rna_qc, scrublet_provenance |
| Private passage/time RNA comparison | 35_crosscheck_existing_rna, 43_plot_private_rna |

Numeric prefixes identify scripts and do not imply that every optional script
is executed in the standard pipeline. Private RNA scripts require owner-supplied
inputs and write outside the public repository.

## Validation

```bash
python run_paper.py --mode figures
python scripts/26_validate_temporal.py --tables results/tables --temporal results/temporal
```

Validation recomputes paired single-peak p values and BH q values, checks fixed
membership counts, compares focal statistics with complete-screen outputs,
checks profile normalization, and reconciles fragment-profile barcode/window
hashes with the current analysis. The temporal validator recomputes programme
statistics from matched-pair scores. Supplying `--h5-root` to that validator
also checks example scores directly against the deposited count matrices.

The GitHub Actions workflow regenerates figures and runs these checks. This
does not execute the full multi-gigabyte input pipeline on every commit.

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
| Expanded middle-region target associations | 44_expand_middle_links, 45_middle_count_models.R |
| Independent MYOD1 annotation and follow-up figure | 46_middle_myod_annotation, 47_report_middle_links |
| Unstratified candidate screen and count models | 48_unstratified_links, 49_unstratified_count_models.R, 50_select_unstratified_followup |
| Unstratified validation and visualization | 51_report_unstratified_links |
| Independent human enhancer intersections and locus annotation | 52_external_enhancer_sets, 53_report_external_enhancers |

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

The external enhancer-region exploration is implemented in
`scripts/54_enhancer_region_tests.py` (external block assignment, exact paired
regional tests, and source-level summaries) and
`scripts/55_external_enhancer_selectivity.py` (eight-reference download,
independent selectivity annotation, and fixed-set tests). Reproduction commands
and inference scopes are in `docs/ENHANCER_REGION_EXPLORATION.md`.

The independent C2C12 reference uses `scripts/56_external_c2c12_timecourse.R`
(DESeq2 across 18 RNA samples), `scripts/57_external_mouse_peak_annotation.py`
(streaming hg38-to-mm10 chain mapping and external ATAC annotation), and
`scripts/58_compare_candidate_evidence.py` (public candidate evidence tables).
Reproduction instructions are in `docs/CANDIDATE_EXPANSION_AND_EXTERNAL_DATA.md`.

The GitHub Actions workflow regenerates figures and runs these checks. This
does not execute the full multi-gigabyte input pipeline on every commit.

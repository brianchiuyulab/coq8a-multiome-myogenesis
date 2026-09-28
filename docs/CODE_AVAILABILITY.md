# Code availability

Analysis code, figure source tables, matched barcodes, reference memberships
and computational checks are available at
[brianchiuyulab/coq8a-multiome-myogenesis](https://github.com/brianchiuyulab/coq8a-multiome-myogenesis).
Specify the repository commit with a manuscript submission. No archival DOI
has been assigned.

## Current entry point

Use `run_fixed_day0.py` for the current Day-0 analysis. The older `run_paper.py`
reproduces the historical temporal workflow.

Python dependencies are pinned in `requirements.txt`. Fragment extraction
uses Linux `pysam` (see `requirements_reference.txt`). The exact audit runtime
is recorded in `results/fixed_day0_D206/audit/audit_report.json`.

### Regenerate figures from deposited source tables

```bash
python -m pip install -r requirements.txt
python run_fixed_day0.py --mode figures
```

### Recompute the fixed comparison from processed H5 matrices

Place the original GSE208248 H5 files in `/data/GSE208248`; choose an empty
writable cache directory. This route uses the repository's versioned peak
mapping and fragment-QC measurements, extracts RNA/ATAC counts afresh, repeats
matching and all peak/regional tests, fits local RNA models, then renders figures.
The 2,000,000 regional permutations are recomputed; no full sensitivity grid
is needed.

```bash
python run_fixed_day0.py --mode analysis --h5 /data/GSE208248 --work /work/day0_cache
```

To additionally rebuild the four-library mapping, RNA/ATAC depth inventory and
joint QC gates, include all four H5 files and their original barcode-metric
files, plus GENCODE v48:

```bash
python run_fixed_day0.py --mode analysis --h5 /data/GSE208248 --work /work/day0_cache --gtf /data/gencode.v48.annotation.gtf.gz --rebuild-upstream
```

This reuses deposited per-nucleus fragment-QC scores. To regenerate those
scores themselves, use `00_fragment_qc_reference.py` for each fragment library
and each of `tss`, `nuc`, and `blacklist`, with `--reference-root reference`,
`--fragment-root /data/fragments`, and `--out-root reference/fragment_qc`.
No FASTQ alignment or new peak calling is performed.

### Regenerate positional fragment profiles

Run in a Python environment with `pysam`, using the original two Day-0
indexed fragment files:

```bash
python run_fixed_day0.py --mode figures --fragments /data/fragments --gtf /data/gencode.v48.annotation.gtf.gz
```

### Audit raw matrices and annotations

After creating the cache, the independent audit checks every peak and local
RNA model against raw H5 values, rebuilds scope from GENCODE/GO annotations,
and verifies joint QC and barcodes. Add `--permutations` to repeat all
2,000,000 swaps and compare every exceedance count.

```bash
python run_fixed_day0.py --mode audit --h5 /data/GSE208248 --work /work/day0_cache --gtf /data/gencode.v48.annotation.gtf.gz --permutations
```

## Code map

| Step | Implementation |
|---|---|
| Matrix depths and COQ8A counts | `scripts/01_extract_nuclei.py` |
| Fragment QC and joint gates | `scripts/00_fragment_qc_reference.py`, `scripts/02_qc_and_matching.py` |
| Four-library peak mapping and gene neighborhoods | `scripts/03_define_regions.py` |
| Count extraction, matching, paired permutations | Shared functions in `scripts/88_day0_sensitivity_grid.py` |
| External 25-gene scope and complete fixed statistics | `scripts/98_recompute_fixed_statistics.py` |
| Region nomination, all local RNA models and RNA contrasts | `scripts/92_fixed_day0_workflow.py` |
| Screening, source comparisons and RNA figures | `scripts/93_plot_fixed_day0.py` |
| Raw fragment profiles | `scripts/94_fixed_fragment_profiles.py` |
| Positional heatmaps and genomic tracks | `scripts/96_plot_fragment_panels.py` |
| Independent matrix/annotation/statistical audit | `scripts/97_audit_fixed_day0.py` |
| Published table and profile validation | `scripts/99_validate_fixed_tables.py` |
| Private C2C12 cross-reference | `scripts/95_private_target_crosscheck.py`, external output only |

The local transcript candidate catalog is archived in
`results/unstratified/candidate_pairs.tsv.gz`; script 48 constructs it from
GENCODE transcript starts and uniquely named deposited RNA features. Script
97 independently reconstructs every local candidate for the current selected
peaks. Script 44 supplies the shared sparse-matrix extraction helper.

[Methods and figure legends](FIXED_DAY0_D206_WORKFLOW.md) describe all statistics,
units and display choices. [The audit report](FIXED_DAY0_AUDIT.md) records the
scope of verification. Historical code maps remain in
[CODE_AVAILABILITY_HISTORY.md](CODE_AVAILABILITY_HISTORY.md).

# Methods and provenance

## Study and inputs

This is a retrospective, exploratory analysis of [GSE208248](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208248), following the [source publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC10123345/). Four same-nucleus RNA/ATAC libraries were used: GSM6339597 (line 1 stem), GSM6339599 (line 1 differentiated), GSM6339601 (line 2 stem), GSM6339603 (line 2 differentiated). The matching ATAC fragment accessions are GSM6339598, GSM6339600, GSM6339602, GSM6339604. Libraries are conditions from **two** source lines; they are not four independent human donors. This analysis never merges nuclei across libraries for matching.

The 221-gene search universe is the measured union of [MSigDB Hallmark Myogenesis](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_MYOGENESIS.html) (198 genes here) and [Reactome Myogenesis](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/REACTOME_MYOGENESIS.html) (29 genes here). `reference/myogenesis_221_gene_sources.tsv` records the exact membership of each. This biological scope was defined for the final analysis; it is not a genome-wide all-pathway screen or a preregistered hypothesis test. Named MRF and MEF2 locus summaries are contextual secondary views and do not restrict region discovery.

## 1. Joint nucleus QC and matching

`01_extract_nuclei.py` reads each 10x filtered H5 matrix and barcode metrics. Its COQ8A, total RNA UMI and open-peak counts were extracted for all 41,622 filtered barcodes. `00_fragment_qc_reference.py` computes fragment-based TSS enrichment, nucleosome-free/mononucleosomal fragment ratio, and blacklist fraction. The TSS score uses per-base insertions within ±500 bp of GENCODE v48 gene TSS, divided by the per-base signal in the two outer 100-bp flanks of a ±1 kb window. This is a documented Signac-style geometry, implemented here directly rather than by Signac. Blacklist is ENCODE hg38 v2. FRiP is 10x peak-region fragments / ATAC fragments.

`02_qc_and_matching.py` requires at least 500 RNA UMI and 500 open ATAC peaks, removes the top 5% of either depth within each library, then requires FRiP≥0.25, nucleosome signal<4, blacklist fraction<0.05. TSS≥3 is primary and TSS≥2 sensitivity. Nuclei are contrasted as COQ8A ≥2 UMI versus exactly 1 UMI; ≥3 versus 1 is the stronger-count sensitivity. COQ8A-zero nuclei are not used because zero is ambiguous under sparse detection. Within each library and condition, high and low nuclei are matched 1:1 without replacement on log1p RNA UMI and log1p open-peak counts, with each coordinate difference ≤0.10. The matching is independent of the tested gene and peak outcomes. The primary set has 958 pairs (211/326/296/125 across the four libraries); the ≥3 set has 201 (46/74/58/23). **All four gate/count combinations are reported**.

## 2. Label-blind region search

`03_define_regions.py` restricts peaks to canonical hg38 chromosomes and width 200–2000 bp, removes blacklist-overlapping peaks, and maps library peaks one-to-one to the first library's peak coordinates using at least 50% reciprocal overlap. Candidate regions are within ±100 kb of **any protein-coding GENCODE v48 transcript TSS** for a gene in the 221-gene universe and observed in ≥3 of four libraries. This yields 7,699 candidate peak–gene pairs. A stricter all-four-library common set is used for the gene-level accessibility summary. A nearby peak is a *candidate* regulatory region, not a proven enhancer for that gene. Neither MRF motifs nor COQ8A outcome is used to enter this search space.

## 3. Broad RNA and ATAC tests

`04_gene_programme_effects.py` uses RNA log1p(count / total RNA UMI × 10,000). For each gene, the ATAC value is the fraction of its all-four-library common nearby peaks open in a nucleus. Programme values average available gene values within the stated set. It reports paired high-minus-low means, 95% pair-level t intervals and two-sided one-sample paired t p values, plus number of libraries with a positive effect. BH q values are computed separately across 221 gene-level tests and across four programme summaries, for each gate/count/modality. These p and q values describe nucleus-level association; they do **not** supply independent-source inference with only two source lines.

## 4. Peaks, links, motifs and ranking

`05_peak_effects.py` tests COQ8A high/low accessibility for all candidate peaks and applies BH over the full candidate-peak test family. `06_primary_decision.py` records the complete primary test-family results before any target ranking.

`07_peak_gene_links.py` checks each candidate peak–gene pair in the primary TSS≥3 matched nuclei and repeats link discovery at TSS≥2 for sensitivity. Within each library, it correlates binary peak accessibility with normalized RNA after residualizing both on an intercept, COQ8A group, log RNA depth and log ATAC depth. A peak must be open in ≥10 matched nuclei and a gene expressed in ≥20. Fisher-z values are weighted across ≥3 evaluable libraries and BH corrected across all 4,925 evaluable candidate links at TSS≥3 or 5,054 at TSS≥2. This peak-to-gene *association* can indicate a functional hypothesis; it is not a causal or physical link. JASPAR 2024 MRF/MEF2 PWM scores are attached **after** this screen for interpretation, not used to nominate peaks. The published [`peak_gene_links.tsv`](../results/tables/peak_gene_links.tsv) and [`peak_gene_links_tss2.tsv`](../results/tables/peak_gene_links_tss2.tsv) contain every tested link at each gate.

`08_candidate_ranking.py` combines positive same-data RNA links with the direction of the COQ8A peak effect. That ranking is exploratory and not a held-out confirmation. MYOD1 was selected from this transparent 221-gene ranking, not treated as the only possible target. `09_link_gate_reconciliation.py` records peak-level support at both link-discovery QC gates.

## 5. Region and QC sensitivity

`10_locus_sensitivity.py` measures MYOD1 in five region sets: all 32 nearby candidates, 19 four-library common peaks, 10 peaks with positive RNA links under TSS≥3, five of those without a strong MRF motif score, and **six** weak-MRF linked peaks identified by repeating link discovery at TSS≥2. It crosses these with TSS≥2/≥3 for the *effect test* and COQ8A≥2/≥3 vs 1 UMI: 20 visible tests, BH-corrected as a descriptive family. The TSS≥2-linked six-peak set is evaluated on the **same TSS≥3 201 pairs** as the five-peak set for the 1.280 versus 1.241 comparison. Here `fold_open` = aggregate high open fraction / aggregate low open fraction; `delta_pp` = 100 × their absolute difference. Linked subsets are learned using overlapping nuclei and the ≥3-UMI contrast is low powered. `11_doublet_sensitivity.py` fixes the matched-pair IDs and removes pairs involving nuclei in the top 0, 1, 2.5, 5 and 10% within-library RNA Scrublet score tail. This is a robustness probe, not complete ATAC doublet calling.

## 6. Figures

`12_make_figures.py` produces PDF and PNG from the committed tables. Figure 1 shows design, programme-level effects and all 221 gene effects. Figure 2 displays the complete primary peak and gene-region ATAC test families. The exploratory MYOD1 locus and threshold sensitivity appear in a supplementary figure, alongside a separate QC and depth-matching supplement. The figure structure follows common locus/forest/quality-control conventions of [Signac](https://stuartlab.org/signac/articles/visualization), [ArchR](https://www.archrproject.com/bookdown/peak2genelinkage-with-archr.html) and published [single-cell ATAC analysis guidance](https://www.nature.com/articles/s41467-024-53089-5).

## Regenerating optional reference tables

The small barcode-indexed reference tables are committed so `run_all.py` does not need multi-gigabyte fragments or hg38 2bit. To regenerate them, obtain indexed ATAC fragments and hg38 2bit. `00_prepare_reference.py` derives barcode and TSS reference inputs from the public H5 extraction and GTF. The generated tables were checked against the committed files: the barcode tables match exactly; the 20,044 TSS coordinate sets match exactly, with only a different row order. The fragment script also requires `hg38-blacklist.v2.bed.gz` in `reference/`:

```powershell
python -m pip install -r requirements_reference.txt
python scripts/01_extract_nuclei.py --h5-root C:\path\to\GSE208248_processed --out results\tables\nuclei.tsv.gz
python scripts/00_prepare_reference.py --nuclei results\tables\nuclei.tsv.gz --gtf C:\path\to\gencode.v48.annotation.gtf.gz --out reference
python scripts/00_fragment_qc_reference.py GSM6339598 tss --reference-root reference --fragment-root C:\path\to\fragments --out-root reference\fragment_qc
python scripts/00_fragment_qc_reference.py GSM6339598 nuc --reference-root reference --fragment-root C:\path\to\fragments --out-root reference\fragment_qc
python scripts/00_fragment_qc_reference.py GSM6339598 blacklist --reference-root reference --fragment-root C:\path\to\fragments --out-root reference\fragment_qc
python scripts/scrublet_provenance.py --h5-root C:\path\to\GSE208248_processed --barcodes reference\gse208248_qc_barcodes.tsv --out reference\rna_scrublet_qc.tsv.gz
python scripts/motif_score_provenance.py --consensus results\tables\consensus_peak_map.tsv.gz --matrices reference\jaspar2024_mrf_mef2_matrices.json --genome-2bit C:\path\to\hg38.2bit --out reference\jaspar2024_peak_scores.tsv.gz
```

Repeat the fragment command for the remaining three ATAC GSMs and all three modes. The full downstream pipeline uses the existing fixed reference tables to keep the reproduction command short. All parameter choices above are recorded retrospectively; they should not be described in a manuscript as prospective preregistration.

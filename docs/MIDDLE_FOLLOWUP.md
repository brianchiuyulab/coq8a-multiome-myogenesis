# Middle-opening region follow-up

## Analysis sequence

External GSE109828 timing assigns 191 early, 34 middle and four late peaks. GSE208248 COQ8A associations use the existing same-nucleus RNA/ATAC data, joint QC and within-library depth matching. Both modalities refer to the same barcode within a nucleus; the high-versus-low comparison subsequently matches two distinct nuclei within the same library.

All displayed stages use COQ8A >=3 versus exactly 1 RNA UMI, TSS enrichment >=3, and the same 201 matched pairs. The other three COQ8A/TSS settings remain in the complete results. RNA UMI >=500, open ATAC peaks >=500, within-library top 5% depth exclusions, FRiP >=0.25, nucleosome signal <4, blacklist fraction <0.05 and matching on log1p RNA/ATAC depths are unchanged; see [Temporal Methods](TEMPORAL_METHODS.md) and original Methods for definitions. No parameter changes between the phase comparison and the 34-peak follow-up.

## Three-phase figure

[Figure](../figures/temporal/Figure_three_temporal_phases.png). Each panel shows the average fraction of the specified peak set accessible per nucleus. Colored lines show the four libraries in accession order (line 1 stem, line 1 differentiated, line 2 stem, line 2 differentiated); black diamonds show the matched-pair pooled estimate. The same nucleus pairs are used in every panel. The paired p values are nucleus-level tests; q values retain the original ten-programme BH family. The figure displays three of those programmes. Two source lines provide the four libraries.

| Phase | Peaks | High/low ratio | Paired p | Programme q |
|---|---:|---:|---:|---:|
| Early | 191 | 1.063824 | 0.068270 | 0.174759 |
| Middle | 34 | 1.222930 | 0.001785 | 0.017848 |
| Late | 4 | 1.034483 | 0.739788 | 0.965541 |

## Middle peak heatmap

[Figure](../figures/temporal/Figure_middle_34_peak_heatmap.png). Rows are all 34 middle peaks, ordered by two-sided exact paired-binomial p. Left: pooled fractions of open nuclei in COQ8A-low and COQ8A-high groups. Middle: high-minus-low percentage-point differences in each library. Right: fold change, raw p and BH q across the 34 middle peaks. Shading marks raw p<0.05. Gene labels denote nearby candidates from the original peak-gene map, not validated regulatory assignments.

| Peak | Nearby gene | Ratio | Exact paired p | Middle-only q |
|---|---|---:|---:|---:|
| chr9:98916421-98917309 | COL15A1 | 2.090909 | 0.000936 | 0.031832 |
| chr10:86674717-86675570 | LDB3 | 1.937500 | 0.031539 | 0.536171 |
| chr3:52446715-52447619 | TNNC1 | 1.764706 | 0.065994 | 0.747932 |

COL15A1 has both the largest ratio and smallest p in this 34-peak screen. One LDB3 peak also has raw p<0.05. Only COL15A1 has middle-only q<0.05. A second LDB3 peak has ratio 1.625 and p=0.3833. The two peaks are retained separately; the earlier gene-region score is a different endpoint.

### Statistical record

The middle programme was prioritized after reviewing the temporal-module results. The 34-peak q is an exploratory within-middle adjustment; it does not replace the full 229-peak screening q or account for selecting this phase from the observed module results. COL15A1 has q=0.071465 in the full 229-peak family. Both columns remain in `middle_peak_followup.tsv`; the figure presents only the middle analysis. Programme q, middle-peak q and earlier gene-region q have different test units and families.

## Functional interpretation of leading loci

- **COL15A1:** extracellular-matrix collagen XV. Knockout mice develop skeletal myopathy and microvascular defects, supporting muscle/connective-tissue stability. This does not demonstrate that COL15A1 overexpression promotes myoblast differentiation or fusion. [Eklund et al., 2001](https://pmc.ncbi.nlm.nih.gov/articles/PMC14731/).
- **LDB3 / ZASP / Cypher:** Z-disc scaffold and structural integrity during muscle contraction. Cypher ablation causes severe congenital myopathy and disrupted Z-lines. The functional evidence is strongest for muscle structure and maintenance, rather than initiating myogenic differentiation. [Zhou et al., 2001](https://pmc.ncbi.nlm.nih.gov/articles/PMC2198871/).
- **TNNC1:** calcium-sensitive troponin component expressed in slow skeletal and cardiac muscle; associated with contractile function. It does not reach raw p<0.05 in this single-peak test. [Functional study](https://pubmed.ncbi.nlm.nih.gov/28473771/).
- **ITGB1**, also present among the 30 nearby genes, has direct experimental evidence for myoblast fusion and sarcomere assembly, but its middle ATAC peak has p=0.6076 and q=1. The functional literature alone does not make it a supported COQ8A-associated peak. [Schwander et al., 2003](https://pubmed.ncbi.nlm.nih.gov/12737803/).

These are literature-based functional descriptions, not a newly computed pathway-enrichment result.

## Reproduction

```bash
python scripts/29_plot_middle_followup.py --temporal results/temporal --figures figures/temporal
```

The step is included in `run_temporal.py`. The released table includes every middle peak in all four target settings. Private C2C12 expression data and its local cross-check are kept outside this public repository.

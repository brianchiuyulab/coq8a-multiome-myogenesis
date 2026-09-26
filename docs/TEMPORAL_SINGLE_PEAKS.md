# COQ8A association at individual externally timed peaks

## Question and analysis population

Does the modest average change across early-opening regions conceal larger changes at individual sites? All **191 early, 34 middle and 4 late peaks** are tested together. Membership comes from the independent GSE109828 differentiation time course, using the fixed 221-gene candidate universe, 2-kb external marker-promoter definition and 50% reciprocal overlap. No COQ8A result is used to remove a peak from this screen. External region selection is detailed in [Temporal Methods](TEMPORAL_METHODS.md).

GSE208248 provides paired RNA/ATAC nuclei in four libraries from two source lines, each measured in stem and differentiated conditions. Existing joint QC, within-library matching and peak mappings are unchanged. Each of the four COQ8A/TSS settings tests the identical 229 sites. The main display retains COQ8A >=3 versus exactly 1 RNA UMI and TSS enrichment >=3: 201 high/low pairs (46, 74, 58 and 23 pairs in GSM6339597, GSM6339599, GSM6339601 and GSM6339603).

## Endpoint and statistics

For each peak, accessibility is the binary presence of ATAC counts in each nucleus. High and low open fractions are the numbers of open nuclei divided by their respective equal group sizes. Fold change is the ratio of those fractions, with no pseudocount. A zero denominator yields infinity when the numerator is positive, or an undefined value when both counts are zero; counts and percentage-point differences remain available.

Each matched pair contributes to high-only, low-only, both-open or neither-open counts. A two-sided exact paired binomial test (exact McNemar test) compares high-only and low-only counts. All 229 p values, including p=1 for peaks without discordant pairs, undergo Benjamini-Hochberg adjustment together **within each setting**. There is no minimum target-open-count filter and no adjustment across sensitivity settings. Early, middle and late labels do not divide the test family. Four-library effect directions and all raw p values are retained.

These are matched-nucleus association tests in the observed two-line dataset. The 201 pairs are not 201 independent source donors. Gene labels denote nearby genes in the original candidate map, not demonstrated regulatory targets. This analysis stops at testing peak accessibility against COQ8A expression.

## Results

COQ8A >=3 versus 1 UMI, TSS >=3:

| Peak (hg38) | Nearby candidate gene | External phase | Low open | High open | Ratio | Exact paired p | BH q, 229 peaks | Positive libraries |
|---|---|---|---:|---:|---:|---:|---:|---:|
| chr16:56341742-56342652 | GNAO1 | Early | 23/201 (11.44%) | 48/201 (23.88%) | 2.087 | 0.000346 | 0.0715 | 4/4 |
| chr16:1152563-1153401 | CACNA1H | Early | 22/201 (10.95%) | 5/201 (2.49%) | 0.227 | 0.000911 | 0.0715 | 0/4 |
| chr9:98916421-98917309 | COL15A1 | Middle | 22/201 (10.95%) | 46/201 (22.89%) | 2.091 | 0.000936 | 0.0715 | 4/4 |
| chr17:39704708-39705562 | TCAP | Early | 9/201 (4.48%) | 21/201 (10.45%) | 2.333 | 0.028959 | 0.8130 | 4/4 |
| chr3:138356441-138357333 | MRAS | Early | 43/201 (21.39%) | 27/201 (13.43%) | 0.628 | 0.029305 | 0.8130 | 0/4 |
| chr11:19201752-19202603 | CSRP3 | Early | 7/201 (3.48%) | 17/201 (8.46%) | 2.429 | 0.030884 | 0.8130 | 2/4 |
| chr10:86674717-86675570 | LDB3 | Middle | 16/201 (7.96%) | 31/201 (15.42%) | 1.938 | 0.031539 | 0.8130 | 3/4 |
| chr14:103595381-103596012 | CKB | Early | 12/201 (5.97%) | 3/201 (1.49%) | 0.250 | 0.035156 | 0.8130 | 1/4 |

There are ten peaks with nominal p<0.05: eight early (five increases, three decreases) and two middle (both increases). None reaches q<0.05 across the 229-site family. The three lowest-p sites have q=0.0715; this numerical result is reported without replacing the q<0.05 criterion. No setting yields a q<0.05 single-peak result in this family.

GNAO1 and COL15A1 are the leading positive individual-site associations. Both have positive effects in all four libraries under the main setting and all four target settings, although effect magnitude depends on the COQ8A contrast:

| Nearby gene | >=2 vs 1, TSS>=2 | >=2 vs 1, TSS>=3 | >=3 vs 1, TSS>=2 | >=3 vs 1, TSS>=3 |
|---|---:|---:|---:|---:|
| GNAO1 | 1.259 | 1.229 | 2.042 | 2.087 |
| COL15A1 | 1.164 | 1.155 | 2.043 | 2.091 |

**The early module average of 1.064 does not imply that all early sites change by 1.064.** Site-level increases and decreases coexist. Focusing exclusively on the 34 middle sites would miss the early GNAO1 candidate. These observations establish local heterogeneity; they do not establish that the remaining early sites are background noise.

### Relation to previous summaries

The earlier early-module 1.064 and middle-module 1.223 ratios are unchanged. Summing the individual-peak counts reproduces all 16 module ratios (four phases/sets across four settings). Earlier gene-neighborhood results used paired t tests of region scores and separate gene-level families (117 early genes or 30 middle genes). This screen uses exact binary-peak tests and one 229-peak family. In particular, the previous COL15A1 gene-screen p=0.000600, q=0.0180 and the present single-peak p=0.000936, q=0.0715 describe different tests/families. Their underlying high/low counts and 2.091-fold ratio agree.

## Figure

[Individual temporal peak screen](../figures/temporal/Figure_temporal_single_peaks.png) ([vector PDF](../figures/temporal/Figure_temporal_single_peaks.pdf)):

- **A:** All 229 peaks; horizontal position is the high-minus-low open fraction in percentage points, vertical position is the raw exact p value on a negative log10 scale. Colors indicate external timing. Dashed line: nominal p=0.05. Overlapping dots can occupy the same coordinate.
- **B:** Eight lowest raw p values in the main setting, retaining both effect directions. Points show actual low/high open fractions; text shows ratio, p and q. This is a descriptive display of the complete screen, with no retesting of the displayed subset.
- **C:** The same eight sites under all four COQ8A/TSS settings. Cell text is the high/low ratio; color represents log2 ratio centered at zero.
- **D:** Site-level percentage-point differences for each library at the main setting. Red denotes higher accessibility in COQ8A-high nuclei; blue denotes lower accessibility. Four libraries represent two source lines measured in two states.

## Reproduction and source data

```bash
python scripts/27_temporal_single_peaks.py --tables results/tables --temporal results/temporal
python scripts/28_plot_temporal_single_peaks.py --temporal results/temporal --out figures/temporal
```

Both steps are included in `run_temporal.py`. The single-peak screen reuses the paired binary counts produced from the target H5 matrices by `05_peak_effects.py`; it does not reuse that script's count-filtered q values. The script checks unique memberships, all library/setting combinations, pair counts, discordant-count consistency and exact reconciliation to the previously generated module fractions and ratios. GitHub Actions reruns the screen.

Files in `results/temporal/`:

- `single_peak_effects.tsv`: all 916 peak/setting results, raw counts, fractions, fold changes, p and q.
- `single_peak_by_library.tsv.gz`: all 3,664 peak/library/setting rows and source/state identities.
- `single_peak_screen_summary.tsv`: phase-level counts of nominal and adjusted results in every setting.
- `single_peak_audit.json`: coverage and reconciliation checks.

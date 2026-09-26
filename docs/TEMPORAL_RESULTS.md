# Results: externally timed myogenic accessibility

## Individual sites across all opening phases

The [complete single-peak screen](TEMPORAL_SINGLE_PEAKS.md) tests all 191 early, 34 middle and 4 late peaks together. At COQ8A>=3 versus 1 UMI/TSS>=3, an early GNAO1-neighborhood site is 2.087-fold more accessible (exact paired p=0.000346), and a middle COL15A1-neighborhood site is 2.091-fold (p=0.000936). Both increase in all four libraries. Across all 229 peaks, both have BH q=0.0715; no site reaches q<0.05. An early CACNA1H-neighborhood site decreases to 0.227-fold (p=0.000911, q=0.0715). Thus the modest early-module average conceals larger local changes of both directions. The module summaries below are unchanged.

## Main result

Within the fixed 221-gene universe, independent external differentiation timing defines a **34-peak, 30-gene middle-opening module**. In GSE208248, COQ8A>=3 versus exactly 1 UMI at TSS>=3 gives **1.222930-fold accessibility**, from **6.8920% to 8.4284%** open peaks per nucleus (+1.5364 percentage points), paired **p=0.001785**, within-setting programme **BH q=0.017848**. There are 201 matched pairs, four libraries and two source lines; three libraries have positive effects.

This module reaches half its external opening gain between **24 and 48 h**. It is not the earliest-opening programme. Early-by-24-h regions are weaker; earliest timing quantiles show no amplified association.

## Complete module comparison

Primary external definition: marker-promoter radius 2 kb, reciprocal overlap>=50%. Target: COQ8A>=3 vs 1, TSS>=3.

| External module | Peaks | High/low ratio | Paired p |
|---|---:|---:|---:|
| All mapped | 3,324 | 1.0061 | 0.3287 |
| All opening | 229 | 1.0808 | 0.0180 |
| Early, by 24 h | 191 | 1.0638 | 0.0683 |
| Middle, 24-48 h | 34 | **1.2229** | **0.001785** |
| Late, after 48 h | 4 | 1.0345 | 0.7398 |
| Low-change profile | 470 | 1.0005 | 0.9655 |
| Closing | 181 | 0.9608 | 0.0699 |
| Earliest 20% | 48 | 1.0031 | 0.9358 |
| Earliest 25% | 58 | 1.0067 | 0.8570 |
| Earliest 33% | 76 | 1.0148 | 0.6816 |

The middle-minus-low-change fold difference is 0.2225 (paired bootstrap CI 0.0853-0.3720), scaled paired p=0.00197. Middle-minus-early is 0.1591 (CI 0.0280-0.3042), p=0.01817. This directly tests differing module effects. The four-peak late set cannot represent all late myogenic regulation.

## Sensitivity

Identical 34-peak membership is used in all four target settings:

| COQ8A contrast | TSS | Pairs | Ratio | Paired p | Positive libraries |
|---|---|---:|---:|---:|---:|
| >=2 vs 1 | >=2 | 1,020 | 1.1014 | 0.000791 | 4/4 |
| >=2 vs 1 | >=3 | 958 | 1.1065 | 0.000676 | 4/4 |
| >=3 vs 1 | >=2 | 212 | 1.2140 | 0.001550 | 3/4 |
| >=3 vs 1 | >=3 | 201 | 1.2229 | 0.001785 | 3/4 |

At target >=3 vs 1/TSS>=3, changing the external definition gives:

| External promoter radius | Reciprocal overlap | Middle peaks | Ratio | Paired p |
|---|---|---:|---:|---:|
| 2 kb | 50% | 34 | 1.2229 | 0.001785 |
| 2 kb | 25% | 45 | 1.1732 | 0.004629 |
| 1 kb | 50% | 33 | 1.1481 | 0.023937 |
| 1 kb | 25% | 44 | 1.1185 | 0.039786 |

Direction persists, but **1.22 is not invariant to external classification**. Complete q values are retained in `temporal_effect_grid.tsv`.

## Cell state

| State/library | Pairs | Ratio | Paired p |
|---|---:|---:|---:|
| Undifferentiated, combined | 104 | **1.3440** | **0.001190** |
| Differentiated, combined | 97 | 1.0860 | 0.3560 |
| Line 1 undifferentiated | 46 | 1.2857 | 0.04085 |
| Line 2 undifferentiated | 58 | 1.4103 | 0.01250 |
| Line 1 differentiated | 74 | 1.1611 | 0.1528 |
| Line 2 differentiated | 23 | 0.9306 | 0.6811 |

The source-adjusted state-interaction p=0.0556 (ten-module q=0.556). The larger undifferentiated estimate is exploratory, not a demonstrated state difference. These are culture states, not young/old groups.

## Candidate genes

The 34 middle-opening peaks map to these 30 candidate neighborhoods:

AK1, ANKRD2, APOD, CACNA1H, CAMK2B, CASQ1, COL15A1, ENO3, HRC, ITGB1, LARGE1, LDB3, MEF2B, MRAS, MYH2, MYH3, MYL7, MYO1C, MYOM1, PSEN2, PVALB, RB1, REEP1, RIT1, RYR1, SLC6A8, SVIL, TNNC1, TNNT2, TPD52L1.

There is no MYOD1 neighborhood in this module. Within the complete 30-gene screen, COL15A1's one peak has ratio 2.0909, p=0.000600, q=0.01801; LDB3's two peaks have ratio 1.8333, p=0.01970, q=0.29545. All genes, including negative effects, remain in the output. Proximity assignments are candidate gene neighborhoods, not established regulatory targets.

## Additional RNA QC

All matched nuclei have >2,300 detected RNA genes. Mitochondrial UMI medians vary by library from approximately 10.2% to 17.7%; the inherited baseline did not impose a mitochondrial-percentage cutoff. A universal 5% gate therefore cannot be described as already passed.

Removing Scrublet predicted doublets retains 192 primary pairs: **1.2045**, p=0.00468. Removing the highest 1%, 2.5%, 5% or 10% Scrublet-score tails gives **1.229-1.235**, p=0.00157-0.00682. The signal is not confined to the high-doublet-score tail.

Mitochondrial 5/10/15/20% gates are explicit additional sensitivities, preserving pair identities. At 5%, no primary pairs remain; at 10%, only 33 remain (1.2273, p=0.318). At 15%, 153 pairs remain (1.1909, p=0.01375); at 20%, 190 remain (1.2374, p=0.00120). Complete counts and outcomes are in `rna_qc_sensitivity.tsv`; the 5/10% branches do not establish statistical persistence.

## Figures

1. **Figure_external_timing**: rows are external opening peaks, columns actual sampling times. Color is scaled within each row and represents timing rather than absolute accessibility; median profiles explain the classification.
2. **Figure_COQ8A_temporal_sensitivity**: fixed temporal modules by four target settings, displaying actual ratios and paired p values. The middle-opening row is strongest; earliest-quantile rows remain near one.
3. **Figure_temporal_states_and_genes**: paired differences/95% intervals by cell state, followed by every gene in the middle module. Dot size is peak count; teal gene dots have within-module q<0.05.

## Interpretation

The results support exploring high COQ8A in association with a **specific differentiation-associated temporal accessibility module**, with a larger estimate in undifferentiated cultures. They do not support preferential association with the very earliest-opening elements. This independent external-definition route is distinct from the previous six-peak MYOD1 selection.

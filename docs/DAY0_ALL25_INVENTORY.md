# All 25 Day-0 neighborhoods: sensitivity inventory

All 336 configurations were attempted; 312 were evaluable. This report summarizes existing results, without rerunning or changing the grid.

For each gene and direction, the representative is the smallest regional p among settings whose top peak changes in that direction in both sources. FC and peak p refer to that top peak, not the whole neighborhood. Regional p accounts for searching peaks within the neighborhood; q25 is BH across the 25 neighborhoods in that setting and direction. The separately reported minimum q can come from a different setting.

Counts of settings are sensitivity descriptions, not independent replications. Settings reuse nuclei. Best-setting p/q do not account for selection across the grid. Missing FC means a zero low-group detection denominator, not a measured infinite biological effect. Different directions or settings may select different peaks. Overlapping neighborhoods (including MYF5/MYF6) can share a peak.

Full machine-readable results, including both peak and regional statistics, are in `results/day0_sensitivity_grid/all25_directional_summary.tsv`.

## Opening

| Neighborhood | Config / population / COQ8A contrast / TSS / caliper / matching | Pairs | Top peak FC | Peak p | Region p | Region q25 | Settings p<.05 / q<.1 |
|---|---|---:|---:|---:|---:|---:|---:|
| ADAM12 | D064; Myogenic-marker+; 3plus_vs_1; 2; 0.05; depth | 80 | 4.29 | 3.395e-05 | 9.3e-05 | 0.002325 | 40 / 16 |
| MYOD1 | D036; All D0; 2plus_vs_zero; 2; 0.2; depth | 528 | 1.51 | 5.365e-05 | 0.00025 | 0.00625 | 70 / 44 |
| CSRP3 | D042; All D0; 2plus_vs_le1; 2; 0.1; depth | 525 | 1.75 | 0.0003096 | 0.00055 | 0.01375 | 50 / 16 |
| MYH9 | D097; Myogenic-marker+; 2plus_vs_le1; 2; 0.05; depth_state | 324 | 1.95 | 9.264e-05 | 0.001 | 0.025 | 37 / 9 |
| CACNA1H | D034; All D0; 2plus_vs_zero; 2; 0.1; depth | 523 | 2.47 | 0.0003069 | 0.00205 | 0.01708 | 18 / 8 |
| BOC | D239; Myogenic-marker+; 3plus_vs_1; 3; 0.3; depth_state | 91 | NA (low=0) | 0.0009766 | 0.00235 | 0.05875 | 16 / 4 |
| TGFB1 | D108; Myogenic-marker+; positive_CP10k_Q25; 2; 0.2; depth | 88 | 15 | 0.0005188 | 0.0026 | 0.065 | 32 / 4 |
| CAV3 | D224; Myogenic-marker+; 2plus_vs_1; 3; 0.05; depth | 373 | 4.8 | 0.0005461 | 0.003021 | 0.07552 | 37 / 14 |
| CDON | D222; All D0; positive_CP10k_Q25; 3; 0.3; depth | 128 | 14 | 0.0009766 | 0.0044 | 0.11 | 12 / 0 |
| IGF1 | D015; All D0; 3plus_vs_1; 2; 0.3; depth_state | 97 | 2.67 | 0.004077 | 0.0066 | 0.0825 | 10 / 2 |
| RB1 | D081; Myogenic-marker+; detected_vs_zero; 2; 0.05; depth_state | 1758 | 1.71 | 0.0006398 | 0.00685 | 0.1712 | 10 / 0 |
| NOTCH1 | D122; PAX7+; 3plus_vs_1; 2; 0.1; depth | 13 | 4.5 | 0.01562 | 0.0075 | 0.1875 | 16 / 0 |
| MAPK14 | D316; PAX7+; 2plus_vs_zero; 3; 0.2; depth | 89 | 5.5 | 0.003906 | 0.00865 | 0.2162 | 13 / 0 |
| NOS1 | D034; All D0; 2plus_vs_zero; 2; 0.1; depth | 523 | 11 | 0.006348 | 0.01235 | 0.06175 | 28 / 8 |
| MYF5 | D287; PAX7+; 2plus_vs_1; 3; 0.3; depth_state | 71 | NA (low=0) | 0.007812 | 0.01273 | 0.1592 | 27 / 8 |
| MYF6 | D287; PAX7+; 2plus_vs_1; 3; 0.3; depth_state | 71 | NA (low=0) | 0.007812 | 0.01273 | 0.1592 | 27 / 8 |
| NEO1 | D173; All D0; 2plus_vs_1; 3; 0.2; depth_state | 449 | 3.29 | 0.003719 | 0.01545 | 0.2094 | 6 / 0 |
| IFRD1 | D282; PAX7+; 2plus_vs_1; 3; 0.1; depth | 73 | 2.14 | 0.007 | 0.01725 | 0.4312 | 4 / 0 |
| MAPK12 | D211; All D0; 2plus_vs_le1; 3; 0.1; depth_state | 458 | 1.28 | 0.001935 | 0.0189 | 0.4725 | 6 / 0 |
| MYOG | D221; All D0; positive_CP10k_Q25; 3; 0.2; depth_state | 65 | 1.67 | 0.005223 | 0.02405 | 0.3006 | 12 / 0 |
| KCNH1 | D120; PAX7+; 3plus_vs_1; 2; 0.05; depth | 10 | NA (low=0) | 0.0625 | 0.0296 | 0.37 | 4 / 0 |
| IGFBP3 | D183; All D0; 3plus_vs_1; 3; 0.3; depth_state | 97 | 3.25 | 0.01172 | 0.03995 | 0.2497 | 2 / 0 |
| ITGB1 | D107; Myogenic-marker+; positive_CP10k_Q25; 2; 0.1; depth_state | 28 | NA (low=0) | 0.03125 | 0.06835 | 0.8392 | 0 / 0 |
| MEF2C | D234; Myogenic-marker+; 3plus_vs_1; 3; 0.1; depth | 96 | NA (low=0) | 0.01562 | 0.0691 | 0.4645 | 0 / 0 |
| ANKRD2 | D030; All D0; detected_vs_zero; 2; 0.3; depth | 3008 | 1.17 | 0.005705 | 0.0697 | 0.2904 | 0 / 0 |

## Closing

| Neighborhood | Config / population / COQ8A contrast / TSS / caliper / matching | Pairs | Top peak FC | Peak p | Region p | Region q25 | Settings p<.05 / q<.1 |
|---|---|---:|---:|---:|---:|---:|---:|
| NOTCH1 | D037; All D0; 2plus_vs_zero; 2; 0.2; depth_state | 494 | 0.397 | 8.14e-06 | 6.45e-05 | 0.001612 | 85 / 34 |
| MYOG | D325; PAX7+; 2plus_vs_le1; 3; 0.2; depth_state | 80 | 0.143 | 4.005e-05 | 0.000132 | 0.0033 | 15 / 12 |
| ITGB1 | D248; Myogenic-marker+; detected_vs_zero; 3; 0.05; depth | 2429 | 0.62 | 3.902e-05 | 0.00025 | 0.00625 | 36 / 24 |
| CAV3 | D224; Myogenic-marker+; 2plus_vs_1; 3; 0.05; depth | 373 | 0.45 | 0.0001123 | 0.0004545 | 0.01136 | 16 / 8 |
| MYH9 | D004; All D0; 2plus_vs_1; 2; 0.2; depth | 526 | 0.39 | 5.571e-05 | 0.00065 | 0.01625 | 56 / 22 |
| KCNH1 | D083; Myogenic-marker+; detected_vs_zero; 2; 0.1; depth_state | 2223 | 0.737 | 0.0001633 | 0.00105 | 0.02625 | 28 / 10 |
| BOC | D038; All D0; 2plus_vs_zero; 2; 0.3; depth | 529 | 0.233 | 0.0001911 | 0.00125 | 0.02812 | 47 / 16 |
| CSRP3 | D248; Myogenic-marker+; detected_vs_zero; 3; 0.05; depth | 2429 | 0.854 | 0.0005406 | 0.00175 | 0.02187 | 20 / 16 |
| IGFBP3 | D222; All D0; positive_CP10k_Q25; 3; 0.3; depth | 128 | 0.4 | 0.001431 | 0.00405 | 0.08062 | 13 / 2 |
| IFRD1 | D097; Myogenic-marker+; 2plus_vs_le1; 2; 0.05; depth_state | 324 | 0.766 | 0.001091 | 0.005 | 0.125 | 21 / 0 |
| MYOD1 | D325; PAX7+; 2plus_vs_le1; 3; 0.2; depth_state | 80 | 0.0909 | 0.001953 | 0.005649 | 0.07061 | 12 / 4 |
| MAPK14 | D283; PAX7+; 2plus_vs_1; 3; 0.1; depth_state | 33 | 0.474 | 0.01294 | 0.0066 | 0.165 | 8 / 4 |
| MAPK12 | D221; All D0; positive_CP10k_Q25; 3; 0.2; depth_state | 65 | 0.286 | 0.001953 | 0.00855 | 0.2137 | 14 / 0 |
| TGFB1 | D030; All D0; detected_vs_zero; 2; 0.3; depth | 3008 | 0.582 | 0.0009062 | 0.00865 | 0.07208 | 26 / 4 |
| CACNA1H | D255; Myogenic-marker+; detected_vs_zero; 3; 0.3; depth_state | 2494 | 0.765 | 0.0006369 | 0.00895 | 0.2237 | 17 / 0 |
| MEF2C | D232; Myogenic-marker+; 3plus_vs_1; 3; 0.05; depth | 80 | 0.167 | 0.006348 | 0.01536 | 0.3839 | 8 / 0 |
| CDON | D153; PAX7+; 2plus_vs_le1; 2; 0.05; depth_state | 46 | 0.214 | 0.007385 | 0.0155 | 0.3875 | 10 / 0 |
| MYF6 | D119; PAX7+; 2plus_vs_1; 2; 0.3; depth_state | 71 | 0.125 | 0.01562 | 0.0244 | 0.2033 | 10 / 0 |
| MYF5 | D119; PAX7+; 2plus_vs_1; 2; 0.3; depth_state | 71 | 0.125 | 0.01562 | 0.0244 | 0.2033 | 10 / 0 |
| ADAM12 | D171; All D0; 2plus_vs_1; 3; 0.1; depth_state | 347 | 0.619 | 0.001993 | 0.02485 | 0.3537 | 2 / 0 |
| NOS1 | D049; All D0; positive_CP10k_Q25; 2; 0.05; depth_state | 14 | 0 | 0.0625 | 0.0299 | 0.7112 | 2 / 0 |
| NEO1 | D033; All D0; 2plus_vs_zero; 2; 0.05; depth_state | 351 | 0.154 | 0.007385 | 0.03205 | 0.2671 | 2 / 0 |
| RB1 | D024; All D0; detected_vs_zero; 2; 0.05; depth | 2876 | 0.5 | 0.003836 | 0.0382 | 0.191 | 6 / 0 |
| ANKRD2 | D157; PAX7+; 2plus_vs_le1; 2; 0.2; depth_state | 80 | 0.2 | 0.007812 | 0.0406 | 0.2537 | 2 / 0 |
| IGF1 | D104; Myogenic-marker+; positive_CP10k_Q25; 2; 0.05; depth | 36 | 0 | 0.0625 | 0.0477 | 0.7245 | 2 / 0 |

## MYOD1 interpretation

At TSS>=3 in all Day-0 cultures, COQ8A>=2 versus zero and depth-only matching gives the same peak chr11:17653349–17654252 at every caliper (0.05/0.10/0.20/0.30), FC 1.495/1.490/1.510/1.510 and regional q25 0.0100/0.01437/0.00875/0.007962. Additional state matching weakens regional q to 0.7978–0.9028. This is a worthwhile state-associated opening candidate, not proof of a state-independent effect or a confirmed MYOD1 RNA target. The same coordinates are not the historical chr11:17649919–17650798 focal peak.

Opening q<0.1 occurs in at least one concordant setting for 12 neighborhood labels: ADAM12, MYOD1, CSRP3, CACNA1H, MYH9, CAV3, BOC, NOS1, TGFB1, IGF1, MYF5 and MYF6. These are not necessarily 12 independent loci. Large FC values should be read alongside the detection counts in the complete peak table.

ANKRD2, MEF2C and ITGB1 have no concordant regional opening p<0.05 in this grid. That does not establish no chromatin association: ITGB1 has closing evidence. No gene should be described as universally unchanged simply from failure of its opening test.

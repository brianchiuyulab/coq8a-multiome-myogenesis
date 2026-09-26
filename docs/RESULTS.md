# Results

## Analysis population

The main comparison uses COQ8A RNA ≥3 versus 1 UMI, TSS enrichment ≥3, and
201 pairs matched within four libraries from two source lines. The biological
universe, external temporal sets and functional membership are fixed across
all four sensitivity settings.

## Temporal modules

| External class | Peaks | Accessibility FC, high/low | Paired p |
|---|---:|---:|---:|
| Early | 191 | 1.064 | 0.0683 |
| Middle | 34 | 1.223 | 0.00178 |
| Late | 4 | 1.034 | 0.740 |
| Closing | 181 | 0.961 | 0.0699 |

The middle module has BH q = 0.0178 within the ten temporal programme tests.
The module effect is an average across all 34 member regions; it does not
mean that each region individually passes FDR. Figure 1 shows the external
definition and the matched target comparison. Supplementary Figure S2 shows
every TSS/count setting.

## Functional dynamic regions

External differentiation, fusion and positive/negative regulation annotations
define 54 distinct regions: 24 early, five middle and 25 closing. These contain
both positive and negative COQ8A associations. None of the eight aggregate
functional modules has nominal p < 0.05 in the four evaluated settings.
Figure 2 displays the entire set without target-effect filtering.

## Follow-up loci

ATAC p values are exact paired tests and q values are BH adjusted across 54
functional peaks. RNA p values are paired t tests and q values are BH adjusted
across 221 genes. RNA Δ is on the loge(1 + CP10k) scale.

| Nearby gene | ATAC FC | ATAC p | ATAC q | RNA Δ | RNA p | RNA q |
|---|---:|---:|---:|---:|---:|---:|
| CSRP3 | 2.429 | 0.0309 | 0.711 | +0.0716 | 0.0121 | 0.0743 |
| CAV3 | 1.765 | 0.0533 | 0.711 | +0.0838 | 0.0432 | 0.163 |
| MYOD1 | 1.706 | 0.0730 | 0.711 | +0.0782 | 0.0510 | 0.185 |
| CACNA1H | 0.227 | 0.000911 | 0.0492 | +0.00316 | 0.196 | 0.402 |

CSRP3, CAV3 and MYOD1 have positive pooled ATAC and RNA directions, with
different consistency across libraries. Their displayed ATAC peaks do not
pass the 54-peak FDR threshold. CACNA1H has decreased accessibility in all four
libraries and passes that threshold in the main setting. These observations
nominate locus-specific follow-up; they are not four confirmed regulatory links.

The expanded peak–RNA analysis yields no FDR-supported link to these four
genes in the main setting. CAV3 and MYOD1 have eligible positive but
nonsignificant partial correlations; CSRP3 and CACNA1H lack sufficient RNA
detection for the required per-library link tests. Genomic neighborhood labels
in figures therefore remain candidate assignments.

## Interpretation

The results support a positive COQ8A association with the externally defined
middle-opening module and heterogeneous changes at individual myogenic
regions. They do not show uniform opening of the functional set. The public
comparison supplies association evidence for a metabolism–chromatin hypothesis;
it does not measure the intervening metabolites or perturb COQ8A.

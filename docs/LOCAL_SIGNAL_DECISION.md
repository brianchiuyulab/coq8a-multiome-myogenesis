# Local-signal analysis and candidate decision

## Decision

The closest coherent exploratory route is **CSRP3**. The old culture data have
a direction-specific gene-neighborhood opening signal in the COQ8A>=2 vs 1
setting; RNA also increases under that same setting. The new adult-muscle data
provide a weaker CSRP3 promoter/RNA association under a different COQ8A
definition. This is partial support, not a fully significant chain or a
same-setting independent replication. MYOD1's large sensitivity FC is an
ATAC-only lead with sparse nucleus support; it is not the strongest integrated
candidate simply because its FC is the largest.

## Task 1: test local signals within each of the 25 neighborhoods

Script 84 tests whether any peak in each frozen gene neighborhood becomes more
accessible, without averaging all peaks. The 25 memberships are unchanged.
For each peak, D is the vector of within-pair differences in binary accessibility.
Its signed statistic is sum(D)/sqrt(sum(D^2)), or zero when all differences are
zero. A gene's opening statistic is the maximum signed statistic across its
peaks, truncated at zero. The two-sided sensitivity uses the maximum absolute
statistic. All eligible peaks participate, including peaks without differences.

High/low labels are swapped independently within matched nucleus pairs. The same
swap is applied across all peaks, preserving peak dependence and library pairing.
The analysis uses 2,000,000 Monte Carlo permutations with seed 20260928;
p=(exceedances+1)/(B+1), followed by BH over all 25 gene-neighborhood hypotheses
within each contrast/direction. An initial 20,000-permutation estimate was near
q=0.10; the final higher precision was used to quantify Monte Carlo uncertainty,
not to choose a different seed. Approximate Monte Carlo uncertainty intervals
are saved in the source table.

This changes the tested hypothesis. A gene-neighborhood q is not a single-peak
q and cannot be attached to the most favorable peak as if that peak passed
correction. These remain paired-nucleus association tests from two biological
source lines, not 958 independent donors. The method and settings are a
follow-up after prior exploration, not a claim of preregistration.

### Best-supported exploratory setting

**TSS>=3; COQ8A>=2 versus exactly 1 UMI; 958 matched pairs.**

CSRP3 has 10 candidate peaks. Its opening-region test gives **p=0.003825,
q25=0.095625**. The Monte Carlo p interval is approximately 0.003740–0.003911.
The two-sided region test gives p=0.007730, q25=0.193250. Thus the finding
depends on testing the directional opening hypothesis. In the retained >=3 vs 1
contrast with 201 pairs, its opening-region q is 0.6718. Both contrasts and
both direction definitions are reported, without relabeling the favorable
sensitivity as the original primary analysis.

The strongest CSRP3 peak under the >=2 vs 1 contrast is
**chr11:19218592–19219518**: high 139/958, low 92/958, **FC=1.51087**,
two-sided paired-binomial p=0.001190, single-peak q565=0.3361.
All four libraries have higher accessibility:

| Library | High open | Low open | Matched pairs |
|---|---:|---:|---:|
| GSM6339597 | 32 | 23 | 211 |
| GSM6339599 | 42 | 29 | 326 |
| GSM6339601 | 58 | 38 | 296 |
| GSM6339603 | 7 | 2 | 125 |

In the same 958 pairs, CSRP3 RNA mean CP10k is 0.142624 versus 0.058631:
**FC=2.43259**, paired log1p-CP10k p=0.012362, q25=0.104054. RNA direction
is positive in two of the four individual libraries. The sparse-expression
1% link analysis gives partial r=0.0226 for this peak and CSRP3; after additionally
conditioning on the myogenesis state score, r is 0.00248 and p=0.695. This is
not evidence for a strong state-independent cis relationship. The primary 5%
detection criterion does not provide an evaluable CSRP3 link.

Owner-provided C2C12 time-course results were cross-referenced privately in
script 83. Their direction supports retaining CSRP3 for exploration; those
time contrasts are not COQ8A overexpression comparisons and are not published
here.

## Task 2: retain partial chains and compare every candidate consistently

Script 85 evaluates all matched settings from the complete native 645-peak
analysis. The descriptive evidence steps are:

| Evidence | Unique peaks | Unique scope genes |
|---|---:|---:|
| ATAC increases with raw p<0.05 in at least one matched setting | 310 | 25 |
| Also RNA increases with raw p<0.05 in that same setting | 20 | 10 |
| Also the peak–RNA correlation has positive direction | 9 | 6 |
| All three individual p values <0.05 | 0 | 0 |

These settings are correlated analyses of the same data, not independent
replications. The final row measures completeness of nominal statistical
support; it is **not an eligibility gate for exploratory candidates**. Missing
RNA/link estimates are retained as missing rather than interpreted as zero.
All 36 ATAC-plus-RNA nominal setting/peak/gene rows are supplied, along with
unrestricted transcript annotation. Gene-neighborhood labels alone do not
identify regulated RNA.

### CSRP3: coherent direction, modest new-data effect

In Fast Pre nuclei, author QC, COQ8A detected versus zero, caliper 0.20
(C138), six donors and 3,417 pairs:

- Native promoter chr11:19201660–19202638 overlaps a CSRP3 TSS.
- ATAC **FC=1.13624, p=0.02485, native-union q=0.92092**.
- RNA **FC=1.13710, p=0.002614, q25=0.05751**.
- Donor-combined partial peak–RNA **r=0.02788, p=0.1287**.

With the 95th-percentile depth cap and the same detected/zero definition
(C166), ATAC FC=1.12624, p=0.04637, q=0.98432; RNA FC=1.14613,
p=0.001887, q=0.04151. Under author QC, >=2 vs 1, caliper 0.10
(C113), ATAC FC=1.0806, p=0.3224 and RNA FC=1.0195, p=0.7229.
The new dataset therefore offers a weak same-gene association under a different
expression definition, not successful replication of the old >=2 vs 1 setting.
The new promoter is not the old chr11:19218592–19219518 region-test maximum.

Other nominal ATAC/RNA pairs include MYH9, MEF2C, ANKRD2, MAPK12 and RB1.
Their complete results and nearest-transcript annotations are in the candidate
table. For example, one ANKRD2-neighborhood signal actually overlaps a HOGA1
TSS, and MYH9's largest nominal effect has a near-zero peak–RNA correlation.
Such findings are retained but do not become established target-gene links.

### MYOD1: two q values correspond to two settings

For Satellite Post, COQ8A>=2 vs 1:

| Matching caliper | Donors / pairs | ATAC FC | Raw p | Native-union q |
|---|---|---:|---:|---:|
| 0.20 (C058) | 4 / 11 | 10.4266 | 0.001506 | 0.03765 |
| 0.30 (C059) | 6 / 16 | 4.4280 | 0.01598 | 0.5034 |

The first result is not invalidated by the second; the comparison shows
sensitivity to matching and available counts. The 22 C058 nuclei have no
detected MYOD1 RNA. It remains a sparse promoter-accessibility lead, rather
than the best integrated explanation of the user's phenotype. The chr16
UBE2I-promoter opening result likewise has stronger ATAC statistics but no
established matching RNA response.

## Reproduction

```bash
python scripts/84_gene_region_permutation.py --work GSE208248_LINK_WORK
python scripts/85_local_candidate_decision.py
```

Outputs are in `results/local_signal_decision/`. This report supplements the
complete 565/645 individual-peak screen; those previously reported peak
statistics remain unchanged.

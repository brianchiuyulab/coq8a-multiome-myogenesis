# Day-0 ADAM12-neighborhood candidate and local RNA assessment

## Analysis sequence

1. The external Hallmark/Reactome myogenesis union contains 221 genes. Its intersection with GO myoblast differentiation or myoblast fusion contains 25 genes. ADAM12 belongs to the fusion set in the frozen `results/differentiation_fusion25/genes.tsv`; it was not added after the current peak screen.
2. The existing 565 common candidate peaks in those neighborhoods form the testing scope. This is a pathway-directed screen, not a whole-genome screen.
3. Use only undifferentiated GSE208248 libraries GSM6339597 and GSM6339601. Apply joint RNA/ATAC QC and compare COQ8A groups within each library. The complete sensitivity grid and all outcomes are archived in `DAY0_SENSITIVITY_GRID.md`.
4. D232 uses TSS enrichment >=3, at least one detected UMI in PAX7/MYF5/MYOD1/MYOG, COQ8A >=3 versus exactly 1 UMI, and within-library RNA/ATAC log-depth matching with caliper 0.05. There are 35 and 45 pairs, respectively. This marker-defined population is not an independently validated pure MuSC annotation.
5. All 565 peaks and all 25 regions receive the same tests. ADAM12 emerges as an opening candidate. Its focal peak is chr10:126373783–126374503: high 30/80 versus low 7/80 open nuclei, FC 4.2857, peak p approximately 0.0000340, BH q565=0.01918. The regional maximum-statistic permutation test across its 48 peaks gives opening p=0.0000965 and BH q25=0.0024125. Both biological sources increase. These are nucleus-pair tests from two sources, not 80 biological replicates.
6. Having selected the ATAC locus, assess all deposited RNA features with an annotated transcript TSS within 500 kb. Do not require membership in the original 25-gene list or assume the neighborhood label identifies its target.
7. Cross-reference measurable mouse homologs against the private C2C12 workbook. Private expression values are not published in this repository.

The external gene scope is independent of the current COQ8A effects. D232 itself was selected by sensitivity exploration, not prospectively or blindly selected. Nearby depth-only calipers also support the regional signal; additional state matching attenuates it. The within-setting q values do not adjust for selecting a setting from the grid.

## Local RNA follow-up

Script `90_adam12_day0_targets.py` re-extracts the six local RNA features and focal ATAC peak from deposited H5 files. It confirms exactly 80 pairs and the 30/7 opening counts. Transcript-TSS distances from the peak midpoint are:

| RNA | Type | Distance (bp) |
|---|---|---:|
| ADAM12 | protein coding | 14,311 |
| LINC00601 | lncRNA | 45,098 |
| C10orf90 | protein coding | 89,232 |
| FANK1 | protein coding | 369,996 |
| FANK1-AS1 | lncRNA | 401,018 |
| DHX32 | protein coding | 477,708 |

Links are partial correlations between binary peak accessibility and log1p RNA CP10k, estimated separately by source, controlling log RNA depth, log open-peak depth and log1p COQ8A CP10k. A second model adds the grid's myogenesis-state score. Fisher-z combinations use residual degrees of freedom; each source needs at least ten peak-positive and ten RNA-positive nuclei. These asymptotic nucleus-level tests are exploratory target assessments. Both the full eligible population (13,434 nuclei) and the exact matched population (160 nuclei) are retained. Local BH includes all six candidate hypotheses, assigning missing tests p=1 for adjustment and displaying their p/q as unavailable.

No local RNA has a convincing positive link in these assessments. ADAM12 in the eligible population has r=-0.00614, p=0.4766; with state adjustment, r=0.00578, p=0.5032. C10orf90 has r=0.01294, p=0.1338; with state adjustment, r=0.01528, p=0.07676. These correlations are small. FANK1-AS1 is undetected; LINC00601 has only 55 detected nuclei across the full eligible population.

Separately, the same 80 pairs give ADAM12 RNA FC=0.8326, paired log1p-CP10k t-test p=0.3253. C10orf90 RNA FC=4.9076, p=0.05220, but only nine matched nuclei detect it. This does not establish C10orf90 as the target. RNA high/low comparison and peak–RNA target identification are distinct analyses; neither is a mandatory gate for retaining the accessibility candidate.

**Current conclusion:** a COQ8A-associated opening signal at an ADAM12-neighborhood locus is supported in this exploratory setting. The regulated RNA remains unresolved. It should not yet be called an ADAM12 transcriptional activation mechanism.

## Reproduction

```powershell
python scripts/90_adam12_day0_targets.py --work <day0_sensitivity_grid_cache> --h5 <GSE208248_processed>
```

Outputs are in `results/adam12_day0_targets/`: candidate identities, source-specific and combined links, and paired RNA comparisons. All local candidates, including undetected features, remain in the output.

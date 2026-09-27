# Task 1: undifferentiated-culture analysis

This report covers the restricted two-contrast run. The subsequently requested
[expanded sensitivity grid](DAY0_SENSITIVITY_GRID.md) identifies additional
opening and closing candidates under other definitions. The negative opening
FDR result below must not be generalized to that expanded grid.

## Main result

The stage-specific analysis is complete. In undifferentiated cultures, CSRP3
has the smallest regional opening p value among the 25 functional neighborhoods
under the fixed primary setting (TSS >=3, COQ8A >=2 vs exactly 1 UMI). Its
regional p=0.03684 and BH q25=0.51367. Thus the stage-specific opening screen
does not retain the earlier pooled regional q<0.10 result. CSRP3 remains a
nominal candidate, not an FDR-supported Day-0 opening discovery.

The primary Day-0 focal peak is chr11:19201752-19202603, proximal to a CSRP3
transcript TSS. Its open-fraction ratio is 1.85, with positive direction in
both sources. This is a different peak from the previously highlighted
chr11:19218592-19219518, whose Day-0 ratio is 1.4754. Both are retained.

## Design fixed for this run

- Dataset: GSE208248, paired RNA/ATAC in the same nuclei.
- Primary libraries: GSM6339597 and GSM6339601, undifferentiated cultures.
- Separate Day-7 libraries: GSM6339599 and GSM6339603.
- Primary contrast: raw COQ8A UMI >=2 versus exactly 1; zeros excluded.
- Retained contrast sensitivity: >=3 versus exactly 1.
- Joint QC: RNA UMI >=500, detected ATAC peaks >=500, within-library 95th
  percentile depth caps, FRiP >=0.25, nucleosome signal <4, blacklist fragment
  fraction <0.05, custom archived TSS enrichment >=3.
- Matching: within-library 1:1 without replacement, both log-depth calipers
  <=0.10, the same implementation as the preceding analysis.
- Scope: reconstruct the 221-gene Hallmark/Reactome union intersected with the
  GO differentiation/fusion union, yielding 25 genes. Reconstruct the 565
  distinct peaks from the frozen 5,097 four-library-common candidate map.
- Coordinate scope remains the original canonical-chromosome, 200-2,000-bp,
  nonblacklist, reciprocal-overlap >=50%, midpoint-to-transcript-TSS <=100-kb
  definition. The peak catalog is kept common to four libraries for direct
  comparison; it is not a new Day-0-only peak-calling catalog.
- No early/middle/late, motif, target identity or observed RNA difference is
  an eligibility requirement.

The saved fragment QC scores and processed author peak map were reused; QC
thresholds and depth matching were reapplied. Rebuilt barcode pairs exactly
match the archived within-library pairs. All downstream stage-specific peak,
regional and RNA comparisons were newly computed from the extracted matrices.
No data were downloaded again, and alignment/peak calling were not rerun.

Undifferentiated is the culture condition. This run does not assign a new cell
type to every nucleus or isolate a marker-defined pure MuSC subset. Depth,
myogenesis-score and mitochondrial-fraction summaries for high/low groups are
reported in depth_and_state_balance.tsv; they are descriptive, not new gates.
The biological source count is two. Tests below use nuclei, not independent
donor replication.

| Stage | Contrast | Source 1 pairs | Source 2 pairs | Total pairs |
|---|---|---:|---:|---:|
| Undifferentiated | >=2 vs 1 | 211 | 296 | 507 |
| Undifferentiated | >=3 vs 1 | 46 | 58 | 104 |
| Day 7 | >=2 vs 1 | 326 | 125 | 451 |
| Day 7 | >=3 vs 1 | 74 | 23 | 97 |

Nuclei are never paired across stages. These contrasts do not track an
individual Day-0 nucleus into Day 7.

## Regional and individual-peak statistics

Use the previously specified maximum standardized paired-difference statistic
within each gene neighborhood. Run 2,000,000 within-pair permutations for every
stage/contrast with seed 20260928; apply BH across all 25 regions separately
for opening and two-sided hypotheses. Individual two-sided discordant-pair
binomial tests retain BH over all 565 peaks in each stage/contrast. Permutation
statistics, counts, raw p and q are retained for every candidate, not only CSRP3.

| Stage / contrast | First opening region | Regional p | Regional q25 | CSRP3 p | CSRP3 q25 |
|---|---|---:|---:|---:|---:|
| Undifferentiated, >=2 vs 1 | CSRP3 | 0.03684 | 0.51367 | 0.03684 | 0.51367 |
| Undifferentiated, >=3 vs 1 | ADAM12 | 0.005005 | 0.12512 | 0.02054 | 0.25674 |
| Day 7, >=2 vs 1 | MAPK12 | 0.14035 | 0.74308 | 0.15927 | 0.74308 |
| Day 7, >=3 vs 1 | MYOG | 0.04789 | 0.72119 | 0.48489 | 0.95337 |

No opening region reaches q<0.10 in these four comparisons. Sensitivity results
remain separate; the >=3 contrast is not substituted for the primary contrast.

### Day-0 CSRP3 peaks under the primary setting

| Peak | High / 507 | Low / 507 | Open-fraction FC | Individual p | Individual q565 |
|---|---:|---:|---:|---:|---:|
| chr11:19201752-19202603 | 37 | 20 | 1.8500 | 0.01151 | 1.000 |
| chr11:19218592-19219518 | 90 | 61 | 1.4754 | 0.01412 | 1.000 |

The proximal peak increases from 3.94% to 7.30% (3.35 percentage points).
Source 1 counts are 17 high versus 11 low out of 211; source 2 counts are
20 versus nine out of 296. A region's q is not an individual-peak q.

### Opposite-direction finding retained

The primary Day-0 peak chr22:36352375-36353254, in the MYH9 search neighborhood,
is less accessible: 23/507 versus 58/507, FC=0.39655, p=0.0001026,
q565=0.05797. Both sources decrease (10 vs 27; 13 vs 31). The two-sided MYH9
regional test gives p=0.003068, q25=0.07670. This is an accessibility-closing
candidate and does not establish MYH9 as its RNA target or support the requested
opening mechanism. It is preserved rather than filtered out.

## RNA response and target assignment are separate

CSRP3 RNA in the same Day-0 507 pairs has mean-CP10k FC=2.67624,
paired-log1p p=0.007305 and q25=0.09131. Both sources have higher mean expression
(FC 3.72 and 1.71). At Day 7, combined RNA FC=0.8381, p=0.6455, q25=0.7335;
the second Day-7 source has zero detected CSRP3 RNA in both matched groups.
These RNA contrasts are supplementary and are not ATAC candidate-selection gates.
They do not demonstrate a delayed rise in RNA after differentiation or an
increased response speed.

All eligible local RNA links near all 565 peaks were recombined by stage from
the archived within-library fits. Both sources must meet the detection rule;
5% and 1% peak/RNA detection are reported separately. This uses the broader
QC-passing nuclei, not only the 507 pairs. BH is within each stage, detection
criterion and model over all evaluable links; missing links remain unavailable.

For the newly leading proximal peak, the CSRP3 association at 1% detection is
positive in both Day-0 sources. Technical/COQ8A-adjusted r=0.05863,
p=2.47e-14, q=1.70e-12 across 5,656 links. Adding myogenesis state reduces it to
r=0.01911, p=0.01304, q=0.16769. The 5% analysis has no evaluable link for this
peak. The positive association and proximal position make CSRP3 a plausible
target candidate but do not establish direct regulation. All alternative RNA
targets remain in stage_specific_local_RNA_links.tsv.gz.

## Figures and reading guide

### Figure 1: complete Day-0 screen

Panel A ranks all 25 neighborhoods by regional opening p. Bar length is
-log10(p); the right column contains BH q over 25 regional tests. Panel B shows
all 565 finite positive peak ratios against individual two-sided -log10(p).
The orange points identify the two discussed CSRP3 peaks; the blue point is the
MYH9-neighborhood closing candidate. Color alone does not encode significance.
Panel A and B use different test units and are not interchangeable p values.

### Figure 2: proximal CSRP3 peak by source and stage

Each line joins the low- and high-group open fractions within one library.
Panel A contains only the two undifferentiated libraries, and panel B only the
two Day-7 libraries. These lines do not represent single-cell trajectories or
measure the same nucleus twice. Both panels use the same percentage scale.

Figures are in figures/task1_stage_specific/ as PNG and vector PDF. Source
tables are in results/task1_stage_specific/. Earlier pooled and temporal
figures remain separate historical results.

## Reproduction and verification

```bash
python scripts/86_stage_specific_task1.py --work GSE208248_LINK_WORK
python scripts/87_report_stage_specific_task1.py
```

The work directory is the existing script-48 extraction, containing peaks.txt,
genes.txt and per-library metadata/ATAC/RNA matrices. Script 86 also consumes
the committed processed nucleus inventory, fragment-QC references and common
peak map. Script 87 independently reconciles all stage/library peak counts
against the previous library tables, recomputes all single-peak p/q values,
checks regional permutation p values and BH families, and renders figures.
validation.json records successful checks; provenance.json records input hashes.

The current result is a nominal Day-0 CSRP3 accessibility candidate with
concordant RNA, not a successful q<0.10 opening-region screen. The earlier
pooled q=0.0956 is unchanged but answers a different, all-stage question.

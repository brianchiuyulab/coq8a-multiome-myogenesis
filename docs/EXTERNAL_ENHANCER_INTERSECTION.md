# Independent enhancer annotation of the MYOD1 candidate

## Temporal membership and candidate discovery answer different questions

The highlighted GRCh38 interval, **chr11:17649919–17650798**, was recovered from
the full 5,097-peak COQ8A screen and linked to MYOD1 RNA. This is a candidate
association, not a newly calculated gene-set or motif enrichment result.

Three other MYOD1-neighborhood peaks are in the external early-opening set:
chr11:17696914–17697783, chr11:17714073–17714965 and
chr11:17806817–17807701. Thus MYOD1-associated regions are not altogether absent
from the temporal analysis. Each genomic interval has its own classification.

For the highlighted interval, depth-standardized detection in GSE109828 experiment
1 is 30.41%, 25.30%, 24.28% and 25.88% at 0, 24, 48 and 72 hours. Experiment 2
has 21.97% at 0 hours and 20.04% at 72 hours. The external dynamic test has
p=0.0542 and q=0.1823. Neither experiment shows a positive endpoint change.
Consequently this region does not meet the opening criteria used before assigning
early/middle/late onset. It is already accessible at the first measured time;
these data do not establish statistical equivalence over time.

The COQ8A comparison, in contrast, asks about high/low nuclei within the target
GSE208248 libraries, not progression along GSE109828 differentiation time.

## External sources and fixed intersections

Human HSMM chromatin-state segmentation was downloaded from the ENCODE/Broad
ChromHMM track `wgEncodeBroadHmmHsmmHMM` in hg19. States **4_Strong_Enhancer**
and **5_Strong_Enhancer** define the enhancer set. These are external chromatin-state
annotations, not individual perturbation validation of every element. The hg19
track is the published UCSC-distributed coordinate version; no new state model
was fitted to the target data.

Human MYOD1 ChIP peaks were taken from GSE50413 myoblast and 72-hour myotube
samples, GSM1218849 and GSM1218850. A candidate peak is MYOD-bound if it overlaps
at least one author peak call from either condition. Both source files and the
ChromHMM file are identified by SHA256 in the reference manifest.

Of 5,097 starting intervals, 5,095 have valid unique hg19 mappings. The two
unmapped intervals are not used for external annotation. Any positive overlap
with an external interval defines membership; covered fractions are retained.
The rules do not use COQ8A effects, target p values, or peak–RNA correlations.

| Intersection | Number of target peaks |
|---|---:|
| Existing four-library myogenesis neighborhoods | 5,097 |
| Overlap HSMM strong-enhancer states | 1,777 |
| Additionally overlap independent MYOD1 ChIP | 942 |
| Additionally belong to original candidate neighborhoods of MYOD1, MYOG, MYF5 or MYF6 | 23 |

The four MRF genes are the externally established core regulator family discussed
in the human epigenome study below. The final 23 comprise ten MYOD1-neighborhood
peaks, eleven MYOG-neighborhood peaks and two peaks shared by the MYF5/MYF6
neighborhoods. The shared pair is one regulatory set, not two independent tests.
Neighborhood assignment alone is not proof of a target gene.

## Focal coordinate verification

The highlighted target maps uniquely to **hg19 chr11:17671466–17672345**. All
879 bp are inside the independently called HSMM segment
**chr11:17671424–17673224, state 4_Strong_Enhancer**. The same target interval
overlaps MYOD1 ChIP calls in both myoblast and myotube samples.

Human MYOD1 neighborhood work describes an approximately 6-kb enhancer-chromatin
region centered roughly 67 kb upstream of MYOD1, with myogenic DNase sensitivity.
The position of the present interval is consistent with that reported far-distal
enhancer-chromatin neighborhood. Exact overlap with the downloaded strong-enhancer
segment is verified independently of the paper's rounded distance designation.
This should not be relabelled as a previously proven reporter element solely
from that positional agreement.

## Exploratory tests after external annotation

Existing TSS ≥3 QC, COQ8A ≥3-versus-1 UMI definitions and 201 matched pairs remain
fixed. Within each nucleus, a set score is the mean binary detection across its
member peaks. Paired score differences are tested with a two-sided t test across
the 201 pairs. All six specified, distinct sets are reported; module q values
adjust across those six tests. Per-peak q values within each external set are
provided separately and are not mixed with module q values. These newly examined
sets remain exploratory rather than retrospectively preregistered hypotheses.

| Set | Peaks | Module high/low FC | Raw p |
|---|---:|---:|---:|
| HSMM strong enhancers | 1,777 | 1.004 | 0.604 |
| HSMM strong, MYOD-bound enhancers | 942 | 1.017 | 0.145 |
| Core MRF neighborhood union | 23 | 1.090 | 0.123 |
| MYOD1 neighborhood | 10 | 1.126 | 0.155 |
| MYOG neighborhood | 11 | 1.029 | 0.708 |
| Shared MYF5/MYF6 neighborhood | 2 | 1.560 | 0.0566 |

The highlighted MYOD1 interval ranks first by nominal ATAC p among the 23 core
MRF peaks. Its original effect is still FC=1.630, p=0.02129; BH q is 0.4898
across the 23 core MRF peaks or 0.2129 within the ten MYOD1-neighborhood peaks.
The external intersections therefore strengthen its regulatory annotation and
make the candidate prioritization explicit, but do not produce an FDR-significant
individual COQ8A accessibility effect. No set gives a significant module result.

The two MYF5/MYF6 intervals have a combined positive direction in four libraries,
but one has only nine versus two detected matched nuclei (FC=4.5, p=0.0654).
Its large ratio does not by itself establish a stronger candidate than MYOD1.

## Candidate framing and other references

The current priority is a **pre-existing myogenic enhancer-chromatin interval
associated with MYOD1 RNA and COQ8A status**, rather than a claim that COQ8A
initiates a time-specific opening event. Independent chromatin-state annotation,
MYOD1 occupancy, within-library RNA coupling and the target accessibility
contrast represent distinct evidence components.

GSE37525 supplies C2C12 myoblast/myotube p300 and H3K27ac peak sets and is a
relevant second reference for enhancer activity. Mouse peaks require assembly
verification and cross-species mapping before intersection with the human target
intervals; no mouse intersection is claimed here. A July 2026 preprint describes
reporter-active mouse MyoD enhancers near −36 and −60 kb. It has not been used
as a validated human-coordinate reference or as a filter in this analysis.

## Figure legend

**Independent annotation and two distinct accessibility comparisons at the
MYOD1 distal candidate.** A, hg19 genomic intervals from HSMM strong-enhancer
annotation and independent MYOD1 ChIP peak calls, with the mapped target interval
shaded. Rectangles indicate interval calls, not quantitative ChIP signal heights.
B, the geometrically matched external GSE109828 accessibility profile; experiment
2 has endpoint measurements only. Lines connect measured time points and are
not fitted trajectories. C, peak-detection fractions in 201 matched COQ8A-high
and low nuclei in GSE208248, with the original paired nominal p and fold change.
The figure does not imply that the external time-course experiment has measured
COQ8A or that the target experiment follows the same nuclei through time.

## Reproduction and sources

```bash
python scripts/52_external_enhancer_sets.py --hmm /data/wgEncodeBroadHmmHsmmHMM.bed.gz --chip /work/middle_links/external --cache /work/unstratified --out /work/external_enhancers
python scripts/53_report_external_enhancers.py --work /work/external_enhancers --hmm /data/wgEncodeBroadHmmHsmmHMM.bed.gz --chip /work/middle_links/external
```

The reference-manifest URL provides the exact ChromHMM download. The MYOD1
source downloads are described in the middle-region follow-up. Render the figure
from committed tables with `python scripts/53_report_external_enhancers.py`.

- [Human MRF neighborhoods and far-distal enhancer chromatin, 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4512632/)
- [ENCODE/Broad ChromHMM track documentation](https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&g=wgEncodeBroadHmm)
- [Human MYOD1 ChIP, GSE50413](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50413)
- [C2C12 enhancer study, 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3533080/)
- [C2C12 enhancer source data, GSE37525](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37525)
- [Mouse MyoD enhancer preprint, July 2026](https://pubmed.ncbi.nlm.nih.gov/42619749/)

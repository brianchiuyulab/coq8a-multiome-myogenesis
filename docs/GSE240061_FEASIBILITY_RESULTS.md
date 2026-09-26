# GSE240061: measured feasibility results

## Input and pairing

The deposited `GSE240061_integrated11723.rds.gz` was downloaded (3,736,503,149
bytes; SHA-256 `556f19882a6797ed46344d54c8c1ec948a7491270522bd43f1d45614b6239fad`).
Its outer gzip layer was decompressed with CRC verification, yielding the
gzip-compressed author RDS. The object contains raw sparse RNA counts for
36,601 genes and ATAC counts for 144,663 peaks, across 37,154 nuclei.
RNA, ATAC and metadata barcode sets were verified identical and aligned by name.
Genome assembly hg38 is specified in GSM7680744 GEO processing metadata.

The author's `refined_annotations_wknn_0.8` field supplies fourteen cell types.
It contains 1,502 Satellite Cells (984 Pre and 518 Post). No new clustering or
cell-type relabeling was performed. Sample labels were reconciled against the
GEO manifest: E/G/I/J are exercise participants and L/N are controls, each with
M1/Pre and M2/Post. GEO's PBMC tissue labels conflict with the paper and actual
muscle cell-type annotations; they remain verbatim in the source manifest.

## COQ8A group availability before matching

| Baseline cell type | All nuclei | COQ8A >=3 | COQ8A >=2 | COQ8A =1 |
|---|---:|---:|---:|---:|
| Satellite Cells | 984 | 16 | 38 | 131 |
| Fast myonuclei | 8,735 | 1,045 | 2,240 | 2,549 |
| Slow myonuclei | 7,332 | 1,446 | 2,744 | 2,044 |

Sixteen baseline satellite high nuclei occur in four participants. The baseline
satellite counts by participant (all/high>=3/low=1) are E 255/4/29, G 49/0/5,
I 148/3/14, J 363/7/54, L 12/0/2 and N 157/2/27. Absence of a high group, rather
than an arbitrary per-donor nucleus cutoff, prevents some paired comparisons.

## Transfer of the existing matching rule

The original matching function was reused without changing its algorithm:
log1p total RNA UMI and log1p detected ATAC peaks, maximum absolute mismatch
0.10 in each coordinate, without replacement. Matching is now within both
sample and author cell type. It uses author-retained nuclei; full old fragment
QC and the old library-specific upper-depth caps have not been reconstructed.

| Baseline cell type | >=3 versus 1 matched pairs | Donors with pairs | >=2 versus 1 matched pairs | Donors with pairs |
|---|---:|---:|---:|---:|
| Satellite Cells | 3 | 1 | 12 | 3 |
| Fast myonuclei | 846 | 6 | 1,707 | 6 |
| Slow myonuclei | 1,006 | 5 | 1,630 | 5 |

Satellite matching sensitivity was evaluated without inspecting ATAC outcomes:

| Log-depth caliper | >=3 versus 1 pairs / donors | >=2 versus 1 pairs / donors |
|---|---|---|
| 0.10 | 3 / 1 | 12 / 3 |
| 0.20 | 4 / 1 | 17 / 4 |
| 0.30 | 6 / 3 | 21 / 4 |

These are feasibility counts, not differential-accessibility significance tests.
Increasing the caliper permits more depth imbalance and does not create new
biological replicates. Strictly copying the current high/low matching design
therefore gives poor satellite-cell replication capacity despite six recruited
participants. Myonuclei have substantially more usable pairs, but test a mature
muscle context, not myoblast differentiation. Baseline and post samples are
not independent cohorts; post results must preserve within-participant dependence.

## Regional correspondence

| Original scope | Old peaks with any overlap | Old peaks with >=50% reciprocal overlap | Distinct author peaks at >=50% reciprocal overlap |
|---|---:|---:|---:|
| 5,097 candidates | 3,212 | 3,071 | 3,071 |
| 1,777 HSMM enhancers | 1,087 | 1,026 | 1,026 |
| 401 muscle-selective enhancers | 209 | 191 | 191 |
| 31 fusion-promoter candidates | 24 | 23 | 23 |

CAV3's original interval chr3:8733438-8733968 is completely contained in the
author peak chr3:8733222-8734203. The original width is 530 bp and the author
width is 981 bp; reciprocal overlap is 1.00 and 0.5403. Thus this region can be
assessed with the author peak matrix, while exact-coordinate replication requires
fragment recounting. Original intervals without overlap remain unmeasured rather
than being assigned zero accessibility.

## Remaining processing requirements

The object includes RNA/ATAC count and feature totals, mitochondrial percentage
and author labels, but lacks per-nucleus TSS enrichment, FRiP, blacklist ratio
and nucleosome signal; the ATAC position-enrichment slot is empty. Consequently
no TSS>=3 replication or TSS sensitivity has yet been performed. Those metrics
require deposited fragments or a separate author QC table. Author retained-cell
analysis is possible as a clearly defined initial scope; it is not numerically
identical preprocessing to the old dataset.

No COQ8A differential-accessibility p or q has been computed for GSE240061 in
this feasibility run. The useful next analyses are (i) source-aware associations
in separately analyzed myonuclear subtypes, and (ii) a distinct satellite
analysis assessing whether a continuous/weighted model has sufficient expression
support without pretending that three matched pairs represent six donors.
Neither is labeled a confirmed CAV3 effect before being tested.

## Reproduction

After the author-object inspection in script 61:

```bash
python scripts/63_replication_feasibility.py --work /work/GSE240061_inventory --out results/replication_GSE240061
```

The public directory contains cell-type/sample counts, matching summaries,
caliper sensitivity, coordinate coverage and the GEO sample manifest. Large
downloaded inputs, raw count subsets and nucleus-level matching files remain
outside the release repository.

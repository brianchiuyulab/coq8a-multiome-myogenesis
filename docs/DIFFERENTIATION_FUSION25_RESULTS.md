# Differentiation/fusion union: completed analysis

Analysis date: 2026-09-27. Public analysis source code: scripts 76–82. All numerical source
tables are in `results/differentiation_fusion25/`.

## Findings

The 25-gene functional union defines **565 distinct GSE208248 peaks** and
**645 native GSE240061 peaks**. The new dataset includes 335 peaks mapped to
old candidates plus 310 additional peaks, independently selected from its
complete author peak set using the same gene-neighborhood and coordinate QC
rules. Final Task 2 results are in `native_GSE240061/`. The 335-peak mapped-only
tables in the parent directory are retained as a scope comparison.

In GSE208248, the complete union does not show increased average accessibility:
COQ8A-high 10.483%, low 10.519%, FC **0.99657**, matched-nucleus paired t
**p=0.7831**. Individual opening candidates remain, but no individual peak
passes union-wide q<0.10 under any of the four retained QC/COQ8A settings.

In GSE240061, none of the primary cell-type/time comparisons yields an opening
peak at q<0.10. One closing interval in Slow Post passes q<0.10. Sensitivity
analyses identify additional opening candidates, including a sparse MYOD1
neighborhood interval; their complete sample support and RNA results are below.

## 1. Scope and computation

The frozen MSigDB GO myoblast differentiation and myoblast fusion memberships
were intersected with the existing Hallmark/Reactome 221-gene universe. The
17 differentiation and 12 fusion genes share four members, yielding 25 distinct
genes. Membership is recorded in `genes.tsv`. No temporal class, motif,
promoter-only window, muscle-selective enhancer annotation, or observed effect
was required for this union analysis.

The existing GSE208248 candidate map supplies peaks within the original 100-kb
TSS-neighborhood rule and present in all four libraries. Each interval is
tested once even if it belongs to multiple gene neighborhoods. These gene
labels define inclusion, not experimentally established target genes.

The original extraction, QC, pairing, normalization, and per-peak model inputs
are retained. GSE208248 paired binomial p values were recomputed from high-only
and low-only accessible-pair counts. All four count fields for all 565 primary
peaks were independently checked against the extracted nucleus matrices.
The initial GSE240061 mapped-only comparison retains the original model fits
and changes only BH families. The final native-peak analysis refits edgeR over
the original candidate background plus the additional native union intervals,
retaining genome-wide TMM factors, raw library sizes, nucleus membership and
design matrices. Dispersion is estimated over this expanded background.

BH correction is performed over all 565 GSE208248 peaks per setting, and over
all count-supported native union peaks per GSE240061 configuration. RNA BH
families contain the evaluable union genes. Sensitivity settings are retained
separately, not pooled into a new BH family. Searching those settings is
exploratory. The union was agreed after prior analyses; it is not a retrospective
claim of preregistration or an independent validation cohort.

## 2. Task 1: GSE208248

### Design and primary setting

Two donor-derived cell lines, each undifferentiated and day-7 differentiated,
provide four libraries (GSM6339597, GSM6339599, GSM6339601, GSM6339603).
Same-nucleus RNA/ATAC measurements are compared within library.

The retained primary setting is TSS enrichment >=3, COQ8A >=3 versus exactly
1 RNA UMI, with 201 within-library depth-matched pairs. Zero-UMI nuclei are
not the low group. Both RNA and ATAC pairing/QC are unchanged. Pair-level
tests quantify nucleus associations; the biological source count remains two.
Per-library effects accompany pooled results.

### Whole union and all-peak screen

Mean per-nucleus accessibility is the fraction of the 565 peaks with nonzero
ATAC signal. Its FC is the high/low mean ratio, not an average of individual
peak FCs. Library-specific FCs are 0.9904, 0.9693, 1.0448 and 0.9930, respectively.

| Setting | N matched pairs | Raw p<0.05 opening peaks | Raw p<0.05 closing peaks | Any peak q<0.10 |
|---|---:|---:|---:|---:|
| TSS>=3, COQ8A>=3 vs 1 | 201 | 12 | 7 | 0 |
| TSS>=2, COQ8A>=3 vs 1 | 212 | 13 | 8 | 0 |
| TSS>=3, COQ8A>=2 vs 1 | 958 | 15 | 16 | 0 |
| TSS>=2, COQ8A>=2 vs 1 | 1020 | 18 | 13 | 0 |

The opening candidates below are selected for discussion after the complete
565-peak screen. Their q values refer to that screen, not a later smaller subset.

| Neighborhood gene / interval (hg38) | ATAC FC | Raw p | q, 565 peaks |
|---|---:|---:|---:|
| CSRP3; chr11:19201752–19202603 | 2.4286 | 0.03088 | 1.000 |
| CSRP3; chr11:19218592–19219518 | 1.9412 | 0.01659 | 1.000 |
| MYOD1; chr11:17649919–17650798 | 1.6296 | 0.02129 | 1.000 |
| CAV3; chr3:8733438–8733968 | 1.7647 | 0.05325 | 1.000 |
| MYF5; chr12:80679105–80680024 | 4.5000 | 0.06543 | 1.000 |

MYF5's 4.5 ratio comes from 9 versus 2 accessible nuclei, not high-count evidence.
The smallest primary p is a closing CACNA1H-neighborhood peak at
chr16:1152563–1153401: FC=0.2273, p=0.000911, q=0.5144.

### RNA integration

RNA high/low tests use the same 201 matched pairs. The reported RNA FC below is
the ratio of mean CP10k expression. The paired p is calculated on log1p CP10k;
it is not a donor-level p or a count-model FC. RNA q is across the 25 genes.

| Gene | RNA mean-CP10k FC | Raw p | RNA q |
|---|---:|---:|---:|
| CSRP3 | 5.2248 | 0.01210 | 0.1512 |
| MYOD1 | 1.3021 | 0.05104 | 0.1823 |
| CAV3 | 1.2119 | 0.04316 | 0.1823 |
| MYF5 | 1.0726 | 0.5550 | 0.6514 |
| MYOG | 1.5035 | 0.04567 | 0.1823 |
| MEF2C | 1.0046 | 0.6510 | 0.7076 |

CSRP3 is sparsely detected; its RNA ratio should be read with the absolute
expression and detection-supported links. At the primary >5% RNA-detection
criterion, its focal peak–RNA link is not evaluable. At the retained 1% sensitivity,
the chr11:19201752–19202603 link has partial r=0.050. Primary MYOD1 and CAV3
focal links have partial r=0.135 and 0.074, respectively. These are small positive
associations in the QC-passing culture nuclei, adjusted for depth and COQ8A;
they are not correlation estimates restricted to the 201 high/low pairs.
Complete local-target results, including alternative RNAs and the two detection
thresholds, are retained. No union RNA passes q<0.10.

## 3. Task 2: GSE240061

### Design and primary setting

Six young adults (exercise E/G/I/J; rest L/N), each with Pre and Post samples,
provide paired RNA/ATAC. This adult-muscle study is not a cultured-myoblast
differentiation time course. Original cell annotations are retained. Satellite
Cells, Fast, Slow and Intermediate nuclei, and Pre/Post, are analyzed separately.

The common primary comparison uses author QC, COQ8A >=3 versus 1 UMI and depth
matching caliper 0.10. The author object lacks per-nucleus TSS/FRiP; these
metrics have not been reconstructed or claimed equivalent to GSE208248 TSS>=3.
Donor pseudobulk models use genome-wide ATAC TMM normalization and ~donor+group.
The outcome-blind count filter is total counts>=10 and presence in >=2 donors.

| Context | Donors | Matched pairs | Tested union peaks | Opening q<0.10 | Closing q<0.10 |
|---|---:|---:|---:|---:|---:|
| Satellite Cells Pre | 1 | 3 | Not estimable | — | — |
| Satellite Cells Post | 1 | 1 | Not estimable | — | — |
| Fast Pre | 6 | 846 | 548 | 0 | 0 |
| Fast Post | 6 | 149 | 387 | 0 | 0 |
| Slow Pre | 5 | 1006 | 568 | 0 | 0 |
| Slow Post | 6 | 307 | 487 | 0 | 1 |
| Intermediate Pre | 5 | 94 | 221 | 0 | 0 |
| Intermediate Post | 4 | 8 | 10 | 0 | 0 |

The complete 645-peak module also shows no primary increase. Equal-donor mean
accessibility FCs for Fast Pre/Post, Slow Pre/Post and Intermediate Pre/Post
are 0.9532/1.0069, 0.9889/1.0028 and 0.8921/0.9160. Their donor-paired raw
p values are 0.346/0.651, 0.510/0.756 and 0.273/0.101. BH across the six
evaluable primary contexts gives q>=0.605. Module scores use all 645 peaks,
not only count-supported individual-peak tests.

### Primary closing interval

Slow Post **chr16:1308794–1309860**: FC=**0.6564**, p=**0.0001789**,
q=**0.08713** over 487 tested native union peaks. Normalized counts decrease in all
six donors. This is a different interval from the previously discussed
chr16:1311478–1312392 opening peak.

It enters through the CACNA1H neighborhood but lies within the UBE2I genomic
span. Nearby transcript annotation shows that it includes an annotated lncRNA
ENSG00000274751 TSS and several UBE2I transcript TSSs. This
does not establish which transcript is regulated. The annotation comparison
uses unrestricted GENCODE v48 annotation of every native candidate interval.
At matching calipers 0.20/0.30, FC=0.714/0.741 and raw p=0.00132/0.00357,
but union q=0.671/0.889. Thus the closing direction persists, while adjusted
significance is not invariant to matching.

### Opening sensitivity findings

All 448 existing configurations were refitted for the native union; 392 yield
count-model results. The remaining configurations have fewer than two donors
(38), insufficient count-supported background features (15), a non-estimable
design (1), or numerical fitting failures (2). Every status is retained; a
failed fit is not recorded as a null effect.

1. **Fast Post chr16:1311478–1312392**, trim95, COQ8A>=2 vs 1,
   caliper 0.20 (C198): six donors, 490 pairs, FC=3.0067, p=1.006e-5,
   **native-union q=0.005222 over 519 peaks**. The preceding mapped-only union
   gave q=0.002245 over 274 peaks; the original broad candidate scan gave
   q=0.0224. The native fit changes the dispersion background and test family,
   with essentially unchanged effect size.
   The preceding local-RNA follow-up still applies: UBE2I RNA FC=1.005,
   p=0.983 and peak–RNA r=0.006. The narrower q does not supply an RNA target.

2. **Satellite Cells Post chr11:17719024–17719942**, MYOD1 neighborhood,
   author QC, COQ8A>=2 vs 1, caliper 0.20 (C058): four donors, **11 pairs**,
   FC=10.4266, p=0.001506, **q=0.03765 over 25 count-supported native union peaks**.
   These are 25 peaks passing the count filter, not the 25 input genes.
   Unrestricted annotation confirms overlap with a protein-coding MYOD1
   transcript TSS (87 bp from the peak midpoint).
   Raw high/low counts total 15/1. Per-donor high/low counts are E 6/0,
   G 7/1, J 2/0, L 0/0, from 5, 4, 1 and 1 pairs, respectively.
   At caliper 0.30 there are six donors and 16 pairs: FC=4.428,
   p=0.01598, q=0.5034. MYOD1 RNA has **zero counts in all 22 C058 nuclei**;
   the broader Satellite Post peak–MYOD1 RNA partial r=-0.047, p=0.175.
   This is an ATAC exploration lead, not a matched RNA-supported MYOD1 mechanism.
   Its corresponding old-data interval chr11:17719034–17719916 has FC=0.9362,
   p=0.7877, q=1 in the old primary analysis. It is not the old nominal MYOD1
   interval chr11:17649919–17650798 and does not replicate that interval's effect.

No newly added native-only interval yields a matched opening q<0.10. The
mapped-only C327 signal at chr16:1105996–1106999 does not pass native-union
q<0.10. All native results, including unmatched models, remain available.

The CDON-neighborhood C057 result has q=0.0654 because only one union peak
passes the count filter, but raw p=0.0654 and only four matched pairs from
three donors. It is retained in the full sensitivity table, not called p<0.05.
Unmatched model results are also retained; the main conclusions use matched
comparisons. No claim is made that the minimum p across settings is a calibrated
single confirmatory test.

## 4. Figure and interpretation

`figures/differentiation_fusion25/union_peak_screen.png` and `.pdf` show the
complete old primary screen and the Slow Post primary screen. Each dot is
one peak; x is log2(high/low FC) and y is -log10(raw p). Grey, blue and magenta
denote p>=0.05, nominal p<0.05, and q<0.10, respectively, with q taking priority.
Old peaks with undefined or zero FC are excluded from the log-axis drawing
but retained in statistical tables and the BH family. Slow Post was chosen
for this display after observing the primary closing hit; all eight contexts
are reported above. Old FC is a binary-open proportion ratio; new FC is the
edgeR normalized-count ratio, so their exact magnitudes are not identical
measurement scales.

`native_primary_contexts.png` and `.pdf` show all six estimable native primary
contexts. Script 77 separately renders the mapped-only comparison; script 81
renders the final native-union figures. Script 82 annotates all 645 intervals
against unrestricted transcript TSS and gene bodies, independently of the
25-gene inclusion labels.

The agreed broader functional scope retains plausible individual candidates,
but does not establish a coherent opening-and-RNA-up programme. The old data
retain nominal CSRP3/MYOD1/CAV3 signals. The new data add a closing locus and
parameter-sensitive opening loci, without a completed matching RNA chain.
These results guide exploratory priority; they do not support presenting a
uniform COQ8A-driven chromatin-opening mechanism.

## Reproduction

```bash
python scripts/76_differentiation_fusion_analysis.py --work GSE240061_WORK
python scripts/77_validate_union_and_plot.py --old-work GSE208248_LINK_WORK --new-work GSE240061_WORK
python scripts/78_define_native_union_peaks.py --work GSE240061_WORK --gtf GENCODE_V48_GTF_GZ
Rscript scripts/79_native_union_donor_models.R . GSE240061_WORK results/differentiation_fusion25/native_GSE240061
python scripts/80_native_union_RNA_links.py --work GSE240061_WORK
python scripts/81_summarize_native_union.py --work GSE240061_WORK
python scripts/82_annotate_native_loci.py --gtf GENCODE_V48_GTF_GZ
```

The work directories contain the previously extracted matrices and metadata
from scripts 48, 61, 64 and 66. Underlying per-peak models are produced by
scripts 05, 48, 51 and 67–69. `manifest.json` specifies scope and inference;
`validation.json` records matrix checks. Private C2C12 data are not published
in this report or its source tables.

Optional script 83 cross-references the 25 symbols against owner-provided RNA
contrast results, with outputs required to be outside the public repository.
Day-versus-baseline contrasts retain their original meaning and are not renamed
as COQ8A overexpression effects. Rounded zero p/q values in the source are
flagged rather than interpreted as exact zero probabilities.

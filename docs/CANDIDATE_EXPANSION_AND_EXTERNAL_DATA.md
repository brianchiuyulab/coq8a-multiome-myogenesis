# Candidate comparison and independent differentiation references

## Objective and statistical units

The objective is a COQ8A-associated chromatin feature with a plausible
relationship to myogenic differentiation; MYOD1 is not a required outcome.
Opening and closing associations are both retained. The external 401-peak
selective-enhancer set is one biological scope, not a requirement that every
useful myogenic enhancer be exclusive to muscle in a reference panel.

The 401 peaks occupy 370 externally segmented enhancer blocks. For example,
three distinct target ATAC peaks can fall inside the same continuous external
state-4/state-5 enhancer block and contribute to one regional score. This changes
the feature being tested, not the number of nuclei, donors, or libraries.
It does not assert that gaps between measured ATAC peaks are accessible.

The target experiment still has two donor-derived lines, four libraries and
201 matched high/low pairs. The original p values test conditional differences
between paired nuclei, not differences across two independent donor means.
Biological source N is two. COQ8A >=3 UMI occurs in 206 of 32,977 QC-pass nuclei;
201 are matched. The remaining exactly-1-UMI pool contains 4,602 nuclei.
Thus the scarcity of the high-expression group, rather than extensive matching
loss, determines the present high/low comparison size. The source-level ratios
below give equal weight to the two stages within each source and are descriptive.

## Existing human candidates compared without changing the contrast

| Candidate and GRCh38 peak | High/low ATAC FC | Nominal ATAC p | Equal-stage FC, source 1 | Equal-stage FC, source 2 | Adjusted peak-open RNA count ratio | COQ8A high/low RNA p |
|---|---:|---:|---:|---:|---:|---:|
| MYOD1, chr11:17649919-17650798 | 1.630 | 0.02129 | 1.701 | 1.424 | 1.461 | 0.05104 |
| CAV3, chr3:8733438-8733968 | 1.765 | 0.05325 | 1.209 | 2.800 | 1.185 | 0.04316 |
| DMD, chrX:31220862-31221781 | 2.400 | 0.02006 | 2.589 | 1.822 | 1.118 | 0.000469 |

Peak-open RNA ratios compare RNA counts in peak-detected versus peak-undetected
nuclei, adjusting for technical depth, COQ8A expression and the fixed myogenesis
score. These are different contrasts from COQ8A high/low ATAC FC. Original count
models, RNA-link eligibility and full-screen ATAC q values are unchanged.
The focal ATAC results are exploratory nominal associations, not new FDR-positive
discoveries. CAV3 p=0.05325 is reported as such.

CAV3 and MYOD1 overlap HSMM strong enhancers without strong-state overlap in
the eight comparator tracks. The DMD peak also overlaps a HSMM strong enhancer,
but four comparator types have a strong state at that peak. Consequently it
remains in the 1,777-peak enhancer scope and is absent from the strict 401-peak
scope. This external annotation does not establish that DMD lacks a myogenic role.

The code preserves the complete existing 284 count-model follow-ups with new
annotations; it does not present that selected follow-up set as a new unbiased
whole-genome screen. It also preserves all five previously identified nominal
positive chains, including CREBRF, DMD, PFN2, HSPB8 and RTN2.

### Closing-region reassessment

The two peaks in the previously highlighted LSP1/TNNI2/TNNT3 neighborhood were
checked against the full previously computed cis-RNA candidate space. At the 5%
primary detection gate, fully adjusted correlations with LSP1 and TNNT3 are
small and nonsignificant; TNNI2 is not eligible at that gate. Its 1% sensitivity
does not provide consistent positive coupling either. The ATAC block-level
q=0.0130 remains numerically unchanged, but the nearest-gene labels do not yet
provide a persuasive functional target. This closing signal is retained without
promoting its nearby genes to established enhancer targets.

## New independent reference: GSE224489

GSE224489 provides C2C12 RNA-seq at GM/0, 12, 24, 48, 60 and 96 hours, with
three RNA samples at each time. It also provides two ATAC samples per condition
for GM and 60 hours, and three Hi-C samples per condition. These are separate
bulk assays, not paired RNA/ATAC measurements from single nuclei.

The GEO RNA count matrix was downloaded and analyzed with DESeq2 using
`design = ~time_h`. Genes with total count >=10 were retained; integer raw
counts, median-ratio size-factor estimation, dispersion estimation and Wald
contrasts for each later time versus GM were used. Default independent filtering
and BH adjustment are applied across genes within each time contrast. All-gene
results, normalized counts, sample design, size factors and R session information
are supplied. No outcome-dependent selection of replicates or time points was
performed. Mouse symbols are matched to unique uppercase symbols for annotation;
this is not a claim of verified one-to-one orthology for every gene in the merged
table. The focal CAV3/Cav3, DMD/Dmd and MYOD1/Myod1 comparisons are named homologs.

| Gene | 24h/GM RNA FC | 24h q | 48h/GM RNA FC | 48h q |
|---|---:|---:|---:|---:|
| Coq8a | 2.183 | 1.23e-6 | 5.885 | 3.26e-34 |
| Coq9 | 1.076 | 0.129 | 1.061 | 0.241 |
| Cav3 | 2.218 | 0.0919 | 16.833 | 1.17e-14 |
| Dmd | 5.174 | 0.175 | 33.262 | 0.000184 |
| Myod1 | 0.836 | 0.00149 | 1.539 | 1.57e-16 |
| Csrp3 | 4.980 | 6.78e-8 | 22.383 | 3.89e-30 |

Dmd has zero raw counts in all three GM samples and only 6, 4 and 7 raw counts
at 48h, rising to 101, 82 and 85 at 96h. Its modeled fold change is highly
baseline-sensitive and should not be used as a large-effect headline. Cav3 is
also low at baseline (1, 6 and 2 counts), with 49, 45 and 46 at 48h and 685, 645
and 681 at 96h. These patterns support later differentiation association, not a
COQ8A perturbation mechanism. Coq9 does not significantly increase at 24 or 48h;
the metabolic intermediate in the proposed pathway remains unmeasured here.

### Cross-species ATAC annotation

Author ATAC metadata explicitly specify mm10 and describe merged replicates
followed by MACS2 broad-peak calling. All 5,097 human target peaks were evaluated
using the UCSC hg38-to-mm10 chain. Start, midpoint and end-1 must map uniquely to
the same chain and strand, with a mapped span 0.5 to 2 times the human interval.
There are 1,367 consistent interval mappings. Unmapped or ambiguous intervals
remain unannotated; they are not treated as closed in mouse.

The CAV3 candidate maps to **mm10 chr6:112459366-112459903**. It overlaps author
ATAC peak calls in both GM and DM60h, with coverage fractions 0.708 and 1.000.
The earlier CSRP3 candidate also maps and overlaps both conditions. The focal
MYOD1 and DMD intervals do not pass these strict interval-mapping rules.
Coverage fractions describe geometric overlap with called peaks; they are not
accessibility fold changes or differential-ATAC p values. No stronger opening
claim is made from the change in peak-call coverage alone.

The mapper streams chain blocks to retain requested positions rather than
indexing the whole mammalian alignment in memory. Source coordinates are 0-based,
half-open; minus-strand target coordinates use qSize-qPosition-1. Its results were
checked against pyliftover on positive/negative chains with alignment gaps.
Input hashes are recorded in `mouse_ATAC_reference_manifest.json`.

## Candidate priority after adding external evidence

**CAV3 is the leading differentiation/fusion-oriented candidate for follow-up.**
Its combination of source-consistent ATAC direction, positive RNA coupling,
independent muscle enhancer annotation, mapped C2C12 ATAC interval, independent
differentiation-associated RNA induction and prior functional literature provides
a coherent reason to pursue it despite borderline nominal human ATAC p.
The literature reports induction during C2C12 differentiation and inhibition of
myotube formation after reducing caveolin-3. It does not establish that arbitrary
CAV3 overexpression will always enhance fusion or that COQ8A acts through CAV3.

**DMD is a complementary maturation-associated candidate** with a stronger
nominal high/low ATAC contrast and concordant RNA association. Its cross-species
regulatory interval has not been established, and low baseline RNA prevents
interpreting a very large external RNA ratio as a precise effect size.

MYOD1 remains a regulatory candidate; it is not the mandatory endpoint. CSRP3
remains differentiation-associated, but the earlier raw-count peak-RNA association
attenuated after state adjustment. No new composite p value is constructed by
combining these different assays, candidate selection steps or public/private
experiments. Private passage/time RNA comparisons remain outside this repository
and are not relabelled as COQ8A OE/KD experiments.

## Additional dataset inventory

| Dataset | Design | Appropriate role | Work completed |
|---|---|---|---|
| GSE240061 | Same-nucleus human RNA/ATAC; 6 participants, 12 biopsies; 4 exercise and 2 resting controls, each pre/post | Independent COQ8A-accessibility assessment within cell type, donor and time | RDS downloaded; paired counts, cell types, COQ8A groups and matching audited; differential-accessibility test not yet run |
| GSE224489 | Bulk C2C12 differentiation RNA time course, GM/60h ATAC and Hi-C | Independent differentiation relevance and orthologous-region annotation | RNA analyzed and author ATAC calls intersected as above |
| GSE191190, subseries GSE191188/191189; GSE145297 | Regeneration scRNA, separate MuSC ATAC after PGE2/aging and bulk RNA | Potential external functional-response reference | Metadata audited; not a same-nucleus COQ8A high/low test |
| GSE221736 | Muscle denervation, including limited multiome samples | Lower-priority alternative context | GEO reviewed; no new analysis |
| GSE327285 | Human chondrocyte dedifferentiation | Outside skeletal myogenesis scope | Excluded by tissue context |

GSE240061's GEO overall-design text mentions PBMCs, but the paper and series
summary identify vastus lateralis skeletal muscle. The full paper describes
rare satellite cells and no significant exercise-induced DARs in that population;
this is not a result about COQ8A. Actual cell-type counts, COQ8A detection and
paired matrices are now checked in `GSE240061_FEASIBILITY_RESULTS.md`; the old
strict high/low matching rule leaves only three baseline satellite pairs from
one donor. A six-donor MuSC replication claim is therefore not justified. Twelve
biopsies and multiple sequencing libraries must not be counted as twelve
independent donors. New cohorts are assessed for relevance and usable design,
not chosen on the basis of a favorable COQ8A association.

## Reproduction

```bash
Rscript scripts/56_external_c2c12_timecourse.R /data/GSE224489_C2C12_differentiation_all_samples_reads_count.txt.gz results/external_c2c12
python scripts/57_external_mouse_peak_annotation.py --reference /data/GSE224489 --out results/external_c2c12
python scripts/58_compare_candidate_evidence.py --cache /work/unstratified --out results/candidate_comparison
```

The RNA and two broadPeak files are in the [GSE224489 supplementary directory](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE224nnn/GSE224489/suppl/).
The mapping file is [UCSC hg38ToMm10.over.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToMm10.over.chain.gz).
The focal comparison and full existing follow-up inventory are in
`results/candidate_comparison`. RNA reference results and interval annotations
are in `results/external_c2c12`.

## Primary sources

- [GSE224489: C2C12 differentiation RNA, ATAC and Hi-C](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224489)
- [Rubenstein et al., Genome Research 2025: same-cell muscle multiome, GSE240061](https://pubmed.ncbi.nlm.nih.gov/40393809/)
- [Wang et al., Cell Stem Cell 2025: PGE2 and MuSC chromatin/regeneration](https://pubmed.ncbi.nlm.nih.gov/40513560/)
- [Targeted caveolin-3 reduction inhibits myotube formation in differentiating C2C12](https://www.sciencedirect.com/science/article/pii/S0021925819519007)
- [FAK signaling, caveolin-3 and normal myoblast fusion](https://pmc.ncbi.nlm.nih.gov/articles/PMC2710835/)
- [Differentiation-dependent expression of muscle adhesion-complex components](https://link.springer.com/article/10.1186/s13395-019-0218-x)

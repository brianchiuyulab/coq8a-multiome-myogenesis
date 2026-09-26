# Linking temporal accessibility to candidate gene regulation

## Scope

This extension evaluates all 410 externally defined dynamic regions in the
existing 221-gene search space: 191 early-opening, 34 middle-opening, four
late-opening and 181 closing peaks. Closing refers to the external differentiation
time course, independently of the direction in COQ8A-high nuclei. A region near
a gene is a candidate regulatory site, not a confirmed target assignment.

The main comparison remains TSS enrichment ≥3 and COQ8A ≥3 versus exactly
1 UMI: 201 within-library depth-matched pairs across four libraries from two
source lines. All four existing COQ8A/TSS combinations are retained. No new
COQ8A threshold, nucleus matching, peak geometry or annotation was introduced.

## Reproduction

```powershell
python scripts/32_bidirectional_peak_screen.py --tables results/tables --temporal results/temporal
python scripts/30_middle_cis_links.py --tables results/tables --temporal results/temporal --gtf ../data/public/gencode.v48.annotation.gtf.gz --h5-root ../data/public/GSE208248_processed --region-family opening-closing
python scripts/33_dynamic_evidence.py --tables results/tables --temporal results/temporal
```

The default `middle` mode of script 30 remains available. It writes the original
middle-specific filenames; `opening-closing` writes separate `dynamic_*` files.

Peak differences use the same two-sided exact paired binomial test. The expanded
410-peak family and phase-specific BH values are exported for each comparison.
The original 229-peak q values remain unchanged in a separate column. No
correction is applied across the sensitivity parameter combinations. All 916
original opening-peak tests and 16 module summaries were reconciled numerically.

Candidate cis links include all protein-coding genes with a transcript TSS within
500 kb, yielding 5,734 peak–gene pairs. The inherited RNA normalization,
within-library adjustment and eligibility rules are described in
[Middle cis links](MIDDLE_CIS_LINKS.md). Of these pairs, 1,022 were testable in the
main comparison: 492 early, 38 middle, nine late and 483 closing. Seventeen links
had BH q<0.05 across the 1,022 tests, all positive. Link q values and differential
accessibility q values test different hypotheses and are never substituted.

The gene-level RNA comparison comes from the original fixed 221-gene analysis.
Its statistic is a high-minus-low mean log1p(CP10K) difference, not a fold change.
Genes beyond that screen retain missing RNA differential statistics rather than
being assigned negative findings. Statistical inference remains exploratory at
the nucleus level; four libraries are not four independent donor replicates.

## Main candidate assessment

These illustrative loci were prioritized by biological relevance and concordance,
not by a new composite p value. Every tested locus and every cis candidate remain
available in `dynamic_evidence.tsv.gz`.

| Gene / phase | hg38 locus | ATAC ratio | ATAC p | ATAC q, 410 peaks | RNA difference | RNA p | RNA q, 221 genes |
|---|---|---:|---:|---:|---:|---:|---:|
| CSRP3 / early | chr11:19201752–19202603 | 2.4286 | 0.03088 | 0.86864 | +0.07162 | 0.01210 | 0.07427 |
| CAV3 / early | chr3:8733438–8733968 | 1.7647 | 0.05325 | 0.86864 | +0.08381 | 0.04316 | 0.16298 |
| LDB3 / middle | chr10:86674717–86675570 | 1.9375 | 0.03154 | 0.86864 | +0.10858 | 0.01205 | 0.07427 |
| TNNC1 / middle | chr3:52446715–52447619 | 1.7647 | 0.06599 | 0.86864 | +0.16321 | 0.06826 | 0.22085 |
| SPARC / closing | chr5:151698875–151699764 | 0.6667 | 0.00404 | 0.41402 | +0.00190 | 0.96389 | 0.99542 |
| GNAO1 / early | chr16:56341742–56342652 | 2.0870 | 0.000346 | 0.12795 | Not detected | — | — |

### Differentiation cofactor: CSRP3

The peak midpoint is 148.5 bp from an annotated CSRP3 transcript TSS. The 2.43
ratio comes from 17 versus seven accessible nuclei; library directions are two
positive, one negative and one tied. CSRP3 RNA was detected in 37/402 nuclei;
no library met the existing eligibility rules for this peak–RNA link. Missing
link statistics are therefore uninformative, not evidence of no connection.
Across the four settings ATAC ratios are 1.61–2.43, with three nominal p<0.05.

CSRP3/MLP has direct C2C12 evidence for enhancing differentiation and cooperating
with MyoD ([Kong et al.](https://pubmed.ncbi.nlm.nih.gov/9234731/)). This makes it a
functional candidate, with sparse public ATAC support. Isoforms require care:
the MLP-b variant has an opposing differentiation phenotype
([Vafiadaki et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC4416226/)).

### Fusion-related candidate: CAV3

The candidate midpoint is 98 bp from an annotated CAV3 transcript TSS. Main
ATAC counts are 30 versus 17, and the peak–RNA link is r=0.10168, p=0.05955,
expanded-link q=0.51575. Main gene RNA is nominally higher, but neither main
ATAC nor RNA passes its BH threshold.

Under the separately reported TSS≥3, COQ8A≥2 comparison, the peak–RNA link is
r=0.09142, p=0.000195, expanded-link q=0.00720, while the ATAC high/low ratio is
1.2093 with p=0.16772. This supports a regulatory association but not a significant
COQ8A accessibility effect in that setting. Across all four settings the ATAC
direction is positive; ratios span 1.21–1.82.

In C2C12, suppressing caveolin-3 prevented fusion despite retained differentiation
marker expression ([Galbiati et al.](https://pubmed.ncbi.nlm.nih.gov/10514527/)).
This is directly relevant to fusion; it does not establish that increasing CAV3
above its normal range is sufficient to improve fusion.

### Middle module

The pre-existing 34-peak aggregate remains 1.22293 (p=0.001785, ten-module
q=0.01785). LDB3 has a positive main peak–RNA association (r=0.13081, p=0.01512),
which is nominal: q=0.32193 in the expanded link family, and the previously
reported middle-only q=0.19152 is preserved. TNNC1's middle peak has only two
eligible libraries for linking; an early TNNC1 peak links strongly to RNA but
has nearly unchanged COQ8A accessibility (1.0588, p=0.8877).

### Closing branch: SPARC

The leading closing peak is less accessible in every library under the main
comparison (54 versus 81 nuclei; ratio 0.6667). However, its SPARC RNA link is
r=0.01310, p=0.7980, and SPARC RNA does not decrease in COQ8A-high nuclei.
Therefore the current data do not support the complete chain “peak closes →
SPARC transcription falls → fusion improves”.

Functional evidence is context-dependent: SPARC overexpression inhibited C2C12
differentiation in [Petersson et al.](https://pubmed.ncbi.nlm.nih.gov/23670848/),
whereas [Munk et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC6440632/) reported
evidence for SPARC supporting myogenesis. It is not classified as a universal
anti-myogenic factor.

### Core-TF and proximity checks

An early MYOG-linked region, chr1:203115894–203116803, has r=0.29097 and link
q=1.609×10−5, but its COQ8A accessibility ratio is only 1.0714 (p=0.7982).
Thus a reliable peak–RNA link does not imply COQ8A-dependent opening.
GNAO1 RNA is undetected in all 402 main nuclei. Its strong-looking neighboring
ATAC effect cannot establish GNAO1 as the expressed regulatory target.

## Decision

For mechanistic prioritization, CSRP3 addresses differentiation cofactor activity
and CAV3 addresses fusion more directly. LDB3/TNNC1 remain downstream muscle
programme candidates. No selected locus yet combines an FDR-positive COQ8A
ATAC difference, an FDR-positive cis link and an FDR-positive same-nucleus gene
RNA difference. The strongest aggregate result remains the middle-opening
module; the most phenotype-relevant early candidates have weaker statistical
support. Opening times describe the external regions, not a demonstrated
CSRP3→CAV3→LDB3 temporal or causal sequence.

## Optional private RNA cross-check and figure

`35_crosscheck_existing_rna.py` joins an existing RNA differential table without
refitting it, using case-insensitive gene symbols. This is not human–mouse peak
orthology. Pass a private output directory. The original passage experiment is
not an OE experiment, and unchanged original contrast FDRs are retained.

```powershell
python scripts/35_crosscheck_existing_rna.py --evidence results/temporal/dynamic_neighborhood_evidence_main.tsv --rna-contrasts PATH_TO_PRIVATE_RNA_RESULTS.csv.gz --out PRIVATE_OUTPUT
python scripts/34_plot_candidate_bridge.py --evidence PRIVATE_OUTPUT/dynamic_candidate_evidence.tsv --tables results/tables --out PRIVATE_OUTPUT
```

Figure legend: Selected illustrative loci and matched gene-level RNA associations.
Panel A shows low/high accessible-nucleus percentages, exact paired p and BH q
across all 410 dynamic peaks. Panel B shows mean paired log1p(CP10K) differences
and 95% t intervals, with BH across the existing 221 RNA genes; undetected RNA
is labeled explicitly. Panel C shows existing P33/P11 RNA ratios at D2 and D6
on a log2 axis; these are 2^model log2FC from log2(normalized count+1) data.
Peak counts in A and gene RNA in B are different endpoints. The plot is an
evidence comparison, not proof that each peak controls the named gene.

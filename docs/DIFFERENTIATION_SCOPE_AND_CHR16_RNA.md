# Differentiation scope and chr16 RNA-target follow-up

## External gene definitions

The frozen MSigDB set `GOBP_MYOBLAST_DIFFERENTIATION` (GO:0045445)
intersects the existing 221-gene universe in 17 genes:

ANKRD2, BOC, CSRP3, IFRD1, IGF1, IGFBP3, ITGB1, MAPK12, MAPK14,
MEF2C, MYF5, MYF6, MYOD1, MYOG, NOTCH1, RB1, TGFB1.

These genes define 345 distinct GSE208248 peaks shared across all four libraries
under the existing 100-kb neighborhood rule. No temporal-phase filter is required.
The differentiation and fusion intersections share ITGB1, MAPK14, MYOD1 and MYOG;
their union contains 25 genes. MYF5 is included by the differentiation definition,
while CAV3 is retained by the fusion definition. The union is an external functional
scope, not a newly discovered list of statistically significant genes.

Source memberships are frozen in `reference/muscle_subprogrammes.json` and
`results/functional/functional_genes.tsv`. The union table is
`results/chr16_RNA_followup/differentiation_fusion_membership.tsv`.

## chr16:1311478–1312392: what RNA does it track?

This is the GSE240061 interval with normalized ATAC FC=3.007 and candidate-wide
q=0.0224 under configuration C198. It entered the screen through the CACNA1H
100-kb neighborhood. Unrestricted GENCODE v48 annotation places it within UBE2I
and across the TSS of a protein-coding UBE2I transcript. That positional result
does not by itself establish a target RNA.

We therefore exported all author RNA features corresponding by exact gene name
to genes with an annotated transcript TSS within 100 kb of the peak midpoint.
Ten genes were present as RNA features: UBE2I, BAIAP3, TSR3, GNPTG, UNKL,
TPSD1, TPSAB1, TPSB2, TPSG1, CACNA1H. Missing exact-name features are recorded,
not treated as zero expression; cross-version remapping of unnamed lncRNAs is
outside this follow-up.

### Same-nucleus link analysis

The primary context is the same 490 high/low pairs in C198, Fast myonuclei at
post sampling, six donors. C197/C199 repeat matching at calipers0.10/0.30.
All author-retained Fast Post nuclei form an additional context.

Within each donor we correlate binary accessibility with log1p RNA CP10k,
controlling log1p RNA depth, log1p detected ATAC peaks, mitochondrial fraction
and log1p COQ8A CP10k. A second model adds the fixed myogenesis RNA score.
None of the ten local genes belongs to that score's gene list.

An eligible donor-target comparison requires at least eight nuclei, at least
three accessible and three inaccessible nuclei, and target RNA detected in at
least max(3, 1% of nuclei). Donor correlations are combined using equal-weight
Fisher-z means and a two-sided t test across donors. BH is applied to all
evaluable local genes within each context/model. Four genes have sufficient
information for a multi-donor combined test; their donor counts are retained.

For C198, the UBE2I link is **r=0.00562, p=0.884, q=0.884, N=6 donors**.
Adding the programme score yields r=0.00558, p=0.886. Other matched calipers
and the complete Fast Post context also show no convincing UBE2I link.
TSR3, GNPTG and UNKL do not show a significant local-target-adjusted link.
CACNA1H does not meet the RNA detection support requirement.

### COQ8A-high versus low RNA

RNA counts for the same matched nuclei were aggregated by donor and group.
The model follows the preceding RNA analysis: edgeR robust QL, `~donor+group`,
raw whole-transcriptome library-size offsets, and an outcome-blind filter of
total counts>=10 and counts in at least two donors. Dispersion estimation uses
the original RNA candidate background plus additional local genes. BH for the
local follow-up is over all count-supported local targets (four in C198).

| Local gene | C198 RNA FC | p | Local-target q |
|---|---:|---:|---:|
| UBE2I | 1.0053 | 0.9827 | 0.9827 |
| TSR3 | 0.6824 | 0.5315 | 0.9827 |
| GNPTG | 0.9408 | 0.7560 | 0.9827 |
| UNKL | 0.8473 | 0.4864 | 0.9827 |

UBE2I RNA FC is 0.989–1.005 across the three depth-cap95 matched settings.
The ATAC difference has therefore not been connected to increased expression
of an identifiable local RNA target by these analyses.

## Biological relevance and decision

UBE2I encodes UBC9. A primary C2C12 study reported impaired terminal
differentiation after Ubc9 knockdown, without corresponding loss of MyoD or
myogenin expression/activity:
[Riquelme et al., 2006](https://pubmed.ncbi.nlm.nih.gov/16631162/).
This supplies a plausible biological reason to inspect UBE2I, but it does not
replace the currently unsupported link between this peak and UBE2I RNA.

For the present COQ8A-to-myoblast-differentiation question, the chr16 interval
remains an ATAC observation but is **deprioritized as the mechanistic lead**.
No claim is made that CACNA1H or UBE2I is activated downstream of COQ8A.
Owner-provided C2C12 expression checks remain in private files outside this
repository. Their findings are not public study results.

## Reproduction

```bash
Rscript scripts/73_export_locus_RNA.R AUTHOR_RDS . LOCUS_WORK
python scripts/74_test_locus_RNA_links.py --work MULTIOME_WORK --locus LOCUS_WORK --out results/chr16_RNA_followup
Rscript scripts/75_locus_RNA_group_models.R . MULTIOME_WORK LOCUS_WORK results/chr16_RNA_followup
```

`MULTIOME_WORK` contains the aligned counts and matching outputs from scripts61–66;
`LOCUS_WORK` is outside the public repository. Public source tables contain
gene availability, per-donor detection, all evaluated correlations and RNA-group
models. These are targeted follow-up tests after the ATAC screen, not an
independent validation cohort.

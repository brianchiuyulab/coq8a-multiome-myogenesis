# Search-space and target-selection audit

## Biological search space

The analysis starts from two independent knowledge-based gene sets: [MSigDB Hallmark Myogenesis](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_MYOGENESIS.html) and [Reactome Myogenesis](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/REACTOME_MYOGENESIS.html). Membership was checked against the public MSigDB pages on 2026-09-26 and is fixed in [`myogenesis_221_gene_sources.tsv`](../reference/myogenesis_221_gene_sources.tsv).

| Set operation | Genes |
|---|---:|
| Hallmark Myogenesis source set | 200 |
| Hallmark genes present in the GSE208248 RNA matrix | 198 |
| Reactome Myogenesis genes present in the RNA matrix | 29 |
| Genes shared by the two measured sets | 6 |
| Unique genes entering the analysis | **221** |

The two Hallmark genes absent from the input matrix are **DENND2B** and **MYL11**. The six genes shared by the measured sets are MAPK12, MEF2A, MEF2C, MEF2D, MYF6 and MYOG. The `named_core_TF` column labels seven familiar regulators (MYOD1, MYOG, MYF5, MYF6, MEF2A, MEF2C and MEF2D); it does **not** filter the 221-gene search. MRF and MEF2 summaries are secondary contextual results in the tables. In particular, MYF5 enters through Reactome and is part of the complete search even though it is absent from the Hallmark set.

## Primary ATAC test family

The search uses hg38 canonical chromosomes, anchor-peak width 200–2000 bp, blacklist exclusion, ≥50% reciprocal overlap to corresponding library peaks, occurrence in ≥3 of four libraries, and position within ±100 kb of a transcript TSS of any of the 221 genes. The definition does not use motif scores or COQ8A high/low effects. It yields **7,699 peak–gene pairs** involving **7,132 distinct anchor peaks**. The gene-region summary uses only corresponding peaks present in all four libraries.

Of the 7,132 candidate peaks, **6,871** have at least 20 pooled open observations in the primary matched nuclei and receive a paired-test BH q value. The other 261 remain in the output table with no q value. All 221 gene-region summaries and all 6,871 eligible peaks are shown in the primary results before ranking a target. None has ATAC q<0.05.

## Exploratory MYOD1 follow-up

MYOD1 was nominated only after the complete search. Its 32 nearby candidate peaks narrow to 19 corresponding peaks present in all four libraries. Among all tested peak–RNA links, 10 MYOD1 peaks meet the positive-link, q<0.05, four-evaluable-library rule when links are learned at TSS≥3. Five of these have normalized maximum MRF PWM score <0.95. Repeating link learning at TSS≥2 yields 11 positive four-library links and **six** below that score cutoff.

The <0.95 motif-score rule is an exploratory, uncalibrated threshold; it does not establish that MRF proteins cannot bind those regions. The six-peak set is also downstream of a target ranking that used COQ8A effect direction. It therefore cannot serve as an outcome-independent confirmation. The complete gate/count/region-set grid is in [`myod1_locus_summary.tsv`](../results/tables/myod1_locus_summary.tsv), and each tested link is in [`peak_gene_links_tss2.tsv`](../results/tables/peak_gene_links_tss2.tsv) or [`peak_gene_links.tsv`](../results/tables/peak_gene_links.tsv).

The six-peak TSS≥2-link set gives a high/low open-fraction ratio of **1.280** in the TSS≥3, COQ8A≥3 versus 1 UMI effect test (201 matched pairs). It gives **1.110** in the primary COQ8A≥2 versus 1 contrast, and falls to about **1.192** after excluding pairs that touch the highest 2.5% RNA doublet-score tail. These are exploratory sensitivity estimates; the complete primary ATAC families remain without an FDR-positive result.

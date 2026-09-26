# Extended exploratory screen

## Scope and methods

The original 20-setting table covered five MYOD1 region definitions, two TSS-enrichment gates and two COQ8A contrasts. It did not exhaust the gene-wide search, link-discovery contrasts or TSS plotting windows.

This extension applies each region definition to every gene in the fixed 221-gene Hallmark/Reactome union. Effect tests use existing depth-matched nuclei: COQ8A >=2 versus 1 or >=3 versus 1 UMI, each at TSS enrichment >=2 or >=3. Link discovery is now evaluated under all four settings. The original default link tables retain their original settings.

Eight position-based definitions are evaluated: all four-library-common peaks; promoter-proximal peaks within 0.5, 1, 2 or 5 kb of the nearest annotated transcript TSS; and distal peaks beyond 1, 2 or 5 kb, within the original 100-kb search. Four link-based definitions are evaluated at each link-discovery setting: all positive links, links with maximum MRF PWM score <0.95, links with score >=0.95, and distal positive links. Positive links require partial r>0, link-family q<0.05 and four evaluable libraries. The motif split is a numerical score split, not a measurement of protein occupancy.

Each definition is evaluated by gene and over six unions: all 221 genes, Hallmark, Reactome, MRF, MEF2 and the seven MRF/MEF2 regulators. Programme scores measure the fraction of distinct selected peaks open per nucleus. Every result records the effective peak and gene counts. A programme label identifies the starting search space, not a guarantee that every starting gene survives filtering.

Thirteen additional external gene sets from MSigDB are frozen in `reference/muscle_subprogrammes.json` with source URLs and download hashes. Their intersections with the original 221 genes are documented in `subprogramme_sources.tsv`. They cover GO muscle fate, differentiation, fusion, assembly and contraction, MYOD perturbation signatures, SIN3A and MEF2 target sets. These are functional definitions; they are not automatically early-stage definitions. MSigDB gene sets are attributed to the Broad Institute, MIT and Regents of the University of California under [CC BY 4.0](https://www.gsea-msigdb.org/gsea/msigdb/license_terms.jsp). de la Serna signatures derive from mouse fibroblast MYOD-induction experiments represented as human symbols in MSigDB.

Seven symmetric TSS windows (±100, ±200, ±500, ±1000, ±2000, ±3000 and ±5000 bp) are evaluated from existing 100-bp insertion profiles in the 201 main pairs. These pooled counts are descriptive; no nucleus-level p value is inferred from them.

## Exploratory route to MYOD1

The new figure fixes effect estimation at COQ8A >=3 versus 1 UMI and TSS enrichment >=3 throughout, except the explicitly labelled sensitivity panel F. This comparison contains 201 matched pairs across four libraries from two source lines.

1. **RNA programmes:** higher COQ8A accompanies higher Hallmark Myogenesis, Reactome Myogenesis and MRF RNA scores; four-programme BH q values are 0.00990, 0.03158 and 0.00779.
2. **Complete ATAC screen:** all 221 gene regions are displayed. Their union contains 5,097 distinct common peaks and has a high/low open-fraction ratio of 1.0036 (nominal paired p=0.5109). The previously reported ±5-kb TSS insertion ratio is 1.0162. These are distinct measurements and neither describes a broad 1.2-fold increase.
3. **Candidate ranking:** links learned in the larger COQ8A >=2-versus-1 population at TSS>=2 are joined to the fixed >=3-versus-1 effect table. Candidates are ordered by the count of positive RNA-linked peaks with increased accessibility in at least three libraries, then total positive-link count, then gene-region effect. **MYOD1 ranks first**, with 11 positive links and four concordant accessibility links. TNNT2, TEAD4 and CKB each have three concordant links. The ranking rule is applied across the starting 221 genes.
4. **MYOD1 region:** its 19 common peaks yield 1.0934-fold opening; its 11 positive links yield 1.0985-fold. The six links with MRF score <0.95 yield 13.267% versus 10.365% open: **1.2800-fold**, +2.902 percentage points (95% paired-nucleus CI +0.358 to +5.446 points), nominal p=0.02558, four positive library directions. The five higher-score links yield 0.985-fold.
5. **Fixed-set sensitivity:** the same six peaks yield 1.2769-fold at TSS>=2 with COQ8A>=3 versus 1 (212 pairs, p=0.02522), and 1.1097–1.1207-fold under COQ8A>=2 versus 1.

Link discovery and effect estimation in this route use explicitly different populations. Learning links using only >=3-versus-1 nuclei retains no qualifying four-library MYOD1 link at either TSS gate. The same MYOD1 route is therefore not obtained when both stages are restricted to the smaller population. The alternative link tables and all rankings are included.

## Additional candidates and windows

The MYOD/SMARCA4 signature intersects the 221-gene search at nine genes. Their 13 common promoter-proximal peaks within 2 kb represent seven genes and yield **1.1483-fold** (20.742% versus 18.064%, p=0.00654, four positive libraries). Using 0.5, 1 or 5 kb yields 1.1404, 1.1250 and 1.1437-fold, with nominal p=0.02782, 0.03085 and 0.00241; every window has four positive library directions. This set originates from genes induced by MYOD whose induction depended on functional SMARCA4/BRG1 in the source experiment: [MSigDB source description](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/DELASERNA_TARGETS_OF_MYOD_AND_SMARCA4.html).

The GO positive-regulation-of-muscle-differentiation set intersects at 14 genes. TSS>=3, >=2-versus-1 link discovery retains 18 positive links at six genes. Applying MRF score <0.95 retains seven peaks at **MYOD1, CAV3 and LAMA2**. The main >=3-versus-1 effect is **1.2298-fold**, 14.0725% versus 11.4428% open, p=0.03557, three positive libraries and one exact tie. Numerical direction counts use a 1e-10 percentage-point tolerance to avoid counting floating-point cancellation as an increase. The preceding all-region and positive-link-only ratios are 0.9720 and 1.0824. The seven-peak result should not be attributed to the whole 14-gene set. [GO/MSigDB source](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/GOBP_POSITIVE_REGULATION_OF_MUSCLE_CELL_DIFFERENTIATION.html).

Under COQ8A>=3 versus 1 and TSS>=3, **MYL1** has a 1.5231-fold ratio over its complete seven common peaks (7.036% versus 4.620%, nominal p=0.00554, four positive libraries). Its six distal peaks yield 1.6071-fold (p=0.00357, four positive libraries). Thus 1.28 is not the maximum local effect in this dataset. The original 221-gene q value for the complete MYL1 region is 0.4079.

Across all 221 TSSs, shrinking the half-window from 5 kb to 100 bp gives ratios of 0.9794–1.0162. It does not amplify the overall TSS signal. A strongly filtered programme result reaches 1.2645-fold but represents only three peaks at three genes, not the full 221-gene programme and not MYOD1. Retained-gene/peak counts must accompany these summaries.

All effect-grid p values are nominal paired-nucleus statistics. Candidate and setting selection do not make them independent confirmation. Existing complete-screen multiplicity results are preserved; no correction across sensitivity settings is applied.

## Figure legend

**Figure_exploratory_route.** (A) RNA programme differences with paired-nucleus 95% confidence intervals; q values correct four RNA programmes. (B) All 221 common-peak gene-region effects, log2 open-fraction ratio against nominal paired p; MYOD1 and MYL1 are highlighted. (C) Leading candidates using TSS>=2, COQ8A>=2-versus-1 links and main-contrast effects. (D) Open fractions across nested 19-, 11- and six-peak MYOD1 sets. (E) Six-peak open fractions by library and in pooled matched nuclei. (F) Both COQ8A contrasts and both TSS QC gates applied to the fixed six-peak set. Four libraries represent two source lines.

## Reproduction and outputs

Add `--extended-exploration` to the documented `run_all.py` command. Scripts 07, 15 and 16 generate extra link tables, the complete grid and the PNG/PDF figure.

- `results/exploration/effect_grid.tsv`: all evaluated effects, effective counts, paired confidence intervals and library directions.
- `results/exploration/effects_by_library.tsv.gz`: library-specific values.
- `results/exploration/region_memberships.tsv.gz`: exact peak membership for every scored set.
- `results/exploration/candidate_ranking_grid.tsv`: all link/effect combinations; genes with no qualifying links have no numerical rank.
- `results/exploration/tss_window_grid.tsv`: seven windows per programme.
- `figures/exploration/Figure_exploratory_route.png` and `.pdf`: six-panel exploratory figure.
- `figures/exploration/Figure_subprogrammes.png` and `.pdf`: external subprogramme screen, MYOD/SMARCA4 promoter-window sensitivity and the seven-peak pro-differentiation subset.
- `results/exploration/subprogramme_sources.tsv`: external source-set sizes, exact intersections and references.

## Early differentiation question

The present subprogramme results do not determine when these regions open during differentiation. The next analysis is specified in [EARLY_ACCESSIBILITY_DESIGN.md](EARLY_ACCESSIBILITY_DESIGN.md): use an independent temporal ATAC reference to define early-opening regions, then assess those fixed regions in the COQ8A multiome. It has not yet been executed.

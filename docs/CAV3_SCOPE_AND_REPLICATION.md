# CAV3 prioritization and independent replication design

## Fixed target comparison

The target dataset is GSE208248: two donor-derived human muscle cell lines,
each profiled before and after seven days of differentiation. Four paired
RNA/ATAC libraries are available. The current analysis retains 32,977 nuclei
after QC and compares 201 matched pairs, within libraries, using TSS enrichment
>=3 and COQ8A RNA >=3 versus exactly 1 UMI. These p values quantify matched-nucleus
associations; four libraries are not four independent donors.

## Candidate rationale

The 221-gene myogenesis universe defines 5,097 peaks shared across all four
libraries. Independent HSMM strong-enhancer annotation retains 1,777 peaks.
Requiring no strong-enhancer annotation in eight nonmuscle comparator cell types
defines an optional 401-peak muscle-selective scope. It is not a universal filter:
a myogenic regulatory element can also function in another cell type.

Within these scopes, COQ8A-associated accessibility, effect direction across
sources, RNA coupling and independent differentiation evidence are inspected.
This is exploratory prioritization rather than a sequence of FDR-positive gates.
CAV3 is prioritized because multiple types of evidence support a differentiation
relationship, not because it passed the accessibility FDR threshold.

The focal GRCh38 interval is chr3:8733438-8733968. Its midpoint is 99 bp from the
nearest annotated CAV3 transcript TSS in the existing candidate annotation; the
interval therefore covers that TSS. This supports a promoter-proximal assignment,
without asserting that it regulates every CAV3 isoform. The same interval also
overlaps an external HSMM strong-enhancer state. Chromatin-state and transcript
annotations need not give mutually exclusive labels.

## Accessibility FDR scope audit

The focal peak is detected in 30/201 high and 17/201 low nuclei, giving FC=1.7647
and two-sided paired-binomial p=0.0532519141. Recalculated BH values are:

| Hypothesis family | Number of peaks | Focal q |
|---|---:|---:|
| Shared myogenesis candidates | 5,097 | 1.000 |
| HSMM strong enhancers | 1,777 | 1.000 |
| Muscle-selective enhancers | 401 | 1.000 |
| CAV3-neighboring muscle-selective peaks, post hoc scope audit | 8 | 0.426015 |

The last row is a post hoc diagnostic and is not presented as an independently
prespecified family. No listed scope yields CAV3 accessibility q<0.1. Small
peak-RNA association q values or external differentiation RNA q values answer
different hypotheses and cannot substitute for this accessibility q.

All eight CAV3-neighboring selective peaks combined yield FC=1.036 and
p=0.6797 in the existing region analysis. Thus the focal peak does not imply
uniform opening across the whole CAV3 neighborhood.

Independent C2C12 GSE224489 RNA rises at 48h versus growth medium (Cav3 FC=16.833,
q=1.17e-14). The mapped mouse interval overlaps author ATAC calls in both growth
medium and differentiation conditions. This is interval annotation, not a new
quantitative differential-ATAC result. Private passage/time RNA comparisons are
kept outside this repository and are not COQ8A perturbation contrasts.

## Dataset comparison

| Feature | GSE208248, current target | GSE240061, replication candidate |
|---|---|---|
| Assay | Same-nucleus 10x RNA + ATAC | Same-nucleus 10x RNA + ATAC |
| Material | Cultured human muscle cells | Human vastus lateralis biopsies |
| Independent biological sources | 2 donor-derived lines | 6 participants |
| Design | Each line before/after day-7 differentiation; 4 libraries | 4 exercise and 2 resting participants, each pre/post; 12 biopsies |
| Perturbation/time | In vitro differentiation | 40-min cycling or rest, post biopsy 3.5h later |
| Nuclei | 32,977 after the current analysis QC | 37,154 after author QC |
| Cell composition | Cultured myogenic populations | 14 annotated cell types, including fibers and rare satellite cells |
| Processed inputs | RNA/peak matrices, barcode metrics, fragments | Approximately 3.5-GB compressed RDS; approximately 21-GiB fragment archive |
| COQ8A analysis status | Current matched comparison completed | Metadata reviewed; paired matrices and COQ8A associations not yet inspected |

Counts after different QC procedures are descriptive and are not a direct
comparison of cell recovery. Technical assay accessions and repeated biopsies
must not be counted as additional independent participants. GSE240061 is an
exercise/tissue study, not another cultured differentiation time course.

The next feasibility step is to inspect the RDS RNA counts, ATAC assay and
barcode metadata while preserving author cell-type labels. Report counts and
COQ8A/CAV3 detection by participant, time and cell type before selecting an
analysis population. Baseline samples can provide a six-participant replication
context if enough usable nuclei exist; pre/post changes require repeated-measure
handling. The focal CAV3 interval and additional existing candidates can be
evaluated as an explicitly defined replication family. No favorable result is
assumed. Rare MuSC recovery may limit the usable biological sample count.

## Reproduction

```bash
python scripts/59_candidate_fdr_scope.py --gene CAV3
```

The output `results/candidate_comparison/CAV3_fdr_scope_audit.tsv` retains every
CAV3-neighboring peak in each declared scope. BH is calculated over the complete
scope before extracting CAV3 rows. Region-level results and external RNA results
are generated by the existing scripts documented in
`CANDIDATE_EXPANSION_AND_EXTERNAL_DATA.md`.

## Sources

- [GSE208248](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208248)
- [Human muscle cell model, iScience 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10123345/)
- [GSE240061](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE240061)
- [Rubenstein et al., Genome Research 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12212352/)
- [Caveolin-3 and normal myoblast fusion](https://pmc.ncbi.nlm.nih.gov/articles/PMC2710835/)

# Function-defined myogenic region analysis

## Question and external definition

Within the fixed Hallmark/Reactome myogenesis universe of 221 genes, are
differentiation- or fusion-associated dynamic regions differentially accessible
between COQ8A-high and COQ8A-low nuclei?

Four human MSigDB GO Biological Process sets define membership:

| Functional annotation | MSigDB set | Genes within 221 | Early-opening peaks | Middle-opening peaks | Late-opening peaks | Closing peaks |
|---|---|---:|---:|---:|---:|---:|
| Myoblast differentiation | GOBP_MYOBLAST_DIFFERENTIATION | 17 | 16 | 4 | 0 | 12 |
| Myoblast fusion | GOBP_MYOBLAST_FUSION | 12 | 15 | 2 | 0 | 15 |
| Positive regulation of muscle cell differentiation | GOBP_POSITIVE_REGULATION_OF_MUSCLE_CELL_DIFFERENTIATION | 14 | 13 | 0 | 0 | 8 |
| Negative regulation of muscle cell differentiation | GOBP_NEGATIVE_REGULATION_OF_MUSCLE_CELL_DIFFERENTIATION | 4 | 3 | 0 | 0 | 2 |

The union contains 30 genes. Of these, 22 have one or more of the previously
defined 410 dynamic regions in their fixed candidate neighborhoods. There are
54 distinct functional dynamic peaks: 24 early, five middle, no late, and 25
closing. Memberships overlap and are not additive. No gene is added manually.

Gene annotations and source hashes are archived in `reference/` and
`results/functional/analysis_manifest.json`. Myoblast differentiation is the GO
annotation, not a claim that every member exclusively controls the
myoblast-to-myotube transition. Positive and negative regulation annotations
may overlap and depend on biological context; CAV3 and NOTCH1 occur in both.
These annotations do not determine the sign of a peak's COQ8A effect.

## Region selection

1. Intersect each complete external gene set with the existing 221 genes.
2. Retain existing candidate peaks within 100 kb of annotated protein-coding
   transcript TSSs of these genes, with valid coordinates, blacklist exclusion,
   and the existing four-library correspondence. This is genomic proximity,
   not an inferred or validated enhancer-target assignment.
3. Intersect with the externally defined GSE109828 dynamic regions. Opening
   includes early, middle and late classes; closing is analyzed separately.
4. Write gene/peak memberships and their SHA-256 before loading COQ8A effects.

The analysis was developed after earlier exploration of this target dataset.
External annotations make membership reproducible and independent of current
effect values; they do not make this a retrospectively preregistered study.
Earlier broad GO and temporal analyses remain available. This is an additional
specified exploratory analysis, not a replacement for them.

## Fixed comparison and inference

QC and matched nuclei are reused without modification. The main comparison is
COQ8A RNA >=3 UMI versus exactly 1 UMI, TSS enrichment >=3, with 201 matched
pairs within four libraries from two source cell lines. Existing paired RNA and
ATAC depth matching and joint QC apply. Sensitivities retain TSS >=2/3 and
RNA >=2/3 versus exactly 1, independently of functional membership.

For each nucleus and module, the score is the fraction of the distinct member
peaks detected. FC is mean high-group score divided by mean low-group score.
Module P values are paired t tests on high-minus-low scores, with confidence
intervals for that difference. BH correction includes all eight functional
opening/closing modules within a QC/contrast setting. The sample-level
inferential unit is a matched nucleus pair; this is not a four-donor test.
Library effects are provided separately.

Individual peak P values retain the previous two-sided exact paired discordance
test. BH correction includes all 54 distinct functional dynamic peaks together
within each setting, including increased and decreased effects. Original
410-peak, 229-peak and phase-family q values remain in the output. There is no
correction across the four sensitivity settings; choosing the most favorable
setting would require reporting that choice as exploratory.

## Results

| Module, external opening regions | Peaks | Main FC | P | q, eight modules |
|---|---:|---:|---:|---:|
| Differentiation | 20 | 1.0242 | 0.7664 | 0.8759 |
| Fusion | 17 | 1.0448 | 0.6114 | 0.8759 |
| Positive regulation | 13 | 1.0039 | 0.9651 | 0.9651 |
| Negative regulation | 3 | 1.1333 | 0.4004 | 0.8759 |

No opening or closing module passes nominal P <0.05 in any of the four fixed
settings. The earlier 34-middle-peak module result is unchanged; it represents
a different complete temporal set, not the smaller functional subsets here.

Selected loci for interpretation, with all 54 results retained in the table:

| Nearby gene | Peak | Main FC | P | q, 54 peaks |
|---|---|---:|---:|---:|
| CACNA1H | chr16:1152563-1153401 | 0.2273 | 0.0009105 | 0.04917 |
| CSRP3 | chr11:19201752-19202603 | 2.4286 | 0.030884 | 0.71082 |
| CAV3 | chr3:8733438-8733968 | 1.7647 | 0.053252 | 0.71082 |
| MYOD1 | chr11:17696914-17697783 | 1.7059 | 0.072951 | 0.71082 |

CACNA1H is the only q <0.05 peak in this main 54-peak family. Accessibility is
5/201 versus 22/201 nuclei, with decreased effects in all four libraries. All
four sensitivity settings give decreased accessibility and nominal P <0.05,
but only the main setting gives q <0.05 (other q values: 0.0834, 0.3044,
0.3542). It lies approximately 121 bp from an annotated CACNA1H transcript TSS.
RNA detection is insufficient for a reliable peak-to-CACNA1H RNA link in the
existing link analysis. It is not evidence that reducing CACNA1H promotes fusion.

CSRP3 and CAV3 enter through differentiation and fusion annotations,
respectively; neither passes the 54-peak correction. This external narrowing
does not establish a generalized opening of fusion/differentiation regions.

## Reproduction

From the repository root, with the processed GEO H5 files available:

```bash
python scripts/37_functional_subsets.py freeze
python scripts/37_functional_subsets.py test --h5-root /path/to/GSE208248_processed
```

The script verifies the membership hash before testing and reconciles all 32
module FC values against existing single-peak high/low counts. Outputs include
complete gene memberships, peak memberships, library effects, matched scores,
four-setting module results and four-setting single-peak results.

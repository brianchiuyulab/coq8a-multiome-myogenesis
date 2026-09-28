# COQ8A-associated chromatin accessibility during myogenesis

## Current fixed analysis: Day 0, COQ8A >=2 versus 0

![Position-aligned ATAC profiles](figures/fixed_day0_D206/Figure_4_position_aligned_ATAC.png)
![Local ATAC tracks](figures/fixed_day0_D206/Figure_5_local_ATAC_tracks.png)

These fragment-based panels use the same matched nuclei as the peak tests:
TSS-centered profiles, peak-centered heatmaps, and highlighted locus tracks.
Scripts 94 and 96 regenerate them from indexed fragments. The owner's C2C12
cross-reference is kept outside this public repository.

[Complete methods, results and figure legends](docs/FIXED_DAY0_D206_WORKFLOW.md)
describe the fixed D206 workflow: external 221-to-25-gene scope, 565 peaks,
joint QC, 527 within-source matched pairs, all-region opening/closing tests,
and 170 local RNA target assessments. Six neighborhoods pass regional q<0.05;
four additional labels pass q<0.10. The nominated MYOD1 peak has a positive
MYOD1 RNA association, including in the exact matched population; CAV3 has
concordant ATAC/RNA group changes and a weaker population-dependent link.
The same analysis parameters apply to every candidate.

Code: `scripts/92_fixed_day0_workflow.py`, `scripts/93_plot_fixed_day0.py`.
Source tables: `results/fixed_day0_D206/`. PNG/PDF: `figures/fixed_day0_D206/`.
Earlier settings below are retained as analysis history, not substituted into
this fixed comparison.

![Complete screen](figures/fixed_day0_D206/Figure_1_screen.png)
![Selected loci](figures/fixed_day0_D206/Figure_2_selected_loci.png)
![Local RNA targets](figures/fixed_day0_D206/Figure_3_RNA_targets.png)

[All 25 neighborhoods across the Day-0 sensitivity grid](docs/DAY0_ALL25_INVENTORY.md)
reports opening and closing representatives, exact parameter settings, peak and
regional statistics, and the frequency of signals across evaluated settings.

## ADAM12-neighborhood RNA follow-up

[The complete selection and local RNA assessment](docs/ADAM12_LOCUS_FOLLOWUP.md)
traces the external 25-gene scope to the D232 opening candidate, then tests all
six local RNA features. The 4.2857-fold accessibility contrast is reproduced;
the regulated RNA remains unresolved. Public results and script 90 are included.

## Expanded Day-0 sensitivity grid

The [336-configuration exploration](docs/DAY0_SENSITIVITY_GRID.md) is complete:
312 configurations are evaluable across two sources. It tests opening and
closing, seven COQ8A definitions, two TSS thresholds, four matching calipers,
state matching and three explicitly defined populations. ADAM12-neighborhood
opening and several closing candidates emerge; their sensitivity and exact
population definitions are reported. These are selected exploratory results,
not a claim of across-grid FDR control or purified MuSC annotation.

![Parameter sensitivity](figures/day0_sensitivity_grid/Parameter_sensitivity.png)

## Completed Day-0 primary analysis

[Stage-specific results and methods](docs/TASK1_DAY0_RESULTS.md) now separate
undifferentiated cultures (507 matched pairs) from Day 7. The entire frozen
25-gene / 565-peak scope was retested. CSRP3 ranks first in the primary opening
screen, but regional q=0.514; the pooled q=0.0956 below is not a Day-0 result.
Scripts 86–87, complete source tables and new figures are included.

![Day-0 complete screen](figures/task1_stage_specific/Figure_1_Day0_complete_screen.png)

## Current focus: Task 1, GSE208248

The [Task 1 workflow](docs/TASK1_ANALYSIS_WORKFLOW.md) documents the complete
functional-scope analysis: 221 myogenesis genes -> 25 differentiation/fusion
genes -> 565 candidate peaks -> 25 regional tests -> CSRP3 follow-up. It specifies
joint QC, within-library matching, both multiplicity levels and same-pair RNA
comparisons. The reported exploratory setting is TSS >=3 and COQ8A >=2 vs 1.
GSE240061 exploration is retained below but is not the current main evidence.
Earlier temporal figures are separate analyses, not figures for this regional test.

## Local-signal follow-up and candidate decision

The [local-signal report](docs/LOCAL_SIGNAL_DECISION.md) tests sparse opening
within each of the 25 gene neighborhoods and compares complete versus partial
ATAC/RNA evidence. CSRP3 is the closest integrated exploratory route; a regional
q value is distinguished from individual-peak q. Scripts 84–85 and all candidate
tables are provided.

## Differentiation + fusion union — completed 2026-09-27

The agreed 25-gene union has now been analyzed in both datasets: 565 GSE208248
peaks and 645 native GSE240061 peaks (335 mapped plus 310 additional). See the [complete results and methods](docs/DIFFERENTIATION_FUSION25_RESULTS.md),
[source tables](results/differentiation_fusion25), and [figure](figures/differentiation_fusion25/union_peak_screen.png).
Scripts 76–82 reproduce this analysis. The preceding temporal and narrower
fusion analyses below retain their own scopes and are not the new union results.

## Completed CAV3 and GSE240061 exploration — 2026-09-27

The new [analysis and results](docs/REPLICATION_ANALYSIS_AND_RESULTS.md) include 288 external-scope comparisons,
448 replication settings, donor-level statistics, model stress tests and figures.
GSE240061 differential accessibility has now been analyzed. The sections below
retain the preceding analysis and its specific settings.

New figures: [external CAV3 scope](figures/replication/Figure_1_CAV3_functional_scope.png),
[CAV3 sensitivity heatmap](figures/replication/Figure_2_CAV3_sensitivity.png),
[six-donor follow-up](figures/replication/Figure_3_CAV3_donor_followup.png),
[distal candidate](figures/replication/Supplement_distal_candidate.png).
Regenerate with `python scripts/71_summarize_replication.py`.

Reproducible analysis of paired human RNA/ATAC measurements at regions defined
by an independent myoblast differentiation time course.

## Study design

**Biological scope → external temporal definition → COQ8A association → locus and RNA follow-up**

1. The measured Hallmark/Reactome myogenesis union contains **221 genes**.
2. Their annotated neighborhoods contain **5,097 peaks common to four GSE208248 libraries**.
3. Independent GSE109828 differentiation data define **191 early, 34 middle,
   four late and 181 closing regions**.
4. External GO annotations for differentiation, fusion and regulation define
   **54 functional dynamic regions**, including opening and closing regions.
5. Compare COQ8A-high (≥3 RNA UMI) with low (1 UMI), using joint TSS ≥3 QC and
   **201 depth-matched pairs within libraries**. Repeat the complete comparisons
   with the specified TSS/count sensitivities.

The target data comprise **two donor-derived lines, each undifferentiated and
day-7 differentiated**. The four libraries are not four independent donors.
This repository analyzes association in processed public data.

## Main figures

### Figure 1 — External definition and temporal-module association

The external time course defines the regions. The aligned COQ8A comparison
then tests each temporal class. The middle module has FC **1.223**, paired
p **0.00178**, and BH q **0.0178** across the ten temporal programme tests.

![Figure 1](figures/main/Figure_1_External_definition.png)

### Figure 2 — Complete functional-region screen

All **54 regions** are shown in the same external-time order: matched-nucleus
ATAC insertion profiles, within-library accessibility differences and pooled
fold changes. Each row is one genomic interval; nearby gene names can repeat.

![Figure 2](figures/main/Figure_2_Functional_accessibility.png)

### Figure 3 — Candidate loci and paired RNA

Four follow-up loci are shown with actual fragment profiles, accessibility
fractions and RNA from the same matched nuclei. CSRP3, CAV3 and MYOD1 have
positive pooled accessibility directions; CACNA1H illustrates a decrease.
The figure reports both p and q and does not imply that all four pass FDR.

![Figure 3](figures/main/Figure_3_Candidate_loci.png)

Vector PDFs accompany all PNGs. [Supplementary figures](figures/supplement)
contain the 221-gene TSS heatmap and the complete four-setting module comparison.

## Reproduce the figures

Python 3.13:

```bash
python -m pip install -r requirements.txt
python run_paper.py --mode figures
```

This command uses committed source tables and fragment profiles, regenerates
the figures, and validates their statistical inputs. No private RNA data or
multi-gigabyte downloads are needed for this route.

## Recompute from deposited processed data

Download the GSE208248 H5 matrices and barcode metrics listed in
[Data availability](docs/DATA_AVAILABILITY.md), and GENCODE v48. Then run:

```bash
python run_paper.py --mode analysis --h5-root /data/GSE208248 --gtf /data/gencode.v48.annotation.gtf.gz --external-raw /data/GSE109828 --external-work /work/GSE109828
```

The command downloads the external reference, reconstructs temporal membership,
tests the fixed contrasts and generates figures. It reuses the versioned
fragment-QC tables. Add `--fragments /data/fragments` with `pysam` installed to
rebuild figure fragment profiles. Full QC regeneration and all analytical
definitions are in [Methods](docs/METHODS.md).

## Documentation and source data

- [Candidate provenance and independent replication](docs/CANDIDATE_PROVENANCE_AND_REPLICATION_PLAN.md):
  how CAV3 entered follow-up, identical evidence columns for every candidate,
  and the GSE240061 replication design.
  [Complete candidate audit](results/candidate_provenance).
  [GSE240061 measured feasibility](docs/GSE240061_FEASIBILITY_RESULTS.md)
  reports actual group sizes, matching and regional coverage.

- [Independent enhancer intersections](docs/EXTERNAL_ENHANCER_INTERSECTION.md):
  human HSMM strong-enhancer annotation, MYOD1 ChIP and core MRF neighborhoods.
  [Coordinate/time-course figure](figures/external_enhancers/MYOD1_external_annotation_and_time.png)
  and [source tables](results/external_enhancers).

- [Temporally unstratified exploration](docs/UNSTRATIFIED_LINK_EXPLORATION.md): all
  5,097 shared candidate peaks, nominal effects and peak–RNA count-model follow-up.
  [Exploration figure](figures/unstratified/Unstratified_peak_RNA_exploration.png)
  and [complete source tables](results/unstratified).

- [Middle-region target-link follow-up](docs/MIDDLE_LINK_REASSESSMENT.md): expanded
  nucleus population, count models, stratified bootstrap and independent MYOD1
  annotation. [Follow-up figure](figures/link_followup/Middle_peak_RNA_associations.png)
  and [source tables](results/link_followup) are separate from the main high/low tests.

- [Enhancer-region exploration](docs/ENHANCER_REGION_EXPLORATION.md): sample-size
  audit, external cell-context selectivity, and region-level tests without a MYOD
  occupancy filter. [Broad results](results/enhancer_regions) and
  [selective-region results](results/enhancer_selective_regions).
- [Candidate expansion and external datasets](docs/CANDIDATE_EXPANSION_AND_EXTERNAL_DATA.md):
  source-level comparisons, independent C2C12 RNA/ATAC reference, and additional
  dataset inventory. [Candidate evidence](results/candidate_comparison) and
  [external C2C12 results](results/external_c2c12).
- [Methods](docs/METHODS.md): samples, preprocessing, thresholds, matching, statistics and provenance.
- [Figure legends](docs/FIGURE_LEGENDS.md): panel definitions, signal units and statistical families.
- [Results](docs/RESULTS.md): current numerical findings.
- [Code availability](docs/CODE_AVAILABILITY.md) and [Data availability](docs/DATA_AVAILABILITY.md).
- [Figure source data](results/figure_source), [functional results](results/functional),
  [temporal results](results/temporal) and [sample/scope inventories](results/design).

Private C2C12 passage/time results are analyzed separately and are not deposited
here. Generic code for that comparison is provided with an external-output requirement.

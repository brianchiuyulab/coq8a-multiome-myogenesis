# TSS overview and focused testing

## TSS-aligned figure

`figures/overview/Figure_221_gene_TSS_heatmap.png/pdf` shows all 221 fixed
Hallmark/Reactome genes in the same 201 matched high/low nucleus pairs. Each
row is a gene-feature TSS, oriented by strand, with a ±5 kb window and 100 bp
bins. Rows are sorted by pooled high+low insertion intensity, identically in
both panels. Colors use one shared absolute scale in insertions per 100 nuclei
per 100 bp, with saturation at the pooled 99.5th percentile; there is no per-row
z-score. Gaussian smoothing (sigma one 100 bp bin) is for display only.

The figure is promoter context. It does not display all 5,097 distinct candidate
peaks: many candidate peaks are distal within the larger ±100 kb transcript-TSS
search space. It does not nominate candidates or alter the tests.

Reproduction:

```powershell
python scripts/36_plot_tss_heatmap.py --tables results/tables --out figures/overview
```

### Coordinate correction, 2026-09-27

The previous fragment-profile script applied +4/−5 bp to Cell Ranger ARC
fragment endpoints, although these endpoints are already Tn5-adjusted according
to [10x documentation](https://www.10xgenomics.com/support/cn/software/cell-ranger-arc/latest/analysis/outputs/fragments-file).
The duplicate offset was removed and all four fragment libraries were re-read.
Gene identities, alignment TSSs, bin widths and matched nuclei were unchanged.
The whole-window totals changed from 57,556/56,638 to 57,558/56,634 high/low
endpoints; the descriptive ratio is now 1.016315 (previously 1.016172).
Both the main Figure 1 and the simplified blue heatmap use the corrected table.
Historical exploratory TSS-window summaries remain outputs of the prior
coordinate implementation; they are not used to nominate the current candidates.

Differential peak counts, temporal modules, peak–RNA links and their p/q values
are calculated from the H5 peak matrices and were not recalculated or changed
by this display correction. The release validator passes after the correction.

## Why large ratios can have high q values

CSRP3's selected early peak has 17 versus seven accessible nuclei. Among 201
matched pairs, 14 are high-only and four low-only, giving the exact two-sided
p=0.0308838. It ranks sixth among the original 229 opening peaks. CAV3 has 30
versus 17 accessible nuclei, 26 versus 13 discordant pairs, and p=0.0532519.
Their 229-peak BH q values are both 0.812979; phase-only q values are both
0.782394. The larger 410-peak extension gives q=0.868644 for both.
Large relative effects at low counts need not yield precise estimates. BH
depends on the full ordered p-value distribution, not fold change alone.
Sensitivity settings have never been included as additional BH discoveries.

## Defensible focused analyses

1. Use externally defined functional groups (e.g. direct fusion regulators from
   a perturbation screen) and external promoter/enhancer evidence to specify
   which regions to test. Freeze this rule before examining the new subset's
   COQ8A test outcomes, retain all qualifying regions, and report the subset as
   an exploratory extension because the broader outcomes have already been seen.
2. Use suitable independent filtering of low-information features, with the
   filter/test combination justified under the null; blindly treating every
   pooled-count filter as independent is not sufficient. This may remove sparse
   attractive candidates rather than rescue them. See
   [Bourgon et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC2906865/).
3. Where the scientific question concerns a programme, predefine module-level
   primary tests and adjust across modules; individual regions can be presented
   as exploratory follow-ups. A significant module does not confer significance
   on every constituent region. The existing middle module result is an example
   of a distinct aggregate estimand, not a reclassification of nonsignificant peaks.
4. Raw p values may be reported transparently for exploratory candidates. They
   are not FDR-controlled discoveries. A genuinely prespecified single primary
   test in independent follow-up data can be assessed without a genome-wide
   correction; post hoc selection of two favorable peaks from this same screen
   does not create such a design.

Changing algorithms or repeatedly shrinking gene sets until significance is
obtained is not an independent validation. No reviewer acceptance is guaranteed
by a particular procedure. Functional relevance, effect uncertainty and test
families must remain explicit in Methods and figure legends.

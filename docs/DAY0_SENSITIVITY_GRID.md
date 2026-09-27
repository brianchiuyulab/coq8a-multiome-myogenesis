# Day-0 sensitivity exploration: local opening and closing

## Result and relation to the preceding analysis

The expanded grid identifies both opening and closing candidates in
undifferentiated GSE208248 cultures. The earlier Day-0 analysis tested only two
COQ8A definitions at TSS>=3 and one matching caliper. Its negative regional FDR
result did not establish absence of signals under other definitions.

The strongest reported opening route is an ADAM12-neighborhood peak in
myogenic-marker-positive nuclei: FC=4.2857, exact individual-peak p about
3.40e-5, q565=0.01918; regional opening p=0.0000965, q25=0.0024125. Both sources
increase. This is an exploratory choice from the sensitivity grid, not a
prospectively selected result or an independent validation.

Strict PAX7-detected nuclei also yield closing candidates. No opening region
in that subset has both regional q<0.10 and a top-peak increase in both sources.
Consequently myogenic-marker-positive, PAX7-detected and all-culture results
must not be collapsed into a single claim about purified MuSCs.

## Complete grid, fixed before target tests

The unchanged biological search space is the 25 differentiation/fusion genes
and their 565 distinct four-library-common candidate peaks. Only the two
undifferentiated libraries GSM6339597 and GSM6339601 are tested. Day 7 is not
included. The same raw matrices and archived non-TSS joint-QC filters are used;
the extra TSS>=2 nuclei were re-extracted from the deposited H5 matrices.

| Axis | Values |
|---|---|
| TSS enrichment | >=2, >=3 |
| Population | All undifferentiated QC-passing nuclei; myogenic-marker detected; PAX7 detected |
| COQ8A definition | >=2 vs 1; >=3 vs 1; >=4 vs 1; detected vs zero; >=2 vs zero; >=2 vs <=1; upper vs lower quartiles of CP10k among positive nuclei |
| Matching depth caliper | 0.05, 0.10, 0.20, 0.30 in each log1p depth dimension |
| Matching covariates | RNA/ATAC depth; depth plus myogenesis state |

This is 2 x 3 x 7 x 4 x 2 = **336 configurations**. Of these, 312 retain pairs
in both sources and are evaluated; 24 have no pairs in at least one source and
remain explicitly unavailable. No minimum of 20 nuclei per source is imposed.

Myogenic-marker detected means at least one RNA UMI in PAX7, MYF5, MYOD1 or
MYOG. PAX7 detected means at least one PAX7 RNA UMI. These are transparent
expression-based sensitivity populations, not new validated cell-type labels;
marker non-detection does not prove a nonmyogenic identity. Selection on
expression also changes the population being compared.

State is mean log1p CP10k over the broad myogenesis genes excluding all 25
tested scope genes. State matching requires a within-library standardized
score difference <=0.25 SD in addition to both depth calipers. The 100-neighbor
greedy matching procedure is retained, without replacement and within library.
State matching asks a conditional question and can attenuate either state
confounding or a biological association mediated by state; it is not an
automatic verdict that the depth-only result is invalid.

COQ8A zero means no detected RNA UMI, not proven biological absence. It defines
a different comparison from exactly one UMI and is not substituted silently.

## Statistics and precision

Every evaluable configuration tests all 565 peaks and all 25 neighborhoods.
Exact two-sided paired-binomial p values are BH-adjusted over 565 peaks.
Regional maximum standardized paired-difference statistics are tested for
opening and closing separately, preserving correlation among peaks during
within-pair label swaps. BH is across 25 regions within direction/configuration.

All configurations initially use 20,000 permutations. The top three concordant
configurations per direction, plus the retained baseline, were then rerun at
2,000,000 permutations. Additional precision was applied to the best concordant
setting per population/direction at TSS>=3, for nine refined configurations in
total. Each refinement uses its configuration-specific original seed
20260928 + configuration index; no seed is searched. Full initial and refined
tables, Monte Carlo intervals and chosen-configuration lists are retained.

The reported q values describe within-configuration test families. They do not
control the selection of the best result across this exploratory grid. Repeated
settings reuse nuclei and are not independent replication. N is two source
lines; paired-nucleus tests are not donor-level population tests.

## Representative results after high-precision refinement

| Config | Population | Direction / neighborhood | COQ8A high vs low | Depth caliper / state match | Pairs | Top peak FC | Regional p | Regional q25 | Peak q565 |
|---|---|---|---|---|---:|---:|---:|---:|---:|
| D232 | Myogenic-marker detected | Opening / ADAM12 | >=3 vs 1 | 0.05 / no | 80 | 4.2857 | 0.0000965 | 0.0024125 | 0.01918 |
| D206 | All D0 culture nuclei | Opening / MYOD1 | >=2 vs 0 | 0.30 / no | 527 | 1.5098 | 0.0003185 | 0.0079625 | 0.03031 |
| D205 | All D0 culture nuclei | Closing / NOTCH1 | >=2 vs 0 | 0.20 / yes | 493 | 0.3968 | 0.0000705 | 0.0017625 | 0.00460 |
| D325 | PAX7 detected | Closing / MYOG | >=2 vs <=1 | 0.20 / yes | 80 | 0.1429 | 0.0001320 | 0.0033000 | 0.02263 |

All four use TSS>=3. FC is the high/low fraction of nuclei with a positive
ATAC entry. A value below one indicates closing. Regional q and individual
peak q answer different hypotheses; both are shown here.

| Neighborhood / exact peak | Source 1 high / low open (pairs) | Source 2 high / low open (pairs) |
|---|---|---|
| ADAM12: chr10:126373783-126374503 | 12 / 4 (35) | 18 / 3 (45) |
| MYOD1: chr11:17653349-17654252 | 71 / 52 (226) | 83 / 50 (301) |
| NOTCH1: chr9:136586858-136587780 | 13 / 30 (201) | 12 / 33 (292) |
| MYOG: chr1:203091609-203092395 | 1 / 10 (20) | 2 / 11 (60) |

These are neighborhood labels from the external functional scope. This screen
does not establish that the named gene is the RNA target of the focal peak,
that its RNA must change immediately, or that closing increases differentiation.
No RNA response or link p value was used as an inclusion gate.

## Are these isolated settings?

- **ADAM12 opening:** At TSS>=3, myogenic-marker detected, >=3 vs 1 and
  depth-only matching, calipers 0.05/0.10/0.20/0.30 give top-peak FC
  4.29/3.44/3.44/3.44 and regional q approximately
  0.0024/0.011/0.031/0.025. The same peak increases in both sources throughout.
  TSS>=2 gives the same 80 pairs at caliper 0.05 and a refined q=0.002325.
  Adding state matching weakens most settings; at caliper 0.10, FC=3.67 and
  screening q=0.0687, but the other state-matched settings have q>0.10.
- **MYOD1 opening:** In all D0 nuclei with >=2 vs zero and depth-only matching,
  all four TSS>=3 calipers give FC about 1.49-1.51 and regional q<0.05.
  The corresponding state-matched regional tests have q>0.79. This is not
  established as a state-independent COQ8A association.
- **NOTCH1-neighborhood closing:** In all D0 nuclei with >=2 vs zero, all eight
  TSS>=3 depth/state matching combinations give regional q<0.10 and FC about
  0.40-0.46, with both sources decreasing. This is the more consistent
  conditional closing route among the displayed candidates.
- **PAX7-positive MYOG-neighborhood closing:** State-matched calipers
  0.10/0.20/0.30 give FC 0.19/0.14/0.27 and q about 0.076/0.0033/0.082.
  Depth-only settings do not reach q<0.10. The result is population- and
  matching-dependent and is not a demonstration that MYOG RNA is suppressed.

Unrefined adjacent-setting regional values above are 20,000-permutation
estimates; the focal rows in the representative table use 2,000,000. Exact
individual-peak p values do not depend on Monte Carlo precision.

For context, a PAX7-positive TGFB1-neighborhood opening configuration has
screening q=0.076 and FC=2.67, but source 1 is flat (3 vs 3) while source 2
increases (21 vs 6). It is retained in the full grid rather than reported as
a concordant opening lead.

## Figure

figures/day0_sensitivity_grid/Parameter_sensitivity.png (and vector PDF)
shows four displayed neighborhoods at TSS>=3. Columns vary the depth caliper;
rows add or omit state matching. Color encodes the log2 FC of the region's
leading peak and text gives FC and regional q25. The leading peak can change
between settings; the full tables supply its coordinate. Red denotes opening
and blue closing, not significance. Complete results, not only the four
displayed candidates, remain in results/day0_sensitivity_grid/.

## Reproduction and audit

```bash
python scripts/88_day0_sensitivity_grid.py --work DAY0_GRID_WORK --h5 GSE208248_H5_DIRECTORY
python scripts/89_summarize_day0_grid.py --work DAY0_GRID_WORK
```

The initial grid and design are written before testing. Script 88 re-extracts
the two libraries from H5 for the expanded QC set, saves every pairing and
support count, and fits the screen/refinement. Script 89 completes subgroup
precision, validates all exact peak p/q and all regional BH values, verifies
that baseline D170 reproduces the preceding 507 matched pairs exactly, and
checks that each high/low pairing is disjoint and without replacement.
validation.json records successful checks. This does not replace or invalidate
the earlier restricted-setting results; it answers the broader sensitivity
request that had not previously been completed for Day 0.

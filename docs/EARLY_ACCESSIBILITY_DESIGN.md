# Independent temporal definition of early myogenic accessibility

## Question

Does higher COQ8A expression accompany increased accessibility at regulatory elements that normally open early as myoblasts enter differentiation? The proposed mechanistic chain is COQ8A/CoQ and metabolic changes, altered chromatin regulation, and improved entry into myogenic differentiation. The current observational contrast tests the accessibility association in that chain.

## External reference and rationale

[Pliner et al., Molecular Cell 2018](https://doi.org/10.1016/j.molcel.2018.06.044), the original Cicero study, measured human skeletal muscle myoblast ATAC at 0, 24, 48 and 72 hours after differentiation induction. The study identified dynamically opening elements, ordered them by accessibility timing and examined early-opening MYOD-associated sites and PBX1/MEIS1. [Author data inventory](https://cole-trapnell-lab.github.io/cicero-release/data/) identifies GSE109828 and its processed sci-ATAC datasets. Temporal coverage differs between its experiments and must be read from the sample metadata before assigning replication units.

This reference defines accessibility timing independently of COQ8A expression or accessibility effects in GSE208248. Early-opening elements may regulate genes expressed later. Consequently, an early-RNA-only list would discard some plausible priming elements. The MYOD/SMARCA4 gene set already explored here derives from a 24-hour induction experiment and must not be renamed an early-opening signature.

## Planned analysis

1. Retrieve published element-level dynamic-accessibility and timing annotations where available, plus coordinates, genome build and cell metadata. If timing annotations are unavailable, reproduce time-associated accessibility profiles from processed counts and use the reported temporal ordering, checking it against actual sampling times.
2. Define opening, closing and stable elements using only the external reference. Rank opening elements by onset or half-maximum timing; display continuous timing and early/middle/late groups. Explore the earliest 20%, 25% and one-third as explicitly labelled timing thresholds. These fractions are relative differentiation stages, not measured hour boundaries.
3. Harmonize genome coordinates and intersect external elements with the existing GSE208248 consensus peaks. Retain mapping and reciprocal-overlap provenance. Start with elements linked to or in the fixed candidate region space of the 221 genes; report attrition. Use external Cicero links preferentially for distal target assignment, avoiding selection based on the COQ8A effect.
4. Compare COQ8A high/low at those fixed early-, middle- and late-opening sets using the existing joint-QC, within-library matched nuclei. Evaluate >=2 versus 1 and >=3 versus 1 under TSS>=2 and >=3. Analyse undifferentiated and differentiated libraries separately as well as pooled, retaining source-line identity.
5. Compare the COQ8A-associated effect between early and late sets directly, rather than interpreting one small p value and one large p value as proof that the effects differ. Account for the differing baseline accessibility and peak counts. Report each timing stratum's fold ratio, percentage-point difference, uncertainty and library direction.
6. Within an associated temporal module, rank its externally linked genes and visualize the relevant loci. MYOD1 is one possible target; it is not required to be the winner.

## Figures

External-reference heatmap: peak accessibility across differentiation time, with rows ordered by activation timing. COQ8A comparison: effects on early, middle and late opening sets, separated by undifferentiated/differentiated library state. Locus panels: early-opening elements and external gene links around selected candidates. Effect-size and timing-threshold sensitivity plots accompany these panels.

## Status

Reference study and public accession verified. The external temporal annotations have not yet been imported or mapped, so no COQ8A early-opening result is claimed. Completed functional-subprogramme results remain in `EXTENDED_EXPLORATION.md`.

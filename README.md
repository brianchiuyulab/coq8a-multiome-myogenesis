# COQ8A-associated chromatin accessibility

## Current analysis: undifferentiated GSE208248 cultures

**External myogenesis scope → joint RNA/ATAC QC → COQ8A comparison → local RNA candidates**

- Hallmark/Reactome myogenesis: 221 genes and 5,097 shared candidate peaks.
- GO myoblast differentiation/fusion intersection: 25 genes, 565 distinct peaks.
- Day 0 only; TSS enrichment ≥3; COQ8A ≥2 versus 0 raw RNA UMI.
- Within-source depth matching: 527 pairs (226 + 301), two biological sources.
- Complete peak and regional tests, followed by 170 local peak–RNA assessments.

[Methods, results and figure legends](docs/FIXED_DAY0_D206_WORKFLOW.md) ·
[Computational audit](docs/FIXED_DAY0_AUDIT.md) ·
[Code availability and reproduction](docs/CODE_AVAILABILITY.md) ·
[Data availability](docs/DATA_AVAILABILITY.md)

## Main results

Four neighborhoods show increased accessibility at regional q<0.05: MYOD1,
CSRP3, CACNA1H and CAV3. BOC and NOTCH1 neighborhoods show decreased
accessibility. Names identify genomic neighborhoods; local RNA models assess
candidate targets separately. All comparisons use the same fixed definition.

| Neighborhood | Direction | Representative peak FC | Regional q |
|---|---|---:|---:|
| MYOD1 | Opening | 1.510 | 0.00796 |
| CSRP3 | Opening | 1.717 | 0.0135 |
| CACNA1H | Opening | 2.474 | 0.0229 |
| CAV3 | Opening | 1.316 | 0.0401 |
| BOC | Closing | 0.233 | 0.0302 |
| NOTCH1 | Closing | 0.459 | 0.0302 |

MYOD1 has a positive local RNA association in both eligible and matched nuclei.
CAV3 has concordant ATAC/RNA group changes; its link is weaker in matched nuclei.
The supplied tests concern nucleus-level associations within two sources.

## Figures

### 1. External scope and complete peak screen

![Candidate screen](figures/fixed_day0_D206/Figure_1_screen.png)

### 2. Six primary accessibility loci

![Source-specific accessibility](figures/fixed_day0_D206/Figure_2_selected_loci.png)

### 3. Local RNA candidate assessment

![Local RNA associations and group changes](figures/fixed_day0_D206/Figure_3_RNA_targets.png)

### 4. TSS- and peak-aligned ATAC profiles

![Positional ATAC heatmaps](figures/fixed_day0_D206/Figure_4_position_aligned_ATAC.png)

### 5. Representative genomic tracks

![Local ATAC tracks](figures/fixed_day0_D206/Figure_5_local_ATAC_tracks.png)

[All 25 regions](figures/fixed_day0_D206/Figure_S1_all_regions.png) and
[secondary loci](figures/fixed_day0_D206/Figure_S2_secondary_loci.png) accompany
the main panels. Editable PDF/SVG and 400-dpi PNG versions are provided.

## Reproduce

```bash
python -m pip install -r requirements.txt
python run_fixed_day0.py --mode figures
```

This regenerates the current figure set from committed source tables.
Full matrix-based recomputation and raw-input audit commands are in
[Code availability](docs/CODE_AVAILABILITY.md).

- Code: `run_fixed_day0.py`; scripts 92–98 and their documented upstream dependencies.
- Data: `results/fixed_day0_D206/`.
- Figures: `figures/fixed_day0_D206/`.
- Audits and checksums: `results/fixed_day0_D206/audit/`.

[Parameter sensitivity](docs/DAY0_SENSITIVITY_GRID.md) and
[earlier analysis routes](ANALYSIS_HISTORY.md) are retained separately.
Private C2C12 data and derived figures are outside this repository.

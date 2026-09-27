"""Summarize all 25 neighborhoods across the completed Day-0 grid."""
from pathlib import Path
import pandas as pd

root = Path(__file__).resolve().parents[1]
out = root / 'results/day0_sensitivity_grid'
x = pd.read_csv(out/'regional_results.tsv.gz', sep='\t')
rows = []
for (gene, mode), allrows in x.groupby(['gene','mode']):
    d = allrows[allrows.both_sources_same_direction]
    best = d.sort_values(['p','q25','config']).iloc[0].to_dict()
    best.update(n_evaluable=len(allrows), n_concordant=len(d),
                n_concordant_p05=int((d.p<.05).sum()),
                n_concordant_q01=int((d.q25<.1).sum()),
                minimum_q25=d.q25.min(),
                minimum_q_config=d.sort_values(['q25','p','config']).iloc[0]['config'])
    rows.append(best)
summary = pd.DataFrame(rows)
assert len(summary)==50 and summary.gene.nunique()==25
summary.to_csv(out/'all25_directional_summary.tsv',sep='\t',index=False,na_rep='NA')
pop={'all_culture':'All D0','myogenic_detected':'Myogenic-marker+','PAX7_detected':'PAX7+'}
parts=['# All 25 Day-0 neighborhoods: sensitivity inventory',
       'All 336 configurations were attempted; 312 were evaluable. This report summarizes existing results, without rerunning or changing the grid.',
       'For each gene and direction, the representative is the smallest regional p among settings whose top peak changes in that direction in both sources. FC and peak p refer to that top peak, not the whole neighborhood. Regional p accounts for searching peaks within the neighborhood; q25 is BH across the 25 neighborhoods in that setting and direction. The separately reported minimum q can come from a different setting.',
       'Counts of settings are sensitivity descriptions, not independent replications. Settings reuse nuclei. Best-setting p/q do not account for selection across the grid. Missing FC means a zero low-group detection denominator, not a measured infinite biological effect. Different directions or settings may select different peaks. Overlapping neighborhoods (including MYF5/MYF6) can share a peak.',
       'Full machine-readable results, including both peak and regional statistics, are in `results/day0_sensitivity_grid/all25_directional_summary.tsv`.']
for mode in ['opening','closing']:
    parts += ['## '+mode.title(),
    '| Neighborhood | Config / population / COQ8A contrast / TSS / caliper / matching | Pairs | Top peak FC | Peak p | Region p | Region q25 | Settings p<.05 / q<.1 |',
    '|---|---|---:|---:|---:|---:|---:|---:|']
    for _,r in summary[summary['mode']==mode].sort_values('p').iterrows():
        fc=f'{r.top_peak_FC:.3g}' if pd.notna(r.top_peak_FC) else 'NA (low=0)'
        setting=f'{r.config}; {pop[r.population]}; {r.contrast}; {r.TSS}; {r.caliper}; {r.matching}'
        parts.append(f'| {r.gene} | {setting} | {r.n_pairs} | {fc} | {r.top_peak_p:.4g} | {r.p:.4g} | {r.q25:.4g} | {r.n_concordant_p05} / {r.n_concordant_q01} |')
parts += ['## MYOD1 interpretation',
          'At TSS>=3 in all Day-0 cultures, COQ8A>=2 versus zero and depth-only matching gives the same peak chr11:17653349–17654252 at every caliper (0.05/0.10/0.20/0.30), FC 1.495/1.490/1.510/1.510 and regional q25 0.0100/0.01437/0.00875/0.007962. Additional state matching weakens regional q to 0.7978–0.9028. This is a worthwhile state-associated opening candidate, not proof of a state-independent effect or a confirmed MYOD1 RNA target. The same coordinates are not the historical chr11:17649919–17650798 focal peak.',
          'Opening q<0.1 occurs in at least one concordant setting for 12 neighborhood labels: ADAM12, MYOD1, CSRP3, CACNA1H, MYH9, CAV3, BOC, NOS1, TGFB1, IGF1, MYF5 and MYF6. These are not necessarily 12 independent loci. Large FC values should be read alongside the detection counts in the complete peak table.',
          'ANKRD2, MEF2C and ITGB1 have no concordant regional opening p<0.05 in this grid. That does not establish no chromatin association: ITGB1 has closing evidence. No gene should be described as universally unchanged simply from failure of its opening test.']
(root/'docs/DAY0_ALL25_INVENTORY.md').write_text(('\n\n'.join(parts)+'\n').replace('|\n\n|','|\n|'),encoding='utf-8')
print(summary[['gene','mode','n_concordant_p05','n_concordant_q01']].to_string(index=False))

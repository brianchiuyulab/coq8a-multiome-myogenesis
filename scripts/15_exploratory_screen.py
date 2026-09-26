"""Evaluate a finite exploratory grid over all 221 myogenesis genes.

This extends the original MYOD1-only grid. Region definitions are applied to
every gene before effects are compared. Outputs retain every tested setting.
Programme scores are the fraction of unique selected peaks open per nucleus.
P values describe paired nuclei; four libraries derive from two source lines.
"""

import argparse
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse, stats

from multiome_core import h5_for, read_barcodes


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--h5-root', type=Path, required=True)
    p.add_argument('--tables', type=Path, required=True)
    p.add_argument('--genes', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--extra-gene-sets', type=Path)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    genes = pd.read_csv(a.genes, sep='\t')
    candidates = pd.read_csv(a.tables / 'candidate_peak_gene.tsv.gz', sep='\t')
    common = candidates[candidates.n_libraries == 4].copy()
    mapping = pd.read_csv(a.tables / 'consensus_peak_map.tsv.gz', sep='\t').set_index('peak')
    pairs = pd.read_csv(a.tables / 'matched_pairs.tsv.gz', sep='\t')
    groups = {
        'All_221': set(genes.gene),
        'Hallmark': set(genes.loc[genes.Hallmark_Myogenesis, 'gene']),
        'Reactome': set(genes.loc[genes.Reactome_Myogenesis, 'gene']),
        'MRF': {'MYOD1', 'MYOG', 'MYF5', 'MYF6'},
        'MEF2': {'MEF2A', 'MEF2C', 'MEF2D'},
        'Core_TF7': set(genes.loc[genes.named_core_TF, 'gene']) if 'named_core_TF' in genes else {'MYOD1','MYOG','MYF5','MYF6','MEF2A','MEF2C','MEF2D'},
    }
    if a.extra_gene_sets:
        external=json.loads(a.extra_gene_sets.read_text(encoding='utf-8'))
        provenance=[]
        for name,source in external.items():
            available=set(source['geneSymbols']) & set(genes.gene)
            if available:
                groups[name]=available
            provenance.append(dict(programme=name,n_source_genes=len(source['geneSymbols']),
                                   n_genes_in_221=len(available),genes_in_221=';'.join(sorted(available)),
                                   source_url=source.get('msigdbURL',''),pmid=source.get('pmid','')))
        pd.DataFrame(provenance).to_csv(a.out/'subprogramme_sources.tsv',sep='\t',index=False)
    regions = [('all_common', 'none', 'none', common)]
    for width,label in [(500,'05kb'),(1000,'1kb'),(2000,'2kb'),(5000,'5kb')]:
        regions.append(('promoter_'+label,'none','none',common[common.nearest_tss_bp<=width]))
    for width,label in [(1000,'1'),(2000,'2'),(5000,'5')]:
        regions.append(('distal_'+label+'to100kb','none','none',common[common.nearest_tss_bp>width]))
    link_tables = {}
    for gate, suffix in [('TSS_ge_2', '_tss2'), ('TSS_ge_3', '')]:
        for contrast in ['2plus_vs_1', '3plus_vs_1']:
            cs = '' if contrast == '2plus_vs_1' else '_' + contrast
            links = pd.read_csv(a.tables / f'peak_gene_links{suffix}{cs}.tsv', sep='\t')
            link_tables[gate, contrast] = links
            selected = links[(links.q_all_links < .05) & (links.partial_r > 0) & (links.n_libraries == 4)]
            selected = selected.merge(common[['peak','gene']], on=['peak','gene'], validate='one_to_one')
            for rule, frame in [
                ('positive_link', selected),
                ('linked_MRFscore_lt095', selected[selected.mrf_max_score < .95]),
                ('linked_MRFscore_ge095', selected[selected.mrf_max_score >= .95]),
                ('linked_distal_2to100kb', selected[selected.nearest_tss_bp > 2000]),
            ]:
                regions.append((rule, gate, contrast, frame))
    metadata, memberships = [], []
    for rule, gate, contrast, frame in regions:
        selections = [('gene', name, sub) for name, sub in frame.groupby('gene')]
        selections += [('programme', name, frame[frame.gene.isin(members)]) for name, members in groups.items()]
        for scope, name, sub in selections:
            peaks = sorted(set(sub.peak))
            if not peaks:
                continue
            col = len(metadata)
            metadata.append(dict(set_id=col, rule=rule, link_gate=gate, link_contrast=contrast,
                                 scope=scope, name=name, n_peaks=len(peaks), n_genes=sub.gene.nunique(),
                                 represented_genes=';'.join(sorted(set(sub.gene)))))
            memberships.extend(dict(set_id=col, peak=peak) for peak in peaks)
    meta = pd.DataFrame(metadata)
    members = pd.DataFrame(memberships)
    members.to_csv(a.out / 'region_memberships.tsv.gz', sep='\t', index=False)
    pooled, library_rows = {}, []
    for gsm, libpairs in pairs.groupby('gsm'):
        barcodes = sorted(set(libpairs.high_barcode) | set(libpairs.low_barcode))
        with h5py.File(h5_for(a.h5_root, gsm)) as h5:
            _, P, _, peak_names = read_barcodes(h5, barcodes)
        bidx = {b:i for i,b in enumerate(barcodes)}
        pidx = {str(p):i for i,p in enumerate(peak_names)}
        local = members.copy()
        target = local.peak if gsm == 'GSM6339597' else local.peak.map(mapping[gsm + '_peak'])
        local['local_peak'] = target.map(pidx)
        assert local.local_peak.notna().all()
        local = local.drop_duplicates(['set_id', 'local_peak'])
        denom = local.groupby('set_id').size()
        matrix = sparse.csr_matrix((1 / local.set_id.map(denom).to_numpy(),
                                   (local.local_peak.astype(int), local.set_id)),
                                  shape=(len(peak_names), len(meta)))
        values = (P @ matrix).toarray()
        for (gate, contrast), ps in libpairs.groupby(['gate','contrast']):
            hi = values[[bidx[b] for b in ps.high_barcode]]
            lo = values[[bidx[b] for b in ps.low_barcode]]
            pooled.setdefault((gate,contrast), []).append((hi,lo))
            summary = meta.copy()
            summary['gsm'], summary['gate'], summary['contrast'] = gsm, gate, contrast
            summary['n_pairs'] = len(ps)
            summary['high_open_pct'], summary['low_open_pct'] = 100*hi.mean(0), 100*lo.mean(0)
            summary['delta_pp'] = 100*(hi-lo).mean(0)
            summary.loc[summary.delta_pp.abs()<1e-10,'delta_pp']=0.0
            summary['fold_open'] = np.divide(hi.mean(0),lo.mean(0),out=np.full(len(meta),np.nan),where=lo.mean(0)>0)
            library_rows.append(summary)
        print(gsm, 'exploratory sets', len(meta), flush=True)
    libraries = pd.concat(library_rows, ignore_index=True)
    libraries.to_csv(a.out/'effects_by_library.tsv.gz',sep='\t',index=False)
    summaries=[]
    for (gate,contrast), arrays in pooled.items():
        hi=np.concatenate([x[0] for x in arrays]); lo=np.concatenate([x[1] for x in arrays])
        d=hi-lo; n=len(d); h=hi.mean(0); l=lo.mean(0)
        summary=meta.copy()
        summary['gate'],summary['contrast'],summary['n_pairs']=gate,contrast,n
        summary['high_open_pct'],summary['low_open_pct']=100*h,100*l
        summary['fold_open']=np.divide(h,l,out=np.full(len(meta),np.nan),where=l>0)
        summary['delta_pp']=100*d.mean(0)
        se=stats.sem(d,axis=0); margin=stats.t.ppf(.975,n-1)*se
        summary['delta_ci_low_pp'],summary['delta_ci_high_pp']=100*(d.mean(0)-margin),100*(d.mean(0)+margin)
        summary['p_pair']=stats.ttest_1samp(d,0,axis=0).pvalue
        ss=libraries[(libraries.gate==gate)&(libraries.contrast==contrast)]
        summary['positive_libraries']=summary.set_id.map(ss.groupby('set_id').delta_pp.apply(lambda x: int((x>1e-10).sum())))
        summaries.append(summary)
    pd.concat(summaries,ignore_index=True).to_csv(a.out/'effect_grid.tsv',sep='\t',index=False)
    # Retain the original ranking rule and apply it to every link/effect setting.
    effects=pd.read_csv(a.tables/'candidate_peak_effects_pooled.tsv.gz',sep='\t')
    gene_effects=pd.read_csv(a.tables/'gene_effects_pooled.tsv',sep='\t')
    ranks=[]
    for (lg,lc),links in link_tables.items():
        links=links[(links.q_all_links<.05)&(links.partial_r>0)&(links.n_libraries==4)]
        for (eg,ec),eff in effects.groupby(['gate','contrast']):
            joined=links.merge(eff[['peak','delta_pp','positive_libraries']],on='peak',how='left')
            joined['concordant']=(joined.delta_pp>0)&(joined.positive_libraries>=3)
            rank=joined.groupby('gene').agg(n_linked=('peak','size'), n_concordant=('concordant','sum'),mean_delta_pp=('delta_pp','mean')).reset_index()
            rank=genes[['gene']].merge(rank,on='gene',how='left').fillna(0)
            ge=gene_effects[(gene_effects.gate==eg)&(gene_effects.contrast==ec)&(gene_effects.modality=='ATAC')]
            rank=rank.merge(ge[['gene','difference']],on='gene',how='left',validate='one_to_one')
            rank=rank.sort_values(['n_concordant','n_linked','difference','gene'],ascending=[False,False,False,True])
            rank['rank']=np.where(rank.n_linked>0,np.arange(1,len(rank)+1),np.nan)
            rank['link_gate'],rank['link_contrast'],rank['gate'],rank['contrast']=lg,lc,eg,ec
            ranks.append(rank)
    pd.concat(ranks,ignore_index=True).to_csv(a.out/'candidate_ranking_grid.tsv',sep='\t',index=False)
    # Existing 100-bp insertion bins suffice for seven exact symmetric windows.
    profile=pd.read_csv(a.tables/'tss_fragment_profile_221.tsv.gz',sep='\t')
    windows=[]
    for width in [100,200,500,1000,2000,3000,5000]:
        sub=profile[(profile.bin_start_bp>=-width)&(profile.bin_end_bp<=width)]
        for name, names in groups.items():
            s=sub[sub.gene.isin(names)]
            high,low=int(s.high_cuts.sum()),int(s.low_cuts.sum())
            windows.append(dict(programme=name,half_window_bp=width,high_cuts=high,low_cuts=low,
                                fold_cuts=high/low if low else np.nan,n_pairs=201,
                                gate='TSS_ge_3',contrast='3plus_vs_1'))
    pd.DataFrame(windows).to_csv(a.out/'tss_window_grid.tsv',sep='\t',index=False)
    print('Wrote complete exploratory grid to',a.out,flush=True)


if __name__ == '__main__':
    main()

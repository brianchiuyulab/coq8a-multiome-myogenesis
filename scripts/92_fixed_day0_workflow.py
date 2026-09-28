"""Fixed D206 screen, source-level summaries, and unrestricted local RNA models."""
import argparse
import importlib.util
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests


def save(d,p):
    d.to_csv(p,sep='\t',index=False,na_rep='NA')


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--work',type=Path,required=True)
    ap.add_argument('--h5',type=Path,required=True)
    a=ap.parse_args()
    root=Path(__file__).resolve().parents[1]
    out=root/'results/fixed_day0_D206';out.mkdir(exist_ok=True)
    grid=root/'results/day0_sensitivity_grid'
    config=pd.read_csv(grid/'configurations.tsv',sep='\t').query("config=='D206'").iloc[0].to_dict()
    config.update(region_discovery_q=0.05,region_exploratory_q=0.10,local_RNA_window_bp=500000,
                  selection_history='Fixed after sensitivity exploration',RNA_required_for_ATAC_selection=False)
    (out/'design.json').write_text(json.dumps(config,indent=2)+'\n')
    regions=pd.read_csv(grid/'regional_results.tsv.gz',sep='\t').query("config=='D206'").copy()
    peaks=pd.read_csv(grid/'peak_screen.tsv.gz',sep='\t').query("config=='D206'").copy()
    pairs=pd.read_csv(grid/'matched_pairs.tsv.gz',sep='\t').query("config=='D206'").copy()
    assert len(peaks)==565 and len(regions)==50 and len(pairs)==527
    assert np.allclose(multipletests(peaks.p,method='fdr_bh')[1],peaks.q565)
    for _,d in regions.groupby('mode'):
        assert np.allclose(multipletests(d.p,method='fdr_bh')[1],d.q25)
    hits=regions[regions.q25<.1].sort_values(['mode','q25']).copy()
    hits['tier']=np.where(hits.q25<.05,'q<0.05','0.05<=q<0.10')
    save(regions,out/'all_regions.tsv');save(peaks,out/'all_peaks.tsv');save(pairs,out/'matched_pairs.tsv')
    save(hits,out/'selected_regions.tsv')
    selected=sorted(hits.top_peak.unique())
    candidates=pd.read_csv(root/'results/unstratified/candidate_pairs.tsv.gz',sep='\t')
    candidates=candidates[candidates.peak.isin(selected)].reset_index(drop=True)
    assert set(candidates.peak)==set(selected)
    save(candidates,out/'local_RNA_candidates.tsv')
    genes=sorted(candidates.gene.unique());gi={g:i for i,g in enumerate(genes)};pi={p:i for i,p in enumerate(selected)}
    spec=importlib.util.spec_from_file_location('extractor',root/'scripts/44_expand_middle_links.py')
    ext=importlib.util.module_from_spec(spec);spec.loader.exec_module(ext)
    mapping=pd.read_csv(root/'results/tables/consensus_peak_map.tsv.gz',sep='\t').set_index('peak')
    rows=[];counts=[];rna_sources=[];hh=[];ll=[];qcs=[]
    for gsm,pr in pairs.groupby('gsm'):
        full=pd.read_csv(a.work/f'{gsm}_meta.tsv',sep='\t')
        meta=full[full.pass_tss3].reset_index(drop=True)
        qcs.append(dict(gsm=gsm,TSS2_n=len(full),TSS3_n=len(meta),high_eligible=int((meta.COQ8A_umi>=2).sum()),low_eligible=int((meta.COQ8A_umi==0).sum()),pairs=len(pr)))
        local=[p if gsm=='GSM6339597' else mapping.loc[p,gsm+'_peak'] for p in selected]
        mat,_=ext.extract(next(a.h5.glob(gsm+'_*filtered_feature_bc_matrix.h5')),genes,local,meta.barcode.tolist())
        raw=mat[:len(genes)].toarray().T.astype(float)
        x=(mat[len(genes):].toarray().T>0).astype(float)
        cp=raw/meta.total_rna_umi.to_numpy()[:,None]*10000;y=np.log1p(cp)
        ix=pd.Index(meta.barcode);hi=ix.get_indexer(pr.high_barcode);lo=ix.get_indexer(pr.low_barcode)
        assert min(hi.min(),lo.min())>=0
        assert (meta.iloc[hi].COQ8A_umi>=2).all() and (meta.iloc[lo].COQ8A_umi==0).all()
        hh.append(cp[hi]);ll.append(cp[lo])
        for j,p in enumerate(selected):
            counts.append(dict(gsm=gsm,peak=p,pairs=len(pr),high=int(x[hi,j].sum()),low=int(x[lo,j].sum())))
        for j,g in enumerate(genes):
            rna_sources.append(dict(gsm=gsm,gene=g,pairs=len(pr),high_mean_CP10k=cp[hi,j].mean(),low_mean_CP10k=cp[lo,j].mean(),high_detected=int((raw[hi,j]>0).sum()),low_detected=int((raw[lo,j]>0).sum())))
        for population,ids in [('eligible_D0',np.arange(len(meta))),('matched',np.r_[hi,lo])]:
            m=meta.iloc[ids]
            for state in [False,True]:
                cov=[np.ones(len(ids)),np.log1p(m.total_rna_umi),np.log1p(m.total_open_peaks),np.log1p(m.coq_CP10k)]
                if state:cov.append(m.state_score.to_numpy())
                c=np.column_stack(cov);df=len(ids)-np.linalg.matrix_rank(c)-1
                xr=x[ids]-c@np.linalg.lstsq(c,x[ids],rcond=None)[0]
                yr=y[ids]-c@np.linalg.lstsq(c,y[ids],rcond=None)[0]
                for _,target in candidates.iterrows():
                    j,k=pi[target.peak],gi[target.gene]
                    den=np.linalg.norm(xr[:,j])*np.linalg.norm(yr[:,k])
                    r=np.dot(xr[:,j],yr[:,k])/den if den>1e-10 else np.nan
                    p=2*stats.t.sf(abs(r)*np.sqrt(df/max(1-r*r,1e-12)),df) if np.isfinite(r) else np.nan
                    rows.append(dict(peak=target.peak,gene=target.gene,population=population,model='depth_COQ_state' if state else 'depth_COQ',gsm=gsm,n=len(ids),df=df,r=r,p=p,peak_detected=int(x[ids,j].sum()),RNA_detected=int((raw[ids,k]>0).sum())))
        print(gsm,len(meta),'eligible',len(pr),'pairs',flush=True)
    counts=pd.DataFrame(counts)
    for p,d in counts.groupby('peak'):
        z=peaks.set_index('peak').loc[p]
        assert d.high.sum()==z.high and d.low.sum()==z.low
    save(counts,out/'selected_ATAC_by_source.tsv');save(pd.DataFrame(qcs),out/'QC_flow.tsv')
    save(pd.DataFrame(rna_sources),out/'RNA_by_source.tsv')
    per=pd.DataFrame(rows);save(per,out/'links_by_source.tsv')
    combined=[]
    for (pop,model,peak,gene),d in per.groupby(['population','model','peak','gene']):
        rec=dict(population=pop,model=model,peak=peak,gene=gene,n=int(d.n.sum()),RNA_detected=int(d.RNA_detected.sum()),positive_sources=int((d.r>0).sum()))
        ok=d.r.notna()&(d.RNA_detected>=10)&(d.peak_detected>=10)
        rec['estimable_sources']=int(ok.sum())
        if ok.all():
            w=d.df-1;z=np.average(np.arctanh(np.clip(d.r,-.999999,.999999)),weights=w)
            rec.update(r=np.tanh(z),p=2*stats.norm.sf(abs(z)*np.sqrt(w.sum())))
        combined.append(rec)
    links=pd.DataFrame(combined);links['q_all_local_links']=np.nan
    for _,d in links.groupby(['population','model']):
        links.loc[d.index,'q_all_local_links']=multipletests(d.p.fillna(1),method='fdr_bh')[1]
    links.loc[links.p.isna(),'q_all_local_links']=np.nan
    save(links,out/'links_combined.tsv')
    h,l=np.vstack(hh),np.vstack(ll);rna=[]
    for j,g in enumerate(genes):
        delta=np.log1p(h[:,j])-np.log1p(l[:,j]);p=stats.ttest_1samp(delta,0).pvalue if np.std(delta)>0 else np.nan
        rna.append(dict(gene=g,n_pairs=len(h),high_mean_CP10k=h[:,j].mean(),low_mean_CP10k=l[:,j].mean(),FC=h[:,j].mean()/l[:,j].mean() if l[:,j].mean()>0 else np.nan,p=p))
    rna=pd.DataFrame(rna);rna['q_local_RNAs']=multipletests(rna.p.fillna(1),method='fdr_bh')[1]
    rna.loc[rna.p.isna(),'q_local_RNAs']=np.nan;save(rna,out/'RNA_high_low.tsv')
    base=links.query("population=='eligible_D0' and model=='depth_COQ'").merge(candidates,on=['peak','gene']).merge(rna,on='gene',suffixes=('_link','_RNA'))
    state=links.query("population=='eligible_D0' and model=='depth_COQ_state'")[['peak','gene','r','p','q_all_local_links']].rename(columns={'r':'state_r','p':'state_p','q_all_local_links':'state_q'})
    base=base.merge(state,on=['peak','gene']);save(base,out/'target_evidence.tsv')
    (out/'validation.json').write_text(json.dumps(dict(config='D206',pairs=527,peak_tests=565,regions_per_direction=25,selected_regions=len(hits),unique_selected_peaks=len(selected),local_RNA_pairs=len(candidates),local_RNAs=len(genes),raw_ATAC_counts_reproduced=True,BH_reproduced=True),indent=2)+'\n')
    print(base.sort_values('p_link')[['peak','gene','r','p_link','q_all_local_links','state_r','state_q','FC','p_RNA']].head(15).to_string(index=False))


if __name__=='__main__':main()

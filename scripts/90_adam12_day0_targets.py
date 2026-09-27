"""Local RNA target assessment for the D232 Day-0 accessibility candidate."""
import argparse
import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--work', type=Path, required=True)
    ap.add_argument('--h5', type=Path, required=True)
    a = ap.parse_args()
    root = Path(__file__).resolve().parents[1]
    out = root / 'results/adam12_day0_targets'
    out.mkdir(exist_ok=True)
    peak = 'chr10:126373783-126374503'
    candidates = pd.read_csv(root/'results/unstratified/candidate_pairs.tsv.gz', sep='\t')
    candidates = candidates[candidates.peak.eq(peak)].copy()
    genes = candidates.gene.tolist()
    candidates.to_csv(out/'candidate_RNAs.tsv', sep='\t', index=False)
    pairs = pd.read_csv(root/'results/day0_sensitivity_grid/matched_pairs.tsv.gz', sep='\t').query("config=='D232'")
    assert len(pairs) == 80
    spec = importlib.util.spec_from_file_location('extractor', root/'scripts/44_expand_middle_links.py')
    ext = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ext)
    mapping = pd.read_csv(root/'results/tables/consensus_peak_map.tsv.gz', sep='\t').set_index('peak')
    rows, contrasts, high_all, low_all = [], [], [], []
    opens = [0, 0]
    for gsm, pr in pairs.groupby('gsm'):
        meta = pd.read_csv(a.work/f'{gsm}_meta.tsv', sep='\t')
        meta = meta[meta.pass_tss3 & meta.myogenic_detected].reset_index(drop=True)
        local = peak if gsm == 'GSM6339597' else mapping.loc[peak, gsm+'_peak']
        mat, _ = ext.extract(next(a.h5.glob(gsm+'_*filtered_feature_bc_matrix.h5')), genes, [local], meta.barcode.tolist())
        raw = mat[:len(genes)].toarray().T.astype(float)
        x = (mat[len(genes)].toarray().ravel() > 0).astype(float)
        cp = raw / meta.total_rna_umi.to_numpy()[:,None] * 10000
        y = np.log1p(cp)
        ix = pd.Index(meta.barcode)
        hi, lo = ix.get_indexer(pr.high_barcode), ix.get_indexer(pr.low_barcode)
        assert min(hi.min(), lo.min()) >= 0
        opens[0] += int(x[hi].sum()); opens[1] += int(x[lo].sum())
        high_all.append(cp[hi]); low_all.append(cp[lo])
        for j,g in enumerate(genes):
            contrasts.append(dict(gsm=gsm,gene=g,n_pairs=len(pr),high_mean_CP10k=cp[hi,j].mean(),low_mean_CP10k=cp[lo,j].mean(),high_detected=int((raw[hi,j]>0).sum()),low_detected=int((raw[lo,j]>0).sum())))
        for population, ids in [('D232_eligible',np.arange(len(meta))),('D232_matched',np.r_[hi,lo])]:
            m = meta.iloc[ids]
            for state in [False,True]:
                cov = [np.ones(len(ids)),np.log1p(m.total_rna_umi),np.log1p(m.total_open_peaks),np.log1p(m.coq_CP10k)]
                if state: cov.append(m.state_score.to_numpy())
                c = np.column_stack(cov)
                xr = x[ids]-c@np.linalg.lstsq(c,x[ids],rcond=None)[0]
                yr = y[ids]-c@np.linalg.lstsq(c,y[ids],rcond=None)[0]
                df = len(ids)-np.linalg.matrix_rank(c)-1
                for j,g in enumerate(genes):
                    den = np.linalg.norm(xr)*np.linalg.norm(yr[:,j])
                    r = np.dot(xr,yr[:,j])/den if den>1e-10 else np.nan
                    p = 2*stats.t.sf(abs(r)*np.sqrt(df/max(1-r*r,1e-12)),df) if np.isfinite(r) else np.nan
                    rows.append(dict(population=population,model='technical_coq_state' if state else 'technical_coq',gsm=gsm,gene=g,n=len(ids),df=df,r=r,p=p,gene_detected=int((raw[ids,j]>0).sum()),peak_detected=int(x[ids].sum())))
        print(gsm,len(meta),'eligible',len(pr),'pairs',flush=True)
    assert opens == [30,7], opens
    per = pd.DataFrame(rows)
    per.to_csv(out/'links_by_source.tsv',sep='\t',index=False)
    combined=[]
    for (pop,model,gene), d in per.groupby(['population','model','gene']):
        rec=dict(population=pop,model=model,gene=gene,n=int(d.n.sum()),gene_detected=int(d.gene_detected.sum()))
        # Both sources must have enough detected observations to estimate a link.
        ok=d.r.notna() & (d.gene_detected>=10) & (d.peak_detected>=10)
        rec['estimable_sources']=int(ok.sum())
        if ok.all():
            w=d.df-1
            z=np.average(np.arctanh(np.clip(d.r,-.999999,.999999)),weights=w)
            rec.update(r=np.tanh(z),p=2*stats.norm.sf(abs(z)*np.sqrt(w.sum())))
        combined.append(rec)
    result=pd.DataFrame(combined)
    result['q_local6']=np.nan
    for _,d in result.groupby(['population','model']):
        result.loc[d.index,'q_local6']=multipletests(d.p.fillna(1),method='fdr_bh')[1]
        result.loc[d.index[d.p.isna()],'q_local6']=np.nan
    result.to_csv(out/'links_combined.tsv',sep='\t',index=False,na_rep='NA')
    h,l=np.vstack(high_all),np.vstack(low_all)
    rna=[]
    for j,g in enumerate(genes):
        delta=np.log1p(h[:,j])-np.log1p(l[:,j])
        p=stats.ttest_1samp(delta,0).pvalue if np.std(delta)>0 else np.nan
        rna.append(dict(gene=g,n_pairs=80,high_mean_CP10k=h[:,j].mean(),low_mean_CP10k=l[:,j].mean(),FC=h[:,j].mean()/l[:,j].mean() if l[:,j].mean()>0 else np.nan,p=p))
    rna=pd.DataFrame(rna)
    rna['q_local6']=multipletests(rna.p.fillna(1),method='fdr_bh')[1]
    rna.loc[rna.p.isna(),'q_local6']=np.nan
    rna.to_csv(out/'RNA_high_low.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(contrasts).to_csv(out/'RNA_high_low_by_source.tsv',sep='\t',index=False)
    print(result.to_string(index=False));print(rna.to_string(index=False))


if __name__=='__main__':
    main()

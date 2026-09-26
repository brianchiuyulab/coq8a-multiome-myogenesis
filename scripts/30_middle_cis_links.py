"""Explore cis RNA links for all 34 middle peaks using unchanged matched nuclei.

Candidates include all protein-coding genes with a transcript TSS within 500 kb
of the peak midpoint. Within-library partial correlations adjust COQ8A group and
RNA/ATAC depth, then use Fisher-z meta-analysis as in the existing link screen.
This is not Signac LinkPeaks: no GC/background-peak null is estimated here.
"""

import argparse
import gzip
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from multiome_core import h5_for, read_barcodes


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['tables','temporal','gtf','h5-root']:
        p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    members=pd.read_csv(a.temporal/'temporal_set_memberships.tsv.gz',sep='\t')
    peaks=sorted(members.loc[members.set_id.eq('P2_O50_Middle_24to48h'),'peak'])
    assert len(peaks)==34
    positions={s:(s.split(':')[0],np.mean(list(map(int,s.split(':')[1].split('-'))))) for s in peaks}
    chromosomes={x[0] for x in positions.values()}
    rows=[]
    with gzip.open(a.gtf,'rt') as f:
        for line in f:
            if line.startswith('#'):continue
            q=line.rstrip().split('\t')
            if q[2]!='transcript' or q[0] not in chromosomes or 'gene_type "protein_coding"' not in q[8]:continue
            name=re.search(r'gene_name "([^"]+)"',q[8]).group(1)
            tss=int(q[3] if q[6]=='+' else q[4])-1
            for peak,(chrom,center) in positions.items():
                distance=abs(center-tss)
                if q[0]==chrom and distance<=500000:
                    rows.append((peak,name,distance))
    candidates=pd.DataFrame(rows,columns=['peak','gene','tss_distance_bp']).groupby(['peak','gene'],as_index=False).tss_distance_bp.min()
    candidates.to_csv(a.temporal/'middle_cis_candidates_500kb.tsv',sep='\t',index=False)
    pairs=pd.read_csv(a.tables/'matched_pairs.tsv.gz',sep='\t')
    qc=pd.read_csv(a.tables/'qc_nuclei.tsv.gz',sep='\t')
    mapping=pd.read_csv(a.tables/'consensus_peak_map.tsv.gz',sep='\t').set_index('peak')
    output=[]
    for gsm,ps in pairs.groupby('gsm'):
        barcodes=sorted(set(ps.high_barcode)|set(ps.low_barcode))
        with h5py.File(h5_for(a.h5_root,gsm)) as h5:
            R,P,names,localpeaks=read_barcodes(h5,barcodes)
        ri={g:i for i,g in enumerate(names)};pi={str(g):i for i,g in enumerate(localpeaks)}
        present=candidates[candidates.gene.isin(ri)].copy()
        genes=sorted(present.gene.unique());gi={g:i for i,g in enumerate(genes)}
        peakids=[pi[s if gsm=='GSM6339597' else mapping.loc[s,gsm+'_peak']] for s in peaks]
        X=P[:,peakids].toarray().astype(float)
        raw=R[:,[ri[g] for g in genes]].toarray().astype(float)
        bmap={b:i for i,b in enumerate(barcodes)}
        quality=qc[qc.gsm.eq(gsm)].set_index('barcode').loc[barcodes]
        peaklookup={s:i for i,s in enumerate(peaks)}
        pix=present.peak.map(peaklookup).to_numpy();gix=present.gene.map(gi).to_numpy()
        for (gate,contrast),s in ps.groupby(['gate','contrast']):
            ids=[bmap[b] for b in s.high_barcode]+[bmap[b] for b in s.low_barcode]
            depth=quality.iloc[ids].total_rna_umi.to_numpy();atac=quality.iloc[ids].total_open_peaks.to_numpy()
            xx=X[ids];rr=raw[ids];yy=np.log1p(10000*rr/depth[:,None])
            cov=np.column_stack([np.ones(len(ids)),np.r_[np.ones(len(s)),np.zeros(len(s))],np.log1p(depth),np.log1p(atac)])
            cov[:,2:]=(cov[:,2:]-cov[:,2:].mean(axis=0))/cov[:,2:].std(axis=0)
            xr=xx-cov@np.linalg.lstsq(cov,xx,rcond=None)[0]
            yr=yy-cov@np.linalg.lstsq(cov,yy,rcond=None)[0]
            den=np.linalg.norm(xr[:,pix],axis=0)*np.linalg.norm(yr[:,gix],axis=0)
            r=np.divide(np.einsum('ij,ij->j',xr[:,pix],yr[:,gix]),den,out=np.full(len(present),np.nan),where=den>1e-12)
            out=present.copy();out['gsm']=gsm;out['gate']=gate;out['contrast']=contrast;out['n_nuclei']=len(ids)
            out['peak_open_nuclei']=xx[:,pix].sum(axis=0).astype(int)
            out['gene_detected_nuclei']=(rr[:,gix]>0).sum(axis=0)
            out['partial_r_descriptive']=r
            out['eligible']=(out.peak_open_nuclei>=10)&(out.gene_detected_nuclei>=20)&np.isfinite(r)
            out['partial_r']=out.partial_r_descriptive.where(out.eligible)
            output.append(out)
        print(gsm,'candidate links',len(present),flush=True)
    lib=pd.concat(output,ignore_index=True)
    lib.to_csv(a.temporal/'middle_cis_links_by_library.tsv.gz',sep='\t',index=False,na_rep='NA',compression={'method':'gzip','mtime':0})
    rows=[]
    for (gate,contrast,peak,gene),s in lib.groupby(['gate','contrast','peak','gene']):
        valid=s[s.eligible]
        record=dict(gate=gate,contrast=contrast,peak=peak,gene=gene,tss_distance_bp=s.tss_distance_bp.iloc[0],n_eligible_libraries=len(valid),n_gene_detected=int(s.gene_detected_nuclei.sum()),n_peak_open=int(s.peak_open_nuclei.sum()))
        if len(valid)>=3:
            z=np.arctanh(np.clip(valid.partial_r.to_numpy(),-.999999,.999999));w=valid.n_nuclei.to_numpy()-5
            mz=np.average(z,weights=w);se=1/np.sqrt(w.sum())
            record.update(partial_r=np.tanh(mz),r_ci_low=np.tanh(mz-1.96*se),r_ci_high=np.tanh(mz+1.96*se),p_link=2*stats.norm.sf(abs(mz/se)),positive_libraries=int((valid.partial_r>0).sum()))
        rows.append(record)
    pooled=pd.DataFrame(rows)
    for _,s in pooled.groupby(['gate','contrast']):
        valid=s.p_link.notna()
        pooled.loc[s.index[valid],'q_all_testable_middle_cis_links']=multipletests(s.loc[valid,'p_link'],method='fdr_bh')[1]
    pooled.sort_values(['gate','contrast','p_link','peak','gene']).to_csv(a.temporal/'middle_cis_links.tsv',sep='\t',index=False,na_rep='NA')
    # Regression check: shared peak-gene correlations must reproduce the prior screen.
    check=pooled[(pooled.gate=='TSS_ge_3')&(pooled.contrast=='3plus_vs_1')].merge(pd.read_csv(a.tables/'peak_gene_links_3plus_vs_1.tsv',sep='\t'),on=['peak','gene'],suffixes=('','_old'))
    assert np.allclose(check.partial_r,check.partial_r_old,atol=2e-7)
    for key,s in pooled.groupby(['gate','contrast']):
        print(key,'testable',int(s.p_link.notna().sum()),'positive q<.05',int(((s.partial_r>0)&(s.q_all_testable_middle_cis_links<.05)).sum()),flush=True)
        print(s.sort_values('p_link').head(8)[['peak','gene','partial_r','p_link','q_all_testable_middle_cis_links']].to_string(index=False),flush=True)


if __name__=='__main__':main()

"""Extract actual fragment-end profiles for the fixed D206 nuclei (Linux pysam)."""
import argparse,csv,gzip,importlib.util
from pathlib import Path
import pysam

ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--fragments',type=Path,required=True);ap.add_argument('--gtf',type=Path,required=True);a=ap.parse_args()
out=a.root/'results/fixed_day0_D206'
def read(p):
    with open(p) as f:return list(csv.DictReader(f,delimiter='\t'))
pairs=read(out/'matched_pairs.tsv')
spec=importlib.util.spec_from_file_location('tss',a.root/'scripts/14_tss_fragment_profiles.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
genes={r['gene'] for r in read(a.root/'reference/myogenesis_221_gene_sources.tsv')};tss=mod.gene_tss(a.gtf,genes)
windows=[]
for gene,(chrom,pos,strand) in tss.items():windows.append(dict(id=gene,kind='TSS',chrom=chrom,center=pos,start=pos-5000,end=pos+5000,strand=strand,bin_bp=100))
for r in read(out/'all_peaks.tsv'):
    chrom,coord=r['peak'].split(':');lo,hi=map(int,coord.split('-'));pos=(lo+hi)//2
    windows.append(dict(id=r['peak'],kind='peak',chrom=chrom,center=pos,start=pos-2000,end=pos+2000,strand='+',bin_bp=50))
for g,pk in [('MYOD1','chr11:17653349-17654252'),('CAV3','chr3:8757835-8758737'),('CSRP3','chr11:19218592-19219518')]:
    chrom,coord=pk.split(':');lo,hi=map(int,coord.split('-'));pos=tss[g][1]
    windows.append(dict(id=g,kind='locus',chrom=chrom,center=pos,start=min(lo,pos)-5000,end=max(hi,pos)+5000,strand='+',bin_bp=100))
with open(out/'profile_windows.tsv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(windows[0]),delimiter='\t');w.writeheader();w.writerows(windows)
fields=['id','kind','gsm','group','bin','relative_bp','insertions','per100_nuclei_per100bp']
with gzip.open(out/'fragment_profiles.tsv.gz','wt',newline='') as f:
    wr=csv.DictWriter(f,fieldnames=fields,delimiter='\t');wr.writeheader()
    for gsm,fgsm in [('GSM6339597','GSM6339598'),('GSM6339601','GSM6339602')]:
        ps=[p for p in pairs if p['gsm']==gsm];groups={p[g+'_barcode']:g for p in ps for g in ['high','low']}
        path=next(a.fragments.glob(fgsm+'*fragments.tsv.gz'))
        with pysam.TabixFile(str(path)) as tb:
            for k,win in enumerate(windows):
                start,end,bp=win['start'],win['end'],win['bin_bp'];n=(end-start+bp-1)//bp;c={g:[0]*n for g in ['high','low']}
                for line in tb.fetch(win['chrom'],max(0,start),end):
                    z=line.split('\t');g=groups.get(z[3])
                    if g is None:continue
                    for cut in [int(z[1]),int(z[2])]:
                        if start<=cut<end:c[g][(cut-start)//bp]+=1
                for g in ['low','high']:
                    vals=c[g][::-1] if win['strand']=='-' else c[g]
                    for j,v in enumerate(vals):
                        original=n-1-j if win['strand']=='-' else j;width=min(bp,end-start-original*bp)
                        wr.writerow(dict(id=win['id'],kind=win['kind'],gsm=gsm,group=g,bin=j,relative_bp=start+j*bp-win['center'],insertions=v,per100_nuclei_per100bp=v/len(ps)*100*100/width))
                if k%200==0:print(gsm,k,'/',len(windows),flush=True)
print('Completed',len(windows),'windows in both sources')

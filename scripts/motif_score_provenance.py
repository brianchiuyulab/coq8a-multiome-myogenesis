"""Score JASPAR 2024 MRF/MEF2 matrices on 500-bp peak-centred hg38 windows.

Outputs normalized best PWM score (0..1) for either strand, per consensus peak.
Numerical score cutoffs are examined downstream as sensitivity settings.
"""
from pathlib import Path
import argparse
import json,time
import numpy as np
import pandas as pd
from twobitreader import TwoBitFile
from numba import njit,prange,set_num_threads

NAMES=['MYOD1','MYOG','MYF5','MYF6','MEF2A','MEF2B','MEF2C','MEF2D']
W=500

@njit(parallel=True)
def max_scores(seq,pwm,rc,lengths,lo,hi):
 n=seq.shape[0];k=pwm.shape[2];out=np.zeros((n,k),dtype=np.float32)
 for i in prange(n):
  for m in range(k):
   l=lengths[m]
   best=-1e8
   for start in range(seq.shape[1]-l+1):
    s=0.;t=0.;ok=True
    for j in range(l):
     b=seq[i,start+j]
     if b>3:
      ok=False;break
     s+=pwm[j,b,m];t+=rc[j,b,m]
    if ok:
     if s>best:best=s
     if t>best:best=t
   out[i,m]=(best-lo[m])/(hi[m]-lo[m]) if best>-1e7 else 0.
 return out

def main():
 p=argparse.ArgumentParser()
 p.add_argument('--consensus',type=Path,required=True)
 p.add_argument('--matrices',type=Path,required=True)
 p.add_argument('--genome-2bit',type=Path,required=True)
 p.add_argument('--out',type=Path,required=True)
 a=p.parse_args()
 d=pd.read_csv(a.consensus,sep='\t')
 if 'n_matched_samples' not in d: d['n_matched_samples']=d[[c for c in d if c.endswith('_peak')]].notna().sum(axis=1)+1
 d=d[(d.n_matched_samples>=3)&d.peak.str.match(r'^chr(?:[1-9]|1[0-9]|2[0-2]|X|Y):')].copy().reset_index(drop=True)
 coord=d.peak.str.extract(r'^(chr[^:]+):(\d+)-(\d+)$')
 chrom=coord[0].to_numpy();start=coord[1].astype(int).to_numpy();end=coord[2].astype(int).to_numpy()
 mid=(start+end)//2
 j=json.loads(a.matrices.read_text())
 lengths=np.array([len(j[n]['pfm']['A']) for n in NAMES],dtype=np.int64)
 L=int(lengths.max());pwm=np.zeros((L,4,8),dtype=np.float32)
 for i,name in enumerate(NAMES):
  pfm=np.stack([j[name]['pfm'][b] for b in 'ACGT'],axis=1).astype(float)
  p=(pfm+.25)/(pfm.sum(axis=1,keepdims=True)+1.)
  pwm[:lengths[i],:,i]=np.log2(p/.25)
 rc=np.zeros_like(pwm)
 for i in range(8):rc[:lengths[i],:,i]=pwm[:lengths[i],:,i][::-1,::-1]
 lows=np.zeros(8,dtype=np.float32);highs=np.zeros(8,dtype=np.float32)
 for i in range(8):
  lows[i]=np.min(pwm[:lengths[i],:,i],axis=1).sum()
  highs[i]=np.max(pwm[:lengths[i],:,i],axis=1).sum()
 tr=bytearray([4]*256)
 for c,v in [('A',0),('C',1),('G',2),('T',3),('a',0),('c',1),('g',2),('t',3)]:tr[ord(c)]=v
 seq=np.full((len(d),W),4,dtype=np.uint8)
 genome=TwoBitFile(str(a.genome_2bit))
 t=time.monotonic()
 for c in sorted(set(chrom)):
  ix=np.flatnonzero(chrom==c)
  if c not in genome:continue
  reference=genome[c]
  for i in ix:
   left=max(0,int(mid[i]-W//2));right=min(len(reference),left+W)
   s=reference[left:right]
   arr=np.frombuffer(s.encode('ascii').translate(tr),dtype=np.uint8)
   seq[i,:len(arr)]=arr
  print('sequence',c,len(ix),'seconds',round(time.monotonic()-t,1),flush=True)
 gc=np.mean((seq==1)|(seq==2),axis=1).astype(np.float32)
 print('sequence matrix',seq.shape,'GC mean',gc.mean(),flush=True)
 set_num_threads(6)
 scores=max_scores(seq[:10],pwm,rc,lengths,lows,highs)
 t=time.monotonic();scores=max_scores(seq,pwm,rc,lengths,lows,highs)
 print('motif scan seconds',round(time.monotonic()-t,1),flush=True)
 out=pd.DataFrame({'peak':d.peak,'gc_500bp':gc,'width':end-start,'n_matched_samples':d.n_matched_samples})
 for i,name in enumerate(NAMES):out[name+'_score']=scores[:,i]
 a.out.parent.mkdir(parents=True,exist_ok=True)
 out.to_csv(a.out,sep='\t',index=False,compression='gzip')
 print(out[[x+'_score' for x in NAMES]].describe().to_string(),flush=True)

if __name__=='__main__':main()

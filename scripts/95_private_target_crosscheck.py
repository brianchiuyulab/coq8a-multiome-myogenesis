"""Read-only C2C12 expression cross-reference; outputs must remain private."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import openpyxl

ap=argparse.ArgumentParser();ap.add_argument('--counts',type=Path,required=True);ap.add_argument('--de',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
root=Path(__file__).resolve().parents[1]
assert root not in a.out.resolve().parents and a.out.resolve()!=root
a.out.mkdir(parents=True,exist_ok=True)
scope=pd.read_csv(root/'results/fixed_day0_D206/local_RNA_candidates.tsv',sep='\t');wanted=set(scope.gene)
w=openpyxl.load_workbook(a.counts,read_only=True,data_only=True);sheet=w.worksheets[0];rows=sheet.iter_rows(values_only=True);header=next(rows);found=[];wide=[]
for ix,row in enumerate(rows,2):
    symbol=str(row[1]);human=symbol.upper()
    if human not in wanted:continue
    v=np.asarray(row[4:],float)
    found.append(dict(human_symbol=human,mouse_symbol=symbol,mouse_ensembl=row[0],sheet=sheet.title,row=ix,n_samples=len(v),detected=int((v>0).sum()),at_least10=int((v>=10).sum()),minimum=v.min(),maximum=v.max(),mapping='case-normalized symbol match; not a full orthology inference'))
    for sample,value in zip(header[4:],v):wide.append(dict(human_symbol=human,mouse_symbol=symbol,sample=sample,normalized_count=value))
w.close()
det=pd.DataFrame(found);det.to_csv(a.out/'detection.tsv',sep='\t',index=False)
pd.DataFrame(wide).to_csv(a.out/'normalized_counts.tsv',sep='\t',index=False)
pd.DataFrame({'human_symbol':sorted(wanted-set(det.human_symbol)),'status':'No same-symbol match; not proof of absent expression or absent mouse ortholog'}).to_csv(a.out/'unresolved_symbols.tsv',sep='\t',index=False)
de=pd.read_csv(a.de);de['human_symbol']=de.symbol.str.upper();de=de[de.human_symbol.isin(wanted)].copy();de['FC']=2**de.log2FC
de['p_rounded_to_zero_in_source']=de.p_value.eq(0);de['fdr_rounded_to_zero_in_source']=de.fdr.eq(0)
de.to_csv(a.out/'all_contrasts.tsv',sep='\t',index=False)
sub=de[de.contrast.str.match(r'P(11|22|33)_D[12]_vs_D0$')];sub.to_csv(a.out/'D1_D2_vs_D0.tsv',sep='\t',index=False)
focus=['MYOD1','CAV3','CSRP3','SOX8','LMF1','NOTCH1','BOC','CCDC80','MYH9','CACNA1H','MYF5','MYF6','NOS1']
cols=[f'P{p}_D{day}_vs_D0' for p in [11,22,33] for day in [1,2]]
v=sub.pivot(index='human_symbol',columns='contrast',values='log2FC').reindex(index=focus,columns=cols)
q=sub.pivot(index='human_symbol',columns='contrast',values='fdr').reindex(index=focus,columns=cols)
print(det[det.human_symbol.isin(focus)][['human_symbol','row','detected','at_least10']].to_string(index=False))
print(sub[sub.human_symbol.isin(focus)&sub.contrast.str.contains('_D2_')][['symbol','contrast','FC','fdr']].to_string(index=False))

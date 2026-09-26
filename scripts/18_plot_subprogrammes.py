"""Show external muscle subprogramme screens and two selected exploratory routes."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd

SMARCA='DELASERNA_TARGETS_OF_MYOD_AND_SMARCA4'
POS='GOBP_POSITIVE_REGULATION_OF_MUSCLE_CELL_DIFFERENTIATION'
SETS=[
    (SMARCA,'MYOD / SMARCA4-dependent'),
    ('DELASERNA_MYOD_TARGETS_UP','MYOD-induced'),
    (POS,'Promote muscle differentiation'),
    ('GOBP_MUSCLE_CELL_FATE_COMMITMENT','Muscle fate commitment'),
    ('GOBP_MYOBLAST_DIFFERENTIATION','Myoblast differentiation'),
    ('GOBP_MYOTUBE_DIFFERENTIATION','Myotube differentiation'),
    ('GOBP_MYOBLAST_FUSION','Myoblast fusion'),
    ('GOBP_MYOFIBRIL_ASSEMBLY','Myofibril assembly'),
    ('GOBP_MUSCLE_CONTRACTION','Muscle contraction'),
    ('GOBP_MUSCLE_CELL_DIFFERENTIATION','Muscle differentiation'),
    ('MEF2C_TARGET_GENES','MEF2C target set'),
    ('MEF2D_TARGET_GENES','MEF2D target set'),
]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exploration',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    d=pd.read_csv(a.exploration/'effect_grid.tsv',sep='\t')
    m=d[(d.scope=='programme')&(d.gate=='TSS_ge_3')&(d.contrast=='3plus_vs_1')]
    sources=pd.read_csv(a.exploration/'subprogramme_sources.tsv',sep='\t').set_index('programme')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
    fig=plt.figure(figsize=(14,10));gs=fig.add_gridspec(2,2,width_ratios=[1.1,1],left=.26,right=.96,top=.88,bottom=.16,wspace=.7,hspace=.62)
    fig.suptitle('Narrowing the 221-gene search with external muscle programmes',x=.05,y=.98,ha='left',fontsize=16,fontweight='bold')
    fig.text(.05,.94,'COQ8A >=3 vs 1 UMI | TSS enrichment >=3 | 201 matched pairs | 4 libraries from 2 source lines',fontsize=10)
    ax=fig.add_subplot(gs[:,0]);matrix=np.full((len(SETS),3),np.nan);texts={};source_rows=[]
    for i,(name,label) in enumerate(SETS):
        for j,rule in enumerate(['all_common','promoter_2kb','linked_MRFscore_lt095']):
            sub=m[(m.name==name)&(m.rule==rule)]
            if j==2:
                sub=sub[(sub.link_gate=='TSS_ge_3')&(sub.link_contrast=='2plus_vs_1')]
            if sub.empty:continue
            row=sub.iloc[0];matrix[i,j]=row.fold_open
            texts[i,j]=f'{row.fold_open:.2f}'+('*' if row.p_pair<.05 else '')+f'\n{int(row.n_genes)}g / {int(row.n_peaks)}p'
            source_rows.append(row)
    cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#EEEEEE')
    image=ax.imshow(matrix,aspect='auto',cmap=cmap,norm=TwoSlopeNorm(vmin=.75,vcenter=1,vmax=1.25))
    for (i,j),label in texts.items():
        ax.text(j,i,label,ha='center',va='center',fontsize=9,color='white' if matrix[i,j]>1.19 or matrix[i,j]<.81 else '#202020')
    for i,j in zip(*np.where(np.isnan(matrix))):ax.text(j,i,'No set',ha='center',va='center',fontsize=8,color='#666666')
    ax.set_yticks(range(len(SETS)),[f'{label}\n({int(sources.loc[name,"n_genes_in_221"])} of 221 genes)' for name,label in SETS],fontsize=9)
    ax.set_xticks(range(3),['All nearby\ncommon peaks','Promoter\n<=2 kb','RNA-linked\nMRF score <0.95'],fontsize=9)
    ax.set_title('A  Externally defined subprogrammes',loc='left',fontsize=12,fontweight='bold',pad=15)
    cb=fig.colorbar(image,ax=ax,orientation='horizontal',fraction=.035,pad=.11);cb.set_label('High / low open-fraction ratio',fontsize=9)
    ax=fig.add_subplot(gs[0,1]);r=m[(m.name==SMARCA)&m.rule.str.startswith('promoter')].set_index('rule').loc[['promoter_05kb','promoter_1kb','promoter_2kb','promoter_5kb']]
    x=np.arange(4)
    ax.plot(x,r.fold_open,'o-',color='#16877C',lw=2,markersize=6)
    ax.axhline(1,color='#AEB8C4',lw=1);ax.set_ylim(.99,1.20);ax.set_xlim(-.3,3.3)
    for i,row in enumerate(r.itertuples()):
        ax.text(i,row.fold_open+.012,f'{row.fold_open:.3f}x',ha='center',fontsize=9)
        ax.text(i,1.012,f'p={row.p_pair:.4f}',ha='center',fontsize=8,rotation=45)
    ax.set_xticks(x,['0.5','1','2','5']);ax.set_xlabel('Maximum distance to transcript TSS (kb)');ax.set_ylabel('Open-fraction ratio')
    ax.set_title('B  MYOD / SMARCA4 programme',loc='left',fontsize=12,fontweight='bold',pad=15)
    ax.text(.02,.93,'4/4 positive libraries at every window',transform=ax.transAxes,fontsize=9)
    ax.spines[['top','right']].set_visible(False)
    ax=fig.add_subplot(gs[1,1]);sel=m[(m.name==POS)&(m.rule=='linked_MRFscore_lt095')&(m.link_gate=='TSS_ge_3')&(m.link_contrast=='2plus_vs_1')].iloc[0]
    assert set(sel.represented_genes.split(';'))=={'MYOD1','CAV3','LAMA2'}
    lib=pd.read_csv(a.exploration/'effects_by_library.tsv.gz',sep='\t')
    lib=lib[(lib.set_id==sel.set_id)&(lib.gate=='TSS_ge_3')&(lib.contrast=='3plus_vs_1')].sort_values('gsm')
    for row,col in zip(lib.itertuples(),['#4169A1','#79A8CD','#915794','#C193B9']):
        ax.plot([0,1],[row.low_open_pct,row.high_open_pct],'-o',color=col,lw=1.6,markersize=5)
    ax.plot([0,1],[sel.low_open_pct,sel.high_open_pct],'-D',color='#222222',lw=3,markersize=6)
    ax.set_xticks([0,1],['COQ8A low','COQ8A high']);ax.set_xlim(-.15,1.15);ax.set_ylim(0,26)
    ax.set_ylabel('Seven-peak open fraction (%)')
    ax.set_title('C  Pro-differentiation linked subset',loc='left',fontsize=12,fontweight='bold',pad=15)
    npos=int((lib.delta_pp>1e-10).sum());nflat=int((lib.delta_pp.abs()<=1e-10).sum())
    ax.text(.03,.95,f'1.230x | p = {sel.p_pair:.4f}\nMYOD1, CAV3, LAMA2\n{npos} libraries higher; {nflat} unchanged',transform=ax.transAxes,va='top',fontsize=10)
    ax.spines[['top','right']].set_visible(False)
    fig.text(.05,.06,'A: ratios and effective genes (g) / peaks (p); * nominal paired p<0.05. Grey cells: no qualifying peak set.\nRNA links: COQ8A >=2 vs 1, TSS>=3, partial r>0, link q<0.05 and four evaluable libraries.\nC: 221 starting genes -> 14 GO pro-differentiation genes -> 6 genes with 18 positive links -> 3 genes with 7 lower-MRF-score links.\nSIN3A signature overlaps only one of the 221 genes; its results remain in the complete tables.',fontsize=9,color='#555555',linespacing=1.5)
    pd.DataFrame(source_rows).to_csv(a.exploration/'subprogramme_figure_source.tsv',sep='\t',index=False)
    for ext in ['png','pdf']:fig.savefig(a.out/f'Figure_subprogrammes.{ext}',dpi=190,facecolor='white')
    plt.close(fig)
    print(a.out/'Figure_subprogrammes.png')


if __name__=='__main__':main()

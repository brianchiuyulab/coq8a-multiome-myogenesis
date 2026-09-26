"""Plot the exploratory route from RNA programmes to ranked ATAC loci."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


RED, BLUE, GREY, TEAL = '#B83745', '#3477A9', '#BAC2CB', '#16877C'


def style(ax, title):
    ax.set_title(title, loc='left', fontsize=12, fontweight='bold', pad=15)
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(labelsize=10)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--tables',type=Path,required=True)
    p.add_argument('--exploration',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args(); a.out.mkdir(parents=True,exist_ok=True)
    grid=pd.read_csv(a.exploration/'effect_grid.tsv',sep='\t')
    ranks=pd.read_csv(a.exploration/'candidate_ranking_grid.tsv',sep='\t')
    main=grid[(grid.gate=='TSS_ge_3')&(grid.contrast=='3plus_vs_1')]
    focal=main[(main.name=='MYOD1')&(main.link_gate=='TSS_ge_2')&(main.link_contrast=='2plus_vs_1')&(main.rule=='linked_MRFscore_lt095')].iloc[0]
    assert np.isclose(focal.fold_open,1.28,atol=1e-5) and focal.n_pairs==201
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'ps.fonttype':42})
    fig,axes=plt.subplots(3,2,figsize=(12.5,13.2))
    fig.subplots_adjust(left=.12,right=.96,top=.91,bottom=.10,hspace=.66,wspace=.46)
    fig.suptitle('COQ8A-associated myogenesis: from the 221-gene screen to MYOD1',fontsize=16,fontweight='bold',x=.12,ha='left',y=.975)
    fig.text(.12,.943,'Main effect comparison: COQ8A >=3 vs 1 UMI | TSS enrichment >=3 | 201 matched pairs',fontsize=10,color='#444444')
    ax=axes[0,0]
    programme=pd.read_csv(a.tables/'programme_effects_pooled.tsv',sep='\t')
    programme=programme[(programme.gate=='TSS_ge_3')&(programme.contrast=='3plus_vs_1')&(programme.modality=='RNA')].set_index('programme')
    ps=programme.loc[['Hallmark_myogenesis','Reactome_myogenesis','MRF_loci']]
    y=np.arange(3)
    ax.errorbar(ps.difference,y,xerr=[ps.difference-ps.ci_low,ps.ci_high-ps.difference],fmt='o',color=RED,capsize=4,markersize=7)
    for i,row in enumerate(ps.itertuples()):
        ax.text(.102,i,f'q = {row.q_4_programmes:.3f}',va='center',fontsize=10)
    ax.axvline(0,c=GREY,lw=1);ax.set_xlim(-.003,.16)
    ax.set_yticks(y,['Hallmark\nmyogenesis','Reactome\nmyogenesis','MRF genes']);ax.invert_yaxis()
    ax.set_xlabel('RNA programme difference (high - low)\nMean log1p(CP10K) per gene; 95% CI')
    style(ax,'A  RNA programme association')
    ax=axes[0,1]
    allgenes=main[(main.scope=='gene')&(main.rule=='all_common')].copy()
    ax.scatter(np.log2(allgenes.fold_open),-np.log10(allgenes.p_pair),s=15,c=GREY,alpha=.8)
    for name,col,offset in [('MYOD1',RED,(8,4)),('MYL1',TEAL,(-40,8))]:
        row=allgenes[allgenes.name==name].iloc[0]
        xy=(np.log2(row.fold_open),-np.log10(row.p_pair))
        ax.scatter(*xy,s=50,c=col,zorder=3);ax.annotate(name,xy,xytext=offset,textcoords='offset points',color=col,fontweight='bold')
    ax.axvline(0,c=GREY,lw=1)
    ax.set_xlabel('log2 accessibility ratio (high / low)')
    ax.set_ylabel('-log10 nominal paired p')
    style(ax,'B  All 221 ATAC gene regions')
    ax=axes[1,0]
    rank=ranks[(ranks.link_gate=='TSS_ge_2')&(ranks.link_contrast=='2plus_vs_1')&(ranks.gate=='TSS_ge_3')&(ranks.contrast=='3plus_vs_1')].head(6)
    assert rank.iloc[0].gene=='MYOD1'
    ax.barh(rank.gene,rank.n_concordant,color=[RED if g=='MYOD1' else GREY for g in rank.gene],height=.65)
    ax.invert_yaxis();ax.set_xticks(range(5));ax.set_xlim(0,4.5)
    ax.set_xlabel('Positive RNA links with increased ATAC\nin at least 3 of 4 libraries')
    style(ax,'C  Rank candidates across 221 genes')
    ax=axes[1,1]
    rows=[allgenes[allgenes.name=='MYOD1'].iloc[0],
          main[(main.name=='MYOD1')&(main.rule=='positive_link')&(main.link_gate=='TSS_ge_2')&(main.link_contrast=='2plus_vs_1')].iloc[0],focal]
    for i,row in enumerate(rows):
        ax.plot([row.low_open_pct,row.high_open_pct],[i,i],color=GREY,lw=2)
        ax.scatter(row.low_open_pct,i,c=BLUE,s=55,label='Low' if i==0 else None,zorder=3)
        ax.scatter(row.high_open_pct,i,c=RED,s=55,label='High' if i==0 else None,zorder=3)
        ax.text(18,i,f'{row.fold_open:.2f}x',va='center',fontweight='bold',color=RED if i==2 else '#333333')
    ax.set_yticks(range(3),['19 common peaks','11 positive RNA links','6 links with MRF\nscore <0.95'])
    ax.set_ylim(2.7,-.4);ax.set_xlim(0,23);ax.set_xlabel('Mean fraction of peaks open (%)')
    ax.legend(frameon=False,ncol=2,loc='lower right',fontsize=9)
    style(ax,'D  Focus the MYOD1 region')
    ax=axes[2,0]
    libs=pd.read_csv(a.exploration/'effects_by_library.tsv.gz',sep='\t')
    libs=libs[(libs.set_id==focal.set_id)&(libs.gate=='TSS_ge_3')&(libs.contrast=='3plus_vs_1')].sort_values('gsm')
    labels=['Line 1 stem','Line 1 diff.','Line 2 stem','Line 2 diff.']
    colors=['#4169A1','#79A8CD','#915794','#C193B9']
    for row,label,color in zip(libs.itertuples(),labels,colors):
        ax.plot([0,1],[row.low_open_pct,row.high_open_pct],'-o',color=color,lw=1.7,markersize=5,label=label)
        ax.text(1.06,row.high_open_pct,label,va='center',color=color,fontsize=9)
    ax.plot([0,1],[focal.low_open_pct,focal.high_open_pct],'-D',color='#222222',lw=3,markersize=6,label='Pooled')
    ax.text(1.06,focal.high_open_pct,'Pooled',va='center',fontsize=9,fontweight='bold')
    ax.set_xticks([0,1],['COQ8A low','COQ8A high']);ax.set_xlim(-.12,1.55);ax.set_ylim(0,24)
    ax.set_ylabel('Six-peak open fraction (%)')
    ax.text(.04,.95,f'1.28x | p = {focal.p_pair:.4f}\n4/4 library directions',transform=ax.transAxes,va='top',fontsize=10)
    style(ax,'E  The six-peak effect by library')
    ax=axes[2,1]
    s=grid[(grid.name=='MYOD1')&(grid.rule=='linked_MRFscore_lt095')&(grid.link_gate=='TSS_ge_2')&(grid.link_contrast=='2plus_vs_1')]
    order=[(g,c) for g in ['TSS_ge_2','TSS_ge_3'] for c in ['2plus_vs_1','3plus_vs_1']]
    s=s.set_index(['gate','contrast']).loc[order]
    matrix=s.fold_open.to_numpy().reshape(2,2)
    im=ax.imshow(matrix,cmap='Reds',vmin=1,vmax=1.30,aspect='auto')
    for i,row in enumerate(s.itertuples()):
        ax.text(i%2,i//2,f'{row.fold_open:.3f}x\np = {row.p_pair:.4f}\n{row.n_pairs} pairs',ha='center',va='center',color='white' if row.fold_open>1.22 else '#222222',fontsize=11)
    ax.set_xticks([0,1],['>=2 vs 1 UMI','>=3 vs 1 UMI']);ax.set_yticks([0,1],['TSS >=2','TSS >=3'])
    ax.set_xlabel('COQ8A effect contrast');ax.set_ylabel('ATAC QC threshold')
    style(ax,'F  Fixed six peaks: sensitivity')
    fig.text(.12,.04,'Links in C-F: COQ8A >=2 vs 1, TSS >=2; positive partial correlation, link q<0.05, four evaluable libraries.\nB, E and F show nominal paired-nucleus p values. Four libraries represent two source lines. Exploratory analysis.',fontsize=9,color='#555555',linespacing=1.5)
    for suffix in ['png','pdf']:
        fig.savefig(a.out/f'Figure_exploratory_route.{suffix}',dpi=190,facecolor='white')
    plt.close(fig)
    print(a.out/'Figure_exploratory_route.png')


if __name__=='__main__':
    main()

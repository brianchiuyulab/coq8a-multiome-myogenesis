"""Publication panels for the fixed D206 candidate screen and RNA follow-up."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

root=Path(__file__).resolve().parents[1]
src=root/'results/fixed_day0_D206';out=root/'figures/fixed_day0_D206';out.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':220})
blue='#326b9b';red='#bf5145';gray='#b8bdc3'
p=pd.read_csv(src/'all_peaks.tsv',sep='\t');r=pd.read_csv(src/'all_regions.tsv',sep='\t');hits=pd.read_csv(src/'selected_regions.tsv',sep='\t')
labels=hits.groupby('top_peak',sort=False).gene.apply(lambda z:'/'.join(z)).to_dict()
order=hits.sort_values(['mode','q25']).top_peak.drop_duplicates().tolist()
order=hits.sort_values('q25').top_peak.drop_duplicates().tolist()


def output(fig,name):
    fig.savefig(out/f'{name}.png',bbox_inches='tight',facecolor='white')
    fig.savefig(out/f'{name}.pdf',bbox_inches='tight',facecolor='white')
    plt.close(fig)


fig=plt.figure(figsize=(13,8),layout='constrained')
gs=fig.add_gridspec(2,2,height_ratios=[.35,1],width_ratios=[1,1.35])
ax=fig.add_subplot(gs[0,:]);ax.axis('off')
steps=['221 myogenesis genes','25 differentiation /\nfusion genes','565 candidate peaks','D0: 527 matched pairs\nCOQ8A >=2 vs 0','6 regions q < 0.05\n+ 4 regions q < 0.10']
for i,t in enumerate(steps):
    ax.text((i+.5)/5,.5,t,ha='center',va='center',fontsize=10,bbox=dict(boxstyle='round,pad=.65',facecolor='#edf2f6',edgecolor='#b8c7d4'),transform=ax.transAxes)
    if i<4:ax.annotate('',xy=((i+1.08)/5,.5),xytext=((i+.93)/5,.5),xycoords='axes fraction',arrowprops=dict(arrowstyle='->',color=gray))
ax.set_title('A   External scope and fixed Day-0 comparison',loc='left',fontweight='bold')
ax=fig.add_subplot(gs[1,0])
vals=np.column_stack([p.source1_low/p.source1_pairs,p.source1_high/p.source1_pairs,p.source2_low/p.source2_pairs,p.source2_high/p.source2_pairs])*100
idx=np.argsort(vals.mean(1))[::-1]
im=ax.imshow(vals[idx],aspect='auto',cmap='Blues',vmin=0,vmax=100,interpolation='nearest')
ax.set_xticks(range(4),['Source 1\nLow','Source 1\nHigh','Source 2\nLow','Source 2\nHigh'])
ax.set_ylabel('565 peaks, ordered by mean accessibility');ax.set_yticks([0,564],[1,565])
ax.set_title('B   Complete candidate-peak landscape',loc='left',fontweight='bold')
fig.colorbar(im,ax=ax,shrink=.65,label='Nuclei with accessible peak (%)')
ax=fig.add_subplot(gs[1,1]);xx=np.log2((p.high+.5)/(p.low+.5));yy=-np.log10(p.p)
ax.scatter(xx,yy,s=15,c=gray,alpha=.6,rasterized=True)
for direction,color,mask in [('Opening',red,(p.q565<.05)&(p.FC>1)),('Closing',blue,(p.q565<.05)&(p.FC<1))]:
    ax.scatter(xx[mask],yy[mask],s=34,c=color,label='Peak q < 0.05: '+direction)
for peak,offset in [('chr11:17653349-17654252',(8,15)),('chr11:19218592-19219518',(15,-6)),('chr16:1078248-1079160',(15,10)),('chr3:113117234-113118143',(-5,15)),('chr9:136586858-136587780',(-85,-18)),('chr3:8757835-8758737',(-15,-25))]:
    j=p.index[p.peak.eq(peak)][0];ax.annotate(labels[peak],(xx[j],yy[j]),xytext=offset,textcoords='offset points',fontsize=9,arrowprops=dict(arrowstyle='-',color='#666666',lw=.6))
ax.axvline(0,color='#777777',lw=.8);ax.set_xlabel('log2 accessibility ratio (0.5-count display offset)');ax.set_ylabel('-log10 paired peak p')
ax.set_title('C   All 565 peak tests',loc='left',fontweight='bold');ax.legend(loc='lower left',fontsize=8,frameon=False)
output(fig,'Figure_1_screen')

fig,axs=plt.subplots(3,3,figsize=(12,10),layout='constrained')
c=pd.read_csv(src/'selected_ATAC_by_source.tsv',sep='\t')
for ax,peak in zip(axs.flat,order):
    d=c[c.peak.eq(peak)].sort_values('gsm');z=p.set_index('peak').loc[peak];region=hits[hits.top_peak.eq(peak)]
    for i,(_,row) in enumerate(d.iterrows()):
        yy=np.array([row.low,row.high])/row.pairs*100
        ax.plot([0,1],yy,'o-',color=[blue,red][i],lw=1.8,label=f'Source {i+1} (n={row.pairs} pairs)')
    ax.set_xticks([0,1],['COQ8A low\n0 UMI','COQ8A high\n>=2 UMI']);ax.set_xlim(-.25,1.25)
    ax.set_ylim(bottom=0);ax.set_ylabel('Accessible nuclei (%)')
    ax.set_title(f'{labels[peak]} neighborhood\nFC {z.FC:.2f} | regional q {region.q25.min():.3g}\n{peak}',fontsize=9,pad=10)
    ax.legend(fontsize=7,frameon=False,loc='best')
fig.suptitle('Selected local accessibility signals | Same 527 pairs throughout',fontsize=15,fontweight='bold')
output(fig,'Figure_2_selected_loci')

e=pd.read_csv(src/'target_evidence.tsv',sep='\t')
shown=[]
for peak in order:
    d=e[e.peak.eq(peak)].copy()
    scope=set(labels[peak].split('/'))
    extra=d[d.r.notna()].assign(abs_r=lambda z:z.r.abs()).sort_values('abs_r',ascending=False).head(2).gene.tolist()
    wanted=scope|set(extra)
    shown.append(d[d.gene.isin(wanted)].sort_values('gene'))
d=pd.concat(shown,ignore_index=True)
savecols=['peak','gene','r','q_all_local_links','state_r','state_q','FC','p_RNA','q_local_RNAs']
d[savecols].to_csv(src/'figure3_display_rows.tsv',sep='\t',index=False,na_rep='NA')
fig,axs=plt.subplots(1,3,figsize=(13,max(7,len(d)*.34)),layout='constrained',gridspec_kw={'width_ratios':[1,1,.9]})
yt=[f'{labels[row.peak]} region  >  {row.gene}' for _,row in d.iterrows()]
for ax,cols,title in [(axs[0],['r','state_r'],'A   Peak-RNA partial correlation'),(axs[1],['q_all_local_links','state_q'],'B   Peak-RNA evidence')]:
    v=d[cols].to_numpy(float)
    if ax==axs[0]:im=ax.imshow(np.ma.masked_invalid(v),cmap='RdBu_r',vmin=-.16,vmax=.16,aspect='auto');fmt=lambda z:f'{z:.3f}'
    else:im=ax.imshow(np.ma.masked_invalid(-np.log10(np.maximum(v,1e-12))),cmap='Blues',vmin=0,vmax=6,aspect='auto');fmt=lambda z:f'{z:.1e}' if z<.001 else f'{z:.3f}'
    for i in range(len(d)):
        for j in range(2):ax.text(j,i,fmt(v[i,j]) if np.isfinite(v[i,j]) else 'NA',ha='center',va='center',fontsize=8,color='white' if (abs(v[i,j])>.09 if ax==axs[0] else v[i,j]<.001) else 'black')
    ax.set_xticks([0,1],['Depth + COQ8A','+ cell state'],rotation=20,ha='right');ax.set_yticks(range(len(d)),yt if ax==axs[0] else [])
    ax.set_title(title,loc='left',fontsize=11,fontweight='bold')
    fig.colorbar(im,ax=ax,shrink=.35,label='Partial r' if ax==axs[0] else '-log10 link q (color capped at 6)')
ax=axs[2];fc=d.FC.to_numpy();ys=np.arange(len(d));valid=np.isfinite(fc)&(fc>0)
ax.scatter(np.log2(fc[valid]),ys[valid],c=[red if v>1 else blue for v in fc[valid]],s=32)
for i in range(len(d)):
    if valid[i]:ax.annotate(f'{fc[i]:.2f}',(np.log2(fc[i]),i),xytext=(6,0),textcoords='offset points',va='center',fontsize=8)
ax.axvline(0,c=gray,lw=1);ax.set_ylim(len(d)-.5,-.5);ax.set_yticks([]);ax.set_xlabel('RNA log2 FC (high / low)');ax.set_title('C   Same-pair RNA response',loc='left',fontsize=11,fontweight='bold')
ax.margins(x=.4)
output(fig,'Figure_3_RNA_targets')

fig,ax=plt.subplots(figsize=(8,9),layout='constrained')
v=r.pivot(index='gene',columns='mode',values='q25')[['opening','closing']]
v=v.loc[v.min(1).sort_values().index]
im=ax.imshow(-np.log10(v),cmap='Blues',vmin=0,vmax=3,aspect='auto')
for i in range(len(v)):
    for j in range(2):ax.text(j,i,f'{v.iloc[i,j]:.3f}',ha='center',va='center',color='white' if v.iloc[i,j]<.03 else 'black')
ax.set_yticks(range(len(v)),v.index);ax.set_xticks([0,1],['Opening','Closing']);ax.set_title('All 25 candidate neighborhoods\nRegional q values | Fixed D206 comparison',fontweight='bold')
fig.colorbar(im,ax=ax,label='-log10 regional q25',shrink=.5)
output(fig,'Figure_S1_all_regions')
print('Four PNG/PDF figure pairs saved to',out)

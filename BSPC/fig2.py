# Fig. 2: quantiles of |timing error| of matched beats (from bspc_incart_diag.m output)
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'serif','font.size':8,'axes.linewidth':0.6})
q=[5,10,25,50,75,90,95,99]
D={'MIT-BIH (2.8 ms/sample)':{'L (proposed)':[0,0,0,2.8,2.8,8.3,22.2,55.6],'G (D-TMA-I)':[0,0,0,2.8,2.8,8.3,25.0,47.2],
   'Elgendi':[0,0,0,0,0,19.4,38.9,72.2],'Pan-Tompkins':[0,0,0,2.8,2.8,8.3,38.9,52.8]},
 'SVDB (7.8 ms/sample)':{'L (proposed)':[7.8,7.8,7.8,15.6,23.4,31.2,31.2,62.5],'G (D-TMA-I)':[0,0,7.8,15.6,15.6,23.4,31.2,54.7],
   'Elgendi':[7.8,7.8,7.8,15.6,23.4,31.2,46.9,62.5],'Pan-Tompkins':[7.8,7.8,15.6,15.6,23.4,46.9,54.7,62.5]},
 'INCART (3.9 ms/sample)':{'L (proposed)':[0,0,3.9,3.9,11.7,19.5,23.3,31.1],'G (D-TMA-I)':[0,0,3.9,7.8,11.7,19.5,23.3,35.0],
   'Elgendi':[0,0,3.9,7.8,15.6,27.2,31.1,62.3],'Pan-Tompkins':[0,0,3.9,7.8,15.6,23.3,50.6,70.0]}}
sty={'L (proposed)':dict(color='k',marker='o',lw=1.4),'G (D-TMA-I)':dict(color='0.55',marker='s',lw=0.9,ls='--'),
     'Elgendi':dict(color='#c0392b',marker='^',lw=1.0),'Pan-Tompkins':dict(color='#1f4e9c',marker='v',lw=1.0,ls=':')}
fig,ax=plt.subplots(1,3,figsize=(7.0,2.3),sharey=True)
for a,(t,d) in zip(ax,D.items()):
    for m,v in d.items(): a.plot(q,v,ms=3,label=m,**sty[m])
    a.set_title(t,fontsize=8); a.set_xlabel('percentile of matched beats'); a.set_xticks([5,25,50,75,95]); a.grid(alpha=0.3,lw=0.4)
ax[0].set_ylabel('|timing error| (ms)')
h,l=ax[0].get_legend_handles_labels(); fig.legend(h,l,loc='lower center',ncol=4,fontsize=7,frameon=False,bbox_to_anchor=(0.5,-0.02))
fig.tight_layout(rect=(0,0.09,1,1)); fig.savefig('fig2.pdf'); fig.savefig('fig2.png',dpi=150); print('ok')

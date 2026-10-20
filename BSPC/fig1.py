# Fig. 1 of the BSPC manuscript (filter: M. Karimi). Exact responses, no data.
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'serif','font.size':8,'axes.linewidth':0.6,'lines.linewidth':1.0})
pi=np.pi
def HK(K,f):
    s=np.sin(pi*f); out=np.ones_like(f); nz=s!=0
    out[nz]=(np.sin(pi*K*f[nz])/(K*s[nz]))**2; return out
def GK(K,f): return np.exp(-2*pi**2*(K*K-1)/6*f*f)
def dog_real(K1,K2,f):
    v1=(K1*K1-1)/6; v2=(K2*K2-1)/6; L=int(np.ceil(4*np.sqrt(v1+v2))); t=np.arange(-L,L+1)
    g2=np.exp(-t**2/(2*v2)); g2/=g2.sum(); g1=np.exp(-t**2/(2*(v1+v2))); g1/=g1.sum()
    return ((g2-g1)[None,:]*np.cos(2*pi*f[:,None]*t[None,:])).sum(1)
fig,ax=plt.subplots(1,2,figsize=(7.0,2.4))
for a,(K1,K2,fs,title) in zip(ax,[(23,2,360,'(a) MIT-BIH median, $K_1=23$ (odd), $K_2=2$'),
                                 (10,2,128,'(b) SVDB median, $K_1=10$ (even), $K_2=2$')]):
    f=np.linspace(0,0.5,4001)
    Hz=HK(K2,f)*(1-HK(K1,f)); D=GK(K2,f)*(1-GK(K1,f))
    a.plot(f*fs,Hz,'k-',label=r'$H_{\rm tot}$ (D-TMA-Z)')
    if K1%2==0:
        Hi=np.abs(HK(K2,f)*(1-np.exp(1j*2*pi*f)*HK(K1,f)))
        a.plot(f*fs,Hi,color='#c0392b',lw=0.9,label=r'$|H_{\rm impl}|$ (D-TMA-I)')
    a.plot(f*fs,D,'--',color='#1f4e9c',label=r'DoG target $D$')
    a.plot(f*fs,dog_real(K1,K2,f),':',color='#1f4e9c',lw=1.1,label='sampled DoG (77/35 taps)')
    a.set_xlim(0,fs/2); a.set_ylim(0,1.05); a.set_xlabel('frequency (Hz)'); a.set_title(title,fontsize=8)
ax[0].set_ylabel('magnitude')
h,l=[],[]
for a in ax:
    for hh,ll in zip(*a.get_legend_handles_labels()):
        if ll not in l: h.append(hh); l.append(ll)
fig.legend(h,l,loc='lower center',ncol=4,fontsize=7,frameon=False,bbox_to_anchor=(0.5,-0.02))
fig.tight_layout(rect=(0,0.08,1,1)); fig.savefig('fig1.pdf'); fig.savefig('fig1.png',dpi=150)
print('ok')

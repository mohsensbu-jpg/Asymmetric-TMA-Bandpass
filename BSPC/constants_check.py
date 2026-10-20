import numpy as np
from scipy.optimize import brentq
pi=np.pi
def HK(K,f):
    f=np.asarray(f,float); out=np.ones_like(f); nz=np.abs(f)>1e-15
    out[nz]=(np.sin(pi*K*f[nz])/(K*np.sin(pi*f[nz])))**2; return out
def GK(K,f,m=2): s2=m*(K**2-1)/12; return np.exp(-2*pi**2*s2*f**2)
gam=lambda t: np.log(np.sin(t)/t)+t**2/6
# 1) 1% thresholds (K-independent bound)
tH=brentq(lambda t: np.exp(-2*gam(t))-1-0.01,0.1,2)/pi   # (G-H)/H
tG=brentq(lambda t: 1-np.exp(2*gam(t))-0.01,0.1,2)/pi     # (G-H)/G
print("1%% thresholds |Kf|: rel-to-H %.5f rel-to-G %.5f"%(tH,tG))
# check actual rel errors <1% at |Kf|<=0.30 for K=2..500
worst=0
for K in range(2,501):
    f=np.linspace(1e-9,0.30/K,2000); r=(GK(K,f)-HK(K,f))/HK(K,f); worst=max(worst,r.max())
print("max rel err (G-H)/H on |Kf|<=0.30, K=2..500: %.5f"%worst)
# 2) certified eps_K: grid max + L*dx/2, L=2pi(K-1)+2pi*sigma/sqrt(e)
def epsK(K,n=2_000_001,m=2):
    f=np.linspace(0,0.5,n); d=np.abs(HK(K,f) if m==2 else 0)
    e=np.max(np.abs(HK(K,f)-GK(K,f))); s=np.sqrt((K**2-1)/6)
    L=2*pi*(K-1)+2*pi*s/np.sqrt(np.e); return e, e+L*(0.5/(n-1))/2
res={}
for K in [2,3,4,5,10,23,50,100,165,500]:
    e,c=epsK(K); res[K]=c; print("K=%d eps_hat=%.4f certified<=%.4f"%(K,e,c))
mx=0;arg=0
for K in range(2,501):
    e,c=epsK(K,400_001)
    if c>mx: mx,arg=c,K
print("max certified over K=2..500: %.4f at K=%d"%(mx,arg))
u=np.linspace(1e-9,40,4_000_001); lim=np.max(np.abs(np.sinc(u)**2-np.exp(-pi**2*u**2/3)))
print("K->inf limit sup|sinc^2(u)-exp(-pi^2u^2/3)| = %.4f"%lim)
# 3) operating point (23,2): bounds
K1,K2=23,2; f=np.linspace(0,0.5,2_000_001)
H1,H2,G1,G2=HK(K1,f),HK(K2,f),GK(K1,f),GK(K2,f)
Ht=H2*(1-H1); D=G2*(1-G1); err=np.abs(Ht-D)
print("(23,2) sup err %.4f ; sum bound %.4f ; relL2 %.4f"%(err.max(),res[23]+res[2],np.sqrt(np.trapezoid((Ht-D)**2,f)/np.trapezoid(D**2,f))))
band=f<1/K1
env=np.maximum(G2*(1-np.exp(2*gam(pi*K2*np.maximum(f,1e-12))))*(1-H1), G2*G1*(1-np.exp(2*gam(pi*K1*np.maximum(f,1e-12)))))
print("common main lobe |f|<1/K1: actual %.4f envelope %.4f, sum-bound %.4f"%(err[band].max(),env[band].max(),res[23]+res[2]))
print("  check envelope>=err on band:",np.all(env[band]+1e-12>=err[band]))
# where is sup attained
print("  argmax f=%.4f (fW-ish), D there=%.3f"%(f[err.argmax()],D[err.argmax()]))
# 4) misalignment: implemented filter for even K1 => H1 advanced by 1 sample
for K1 in [22,23,24,40]:
    H1=HK(K1,f); th=2*pi*f*(K1%2==0)
    Himp=np.abs(HK(2,f)*(1-np.exp(1j*th)*H1)); Hz=HK(2,f)*(1-H1)
    print("K1=%d |Himpl|-|H0| sup %.4f ; bound sup 2 sin(pi f) H1 = %.4f"%(K1,np.max(np.abs(Himp-Hz)),np.max(2*np.sin(pi*f)*H1)))
# 5) fixed-target design for wide factor, even m
for s2 in [(23**2-1)/6,(165**2-1)/6]:
    print("target sigma^2=%.1f"%s2)
    for m in [2,4,6,8]:
        L=max(1,int(round(np.sqrt(12*s2/m+1))))
        Hm=(np.sin(pi*L*f)/(L*np.sin(pi*np.maximum(f,1e-12))))**m; Hm[0]=1
        G=np.exp(-2*pi**2*s2*f**2)
        print("  m=%d L=%d var=%.1f adds/sample=%d sup err=%.4f"%(m,L,m*(L*L-1)/12,2*m,np.max(np.abs(Hm-G))))
# 6) ablation DoG kernel vs target D at (23,2)
K1,K2=23,2; s1=K1/np.sqrt(6); s2=max(K2/np.sqrt(6),0.5); Lk=int(np.ceil(4*s1)); t=np.arange(-Lk,Lk+1)
g2=np.exp(-t**2/(2*s2**2)); g2/=g2.sum(); g1=np.exp(-t**2/(2*s1**2)); g1/=g1.sum(); g=g2-g1
ff=np.linspace(0,0.5,200001); Hg=(g[None,:]*np.cos(2*pi*ff[:,None]*t[None,:])).sum(1)
Dt=GK(K2,ff)*(1-GK(K1,ff)); Htma=HK(K2,ff)*(1-HK(K1,ff))
print("ablation DoG taps=%d sup|H_abl-D|=%.4f sup|H_abl-H_tma|=%.4f"%(len(t),np.max(np.abs(Hg-Dt)),np.max(np.abs(Hg-Htma))))

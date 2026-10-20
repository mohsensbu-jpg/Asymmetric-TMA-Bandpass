import numpy as np, math
pi=np.pi; n=2_000_001; f=np.linspace(0,0.5,n); D=0.5/(n-1)
sf=np.sin(pi*f); sf[0]=1
EVAL=1e-12  # allowance for floating-point evaluation error
def up(v,d=4): return math.ceil(v*10**d)/10**d
lo={};hi={}
for K in range(2,501):
    H=(np.sin(pi*K*f)/(K*sf))**2; H[0]=1
    s2=(K*K-1)/6; e=np.max(np.abs(H-np.exp(-2*pi**2*s2*f*f)))
    lo[K]=e; hi[K]=e+EVAL+D/2*(2*pi*(K-1)+2*pi*np.sqrt(s2)/np.sqrt(np.e))
for K in [2,3,4,5,10,23,50,100,165,500]:
    print(K,"grid %.6f  certified upper %.6f -> %.4f"%(lo[K],hi[K],up(hi[K])))
print("max upper K=3..500: %.6f ; lower K=2: %.6f"%(max(hi[k] for k in range(3,501)),lo[2]))
print("upper 11..500 %.6f ; 50..500 %.6f"%(max(hi[k] for k in range(11,501)),max(hi[k] for k in range(50,501))))
# DoG theorem-matched realisation error
for K1 in [10,23]:
    K2=2; v1=(K1*K1-1)/6; v2=0.5; L=int(np.ceil(4*np.sqrt(v1+v2))); t=np.arange(-L,L+1)
    g2=np.exp(-t**2/(2*v2)); g2/=g2.sum(); g1=np.exp(-t**2/(2*(v1+v2))); g1/=g1.sum(); g=g2-g1
    ff=np.linspace(0,0.5,100001); Hg=(g[None,:]*np.cos(2*pi*ff[:,None]*t[None,:])).sum(1)
    Dt=np.exp(-2*pi**2*v2*ff**2)-np.exp(-2*pi**2*(v1+v2)*ff**2)
    print("K1=%d taps %d sup|H_real - D| = %.4f"%(K1,len(t),np.max(np.abs(Hg-Dt))))
# main table for Z
for name,B,FP,FN in [("MIT",109494,794,326),("SVDB",184582,1344,1050),("MIT-I",109494,799,299),("SVDB-I",184582,2202,1020)]:
    TP=B-FN; print(name,"Se %.4f +P %.4f F1 %.5f"%(TP/B,TP/(TP+FP),2*TP/(2*TP+FP+FN)))
# envelope more digits
K1,K2=23,2; ff=np.linspace(1e-9,1/23-1e-9,400001)
H=lambda K:(np.sin(pi*K*ff)/(K*np.sin(pi*ff)))**2; G=lambda K:np.exp(-2*pi**2*(K*K-1)/6*ff**2)
gam=lambda t:np.log(np.sin(t)/t)+t*t/6
err=np.abs(H(2)*(1-H(23))-G(2)*(1-G(23)))
env=np.maximum(G(2)*(1-np.exp(2*gam(pi*2*ff)))*(1-H(23)),G(2)*G(23)*(1-np.exp(2*gam(pi*23*ff))))
print("band max err %.6f env %.6f; pointwise min ratio err/env on f>0.002: %.3f"%(err.max(),env.max(),np.min(err[ff>0.002]/env[ff>0.002])))

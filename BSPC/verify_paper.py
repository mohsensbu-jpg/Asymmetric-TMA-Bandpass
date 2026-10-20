#!/usr/bin/env python3
"""
verify_paper.py -- independent check of the numbers in the manuscript
"QRS detection with a running-sum triangular-filter cascade and a local
decision rule without learned parameters" (M. Karimi, BSPC submission).

Part A re-derives every theoretical constant from first principles (needs only
numpy).  Part B checks the internal arithmetic of every results table and of the
numbers quoted in the text/abstract against the raw counts printed by the
MATLAB scripts.  Each check prints PASS or FAIL; the summary at the end counts
failures.  Run:  python3 verify_paper.py
"""
import math
import numpy as np

pi = np.pi
FAIL = []
# np.trapezoid exists from numpy 2.0; older numpy has np.trapz
trapz = getattr(np, "trapezoid", None) or np.trapz


def check(name, ok, detail=""):
    print(("PASS  " if ok else "FAIL  ") + name + (("   " + detail) if detail else ""))
    if not ok:
        FAIL.append(name)


def close(a, b, tol):
    return abs(a - b) <= tol


# ----------------------------------------------------------------------------
# PART A  theory
# ----------------------------------------------------------------------------
print("\n=== PART A: theory ===")


def HK(K, f):
    """exact response of the centred triangular kernel (two length-K boxcars)"""
    f = np.asarray(f, float)
    s = np.sin(pi * f)
    out = np.ones_like(f)
    nz = np.abs(s) > 1e-15
    out[nz] = (np.sin(pi * K * f[nz]) / (K * s[nz])) ** 2
    return out


def GK(K, f, var=None):
    v = (K * K - 1) / 6.0 if var is None else var
    return np.exp(-2 * pi ** 2 * v * np.asarray(f, float) ** 2)


def gam(t):
    t = np.asarray(t, float)
    return np.log(np.sin(t) / t) + t ** 2 / 6


# A1: H_K equals the DTFT of the triangular kernel with weights (K-|n|)/K^2
f = np.linspace(0, 0.5, 2001)
ok = True
for K in [2, 3, 4, 7, 10, 23]:
    n = np.arange(-(K - 1), K)
    h = (K - np.abs(n)) / K ** 2
    H = (h[None, :] * np.cos(2 * pi * f[:, None] * n[None, :])).sum(1)
    ok &= np.max(np.abs(H - HK(K, f))) < 1e-12
    ok &= close((h * n ** 2).sum(), (K * K - 1) / 6, 1e-12)   # variance of the kernel
check("A1 squared-Dirichlet response and variance (K^2-1)/6 of the triangular kernel", ok)

# A2: gamma series = -sum_{n>=2} zeta(2n) t^2n / (n pi^2n)
t = np.linspace(0.05, 2.5, 50)
zeta = lambda s: float(np.sum(np.arange(1, 20000, dtype=float) ** (-float(s))))
series = -sum(zeta(2 * m) * t ** (2 * m) / (m * pi ** (2 * m)) for m in range(2, 60))
check("A2 zeta-series of gamma(t)", np.max(np.abs(series - gam(t))) < 1e-8)

# A3: Lemma 1 sandwich G e^{2 gamma(pi K f)} <= H_K <= G on 0<|f|<1/K
worst = 0.0
for K in range(2, 201):
    ff = np.linspace(1e-6, min(0.5, 1.0 / K - 1e-6), 4000)
    H, G = HK(K, ff), GK(K, ff)
    lo = G * np.exp(2 * gam(pi * K * ff))
    worst = max(worst, np.max(H - G), np.max(lo - H))
check("A3 Lemma 1 two-sided bound, K=2..200", worst < 1e-12, "max violation %.1e" % worst)

# A4: Corollary 1 constants
v1 = np.exp(-2 * gam(0.3 * pi)) - 1
v2 = 1 - np.exp(2 * gam(0.3 * pi))
check("A4 Corollary 1: 0.934 % (rel. to H) and 0.925 % (rel. to G)",
      v1 < 0.00935 and v2 < 0.00926 and close(v1, 0.00934, 1e-4) and close(v2, 0.00925, 1e-4),
      "%.5f %.5f" % (v1, v2))
lo_, hi_ = 0.1, 1.5
for _ in range(100):
    mid = (lo_ + hi_) / 2
    lo_, hi_ = (mid, hi_) if np.exp(-2 * gam(mid)) - 1 < 0.01 else (lo_, mid)
t1 = lo_ / pi
lo_, hi_ = 0.1, 1.5
for _ in range(100):
    mid = (lo_ + hi_) / 2
    lo_, hi_ = (mid, hi_) if 1 - np.exp(2 * gam(mid)) < 0.01 else (lo_, mid)
t2 = lo_ / pi
check("A4b exact 1 % limits |Kf| = 0.30498 and 0.30572", close(t1, 0.30498, 1e-4) and close(t2, 0.30572, 1e-4),
      "%.5f %.5f" % (t1, t2))

# A5: whole-band constants (Table of eps_K) with Bernstein certification
N = 2_000_001
ff = np.linspace(0, 0.5, N)
D = 0.5 / (N - 1)
table = {2: (0.11370, 0.1138), 3: (0.10972, 0.1098), 4: (0.07645, 0.0765), 5: (0.07306, 0.0731),
         10: (0.06882, 0.0689), 23: (0.06773, 0.0678), 100: (0.06749, 0.0676), 500: (0.06748, 0.0680)}
ok = True
for K, (g_rep, u_rep) in table.items():
    e = np.max(np.abs(HK(K, ff) - GK(K, ff)))
    up = e + 1e-12 + D / 2 * (2 * pi * (K - 1) + 2 * pi * math.sqrt((K * K - 1) / 6) / math.sqrt(math.e))
    # the paper reports the certified bound rounded UP to 4 decimals: it must be >= up and within 1e-4
    ok &= close(e, g_rep, 1e-5) and (u_rep >= up - 1e-12) and (u_rep - up < 1e-4 + 1e-12)
check("A5 eps_K table (grid value and outward-rounded certified bound)", ok)
up3 = np.max(np.abs(HK(3, ff) - GK(3, ff))) + 1e-12 + D / 2 * (2 * pi * 2 + 2 * pi * math.sqrt(8 / 6) / math.sqrt(math.e))
check("A5b text: certified bound for K=3 is at most 0.10973 (5 decimals, rounded up)", up3 <= 0.10973 and 0.10973 - up3 < 1e-5,
      "%.6f" % up3)

# A6: operating point (K1,K2) = (23,2): sup|Htot-D| 0.114, rel L2 0.116, bound 0.182, main lobe 0.06695/0.06697
H1, H2 = HK(23, ff), HK(2, ff)
Ht = H2 * (1 - H1)
Dd = GK(2, ff) * (1 - GK(23, ff))
err = np.abs(Ht - Dd)
l2 = math.sqrt(trapz((Ht - Dd) ** 2, ff) / trapz(Dd ** 2, ff))
check("A6 (23,2): sup error 0.114 near f=0.43, rel. L2 0.116, bound 0.1138+0.0678=0.182",
      close(err.max(), 0.1137, 2e-4) and close(ff[err.argmax()], 0.43, 0.01) and close(l2, 0.116, 1e-3)
      and close(0.1138 + 0.0678, 0.1816, 1e-9))
band = (ff > 0) & (ff < 1 / 23)
fb = ff[band]
env = np.maximum(GK(2, fb) * (1 - np.exp(2 * gam(pi * 2 * fb))) * (1 - HK(23, fb)),
                 GK(2, fb) * GK(23, fb) * (1 - np.exp(2 * gam(pi * 23 * fb))))
check("A6b main-lobe bound 0.06697 vs measured 0.06695; bound >= error everywhere",
      close(env.max(), 0.06697, 1e-4) and close(err[band].max(), 0.06695, 1e-4) and np.all(env + 1e-12 >= err[band]))
check("A6c Prop. 1 bound <= 0.228 for K<=500 and <= 0.138 for K>10",
      0.1138 * 2 <= 0.228 and 2 * 0.0686 <= 0.138)

# A7: offset (D-TMA-I) terms delta_K and measured magnitude changes
ok = True
for K1, dl, ms in [(22, 0.066, 0.011), (10, 0.145, 0.047), (8, 0.182, 0.070)]:
    H1 = HK(K1, ff)
    delta = np.max(2 * np.abs(np.sin(pi * ff)) * H1)
    Hi = np.abs(HK(2, ff) * (1 - np.exp(1j * 2 * pi * ff) * H1))
    meas = np.max(np.abs(Hi - HK(2, ff) * (1 - H1)))
    ok &= close(delta, dl, 1.5e-3) and close(meas, ms, 1.5e-3)
check("A7 offset terms delta = 0.066/0.145/0.182, measured 0.011/0.047/0.070 (K1 = 22/10/8)", ok)

# A8: accuracy-cost table (m boxcars of integer length L vs fixed Gaussian)
ok = True
rep = {88.0: {2: (23, 0.0677), 4: (16, 0.0215), 6: (13, 0.0093), 8: (12, 0.0409)},
       (165 ** 2 - 1) / 6: {2: (165, 0.0675), 4: (117, 0.0314), 6: (95, 0.0178), 8: (83, 0.0176)}}
fm = np.linspace(1e-9, 0.5, 400001)
for s2, row in rep.items():
    for m, (Lr, er) in row.items():
        L = int(round(math.sqrt(12 * s2 / m + 1)))
        Hm = (np.sin(pi * L * fm) / (L * np.sin(pi * fm))) ** m
        e = np.max(np.abs(Hm - GK(0, fm, var=s2)))
        ok &= (L == Lr) and close(e, er, 1.5e-4)
check("A8 accuracy-cost table", ok)

# A9: integer realisation: 32-bit range and the 529 = 2^9+2^4+1 example.
# The filter is LINEAR in the samples (no squaring): with |x| <= 2^11 (12-bit
# signed), |S_K1 x| <= K1^2 2^11, so |K1^2 x - S_K1 x| <= 2 K1^2 2^11, and S_K2
# (weights sum K2^2 = 4) multiplies the bound by 4: total <= K1^2 * 2^14.
check("A9 32-bit accumulator holds for K1<=362 (K2=2, 12-bit samples); 529 = 512+16+1",
      362 ** 2 * 2 ** 14 < 2 ** 31 <= 363 ** 2 * 2 ** 14 and 23 ** 2 == 2 ** 9 + 2 ** 4 + 1)

# ----------------------------------------------------------------------------
# PART B  arithmetic of the reported results
# ----------------------------------------------------------------------------
print("\n=== PART B: reported results ===")
F1 = lambda tp, fp, fn: 2 * tp / (2 * tp + fp + fn)

# B1: Table 'test' -- (beats, Se, +P, F1, FP, FN) as printed, recomputed from FP, FN
rows = {
    ("MIT-BIH", "L"): (109494, .9980, .9949, .9965, 559, 219),
    ("MIT-BIH", "G"): (109494, .9973, .9927, .9950, 799, 299),
    ("MIT-BIH", "E"): (109494, .9987, .9957, .9972, 468, 146),
    ("MIT-BIH", "PT"): (109494, .9916, .9919, .9918, 885, 915),
    ("SVDB", "L"): (184582, .9949, .9943, .9946, 1048, 935),
    ("SVDB", "G"): (184582, .9945, .9882, .9913, 2201, 1019),
    ("SVDB", "E"): (184582, .9966, .9966, .9966, 625, 627),
    ("SVDB", "PT"): (184582, .9849, .9851, .9850, 2742, 2788),
    ("INCART", "L"): (175900, .9862, .9819, .9840, 3201, 2425),
    ("INCART", "G"): (175900, .9878, .9694, .9785, 5491, 2139),
    ("INCART", "E"): (175900, .9981, .9929, .9955, 1263, 340),
    ("INCART", "PT"): (175900, .9784, .9838, .9811, 2835, 3795),
    ("EDB", "L"): (790565, .9874, .9823, .9849, 14047, 9953),
    ("EDB", "E"): (790565, .9933, .9880, .9907, 9532, 5272),
    ("EDB", "PT"): (790565, .9726, .9636, .9681, 29087, 21648),
}
bad = []
for (db, m), (B, se, pp, f1, fp, fn) in rows.items():
    tp = B - fn
    if not (close(tp / B, se, 1e-4) and close(tp / (tp + fp), pp, 1e-4) and close(F1(tp, fp, fn), f1, 1e-4)):
        bad.append(f"{db}/{m}")
check("B1 Se, +P, F1 consistent with beats, FP, FN (Table test)", not bad, ", ".join(bad))

# B2: matched beats in the localisation table = beats - FN
matched = {("MIT-BIH", "L"): 109275, ("MIT-BIH", "G"): 109195, ("MIT-BIH", "E"): 109348, ("MIT-BIH", "PT"): 108579,
           ("SVDB", "L"): 183647, ("SVDB", "G"): 183563, ("SVDB", "E"): 183955, ("SVDB", "PT"): 181794,
           ("INCART", "L"): 173475, ("INCART", "G"): 173761, ("INCART", "E"): 175560, ("INCART", "PT"): 172105,
           ("EDB", "L"): 780612, ("EDB", "E"): 785293, ("EDB", "PT"): 768917}
bad = [f"{k}" for k, v in matched.items() if rows[k][0] - rows[k][5] != v]
check("B2 matched beats = beats - FN (Table loc)", not bad, ", ".join(bad))

# B3: differences quoted in the text
check("B3 Elgendi - L gross F1: 0.002 (SVDB), 0.011 (INCART), 0.006 (EDB)",
      close(.9966 - .9946, 0.002, 6e-4) and close(.9955 - .9840, 0.011, 6e-4) and close(.9907 - .9849, 0.0058, 1e-4))
# B5: NSTDB table F1 recomputed from counts (1856 beats per SNR, 1785 at -6 dB)
nst = {  # SNR: method: (F1, FP, FN)
    24: {"L": (.9997, 1, 0), "G": (.9995, 2, 0), "E": (.9992, 3, 0), "PT": (.9995, 2, 0)},
    18: {"L": (.9960, 15, 0), "G": (.9938, 23, 0), "E": (.9745, 97, 0), "PT": (.9989, 4, 0)},
    12: {"L": (.9170, 335, 1), "G": (.8777, 517, 0), "E": (.9028, 397, 2), "PT": (.9057, 384, 2)},
    6: {"L": (.7862, 923, 56), "G": (.7492, 1171, 43), "E": (.7676, 1016, 67), "PT": (.7204, 1265, 99)},
}
bad = []
for snr, d in nst.items():
    for m, (f1, fp, fn) in d.items():
        if not close(F1(1856 - fn, fp, fn), f1, 1e-4):
            bad.append(f"{snr}dB/{m}")
check("B5 NSTDB F1 consistent with FP, FN and 1856 beats", not bad, ", ".join(bad))
# B6: INCART missed-beat diagnosis
fnc = {"N": (153672, 2186, 197), "SVEB": (1992, 20, 4), "VEB": (20011, 214, 139), "F": (219, 5, 0), "Q": (6, 0, 0)}
check("B6 AAMI classes sum to 175900 beats, FN 2425 (L) and 340 (E)",
      sum(v[0] for v in fnc.values()) == 175900 and sum(v[1] for v in fnc.values()) == 2425
      and sum(v[2] for v in fnc.values()) == 340)
check("B6b miss rates 1.42/1.00/1.07 % (L) and 0.13/0.20/0.69 % (E)",
      close(100 * 2186 / 153672, 1.42, 0.005) and close(100 * 20 / 1992, 1.00, 0.005) and close(100 * 214 / 20011, 1.07, 0.005)
      and close(100 * 197 / 153672, 0.13, 0.005) and close(100 * 4 / 1992, 0.20, 0.005) and close(100 * 139 / 20011, 0.69, 0.005))
check("B6c causes 1928+462+35+0 = 2425 (80 %, 19 %, 1.4 %); 5-min blocks sum to 2425; I03 1209 = 50 %",
      1928 + 462 + 35 == 2425 and close(1928 / 2425, 0.80, 0.006) and close(462 / 2425, 0.19, 0.006)
      and close(35 / 2425, 0.014, 0.001) and sum([428, 546, 275, 403, 384, 389]) == 2425 and close(1209 / 2425, 0.50, 0.005))
check("B6d I03: W = 30 samples = 117 ms, w1 = 99 = 385 ms at 257 Hz",
      close(30 / 257 * 1000, 117, 0.5) and close(3.3 * 30, 99.0, 0.5) and close(99 / 257 * 1000, 385, 0.5))

# B7: width statistics
check("B7 W medians: 14 samples = 38.9 ms at 360 Hz; 5.5 samples = 43.0 ms at 128 Hz",
      close(14 / 360 * 1000, 38.9, 0.05) and close(5.5 / 128 * 1000, 43.0, 0.05))

# B8: BUT QDB proportions
check("B8 BUT QDB classes 49.3 + 34.5 + 16.3 = 100 %", close(49.3 + 34.5 + 16.3, 100, 0.15))

# B9: EDB offsets and sampling intervals
check("B9 sample intervals 2.8 / 7.8 / 3.9 / 4.0 ms", close(1000 / 360, 2.8, 0.05) and close(1000 / 128, 7.8, 0.05)
      and close(1000 / 257, 3.9, 0.05) and close(1000 / 250, 4.0, 0.05))

print("\nSUMMARY: %d failed check(s)%s" % (len(FAIL), (": " + "; ".join(FAIL)) if FAIL else ""))
print("numpy version:", np.__version__)

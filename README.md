
# Asymmetric-TMA-Bandpass

Running-sum triangular moving-average (TMA) difference filter by Mohsen Karimi
(ORCID 0009-0004-1055-8267).

    x_p = TMA( x - TMA(x, w1), w2 )

An exact band-pass filter with DC rejection and independently controlled band
edges, realisable with integer additions/subtractions only, and close to a
difference of Gaussians (DoG) with certified error bounds.

## Contents

- `tma_detect.m`, `reproduce_tables.m`, `references.bib` — code of the
  IEEE Signal Processing Letters version (release `v1.0-SPL`).
- `BSPC/` — code for the paper "QRS detection with a running-sum
  triangular-filter cascade and a local decision rule without learned
  parameters" (Biomedical Signal Processing and Control, 2026): QRS detector,
  experiments on MIT-BIH, SVDB, INCART, EDB and NSTDB, and Python scripts that
  verify the theory constants. See `BSPC/README.md` for run order.

No PhysioNet data are stored in this repository; the scripts download them.

## Third-party code

Pan-Tompkins comparator: H. Sedghamiz, MATLAB File Exchange 45840 (not included).

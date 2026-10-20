# BSPC paper — reproduction code

Code for the paper *"QRS detection with a running-sum triangular-filter cascade
and a local decision rule without learned parameters"* (Mohsen Karimi,
submitted to Biomedical Signal Processing and Control, 2026).

Author: Mohsen Karimi (ORCID 0009-0004-1055-8267). The TMA-difference filter
and the detector design are the author's work; see the paper's declaration for
the role of AI tools.

## Filter

    x_p = TMA( x - TMA(x, w1), w2 ),   w1 = 3.3 W,  w2 = 0.06 W  (each >= 3 samples)

TMA = two cascaded length-K boxcars (running sums), K = ceil(w/2).
D-TMA-I is the original realisation; D-TMA-Z is the zero-phase realisation.

## Proposed detector (frozen configuration, selected on MIT-BIH only)

- D-TMA-Z filter
- local level mu = 10 s moving median, local scale s = moving std of (x_p - mu)
- candidates: local maxima of |x_p - mu| above 1.75 s
- chain merge within w1, keep the largest e/s
- within 200 ms keep the larger peak; 200-360 ms reject if e/s < 0.4 x previous beat

## Data

No data are stored here. The `get_*.m` scripts download the PhysioNet records
(MIT-BIH, SVDB, INCART, EDB, NSTDB) into the MATLAB working folder.
`read_mitdb.m` reads format-212 records (0-based sample indices).

## Run order (MATLAB)

1. `get_mitdb_full`, `get_svdb`, `get_nstdb`, `get_incart`, `get_edb`
2. `bspc_dev_test` — 48 configurations on MIT-BIH (development), test on the others
3. `bspc_fix` — merge variants, EDB, INCART record I03
4. `bspc_spread`, `bspc_timing_table`, `bspc_bootstrap` — timing and confidence intervals
5. `bspc_runs`, `bspc_runs2`, `bspc_incart_diag`, `ablation_test`, `butqdb_quality` — remaining tables

All outputs are printed as text.

## Theory check (Python 3, numpy, matplotlib)

    python constants_check.py   # constants of Lemma 1 / Corollary 1
    python certify_eps.py       # certified eps_K bounds
    python verify_paper.py      # checks all theory constants and reported numbers
    python fig1.py; python fig2.py

## Comparator

Pan-Tompkins: H. Sedghamiz, MATLAB File Exchange 45840 (not included; download
`pan_tompkin.m` from File Exchange and put it on the MATLAB path).
The Elgendi comparator is implemented inside the scripts.

## Earlier version

The IEEE SPL version of this code is the release tagged `v1.0-SPL`.

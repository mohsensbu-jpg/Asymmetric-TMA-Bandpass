function bspc_dev_test()
%BSPC_DEV_TEST  Clean development/test protocol for a training-free local
%   decision stage on the D-TMA-Z filter (filter, detector and the ideas of a
%   local noise estimate and physiological QRS/T constraints: M. Karimi).
%
%   STEP 1 (development, MIT-BIH only): 48 prespecified configurations of the
%     decision stage are scored on the 48 MIT-BIH records. The configuration
%     with the highest gross F1 (+-100 ms) is selected automatically; ties
%     within 0.0002 go to the configuration with fewer optional rules, then to
%     the longer window. Nothing from SVDB, INCART or NSTDB is used.
%   STEP 2 (one-shot test): the frozen configuration, D-TMA-I (global, frozen
%     earlier), Elgendi (our implementation) and Pan-Tompkins on SVDB,
%     INCART (lead II, 257 Hz, never used before) and NSTDB noisy segments.
%   Needs: read_mitdb.m, pan_tompkin.m, folders mitdb_full, svdb, nstdb and
%   incart (run get_incart first).
%
%   Decision stage (all options fixed in advance):
%     e = |xp - mu|, mu = moving median of xp over Lw seconds,
%     scale s = moving std  (alpha in {1.5, 1.75, 2.0})
%            or 1.4826 * moving median |xp - mu|  (alpha in {4, 6, 8}),
%     Lw in {5, 10} s;  candidates = local maxima of e above alpha*s,
%     chain merge within w1 keeping the largest e/s;
%     option WIDTH: reject events whose half-amplitude width of e is outside
%       [0.3W, 3W];
%     option RR/T: within 200 ms of the previous accepted beat keep only the
%       larger; between 200 and 360 ms reject if e/s < 0.4 of the previous beat
%       (a T wave or artefact should not be counted as a second beat).
root = userpath;  if isempty(root), root = pwd; end
if exist('pan_tompkin', 'file') ~= 2
    f = dir(fullfile(root, '**', 'pan_tompkin.m'));  addpath(f(1).folder);
end
tic;
CFG = [];
for Lw = [10 5]
    for sc = 1:2
        al = [1.5 1.75 2.0];  if sc == 2, al = [4 6 8]; end
        for a = al
            for wd = 0:1
                for rr = 0:1
                    CFG = [CFG; Lw sc a wd rr];   %#ok<AGROW>
                end
            end
        end
    end
end
nC = size(CFG, 1);
mrec = [100:109 111:119 121:124 200:203 205 207:210 212:215 217 219:223 228 230:234];
C = zeros(nC, 3);
fprintf('\nSTEP 1  development on MIT-BIH (%d configurations)\n', nC);
for k = 1:numel(mrec)
    [x, gt] = read_mitdb(fullfile(root, 'mitdb_full'), mrec(k));  x = x - median(x);
    F = front(x, 360);
    P = cell(2, 2);   % (window, scale) statistics computed once per record
    for c = 1:nC
        iw = 1 + (CFG(c,1) == 5);  is = CFG(c,2);
        if isempty(P{iw,is}), P{iw,is} = dstats(F, CFG(c,1), is, 360); end
        d = decide(F, P{iw,is}, CFG(c,:), 360);
        [tp, fp, fn] = score_j(d - 1, gt, round(0.1*360));  C(c,:) = C(c,:) + [tp fp fn];
    end
end
F1 = 2*C(:,1) ./ (2*C(:,1) + C(:,2) + C(:,3));
nopt = CFG(:,4) + CFG(:,5);
fprintf('  Lw  scale  alpha  width  RR/T    F1       FP     FN\n');
sn = {'std', 'MAD'};
for c = 1:nC
    fprintf('  %2d  %-4s  %5.2f   %d      %d     %.5f  %5d  %5d\n', CFG(c,1), sn{CFG(c,2)}, CFG(c,3), CFG(c,4), CFG(c,5), F1(c), C(c,2), C(c,3));
end
best = max(F1);  cand = find(F1 >= best - 0.0002);
[~, o] = sortrows([nopt(cand) -CFG(cand,1) -F1(cand)]);  sel = cand(o(1));
S = CFG(sel,:);
fprintf('  SELECTED: Lw %d s, scale %s, alpha %.2f, width %d, RR/T %d  (MIT-BIH F1 %.5f)\n', ...
    S(1), sn{S(2)}, S(3), S(4), S(5), F1(sel));
fprintf('  development time %.0f s\n', toc);

%% STEP 2  one-shot test
NM = {'TMA-Z local (sel.)', 'TMA-I global', 'Elgendi', 'Pan-Tompkins'};
DBS = {'MIT-BIH (dev)', 'SVDB', 'INCART'};
for b = 1:3
    switch b
        case 1, recs = num2cell(mrec);  fs = 360;
        case 2, recs = num2cell([800:812 820:829 840:894]);  fs = 128;
        case 3   % INCART: every record that get_incart has finished (order I01..I75)
            have = arrayfun(@(r) exist(fullfile(root, 'incart', sprintf('I%02d.mat', r)), 'file') == 2, 1:75);
            recs = num2cell(find(have));  fs = 257;
            if isempty(recs), fprintf('\nINCART: no records yet, skipped\n'); continue; end
    end
    m = numel(recs);  Cg = zeros(m, 4, 2, 3);  Fr = zeros(m, 4);  J = cell(4,1);
    for k = 1:m
        switch b
            case 1, [x, gt] = read_mitdb(fullfile(root, 'mitdb_full'), recs{k});
            case 2, [x, gt] = read_mitdb(fullfile(root, 'svdb'), recs{k});
            case 3, Sx = load(fullfile(root, 'incart', sprintf('I%02d.mat', recs{k})));
                    x = double(Sx.x);  gt = Sx.beats;  fs = Sx.fs;
        end
        x = x - median(x);
        for j = 1:4
            switch j
                case 1, Fx = front(x, fs);  d = decide(Fx, dstats(Fx, S(1), S(2), fs), S, fs);
                case 2, d = detect(x, fs, 1);
                case 3, d = elgendi(x, fs);
                case 4, d = pt_orig(x, fs);
            end
            d = d - 1;   % MATLAB 1-based index -> 0-based WFDB sample number
            for t = 1:2
                tol = round([0.100 0.050]*fs);
                [tp, fp, fn, dt] = score_j(d, gt, tol(t));  Cg(k,j,t,:) = [tp fp fn];
                if t == 1, Fr(k,j) = 2*tp/max(2*tp+fp+fn,1);  J{j} = [J{j}; 1000*dt/fs]; end
            end
        end
    end
    fprintf('\n%s (%d records, %d beats, %d Hz)\n', DBS{b}, m, sum(Cg(:,1,1,1)+Cg(:,1,1,3)), fs);
    fprintf('  method               Se       +P       F1(100)  F1(50)   FP     FN   medF1   |err| med/p95 ms   vs TMA-Z local\n');
    for j = 1:4
        c = squeeze(sum(Cg(:,j,1,:),1));  c2 = squeeze(sum(Cg(:,j,2,:),1));
        s = '-';
        if j > 1
            dd = Fr(:,j) - Fr(:,1);  nb = sum(dd > 1e-9);  nw = sum(dd < -1e-9);
            s = sprintf('%d/%d p=%.3f', nb, nw, sign_test(nb, nw));
        end
        fprintf('  %-19s  %.5f  %.5f  %.5f  %.5f  %5d  %5d  %.4f  %5.1f / %5.1f   %s\n', NM{j}, ...
            c(1)/(c(1)+c(3)), c(1)/(c(1)+c(2)), 2*c(1)/(2*c(1)+c(2)+c(3)), ...
            2*c2(1)/(2*c2(1)+c2(2)+c2(3)), c(2), c(3), median(Fr(:,j)), median(J{j}), qtl(J{j}, 0.95), s);
    end
end
fN = fullfile(root, 'nstdb');  fs = 360;  tol = round(0.100*fs);
snr = {'24','18','12','06','00','_6'};  lab = [24 18 12 6 0 -6];
fprintf('\nNSTDB noisy segments, 118+119 pooled, +-100 ms\n  SNR'); fprintf('%20s', NM{:}); fprintf('\n');
for s = 1:numel(snr)
    Cn = zeros(4, 3);
    for r = {'118','119'}
        [x, gt] = read_mitdb(fN, [r{1} 'e' snr{s}]);  x = x - median(x);
        mask = noisy_mask(numel(x), fs);
        for j = 1:4
            switch j
                case 1, Fx = front(x, fs);  d = decide(Fx, dstats(Fx, S(1), S(2), fs), S, fs);
                case 2, d = detect(x, fs, 1);
                case 3, d = elgendi(x, fs);
                case 4, d = pt_orig(x, fs);
            end
            [tp, fp, fn] = masked_counts(d - 1, gt, tol, mask);  Cn(j,:) = Cn(j,:) + [tp fp fn];
        end
    end
    fprintf('  %3d', lab(s)); fprintf('%20.4f', 2*Cn(:,1)./(2*Cn(:,1)+Cn(:,2)+Cn(:,3))); fprintf('\n');
    fprintf('   FP'); fprintf('%20d', Cn(:,2)); fprintf('\n');
    fprintf('   FN'); fprintf('%20d', Cn(:,3)); fprintf('\n');
end
fprintf('\n  total time %.0f s, MATLAB %s\n\n', toc, version);
end

function F = front(x, fs)
% D-TMA-Z filter and the same automatic width estimate as before
x = double(x(:));  N = numel(x);
W  = estimate_width(x, fs);
W  = min(max(3, round(W)), max(3, round(2*fs)));
F.W  = W;
F.w1 = min(max(3, round(3.3*W)),  max(3, round(2*fs)));
w2   = min(max(3, round(0.06*W)), max(3, round(2*fs)));
F.xp = tma0(x - tma0(x, F.w1), w2);
F.edge = min(max(round(0.01*fs), 2*F.w1), floor((N-3)/2));
end

function P = dstats(F, Lw, sc, fs)
xp = F.xp;  L = min(round(Lw*fs), numel(xp));
P.mu = movmedian(xp, L, 'Endpoints', 'shrink');
P.e  = abs(xp - P.mu);
if sc == 1
    m1  = movmean(xp - P.mu, L, 'Endpoints', 'shrink');
    P.s = sqrt(max(movmean((xp - P.mu).^2, L, 'Endpoints', 'shrink') - m1.^2, eps));
else
    P.s = max(1.4826 * movmedian(P.e, L, 'Endpoints', 'shrink'), eps);
end
end

function d = decide(F, P, cfg, fs)
al = cfg(3);  wd = cfg(4);  rr = cfg(5);
N = numel(F.xp);  e = P.e;  s = P.s;
i = (F.edge+1 : N-F.edge)';
c = i(e(i) > e(i-1) & e(i) >= e(i+1) & e(i) > al*s(i));
a = e ./ s;
d = merge_chain(a, c, F.w1);
if isempty(d), return; end
if wd
    keep = true(numel(d), 1);
    for k = 1:numel(d)
        h = 0.5*e(d(k));  l = d(k);  r = d(k);
        while l > 1 && e(l) > h, l = l - 1; end
        while r < N && e(r) > h, r = r + 1; end
        wk = r - l;  keep(k) = wk >= 0.3*F.W && wk <= 3*F.W;
    end
    d = d(keep);
end
if rr && numel(d) > 1
    r1 = round(0.200*fs);  r2 = round(0.360*fs);
    out = d(1);  ao = a(d(1));
    for k = 2:numel(d)
        g = d(k) - out(end);  ak = a(d(k));
        if g < r1
            if ak > ao, out(end) = d(k); ao = ak; end
        elseif g < r2 && ak < 0.4*ao
            continue
        else
            out(end+1, 1) = d(k);  ao = ak;   %#ok<AGROW>
        end
    end
    d = out;
end
end

function q = qtl(v, p)
v = sort(v(:));  q = v(max(1, round(p*numel(v))));
end

function d = run_method(j, x, fs)
switch j
    case 1, d = detect(x, fs, 1);
    case 2, d = detect_local(x, fs);
    case 3, d = elgendi(x, fs);
end
end

function dets = detect_local(x, fs)
% same filter and width estimate as detect(...,1); threshold mu +- 1.75 sigma
% with mu, sigma over a sliding 10-s window (V5r-loc rule, fixed beforehand)
x = double(x(:));  N = numel(x);
W  = estimate_width(x, fs);
W  = min(max(3, round(W)), max(3, round(2*fs)));
w1 = min(max(3, round(3.3*W)),  max(3, round(2*fs)));
w2 = min(max(3, round(0.06*W)), max(3, round(2*fs)));
xp = tma(x - tma(x, w1), w2);
edge = min(max(round(0.01*fs), 2*w1), floor((N-3)/2));
Lw = min(round(10*fs), N);
mu = movmean(xp, Lw, 'Endpoints', 'shrink');
sg = sqrt(max(movmean(xp.^2, Lw, 'Endpoints', 'shrink') - mu.^2, eps));
i  = (edge+1 : N-edge)';
xc = xp(i);  xl = xp(i-1);  xr = xp(i+1);
cand = sort([ i(xc > xl & xc >= xr & xc > mu(i) + 1.75*sg(i)) ; ...
              i(xc < xl & xc <= xr & xc < mu(i) - 1.75*sg(i)) ]);
dets = merge_chain(abs(xp - mu), cand, w1);
end

function d = elgendi(x, fs)
% our implementation of Elgendi, PLoS ONE 8(9):e73557, 2013
x = double(x(:));
[bb, aa] = butter(3, [8 20]/(fs/2));  y = filtfilt(bb, aa, x);
z = y.^2;
w1 = round(0.097*fs);  w2 = round(0.611*fs);
maq = movmean(z, w1);  mab = movmean(z, w2);
blk = maq > mab + 0.08*mean(z);
e = diff([0; blk; 0]);  on = find(e == 1);  off = find(e == -1) - 1;
d = zeros(numel(on), 1);  n = 0;
for i = 1:numel(on)
    if off(i) - on(i) + 1 >= w1
        [~, k] = max(abs(y(on(i):off(i))));  n = n + 1;  d(n) = on(i) + k - 1;
    end
end
d = d(1:n);
end

function m = noisy_mask(N, fs)
t = (0:N-1)' / fs / 60;
m = t >= 5 & mod(floor((t - 5) / 2), 2) == 0;
end

function [tp, fp, fn] = masked_counts(d, gt, tol, mask)
N = numel(mask);  gt = gt(:);  d = round(d(:));
gi = gt(gt >= 0 & gt < N);  g = gi(mask(gi + 1));
dd = d(d >= 0 & d < N);  keep = mask(dd + 1);
for i = find(~keep)'
    keep(i) = any(abs(g - dd(i)) <= tol);
end
[tp, fp, fn] = score_j(dd(keep), g, tol);
end

function y = tma0(x, w)
% zero-phase TMA: the second pass is shifted so that even K is also centred
x = double(x(:));
K = ceil(max(3, round(w)) / 2);
pad = min(2*K, numel(x)-1);
if pad < 1, y = x; return; end
xp = [flipud(x(1:pad)); x; flipud(x(end-pad+1:end))];
y0 = runmean_o(runmean_o(xp, K, floor(K/2)), K, ceil(K/2)-1);
y  = y0(pad+1 : pad+numel(x));
end

function y = runmean_o(v, K, o)
v = v(:);  M = numel(v);
cs = zeros(M + 2*K + 2, 1);
cs(K+2 : K+M+1) = cumsum(v);
cs(K+M+2 : end) = cs(K+M+1);
idx = (1:M)' + K + 1;
y = (cs(idx + o) - cs(idx + o - K)) / K;
end

function dets = detect(x, fs, kind)
x = double(x(:));  N = numel(x);
W  = estimate_width(x, fs);
W  = min(max(3, round(W)), max(3, round(2*fs)));
w1 = min(max(3, round(3.3*W)),  max(3, round(2*fs)));
w2 = min(max(3, round(0.06*W)), max(3, round(2*fs)));
switch kind
    case 1, xp = tma(x - tma(x, w1), w2);
    case 2, xp = x - tma(x, w1);
    case 3, g = ricker_fit(w1, w2);  xp = conv(x, g, 'same');
    case 4
        K1 = ceil(w1/2);  K2 = ceil(w2/2);  s1 = K1/sqrt(6);  s2 = max(K2/sqrt(6), 0.5);
        L = ceil(4*s1);  t = (-L:L)';
        g2 = exp(-t.^2/(2*s2^2));  g2 = g2/sum(g2);
        g1 = exp(-t.^2/(2*s1^2));  g1 = g1/sum(g1);
        xp = conv(x, g2 - g1, 'same');
    case 5
        [bb, aa] = butter(2, [5 15]/(fs/2));  xp = filtfilt(bb, aa, x);
    case 6   % same cascade, both TMAs exactly centred (zero phase)
        xp = tma0(x - tma0(x, w1), w2);
    case 7   % exact DoG target of Proposition 1: variance-matched, summed variance
        K1 = ceil(max(3,round(w1))/2);  K2 = ceil(max(3,round(w2))/2);
        v1 = (K1^2-1)/6;  v2 = (K2^2-1)/6;
        L = ceil(4*sqrt(v1+v2));  t = (-L:L)';
        g2 = exp(-t.^2/(2*v2));  g2 = g2/sum(g2);
        g1 = exp(-t.^2/(2*(v1+v2)));  g1 = g1/sum(g1);
        xp = conv(x, g2 - g1, 'same');
end
edge = min(max(round(0.01*fs), 2*w1), floor((N-3)/2));
z = xp(edge+1 : N-edge);  mu = median(z);  sigma = max(std(z), eps);
i  = (edge+1 : N-edge)';
xc = xp(i);  xl = xp(i-1);  xr = xp(i+1);
cand = sort([ i(xc > xl & xc >= xr & xc > mu + 1.75*sigma) ; ...
              i(xc < xl & xc <= xr & xc < mu - 1.75*sigma) ]);
dets = merge_chain(abs(xp), cand, w1);
end

function g = ricker_fit(w1, w2)
% Ricker wavelet whose scale maximises correlation with the TMA-diff kernel
L = 8*w1 + 1;  dlt = zeros(L,1);  mid = (L+1)/2;  dlt(mid) = 1;
h = tma(dlt - tma(dlt, w1), w2);  h = h - mean(h);
t = (1:L)' - mid;  best = -inf;  sb = 1;
for s = linspace(0.2*w2 + 0.5, 1.5*w1, 400)
    r = (1 - (t/s).^2) .* exp(-t.^2/(2*s^2));  r = r - mean(r);
    c = abs(h'*r) / (norm(h)*norm(r));
    if c > best, best = c;  sb = s; end
end
Lg = ceil(5*sb);  t = (-Lg:Lg)';
g = (1 - (t/sb).^2) .* exp(-t.^2/(2*sb^2));  g = g - mean(g);
end

function [tp, fp, fn, dt] = score_j(dets, gt, tol)
gt = sort(gt(:));  dets = sort(round(double(dets(:))));
used = false(numel(gt), 1);  tp = 0;  p = 1;  ng = numel(gt);  dt = zeros(numel(dets),1);
for i = 1:numel(dets)
    while p <= ng && gt(p) < dets(i) - tol, p = p + 1; end
    j = p;  bj = 0;  bd = inf;
    while j <= ng && gt(j) <= dets(i) + tol
        if ~used(j), dj = abs(gt(j) - dets(i)); if dj < bd, bd = dj; bj = j; end, end
        j = j + 1;
    end
    if bj > 0, used(bj) = true;  tp = tp + 1;  dt(tp) = bd; end
end
dt = dt(1:tp);  fp = numel(dets) - tp;  fn = ng - tp;
end

function d = pt_orig(x, fs)
% pan_tompkin.m (Sedghamiz, File Exchange 45840): 2nd output = R indices
try
    [~, d] = pan_tompkin(x, fs, 0);
catch err
    error('pan_tompkin(x, fs, 0) failed: %s', err.message);
end
d = round(d(:));
end

function out = merge_chain(a, d, radius)
d = sort(d(:));  out = [];
if isempty(d), return; end
g0 = 1;
for i = 2:numel(d)+1
    if i > numel(d) || d(i) - d(i-1) > radius
        g = d(g0:i-1);  [~, b] = max(a(g));  out(end+1,1) = g(b);   %#ok<AGROW>
        g0 = i;
    end
end
end

function p = sign_test(nb, nw)
% exact two-sided sign test (ties dropped)
n = nb + nw;
if n == 0, p = 1; return; end
k = min(nb, nw);  pk = 0;
for i = 0:k
    pk = pk + exp(gammaln(n+1) - gammaln(i+1) - gammaln(n-i+1) - n*log(2));
end
p = min(1, 2*pk);
end

function y = tma(x, w)
x = double(x(:));
K = ceil(max(3, round(w)) / 2);
pad = min(2*K, numel(x)-1);
if pad < 1, y = x; return; end
xp = [flipud(x(1:pad)); x; flipud(x(end-pad+1:end))];
y0 = runmean(runmean(xp, K), K);
y  = y0(pad+1 : pad+numel(x));
end

function y = runmean(v, K)
v = v(:);  M = numel(v);
o  = floor(K/2);
cs = zeros(M + 2*K + 2, 1);
cs(K+2 : K+M+1) = cumsum(v);
cs(K+M+2 : end) = cs(K+M+1);
idx = (1:M)' + K + 1;
y = (cs(idx + o) - cs(idx + o - K)) / K;
end

function W = estimate_width(x, fs)
x = double(x(:));  N = numel(x);
w1 = max(3, min(round(fs), max(3, round(N/10))));
xp = tma(x - tma(x, w1), max(3, round(0.05*fs)));
med = median(xp);
sig = max([1.4826*median(abs(xp-med)), std(xp), eps]);
d = max(3, round(0.1*fs));
c = local_peaks( xp,  med + 1.5*sig, d);
if numel(c) < 3
    c = local_peaks(-xp, -med + 1.5*sig, d);
end
if isempty(c), W = max(3, round(0.1*fs)); return; end
a = abs(xp);  widths = zeros(1, numel(c));
for i = 1:numel(c)
    h = 0.5*a(c(i));  l = c(i);  r = c(i);
    while l > 1 && a(l) > h, l = l - 1; end
    while r < N && a(r) > h, r = r + 1; end
    widths(i) = r - l;
end
widths = widths(widths >= 3 & widths <= 2*fs);
if isempty(widths), W = round(0.1*fs); else, W = round(median(widths)); end
W = max(3, min(W, round(2*fs)));
end

function idx = local_peaks(y, minheight, mindist)
y = y(:);
cand = find(y(2:end-1) > y(1:end-2) & y(2:end-1) >= y(3:end)) + 1;
cand = cand(y(cand) >= minheight);
if isempty(cand), idx = []; return; end
[~, ord] = sort(y(cand), 'descend');
cand = cand(ord);
keep = zeros(numel(cand), 1);  k = 0;
for i = 1:numel(cand)
    if k == 0 || all(abs(cand(i) - keep(1:k)) >= mindist)
        k = k + 1;  keep(k) = cand(i);
    end
end
idx = sort(keep(1:k))';
end


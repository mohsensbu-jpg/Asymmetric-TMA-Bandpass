function bspc_runs2()
%BSPC_RUNS2  Second set of checks for the BSPC paper (filter/detector: M. Karimi).
%   Part D  decision-stage ablation: SAME filter (D-TMA-I) and width estimate,
%           global threshold vs prespecified local threshold (mu, sigma over a
%           sliding 10-s window, alpha 1.75 - the V5r-loc rule fixed earlier).
%   Part E  close moving-average baseline: our implementation of Elgendi (2013)
%           (Butterworth 8-20 Hz order 3 zero phase, squaring, MA 97 ms and
%           611 ms, beta 0.08, blocks >= 97 ms, R = max |filtered| in block).
%   Part F  timing over all 48 MIT-BIH records (median of 3 runs per record).
%   Databases: MIT-BIH, SVDB (gross, +-100 / +-50 ms), NSTDB noisy segments.
root = userpath;  if isempty(root), root = pwd; end
if exist('pan_tompkin', 'file') ~= 2
    f = dir(fullfile(root, '**', 'pan_tompkin.m'));  addpath(f(1).folder);
end
DB = struct('name', {'MIT-BIH','SVDB'}, 'folder', {fullfile(root,'mitdb_full'), fullfile(root,'svdb')}, ...
    'fs', {360, 128}, 'recs', {[100:109 111:119 121:124 200:203 205 207:210 212:215 217 219:223 228 230:234], ...
    [800:812 820:829 840:894]});
NM = {'TMA-I global','TMA-I local','Elgendi'};
tic;
fprintf('\nPARTS D+E  gross results (Se, +P, F1 at 100 ms; F1 at 50 ms; FP; FN)\n');
for b = 1:2
    fs = DB(b).fs;  recs = DB(b).recs;  m = numel(recs);
    C = zeros(m, 3, 2, 3);  Fr = zeros(m, 3);
    for k = 1:m
        [x, gt] = read_mitdb(DB(b).folder, recs(k));  x = x - median(x);
        for j = 1:3
            d = run_method(j, x, fs);
            for t = 1:2
                tol = round([0.100 0.050]*fs);
                [tp, fp, fn] = score_j(d, gt, tol(t));  C(k,j,t,:) = [tp fp fn];
                if t == 1, Fr(k,j) = 2*tp/max(2*tp+fp+fn,1); end
            end
        end
    end
    fprintf('%s\n  method          Se       +P       F1(100)  F1(50)   FP     FN    rec. better/worse vs TMA-I global\n', DB(b).name);
    for j = 1:3
        c = squeeze(sum(C(:,j,1,:),1));  c2 = squeeze(sum(C(:,j,2,:),1));
        s = '-';
        if j > 1
            dd = Fr(:,j) - Fr(:,1);  nb = sum(dd > 1e-9);  nw = sum(dd < -1e-9);
            s = sprintf('%d/%d p=%.3f', nb, nw, sign_test(nb, nw));
        end
        fprintf('  %-14s  %.5f  %.5f  %.5f  %.5f  %5d  %5d   %s\n', NM{j}, c(1)/(c(1)+c(3)), ...
            c(1)/(c(1)+c(2)), 2*c(1)/(2*c(1)+c(2)+c(3)), 2*c2(1)/(2*c2(1)+c2(2)+c2(3)), c(2), c(3), s);
    end
    fprintf('  per-record median F1: %s\n', sprintf('%.4f  ', median(Fr)));
end
fN = fullfile(root, 'nstdb');  fs = 360;  tol = round(0.100*fs);
snr = {'24','18','12','06','00','_6'};  lab = [24 18 12 6 0 -6];
fprintf('\nNSTDB noisy segments, 118+119 pooled, +-100 ms (W, mu, sigma from the whole record)\n');
fprintf('  SNR'); fprintf('%16s', NM{:}); fprintf('\n');
for s = 1:numel(snr)
    C = zeros(3, 3);
    for r = {'118','119'}
        [x, gt] = read_mitdb(fN, [r{1} 'e' snr{s}]);  x = x - median(x);
        mask = noisy_mask(numel(x), fs);
        for j = 1:3
            d = run_method(j, x, fs);
            [tp, fp, fn] = masked_counts(d, gt, tol, mask);  C(j,:) = C(j,:) + [tp fp fn];
        end
    end
    fprintf('  %3d', lab(s)); fprintf('%16.4f', 2*C(:,1)./(2*C(:,1)+C(:,2)+C(:,3))); fprintf('\n');
    fprintf('   FP'); fprintf('%16d', C(:,2)); fprintf('\n');
    fprintf('   FN'); fprintf('%16d', C(:,3)); fprintf('   (beats %d)\n', C(1,1)+C(1,3));
end
%% Part F timing
fprintf('\nPART F  timing, 48 MIT-BIH records, median of 3 runs each; ns per sample\n');
recs = DB(1).recs;  T = zeros(numel(recs), 6);
for k = 1:numel(recs)
    x = read_mitdb(DB(1).folder, recs(k));  x = x - median(x);  N = numel(x);  fs = 360;
    W = estimate_width(x, fs);  W = min(max(3, round(W)), max(3, round(2*fs)));
    w1 = max(3, round(3.3*W));  w2 = max(3, round(0.06*W));
    K1 = ceil(w1/2);  K2 = ceil(w2/2);  v1 = (K1^2-1)/6;  v2 = (K2^2-1)/6;
    L = ceil(4*sqrt(v1+v2));  t = (-L:L)';
    g = exp(-t.^2/(2*v2))/sum(exp(-t.^2/(2*v2))) - exp(-t.^2/(2*(v1+v2)))/sum(exp(-t.^2/(2*(v1+v2))));
    fns = {@() tma(x - tma(x, w1), w2), @() conv(x, g, 'same'), ...
           @() run_method(1, x, fs), @() run_method(2, x, fs), @() run_method(3, x, fs), @() pt_orig(x, fs)};
    for j = 1:6
        tt = zeros(3,1);
        for rep = 1:3, t0 = tic; fns{j}(); tt(rep) = toc(t0); end
        T(k,j) = 1e9*median(tt)/N;
    end
end
lb = {'filter TMA-I','filter DoG conv','detector global','detector local','Elgendi','Pan-Tompkins'};
for j = 1:6
    fprintf('  %-16s median %7.1f ns/sample  (IQR %6.1f-%6.1f)\n', lb{j}, median(T(:,j)), qtl(T(:,j),0.25), qtl(T(:,j),0.75));
end
fprintf('  DoG taps: %d (record %d)\n', numel(g), recs(end));
fprintf('  computer %s, MATLAB %s, total time %.0f s\n\n', computer, version, toc);
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


function bspc_runs()
%BSPC_RUNS  Extra experiments for the BSPC version of the TMA-diff paper
%   (filter and detector: M. Karimi). Text output only; runs in MATLAB Online.
%   Part A  MIT-BIH (48) and SVDB (78): published TMA-diff vs the zero-phase
%           (exactly centred) cascade and vs the exact DoG target of
%           Proposition 1. Same threshold/search/merge for all.
%   Part B  NSTDB (118, 119 at 24..-6 dB): F1 on the NOISY segments only,
%           all filters + Pan-Tompkins, +-100 ms.
%   Needs: read_mitdb.m, pan_tompkin.m, folders mitdb_full, svdb, nstdb.
root = userpath;  if isempty(root), root = pwd; end
if exist('pan_tompkin', 'file') ~= 2
    f = dir(fullfile(root, '**', 'pan_tompkin.m'));  addpath(f(1).folder);
end
NAMES = {'TMA-diff','HP-only','Ricker','DoG-abl','BP 5-15','TMA-zero','DoG-exact','Pan-Tompkins'};
tic;
%% Part A
DB = struct('name', {'MIT-BIH','SVDB'}, 'folder', {fullfile(root,'mitdb_full'), fullfile(root,'svdb')}, ...
    'fs', {360, 128}, 'recs', {[100:109 111:119 121:124 200:203 205 207:210 212:215 217 219:223 228 230:234], ...
    [800:812 820:829 840:894]});
JA = [1 6 7];
fprintf('\nPART C  estimated width W and window halves K1, K2 (samples)\n');
for b = 1:2
    fs = DB(b).fs;  Wv = zeros(numel(DB(b).recs),1);
    for k = 1:numel(DB(b).recs)
        x = read_mitdb(DB(b).folder, DB(b).recs(k));  x = x - median(x);
        W = estimate_width(double(x(:)), fs);  Wv(k) = min(max(3, round(W)), max(3, round(2*fs)));
    end
    K1 = ceil(max(3,round(3.3*Wv))/2);  K2 = ceil(max(3,round(0.06*Wv))/2);
    fprintf('  %s: W median %d (min %d, max %d) = %.1f ms; K1 median %d (%d-%d); K2 values: %s\n', ...
        DB(b).name, round(median(Wv)), min(Wv), max(Wv), 1000*median(Wv)/fs, ...
        round(median(K1)), min(K1), max(K1), mat2str(unique(K2)'));
    fprintf('     records with even K1: %d of %d\n', sum(mod(K1,2)==0), numel(K1));
    for r2 = [0.06 0.10 0.15 0.30 0.50 0.70]
        K2r = ceil(max(3,round(r2*Wv))/2);
        fprintf('     w2 = %.2fW -> K2 median %d (%d-%d)\n', r2, round(median(K2r)), min(K2r), max(K2r));
    end
end
fprintf('\nPART A  (gross F1; record wins/losses vs TMA-diff at 100 ms, exact sign test)\n');
for b = 1:2
    fs = DB(b).fs;  recs = DB(b).recs;  m = numel(recs);
    C = zeros(m, numel(JA), 2, 3);  Fr = zeros(m, numel(JA));
    for k = 1:m
        [x, gt] = read_mitdb(DB(b).folder, recs(k));  x = x - median(x);
        for jj = 1:numel(JA)
            d = detect(x, fs, JA(jj));
            for t = 1:2
                tol = round([0.100 0.050]*fs);
                [tp, fp, fn] = score_j(d, gt, tol(t));  C(k,jj,t,:) = [tp fp fn];
                if t == 1, Fr(k,jj) = 2*tp/max(2*tp+fp+fn,1); end
            end
        end
    end
    fprintf('%s\n  filter        F1(100)   F1(50)    FP     FN    wins/losses\n', DB(b).name);
    for jj = 1:numel(JA)
        c1 = squeeze(sum(C(:,jj,1,:),1));  c2 = squeeze(sum(C(:,jj,2,:),1));
        s = '-';
        if jj > 1
            dd = Fr(:,jj) - Fr(:,1);  nb = sum(dd > 1e-9);  nw = sum(dd < -1e-9);
            s = sprintf('%d/%d p=%.3f', nb, nw, sign_test(nb, nw));
        end
        fprintf('  %-11s  %.5f  %.5f  %5d  %5d   %s\n', NAMES{JA(jj)}, ...
            2*c1(1)/(2*c1(1)+c1(2)+c1(3)), 2*c2(1)/(2*c2(1)+c2(2)+c2(3)), c1(2), c1(3), s);
    end
end
%% Part B
fN = fullfile(root, 'nstdb');  fs = 360;  tol = round(0.100*fs);
snr = {'24','18','12','06','00','_6'};  lab = [24 18 12 6 0 -6];
fprintf('\nPART B  NSTDB, noisy segments only, records 118+119 pooled, +-100 ms\n');
fprintf('  SNR ');  fprintf('%13s', NAMES{:});  fprintf('\n');
for s = 1:numel(snr)
    C = zeros(numel(NAMES), 3);
    for r = {'118','119'}
        [x, gt] = read_mitdb(fN, [r{1} 'e' snr{s}]);  x = x - median(x);
        mask = noisy_mask(numel(x), fs);
        for j = 1:numel(NAMES)
            if j == numel(NAMES), d = pt_orig(x, fs); else, d = detect(x, fs, j); end
            [tp, fp, fn] = masked_counts(d, gt, tol, mask);  C(j,:) = C(j,:) + [tp fp fn];
        end
    end
    fprintf('  %3d ', lab(s));
    fprintf('%13.4f', 2*C(:,1)./(2*C(:,1)+C(:,2)+C(:,3)));  fprintf('\n');
    fprintf('   FP');  fprintf('%13d', C(:,2));  fprintf('\n');
    fprintf('   FN');  fprintf('%13d', C(:,3));  fprintf('\n');
end
fprintf('\n  time %.0f s,  MATLAB %s\n\n', toc, version);
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


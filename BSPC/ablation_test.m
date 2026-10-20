function ablation_test()
%ABLATION_TEST  Same detector, only the FILTER changes (stages 2-4 identical:
%   global threshold mu +- 1.75 sigma, dual-polarity search, chain merge at
%   radius w1 = 3.3W, same automatic width estimate W). Text output only.
%     TMA-diff  proposed: TMA(x - TMA(x, w1), w2), w2 = 0.06W   (M. Karimi)
%     HP-only   x - TMA(x, w1)            (baseline removal, no low-pass)
%     Ricker    convolution with a Ricker wavelet, scale fitted to TMA-diff
%     DoG       exact Gaussian DoG with the same scales (sigma_i = K_i/sqrt 6)
%     BP 5-15   2nd-order Butterworth 5-15 Hz, zero phase (filtfilt)
%   Reference: Pan-Tompkins (Sedghamiz) on the same data.
%   Tolerances +-100 ms and +-50 ms; jitter = median |detection - beat|.
%   Databases: MIT-BIH full (48) and SVDB (78).

root = userpath;  if isempty(root), root = pwd; end
if exist('pan_tompkin', 'file') ~= 2
    f = dir(fullfile(root, '**', 'pan_tompkin.m'));  addpath(f(1).folder);
end
DB = struct('name', {'MIT-BIH', 'SVDB'}, 'folder', {fullfile(root,'mitdb_full'), fullfile(root,'svdb')}, ...
    'fs', {360, 128}, 'recs', {[100:109 111:119 121:124 200:203 205 207:210 212:215 217 219:223 228 230:234], ...
    [800:812 820:829 840:894]});
NAMES = {'TMA-diff', 'HP-only', 'Ricker', 'DoG', 'BP 5-15', 'Pan-Tompkins'};  nF = numel(NAMES);
tic;
for b = 1:2
    fs = DB(b).fs;  recs = DB(b).recs;  m = numel(recs);
    C = zeros(m, nF, 2, 3);  Fr = zeros(m, nF, 2);  J = cell(nF, 1);
    for k = 1:m
        [x, gt] = read_mitdb(DB(b).folder, recs(k));  x = x - median(x);
        for j = 1:nF
            if j == nF, d = pt_orig(x, fs); else, d = detect(x, fs, j); end
            for t = 1:2
                tol = round([0.100 0.050]*fs);  tol = tol(t);
                [tp, fp, fn, dt] = score_j(d, gt, tol);
                C(k,j,t,:) = [tp fp fn];  Fr(k,j,t) = 2*tp / max(2*tp+fp+fn, 1);
                if t == 1, J{j} = [J{j}; dt]; end
            end
        end
    end
    fprintf('\n%s (%d records, %d beats, %d Hz)\n', DB(b).name, m, sum(C(:,1,1,1)+C(:,1,1,3)), fs);
    fprintf('  filter         F1 +-100ms  F1 +-50ms   FP(100)  FN(100)  jitter(ms)  records better/worse than TMA-diff (100 ms)\n');
    for j = 1:nF
        c1 = squeeze(sum(C(:,j,1,:),1));  c2 = squeeze(sum(C(:,j,2,:),1));
        f1a = 2*c1(1)/(2*c1(1)+c1(2)+c1(3));  f1b = 2*c2(1)/(2*c2(1)+c2(2)+c2(3));
        if j == 1, cmp = '-';
        else
            dd = Fr(:,j,1) - Fr(:,1,1);  nb = sum(dd > 1e-9);  nw = sum(dd < -1e-9);
            cmp = sprintf('%d / %d (p = %.3f)', nb, nw, sign_test(nb, nw));
        end
        fprintf('  %-12s   %.5f     %.5f    %6d   %6d   %6.1f      %s\n', NAMES{j}, f1a, f1b, ...
                c1(2), c1(3), 1000*median(J{j})/fs, cmp);
    end
end
fprintf('\n  time %.0f s,  MATLAB %s\n\n', toc, version);
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


function bspc_timing_table()
%BSPC_TIMING_TABLE  Localisation on matched beats only (+-100 ms matching),
%   frozen methods, no change: n matched, median signed error, median and 95th
%   percentile absolute error, and % of matched beats with |error| < 20 / 50 ms.
%   Sign: positive = detection AFTER the annotation.
root = userpath;  if isempty(root), root = pwd; end
if exist('pan_tompkin', 'file') ~= 2
    f = dir(fullfile(root, '**', 'pan_tompkin.m'));  addpath(f(1).folder);
end
S = [10 1 1.75 0 1];
NM = {'L proposed', 'G D-TMA-I', 'Elgendi', 'Pan-Tompkins'};
fprintf('\n  method          n matched   med signed  med |e|  p95 |e|   <20 ms   <50 ms\n');
for b = 1:3
    switch b
        case 1, recs = num2cell([100:109 111:119 121:124 200:203 205 207:210 212:215 217 219:223 228 230:234]); fs = 360; nm = 'MIT-BIH';
        case 2, recs = num2cell([800:812 820:829 840:894]); fs = 128; nm = 'SVDB';
        case 3, recs = num2cell(1:75); fs = 257; nm = 'INCART';
    end
    E = cell(4,1);
    for k = 1:numel(recs)
        switch b
            case 1, [x, gt] = read_mitdb(fullfile(root, 'mitdb_full'), recs{k});
            case 2, [x, gt] = read_mitdb(fullfile(root, 'svdb'), recs{k});
            case 3, Sx = load(fullfile(root, 'incart', sprintf('I%02d.mat', recs{k})));  x = double(Sx.x);  gt = Sx.beats;
        end
        x = x - median(x);
        for j = 1:4
            switch j
                case 1, Fx = front(x, fs);  d = decide(Fx, dstats(Fx, S(1), S(2), fs), S, fs);
                case 2, d = detect(x, fs, 1);
                case 3, d = elgendi(x, fs);
                case 4, d = pt_orig(x, fs);
            end
            [~, sd] = match_flags(d - 1, gt, round(0.1*fs));
            E{j} = [E{j}; 1000*sd/fs];
        end
    end
    fprintf('%s (1 sample = %.1f ms)\n', nm, 1000/fs);
    for j = 1:4
        e = E{j};  a = abs(e);
        fprintf('  %-13s  %9d   %+8.1f   %6.1f   %6.1f   %6.2f%%  %6.2f%%\n', NM{j}, numel(e), median(e), median(a), qv(a, 0.95), 100*mean(a < 20), 100*mean(a < 50));
    end
end
fprintf('\n');
end

function o = weoptions_()
o = weboptions('Timeout', 60);
end

function v = qv(x, q)
x = sort(x(:));  v = zeros(size(q));
if isempty(x), return; end
for i = 1:numel(q), v(i) = x(max(1, min(numel(x), round(q(i)*numel(x))))); end
end

function c = aami(t)
c = zeros(size(t));
c(ismember(t, [1 2 3 11 25 34])) = 1;      % N, L, R, j, B, e
c(ismember(t, [4 7 8 9 35])) = 2;          % a, J, A, S, n
c(ismember(t, [5 10 30])) = 3;             % V, E, r
c(t == 6) = 4;                             % F
c(ismember(t, [12 13 38 41])) = 5;         % /, Q, f, other
end

function [hit, sd] = match_flags(d, gt, tol)
% same greedy one-to-one matching as score_j; hit(i) true if beat i matched;
% sd = signed error (detection - beat, samples) of matched pairs
gt = gt(:);  [gs, og] = sort(gt);  d = sort(round(double(d(:))));
used = false(numel(gs),1);  p = 1;  ng = numel(gs);  sd = zeros(numel(d),1);  n = 0;
for i = 1:numel(d)
    while p <= ng && gs(p) < d(i) - tol, p = p + 1; end
    j = p;  bj = 0;  bd = inf;
    while j <= ng && gs(j) <= d(i) + tol
        if ~used(j), dj = abs(gs(j) - d(i)); if dj < bd, bd = dj; bj = j; end, end
        j = j + 1;
    end
    if bj > 0, used(bj) = true;  n = n + 1;  sd(n) = d(i) - gs(bj); end
end
hit = false(numel(gt),1);  hit(og(used)) = true;  sd = sd(1:n);
end

function [pos, typ] = read_atr_types(file, N)
fid = fopen(file, 'r');  a = fread(fid, inf, 'uint8=>double');  fclose(fid);
w = a(1:2:end-1) + 256*a(2:2:end);
t = 0;  k = 1;  pos = zeros(numel(w),1);  typ = zeros(numel(w),1);  c = 0;
while k <= numel(w)
    A = floor(w(k)/1024);  I = mod(w(k), 1024);
    if A == 0 && I == 0, break; end
    if A == 59
        v = w(k+1)*65536 + w(k+2);  if v >= 2^31, v = v - 2^32; end
        t = t + v;  k = k + 3;  continue
    elseif A == 63, k = k + 1 + ceil(I/2);  continue
    elseif A >= 60, k = k + 1;  continue
    end
    t = t + I;  c = c + 1;  pos(c) = t;  typ(c) = A;  k = k + 1;
end
pos = pos(1:c);  typ = typ(1:c);
keep = ismember(typ, [1:13 25 30 34 35 38 41]) & pos >= 0 & pos < N;
pos = pos(keep);  typ = typ(keep);
end

function [d, dM, c] = decide_trace(F, P, cfg, fs)
% identical to decide(), but also returns the candidates c and the
% detections after the chain merge dM (before the RR/T rule)
al = cfg(3);  rr = cfg(5);
N = numel(F.xp);  e = P.e;  s = P.s;
i = (F.edge+1 : N-F.edge)';
c = i(e(i) > e(i-1) & e(i) >= e(i+1) & e(i) > al*s(i));
a = e ./ s;
dM = merge_chain(a, c, F.w1);  d = dM;
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


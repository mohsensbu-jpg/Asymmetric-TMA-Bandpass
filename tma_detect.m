function [dets, xp] = tma_detect(x, fs, varargin)
%TMA_DETECT  O(N) transient detector using a difference of triangular
%            moving averages.
%
%   DETS = TMA_DETECT(X, FS) returns the sample indices of detected events
%   in signal X sampled at FS hertz.
%
%   [DETS, XP] = TMA_DETECT(...) also returns the band-pass filtered signal.
%
%   Name-value options:
%     'W'      nominal event width in samples. Estimated from X if omitted.
%     'alpha'  threshold multiplier (default 1.75).
%     'w1'     wide window   (default 3.3*W)
%     'w2'     narrow window (default 0.06*W)
%
%   The four stages correspond directly to Section IV of the paper:
%
%     1. filter     xp = TMA(x - TMA(x, w1), w2)
%     2. threshold  mu +/- alpha*sigma over the record interior
%     3. search     local maxima above and local minima below
%     4. merge      candidates within w1 form one event; keep max |xp|
%
%   Each TMA is a pair of cascaded running-sum moving averages, so the cost
%   is one addition and one subtraction per output sample, independent of
%   window length, and the working memory is a single accumulator rather
%   than a window buffer.
%
%   Example:
%     T = readtable('100.csv');
%     x = T.signal_value - median(T.signal_value);
%     dets = tma_detect(x, 360);
%
%   See also REPRODUCE_TABLES.

p = inputParser;
addParameter(p, 'W',     []);
addParameter(p, 'alpha', 1.75);
addParameter(p, 'w1',    []);
addParameter(p, 'w2',    []);
parse(p, varargin{:});

x = double(x(:));
N = numel(x);
alpha = p.Results.alpha;

W = p.Results.W;
if isempty(W), W = estimate_width(x, fs); end
W = min(max(3, round(W)), max(3, round(2*fs)));

w1 = p.Results.w1;  if isempty(w1), w1 = 3.3*W;  end
w2 = p.Results.w2;  if isempty(w2), w2 = 0.06*W; end
w1 = min(max(3, round(w1)), max(3, round(2*fs)));
w2 = min(max(3, round(w2)), max(3, round(2*fs)));

% --- stage 1: filter -----------------------------------------------------
xp = tma(x - tma(x, w1), w2);

% --- stage 2: threshold --------------------------------------------------
edge = min(max(round(0.01*fs), 2*w1), floor((N-3)/2));
z     = xp(edge+1 : N-edge);
mu    = median(z);
sigma = max(std(z), eps);
thr_hi = mu + alpha*sigma;
thr_lo = mu - alpha*sigma;

% --- stage 3: dual-polarity extremum search ------------------------------
i  = (edge+1 : N-edge)';
xc = xp(i);  xl = xp(i-1);  xr = xp(i+1);
cand = sort([ i(xc > xl & xc >= xr & xc > thr_hi) ; ...
              i(xc < xl & xc <= xr & xc < thr_lo) ]);

% --- stage 4: chain merge ------------------------------------------------
if isempty(cand)
    dets = [];
else
    dets = chain_merge(abs(xp), cand', w1);
end
end

% =========================================================================
function out = chain_merge(a, dets, radius)
%CHAIN_MERGE  Consecutive detections separated by <= RADIUS form one group;
%             the member with the largest A is kept.
dets = sort(dets(:))';
out = zeros(1, numel(dets));
k = 0;  gstart = 1;
for i = 2 : numel(dets)+1
    if i > numel(dets) || dets(i) - dets(i-1) > radius
        g = dets(gstart : i-1);
        [~, b] = max(a(g));
        k = k + 1;  out(k) = g(b);
        gstart = i;
    end
end
out = out(1:k);
end

% =========================================================================
function y = tma(x, w)
%TMA  Triangular moving average: two cascaded running-sum moving averages
%     of length K = ceil(w/2), with reflective padding at the boundaries.
x = double(x(:));
K = ceil(max(3, round(w)) / 2);
pad = min(2*K, numel(x)-1);
if pad < 1, y = x; return; end
xp = [flipud(x(1:pad)); x; flipud(x(end-pad+1:end))];
y0 = runmean(runmean(xp, K), K);
y  = y0(pad+1 : pad+numel(x));
end

function y = runmean(v, K)
%RUNMEAN  Moving average of length K by cumulative sums: one subtraction
%         per output sample, independent of K. Matches
%         conv(v, ones(K,1)/K, 'same') exactly.
v = v(:);  M = numel(v);
o  = floor(K/2);
cs = zeros(M + 2*K + 2, 1);
cs(K+2 : K+M+1) = cumsum(v);
cs(K+M+2 : end) = cs(K+M+1);
idx = (1:M)' + K + 1;
y = (cs(idx + o) - cs(idx + o - K)) / K;
end

% =========================================================================
function W = estimate_width(x, fs)
%ESTIMATE_WIDTH  Median half-amplitude width of prominent events, used to
%                set w1 and w2 when W is not supplied.
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

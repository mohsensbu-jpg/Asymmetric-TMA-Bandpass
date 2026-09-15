function reproduce_tables(datadir)
%REPRODUCE_TABLES  Regenerates the numbers reported in the paper.
%
%   REPRODUCE_TABLES(DATADIR) expects, for each of the 48 MIT-BIH records,
%   a signal file <rec>.csv with columns sample_index, signal_value, and an
%   annotation file <rec>annotations.txt with columns sample_index, type.
%
%   Prints per-record results and the summary figures of Section VI-B, and
%   writes results.csv.
%
%   To reproduce the comparison of Table III you also need pan_tompkin.m
%   from MATLAB File Exchange submission 45840; see README.

if nargin < 1 || isempty(datadir), datadir = pwd; end

fs   = 360;
TOL  = round(0.100*fs);          % +-100 ms
recs = [100:109 111:119 121:124 200:203 205 207:210 212:215 217 ...
        219:223 228 230:234];

n = numel(recs);
R = nan(n, 7);

fprintf('\n  Rec  beats   TP   FP   FN  Recall      F1\n');
fprintf('  ------------------------------------------\n');

for k = 1:n
    [x, gt] = load_record(datadir, recs(k));
    if isempty(x)
        fprintf('  %3d   file not found\n', recs(k));
        continue;
    end
    x = x - median(x);

    dets = tma_detect(x, fs);

    [tp, fp, fn, ~, rc, f1] = score(dets, gt, TOL);
    R(k,:) = [recs(k) numel(gt) tp fp fn rc f1];
    fprintf('  %3d  %5d %4d %4d %4d   %.3f  %.4f\n', R(k,:));
end

ok = ~isnan(R(:,1));
F1 = R(ok,7);
tp = sum(R(ok,3));  fp = sum(R(ok,4));  fn = sum(R(ok,5));
pr = tp/(tp+fp);    rc = tp/(tp+fn);

fprintf('  ------------------------------------------\n');
fprintf('  mean F1        %.4f\n', mean(F1));
fprintf('  median F1      %.4f\n', median(F1));
fprintf('  F1 >= 0.95     %d / %d\n', sum(F1>=0.95), numel(F1));
fprintf('  F1 <  0.80     %d / %d\n', sum(F1<0.80),  numel(F1));
fprintf('  totals         TP %d  FP %d  FN %d  of %d beats\n', ...
        tp, fp, fn, sum(R(ok,2)));
fprintf('  aggregate      precision %.4f  recall %.4f  F1 %.4f\n', ...
        pr, rc, 2*pr*rc/(pr+rc));
fprintf('\n  MATLAB %s on %s\n\n', version, computer);

writetable(array2table(R(ok,:), 'VariableNames', ...
    {'Rec','Beats','TP','FP','FN','Recall','F1'}), 'results.csv');
fprintf('  Written: results.csv\n\n');
end

% =========================================================================
function [x, gt] = load_record(datadir, r)
x = [];  gt = [];
fsig = fullfile(datadir, sprintf('%d.csv', r));
fann = fullfile(datadir, sprintf('%dannotations.txt', r));
if exist(fsig,'file') ~= 2 || exist(fann,'file') ~= 2, return; end
Ts = readtable(fsig);
x  = double(Ts.signal_value);
Ta = readtable(fann, 'Delimiter', ',');
typ = Ta.type;
if ~iscell(typ), typ = cellstr(string(typ)); end
gt = double(Ta.sample_index(~strcmp(strtrim(typ), '+')));
gt = gt(gt >= 0 & gt < numel(x));
end

function [tp, fp, fn, pr, rc, f1] = score(dets, gt, tol)
%SCORE  Greedy nearest matching; each annotated beat is claimed at most once.
gt   = sort(gt(:));
dets = sort(round(double(dets(:))));
used = false(numel(gt), 1);
tp = 0;
for i = 1:numel(dets)
    c = find(~used & abs(gt - dets(i)) <= tol);
    if ~isempty(c)
        [~, j] = min(abs(gt(c) - dets(i)));
        used(c(j)) = true;
        tp = tp + 1;
    end
end
fp = numel(dets) - tp;
fn = numel(gt) - tp;
pr = tp / max(tp + fp, 1);
rc = tp / max(tp + fn, 1);
if pr + rc > 0, f1 = 2*pr*rc/(pr+rc); else, f1 = 0; end
end

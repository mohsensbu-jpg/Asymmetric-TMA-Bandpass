function butqdb_quality()
%BUTQDB_QUALITY  Share of time in each expert quality class in real
%   free-living wearable ECG (BUT QDB, Bittium Faros 180, 1000 Hz).
%   Downloads only the small annotation files (*_ANN.csv). Class 1: full
%   waveform analysis possible; class 2: only QRS detection reliable; class 3:
%   unusable. Smital et al. (IEEE TBME 2020) relate these to SNR limits of
%   18 dB and 5 dB. Uses the consensus annotation (columns 10-12).
bases = {'https://physionet-open.s3.amazonaws.com/butqdb/1.0.0/', 'https://physionet.org/files/butqdb/1.0.0/'};
opt = weboptions('Timeout', 120);
recs = {};
for b = 1:2
    try, t = webread([bases{b} 'RECORDS'], opt); recs = strsplit(strtrim(char(t))); base = bases{b}; break; catch, end
end
if isempty(recs), error('could not read RECORDS'); end
ids = unique(regexprep(recs, '[/_].*$', ''));
tmp = tempname;  mkdir(tmp);
T = zeros(numel(ids), 3);
for k = 1:numel(ids)
    id = ids{k};  f = fullfile(tmp, [id '_ANN.csv']);  ok = false;
    for b = 1:2
        try, websave(f, [bases{b} id '/' id '_ANN.csv'], opt); ok = true; break; catch, end
    end
    if ~ok, fprintf('  %s: annotation not found\n', id); continue; end
    A = readmatrix(f);
    if size(A, 2) < 12, fprintf('  %s: unexpected format (%d columns)\n', id, size(A,2)); continue; end
    A = A(:, 10:12);  A = A(all(~isnan(A), 2), :);
    for c = 1:3
        sel = A(:,3) == c;  T(k,c) = sum(A(sel,2) - A(sel,1) + 1) / 1000 / 60;   % minutes
    end
    fprintf('  %s: annotated %.0f min; class 1 %.1f%%, class 2 %.1f%%, class 3 %.1f%%\n', id, sum(T(k,:)), 100*T(k,:)/max(sum(T(k,:)), eps));
end
tot = sum(T, 1);
fprintf('\nALL annotated (%.0f min, %.1f h): class 1 %.1f%%, class 2 %.1f%%, class 3 %.1f%%\n', sum(tot), sum(tot)/60, 100*tot/sum(tot));
full = sum(T, 2) > 20*60;   % records annotated over (almost) their whole length
if any(full)
    tf = sum(T(full,:), 1);
    fprintf('Fully annotated records only (%d records, %.1f h): class 1 %.1f%%, class 2 %.1f%%, class 3 %.1f%%\n', sum(full), sum(tf)/60, 100*tf/sum(tf));
end
fprintf('Note: the partially annotated segments, including extra poor-quality segments, were selected by the database authors and are not a random sample.\n\n');
rmdir(tmp, 's');
end

function get_mitdb_full()
%GET_MITDB_FULL  Download the 48 full-length MIT-BIH Arrhythmia records
%   (30 min each) from PhysioNet straight into MATLAB Drive.
%   The transfer is server-to-server (MathWorks <-> PhysioNet), so a slow
%   home connection does not matter. Safe to re-run: files already present
%   are skipped, so an interrupted run simply continues.
%
%   Files go to <MATLAB Drive>/mitdb_full/ : <rec>.hea, <rec>.dat, <rec>.atr
%   Total size about 100 MB.

recs = [100:109 111:119 121:124 200:203 205 207:210 212:215 217 ...
        219:223 228 230:234];
base = 'https://physionet.org/files/mitdb/1.0.0/';
root = userpath;  if isempty(root), root = pwd; end
out  = fullfile(root, 'mitdb_full');
if ~exist(out, 'dir'), mkdir(out); end
opt = weboptions('Timeout', 120);

fprintf('\nDownloading into %s\n', out);
nOK = 0;
for r = recs
    ok = true;
    for ext = {'.hea', '.dat', '.atr'}
        f = fullfile(out, sprintf('%d%s', r, ext{1}));
        if exist(f, 'file') == 2
            fi = dir(f);  if fi.bytes > 0, continue; end
        end
        done = false;
        for attempt = 1:3
            try
                websave(f, [base sprintf('%d%s', r, ext{1})], opt);
                done = true;  break
            catch err
                pause(2);
            end
        end
        if ~done
            ok = false;
            fprintf('  %d%s FAILED: %s\n', r, ext{1}, err.message);
        end
    end
    if ok, nOK = nOK + 1; fprintf('  %d ok\n', r); end
end
fprintf('%d of %d records complete.\n', nOK, numel(recs));

% sanity check on record 100 (known: 650000 samples, 2273 beat labels)
try
    [x, beats] = read_mitdb(out, 100);
    fprintf('Check record 100: %d samples (expect 650000), %d beats (expect 2273)\n', ...
            numel(x), numel(beats));
catch err
    fprintf('Check of record 100 failed: %s\n', err.message);
end
end

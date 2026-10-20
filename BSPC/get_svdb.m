function get_svdb()
%GET_SVDB  Download the MIT-BIH Supraventricular Arrhythmia Database
%   (78 records, 30 min, 128 Hz, 2 ECG leads) from PhysioNet into MATLAB Drive.
%   The transfer is server-to-server (MathWorks <-> PhysioNet), so a slow
%   home connection does not matter. Safe to re-run: files already present
%   are skipped, so an interrupted run simply continues.
%
%   Files go to <MATLAB Drive>/svdb/ : <rec>.hea, <rec>.dat, <rec>.atr
%   Total size about 55 MB.

recs = [800:812 820:829 840:894];
base = 'https://physionet.org/files/svdb/1.0.0/';
root = userpath;  if isempty(root), root = pwd; end
out  = fullfile(root, 'svdb');
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

try
    [x, beats] = read_mitdb(out, 800);
    fprintf('Check 800: %d samples (expect 230400), %d beats\n', numel(x), numel(beats));
catch err
    fprintf('Check of 800 failed: %s\n', err.message);
end
end

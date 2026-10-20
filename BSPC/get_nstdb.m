function get_nstdb()
%GET_NSTDB  Download the MIT-BIH Noise Stress Test Database (12 records,
%   30 min each: records 118 and 119 with real electrode-motion noise added
%   at SNR 24, 18, 12, 6, 0 and -6 dB) from PhysioNet into MATLAB Drive.
%   The transfer is server-to-server (MathWorks <-> PhysioNet), so a slow
%   home connection does not matter. Safe to re-run: files already present
%   are skipped, so an interrupted run simply continues.
%
%   Files go to <MATLAB Drive>/nstdb/ : <rec>.hea, <rec>.dat, <rec>.atr
%   Total size about 25 MB.

snr  = {'24', '18', '12', '06', '00', '_6'};
recs = [strcat('118e', snr), strcat('119e', snr)];
base = 'https://physionet.org/files/nstdb/1.0.0/';
root = userpath;  if isempty(root), root = pwd; end
out  = fullfile(root, 'nstdb');
if ~exist(out, 'dir'), mkdir(out); end
opt = weboptions('Timeout', 120);

fprintf('\nDownloading into %s\n', out);
nOK = 0;
for r = recs
    r = r{1};
    ok = true;
    for ext = {'.hea', '.dat', '.atr'}
        f = fullfile(out, [r ext{1}]);
        if exist(f, 'file') == 2
            fi = dir(f);  if fi.bytes > 0, continue; end
        end
        done = false;
        for attempt = 1:3
            try
                websave(f, [base [r ext{1}]], opt);
                done = true;  break
            catch err
                pause(2);
            end
        end
        if ~done
            ok = false;
            fprintf('  %s%s FAILED: %s\n', r, ext{1}, err.message);
        end
    end
    if ok, nOK = nOK + 1; fprintf('  %s ok\n', r); end
end
fprintf('%d of %d records complete.\n', nOK, numel(recs));

% sanity check (record 118 has 2278 beats; the noisy copies keep them)
try
    [x, beats] = read_mitdb(out, '118e06');
    fprintf('Check 118e06: %d samples (expect 650000), %d beats (expect about 2278)\n', ...
            numel(x), numel(beats));
catch err
    fprintf('Check of 118e06 failed: %s\n', err.message);
end
end

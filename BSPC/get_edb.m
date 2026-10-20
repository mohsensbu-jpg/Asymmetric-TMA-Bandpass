function get_edb()
%GET_EDB  European ST-T Database (90 records, 2 h, 250 Hz, 2 leads) from the
%   PhysioNet AWS mirror (PhysioNet as fallback), server-to-server. Only
%   channel 1 (int16) and the beat annotations are kept (edb/<rec>.mat),
%   the .dat file is deleted after extraction. Safe to re-run.
bases = {'https://physionet-open.s3.amazonaws.com/edb/1.0.0/', 'https://physionet.org/files/edb/1.0.0/'};
root = userpath;  if isempty(root), root = pwd; end
out = fullfile(root, 'edb');  if ~exist(out, 'dir'), mkdir(out); end
tmp = tempname;  mkdir(tmp);  opt = weboptions('Timeout', 300);
recs = {};
for b = 1:2
    try, t = webread([bases{b} 'RECORDS'], opt); recs = strsplit(strtrim(char(t))); break; catch, end
end
if isempty(recs), error('could not read the RECORDS list'); end
nOK = 0;
for k = 1:numel(recs)
    rec = strtrim(recs{k});  fm = fullfile(out, [rec '.mat']);
    if exist(fm, 'file') == 2, nOK = nOK + 1; continue; end
    ok = false;  t0 = tic;
    for attempt = 1:4
        base = bases{1 + (attempt > 2)};
        try
            for ext = {'.hea', '.atr', '.dat'}
                websave(fullfile(tmp, [rec ext{1}]), [base rec ext{1}], opt);
            end
            ok = true;  break
        catch err
            fprintf('  %s attempt %d failed: %s\n', rec, attempt, err.message);  pause(3);
        end
    end
    if ~ok, continue; end
    try
        [x, beats] = read_mitdb(tmp, rec);   % channel 1, format 212, beat labels
        x = int16(x);  fs = 250;  save(fm, 'x', 'beats', 'fs');
        nOK = nOK + 1;  fprintf('  %s ok (%d samples, %d beats, %.0f s)\n', rec, numel(x), numel(beats), toc(t0));
    catch err
        fprintf('  %s read FAILED: %s\n', rec, err.message);
    end
    delete(fullfile(tmp, [rec '.*']));
end
rmdir(tmp, 's');
fprintf('%d of %d records ready in %s\n', nOK, numel(recs), out);
end

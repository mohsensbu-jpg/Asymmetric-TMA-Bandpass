function get_incart()
%GET_INCART  St Petersburg INCART 12-lead Arrhythmia Database (75 records,
%   30 min, 257 Hz) from PhysioNet into MATLAB Drive. Server-to-server; only
%   lead II and the beat annotations are kept (about 1 MB per record), the
%   12-lead .dat file is deleted after extraction. Safe to re-run.
%   Output: <MATLAB Drive>/incart/Ixx.mat with x (int16, lead II), beats
%   (0-based sample indices of beat labels), fs, lead.
% the AWS open-data mirror of PhysioNet is usually much faster; PhysioNet is the fallback
bases = {'https://physionet-open.s3.amazonaws.com/incartdb/1.0.0/', ...
         'https://physionet.org/files/incartdb/1.0.0/'};
root = userpath;  if isempty(root), root = pwd; end
out  = fullfile(root, 'incart');  if ~exist(out, 'dir'), mkdir(out); end
tmp  = tempname;  mkdir(tmp);
opt  = weboptions('Timeout', 300);
nOK = 0;
for r = 1:75
    rec = sprintf('I%02d', r);  fm = fullfile(out, [rec '.mat']);
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
            fprintf('  %s attempt %d (%s) failed: %s\n', rec, attempt, base(9:22), err.message);  pause(3);
        end
    end
    if ~ok, continue; end
    try
        [x, beats, fs, lead] = read_wfdb_lead(tmp, rec, 'II');
        x = int16(x);  save(fm, 'x', 'beats', 'fs', 'lead');
        nOK = nOK + 1;  fprintf('  %s ok (%d samples, %d beats, lead %s, %.0f s)\n', rec, numel(x), numel(beats), lead, toc(t0));
    catch err
        fprintf('  %s read FAILED: %s\n', rec, err.message);
    end
    delete(fullfile(tmp, [rec '.*']));
end
rmdir(tmp, 's');
fprintf('%d of 75 records ready in %s\n', nOK, out);
end

function [x, beats, fs, lead] = read_wfdb_lead(folder, rec, want)
% format-16 (or 212 for 2 signals) WFDB reader, one lead chosen by name
L = strsplit(strtrim(fileread(fullfile(folder, [rec '.hea']))), newline);
h = strsplit(strtrim(L{1}));  nsig = str2double(h{2});  fs = str2double(h{3});
nsamp = str2double(h{4});
names = cell(nsig, 1);  fmt = zeros(nsig, 1);
for k = 1:nsig
    t = strsplit(strtrim(L{1+k}));  fmt(k) = str2double(regexprep(t{2}, '[^0-9].*', ''));
    names{k} = upper(t{end});
end
c = find(strcmp(names, upper(want)), 1);
if isempty(c), c = min(2, nsig); end
lead = names{c};
fid = fopen(fullfile(folder, [rec '.dat']), 'r');
if all(fmt == 16)
    v = fread(fid, inf, 'int16=>double');  fclose(fid);
    n = floor(numel(v)/nsig);  v = reshape(v(1:n*nsig), nsig, n);  x = v(c, :)';
elseif all(fmt == 212) && nsig == 2
    b = fread(fid, inf, 'uint8=>double');  fclose(fid);
    n = floor(numel(b)/3);  b = reshape(b(1:3*n), 3, n);
    if c == 1, s = b(1,:) + bitand(b(2,:), 15)*256;
    else,      s = b(3,:) + bitshift(b(2,:), -4)*256; end
    s(s > 2047) = s(s > 2047) - 4096;  x = s(:);
else
    fclose(fid);  error('unsupported format %s', mat2str(fmt'));
end
if ~isnan(nsamp) && numel(x) > nsamp, x = x(1:nsamp); end
beats = read_atr(fullfile(folder, [rec '.atr']), numel(x));
end

function beats = read_atr(file, N)
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
beats = pos(ismember(typ, [1:13 25 30 34 35 38 41]));
beats = beats(beats >= 0 & beats < N);
end

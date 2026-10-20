function [x, beats, types] = read_mitdb(folder, rec)
%READ_MITDB  Read one MIT-BIH record without any toolbox.
%   [x, beats, types] = read_mitdb(folder, rec)   rec: 100 or '118e06'
%   x     : first signal (channel 1, usually MLII), raw ADC units, double
%   beats : 0-based sample indices of BEAT annotations (as in the CSV files)
%   types : annotation codes of those beats
%   Reads format 212 (.dat) and the MIT annotation format (.atr).

if isnumeric(rec), rec = sprintf('%d', rec); end

%% ---- header: number of samples
hea = fileread(fullfile(folder, [rec '.hea']));
first = strsplit(strtrim(strtok(hea, sprintf('\n'))));
nsig = str2double(first{2});
nsamp = str2double(first{4});
if nsig ~= 2, error('expected 2 signals, header says %d', nsig); end

%% ---- signal, format 212 (two 12-bit samples in 3 bytes)
fid = fopen(fullfile(folder, [rec '.dat']), 'r');
b = fread(fid, inf, 'uint8=>double');  fclose(fid);
n = floor(numel(b) / 3);
b = reshape(b(1:3*n), 3, n);
s1 = b(1,:) + bitand(b(2,:), 15) * 256;
s1(s1 > 2047) = s1(s1 > 2047) - 4096;
x = s1(:);
if ~isnan(nsamp) && numel(x) > nsamp, x = x(1:nsamp); end

%% ---- annotations (.atr)
fid = fopen(fullfile(folder, [rec '.atr']), 'r');
a = fread(fid, inf, 'uint8=>double');  fclose(fid);
w = a(1:2:end-1) + 256 * a(2:2:end);       % 16-bit little-endian words
t = 0;  k = 1;  pos = zeros(numel(w), 1);  typ = zeros(numel(w), 1);  c = 0;
while k <= numel(w)
    A = floor(w(k) / 1024);  I = mod(w(k), 1024);
    if A == 0 && I == 0, break; end
    if A == 59                                % SKIP: 32-bit interval follows
        hi = w(k+1);  lo = w(k+2);
        v = hi * 65536 + lo;  if v >= 2^31, v = v - 2^32; end
        t = t + v;  k = k + 3;  continue
    elseif A == 63                            % AUX: I bytes of text follow
        k = k + 1 + ceil(I / 2);  continue
    elseif A >= 60                            % NUM / SUB / CHN: ignore
        k = k + 1;  continue
    end
    t = t + I;  c = c + 1;  pos(c) = t;  typ(c) = A;  k = k + 1;
end
pos = pos(1:c);  typ = typ(1:c);
isBeat = ismember(typ, [1:13 25 30 34 35 38 41]);   % standard QRS labels
beats = pos(isBeat);  types = typ(isBeat);
beats = beats(beats >= 0 & beats < numel(x));
end

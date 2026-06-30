function freqFeat = WM_Compute_FrequencyFeatures(ep, t, useChanIdx, cfg)
%
if nargin < 3, useChanIdx = []; end
if nargin < 4, cfg = struct(); end

% ---- struct -> numeric ----
ep = pickEpochIfStruct(ep, cfg);
t  = pickTimeIfStruct(t, cfg);

% ---- validate ----
if ~isnumeric(ep) || ndims(ep)~=3
    error('WM_Compute_FrequencyFeatures:badEp', 'ep must be numeric [chan x time x trial].');
end
if ~isnumeric(t) || ~isvector(t)
    error('WM_Compute_FrequencyFeatures:badT', 't must be numeric time vector.');
end

t = double(t(:)'); % row
[nChan, nTime, nTrial] = size(ep);
if numel(t) ~= nTime
    error('WM_Compute_FrequencyFeatures:timeMismatch', 'Length(t) != size(ep,2).');
end

if isempty(useChanIdx)
    useChanIdx = 1:nChan;
else
    useChanIdx = unique(useChanIdx(:)');
    useChanIdx = useChanIdx(useChanIdx>=1 & useChanIdx<=nChan);
end

% ---- fs from time vector ----
dt = median(diff(t(~isnan(t))));
if isempty(dt) || isnan(dt) || dt<=0
    if isfield(cfg,'fs') && ~isempty(cfg.fs)
        fs = double(cfg.fs);
    else
        fs = 1; % fallback
    end
else
    fs = 1/dt;
end

% ---- bands (Hz) ----
bands = struct();
bands.delta = [1 4];
bands.theta = [4 8];
bands.alpha = [8 13];
bands.beta  = [13 30];
bands.gamma = [30 45];

if isfield(cfg,'freq') && isstruct(cfg.freq) && isfield(cfg.freq,'bands') && isstruct(cfg.freq.bands)
    bands = cfg.freq.bands; % allow override
end

% ---- analysis mask (post-stim by default) ----
idxUse = (t >= 0);
if ~any(idxUse), idxUse = true(size(t)); end

ep = double(ep);

freqFeat = table();

bandNames = fieldnames(bands);

for tr = 1:nTrial
    row = table();
    for ch = useChanIdx
        x = squeeze(ep(ch, idxUse, tr));
        x = x(:)'; 
        x = x(~isnan(x));

        if numel(x) < 8
            % too short for PSD
            for b = 1:numel(bandNames)
                bn = bandNames{b};
                row.(sprintf('BP_%s_ch%02d', bn, ch))   = NaN;
                row.(sprintf('RBP_%s_ch%02d', bn, ch))  = NaN;
            end
            row.(sprintf('BP_total_ch%02d', ch)) = NaN;
            continue;
        end

        % PSD via pwelch (fallback to periodogram if pwelch missing)
        [Pxx, F] = safePSD(x, fs);

        if isempty(Pxx)
            for b = 1:numel(bandNames)
                bn = bandNames{b};
                row.(sprintf('BP_%s_ch%02d', bn, ch))   = NaN;
                row.(sprintf('RBP_%s_ch%02d', bn, ch))  = NaN;
            end
            row.(sprintf('BP_total_ch%02d', ch)) = NaN;
            continue;
        end

        % total power in [1..45] Hz
        totalBand = [1 45];
        BP_total = trapzBand(Pxx, F, totalBand(1), totalBand(2));
        row.(sprintf('BP_total_ch%02d', ch)) = BP_total;

        for b = 1:numel(bandNames)
            bn = bandNames{b};
            fr = bands.(bn);
            bp = trapzBand(Pxx, F, fr(1), fr(2));
            row.(sprintf('BP_%s_ch%02d', bn, ch)) = bp;

            if ~isnan(BP_total) && BP_total > 0
                row.(sprintf('RBP_%s_ch%02d', bn, ch)) = bp / BP_total;
            else
                row.(sprintf('RBP_%s_ch%02d', bn, ch)) = NaN;
            end
        end
    end
    freqFeat = [freqFeat; row]; %#ok<AGROW>
end

end

% ================= helpers =================
function epOut = pickEpochIfStruct(epIn, cfg)
epOut = epIn;
if ~isstruct(epIn), return; end
ph = "";
if isfield(cfg,'phase') && ~isempty(cfg.phase), ph = string(cfg.phase); end
if ph ~= "" && isfield(epIn, ph), epOut = epIn.(ph); return; end
f = fieldnames(epIn);
for i=1:numel(f)
    v = epIn.(f{i});
    if isnumeric(v) && ndims(v)==3 && ~isempty(v), epOut = v; return; end
end
error('WM_Compute_FrequencyFeatures:structEpEmpty','No numeric 3D field in ep struct.');
end

function tOut = pickTimeIfStruct(tIn, cfg)
tOut = tIn;
if ~isstruct(tIn), return; end
ph = "";
if isfield(cfg,'phase') && ~isempty(cfg.phase), ph = string(cfg.phase); end
if ph ~= "" && isfield(tIn, ph), tOut = tIn.(ph); return; end
f = fieldnames(tIn);
for i=1:numel(f)
    v = tIn.(f{i});
    if isnumeric(v) && isvector(v) && ~isempty(v), tOut = v; return; end
end
error('WM_Compute_FrequencyFeatures:structTVecEmpty','No numeric vector field in t struct.');
end

function [Pxx,F] = safePSD(x, fs)
Pxx = []; F = [];
try
    if exist('pwelch','file') == 2
        n = numel(x);
        wlen = max(64, 2^nextpow2(min(n, 512)));
        nover = round(0.5*wlen);
        nfft  = max(wlen, 256);
        [Pxx, F] = pwelch(x, wlen, nover, nfft, fs);
    else
        % fallback: periodogram
        nfft = max(256, 2^nextpow2(numel(x)));
        [Pxx, F] = periodogram(x, [], nfft, fs);
    end
catch
    Pxx = []; F = [];
end
end

function bp = trapzBand(Pxx, F, f1, f2)
mask = (F >= f1) & (F <= f2);
if ~any(mask)
    bp = NaN;
else
    bp = trapz(F(mask), Pxx(mask));
end
end

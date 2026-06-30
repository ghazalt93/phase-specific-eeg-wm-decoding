function tempFeat = WM_Compute_TemporalFeatures(ep, t, useChanIdx, cfgOrTimeWin)

if nargin < 3, useChanIdx = []; end
if nargin < 4, cfgOrTimeWin = []; end

% ---- interpret 4th arg ----
cfg = struct();
timeWin = [];

if isstruct(cfgOrTimeWin)
    cfg = cfgOrTimeWin;
elseif isnumeric(cfgOrTimeWin) && numel(cfgOrTimeWin)==2
    timeWin = double(cfgOrTimeWin(:))';
elseif isempty(cfgOrTimeWin)
    % nothing
else
    error('WM_Compute_TemporalFeatures:badArg','4th argument must be [] or [t1 t2] or cfg struct.');
end

% ---- struct -> numeric ----
ep = pickEpochIfStruct(ep, cfg);
t  = pickTimeIfStruct(t, cfg);

% ---- validate ----
if ~isnumeric(ep) || ndims(ep)~=3
    error('WM_Compute_TemporalFeatures:badEp','ep must be numeric [chan x time x trial].');
end
if ~isnumeric(t) || ~isvector(t)
    error('WM_Compute_TemporalFeatures:badT','t must be numeric time vector.');
end

t = double(t(:)'); % row
[nChan, nTime, nTrial] = size(ep);
if numel(t) ~= nTime
    error('WM_Compute_TemporalFeatures:timeMismatch','Length(t) != size(ep,2).');
end

% ---- channels ----
if isempty(useChanIdx)
    useChanIdx = 1:nChan;
else
    useChanIdx = unique(useChanIdx(:)');
    useChanIdx = useChanIdx(useChanIdx>=1 & useChanIdx<=nChan);
end

% ---- choose time window/index ----
if isfield(cfg,'idxUse') && ~isempty(cfg.idxUse)
    idx = logical(cfg.idxUse);
    if numel(idx) ~= nTime
        error('WM_Compute_TemporalFeatures:badIdxUse','cfg.idxUse must match size(ep,2).');
    end
elseif isfield(cfg,'tRange') && ~isempty(cfg.tRange) && numel(cfg.tRange)==2
    tr = double(cfg.tRange(:))';
    idx = (t >= tr(1)) & (t <= tr(2));
elseif ~isempty(timeWin)
    idx = (t >= timeWin(1)) & (t <= timeWin(2));
else
    % default: post window
    idx = (t >= 0);
end
if ~any(idx), idx = true(1,nTime); end

tUse = t(idx);

% ---- output ----
nCh = numel(useChanIdx);
nF = 4;
out = nan(nTrial, nCh*nF);

varNames = cell(1, nCh*nF);
k = 0;
for ic=1:nCh
    ch = useChanIdx(ic);
    varNames{k+1} = sprintf('RMS_ch%02d',     ch);
    varNames{k+2} = sprintf('PeakAbs_ch%02d', ch);
    varNames{k+3} = sprintf('AUCabs_ch%02d',  ch);
    varNames{k+4} = sprintf('LineLen_ch%02d', ch);
    k = k + nF;
end
varNames = matlab.lang.makeUniqueStrings(varNames);

ep = double(ep);

for tr=1:nTrial
    k = 0;
    for ic=1:nCh
        ch = useChanIdx(ic);
        x = double(squeeze(ep(ch, idx, tr)));
        x = x(:);

        % remove non-finite
        good = isfinite(x);
        x = x(good);
        tt = tUse(:);
        tt = tt(good);

        if isempty(x)
            out(tr, k+(1:4)) = NaN;
        else
            rmsv = sqrt(mean(x.^2));
            peak = max(abs(x));

            if numel(x) >= 2 && numel(tt) == numel(x)
                auc = trapz(tt, abs(x));
            else
                auc = sum(abs(x));
            end

            if numel(x) >= 2
                ll = sum(abs(diff(x)));
            else
                ll = 0;
            end

            out(tr, k+1) = rmsv;
            out(tr, k+2) = peak;
            out(tr, k+3) = auc;
            out(tr, k+4) = ll;
        end

        k = k + nF;
    end
end

tempFeat = array2table(out, 'VariableNames', varNames);
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
error('WM_Compute_TemporalFeatures:structEpEmpty','No numeric 3D field in ep struct.');
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
error('WM_Compute_TemporalFeatures:structTVecEmpty','No numeric vector field in t struct.');
end

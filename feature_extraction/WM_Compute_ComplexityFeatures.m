function compFeat = WM_Compute_ComplexityFeatures(ep, t, useChanIdx, cfgComp)
%
if nargin < 3, useChanIdx = []; end
if nargin < 4 || isempty(cfgComp), cfgComp = struct(); end

% ---- struct -> numeric ----
ep = pickEpochIfStruct(ep, cfgComp);
t  = pickTimeIfStruct(t,  cfgComp);

% ---- validate ----
if ~isnumeric(ep) || ndims(ep) ~= 3
    error('WM_Compute_ComplexityFeatures:badEp','ep must be numeric [chan x time x trial].');
end
if ~isnumeric(t) || ~isvector(t)
    error('WM_Compute_ComplexityFeatures:badT','t must be numeric time vector.');
end

t = double(t(:)'); % row
[nChan, nTime, nTrial] = size(ep);
if numel(t) ~= nTime
    error('WM_Compute_ComplexityFeatures:timeMismatch','Length(t) != size(ep,2).');
end

% ---- defaults ----
if ~isfield(cfgComp,'ds') || isempty(cfgComp.ds), cfgComp.ds = 1; end
if ~isfield(cfgComp,'maxLen') || isempty(cfgComp.maxLen), cfgComp.maxLen = inf; end
if ~isfield(cfgComp,'binarize') || isempty(cfgComp.binarize), cfgComp.binarize = 'median'; end
if ~isfield(cfgComp,'normalize') || isempty(cfgComp.normalize), cfgComp.normalize = true; end

ds = max(1, round(double(cfgComp.ds)));
maxLen = double(cfgComp.maxLen);
if ~isfinite(maxLen) || maxLen < 2, maxLen = inf; end

% ---- channels ----
if isempty(useChanIdx)
    useChanIdx = 1:nChan;
else
    useChanIdx = unique(useChanIdx(:))';
    useChanIdx = useChanIdx(useChanIdx>=1 & useChanIdx<=nChan);
end

% ---- choose time window ----
if isfield(cfgComp,'idxUse') && ~isempty(cfgComp.idxUse)
    idx = logical(cfgComp.idxUse);
    if numel(idx) ~= nTime
        error('WM_Compute_ComplexityFeatures:badIdxUse','cfgComp.idxUse must match size(ep,2).');
    end
elseif isfield(cfgComp,'tRange') && ~isempty(cfgComp.tRange) && numel(cfgComp.tRange)==2
    tr = double(cfgComp.tRange(:))';
    idx = (t >= tr(1)) & (t <= tr(2));
elseif isfield(cfgComp,'time_win') && ~isempty(cfgComp.time_win) && numel(cfgComp.time_win)==2
    tw = double(cfgComp.time_win(:))';
    idx = (t >= tw(1)) & (t <= tw(2));
else
    idx = (t >= 0);
end
if ~any(idx), idx = true(1,nTime); end

% ---- output ----
nCh = numel(useChanIdx);
out = nan(nTrial, nCh);

varNames = cell(1,nCh);
for ic=1:nCh
    varNames{ic} = sprintf('LZC_ch%02d', useChanIdx(ic));
end
varNames = matlab.lang.makeUniqueStrings(varNames);

ep = double(ep);

for tr=1:nTrial
    for ic=1:nCh
        ch = useChanIdx(ic);

        x = double(squeeze(ep(ch, idx, tr)));
        x(~isfinite(x)) = 0;

        x = x(1:ds:end);
        if numel(x) > maxLen
            x = x(1:maxLen);
        end

        if numel(x) < 2
            out(tr, ic) = NaN;
            continue
        end

        switch lower(string(cfgComp.binarize))
            case 'median'
                thr = median(x);
            case 'mean'
                thr = mean(x);
            case 'zero'
                thr = 0;
            otherwise
                thr = median(x);
        end
        b = (x > thr);

        out(tr, ic) = lzc_kaspar(b, cfgComp.normalize);
    end
end

compFeat = array2table(out, 'VariableNames', varNames);
end

% ================= helpers =================
function epOut = pickEpochIfStruct(epIn, cfg)
epOut = epIn;
if ~isstruct(epIn), return; end
ph = "";
if isfield(cfg,'phase') && ~isempty(cfg.phase), ph = string(cfg.phase); end
if ph ~= "" && isfield(epIn, ph)
    epOut = epIn.(ph); return;
end
f = fieldnames(epIn);
for i=1:numel(f)
    v = epIn.(f{i});
    if isnumeric(v) && ndims(v)==3 && ~isempty(v)
        epOut = v; return;
    end
end
error('WM_Compute_ComplexityFeatures:structEpEmpty','No numeric 3D field in ep struct.');
end

function tOut = pickTimeIfStruct(tIn, cfg)
tOut = tIn;
if ~isstruct(tIn), return; end
ph = "";
if isfield(cfg,'phase') && ~isempty(cfg.phase), ph = string(cfg.phase); end
if ph ~= "" && isfield(tIn, ph)
    tOut = tIn.(ph); return;
end
f = fieldnames(tIn);
for i=1:numel(f)
    v = tIn.(f{i});
    if isnumeric(v) && isvector(v) && ~isempty(v)
        tOut = v; return;
    end
end
error('WM_Compute_ComplexityFeatures:structTVecEmpty','No numeric vector field in t struct.');
end

function c = lzc_kaspar(b, doNormalize)
b = logical(b(:))';
n = numel(b);

if n < 2
    c = NaN; return;
end

i = 1; k = 1; l = 1;
cCnt = 1;
k_max = 1;

while true
    if i+k > n || l+k > n
        cCnt = cCnt + 1;
        break;
    end

    if b(i+k) == b(l+k)
        k = k + 1;
        if l+k > n
            cCnt = cCnt + 1;
            break;
        end
    else
        if k > k_max, k_max = k; end
        i = i + 1;
        if i == l
            cCnt = cCnt + 1;
            l = l + k_max;
            if l > n
                break;
            end
            i = 1;
            k = 1;
            k_max = 1;
        else
            k = 1;
        end
    end
end

if doNormalize
    c = cCnt * (log2(n) / n);
else
    c = cCnt;
end
end

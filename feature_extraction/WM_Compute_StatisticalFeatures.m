function statFeat = WM_Compute_StatisticalFeatures(ep, t, useChanIdx, cfg)


if nargin < 3, useChanIdx = []; end
if nargin < 4, cfg = struct(); end

% ---- struct -> numeric ----
ep = pickEpochIfStruct(ep, cfg);
t  = pickTimeIfStruct(t, cfg);

% ---- validate ----
if ~isnumeric(ep) || ndims(ep)~=3
    error('WM_Compute_StatisticalFeatures:badEp', 'ep must be numeric [chan x time x trial].');
end
if ~isnumeric(t) || ~isvector(t)
    error('WM_Compute_StatisticalFeatures:badT', 't must be numeric time vector.');
end

t = double(t(:)'); % row
[nChan, nTime, nTrial] = size(ep);
if numel(t) ~= nTime
    error('WM_Compute_StatisticalFeatures:timeMismatch', 'Length(t) != size(ep,2).');
end

% ---- channels ----
if isempty(useChanIdx)
    useChanIdx = 1:nChan;
else
    useChanIdx = unique(useChanIdx(:)');
    useChanIdx = useChanIdx(useChanIdx>=1 & useChanIdx<=nChan);
end

% ---- choose time window ----
if isfield(cfg,'idxUse') && ~isempty(cfg.idxUse)
    idxUse = logical(cfg.idxUse);
    if numel(idxUse) ~= nTime
        error('WM_Compute_StatisticalFeatures:badIdxUse', 'cfg.idxUse must match size(ep,2).');
    end
elseif isfield(cfg,'tRange') && ~isempty(cfg.tRange) && numel(cfg.tRange)==2
    t1 = cfg.tRange(1); t2 = cfg.tRange(2);
    idxUse = (t >= t1) & (t <= t2);
else
    % default: post (t>=0)
    idxUse = (t >= 0);
end
if ~any(idxUse), idxUse = true(1,nTime); end

nCh = numel(useChanIdx);
out = nan(nTrial, nCh*4);

varNames = cell(1, nCh*4);
k = 0;
for ic=1:nCh
    ch = useChanIdx(ic);
    varNames{k+1} = sprintf('Mean_ch%02d', ch);
    varNames{k+2} = sprintf('Std_ch%02d',  ch);
    varNames{k+3} = sprintf('Skew_ch%02d', ch);
    varNames{k+4} = sprintf('Kurt_ch%02d', ch);
    k = k+4;
end
varNames = matlab.lang.makeUniqueStrings(varNames);

ep = double(ep);

for tr=1:nTrial
    k = 0;
    for ic=1:nCh
        ch = useChanIdx(ic);
        sig = double(squeeze(ep(ch, idxUse, tr)));
        sig = sig(isfinite(sig));

        if numel(sig) < 3
            out(tr, k+(1:4)) = NaN;
        else
            mu = mean(sig);
            sd = std(sig,0);

            if sd == 0 || ~isfinite(sd)
                sk = NaN; ku = NaN;
            else
                z = (sig - mu) ./ sd;
                sk = mean(z.^3);
                ku = mean(z.^4); % raw kurtosis (not excess)
            end

            out(tr, k+1) = mu;
            out(tr, k+2) = sd;
            out(tr, k+3) = sk;
            out(tr, k+4) = ku;
        end

        k = k + 4;
    end
end

statFeat = array2table(out, 'VariableNames', varNames);
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
error('WM_Compute_StatisticalFeatures:structEpEmpty','No numeric 3D field in ep struct.');
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
error('WM_Compute_StatisticalFeatures:structTVecEmpty','No numeric vector field in t struct.');
end

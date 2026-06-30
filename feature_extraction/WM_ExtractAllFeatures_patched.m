function [featAll, trialInfo, meta] = WM_ExtractAllFeatures_patched(EEG, epochs, timeVec, trialInfo, cfg)


if nargin < 5 || isempty(cfg), cfg = struct(); end
if ~isfield(cfg,'phase'), error('cfg.phase is required'); end
phase = char(string(cfg.phase));

% ---------- Validate epochs/timeVec ----------
if ~isstruct(epochs) || ~isfield(epochs, phase)
    error('WM_ExtractAllFeatures:MissingEpochsField', 'epochs.%s not found', phase);
end
if ~isstruct(timeVec) || ~isfield(timeVec, phase)
    error('WM_ExtractAllFeatures:MissingTimeVecField', 'timeVec.%s not found', phase);
end

ep = epochs.(phase);      % [chan x time x trial]
t  = timeVec.(phase);
t  = double(t(:)');       % row

if isempty(ep) || ndims(ep) ~= 3 || size(ep,3) == 0
    featAll = table();
    warning('[SKIP] empty epochs for phase=%s', phase);
    meta = struct('phase',phase);
    return
end

nTrial = size(ep,3);

% ---------- Align trialInfo to epochs ----------
if ~isempty(trialInfo) && istable(trialInfo) && height(trialInfo) ~= nTrial
    warning('[ALIGN] phase=%s: trialInfo=%d but epochs=%d -> trunc to min', ...
        phase, height(trialInfo), nTrial);
    m = min(height(trialInfo), nTrial);
    ep = ep(:,:,1:m);
    trialInfo = trialInfo(1:m,:);
    nTrial = m;
end

% ---------- Ensure trialInfo ----------
if isempty(trialInfo) || ~istable(trialInfo)
    trialInfo = table();
end
if height(trialInfo) ~= nTrial
    trialInfo = resize_trialinfo(trialInfo, nTrial);
end

% ---------- Ensure epochType in trialInfo ----------
if ~ismember('epochType', trialInfo.Properties.VariableNames)
    trialInfo.epochType = repmat(string(phase), nTrial, 1);
end

% ---------- Channels ----------
useChanIdx = 1:size(ep,1);
if isfield(cfg,'useChanIdx') && ~isempty(cfg.useChanIdx)
    useChanIdx = cfg.useChanIdx(:)';
end
useChanIdx = useChanIdx(useChanIdx>=1 & useChanIdx<=size(ep,1));

% ---------- Sampling rate ----------
fs = NaN;
if ~isempty(EEG) && isstruct(EEG) && isfield(EEG,'srate') && ~isempty(EEG.srate)
    fs = double(EEG.srate);
end
if ~(isfinite(fs) && fs>0)
    dt = median(diff(t(~isnan(t))));
    if isfinite(dt) && dt>0
        fs = 1/dt;
    end
end

% ---------- defaults for baseline / windows ----------
if ~isfield(cfg,'baseline_sec') || isempty(cfg.baseline_sec)
    cfg.baseline_sec = 0.5;
end
if ~isfield(cfg,'epsPow') || isempty(cfg.epsPow)
    cfg.epsPow = 1e-12;
end

% Optional: map cfg.freqBands -> cfg.freq.bands for WM_Compute_FrequencyFeatures
if isfield(cfg,'freqBands') && isstruct(cfg.freqBands) && ~isempty(fieldnames(cfg.freqBands))
    if ~isfield(cfg,'freq') || ~isstruct(cfg.freq), cfg.freq = struct(); end
    cfg.freq.bands = cfg.freqBands;
end

% ---------- Post window: use cfg.win.<phase> if present else default (t>=0) ----------
if isfield(cfg,'win') && isstruct(cfg.win) && isfield(cfg.win, phase) && numel(cfg.win.(phase))==2
    tw = double(cfg.win.(phase));
    idxPost = (t >= tw(1)) & (t <= tw(2));
else
    idxPost = (t >= 0);
end
if ~any(idxPost), idxPost = true(size(t)); end

% ---------- Baseline window: prefer pre-zero, else first baseline_sec ----------
idxBase = (t >= -cfg.baseline_sec) & (t < 0);
if ~any(idxBase)
    nBase = max(8, round(cfg.baseline_sec * fs));
    idxBase = false(size(t));
    idxBase(1:min(nBase,numel(t))) = true;
end

epPost = ep(:, idxPost, :);  tPost = t(idxPost);
epBase = ep(:, idxBase, :);  tBase = t(idxBase);

% ---------- Prepare outputs ----------
freqFeat  = emptyTableNRows(nTrial);
statFeat  = emptyTableNRows(nTrial);
tempFeat  = emptyTableNRows(nTrial);
compFeat  = emptyTableNRows(nTrial);
behavFeat = emptyTableNRows(nTrial);

% ===================== Frequency (logbpblwin) =====================
if isfield(cfg,'do_freq') && cfg.do_freq
    if ~(isfinite(fs) && fs>0)
        warning('[FREQ] fs invalid -> skipping frequency features.');
    else
        try
            cfgF = cfg; cfgF.fs = fs; cfgF.phase = phase;
            Fp = WM_Compute_FrequencyFeatures(epPost, tPost, useChanIdx, cfgF);
            Fb = WM_Compute_FrequencyFeatures(epBase, tBase, useChanIdx, cfgF);

            vn = string(Fp.Properties.VariableNames);
            isBP = startsWith(vn,'BP_') & ~contains(vn,'total');
            bpCols = vn(isBP);

            freqFeat = table();
            for c = 1:numel(bpCols)
                col = bpCols(c);
                xPost = double(Fp.(col));
                xBase = double(Fb.(col));
                feat = log10(xPost + cfg.epsPow) - log10(xBase + cfg.epsPow);
                newName = matlab.lang.makeValidName(regexprep(col,'^BP_','logbpblwin_'));
                freqFeat.(newName) = feat;
            end
        catch ME
            warning('[FREQ] failed: %s', ME.message);
        end
    end
end

% ===================== Statistical =====================
if isfield(cfg,'do_stats') && cfg.do_stats && exist('WM_Compute_StatisticalFeatures','file')==2
    try
        statPost = call_stats(epPost, tPost, useChanIdx, phase);
        if ~isfield(cfg,'stat_blwin') || cfg.stat_blwin
            statBase = call_stats(epBase, tBase, useChanIdx, phase);
            statFeat = statPost; statFeat{:,:} = statPost{:,:} - statBase{:,:};
            prefix = "statblwin_";
        else
            statFeat = statPost;
            prefix = "stat_";
        end
        statFeat.Properties.VariableNames = matlab.lang.makeUniqueStrings(cellstr(prefix + string(statFeat.Properties.VariableNames)));
    catch ME
        warning('[STAT] failed: %s', ME.message);
    end
end

% ===================== Temporal =====================
if isfield(cfg,'do_temporal') && cfg.do_temporal && exist('WM_Compute_TemporalFeatures','file')==2
    try
        % most legacy versions accept [t1 t2]
        timeWin = [min(tPost) max(tPost)];
        tempPost = WM_Compute_TemporalFeatures(epPost, tPost, useChanIdx, timeWin);

        if ~isfield(cfg,'temp_blwin') || cfg.temp_blwin
            timeWinB = [min(tBase) max(tBase)];
            tempBase = WM_Compute_TemporalFeatures(epBase, tBase, useChanIdx, timeWinB);
            tempFeat = tempPost; tempFeat{:,:} = tempPost{:,:} - tempBase{:,:};
            prefix = "tempblwin_";
        else
            tempFeat = tempPost;
            prefix = "temp_";
        end

        tempFeat.Properties.VariableNames = matlab.lang.makeUniqueStrings(cellstr(prefix + string(tempFeat.Properties.VariableNames)));
    catch ME
        warning('[TEMP] failed: %s', ME.message);
    end
end

% ===================== Complexity =====================
if isfield(cfg,'do_complexity') && cfg.do_complexity && exist('WM_Compute_ComplexityFeatures','file')==2
    try
        cfgComp = struct();
        if isfield(cfg,'comp') && ~isempty(cfg.comp), cfgComp = cfg.comp; end
        cfgComp.phase = phase;

        % safe defaults
        if ~isfield(cfgComp,'time_win') || isempty(cfgComp.time_win)
            cfgComp.time_win = [min(tPost) max(tPost)];
        end
        if ~isfield(cfgComp,'ds') || isempty(cfgComp.ds), cfgComp.ds = 2; end
        if ~isfield(cfgComp,'maxLen') || isempty(cfgComp.maxLen), cfgComp.maxLen = 1200; end

        compPost = WM_Compute_ComplexityFeatures(epPost, tPost, useChanIdx, cfgComp);

        if ~isfield(cfg,'comp_blwin') || cfg.comp_blwin
            cfgCompB = cfgComp; cfgCompB.time_win = [min(tBase) max(tBase)];
            compBase = WM_Compute_ComplexityFeatures(epBase, tBase, useChanIdx, cfgCompB);
            compFeat = compPost; compFeat{:,:} = compPost{:,:} - compBase{:,:};
            prefix = "compblwin_";
        else
            compFeat = compPost;
            prefix = "comp_";
        end

        compFeat.Properties.VariableNames = matlab.lang.makeUniqueStrings(cellstr(prefix + string(compFeat.Properties.VariableNames)));
    catch ME
        warning('[COMP] failed: %s', ME.message);
    end
end

% ===================== Behavioral =====================
if exist('WM_Compute_BehavioralFeatures','file')==2
    try
        % prefer (trialInfo,cfg) if supported
        behavFeat = call_behav(trialInfo, phase);
        behavFeat.Properties.VariableNames = matlab.lang.makeUniqueStrings(cellstr("beh_" + string(behavFeat.Properties.VariableNames)));
    catch ME
        warning('[BEH] failed: %s', ME.message);
    end
end

% ---------- enforce same height ----------
freqFeat  = enforceNRows(freqFeat,  nTrial);
statFeat  = enforceNRows(statFeat,  nTrial);
tempFeat  = enforceNRows(tempFeat,  nTrial);
compFeat  = enforceNRows(compFeat,  nTrial);
behavFeat = enforceNRows(behavFeat, nTrial);

% ---------- unique names ----------
freqFeat  = makeVarNamesUnique(freqFeat);
statFeat  = makeVarNamesUnique(statFeat);
tempFeat  = makeVarNamesUnique(tempFeat);
compFeat  = makeVarNamesUnique(compFeat);
behavFeat = makeVarNamesUnique(behavFeat);

% ---------- combine ----------
featAll = [freqFeat statFeat tempFeat compFeat behavFeat];

% ---------- meta ----------
meta = struct();
meta.phase = phase;
meta.fs = fs;
meta.baseline_sec = cfg.baseline_sec;
meta.useChanIdx = useChanIdx;
if isfield(cfg,'subjectName'), meta.subjectName = cfg.subjectName; end
if isfield(cfg,'runTag'), meta.runTag = cfg.runTag; end
if isfield(cfg,'freqBands'), meta.freqBands = cfg.freqBands; end

fprintf('[OK] phase=%s | features=%d rows, %d cols\n', phase, height(featAll), width(featAll));

end

% =================== local helpers ===================
function T = emptyTableNRows(n)
T = table('Size',[n 0],'VariableTypes',{},'VariableNames',{});
end

function T = makeVarNamesUnique(T)
if ~istable(T) || width(T)==0, return; end
T.Properties.VariableNames = matlab.lang.makeUniqueStrings(T.Properties.VariableNames);
end

function T = enforceNRows(T, n)
if ~istable(T) || width(T)==0
    T = emptyTableNRows(n); return
end
h = height(T);
if h==n, return; end
if h > n
    T = T(1:n,:);
else
    pad = n-h;
    vars = T.Properties.VariableNames;
    Tpad = table();
    for j=1:numel(vars)
        ref = T.(vars{j});
        Tpad.(vars{j}) = fillerLike(ref, pad);
    end
    T = [T; Tpad];
end
end

function f = fillerLike(refCol, nRows)
if isstring(refCol)
    f = repmat(missing, nRows, 1);
elseif iscell(refCol)
    f = repmat({''}, nRows, 1);
elseif iscategorical(refCol)
    f = categorical(repmat({''}, nRows, 1));
elseif isnumeric(refCol) || islogical(refCol)
    f = nan(nRows, 1);
else
    f = repmat(missing, nRows, 1);
end
end

function trialInfo = resize_trialinfo(trialInfo, nTrial)
if ~istable(trialInfo)
    trialInfo = table(); %#ok<NASGU>
    trialInfo = table();
end
if height(trialInfo)==nTrial, return; end
if height(trialInfo) > nTrial
    trialInfo = trialInfo(1:nTrial,:);
else
    padN = nTrial - height(trialInfo);
    pad = trialInfo([]);
    for v = 1:width(trialInfo)
        pad.(trialInfo.Properties.VariableNames{v}) = fillerLike(trialInfo.(v), padN);
    end
    trialInfo = [trialInfo; pad];
end
end

function stat = call_stats(ep, t, useChanIdx, phase)
% try 4-arg robust then fallback 3-arg legacy
try
    stat = WM_Compute_StatisticalFeatures(ep, t, useChanIdx, struct('phase',phase,'idxUse',true(1,numel(t))));
catch
    stat = WM_Compute_StatisticalFeatures(ep, t, useChanIdx);
end
end

function beh = call_behav(trialInfo, phase)
try
    beh = WM_Compute_BehavioralFeatures(trialInfo, struct('phase',phase));
catch
    beh = WM_Compute_BehavioralFeatures(trialInfo);
end
end


function [epOut, tOut] = coerce_ep_t_local(epIn, tIn)
epOut = epIn;
tOut  = tIn;

try
    if iscell(epOut)
        while iscell(epOut) && numel(epOut)==1
            epOut = epOut{1};
        end
        if iscell(epOut) && ~isempty(epOut) && all(cellfun(@(x) isnumeric(x) && ndims(x)==2, epOut))
            epOut = cat(3, epOut{:});
        end
    end
catch
end

try
    if iscell(tOut)
        while iscell(tOut) && numel(tOut)==1
            tOut = tOut{1};
        end
        if iscell(tOut) && ~isempty(tOut) && all(cellfun(@(x) isnumeric(x) && isvector(x), tOut))
            tOut = tOut{1};
        end
        if iscellstr(tOut)
            tOut = str2double(string(tOut));
        end
    end
    if isstring(tOut), tOut = str2double(tOut); end
    if ischar(tOut),   tOut = str2double(string(tOut)); end
catch
end
end

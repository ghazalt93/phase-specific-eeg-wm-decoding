function [ep, tvec, ti] = get_phase_epochs_trialinfo(varargin)

% ---- parse inputs ----
phase = '';
T_phases = [];
fs = [];
epochSec = [];
XEEG = [];
epochsAll = struct();
timeVecAll = struct();
trialInfoAll = [];

args = varargin;

% identify phase (char/string)
for k = 1:numel(args)
    if (ischar(args{k}) && ~isempty(args{k})) || (isstring(args{k}) && isscalar(args{k}))
        cand = args{k};
        if isstring(cand); cand = char(cand); end
        % accept only known phases if possible
        if any(strcmpi(cand, {'stim','maint','retr'}))
            phase = lower(cand);
            args{k} = [];
            break;
        end
    end
end

% identify table T_phases
for k = 1:numel(args)
    if istable(args{k})
        T_phases = args{k};
        args{k} = [];
        break;
    end
end

% identify fs (scalar numeric)
for k = 1:numel(args)
    if isnumeric(args{k}) && isscalar(args{k}) && isfinite(args{k}) && args{k} > 0
        fs = double(args{k});
        args{k} = [];
        break;
    end
end

% identify XEEG (numeric 2D)
for k = 1:numel(args)
    if isnumeric(args{k}) && ismatrix(args{k}) && ~isscalar(args{k})
        % Heuristic: XEEG is usually [nSamples x nChan]
        if size(args{k},1) > 10 && size(args{k},2) >= 1
            XEEG = args{k};
            args{k} = [];
            break;
        end
    end
end

% identify epochSec (numeric scalar or vector)
for k = 1:numel(args)
    if isnumeric(args{k}) && ~isempty(args{k}) && isvector(args{k})
        epochSec = double(args{k}(:));
        args{k} = [];
        break;
    end
end

% identify epochsAll/timeVecAll (struct)
for k = 1:numel(args)
    if isstruct(args{k})
        if isempty(fieldnames(epochsAll))
            epochsAll = args{k};
            args{k} = [];
        elseif isempty(fieldnames(timeVecAll))
            timeVecAll = args{k};
            args{k} = [];
        end
    end
end

% identify trialInfoAll (table/struct)
for k = 1:numel(args)
    if istable(args{k}) || isstruct(args{k})
        trialInfoAll = args{k};
        args{k} = [];
        break;
    end
end

if isempty(phase)
    % if phase wasn't provided explicitly, try cfg.phase style inside remaining args
    phase = 'stim';
end

% ---- fast path: use provided epochsAll/timeVecAll ----
if isstruct(epochsAll) && isfield(epochsAll, phase) && ~isempty(epochsAll.(phase)) && isnumeric(epochsAll.(phase))
    ep = epochsAll.(phase);
    if isstruct(timeVecAll) && isfield(timeVecAll, phase) && ~isempty(timeVecAll.(phase))
        tvec = timeVecAll.(phase);
    else
        if isempty(fs)
            fs = 250;
        end
        tvec = (0:size(ep,2)-1) ./ fs;
    end

    % trial info
    ti = [];
    if istable(trialInfoAll)
        if ismember('epochType', trialInfoAll.Properties.VariableNames)
            m = strcmpi(trialInfoAll.epochType, phase);
            if any(m)
                ti = trialInfoAll(m,:);
            end
        end
        if isempty(ti)
            ti = trialInfoAll;
        end
    elseif isstruct(trialInfoAll) && isfield(trialInfoAll, phase)
        ti = trialInfoAll.(phase);
    end
    if ~istable(ti)
        ti = table((1:size(ep,3))','VariableNames',{'trial'});
    end
    if ~ismember('epochType', ti.Properties.VariableNames)
        ti.epochType = repmat({phase}, height(ti), 1);
    end
    return;
end

% ---- fallback: cut from T_phases + XEEG ----
if isempty(fs)
    fs = 250;
end
if isempty(epochSec)
    epochSec = 1; % default
end

if isempty(T_phases) || ~istable(T_phases)
    error('get_phase_epochs_trialinfo:MissingTphases', 'T_phases table is required when epochsAll is not provided.');
end
if isempty(XEEG) || ~isnumeric(XEEG)
    error('get_phase_epochs_trialinfo:MissingEEG', 'XEEG numeric matrix is required when epochsAll is not provided.');
end

% map columns
switch lower(phase)
    case 'stim'
        c0 = 'StimSample'; c1 = 'StimEnd';
    case 'maint'
        c0 = 'MaintSample'; c1 = 'MaintEnd';
    case 'retr'
        c0 = 'RetrSample'; c1 = 'RetrEnd';
    otherwise
        error('get_phase_epochs_trialinfo:BadPhase', 'Unknown phase: %s', phase);
end

if ~ismember(c0, T_phases.Properties.VariableNames)
    error('get_phase_epochs_trialinfo:MissingCol', 'T_phases missing column: %s', c0);
end

s0 = double(T_phases.(c0));
if any(~isfinite(s0))
    s0(~isfinite(s0)) = NaN;
end

hasEnd = ismember(c1, T_phases.Properties.VariableNames);
if hasEnd
    sEnd = double(T_phases.(c1));
else
    sEnd = NaN(size(s0));
end

% 0-based detection: if any zeros present, shift to 1-based
if any(s0 == 0, 'all') || (hasEnd && any(sEnd == 0, 'all'))
    s0 = s0 + 1;
    if hasEnd
        sEnd = sEnd + 1;
    end
end

nSamp = size(XEEG,1);

% duration in samples per trial
if numel(epochSec) == 1
    durSamp = repmat(round(epochSec * fs), numel(s0), 1);
else
    durSamp = round(epochSec(:) * fs);
    if numel(durSamp) ~= numel(s0)
        % resize conservatively
        durSamp = repmat(round(median(epochSec(isfinite(epochSec) & epochSec>0)) * fs), numel(s0), 1);
    end
end

% sanitize durations
badDur = ~isfinite(durSamp) | durSamp <= 0;
if any(badDur)
    dmed = median(durSamp(~badDur));
    if ~isfinite(dmed) || dmed <= 0
        dmed = round(1 * fs);
    end
    durSamp(badDur) = dmed;
end

% decide end sample
s1 = sEnd;
useEnd = hasEnd & isfinite(sEnd) & (sEnd >= s0);
if ~any(useEnd)
    s1 = s0 + durSamp - 1;
else
    s1(~useEnd) = s0(~useEnd) + durSamp(~useEnd) - 1;
end

% clamp and integerize
s0 = max(1, min(nSamp, floor(s0)));
s1 = max(1, min(nSamp, floor(s1)));

% ensure s1>=s0
swap = s1 < s0;
if any(swap)
    tmp = s0(swap);
    s0(swap) = s1(swap);
    s1(swap) = tmp;
end

% compute maxLen
lens = s1 - s0 + 1;
maxLen = max(lens(~isnan(lens) & isfinite(lens)));
if isempty(maxLen) || ~isfinite(maxLen) || maxLen < 1
    error('get_phase_epochs_trialinfo:NoValidEpoch', 'No valid epochs could be constructed for phase=%s', phase);
end

nTr = numel(s0);
nCh = size(XEEG,2);
ep = nan(nCh, maxLen, nTr);

for tr = 1:nTr
    a = s0(tr); b = s1(tr);
    if ~isfinite(a) || ~isfinite(b) || a < 1 || b < 1 || a > nSamp || b > nSamp
        continue;
    end
    x = XEEG(a:b, :)'; % [chan x time]
    L = size(x,2);
    ep(:, 1:L, tr) = x;
end

tvec = (0:maxLen-1) ./ fs;

% trial info
if istable(trialInfoAll)
    if ismember('epochType', trialInfoAll.Properties.VariableNames)
        m = strcmpi(trialInfoAll.epochType, phase);
        if any(m)
            ti = trialInfoAll(m,:);
        else
            ti = trialInfoAll;
        end
    else
        ti = trialInfoAll;
    end
else
    ti = table((1:nTr)', 'VariableNames', {'trial'});
end

if height(ti) ~= nTr
    % align to number of trials
    ti = ti(1:min(height(ti),nTr),:);
    if height(ti) < nTr
        ti2 = table((height(ti)+1:nTr)', 'VariableNames', {'trial'});
        ti = [ti; ti2];
    end
end

if ~ismember('epochType', ti.Properties.VariableNames)
    ti.epochType = repmat({phase}, height(ti), 1);
else
    ti.epochType(:) = {phase};
end

end

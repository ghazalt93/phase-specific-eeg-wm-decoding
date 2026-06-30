function [epochsAll, timeVecAll, trialInfoAll] = WM_BuildEpochs_FromTphases(EEG, T_phases, fs)

if nargin<3 || isempty(fs)
    if isfield(EEG,'srate') && ~isempty(EEG.srate), fs = double(EEG.srate);
    else, fs = 250;
    end
end

% Required start-sample columns:
reqStart = {'StimSample','MaintSample','RetrSample'};
for i=1:numel(reqStart)
    if ~ismember(reqStart{i}, T_phases.Properties.VariableNames)
        error('T_phases missing column: %s', reqStart{i});
    end
end

T_phases = ensure_end_column(T_phases, 'Stim', 0.60, fs);
T_phases = ensure_end_column(T_phases, 'Maint', 1.50, fs);
T_phases = ensure_end_column(T_phases, 'Retr', 2.00, fs);


% Condition
if ~ismember('Condition', T_phases.Properties.VariableNames)
    if ismember('TrigCond', T_phases.Properties.VariableNames)
        T_phases.Condition = double(T_phases.TrigCond);
    else
        error('T_phases missing Condition/TrigCond');
    end
end

X = EEG.data;  % EEGLAB style: [nChan x nSamples]
nChan = size(X,1);
nSamp = size(X,2);
nTr   = height(T_phases);

T_phases = ensure_end_column(T_phases, 'Stim', 0.60, fs);
T_phases = ensure_end_column(T_phases, 'Maint', 1.50, fs);
T_phases = ensure_end_column(T_phases, 'Retr', 2.00, fs);


% Condition
if ~ismember('Condition', T_phases.Properties.VariableNames)
    if ismember('TrigCond', T_phases.Properties.VariableNames)
        T_phases.Condition = double(T_phases.TrigCond);
    else
        error('T_phases missing Condition/TrigCond');
    end
end

[epochsAll.stim,  timeVecAll.stim ] = cut_phase('stim',  T_phases.StimSample,  T_phases.StimEnd);
[epochsAll.maint, timeVecAll.maint] = cut_phase('maint', T_phases.MaintSample, T_phases.MaintEnd);
[epochsAll.retr,  timeVecAll.retr ] = cut_phase('retr',  T_phases.RetrSample,  T_phases.RetrEnd);

% build trialInfoAll (stacked)
Tstim       = T_phases;  Tstim.epochType  = repmat("stim",  nTr, 1);
Tmaint      = T_phases;  Tmaint.epochType = repmat("maint", nTr, 1);
Tretr       = T_phases;  Tretr.epochType  = repmat("retr",  nTr, 1);
trialInfoAll = [Tstim; Tmaint; Tretr];

    function [E, tvec] = cut_phase(tag, S0, S1)
        S0 = double(S0(:)); S1 = double(S1(:));

        ok = ~isnan(S0) & ~isnan(S1) & (S0>=1) & (S1>=1) & (S0<=nSamp) & (S1<=nSamp) & (S1>=S0);
        if any(~ok)
           
           
        end

        L = S1 - S0 + 1;
        L(~ok) = NaN;

  
        maxL = max(L(~isnan(L)));
        if isempty(maxL) || maxL<=0
            E = [];
            tvec = [];
            return
        end

        E = nan(nChan, maxL, nTr, 'like', X);

        for ii=1:nTr
            if ~ok(ii), continue; end
            a = S0(ii); b = S1(ii);
            seg = X(:, a:b);
            E(:,1:size(seg,2),ii) = seg;
        end

        tvec = (0:maxL-1)./fs;
    end
end


function T = ensure_end_column(T, prefix, defaultSec, fs)
% Ensure <prefix>End exists in samples.
endName = [prefix 'End'];
startName = [prefix 'Sample'];
if ismember(endName, T.Properties.VariableNames)
    return;
end

% try duration columns
durNames = { [prefix 'DurSec'], [prefix 'Dur'], [prefix 'DurationSec'], [prefix 'Duration'] };
dur = [];
for i=1:numel(durNames)
    if ismember(durNames{i}, T.Properties.VariableNames)
        dur = T.(durNames{i});
        break;
    end
end

if isempty(dur)
    durSec = defaultSec;
    n = height(T);
    durSamp = repmat(round(durSec*fs), n, 1);
else
    if isnumeric(dur)
        durSec = double(dur(:));
    elseif iscell(dur)
        durSec = nan(numel(dur),1);
        for k=1:numel(dur)
            try
                durSec(k) = double(dur{k});
            catch
            end
        end
    else
        durSec = defaultSec;
        durSec = repmat(durSec, height(T), 1);
    end
    durSamp = round(durSec .* fs);
end

s0 = double(T.(startName));
T.(endName) = s0 + durSamp - 1;
end


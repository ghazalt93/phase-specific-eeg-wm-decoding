function [EEG, epochs, timeVec, trialInfo, cleanInfo, missingSubjects] = ...
    WM_LoadEEG_AndBuildEpochs(loaded_ws, edf_file, subj_name, cfg)

EEG = [];
epochs  = struct('stim',[],'maint',[],'retr',[]);
timeVec = struct('stim',[],'maint',[],'retr',[]);
trialInfo = table();
cleanInfo = struct();
missingSubjects = {};

if nargin < 4 || isempty(cfg), cfg = struct(); end

% ---- defaults ----
if ~isfield(cfg,'win') || isempty(cfg.win), cfg.win = struct(); end
if ~isfield(cfg.win,'stim')  || isempty(cfg.win.stim),  cfg.win.stim  = [0 0.5]; end
if ~isfield(cfg.win,'maint') || isempty(cfg.win.maint), cfg.win.maint = [0.5 2.5]; end
if ~isfield(cfg.win,'retr')  || isempty(cfg.win.retr),  cfg.win.retr  = [0 10]; end

if ~isfield(cfg,'bp') || isempty(cfg.bp), cfg.bp = [0.1 40]; end

if ~isfield(cfg,'nEEGChan') || isempty(cfg.nEEGChan), cfg.nEEGChan = 64; end

% bad channel policy
if ~isfield(cfg,'badChan') || isempty(cfg.badChan), cfg.badChan = struct(); end
if ~isfield(cfg.badChan,'mode') || isempty(cfg.badChan.mode), cfg.badChan.mode = "none"; end
if ~isfield(cfg.badChan,'stdFactor') || isempty(cfg.badChan.stdFactor), cfg.badChan.stdFactor = 5; end

% ---- backward compatibility ----
if isfield(cfg,'badChanStdFactor') && (~isfield(cfg.badChan,'stdFactor') || isempty(cfg.badChan.stdFactor))
    cfg.badChan.stdFactor = cfg.badChanStdFactor;
end
if isfield(cfg,'runICA') && cfg.runICA
    cfg.ica.mode = "run";
end

% ICA policy
if ~isfield(cfg,'ica') || isempty(cfg.ica), cfg.ica = struct(); end
if ~isfield(cfg.ica,'mode') || isempty(cfg.ica.mode), cfg.ica.mode = "none"; end
if ~isfield(cfg,'runICA') || isempty(cfg.runICA), cfg.runICA = false; end

% ---- trial->phase settings ----
if ~isfield(cfg,'trialToPhase') || isempty(cfg.trialToPhase), cfg.trialToPhase = struct(); end

%% ----------- LOAD / BUILD T_phases -----------
if isfield(loaded_ws,'T_phases') && ~isempty(loaded_ws.T_phases)
    T = loaded_ws.T_phases;
else
    if ~isfield(loaded_ws,'T_trials_p5') || ~isfield(loaded_ws,'T_triggers') || ...
       isempty(loaded_ws.T_trials_p5) || isempty(loaded_ws.T_triggers)
        warning('[SKIP] %s : neither T_phases nor (T_trials_p5 + T_triggers) exists', subj_name);
        missingSubjects = {subj_name};
        return
    end
    try
        [T, infoPh] = WM_TrialsToPhases(loaded_ws.T_trials_p5, loaded_ws.T_triggers, cfg.trialToPhase);
        cleanInfo.infoPh = infoPh;
    catch ME
        warning('[SKIP] %s : WM_TrialsToPhases failed -> %s', subj_name, ME.message);
        missingSubjects = {subj_name};
        return
    end
end

reqCols = {'StimSample','MaintSample','RetrSample','Condition'};
if ~all(ismember(reqCols, T.Properties.VariableNames))
    warning('[SKIP] %s : T_phases missing required columns', subj_name);
    missingSubjects = {subj_name};
    return
end

%% ----------- EDF check -----------
if ~isfile(edf_file)
    warning('[SKIP] %s : EDF not found -> %s', subj_name, edf_file);
    missingSubjects = {subj_name};
    return
end

try
    hdr = sopen(edf_file);
    nChanFile = hdr.NS;
    sclose(hdr);
catch ME
    warning('[SKIP] %s : cannot read EDF header -> %s', subj_name, ME.message);
    missingSubjects = {subj_name};
    return
end

fprintf('[EDF] %s | hdr.NS = %d channels\n', subj_name, nChanFile);

if nChanFile < cfg.nEEGChan
    warning('[SKIP] %s : EDF has only %d channels (<%d).', subj_name, nChanFile, cfg.nEEGChan);
    missingSubjects = {subj_name};
    return
end

%% ----------- Load EEG channels only (1..64) -----------
try
    try
        EEG = pop_biosig(edf_file, 'channels', 1:cfg.nEEGChan, 'importevent','off');
    catch
        EEG = pop_biosig(edf_file, 'channels', 1:cfg.nEEGChan);
    end
    EEG = eeg_checkset(EEG);
catch ME
    warning('[SKIP] %s : pop_biosig failed -> %s', subj_name, ME.message);
    missingSubjects = {subj_name};
    EEG = [];
    return
end

EEG.nbchan = size(EEG.data,1);
if EEG.nbchan ~= cfg.nEEGChan
    warning('[WARN] %s : loaded %d chans, expected %d (continuing).', subj_name, EEG.nbchan, cfg.nEEGChan);
end

%% ----------- PREPROCESSING -----------
try
    EEG = pop_eegfiltnew(EEG, cfg.bp(1), cfg.bp(2));
catch ME
    warning('[PREPROCESS] %s : filter failed -> %s', subj_name, ME.message);
end

% optional remove bad channels (OFF by default)
if string(cfg.badChan.mode) == "remove"
    try
        chanStd = std(double(EEG.data),0,2);
        thr = cfg.badChan.stdFactor * median(chanStd);
        badCh = find(chanStd > thr);
        if ~isempty(badCh)
            warning('[PREPROCESS] %s : Removing %d bad channels', subj_name, numel(badCh));
            EEG.data(badCh,:) = [];
            EEG.nbchan = size(EEG.data,1);
            if isfield(EEG,'chanlocs') && numel(EEG.chanlocs) >= max(badCh)
                EEG.chanlocs(badCh) = [];
            end
            EEG = eeg_checkset(EEG);
        end
    catch ME
        warning('[PREPROCESS] %s : bad-channel removal failed -> %s', subj_name, ME.message);
    end
end

% ICA: strongly recommended OFF for batch
doICA = (cfg.runICA == true) || (isfield(cfg,'ica') && isfield(cfg.ica,'mode') && string(cfg.ica.mode) ~= "none");
if doICA
    tICA = tic;
    try
        EEG = pop_runica(EEG,'extended',1,'interrupt','off');
        fprintf('[ICA] %s done in %.1f sec\n', subj_name, toc(tICA));
    catch ME
        warning('[ICA] %s failed -> %s (continuing without ICA)', subj_name, ME.message);
    end
end

maxS = max([T.StimSample; T.MaintSample; T.RetrSample], [], 'omitnan');
if ~isnan(maxS) && maxS > EEG.pnts
    warning('[ALIGN] %s : T_phases sample (%d) exceeds EEG.pnts (%d). Epochs may be empty.', ...
        subj_name, round(maxS), EEG.pnts);
end

%% ================= BUILD EPOCHS (ALIGN-SAFE) =================
win = cfg.win;
fs  = EEG.srate;

epochs.stim  = {};
epochs.maint = {};
epochs.retr  = {};

stim_src  = [];
maint_src = [];
retr_src  = [];

nTrials = height(T);

for i = 1:nTrials
    if ~isnan(T.StimSample(i))
        seg = T.StimSample(i) + round(win.stim*fs);
        if seg(1) > 0 && seg(2) <= EEG.pnts
            epochs.stim{end+1} = EEG.data(:, seg(1):seg(2)); %#ok<AGROW>
            stim_src(end+1,1)  = i; %#ok<AGROW>
        end
    end

    if ~isnan(T.MaintSample(i))
        seg = T.MaintSample(i) + round(win.maint*fs);
        if seg(1) > 0 && seg(2) <= EEG.pnts
            epochs.maint{end+1} = EEG.data(:, seg(1):seg(2)); %#ok<AGROW>
            maint_src(end+1,1)  = i; %#ok<AGROW>
        end
    end

    if ~isnan(T.RetrSample(i))
        seg = T.RetrSample(i) + round(win.retr*fs);
        if seg(1) > 0 && seg(2) <= EEG.pnts
            epochs.retr{end+1} = EEG.data(:, seg(1):seg(2)); %#ok<AGROW>
            retr_src(end+1,1)  = i; %#ok<AGROW>
        end
    end
end

if ~isempty(epochs.stim),  epochs.stim  = cat(3, epochs.stim{:});  end
if ~isempty(epochs.maint), epochs.maint = cat(3, epochs.maint{:}); end
if ~isempty(epochs.retr),  epochs.retr  = cat(3, epochs.retr{:});  end

%% ================= BUILD ALIGNED trialInfo =================
trialInfoAll = table();

if ~isempty(stim_src)
    tmp = T(stim_src,:);
    tmp.epochType = repmat("stim", height(tmp), 1);
    tmp.srcRow    = stim_src;
    trialInfoAll  = [trialInfoAll; tmp];
end
if ~isempty(maint_src)
    tmp = T(maint_src,:);
    tmp.epochType = repmat("maint", height(tmp), 1);
    tmp.srcRow    = maint_src;
    trialInfoAll  = [trialInfoAll; tmp];
end
if ~isempty(retr_src)
    tmp = T(retr_src,:);
    tmp.epochType = repmat("retr", height(tmp), 1);
    tmp.srcRow    = retr_src;
    trialInfoAll  = [trialInfoAll; tmp];
end

trialInfo = trialInfoAll;

%% ================= TIME VECTORS =================
if ~isempty(epochs.stim)
    timeVec.stim  = linspace(win.stim(1),  win.stim(2),  size(epochs.stim,2));
end
if ~isempty(epochs.maint)
    timeVec.maint = linspace(win.maint(1), win.maint(2), size(epochs.maint,2));
end
if ~isempty(epochs.retr)
    timeVec.retr  = linspace(win.retr(1),  win.retr(2),  size(epochs.retr,2));
end

%% ================= CLEAN INFO =================
cleanInfo.subject   = subj_name;
cleanInfo.nTrials   = nTrials;
cleanInfo.nChannels = EEG.nbchan;
cleanInfo.fs        = fs;

nStim  = 0; if ~isempty(epochs.stim),  nStim  = size(epochs.stim,3); end
nMaint = 0; if ~isempty(epochs.maint), nMaint = size(epochs.maint,3); end
nRetr  = 0; if ~isempty(epochs.retr),  nRetr  = size(epochs.retr,3); end
cleanInfo.nEpochs = struct('stim',nStim,'maint',nMaint,'retr',nRetr);

fprintf('[OK] %s | stim=%d | maint=%d | retr=%d\n', subj_name, nStim, nMaint, nRetr);

end

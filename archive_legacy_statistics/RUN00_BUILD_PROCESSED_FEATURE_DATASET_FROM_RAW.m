%% RUN00_BUILD_PROCESSED_FEATURE_DATASET_FROM_RAW.m
% 

clear; clc;

cfgBase = get_project_config();

fprintf('\n============================================================\n');
fprintf('RUN00: Build processed EEG feature dataset from raw files\n');
fprintf('Project root: %s\n', cfgBase.projectRoot);
fprintf('Data root   : %s\n', cfgBase.dataRoot);
fprintf('Output root : %s\n', cfgBase.outputRoot);
fprintf('============================================================\n');

%% ---------------- Paths ----------------
rawSubjectsRoot = fullfile(cfgBase.dataRoot, 'Subjects');
outRoot = fullfile(cfgBase.outputRoot, 'new_article_outputs');
outDir  = fullfile(outRoot, '_wm_ml');

if ~exist(outDir, 'dir'), mkdir(outDir); end

if ~exist(rawSubjectsRoot, 'dir')
    error(['Raw subject folder not found:\n%s\n\n', ...
           'Raw EEG/behavioral files are not included in the public repository by default. ', ...
           'Place subject folders under data/Subjects or edit rawSubjectsRoot in this script.'], rawSubjectsRoot);
end

diary(fullfile(outDir, 'RUN00_build_processed_feature_dataset_log.txt'));

%% ---------------- Feature-extraction configuration ----------------
cfgFE = struct();

% Article-level channel and preprocessing assumptions
cfgFE.nEEGChan = 64;
cfgFE.bp = [0.1 40];
cfgFE.runICA = false;
cfgFE.ica.mode = "none";
cfgFE.badChan.mode = "none";

% Phase windows relative to phase anchors
cfgFE.win.stim  = [-0.2 0.6];
cfgFE.win.maint = [-0.2 1.0];
cfgFE.win.retr  = [-0.2 0.6];

% Feature families
cfgFE.do_freq       = true;
cfgFE.do_stats      = true;
cfgFE.do_temporal   = true;
cfgFE.do_complexity = true;

% Frequency bands
cfgFE.freqBands.delta = [1 4];
cfgFE.freqBands.theta = [4 8];
cfgFE.freqBands.alpha = [8 13];
cfgFE.freqBands.beta  = [13 30];
cfgFE.freqBands.gamma = [30 45];

% Complexity defaults
cfgFE.comp.ds = 2;
cfgFE.comp.maxLen = 1200;
cfgFE.comp.binarize = 'median';
cfgFE.comp.normalize = true;

% PSD summaries can generate many files; keep true for QC, false for speed.
cfgFE.psd.enable = true;

phaseList = {'stim','maint','retr'};

%% ---------------- Discover subject folders ----------------
Dsub = dir(rawSubjectsRoot);
Dsub = Dsub([Dsub.isdir]);
Dsub = Dsub(~ismember({Dsub.name}, {'.','..'}));

if isempty(Dsub)
    error('No subject folders found under:\n%s', rawSubjectsRoot);
end

fprintf('\nFound %d subject folders.\n', numel(Dsub));

allRows = {};
statusRows = {};

%% ---------------- Main subject/run loop ----------------
for si = 1:numel(Dsub)

    subjName = string(Dsub(si).name);
    subjDir = fullfile(rawSubjectsRoot, Dsub(si).name);

    fprintf('\n============================================================\n');
    fprintf('Subject %d/%d: %s\n', si, numel(Dsub), subjName);
    fprintf('Folder: %s\n', subjDir);
    fprintf('============================================================\n');

    matFiles = local_find_workspace_files(subjDir);
    edfFiles = local_find_edf_files(subjDir);

    if isempty(matFiles)
        warning('[SKIP] %s: no usable .mat workspace files found.', subjName);
        statusRows{end+1,1} = local_status(subjName, NaN, "", "", "skip_no_workspace", 0, 0); %#ok<AGROW>
        continue;
    end
    if isempty(edfFiles)
        warning('[SKIP] %s: no EDF files found.', subjName);
        statusRows{end+1,1} = local_status(subjName, NaN, "", "", "skip_no_edf", 0, 0); %#ok<AGROW>
        continue;
    end

    for mi = 1:numel(matFiles)

        matPath = fullfile(matFiles(mi).folder, matFiles(mi).name);
        runNum = local_detect_run_number(matFiles(mi).name, mi);
        edfPath = local_match_edf_for_workspace(matFiles(mi), edfFiles, runNum);

        fprintf('\n--- Subject=%s | run=%d ---\n', subjName, runNum);
        fprintf('MAT: %s\n', matPath);
        fprintf('EDF: %s\n', edfPath);

        try
            loaded_ws = load(matPath);

            cfgLoad = cfgFE;
            cfgLoad.subjectName = char(subjName);
            cfgLoad.runTag = sprintf('run%d', runNum);

            [EEG, epochs, timeVec, trialInfo, cleanInfo, missingSubjects] = ...
                WM_LoadEEG_AndBuildEpochs(loaded_ws, edfPath, char(subjName), cfgLoad);

            if ~isempty(missingSubjects) || isempty(EEG) || isempty(fieldnames(epochs))
                warning('[SKIP] %s run %d: loader returned missing/empty output.', subjName, runNum);
                statusRows{end+1,1} = local_status(subjName, runNum, matPath, edfPath, ...
                    "skip_loader_empty", 0, 0); %#ok<AGROW>
                continue;
            end

            subjOutDir = fullfile(outDir, char(subjName), sprintf('run%d', runNum));
            if ~exist(subjOutDir, 'dir'), mkdir(subjOutDir); end

            for pi = 1:numel(phaseList)

                phaseName = phaseList{pi};

                if ~isfield(epochs, phaseName) || isempty(epochs.(phaseName))
                    warning('[SKIP] %s run %d phase %s: empty epochs.', subjName, runNum, phaseName);
                    statusRows{end+1,1} = local_status(subjName, runNum, matPath, edfPath, ...
                        "skip_empty_" + string(phaseName), 0, 0); %#ok<AGROW>
                    continue;
                end
                if ~isfield(timeVec, phaseName) || isempty(timeVec.(phaseName))
                    warning('[SKIP] %s run %d phase %s: empty time vector.', subjName, runNum, phaseName);
                    statusRows{end+1,1} = local_status(subjName, runNum, matPath, edfPath, ...
                        "skip_empty_time_" + string(phaseName), 0, 0); %#ok<AGROW>
                    continue;
                end

                trialInfoPhase = local_trialinfo_for_phase(trialInfo, phaseName, size(epochs.(phaseName),3));

                cfgPhase = cfgFE;
                cfgPhase.phase = phaseName;
                cfgPhase.subjectName = char(subjName);
                cfgPhase.runTag = sprintf('run%d', runNum);

                [featAll, trialInfoOut, meta] = WM_ExtractAllFeatures_patched( ...
                    EEG, epochs, timeVec, trialInfoPhase, cfgPhase); %#ok<NASGU>

                if isempty(featAll) || height(featAll) == 0 || width(featAll) == 0
                    warning('[SKIP] %s run %d phase %s: no features returned.', subjName, runNum, phaseName);
                    statusRows{end+1,1} = local_status(subjName, runNum, matPath, edfPath, ...
                        "skip_no_features_" + string(phaseName), 0, 0); %#ok<AGROW>
                    continue;
                end

                % Align trialInfoOut and features
                n = min(height(featAll), height(trialInfoOut));
                featAll = featAll(1:n,:);
                trialInfoOut = trialInfoOut(1:n,:);

                Meta = local_make_meta_table(trialInfoOut, subjName, runNum, phaseName, n);
                Tphase = [Meta featAll];

                % Ensure unique variable names after concatenation
                Tphase.Properties.VariableNames = matlab.lang.makeUniqueStrings(Tphase.Properties.VariableNames);

                allRows{end+1,1} = Tphase; %#ok<AGROW>

                outCsv = fullfile(subjOutDir, sprintf('%s_run%d_%s_features.csv', subjName, runNum, phaseName));
                writetable(Tphase, outCsv);

                statusRows{end+1,1} = local_status(subjName, runNum, matPath, edfPath, ...
                    "ok_" + string(phaseName), height(Tphase), width(featAll)); %#ok<AGROW>

                fprintf('[OK] %s run %d phase %s | rows=%d | features=%d\n', ...
                    subjName, runNum, phaseName, height(Tphase), width(featAll));
            end

            try
                save(fullfile(subjOutDir, sprintf('%s_run%d_epochs_cleanInfo.mat', subjName, runNum)), ...
                    'cleanInfo', 'timeVec', '-v7.3');
            catch ME
                warning('[SAVE] Could not save cleanInfo/timeVec for %s run %d: %s', subjName, runNum, ME.message);
            end

        catch ME
            warning('[FAILED] %s run %d: %s', subjName, runNum, ME.message);
            statusRows{end+1,1} = local_status(subjName, runNum, matPath, edfPath, ...
                "failed_" + string(ME.message), 0, 0); %#ok<AGROW>
        end
    end
end

%% ---------------- Combine and save final dataset ----------------
if isempty(allRows)
    Status = local_vertcat_or_empty(statusRows);
    writetable(Status, fullfile(outDir, 'RUN00_feature_extraction_status.csv'));
    save(fullfile(outDir, 'RUN00_feature_extraction_status.mat'), 'Status', 'cfgFE', '-v7.3');
    diary off;
    error('No feature rows were created. Check raw-data paths, EDF files, and workspace contents.');
end

allRows = local_harmonize_tables_for_vertcat(allRows);
DS = vertcat(allRows{:});
T = DS; %#ok<NASGU>

Status = local_vertcat_or_empty(statusRows);

writetable(DS, fullfile(outDir, 'dataset.csv'));
writetable(Status, fullfile(outDir, 'RUN00_feature_extraction_status.csv'));

save(fullfile(outDir, 'dataset.mat'), 'DS', 'T', 'Status', 'cfgFE', '-v7.3');

fprintf('\n============================================================\n');
fprintf('RUN00 finished.\n');
fprintf('Final dataset rows: %d\n', height(DS));
fprintf('Final dataset cols: %d\n', width(DS));
fprintf('Saved dataset:\n%s\n', fullfile(outDir, 'dataset.mat'));
fprintf('============================================================\n');

diary off;

%% ========================================================================
%% Local helper functions
%% ========================================================================

function files = local_find_workspace_files(subjDir)
files = dir(fullfile(subjDir, '**', '*.mat'));
if isempty(files), return; end

keep = true(numel(files),1);
for i = 1:numel(files)
    nm = lower(files(i).name);
    fp = lower(fullfile(files(i).folder, files(i).name));

    % Exclude obvious non-task / intermediate files.
    bad = contains(nm, 'eye') || contains(nm, 'open') || contains(nm, 'close') || ...
          contains(nm, 'rest') || contains(nm, 'debug') || contains(nm, 'feature') || ...
          contains(nm, 'dataset') || contains(fp, '_wm_ml') || contains(fp, 'outputs');

    keep(i) = ~bad;
end
files = files(keep);
end

function files = local_find_edf_files(subjDir)
files = dir(fullfile(subjDir, '**', '*.edf'));
if isempty(files)
    files = dir(fullfile(subjDir, '**', '*.EDF'));
end

if isempty(files), return; end

keep = true(numel(files),1);
for i = 1:numel(files)
    nm = lower(files(i).name);
    keep(i) = ~(contains(nm, 'eye') || contains(nm, 'open') || contains(nm, 'close') || contains(nm, 'rest'));
end
files = files(keep);
end

function runNum = local_detect_run_number(name, fallback)
name = lower(char(name));
tok = regexp(name, 'run[_-]?(\d+)', 'tokens', 'once');
if isempty(tok), tok = regexp(name, '(?:^|[_-])r[_-]?(\d+)(?:[_-]|\.|$)', 'tokens', 'once'); end
if isempty(tok), tok = regexp(name, '(\d+)', 'tokens', 'once'); end
if isempty(tok)
    runNum = fallback;
else
    runNum = str2double(tok{1});
    if ~isfinite(runNum) || runNum < 1
        runNum = fallback;
    end
end
end

function edfPath = local_match_edf_for_workspace(matFile, edfFiles, runNum)
edfPath = fullfile(edfFiles(1).folder, edfFiles(1).name);
if numel(edfFiles) == 1, return; end

matName = lower(matFile.name);
matRunStr = sprintf('%d', runNum);

% Prefer EDF in the same folder.
sameFolder = strcmp({edfFiles.folder}', matFile.folder);
candIdx = find(sameFolder);
if isempty(candIdx)
    candIdx = 1:numel(edfFiles);
end

% Prefer run-matched EDF among candidates.
for k = candIdx(:)'
    nm = lower(edfFiles(k).name);
    if contains(nm, ['run' matRunStr]) || contains(nm, ['run_' matRunStr]) || ...
       contains(nm, ['r' matRunStr]) || contains(nm, ['_' matRunStr])
        edfPath = fullfile(edfFiles(k).folder, edfFiles(k).name);
        return;
    end
end

% Prefer any EDF with a token from MAT name.
[~, baseMat] = fileparts(matName);
for k = candIdx(:)'
    [~, baseEdf] = fileparts(lower(edfFiles(k).name));
    if contains(baseEdf, baseMat) || contains(baseMat, baseEdf)
        edfPath = fullfile(edfFiles(k).folder, edfFiles(k).name);
        return;
    end
end

% Fallback: first candidate.
edfPath = fullfile(edfFiles(candIdx(1)).folder, edfFiles(candIdx(1)).name);
end

function ti = local_trialinfo_for_phase(trialInfo, phaseName, nEpochs)
if istable(trialInfo) && ~isempty(trialInfo)
    if ismember('epochType', trialInfo.Properties.VariableNames)
        m = strcmpi(string(trialInfo.epochType), string(phaseName));
        ti = trialInfo(m,:);
    else
        ti = trialInfo;
    end
else
    ti = table();
end

if ~istable(ti) || isempty(ti)
    ti = table((1:nEpochs)', 'VariableNames', {'TrialIndex'});
end

if height(ti) ~= nEpochs
    n = min(height(ti), nEpochs);
    ti = ti(1:n,:);
end
end

function Meta = local_make_meta_table(trialInfo, subjName, runNum, phaseName, n)
Meta = table();
Meta.Subject = repmat(string(subjName), n, 1);
Meta.Run = repmat(runNum, n, 1);
Meta.Phase = repmat(string(phaseName), n, 1);

% Condition is required by downstream STEP100/104.
cond = local_pick_col(trialInfo, {'Condition','condition','yCondition','Label','label','TrigCond'});
if isempty(cond)
    Meta.Condition = repmat(missing, n, 1);
else
    Meta.Condition = local_to_string_column(cond, n);
end

% Useful non-leakage metadata.
trialNum = local_pick_col(trialInfo, {'TrialNum','trialNum','Trial','trial','TrialIndex','trialIndex'});
if isempty(trialNum)
    Meta.TrialNum = (1:n)';
else
    Meta.TrialNum = local_to_numeric_column(trialNum, n);
end

correct = local_pick_col(trialInfo, {'Correct','correct','yCorrect','Accuracy','accuracy'});
if ~isempty(correct)
    Meta.Correct = local_to_numeric_column(correct, n);
end

rt = local_pick_col(trialInfo, {'RT','rt','ReactionTime','responseTime','ResponseTime'});
if ~isempty(rt)
    Meta.RT = local_to_numeric_column(rt, n);
end
end

function col = local_pick_col(T, candidates)
col = [];
if ~istable(T), return; end
names = T.Properties.VariableNames;
for i = 1:numel(candidates)
    idx = strcmpi(names, candidates{i});
    if any(idx)
        col = T.(names{find(idx,1)});
        return;
    end
end
end

function x = local_to_string_column(col, n)
try
    if iscategorical(col), x = string(col); else, x = string(col); end
    x = x(:);
catch
    x = repmat(missing, n, 1);
end
x = local_pad_or_crop_string(x, n);
end

function x = local_to_numeric_column(col, n)
try
    if isnumeric(col) || islogical(col)
        x = double(col(:));
    else
        x = str2double(string(col(:)));
    end
catch
    x = nan(n,1);
end
x = local_pad_or_crop_numeric(x, n);
end

function x = local_pad_or_crop_string(x, n)
if numel(x) >= n
    x = x(1:n);
else
    x(end+1:n,1) = missing;
end
end

function x = local_pad_or_crop_numeric(x, n)
if numel(x) >= n
    x = x(1:n);
else
    x(end+1:n,1) = NaN;
end
end

function S = local_status(subjName, runNum, matPath, edfPath, status, nRows, nFeatures)
S = table(string(subjName), runNum, string(matPath), string(edfPath), string(status), nRows, nFeatures, ...
    'VariableNames', {'Subject','Run','MatFile','EdfFile','Status','N_rows','N_features'});
end

function T = local_vertcat_or_empty(C)
if isempty(C), T = table(); else, T = vertcat(C{:}); end
end

function tables = local_harmonize_tables_for_vertcat(tables)
allVars = {};
for i = 1:numel(tables)
    allVars = union(allVars, tables{i}.Properties.VariableNames, 'stable');
end

for i = 1:numel(tables)
    T = tables{i};
    vars = T.Properties.VariableNames;
    for v = 1:numel(allVars)
        name = allVars{v};
        if ~any(strcmp(vars, name))
            T.(name) = local_missing_column(T, name);
        end
    end
    T = T(:, allVars);
    tables{i} = T;
end
end

function col = local_missing_column(T, name)
n = height(T);
numericNames = {'Run','TrialNum','Correct','RT'};
if any(strcmp(name, numericNames)) || startsWith(string(name), "WM_") || ...
        contains(string(name), "BP_") || contains(string(name), "RBP_") || ...
        contains(string(name), "Mean_") || contains(string(name), "Std_") || ...
        contains(string(name), "RMS_") || contains(string(name), "LZC_")
    col = NaN(n,1);
else
    col = repmat(missing, n, 1);
end
end

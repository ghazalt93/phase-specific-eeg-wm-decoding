%% RUN06_STEP104_collect_loso_decoding_outputs.m

clear; clc;

cfgBase = get_project_config();

ROOT = fullfile(cfgBase.outputRoot, 'new_article_outputs');

step104OutDir = fullfile(ROOT, ...
    '_wm_STEP104_FisherScore_CHANNEL_ROI_1to22_4classifiers_noLeakage');

collectOutDir = fullfile(ROOT, ...
    '_wm_STEP104_FINAL_collected_outputs');

if ~exist(step104OutDir, 'dir')
    error('STEP104 output folder not found:\n%s\nRun RUN05 first.', step104OutDir);
end

if ~exist(collectOutDir, 'dir')
    mkdir(collectOutDir);
end

fprintf('\n=== RUN06 / STEP104 COLLECT OUTPUTS ONLY ===\n');
fprintf('Reading existing STEP104 task outputs from:\n%s\n', step104OutDir);
fprintf('Saving collected outputs to:\n%s\n', collectOutDir);

diaryFile = fullfile(collectOutDir, 'RUN06_collect_loso_decoding_outputs_log.txt');
if exist(diaryFile, 'file')
    delete(diaryFile);
end
diary(diaryFile);

fprintf('\n=== RUN06 / STEP104 COLLECT OUTPUTS ONLY ===\n');
fprintf('Started: %s\n', datestr(now));
fprintf('Input folder : %s\n', step104OutDir);
fprintf('Output folder: %s\n', collectOutDir);

%% ---------------- Scan STEP104 task folders ----------------

D = dir(step104OutDir);
D = D([D.isdir]);
D = D(~ismember({D.name}, {'.','..'}));

allResults = {};
allPredictions = {};
allSelected = {};
statusRows = {};

for i = 1:numel(D)

    subDir = fullfile(step104OutDir, D(i).name);

    resFiles  = dir(fullfile(subDir, 'STEP104_*_results.csv'));
    predFiles = dir(fullfile(subDir, 'STEP104_*_predictions.csv'));
    selFiles  = dir(fullfile(subDir, 'STEP104_*_selected_features.csv'));

    fprintf('\n[%d/%d] %s\n', i, numel(D), D(i).name);
    fprintf('  results=%d | predictions=%d | selected=%d\n', ...
        numel(resFiles), numel(predFiles), numel(selFiles));

    % Results
    for k = 1:numel(resFiles)
        f = fullfile(subDir, resFiles(k).name);
        [T, statusMsg] = local_safe_readtable(f);
        statusRows{end+1,1} = table(string(f), "results", string(statusMsg), height(T), ...
            'VariableNames', {'File','Type','Status','N_rows'}); %#ok<AGROW>
        if ~isempty(T)
            allResults{end+1,1} = T; %#ok<AGROW>
        end
    end

    % Predictions
    for k = 1:numel(predFiles)
        f = fullfile(subDir, predFiles(k).name);
        [T, statusMsg] = local_safe_readtable(f);
        statusRows{end+1,1} = table(string(f), "predictions", string(statusMsg), height(T), ...
            'VariableNames', {'File','Type','Status','N_rows'}); %#ok<AGROW>
        if ~isempty(T)
            allPredictions{end+1,1} = T; %#ok<AGROW>
        end
    end

    % Selected features
    for k = 1:numel(selFiles)
        f = fullfile(subDir, selFiles(k).name);
        [T, statusMsg] = local_safe_readtable(f);
        statusRows{end+1,1} = table(string(f), "selected_features", string(statusMsg), height(T), ...
            'VariableNames', {'File','Type','Status','N_rows'}); %#ok<AGROW>
        if ~isempty(T)
            allSelected{end+1,1} = T; %#ok<AGROW>
        end
    end
end

%% ---------------- Combine tables ----------------

AllResults = local_vertcat_harmonized(allResults);
AllPredictions = local_vertcat_harmonized(allPredictions);
AllSelected = local_vertcat_harmonized(allSelected);

if isempty(statusRows)
    Status = table();
else
    Status = vertcat(statusRows{:});
end

if ~isempty(AllResults) && all(ismember({'BalAcc','AUC','Accuracy'}, AllResults.Properties.VariableNames))
    AllResults = sortrows(AllResults, {'BalAcc','AUC','Accuracy'}, {'descend','descend','descend'});
elseif ~isempty(AllResults) && all(ismember({'BalancedAccuracy','AUC','Accuracy'}, AllResults.Properties.VariableNames))
    AllResults = sortrows(AllResults, {'BalancedAccuracy','AUC','Accuracy'}, {'descend','descend','descend'});
end

%% ---------------- Save final collected outputs ----------------

outResults = fullfile(collectOutDir, ...
    'STEP104_FINAL_FisherScore_CHANNEL_ROI_featureCountCurve_results.csv');

outPredictions = fullfile(collectOutDir, ...
    'STEP104_FINAL_FisherScore_CHANNEL_ROI_featureCountCurve_predictions.csv');

outSelected = fullfile(collectOutDir, ...
    'STEP104_FINAL_FisherScore_CHANNEL_ROI_selected_features_by_fold.csv');

outStatus = fullfile(collectOutDir, ...
    'STEP104_FINAL_collect_status.csv');

outMAT = fullfile(collectOutDir, ...
    'STEP104_FINAL_collected_outputs.mat');

writetable(AllResults, outResults);
writetable(AllPredictions, outPredictions);
writetable(AllSelected, outSelected);
writetable(Status, outStatus);

save(outMAT, 'AllResults', 'AllPredictions', 'AllSelected', 'Status', ...
    'ROOT', 'step104OutDir', 'collectOutDir', '-v7.3');

fprintf('\n=== RUN06 collect finished ===\n');
fprintf('Results    : %s\n', outResults);
fprintf('Predictions: %s\n', outPredictions);
fprintf('Selected   : %s\n', outSelected);
fprintf('Status     : %s\n', outStatus);
fprintf('MAT        : %s\n', outMAT);
fprintf('Log        : %s\n', diaryFile);

if ~isempty(Status)
    fprintf('\nRead status summary:\n');
    try
        disp(groupsummary(Status, {'Type','Status'}));
    catch
        disp(Status);
    end
end

diary off;

%% ===================== LOCAL FUNCTIONS =====================

function [T, statusMsg] = local_safe_readtable(f)
    T = table();
    statusMsg = "ok";
    try
        if exist(f, 'file') ~= 2
            statusMsg = "missing";
            return;
        end
        T = readtable(f, 'PreserveVariableNames', true);
    catch ME
        T = table();
        statusMsg = "failed_" + string(ME.message);
    end
end

function Tout = local_vertcat_harmonized(listTables)
    if isempty(listTables)
        Tout = table();
        return;
    end

    % Remove empty tables.
    keep = true(numel(listTables),1);
    for i = 1:numel(listTables)
        keep(i) = istable(listTables{i}) && width(listTables{i}) > 0;
    end
    listTables = listTables(keep);

    if isempty(listTables)
        Tout = table();
        return;
    end

    % Union of variable names in first-seen order.
    allVars = {};
    for i = 1:numel(listTables)
        v = listTables{i}.Properties.VariableNames;
        for j = 1:numel(v)
            if ~ismember(v{j}, allVars)
                allVars{end+1} = v{j}; %#ok<AGROW>
            end
        end
    end

    % Add missing variables and align order.
    for i = 1:numel(listTables)
        T = listTables{i};

        for j = 1:numel(allVars)
            vn = allVars{j};
            if ~ismember(vn, T.Properties.VariableNames)
                T.(vn) = local_missing_column(height(T));
            end
        end

        T = T(:, allVars);
        T = local_convert_stringlike_columns(T);
        listTables{i} = T;
    end

    try
        Tout = vertcat(listTables{:});
    catch
        % Last-resort robust conversion: convert all nonnumeric columns to string.
        for i = 1:numel(listTables)
            T = listTables{i};
            for j = 1:width(T)
                x = T.(j);
                if ~(isnumeric(x) || islogical(x))
                    try
                        T.(j) = string(x);
                    catch
                        T.(j) = repmat("", height(T), 1);
                    end
                end
            end
            listTables{i} = T;
        end
        Tout = vertcat(listTables{:});
    end
end

function col = local_missing_column(n)
    col = strings(n,1);
end

function T = local_convert_stringlike_columns(T)
    for j = 1:width(T)
        x = T.(j);
        if iscellstr(x) || iscategorical(x) || ischar(x)
            T.(j) = string(x);
        end
    end
end

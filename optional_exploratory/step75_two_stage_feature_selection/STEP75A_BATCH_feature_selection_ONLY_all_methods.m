%% STEP75A_BATCH_feature_selection_ONLY_all_methods.m

cfg = get_project_config();

clear; clc;

%% ===================== BASE CONFIG =====================

baseCfg = WM75_config();

% Correct project root
baseCfg.ROOT = cfg.outputRoot;

% Feature-selection methods to run
featureMethodList = { ...
    'ranksum'
    'fisher'
    'mrmr'
    'relieff'
    'pca'
};

% Parameters
baseCfg.Kfeat = 100;
baseCfg.nPC = 50;
baseCfg.pcaVarianceToKeep = 95;

% Batch dimensions
runList   = [1 2 3];
phaseList = {'stim','maint','retr'};

% Pairwise binary contrasts between the three conditions
contrastList = { ...
    'orientation_vs_conjunction'
    'color_vs_orientation'
    'color_vs_conjunction'
};

% If you also want one-vs-rest contrasts, uncomment this block instead:
% contrastList = { ...
%     'orientation_vs_conjunction'
%     'color_vs_orientation'
%     'color_vs_conjunction'
%     'conjunction_vs_notConjunction'
%     'orientation_vs_notOrientation'
%     'color_vs_notColor'
% };

%% ===================== OUTPUT =====================

batchOutDir = fullfile(baseCfg.ROOT, '_wm_STEP75_two_stage_FS_then_classifier', ...
    '_BATCH_stage1_feature_selection_ONLY_all_methods');

if ~exist(batchOutDir, 'dir')
    mkdir(batchOutDir);
end

diary(fullfile(batchOutDir, 'STEP75A_BATCH_feature_selection_ONLY_all_methods_log.txt'));

fprintf('\n=== STEP75A BATCH: Feature selection ONLY, all methods ===\n');
fprintf('Started: %s\n', datestr(now));
fprintf('ROOT: %s\n', baseCfg.ROOT);
fprintf('Kfeat=%d | nPC=%d | PCA variance=%g%%\n', ...
    baseCfg.Kfeat, baseCfg.nPC, baseCfg.pcaVarianceToKeep);

fprintf('Methods:\n');
disp(featureMethodList(:));

allSummary = {};
allErrors = {};

%% ===================== MAIN LOOP =====================

job = 0;
nJobs = numel(featureMethodList) * numel(runList) * numel(phaseList) * numel(contrastList);

for m = 1:numel(featureMethodList)

    for c = 1:numel(contrastList)

        for r = 1:numel(runList)

            for p = 1:numel(phaseList)

                job = job + 1;

                cfg = baseCfg;
                cfg.featureMethod = featureMethodList{m};
                cfg.contrastName = contrastList{c};
                cfg.run = runList(r);
                cfg.phase = phaseList{p};

                cfg = local_refresh_step75_paths(cfg);

                fprintf('\n\n============================================================\n');
                fprintf('JOB %d/%d\n', job, nJobs);
                fprintf('FS method: %s\n', cfg.featureMethod);
                fprintf('Contrast: %s | Run: %d | Phase: %s\n', ...
                    cfg.contrastName, cfg.run, cfg.phase);
                fprintf('============================================================\n');

                try
                    %% ---------- Stage 1: Feature selection only ----------

                    if ~exist(cfg.fsOutDir, 'dir')
                        mkdir(cfg.fsOutDir);
                    end

                    [X, y, subj, featureNames, DataInfo] = STEP75_load_binary_data(cfg);

                    fprintf('[LOAD] N=%d | N_A=%d | N_B=%d | Subjects=%d | Features=%d\n', ...
                        numel(y), sum(y=='A'), sum(y=='B'), numel(unique(subj)), size(X,2));

                    [Xz, prep] = STEP75_preprocess_all_for_FS(X);

                    [FSmodel, Ranking] = STEP75_feature_selection_full(Xz, y, featureNames, cfg);

                    outRankCSV = fullfile(cfg.fsOutDir, 'STEP75A_feature_ranking.csv');
                    outTopCSV  = fullfile(cfg.fsOutDir, 'STEP75A_top_selected_features.csv');
                    outFSMAT   = cfg.selectedFeatureFile;

                    writetable(Ranking, outRankCSV);

                    if isfield(FSmodel,'selectedIdx') && ~isempty(FSmodel.selectedIdx)
                        nTop = min(height(Ranking), numel(FSmodel.selectedIdx));
                        TopSelected = Ranking(1:nTop, :);
                    else
                        TopSelected = Ranking;
                    end

                    writetable(TopSelected, outTopCSV);

                    save(outFSMAT, 'FSmodel', 'Ranking', 'prep', 'DataInfo', 'cfg', '-v7.3');

                    fprintf('[SAVED FS] %s\n', outFSMAT);

                    %% ---------- Summary row ----------

                    nSelected = NaN;
                    if isfield(FSmodel, 'selectedIdx')
                        nSelected = numel(FSmodel.selectedIdx);
                    elseif isfield(FSmodel, 'nComp')
                        nSelected = FSmodel.nComp;
                    end

                    topName = "";
                    topScore = NaN;

                    if ~isempty(Ranking) && height(Ranking) >= 1
                        varNames = Ranking.Properties.VariableNames;
                        if any(strcmp(varNames, 'Feature'))
                            topName = string(Ranking.Feature(1));
                        elseif any(strcmp(varNames, 'PCName'))
                            topName = string(Ranking.PCName(1));
                        end

                        scoreCols = {'RankAUC_score','FisherScore','MRMRScore','ReliefFWeight','ExplainedVariancePercent'};
                        for sc = 1:numel(scoreCols)
                            if any(strcmp(varNames, scoreCols{sc}))
                                topScore = Ranking.(scoreCols{sc})(1);
                                break;
                            end
                        end
                    end

                    SummaryOne = table( ...
                        string(cfg.featureMethod), string(cfg.contrastName), cfg.run, string(cfg.phase), ...
                        cfg.Kfeat, cfg.nPC, ...
                        DataInfo.N_total, DataInfo.N_A, DataInfo.N_B, DataInfo.N_subjects, DataInfo.N_features, ...
                        nSelected, topName, topScore, ...
                        string(outFSMAT), string('ok'), ...
                        'VariableNames', {'FeatureMethod','Contrast','Run','Phase','Kfeat','nPC', ...
                        'N_total','N_A','N_B','N_subjects','N_features_input', ...
                        'N_selected_or_components','TopFeatureOrComponent','TopScore', ...
                        'SavedFeatureFile','Status'});

                    allSummary{end+1,1} = SummaryOne; %#ok<AGROW>

                    fprintf('[DONE FS] %s | %s | run=%d | phase=%s | selected=%g | top=%s\n', ...
                        cfg.featureMethod, cfg.contrastName, cfg.run, cfg.phase, nSelected, topName);

                catch ME

                    warning('[FAILED FS] method=%s | %s | run=%d | phase=%s | %s', ...
                        cfg.featureMethod, cfg.contrastName, cfg.run, cfg.phase, ME.message);

                    ErrOne = table( ...
                        string(cfg.featureMethod), string(cfg.contrastName), cfg.run, string(cfg.phase), ...
                        string(ME.message), ...
                        'VariableNames', {'FeatureMethod','Contrast','Run','Phase','ErrorMessage'});

                    allErrors{end+1,1} = ErrOne; %#ok<AGROW>

                    SummaryOne = table( ...
                        string(cfg.featureMethod), string(cfg.contrastName), cfg.run, string(cfg.phase), ...
                        cfg.Kfeat, cfg.nPC, ...
                        NaN, NaN, NaN, NaN, NaN, ...
                        NaN, string(""), NaN, ...
                        string(cfg.selectedFeatureFile), string('failed'), ...
                        'VariableNames', {'FeatureMethod','Contrast','Run','Phase','Kfeat','nPC', ...
                        'N_total','N_A','N_B','N_subjects','N_features_input', ...
                        'N_selected_or_components','TopFeatureOrComponent','TopScore', ...
                        'SavedFeatureFile','Status'});

                    allSummary{end+1,1} = SummaryOne; %#ok<AGROW>
                end

                %% ---------- Progressive save ----------

                if ~isempty(allSummary)
                    BatchSummary_partial = vertcat(allSummary{:});
                    writetable(BatchSummary_partial, ...
                        fullfile(batchOutDir, 'STEP75A_BATCH_feature_selection_ONLY_summary_partial.csv'));
                end

                if ~isempty(allErrors)
                    BatchErrors_partial = vertcat(allErrors{:});
                    writetable(BatchErrors_partial, ...
                        fullfile(batchOutDir, 'STEP75A_BATCH_feature_selection_ONLY_errors_partial.csv'));
                end
            end
        end
    end
end

%% ===================== FINAL SAVE =====================

if isempty(allSummary)
    BatchSummary = table();
else
    BatchSummary = vertcat(allSummary{:});
end

if isempty(allErrors)
    BatchErrors = table();
else
    BatchErrors = vertcat(allErrors{:});
end

outSummary = fullfile(batchOutDir, 'STEP75A_BATCH_feature_selection_ONLY_summary_all.csv');
outErrors  = fullfile(batchOutDir, 'STEP75A_BATCH_feature_selection_ONLY_errors_all.csv');
outMAT     = fullfile(batchOutDir, 'STEP75A_BATCH_feature_selection_ONLY_results_all.mat');

writetable(BatchSummary, outSummary);
writetable(BatchErrors, outErrors);
save(outMAT, 'BatchSummary', 'BatchErrors', 'baseCfg', ...
    'featureMethodList', 'runList', 'phaseList', 'contrastList', '-v7.3');

fprintf('\n=== FEATURE-SELECTION BATCH FINISHED ===\n');
fprintf('Saved summary:\n%s\n', outSummary);
fprintf('Saved errors:\n%s\n', outErrors);
fprintf('Saved MAT:\n%s\n', outMAT);
fprintf('Finished: %s\n', datestr(now));

diary off;

%% ===================== LOCAL FUNCTION =====================

function cfg = local_refresh_step75_paths(cfg)

    tagFS = sprintf('%s_%s_run%d_%s_K%d', ...
        cfg.contrastName, cfg.featureMethod, cfg.run, cfg.phase, cfg.Kfeat);

    tagCLF = sprintf('%s_%s_%s_run%d_%s_K%d', ...
        cfg.contrastName, cfg.featureMethod, cfg.classifier, cfg.run, cfg.phase, cfg.Kfeat);

    cfg.outRoot = fullfile(cfg.ROOT, '_wm_STEP75_two_stage_FS_then_classifier');
    cfg.fsOutDir = fullfile(cfg.outRoot, 'stage1_feature_selection_ONLY', tagFS);
    cfg.clfOutDir = fullfile(cfg.outRoot, 'stage2_classifier_fixed_features', tagCLF);

    cfg.selectedFeatureFile = fullfile(cfg.fsOutDir, 'STEP75A_selected_features.mat');
end

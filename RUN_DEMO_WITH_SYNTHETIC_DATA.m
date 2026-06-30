%% RUN_DEMO_WITH_SYNTHETIC_DATA.m

clear; clc;

%% ---------------- Demo options ----------------
rngSeed = 2026;
nSubjects = 22;          % keep 22 so scripts with final-N checks pass
nRuns = 3;
phaseList = {'stim','maint','retr'};
condList = {'color','orientation','conjunction'};
nTrialsPerCondition = 6; % small demo; increase if desired

runFastPermutations = false;  % set true to run RUN07/RUN08 after setting nPerm quick values manually
runFigures = false;           % set true only if manuscript figure scripts have been checked for demo compatibility

%% ---------------- Resolve project config ----------------
cfgReal = get_project_config();

projectRoot = cfgReal.projectRoot;
demoOutputRoot = fullfile(projectRoot, 'outputs_demo');
demoProcessedRoot = fullfile(demoOutputRoot, 'new_article_outputs');
demoDatasetDir = fullfile(demoProcessedRoot, '_wm_ml');

if ~exist(demoDatasetDir, 'dir'), mkdir(demoDatasetDir); end

fprintf('\n============================================================\n');
fprintf('RUN_DEMO_WITH_SYNTHETIC_DATA\n');
fprintf('Project root     : %s\n', projectRoot);
fprintf('Demo output root : %s\n', demoOutputRoot);
fprintf('Demo dataset dir : %s\n', demoDatasetDir);
fprintf('============================================================\n');

%% ---------------- Create synthetic feature dataset ----------------
rng(rngSeed);

fprintf('\nCreating synthetic processed dataset...\n');

subjects = strings(nSubjects,1);
for s = 1:nSubjects
    subjects(s) = sprintf('s%02d', s);
end

% Feature inventory. Names are intentionally compatible with STEP100/104
% parsers: they include channel tags such as _ch10.
featureNames = {};
for ch = 1:64
    featureNames{end+1} = sprintf('BP_delta_ch%02d', ch); %#ok<SAGROW>
    featureNames{end+1} = sprintf('BP_theta_ch%02d', ch); %#ok<SAGROW>
    featureNames{end+1} = sprintf('BP_alpha_ch%02d', ch); %#ok<SAGROW>
    featureNames{end+1} = sprintf('BP_beta_ch%02d', ch); %#ok<SAGROW>
    featureNames{end+1} = sprintf('BP_gamma_ch%02d', ch); %#ok<SAGROW>
    featureNames{end+1} = sprintf('RBP_delta_ch%02d', ch); %#ok<SAGROW>
    featureNames{end+1} = sprintf('RBP_theta_ch%02d', ch); %#ok<SAGROW>
    featureNames{end+1} = sprintf('RBP_alpha_ch%02d', ch); %#ok<SAGROW>
    featureNames{end+1} = sprintf('RBP_beta_ch%02d', ch); %#ok<SAGROW>
    featureNames{end+1} = sprintf('RBP_gamma_ch%02d', ch); %#ok<SAGROW>
    featureNames{end+1} = sprintf('RMS_ch%02d', ch); %#ok<SAGROW>
    featureNames{end+1} = sprintf('PeakAbs_ch%02d', ch); %#ok<SAGROW>
    featureNames{end+1} = sprintf('AUCabs_ch%02d', ch); %#ok<SAGROW>
    featureNames{end+1} = sprintf('LineLen_ch%02d', ch); %#ok<SAGROW>
    featureNames{end+1} = sprintf('Mean_ch%02d', ch); %#ok<SAGROW>
    featureNames{end+1} = sprintf('Std_ch%02d', ch); %#ok<SAGROW>
    featureNames{end+1} = sprintf('LZC_ch%02d', ch); %#ok<SAGROW>
end

nFeat = numel(featureNames);
nRows = nSubjects * nRuns * numel(phaseList) * numel(condList) * nTrialsPerCondition;

Subject = strings(nRows,1);
Run = nan(nRows,1);
Phase = strings(nRows,1);
Condition = strings(nRows,1);
TrialNum = nan(nRows,1);

X = randn(nRows, nFeat) * 0.40;

row = 0;
for s = 1:nSubjects
    subjOffset = randn(1,nFeat) * 0.08;

    for r = 1:nRuns
        for p = 1:numel(phaseList)
            phaseName = phaseList{p};

            for c = 1:numel(condList)
                condName = condList{c};

                for tr = 1:nTrialsPerCondition
                    row = row + 1;

                    Subject(row) = subjects(s);
                    Run(row) = r;
                    Phase(row) = string(phaseName);
                    Condition(row) = string(condName);
                    TrialNum(row) = tr;

                    % Subject-level stable offset
                    X(row,:) = X(row,:) + subjOffset;

                    % Mild general condition effects
                    if strcmp(condName,'color')
                        X(row, local_feat_idx(featureNames, {'BP_alpha_ch10','RBP_alpha_ch10','Mean_ch09'})) = ...
                            X(row, local_feat_idx(featureNames, {'BP_alpha_ch10','RBP_alpha_ch10','Mean_ch09'})) + 0.35;
                    elseif strcmp(condName,'orientation')
                        X(row, local_feat_idx(featureNames, {'RBP_delta_ch09','RBP_delta_ch10','RMS_ch26'})) = ...
                            X(row, local_feat_idx(featureNames, {'RBP_delta_ch09','RBP_delta_ch10','RMS_ch26'})) + 0.45;
                    elseif strcmp(condName,'conjunction')
                        X(row, local_feat_idx(featureNames, {'RBP_delta_ch26','LineLen_ch26','LZC_ch26'})) = ...
                            X(row, local_feat_idx(featureNames, {'RBP_delta_ch26','LineLen_ch26','LZC_ch26'})) + 0.50;
                    end

                    % Make the manuscript-highlighted demo task separable:
                    % orientation vs conjunction, run 2, retrieval, channel level.
                    if r == 2 && strcmp(phaseName,'retr')
                        if strcmp(condName,'orientation')
                            X(row, local_feat_idx(featureNames, {'RBP_delta_ch09','RBP_delta_ch10','RBP_alpha_ch10'})) = ...
                                X(row, local_feat_idx(featureNames, {'RBP_delta_ch09','RBP_delta_ch10','RBP_alpha_ch10'})) + 0.85;
                        elseif strcmp(condName,'conjunction')
                            X(row, local_feat_idx(featureNames, {'RBP_delta_ch26','RBP_alpha_ch26','RMS_ch26'})) = ...
                                X(row, local_feat_idx(featureNames, {'RBP_delta_ch26','RBP_alpha_ch26','RMS_ch26'})) + 0.85;
                        end
                    end

                    % Make a secondary ROI demo task mildly separable:
                    % color vs orientation, run 3, maintenance.
                    if r == 3 && strcmp(phaseName,'maint')
                        if strcmp(condName,'color')
                            X(row, local_feat_idx(featureNames, {'BP_alpha_ch05','RBP_alpha_ch06','Mean_ch07'})) = ...
                                X(row, local_feat_idx(featureNames, {'BP_alpha_ch05','RBP_alpha_ch06','Mean_ch07'})) + 0.65;
                        elseif strcmp(condName,'orientation')
                            X(row, local_feat_idx(featureNames, {'BP_theta_ch15','RBP_theta_ch16','Std_ch17'})) = ...
                                X(row, local_feat_idx(featureNames, {'BP_theta_ch15','RBP_theta_ch16','Std_ch17'})) + 0.65;
                        end
                    end
                end
            end
        end
    end
end

Meta = table(Subject, Run, Phase, Condition, TrialNum);
FeatTable = array2table(X, 'VariableNames', featureNames);
DS = [Meta FeatTable];
T = DS; %#ok<NASGU>

datasetMat = fullfile(demoDatasetDir, 'dataset.mat');
datasetCsv = fullfile(demoDatasetDir, 'dataset.csv');

save(datasetMat, 'DS', 'T', '-v7.3');
writetable(DS, datasetCsv);

fprintf('Synthetic dataset saved:\n%s\n', datasetMat);
fprintf('Rows=%d | Features=%d | Subjects=%d\n', height(DS), nFeat, numel(unique(DS.Subject)));

%% ---------------- Create demo config shadow ----------------
% The main RUN scripts call get_project_config(). To avoid overwriting real
% outputs, this demo temporarily shadows get_project_config so that outputRoot
% points to outputs_demo/ while projectRoot still points to the repository.

shadowDir = fullfile(demoOutputRoot, '_demo_config_shadow');
if ~exist(shadowDir, 'dir'), mkdir(shadowDir); end

shadowConfigPath = fullfile(shadowDir, 'get_project_config.m');
local_write_demo_get_project_config(shadowConfigPath, projectRoot, demoOutputRoot);

addpath(shadowDir, '-begin');
clear get_project_config;

cfgDemo = get_project_config();
if ~strcmp(cfgDemo.outputRoot, demoOutputRoot)
    error('Demo get_project_config shadow did not take effect. Expected outputRoot:\n%s\nGot:\n%s', ...
        demoOutputRoot, cfgDemo.outputRoot);
end

fprintf('\nDemo config active. Output root is now:\n%s\n', cfgDemo.outputRoot);

%% ---------------- Run demo workflow ----------------
scripts = { ...
    local_resolve_script(projectRoot, { ...
        fullfile('preprocessing_qc','RUN04_STEP100_subject_median_channel_roi_statistics.m') ...
    }), ...
    local_resolve_script(projectRoot, { ...
        fullfile('preprocessing_qc','RUN05_STEP104_loso_fisher_channel_roi_decoding_K1to22.m') ...
    }), ...
    local_resolve_script(projectRoot, { ...
        fullfile('preprocessing_qc','RUN06_STEP104_collect_loso_decoding_outputs.m') ...
    }) ...
};

stepNames = { ...
    'Subject-level channel/ROI statistics', ...
    'Main LOSO Fisher-score decoding grid', ...
    'Collect STEP104 decoding outputs' ...
};

if runFastPermutations
    scripts{end+1} = local_resolve_script(projectRoot, { ...
        fullfile('preprocessing_qc','RUN07_STEP106_permutation_highlighted_models.m') ...
    });
    stepNames{end+1} = 'Permutation tests for highlighted models';

    scripts{end+1} = local_resolve_script(projectRoot, { ...
        fullfile('preprocessing_qc','RUN08_STEP107B_maxperm_primary_channel_binary_model.m') ...
    });
    stepNames{end+1} = 'Primary-family max-permutation correction';
end

fprintf('\n============================================================\n');
fprintf('Running demo workflow with synthetic data.\n');
fprintf('Permutation scripts enabled: %d\n', runFastPermutations);
fprintf('============================================================\n');

tAll = tic;
for i = 1:numel(scripts)
    fprintf('\n[%d/%d] %s\n', i, numel(scripts), stepNames{i});
    fprintf('Running: %s\n', scripts{i});
    tStep = tic;
    run(scripts{i});
    fprintf('[%d/%d] Completed in %.2f minutes.\n', i, numel(scripts), toc(tStep)/60);
end

fprintf('\n============================================================\n');
fprintf('Demo workflow completed in %.2f minutes.\n', toc(tAll)/60);
fprintf('Demo outputs are under:\n%s\n', demoProcessedRoot);
fprintf('============================================================\n');

if runFigures
    warning('runFigures=true was requested, but figure scripts are project-specific. Run checked manuscript_figure scripts manually if needed.');
end

%% ========================================================================
%% Local helper functions
%% ========================================================================

function idx = local_feat_idx(featureNames, names)
idx = [];
for i = 1:numel(names)
    hit = find(strcmp(featureNames, names{i}), 1, 'first');
    if ~isempty(hit), idx(end+1) = hit; end %#ok<AGROW>
end
end

function local_write_demo_get_project_config(pathOut, projectRoot, demoOutputRoot)
fid = fopen(pathOut, 'w');
if fid < 0
    error('Could not write demo config shadow:\n%s', pathOut);
end

projectRoot = strrep(projectRoot, '''', '''''');
demoOutputRoot = strrep(demoOutputRoot, '''', '''''');

fprintf(fid, 'function cfg = get_project_config()\n');
fprintf(fid, '%% Auto-generated demo config shadow. Do not commit as the main config.\n');
fprintf(fid, 'cfg.projectRoot = ''%s'';\n', projectRoot);
fprintf(fid, 'cfg.codeRoot = cfg.projectRoot;\n');
fprintf(fid, 'cfg.dataRoot = fullfile(cfg.projectRoot, ''data'', ''demo'');\n');
fprintf(fid, 'cfg.outputRoot = ''%s'';\n', demoOutputRoot);
fprintf(fid, 'if ~exist(cfg.outputRoot, ''dir''), mkdir(cfg.outputRoot); end\n');
fprintf(fid, 'activeFolders = { ...\n');
fprintf(fid, '    ''config'', ...\n');
fprintf(fid, '    ''helper_functions'', ...\n');
fprintf(fid, '    ''preprocessing_qc'', ...\n');
fprintf(fid, '    ''feature_extraction'', ...\n');
fprintf(fid, '    ''manuscript_figures'' ...\n');
fprintf(fid, '};\n');
fprintf(fid, 'for i = 1:numel(activeFolders)\n');
fprintf(fid, '    p = fullfile(cfg.projectRoot, activeFolders{i});\n');
fprintf(fid, '    if exist(p, ''dir''), addpath(p); end\n');
fprintf(fid, 'end\n');
fprintf(fid, 'end\n');

fclose(fid);
end

function scriptPath = local_resolve_script(projectRoot, candidates)
scriptPath = '';
for i = 1:numel(candidates)
    p = fullfile(projectRoot, candidates{i});
    if exist(p, 'file') == 2
        scriptPath = p;
        return;
    end
end

fprintf('Could not find any candidate script:');
for i = 1:numel(candidates)
    fprintf('  - %s', fullfile(projectRoot, candidates{i}));
end
error('Required demo workflow script not found.');
end

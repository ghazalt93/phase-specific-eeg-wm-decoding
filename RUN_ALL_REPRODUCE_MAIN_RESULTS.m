%% RUN_ALL_REPRODUCE_MAIN_RESULTS.m

clear; clc;

cfg = get_project_config();

fprintf('\n============================================================\n');
fprintf('Main reproducibility workflow\n');
fprintf('Project root: %s\n', cfg.projectRoot);
fprintf('Output root : %s\n', cfg.outputRoot);
fprintf('============================================================\n');

%% ---------------- Check processed dataset ----------------
processedRoot = fullfile(cfg.outputRoot, 'new_article_outputs');

datasetCandidates = { ...
    fullfile(processedRoot, '_wm_ml', 'dataset.mat'), ...
    fullfile(processedRoot, 'Subjects', '_wm_ml', 'dataset.mat'), ...
    fullfile(processedRoot, '_wm_dataset', 'dataset.mat'), ...
    fullfile(processedRoot, 'dataset.mat') ...
};

datasetFound = false;
for i = 1:numel(datasetCandidates)
    if exist(datasetCandidates{i}, 'file')
        datasetFound = true;
        fprintf('\nFound processed dataset:\n%s\n', datasetCandidates{i});
        break;
    end
end

if ~datasetFound
    fprintf('\nProcessed dataset.mat was not found in the expected locations:\n');
    for i = 1:numel(datasetCandidates)
        fprintf('  - %s\n', datasetCandidates{i});
    end
    error(['Processed dataset not found. Raw EEG/behavioral files are not included in the public repository. ', ...
           'Please place the processed dataset.mat in one of the expected locations or update get_project_config.m.']);
end

%% ---------------- Main workflow scripts ----------------
scripts = { ...
    fullfile(cfg.projectRoot, 'preprocessing_qc', 'RUN04_STEP100_subject_median_channel_roi_statistics.m'), ...
    fullfile(cfg.projectRoot, 'preprocessing_qc', 'RUN05_STEP104_loso_fisher_channel_roi_decoding_K1to22.m'), ...
    fullfile(cfg.projectRoot, 'preprocessing_qc', 'RUN06_STEP104_collect_loso_decoding_outputs.m'), ...
    fullfile(cfg.projectRoot, 'preprocessing_qc', 'RUN07_STEP106_permutation_highlighted_models.m'), ...
    fullfile(cfg.projectRoot, 'preprocessing_qc', 'RUN08_STEP107B_maxperm_primary_channel_binary_model.m') ...
};

stepNames = { ...
    'Subject-level channel/ROI statistics', ...
    'Main LOSO Fisher-score decoding grid', ...
    'Collect STEP104 decoding outputs', ...
    'Permutation tests for highlighted models', ...
    'Primary-family max-permutation correction' ...
};

% Check that all required scripts exist before starting.
for i = 1:numel(scripts)
    if exist(scripts{i}, 'file') ~= 2
        error('Required workflow script not found:\n%s', scripts{i});
    end
end

%% ---------------- Run workflow ----------------
tStartAll = tic;

for i = 1:numel(scripts)
    fprintf('\n[%d/%d] %s\n', i, numel(scripts), stepNames{i});
    fprintf('Running: %s\n', scripts{i});

    tStart = tic;
    run(scripts{i});
    fprintf('[%d/%d] Completed in %.2f minutes.\n', i, numel(scripts), toc(tStart)/60);
end

fprintf('\n============================================================\n');
fprintf('Main reproducibility workflow completed in %.2f minutes.\n', toc(tStartAll)/60);
fprintf('Outputs are saved under:\n%s\n', processedRoot);
fprintf('============================================================\n');

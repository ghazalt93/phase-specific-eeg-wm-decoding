%% STEP75A_run_feature_selection_ONLY.m

clear; clc;

cfg = WM75_config();

if ~exist(cfg.fsOutDir, 'dir')
    mkdir(cfg.fsOutDir);
end

diary(fullfile(cfg.fsOutDir, 'STEP75A_feature_selection_ONLY_log.txt'));

fprintf('\n=== STEP75A Feature Selection ONLY ===\n');
fprintf('Started: %s\n', datestr(now));
fprintf('Target: %s | run=%d | phase=%s\n', cfg.contrastName, cfg.run, cfg.phase);
fprintf('Feature method: %s | K=%d | nPC=%d\n', cfg.featureMethod, cfg.Kfeat, cfg.nPC);

[X, y, subj, featureNames, DataInfo] = STEP75_load_binary_data(cfg);

fprintf('Data loaded: N=%d | N_A=%d | N_B=%d | Subjects=%d | Features=%d\n', ...
    numel(y), sum(y=='A'), sum(y=='B'), numel(unique(subj)), size(X,2));

% Preprocess using the full selected subset.
% This is intentional for this diagnostic two-stage pipeline.
[Xz, prep] = STEP75_preprocess_all_for_FS(X);

[FSmodel, Ranking] = STEP75_feature_selection_full(Xz, y, featureNames, cfg);

outRankCSV = fullfile(cfg.fsOutDir, 'STEP75A_feature_ranking.csv');
outTopCSV  = fullfile(cfg.fsOutDir, 'STEP75A_top_selected_features.csv');
outMAT     = cfg.selectedFeatureFile;

writetable(Ranking, outRankCSV);

if isfield(FSmodel,'selectedIdx') && ~isempty(FSmodel.selectedIdx)
    TopSelected = Ranking(1:min(height(Ranking), numel(FSmodel.selectedIdx)), :);
else
    TopSelected = Ranking;
end

writetable(TopSelected, outTopCSV);

save(outMAT, 'FSmodel', 'Ranking', 'prep', 'DataInfo', 'cfg', '-v7.3');

fprintf('\nSaved Stage 1 outputs:\n');
fprintf('  Ranking : %s\n', outRankCSV);
fprintf('  Top     : %s\n', outTopCSV);
fprintf('  MAT     : %s\n', outMAT);

fprintf('Finished: %s\n', datestr(now));
diary off;

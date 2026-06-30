%% STEP66_normality_KS_features.m

cfg = get_project_config();

clear; clc;

%% ===================== USER CONFIG =====================

cfg.ROOT = fullfile(cfg.outputRoot, 'new_article_outputs');

% If your dataset.mat is somewhere else, set it manually here:
cfg.datasetPath = '';

% Common candidate paths if cfg.datasetPath is empty
candidatePaths = {
    fullfile(cfg.ROOT, '_wm_ml', 'dataset.mat')
    fullfile(cfg.ROOT, 'Subjects', '_wm_ml', 'dataset.mat')
    fullfile(cfg.ROOT, '_wm_dataset', 'dataset.mat')
    fullfile(cfg.ROOT, 'dataset.mat')
};

cfg.outDir = fullfile(cfg.ROOT, '_normality_QC');

cfg.alpha = 0.05;
cfg.minN = 20;

cfg.groupMode = "run_phase";

cfg.doFDR = true;

cfg.makeSummaryFigure = true;

%% ===================== LOAD DATASET =====================

if ~exist(cfg.outDir, 'dir')
    mkdir(cfg.outDir);
end

if isempty(cfg.datasetPath)
    cfg.datasetPath = '';
    for i = 1:numel(candidatePaths)
        if exist(candidatePaths{i}, 'file')
            cfg.datasetPath = candidatePaths{i};
            break;
        end
    end
end

if isempty(cfg.datasetPath) || ~exist(cfg.datasetPath, 'file')
    error('dataset.mat not found. Please set cfg.datasetPath manually.');
end

fprintf('Loading dataset:\n%s\n', cfg.datasetPath);
S = load(cfg.datasetPath);

T = [];
if isfield(S, 'DS')
    if istable(S.DS)
        T = S.DS;
    elseif isstruct(S.DS)
        fn = fieldnames(S.DS);
        for i = 1:numel(fn)
            if istable(S.DS.(fn{i}))
                T = S.DS.(fn{i});
                break;
            end
        end
    end
elseif isfield(S, 'T') && istable(S.T)
    T = S.T;
end

if isempty(T) || ~istable(T)
    error('Could not find a table variable named DS or T inside dataset.mat.');
end

fprintf('Loaded table: rows = %d, variables = %d\n', height(T), width(T));

%% ===================== FIND FEATURE COLUMNS =====================

allNames = T.Properties.VariableNames;
isNum = false(1, numel(allNames));

for j = 1:numel(allNames)
    x = T.(allNames{j});
    isNum(j) = isnumeric(x) || islogical(x);
end

numNames = allNames(isNum);

metaNames = { ...
    'Subject'; 'subject'; 'SubjectID'; 'subjectID'; 'subj'; 'subjID'; 'Subj'; 'SubjID'; ...
    'Run'; 'run'; 'runNum'; 'RunNum'; 'session'; 'Session'; ...
    'Phase'; 'phase'; ...
    'Condition'; 'condition'; 'yCondition'; 'Label'; 'label'; ...
    'Correct'; 'correct'; 'yCorrect'; ...
    'TrialNum'; 'trialNum'; 'TrialIndex'; 'trialIndex'; 'Trial'; 'trial'; ...
    'PatternID'; 'patternID'; ...
    'StartRow'; 'EndRow'; 'Second10Row'; 'RetrRow'; ...
    'Fold'; 'fold'; 'CVFold' ...
};
isMeta = false(size(numNames));
for j = 1:numel(numNames)
    isMeta(j) = any(strcmpi(numNames{j}, metaNames));
end

featureNames = numNames(~isMeta);

fprintf('Numeric variables found: %d\n', numel(numNames));
fprintf('Feature candidates after metadata exclusion: %d\n', numel(featureNames));

if isempty(featureNames)
    error('No numeric feature columns found after metadata exclusion.');
end

%% ===================== BUILD GROUPS =====================

[groupLabels, groupInfo] = make_groups(T, cfg.groupMode);

uGroups = unique(groupLabels);
uGroups = uGroups(~ismissing(uGroups));

fprintf('Number of groups: %d\n', numel(uGroups));
disp(table(uGroups, 'VariableNames', {'Group'}));

%% ===================== RUN NORMALITY TESTS =====================

results = struct([]);
rr = 0;

for g = 1:numel(uGroups)

    thisGroup = uGroups(g);
    idxG = groupLabels == thisGroup;

    fprintf('\nGroup %d/%d: %s | rows = %d\n', g, numel(uGroups), thisGroup, sum(idxG));

    for f = 1:numel(featureNames)

        fname = featureNames{f};
        x = T.(fname);
        x = x(idxG);

        x = double(x(:));
        x = x(isfinite(x));

        n = numel(x);
        nMissing = sum(idxG) - n;

        rr = rr + 1;

        results(rr).Group = char(thisGroup);
        results(rr).Feature = fname;
        results(rr).N = n;
        results(rr).N_missing_or_nonfinite = nMissing;

        if n < cfg.minN
            results(rr).Status = 'too_few_samples';
            results(rr).Mean = NaN;
            results(rr).SD = NaN;
            results(rr).Skewness = NaN;
            results(rr).Kurtosis = NaN;
            results(rr).KS_h = NaN;
            results(rr).KS_p = NaN;
            results(rr).KS_stat = NaN;
            results(rr).Lillie_h = NaN;
            results(rr).Lillie_p = NaN;
            results(rr).Lillie_stat = NaN;
            continue;
        end

        mu = mean(x, 'omitnan');
        sd = std(x, 0, 'omitnan');

        results(rr).Mean = mu;
        results(rr).SD = sd;
        results(rr).Skewness = skewness(x, 0);
        results(rr).Kurtosis = kurtosis(x, 0);

        if sd <= eps || ~isfinite(sd)
            results(rr).Status = 'constant_or_near_constant';
            results(rr).KS_h = NaN;
            results(rr).KS_p = NaN;
            results(rr).KS_stat = NaN;
            results(rr).Lillie_h = NaN;
            results(rr).Lillie_p = NaN;
            results(rr).Lillie_stat = NaN;
            continue;
        end

        z = (x - mu) ./ sd;

        % KS test against standard normal
        try
            [hKS, pKS, ksstat] = kstest(z, 'Alpha', cfg.alpha);
        catch
            hKS = NaN;
            pKS = NaN;
            ksstat = NaN;
        end

        % Lilliefors test: better when mean/std are estimated
        try
            [hLil, pLil, lilstat] = lillietest(x, 'Alpha', cfg.alpha);
        catch
            hLil = NaN;
            pLil = NaN;
            lilstat = NaN;
        end

        results(rr).Status = 'ok';
        results(rr).KS_h = hKS;
        results(rr).KS_p = pKS;
        results(rr).KS_stat = ksstat;
        results(rr).Lillie_h = hLil;
        results(rr).Lillie_p = pLil;
        results(rr).Lillie_stat = lilstat;

    end
end

R = struct2table(results);

%% ===================== FDR CORRECTION =====================

R.KS_q = NaN(height(R), 1);
R.Lillie_q = NaN(height(R), 1);

if cfg.doFDR
    for g = 1:numel(uGroups)
        idx = strcmp(string(R.Group), uGroups(g));

        R.KS_q(idx) = bh_fdr(R.KS_p(idx));
        R.Lillie_q(idx) = bh_fdr(R.Lillie_p(idx));
    end
end

R.KS_reject_FDR = R.KS_q < cfg.alpha;
R.Lillie_reject_FDR = R.Lillie_q < cfg.alpha;

%% ===================== SUMMARY TABLE =====================

summaryRows = struct([]);

for g = 1:numel(uGroups)
    idx = strcmp(string(R.Group), uGroups(g)) & strcmp(R.Status, 'ok');

    nTested = sum(idx);

    nKS_raw = sum(R.KS_p(idx) < cfg.alpha, 'omitnan');
    nLil_raw = sum(R.Lillie_p(idx) < cfg.alpha, 'omitnan');

    nKS_fdr = sum(R.KS_reject_FDR(idx), 'omitnan');
    nLil_fdr = sum(R.Lillie_reject_FDR(idx), 'omitnan');

    summaryRows(g).Group = char(uGroups(g));
    summaryRows(g).N_features_tested = nTested;

    summaryRows(g).KS_reject_raw = nKS_raw;
    summaryRows(g).KS_reject_raw_percent = 100 * nKS_raw / max(nTested, 1);

    summaryRows(g).Lillie_reject_raw = nLil_raw;
    summaryRows(g).Lillie_reject_raw_percent = 100 * nLil_raw / max(nTested, 1);

    summaryRows(g).KS_reject_FDR = nKS_fdr;
    summaryRows(g).KS_reject_FDR_percent = 100 * nKS_fdr / max(nTested, 1);

    summaryRows(g).Lillie_reject_FDR = nLil_fdr;
    summaryRows(g).Lillie_reject_FDR_percent = 100 * nLil_fdr / max(nTested, 1);
end

Summary = struct2table(summaryRows);

%% ===================== SAVE OUTPUTS =====================

outCSV = fullfile(cfg.outDir, 'STEP66_feature_normality_KS_Lillie_results.csv');
outSummaryCSV = fullfile(cfg.outDir, 'STEP66_feature_normality_KS_Lillie_summary.csv');
outMAT = fullfile(cfg.outDir, 'STEP66_feature_normality_KS_Lillie_results.mat');

writetable(R, outCSV);
writetable(Summary, outSummaryCSV);
save(outMAT, 'R', 'Summary', 'cfg', 'featureNames', 'groupInfo', '-v7.3');

fprintf('\nSaved detailed results:\n%s\n', outCSV);
fprintf('Saved summary:\n%s\n', outSummaryCSV);
fprintf('Saved MAT:\n%s\n', outMAT);

%% ===================== OPTIONAL SUMMARY FIGURE =====================

if cfg.makeSummaryFigure
    try
        fig = figure('Color', 'w', 'Position', [100 100 1100 450]);

        y = Summary.Lillie_reject_FDR_percent;
        bar(y);
        xticks(1:height(Summary));
        xticklabels(Summary.Group);
        xtickangle(45);
        ylabel('% rejected after FDR');
        title('Feature normality QC: Lilliefors rejection rate after FDR');
        grid on;

        figPath = fullfile(cfg.outDir, 'STEP66_Lillie_FDR_rejection_percent.png');
        saveas(fig, figPath);
        close(fig);

        fprintf('Saved summary figure:\n%s\n', figPath);
    catch ME
        warning('Could not save summary figure: %s', ME.message);
    end
end

%% ===================== LOCAL FUNCTIONS =====================

function [groupLabels, groupInfo] = make_groups(T, groupMode)

    n = height(T);
    groupLabels = strings(n, 1);
    groupInfo = struct();
    groupMode = string(groupMode);

    names = T.Properties.VariableNames;

    phaseName = find_var(names, {'phase','Phase'});
    runName = find_var(names, {'runNum','RunNum','run','Run'});
    condName = find_var(names, {'Condition','condition','yCondition','Label','label'});

    switch groupMode

        case "all"
            groupLabels(:) = "all";

        case "phase"
            if isempty(phaseName)
                warning('Phase column not found. Using all rows as one group.');
                groupLabels(:) = "all";
            else
                groupLabels = "phase_" + clean_str(T.(phaseName));
            end

        case "run_phase"
            if isempty(phaseName) || isempty(runName)
                warning('Phase or run column not found. Using all rows as one group.');
                groupLabels(:) = "all";
            else
                groupLabels = "run" + clean_str(T.(runName)) + "_" + clean_str(T.(phaseName));
            end

        case "run_phase_condition"
            if isempty(phaseName) || isempty(runName) || isempty(condName)
                warning('Phase/run/condition column not found. Falling back to run_phase.');
                if isempty(phaseName) || isempty(runName)
                    groupLabels(:) = "all";
                else
                    groupLabels = "run" + clean_str(T.(runName)) + "_" + clean_str(T.(phaseName));
                end
            else
                groupLabels = "run" + clean_str(T.(runName)) + "_" + ...
                              clean_str(T.(phaseName)) + "_" + ...
                              clean_str(T.(condName));
            end

        otherwise
            warning('Unknown groupMode. Using all rows as one group.');
            groupLabels(:) = "all";
    end

    groupInfo.phaseName = phaseName;
    groupInfo.runName = runName;
    groupInfo.condName = condName;
end

function vname = find_var(names, candidates)
    vname = '';
    for i = 1:numel(candidates)
        idx = strcmp(names, candidates{i});
        if any(idx)
            vname = names{find(idx, 1)};
            return;
        end
    end
end

function s = clean_str(x)
    if isnumeric(x) || islogical(x)
        s = string(x);
    elseif iscell(x)
        s = string(x);
    elseif iscategorical(x)
        s = string(x);
    elseif isstring(x)
        s = x;
    elseif ischar(x)
        s = string(cellstr(x));
    else
        s = string(x);
    end

    s = lower(strtrim(s));
    s = regexprep(s, '\s+', '');
    s = regexprep(s, '[^\w]', '');
end

function q = bh_fdr(p)
    p = double(p(:));
    q = NaN(size(p));

    valid = isfinite(p);
    pv = p(valid);

    m = numel(pv);
    if m == 0
        return;
    end

    [ps, order] = sort(pv, 'ascend');
    ranks = (1:m)';

    qs = ps .* m ./ ranks;
    qs = flipud(cummin(flipud(qs)));
    qs(qs > 1) = 1;

    qv = NaN(size(pv));
    qv(order) = qs;

    q(valid) = qv;
end
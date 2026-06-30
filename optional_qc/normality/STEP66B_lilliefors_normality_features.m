%% STEP66B_lilliefors_normality_features.m

cfg = get_project_config();

clear; clc;

%% ===================== USER CONFIG =====================

cfg.ROOT = fullfile(cfg.outputRoot, 'new_article_outputs');

cfg.datasetPath = '';

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

%% ===================== RUN LILLIEFORS TEST =====================

results = struct([]);
rr = 0;

for g = 1:numel(uGroups)

    thisGroup = uGroups(g);
    idxG = groupLabels == thisGroup;

    fprintf('\nGroup %d/%d: %s | rows = %d\n', ...
        g, numel(uGroups), thisGroup, sum(idxG));

    for f = 1:numel(featureNames)

        fname = featureNames{f};
        x = T.(fname);
        x = double(x(idxG));
        x = x(:);
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
            results(rr).Lillie_h = NaN;
            results(rr).Lillie_p = NaN;
            results(rr).Lillie_stat = NaN;
            continue;
        end

        try
            [h, p, stat] = lillietest(x, 'Alpha', cfg.alpha);
        catch ME
            warning('Lilliefors failed for %s | %s | %s', ...
                char(thisGroup), fname, ME.message);

            h = NaN;
            p = NaN;
            stat = NaN;
        end

        results(rr).Status = 'ok';
        results(rr).Lillie_h = h;
        results(rr).Lillie_p = p;
        results(rr).Lillie_stat = stat;
    end
end

R = struct2table(results);

%% ===================== FDR CORRECTION =====================

R.Lillie_q = NaN(height(R), 1);

if cfg.doFDR
    for g = 1:numel(uGroups)
        idx = strcmp(string(R.Group), uGroups(g));
        R.Lillie_q(idx) = bh_fdr(R.Lillie_p(idx));
    end
end

R.Lillie_reject_raw = R.Lillie_p < cfg.alpha;
R.Lillie_reject_FDR = R.Lillie_q < cfg.alpha;

%% ===================== SUMMARY TABLE =====================

summaryRows = struct([]);

for g = 1:numel(uGroups)

    idx = strcmp(string(R.Group), uGroups(g)) & strcmp(R.Status, 'ok');

    nTested = sum(idx);

    nRaw = sum(R.Lillie_reject_raw(idx), 'omitnan');
    nFDR = sum(R.Lillie_reject_FDR(idx), 'omitnan');

    summaryRows(g).Group = char(uGroups(g));
    summaryRows(g).N_features_tested = nTested;

    summaryRows(g).Lillie_reject_raw = nRaw;
    summaryRows(g).Lillie_reject_raw_percent = 100 * nRaw / max(nTested, 1);

    summaryRows(g).Lillie_reject_FDR = nFDR;
    summaryRows(g).Lillie_reject_FDR_percent = 100 * nFDR / max(nTested, 1);

    summaryRows(g).Median_p = median(R.Lillie_p(idx), 'omitnan');
    summaryRows(g).Median_q = median(R.Lillie_q(idx), 'omitnan');

end

Summary = struct2table(summaryRows);

%% ===================== SAVE OUTPUTS =====================

outCSV = fullfile(cfg.outDir, 'STEP66B_lilliefors_feature_normality_results.csv');
outSummaryCSV = fullfile(cfg.outDir, 'STEP66B_lilliefors_feature_normality_summary.csv');
outMAT = fullfile(cfg.outDir, 'STEP66B_lilliefors_feature_normality_results.mat');

writetable(R, outCSV);
writetable(Summary, outSummaryCSV);

save(outMAT, 'R', 'Summary', 'cfg', 'featureNames', 'groupInfo', '-v7.3');

fprintf('\nSaved detailed results:\n%s\n', outCSV);
fprintf('Saved summary:\n%s\n', outSummaryCSV);
fprintf('Saved MAT:\n%s\n', outMAT);

%% ===================== SUMMARY FIGURE =====================

if cfg.makeSummaryFigure
    try
        fig = figure('Color', 'w', 'Position', [100 100 1100 450]);

        bar(Summary.Lillie_reject_FDR_percent);
        xticks(1:height(Summary));
        xticklabels(Summary.Group);
        xtickangle(45);

        ylabel('% rejected after FDR');
        title('Lilliefors normality test: rejection rate after FDR');
        grid on;

        figPath = fullfile(cfg.outDir, 'STEP66B_lilliefors_FDR_rejection_percent.png');
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
    groupMode = string(groupMode);

    names = T.Properties.VariableNames;

    phaseName = find_var(names, {'phase','Phase'});
    runName   = find_var(names, {'runNum','RunNum','run','Run'});
    condName  = find_var(names, {'Condition','condition','yCondition','Label','label'});

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
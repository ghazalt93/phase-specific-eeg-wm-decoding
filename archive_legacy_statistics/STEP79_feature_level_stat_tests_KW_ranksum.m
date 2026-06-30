%% STEP79_feature_level_stat_tests_KW_ranksum.m

cfg = get_project_config();

clear; clc;

%% ===================== CONFIG =====================

cfg = struct();

cfg.ROOT = cfg.outputRoot;
cfg.datasetPath = '';

cfg.candidatePaths = {
    fullfile(cfg.ROOT, '_wm_ml', 'dataset.mat')
    fullfile(cfg.ROOT, 'Subjects', '_wm_ml', 'dataset.mat')
    fullfile(cfg.ROOT, '_wm_dataset', 'dataset.mat')
    fullfile(cfg.ROOT, 'dataset.mat')
};

cfg.outDir = fullfile(cfg.ROOT, '_wm_STEP79_feature_level_stat_tests');

cfg.runList = [1 2 3];
cfg.phaseList = {'stim','maint','retr'};

cfg.alphaFDR = 0.05;

cfg.minN_threeClass_perClass = 10;
cfg.minN_binary_perGroup = 10;

cfg.labelMap.color       = {'1','color'};
cfg.labelMap.orientation = {'2','orientation'};
cfg.labelMap.conjunction = {'3','conjunction'};

cfg.binaryContrasts = struct([]);

cfg.binaryContrasts(1).name   = 'orientation_vs_conjunction';
cfg.binaryContrasts(1).classA = cfg.labelMap.orientation;
cfg.binaryContrasts(1).classB = cfg.labelMap.conjunction;

cfg.binaryContrasts(2).name   = 'color_vs_orientation';
cfg.binaryContrasts(2).classA = cfg.labelMap.color;
cfg.binaryContrasts(2).classB = cfg.labelMap.orientation;

cfg.binaryContrasts(3).name   = 'color_vs_conjunction';
cfg.binaryContrasts(3).classA = cfg.labelMap.color;
cfg.binaryContrasts(3).classB = cfg.labelMap.conjunction;

cfg.binaryContrasts(4).name   = 'conjunction_vs_notConjunction';
cfg.binaryContrasts(4).classA = cfg.labelMap.conjunction;
cfg.binaryContrasts(4).classB = [cfg.labelMap.color, cfg.labelMap.orientation];

cfg.binaryContrasts(5).name   = 'orientation_vs_notOrientation';
cfg.binaryContrasts(5).classA = cfg.labelMap.orientation;
cfg.binaryContrasts(5).classB = [cfg.labelMap.color, cfg.labelMap.conjunction];

cfg.binaryContrasts(6).name   = 'color_vs_notColor';
cfg.binaryContrasts(6).classA = cfg.labelMap.color;
cfg.binaryContrasts(6).classB = [cfg.labelMap.orientation, cfg.labelMap.conjunction];

%% ===================== PREPARE =====================

if ~exist(cfg.outDir, 'dir')
    mkdir(cfg.outDir);
end

diary(fullfile(cfg.outDir, 'STEP79_feature_level_stat_tests_log.txt'));

fprintf('\n=== STEP79 feature-level statistical tests ===\n');
fprintf('Started: %s\n', datestr(now));
fprintf('ROOT: %s\n', cfg.ROOT);
fprintf('FDR alpha: %.3f\n', cfg.alphaFDR);

%% ===================== LOAD DATASET =====================

[T, col, featureNames] = local_load_dataset_and_columns(cfg);

condVals  = local_clean_str(T.(col.cond));
runVals   = double(T.(col.run));
phaseVals = local_clean_str(T.(col.phase));

fprintf('\nDataset loaded: rows=%d | features=%d\n', height(T), numel(featureNames));

%% ===================== MAIN TESTS =====================

allKW = {};
allRS = {};

for r = 1:numel(cfg.runList)

    runNum = cfg.runList(r);

    for p = 1:numel(cfg.phaseList)

        phaseName = cfg.phaseList{p};

        fprintf('\n\n====================================================\n');
        fprintf('Run %d | Phase %s\n', runNum, phaseName);
        fprintf('====================================================\n');

        idxBase = runVals == runNum & phaseVals == local_clean_str({phaseName});

        if sum(idxBase) < 20
            warning('Too few rows for run=%d phase=%s. Skipping.', runNum, phaseName);
            continue;
        end

        %% ---------- THREE-CLASS KRUSKAL-WALLIS ----------

        idxColor = idxBase & ismember(condVals, local_clean_str(cfg.labelMap.color));
        idxOri   = idxBase & ismember(condVals, local_clean_str(cfg.labelMap.orientation));
        idxConj  = idxBase & ismember(condVals, local_clean_str(cfg.labelMap.conjunction));

        idxKW = idxColor | idxOri | idxConj;

        y3 = strings(sum(idxKW),1);
        y3(idxColor(idxKW)) = "color";
        y3(idxOri(idxKW))   = "orientation";
        y3(idxConj(idxKW))  = "conjunction";
        y3 = categorical(y3);

        X3 = table2array(T(idxKW, featureNames));

        nColor = sum(y3 == 'color');
        nOri   = sum(y3 == 'orientation');
        nConj  = sum(y3 == 'conjunction');

        fprintf('[KW] N color=%d | orientation=%d | conjunction=%d\n', ...
            nColor, nOri, nConj);

        if min([nColor, nOri, nConj]) >= cfg.minN_threeClass_perClass

            KW = local_run_kruskalwallis_features(X3, y3, featureNames);

            KW.Run = repmat(runNum, height(KW), 1);
            KW.Phase = repmat(string(phaseName), height(KW), 1);
            KW.Test = repmat("KruskalWallis_3class", height(KW), 1);
            KW.N_color = repmat(nColor, height(KW), 1);
            KW.N_orientation = repmat(nOri, height(KW), 1);
            KW.N_conjunction = repmat(nConj, height(KW), 1);

            KW = movevars(KW, {'Test','Run','Phase','N_color','N_orientation','N_conjunction'}, 'Before', 1);

            allKW{end+1,1} = KW; %#ok<AGROW>

            fprintf('[KW] Done. Min p = %.3g | N FDR sig = %d\n', ...
                min(KW.p_KW, [], 'omitnan'), sum(KW.q_FDR < cfg.alphaFDR, 'omitnan'));
        else
            warning('[KW] Too few samples for run=%d phase=%s. Skipping KW.', runNum, phaseName);
        end

        %% ---------- BINARY RANK-SUM TESTS ----------

        for c = 1:numel(cfg.binaryContrasts)

            C = cfg.binaryContrasts(c);

            idxA = idxBase & ismember(condVals, local_clean_str(C.classA));
            idxB = idxBase & ismember(condVals, local_clean_str(C.classB));

            idxUse = idxA | idxB;

            nA = sum(idxA);
            nB = sum(idxB);

            fprintf('[RS] %s | N_A=%d | N_B=%d\n', C.name, nA, nB);

            if nA < cfg.minN_binary_perGroup || nB < cfg.minN_binary_perGroup
                warning('[RS] Too few samples for %s run=%d phase=%s. Skipping.', ...
                    C.name, runNum, phaseName);
                continue;
            end

            y2 = strings(sum(idxUse),1);
            y2(idxA(idxUse)) = "A";
            y2(idxB(idxUse)) = "B";
            y2 = categorical(y2);

            X2 = table2array(T(idxUse, featureNames));

            RS = local_run_ranksum_features(X2, y2, featureNames);

            RS.Run = repmat(runNum, height(RS), 1);
            RS.Phase = repmat(string(phaseName), height(RS), 1);
            RS.Contrast = repmat(string(C.name), height(RS), 1);
            RS.Test = repmat("Ranksum_binary", height(RS), 1);
            RS.N_A = repmat(nA, height(RS), 1);
            RS.N_B = repmat(nB, height(RS), 1);

            RS = movevars(RS, {'Test','Contrast','Run','Phase','N_A','N_B'}, 'Before', 1);

            allRS{end+1,1} = RS; %#ok<AGROW>

            fprintf('[RS] Done. Min p = %.3g | N FDR sig = %d\n', ...
                min(RS.p_ranksum, [], 'omitnan'), sum(RS.q_FDR < cfg.alphaFDR, 'omitnan'));
        end
    end
end

%% ===================== CONCATENATE =====================

if isempty(allKW)
    KW_All = table();
else
    KW_All = vertcat(allKW{:});
end

if isempty(allRS)
    RS_All = table();
else
    RS_All = vertcat(allRS{:});
end

%% ===================== SUMMARIES =====================

KW_Summary = local_summary_kw(KW_All, cfg.alphaFDR);
RS_Summary = local_summary_rs(RS_All, cfg.alphaFDR);

KW_FamilySummary = local_family_summary(KW_All, 'p_KW', cfg.alphaFDR);
RS_FamilySummary = local_family_summary(RS_All, 'p_ranksum', cfg.alphaFDR);

%% ===================== SAVE =====================

outKW = fullfile(cfg.outDir, 'STEP79_KruskalWallis_3class_feature_pvalues.csv');
outRS = fullfile(cfg.outDir, 'STEP79_Ranksum_binary_feature_pvalues.csv');

outKWSum = fullfile(cfg.outDir, 'STEP79_KruskalWallis_3class_summary.csv');
outRSSum = fullfile(cfg.outDir, 'STEP79_Ranksum_binary_summary.csv');

outKWFamily = fullfile(cfg.outDir, 'STEP79_KruskalWallis_3class_family_summary.csv');
outRSFamily = fullfile(cfg.outDir, 'STEP79_Ranksum_binary_family_summary.csv');

outMAT = fullfile(cfg.outDir, 'STEP79_feature_level_stat_tests_results.mat');

writetable(KW_All, outKW);
writetable(RS_All, outRS);

writetable(KW_Summary, outKWSum);
writetable(RS_Summary, outRSSum);

writetable(KW_FamilySummary, outKWFamily);
writetable(RS_FamilySummary, outRSFamily);

save(outMAT, 'KW_All', 'RS_All', 'KW_Summary', 'RS_Summary', ...
    'KW_FamilySummary', 'RS_FamilySummary', 'cfg', '-v7.3');

fprintf('\nSaved outputs:\n');
fprintf('  KW p-values        : %s\n', outKW);
fprintf('  Ranksum p-values   : %s\n', outRS);
fprintf('  KW summary         : %s\n', outKWSum);
fprintf('  Ranksum summary    : %s\n', outRSSum);
fprintf('  KW family summary  : %s\n', outKWFamily);
fprintf('  RS family summary  : %s\n', outRSFamily);
fprintf('  MAT                : %s\n', outMAT);

fprintf('\nFinished: %s\n', datestr(now));
diary off;

%% ===================== LOCAL FUNCTIONS =====================

function [T, col, featureNames] = local_load_dataset_and_columns(cfg)

    datasetPath = cfg.datasetPath;

    if isempty(datasetPath)
        for i = 1:numel(cfg.candidatePaths)
            if exist(cfg.candidatePaths{i}, 'file')
                datasetPath = cfg.candidatePaths{i};
                break;
            end
        end
    end

    if isempty(datasetPath) || ~exist(datasetPath, 'file')
        error('dataset.mat not found. Set cfg.datasetPath manually.');
    end

    fprintf('Loading dataset:\n%s\n', datasetPath);
    S = load(datasetPath);

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
        error('Could not find table DS or T inside dataset.mat.');
    end

    names = T.Properties.VariableNames;

    col.subj  = local_find_var(names, {'Subject','subject','SubjectID','subjectID','subj','subjID','Subj','SubjID'});
    col.run   = local_find_var(names, {'runNum','RunNum','run','Run'});
    col.phase = local_find_var(names, {'phase','Phase'});
    col.cond  = local_find_var(names, {'Condition','condition','yCondition','Label','label'});

    if isempty(col.run),   error('Run column not found.'); end
    if isempty(col.phase), error('Phase column not found.'); end
    if isempty(col.cond),  error('Condition column not found.'); end

    featureNames = local_detect_feature_columns(T);
end

function KW = local_run_kruskalwallis_features(X, y, featureNames)

    featureNames = string(featureNames(:));
    nFeat = size(X,2);

    pval = NaN(nFeat,1);
    chi2stat = NaN(nFeat,1);
    df = NaN(nFeat,1);
    epsilon2 = NaN(nFeat,1);

    med_color = NaN(nFeat,1);
    med_orientation = NaN(nFeat,1);
    med_conjunction = NaN(nFeat,1);

    meanRank_color = NaN(nFeat,1);
    meanRank_orientation = NaN(nFeat,1);
    meanRank_conjunction = NaN(nFeat,1);

    for j = 1:nFeat

        x = double(X(:,j));
        ok = isfinite(x) & ~isundefined(y);

        x = x(ok);
        yy = y(ok);

        if numel(categories(removecats(yy))) < 3
            continue;
        end

        if std(x,0,'omitnan') <= eps
            continue;
        end

        try
            [p, tbl] = kruskalwallis(x, yy, 'off');

            pval(j) = p;

            chi2 = tbl{2,5};
            dfi  = tbl{2,3};

            chi2stat(j) = chi2;
            df(j) = dfi;

            N = numel(x);
            k = 3;

            epsilon2(j) = max(0, (chi2 - k + 1) / max(N - k, 1));

            med_color(j)       = median(x(yy == 'color'), 'omitnan');
            med_orientation(j) = median(x(yy == 'orientation'), 'omitnan');
            med_conjunction(j) = median(x(yy == 'conjunction'), 'omitnan');

            ranks = tiedrank(x);
            meanRank_color(j)       = mean(ranks(yy == 'color'), 'omitnan');
            meanRank_orientation(j) = mean(ranks(yy == 'orientation'), 'omitnan');
            meanRank_conjunction(j) = mean(ranks(yy == 'conjunction'), 'omitnan');

        catch
            pval(j) = NaN;
        end
    end

    q = local_bh_fdr(pval);

    KW = table( ...
        featureNames, local_feature_family(featureNames), ...
        pval, q, chi2stat, df, epsilon2, ...
        med_color, med_orientation, med_conjunction, ...
        meanRank_color, meanRank_orientation, meanRank_conjunction, ...
        'VariableNames', {'Feature','Family','p_KW','q_FDR','KW_chi2','KW_df','epsilon2_KW', ...
        'Median_color','Median_orientation','Median_conjunction', ...
        'MeanRank_color','MeanRank_orientation','MeanRank_conjunction'});
end

function RS = local_run_ranksum_features(X, y, featureNames)

    featureNames = string(featureNames(:));

    nFeat = size(X,2);

    pval = NaN(nFeat,1);
    zval = NaN(nFeat,1);
    U_A = NaN(nFeat,1);
    aucA = NaN(nFeat,1);
    rankBiserial = NaN(nFeat,1);

    med_A = NaN(nFeat,1);
    med_B = NaN(nFeat,1);

    meanRank_A = NaN(nFeat,1);
    meanRank_B = NaN(nFeat,1);

    nA_vec = NaN(nFeat,1);
    nB_vec = NaN(nFeat,1);

    for j = 1:nFeat

        x = double(X(:,j));
        ok = isfinite(x) & ~isundefined(y);

        xA = x(ok & y == 'A');
        xB = x(ok & y == 'B');

        nA = numel(xA);
        nB = numel(xB);

        nA_vec(j) = nA;
        nB_vec(j) = nB;

        if nA < 5 || nB < 5
            continue;
        end

        if std([xA; xB],0,'omitnan') <= eps
            continue;
        end

        try
            [p, ~, stats] = ranksum(xA, xB);

            pval(j) = p;

            if isfield(stats, 'zval')
                zval(j) = stats.zval;
            end

            xBoth = [xA; xB];
            ranks = tiedrank(xBoth);

            rA = ranks(1:nA);
            rB = ranks(nA+1:end);

            U = sum(rA) - nA*(nA+1)/2;

            U_A(j) = U;

            auc = U / (nA*nB);

            aucA(j) = auc;

            rankBiserial(j) = 2*auc - 1;

            med_A(j) = median(xA, 'omitnan');
            med_B(j) = median(xB, 'omitnan');

            meanRank_A(j) = mean(rA, 'omitnan');
            meanRank_B(j) = mean(rB, 'omitnan');

        catch
            pval(j) = NaN;
        end
    end

    q = local_bh_fdr(pval);

    RS = table( ...
        featureNames, local_feature_family(featureNames), ...
        pval, q, zval, U_A, aucA, rankBiserial, ...
        med_A, med_B, meanRank_A, meanRank_B, nA_vec, nB_vec, ...
        'VariableNames', {'Feature','Family','p_ranksum','q_FDR','z_ranksum','U_A','AUC_A','RankBiserial_A_vs_B', ...
        'Median_A','Median_B','MeanRank_A','MeanRank_B','N_A_feature','N_B_feature'});
end

function q = local_bh_fdr(p)

    p = double(p(:));
    q = NaN(size(p));

    ok = isfinite(p) & p >= 0 & p <= 1;

    p0 = p(ok);
    m = numel(p0);

    if m == 0
        return;
    end

    [ps, ord] = sort(p0, 'ascend');

    qs = ps .* m ./ (1:m)';

    for i = m-1:-1:1
        qs(i) = min(qs(i), qs(i+1));
    end

    qs(qs > 1) = 1;

    tmp = NaN(m,1);
    tmp(ord) = qs;

    q(ok) = tmp;
end

function Summary = local_summary_kw(KW_All, alpha)

    if isempty(KW_All)
        Summary = table();
        return;
    end

    [G, keys] = findgroups(KW_All(:, {'Run','Phase'}));

    N_features = splitapply(@numel, KW_All.p_KW, G);
    N_sig_raw = splitapply(@(p) sum(p < 0.05, 'omitnan'), KW_All.p_KW, G);
    N_sig_FDR = splitapply(@(q) sum(q < alpha, 'omitnan'), KW_All.q_FDR, G);
    Min_p = splitapply(@(p) min(p, [], 'omitnan'), KW_All.p_KW, G);
    Min_q = splitapply(@(q) min(q, [], 'omitnan'), KW_All.q_FDR, G);
    Max_epsilon2 = splitapply(@(x) max(x, [], 'omitnan'), KW_All.epsilon2_KW, G);

    Summary = [keys, table(N_features, N_sig_raw, N_sig_FDR, Min_p, Min_q, Max_epsilon2)];
end

function Summary = local_summary_rs(RS_All, alpha)

    if isempty(RS_All)
        Summary = table();
        return;
    end

    [G, keys] = findgroups(RS_All(:, {'Contrast','Run','Phase'}));

    N_features = splitapply(@numel, RS_All.p_ranksum, G);
    N_sig_raw = splitapply(@(p) sum(p < 0.05, 'omitnan'), RS_All.p_ranksum, G);
    N_sig_FDR = splitapply(@(q) sum(q < alpha, 'omitnan'), RS_All.q_FDR, G);
    Min_p = splitapply(@(p) min(p, [], 'omitnan'), RS_All.p_ranksum, G);
    Min_q = splitapply(@(q) min(q, [], 'omitnan'), RS_All.q_FDR, G);
    MaxAbsRankBiserial = splitapply(@(x) max(abs(x), [], 'omitnan'), RS_All.RankBiserial_A_vs_B, G);

    Summary = [keys, table(N_features, N_sig_raw, N_sig_FDR, Min_p, Min_q, MaxAbsRankBiserial)];
end

function FamilySummary = local_family_summary(T, pcol, alpha)

    if isempty(T)
        FamilySummary = table();
        return;
    end

    if strcmp(pcol, 'p_KW')
        groupCols = {'Run','Phase','Family'};
    else
        groupCols = {'Contrast','Run','Phase','Family'};
    end

    [G, keys] = findgroups(T(:, groupCols));

    p = T.(pcol);
    q = T.q_FDR;

    N_features = splitapply(@numel, p, G);
    N_sig_raw = splitapply(@(x) sum(x < 0.05, 'omitnan'), p, G);
    N_sig_FDR = splitapply(@(x) sum(x < alpha, 'omitnan'), q, G);
    Min_p = splitapply(@(x) min(x, [], 'omitnan'), p, G);
    Min_q = splitapply(@(x) min(x, [], 'omitnan'), q, G);

    FamilySummary = [keys, table(N_features, N_sig_raw, N_sig_FDR, Min_p, Min_q)];

    FamilySummary = sortrows(FamilySummary, {'N_sig_FDR','Min_q'}, {'descend','ascend'});
end

function vname = local_find_var(names, candidates)

    vname = '';

    for i = 1:numel(candidates)
        idx = strcmp(names, candidates{i});

        if any(idx)
            vname = names{find(idx,1)};
            return;
        end
    end
end

function featureNames = local_detect_feature_columns(T)

    allNames = T.Properties.VariableNames;
    isNum = false(1,numel(allNames));

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
end

function s = local_clean_str(x)

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

function fam = local_feature_family(featureNames)

    s = lower(string(featureNames));
    fam = repmat("unknown", size(s));

    isConn = contains(s,"rie") | contains(s,"pli") | contains(s,"plv") | ...
             contains(s,"coh") | contains(s,"conn") | ...
             contains(s,"c01") | contains(s,"c02") | contains(s,"c03") | contains(s,"_c");

    isRBP = contains(s,"rbp") | contains(s,"relative");

    isBP = contains(s,"bp_") | contains(s,"bandpower") | ...
           contains(s,"alpha") | contains(s,"beta") | contains(s,"theta") | ...
           contains(s,"delta") | contains(s,"gamma");

    isTemporal = contains(s,"skew") | contains(s,"kurt") | contains(s,"hjorth") | ...
                 contains(s,"line") | contains(s,"auc") | contains(s,"rms") | ...
                 contains(s,"mean") | contains(s,"std") | contains(s,"var") | ...
                 contains(s,"max") | contains(s,"min");

    isComplexity = contains(s,"lz") | contains(s,"lzc") | contains(s,"entropy") | ...
                   contains(s,"samp") | contains(s,"perm") | contains(s,"higuchi") | ...
                   contains(s,"fractal");

    isTF = contains(s,"morlet") | contains(s,"tf_") | contains(s,"wavelet") | contains(s,"tfr");

    fam(isBP) = "bandpower";
    fam(isRBP) = "relative_bandpower";
    fam(isTemporal) = "temporal_statistical";
    fam(isComplexity) = "complexity";
    fam(isTF) = "time_frequency";
    fam(isConn) = "connectivity";
end

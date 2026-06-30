function RUN10_STEP215_MAKE_DECODING_PRESENTATION_FIGURES()

cfg = get_project_config();

clc; close all;

%% ========================= USER CONFIG =========================
cfg.root = cfg.outputRoot;

cfg.outDir = fullfile(cfg.root, 'DECODING_FIGURES_PRESENTATION');
if ~exist(cfg.outDir, 'dir'); mkdir(cfg.outDir); end

% Main metric for decoding figures
cfg.metricPriority = {'BalAcc','BalancedAccuracy','balanced_accuracy','MeanBalAcc', ...
                      'TestBalAcc','BA','AUC','Accuracy','Acc'};

cfg.primaryRun        = 2;
cfg.primaryPhaseKey   = 'retr';          % retrieval
cfg.primaryContrastA  = 'orientation';
cfg.primaryContrastB  = 'conjunction';
cfg.primaryClassifier = 'LDA';
cfg.primaryK          = 4;

% Chance levels
cfg.chanceBinary = 0.50;
cfg.chanceThreeClass = 1/3;

% Figure style
cfg.fontName = 'Arial';
cfg.fontSize = 12;
cfg.bigFontSize = 16;
cfg.lineWidth = 2.6;
cfg.markerSize = 8;
cfg.pngResolution = 350;

cfg.allowManualSelectedFeatureFallback = false;

%% ========================= MAIN =========================
fprintf('\n=== STEP215: Making decoding presentation figures ===\n');
fprintf('Root: %s\n', cfg.root);
fprintf('Output: %s\n\n', cfg.outDir);

[Perf, perfSource] = localFindPerformanceTable(cfg.root, cfg);
fprintf('Performance table source:\n%s\n\n', perfSource);

cols = localDetectPerformanceColumns(Perf, cfg);
if isempty(cols.metric)
    error('No performance metric column found. Expected something like BalAcc, BalancedAccuracy, AUC, or Accuracy.');
end
fprintf('Using metric column: %s\n', cols.metric);

bestInfo = localChoosePrimaryOrBestRow(Perf, cols, cfg);
fprintf('\nBest/primary row used for detailed figures:\n');
disp(bestInfo.summaryTable);

localPlotAllCellsHeatmap(Perf, cols, cfg);
localPlotKCurveForBestCell(Perf, cols, bestInfo, cfg);

[permVals, obsValue, pValue, permSource] = localFindPermutationVector(cfg.root, cfg, bestInfo);
if ~isempty(permVals)
    fprintf('\nPermutation source:\n%s\n', permSource);
    localPlotPermutationHistogram(permVals, obsValue, pValue, bestInfo, cfg);
else
    localWriteMissingNote(cfg.outDir, '03_missing_permutation_histogram.txt', ...
        ['Permutation vector was not found automatically. ', ...
         'Put the permutation result CSV/MAT in cfg.root or send me the file name.']);
    warning('Permutation vector was not found. A note file was written instead.');
end

[CM, classNames, cmSource] = localFindConfusionMatrix(cfg.root, cfg, bestInfo);
if ~isempty(CM)
    fprintf('\nConfusion matrix source:\n%s\n', cmSource);
    localPlotConfusionMatrix(CM, classNames, bestInfo, cfg);
else
    localWriteMissingNote(cfg.outDir, '04_missing_confusion_matrix.txt', ...
        ['Prediction table or TP/TN/FP/FN columns were not found automatically. ', ...
         'Expected true/predicted labels or confusion counts.']);
    warning('Confusion matrix data was not found. A note file was written instead.');
end

[Feat, featSource] = localFindSelectedFeatures(cfg.root, cfg, bestInfo);
if ~isempty(Feat)
    fprintf('\nSelected-feature source:\n%s\n', featSource);
    localPlotSelectedFeaturesAndChannels(Feat, cfg, bestInfo);
else
    localWriteMissingNote(cfg.outDir, '05_missing_selected_features.txt', ...
        ['Selected-feature/stability table was not found automatically. ', ...
         'Expected columns like FeatureName and Count/Frequency/Stability.']);
    warning('Selected-feature table was not found. A note file was written instead.');
end

fprintf('\nDone. Figures saved in:\n%s\n', cfg.outDir);
fprintf('Recommended files for PowerPoint: the PNG files.\n\n');

end

%% ========================================================================
function [Tbest, sourceFile] = localFindPerformanceTable(rootDir, cfg)

files = localRecursiveFiles(rootDir, {'.csv','.xlsx','.mat'});
bestScore = -Inf;
Tbest = table();
sourceFile = '';

for i = 1:numel(files)
    f = files{i};
    [~, nm, ext] = fileparts(f);
    lowName = lower(nm);

    % Avoid raw/very large subject files when possible
    if contains(lowName, 'subject') && ~contains(lowName, 'result') && ~contains(lowName, 'loso')
        continue;
    end

    try
        tables = localReadTablesFromFile(f, ext);
    catch
        continue;
    end

    for j = 1:numel(tables)
        T = tables{j};
        if ~istable(T) || height(T) < 1 || width(T) < 2
            continue;
        end
        score = localPerformanceTableScore(T, cfg);
        if score > bestScore
            bestScore = score;
            Tbest = T;
            sourceFile = f;
        end
    end
end

if isempty(Tbest) || bestScore < 3
    error(['Could not find a decoding performance table automatically. ', ...
           'Put your LOSO/decoding summary CSV or MAT file inside cfg.root.']);
end

end

%% ========================================================================
function score = localPerformanceTableScore(T, cfg)
vars = lower(string(T.Properties.VariableNames));
clean = regexprep(vars, '[^a-z0-9]', '');

score = 0;

metricPats = cell(size(cfg.metricPriority));
for ii = 1:numel(cfg.metricPriority)
    metricPats{ii} = lower(regexprep(cfg.metricPriority{ii}, '[^a-zA-Z0-9]', ''));
end

if any(localContainsAny(clean, metricPats))
    score = score + 3;
end
if any(strcmp(clean,'k') | contains(clean,'numfeatures') | contains(clean,'nfeatures') | contains(clean,'topk'))
    score = score + 2;
end
if any(contains(clean,'classifier') | contains(clean,'model') | contains(clean,'learner'))
    score = score + 2;
end
if any(contains(clean,'phase') | contains(clean,'epoch') | contains(clean,'window'))
    score = score + 2;
end
if any(strcmp(clean,'run') | contains(clean,'runid'))
    score = score + 2;
end
if any(contains(clean,'contrast') | contains(clean,'comparison') | contains(clean,'task'))
    score = score + 2;
end
if any(contains(clean,'loso') | contains(clean,'fold'))
    score = score + 1;
end
if height(T) > 20
    score = score + 1;
end
end

%% ========================================================================
function cols = localDetectPerformanceColumns(T, cfg)

vars = string(T.Properties.VariableNames);
clean = lower(regexprep(vars, '[^a-z0-9]', ''));

cols.metric = '';
for i = 1:numel(cfg.metricPriority)
    pat = lower(regexprep(cfg.metricPriority{i}, '[^a-z0-9]', ''));
    idx = find(strcmp(clean, pat) | contains(clean, pat), 1, 'first');
    if ~isempty(idx)
        cols.metric = char(vars(idx));
        break;
    end
end

cols.K          = localFindColumn(vars, {'^K$','numfeatures','nfeatures','topk','numselected','nselected'});
cols.classifier = localFindColumn(vars, {'classifier','model','learner','clf'});
cols.phase      = localFindColumn(vars, {'phase','epoch','window'});
cols.run        = localFindColumn(vars, {'^Run$','runid','runnum'});
cols.contrast   = localFindColumn(vars, {'contrast','comparison','taskcontrast','pair','conditionpair'});
cols.auc        = localFindColumn(vars, {'auc'});
cols.acc        = localFindColumn(vars, {'accuracy','acc'});
cols.p          = localFindColumn(vars, {'pvalue','pval','pperm','pfwer','p_fwer'});
cols.true       = localFindColumn(vars, {'ytrue','true','actual','labeltrue','target'});
cols.pred       = localFindColumn(vars, {'ypred','pred','predicted','labelpred','prediction'});

end

%% ========================================================================
function bestInfo = localChoosePrimaryOrBestRow(T, cols, cfg)

metric = cols.metric;
mask = true(height(T),1);

if ~isempty(cols.run)
    runVals = localToDouble(T.(cols.run));
    if any(~isnan(runVals))
        mask = mask & (runVals == cfg.primaryRun);
    end
end

if ~isempty(cols.phase)
    p = lower(localToString(T.(cols.phase)));
    mask = mask & contains(p, lower(cfg.primaryPhaseKey));
end

if ~isempty(cols.contrast)
    c = lower(localToString(T.(cols.contrast)));
    mask = mask & contains(c, lower(cfg.primaryContrastA)) & contains(c, lower(cfg.primaryContrastB));
end

if ~isempty(cols.classifier)
    m = lower(localToString(T.(cols.classifier)));
    mask = mask & contains(m, lower(cfg.primaryClassifier));
end

if ~isempty(cols.K)
    k = localToDouble(T.(cols.K));
    if any(~isnan(k))
        mask = mask & (k == cfg.primaryK);
    end
end

metricVals = localToDouble(T.(metric));

if any(mask & ~isnan(metricVals))
    idxCandidates = find(mask & ~isnan(metricVals));
    [~, jj] = max(metricVals(idxCandidates));
    idx = idxCandidates(jj);
else
    warning('Exact primary model was not found. Using the best available row in the performance table.');
    [~, idx] = max(metricVals);
end

bestInfo.rowIndex = idx;
bestInfo.metricName = metric;
bestInfo.metricValue = metricVals(idx);

bestInfo.run = localGetCellValue(T, cols.run, idx, '');
bestInfo.phase = localGetCellValue(T, cols.phase, idx, '');
bestInfo.contrast = localGetCellValue(T, cols.contrast, idx, '');
bestInfo.classifier = localGetCellValue(T, cols.classifier, idx, '');
bestInfo.K = localGetCellValue(T, cols.K, idx, '');

bestInfo.summaryTable = T(idx, :);

end

%% ========================================================================
function localPlotAllCellsHeatmap(T, cols, cfg)

metric = cols.metric;
metricVals = localToDouble(T.(metric));

if isempty(cols.run) && isempty(cols.phase) && isempty(cols.contrast)
    warning('Cannot make all-cell heatmap because run/phase/contrast columns were not found.');
    return;
end

cellLabels = localMakeCellLabels(T, cols);
[classLabels, ~] = localGetGroupLabels(T, cols.classifier, 'Model');

[uCells, ~, cellId] = unique(cellLabels, 'stable');
[uClasses, ~, classIdStable] = unique(classLabels, 'stable');

M = nan(numel(uCells), numel(uClasses));
for i = 1:numel(uCells)
    for j = 1:numel(uClasses)
        mask = (cellId == i) & (classIdStable == j) & ~isnan(metricVals);
        if any(mask)
            M(i,j) = max(metricVals(mask));
        end
    end
end

% Sort by run/phase/contrast if possible
[sortIdx, sortedLabels] = localSortCellLabels(uCells);
M = M(sortIdx, :);
uCells = sortedLabels;

fig = figure('Color','w','Position',[50 50 1550 1150]);
ax = axes(fig);
imagesc(M);
colormap(ax, localVividColormap(256));
cb = colorbar;
ylabel(cb, strrep(metric, '_', ' '), 'FontWeight','bold');

mn = min(M(:), [], 'omitnan');
mx = max(M(:), [], 'omitnan');
if isempty(mn) || isnan(mn); mn = cfg.chanceBinary; end
if isempty(mx) || isnan(mx); mx = 1; end
caxis([max(0.30, min(cfg.chanceBinary, mn)-0.02), min(1, max(mx, cfg.chanceBinary+0.05))]);

set(ax, 'XTick', 1:numel(uClasses), 'XTickLabel', uClasses, ...
        'YTick', 1:numel(uCells), 'YTickLabel', uCells, ...
        'TickLabelInterpreter','none', 'FontName', cfg.fontName, 'FontSize', 9, ...
        'LineWidth', 1.2, 'Box','on');
xtickangle(35);
title('Decoding summary across all cells', 'FontSize', cfg.bigFontSize, ...
      'FontWeight','bold', 'FontName', cfg.fontName);
xlabel('Classifier', 'FontSize', cfg.fontSize, 'FontWeight','bold');
ylabel('Run | Phase | Contrast', 'FontSize', cfg.fontSize, 'FontWeight','bold');

% Add values inside cells
for i = 1:size(M,1)
    for j = 1:size(M,2)
        if ~isnan(M(i,j))
            txt = sprintf('%.2f', M(i,j));
            if M(i,j) > (mn+mx)/2
                tc = [0 0 0];
            else
                tc = [1 1 1];
            end
            text(j, i, txt, 'HorizontalAlignment','center', 'FontSize', 8, ...
                 'FontWeight','bold', 'Color', tc, 'FontName', cfg.fontName);
        end
    end
end

localSaveFigure(fig, cfg.outDir, '01_decoding_all_cells_heatmap', cfg);
localWriteHeatmapCSV(M, uCells, uClasses, cfg.outDir);
close(fig);

end

%% ========================================================================
function localPlotKCurveForBestCell(T, cols, bestInfo, cfg)

if isempty(cols.K)
    warning('Cannot make K curve because K / number-of-features column was not found.');
    return;
end

metric = cols.metric;
metricVals = localToDouble(T.(metric));
Kvals = localToDouble(T.(cols.K));

mask = true(height(T),1);
mask = localSameValueMask(mask, T, cols.run, bestInfo.run);
mask = localSameValueMask(mask, T, cols.phase, bestInfo.phase);
mask = localSameValueMask(mask, T, cols.contrast, bestInfo.contrast);

if ~any(mask)
    warning('Could not filter data for K curve. Using all rows.');
    mask = true(height(T),1);
end

[classLabels, ~] = localGetGroupLabels(T, cols.classifier, 'Model');
uClasses = unique(classLabels(mask), 'stable');

fig = figure('Color','w','Position',[100 100 1250 760]);
ax = axes(fig); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
colors = lines(max(numel(uClasses), 4));

bestY = -Inf; bestX = NaN; bestClass = '';
for j = 1:numel(uClasses)
    thisMask = mask & strcmp(classLabels, uClasses{j}) & ~isnan(Kvals) & ~isnan(metricVals);
    if ~any(thisMask); continue; end

    ku = unique(Kvals(thisMask));
    ku = sort(ku(:));
    y = nan(size(ku));
    for ii = 1:numel(ku)
        y(ii) = max(metricVals(thisMask & Kvals == ku(ii)));
    end

    plot(ax, ku, y, '-o', 'LineWidth', cfg.lineWidth, ...
         'MarkerSize', cfg.markerSize, 'Color', colors(j,:), ...
         'MarkerFaceColor', colors(j,:), 'DisplayName', uClasses{j});

    [yj, ij] = max(y);
    if yj > bestY
        bestY = yj; bestX = ku(ij); bestClass = uClasses{j};
    end
end

yline(cfg.chanceBinary, '--', 'Chance', 'LineWidth', 1.8, ...
      'LabelHorizontalAlignment','left', 'FontWeight','bold');

if ~isnan(bestX)
    scatter(ax, bestX, bestY, 170, 'p', 'filled', 'MarkerEdgeColor','k', ...
            'MarkerFaceColor',[1 0.85 0.1], 'DisplayName','Best point');
    text(bestX, bestY+0.015, sprintf('  Best: %s, K=%g, %.3f', bestClass, bestX, bestY), ...
         'FontSize', cfg.fontSize, 'FontWeight','bold', 'FontName', cfg.fontName);
end

xlabel('Number of selected features (K)', 'FontSize', cfg.fontSize, 'FontWeight','bold');
ylabel(strrep(metric, '_', ' '), 'FontSize', cfg.fontSize, 'FontWeight','bold');
title(sprintf('Feature-count curve for selected decoding cell\n%s', localBestCellText(bestInfo)), ...
      'FontSize', cfg.bigFontSize, 'FontWeight','bold', 'Interpreter','none');
legend('Location','bestoutside', 'Interpreter','none');
set(ax, 'FontName', cfg.fontName, 'FontSize', cfg.fontSize, 'LineWidth', 1.2);
ylow = min(metricVals(mask),[],'omitnan');
yhigh = max(metricVals(mask),[],'omitnan');
if isempty(ylow) || isnan(ylow); ylow = 0.3; end
if isempty(yhigh) || isnan(yhigh); yhigh = 0.9; end
ylim([max(0.25, ylow-0.05), min(1, yhigh+0.08)]);

localSaveFigure(fig, cfg.outDir, '02_K_vs_balanced_accuracy_best_cell', cfg);
close(fig);

end

%% ========================================================================
function [permVals, obsValue, pValue, sourceFile] = localFindPermutationVector(rootDir, cfg, bestInfo)

permVals = [];
sourceFile = '';
obsValue = bestInfo.metricValue;
pValue = NaN;

files = localRecursiveFiles(rootDir, {'.csv','.mat','.xlsx'});
bestScore = -Inf;
bestVals = [];
bestFile = '';

for i = 1:numel(files)
    f = files{i};
    [~, nm, ext] = fileparts(f);
    lowName = lower(nm);
    if ~(contains(lowName, 'perm') || contains(lowName, 'shuffle') || contains(lowName, 'null'))
        continue;
    end

    try
        tables = localReadTablesFromFile(f, ext);
    catch
        tables = {};
    end

    % Tables
    for j = 1:numel(tables)
        T = tables{j};
        if ~istable(T) || height(T) < 20; continue; end
        vars = string(T.Properties.VariableNames);
        clean = lower(regexprep(vars, '[^a-z0-9]', ''));

        for c = 1:numel(vars)
            vals = localToDouble(T.(vars(c)));
            vals = vals(isfinite(vals));
            if numel(vals) < 20; continue; end
            medv = median(vals, 'omitnan');
            if medv < 0 || medv > 1.2; continue; end

            nameScore = 0;
            if contains(clean(c),'perm') || contains(clean(c),'null') || contains(clean(c),'shuffle'); nameScore = nameScore + 3; end
            if contains(clean(c),'balacc') || contains(clean(c),'balanced') || contains(clean(c),'auc') || contains(clean(c),'accuracy'); nameScore = nameScore + 2; end
            score = nameScore + log(numel(vals));

            if score > bestScore
                bestScore = score;
                bestVals = vals(:);
                bestFile = f;
            end
        end
    end

    % MAT numeric vectors
    if strcmpi(ext, '.mat')
        try
            S = load(f);
            names = fieldnames(S);
            for j = 1:numel(names)
                v = S.(names{j});
                if isnumeric(v) && numel(v) >= 20
                    vals = v(:);
                    vals = vals(isfinite(vals));
                    if numel(vals) < 20; continue; end
                    medv = median(vals, 'omitnan');
                    if medv < 0 || medv > 1.2; continue; end

                    nm2 = lower(regexprep(names{j}, '[^a-z0-9]', ''));
                    nameScore = 0;
                    if contains(nm2,'perm') || contains(nm2,'null') || contains(nm2,'shuffle'); nameScore = nameScore + 3; end
                    if contains(nm2,'balacc') || contains(nm2,'balanced') || contains(nm2,'auc') || contains(nm2,'accuracy'); nameScore = nameScore + 2; end
                    score = nameScore + log(numel(vals));
                    if score > bestScore
                        bestScore = score;
                        bestVals = vals(:);
                        bestFile = [f ' :: ' names{j}];
                    end
                end
            end
        catch
        end
    end
end

if ~isempty(bestVals)
    permVals = bestVals(:);
    sourceFile = bestFile;
    pValue = (sum(permVals >= obsValue) + 1) / (numel(permVals) + 1);
end

end

%% ========================================================================
function localPlotPermutationHistogram(permVals, obsValue, pValue, bestInfo, cfg)

fig = figure('Color','w','Position',[120 120 1150 720]);
ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');

histogram(ax, permVals, 35, 'FaceColor',[0.20 0.45 0.95], ...
          'EdgeColor','none', 'FaceAlpha',0.88);
yl = ylim(ax);

xline(obsValue, '-', sprintf(' Observed = %.3f', obsValue), ...
      'LineWidth', 3.2, 'Color',[0.90 0.15 0.15], ...
      'LabelVerticalAlignment','middle', 'FontWeight','bold');

xline(cfg.chanceBinary, '--', ' Chance', 'LineWidth', 2.0, 'Color',[0.1 0.1 0.1], ...
      'LabelVerticalAlignment','bottom', 'FontWeight','bold');

text(ax, 0.98, 0.92, sprintf('Permutation p = %.4f\nN = %d', pValue, numel(permVals)), ...
     'Units','normalized', 'HorizontalAlignment','right', 'VerticalAlignment','top', ...
     'FontSize', cfg.bigFontSize, 'FontWeight','bold', 'FontName', cfg.fontName, ...
     'BackgroundColor','w', 'EdgeColor',[0.8 0.8 0.8], 'Margin',10);

xlabel('Permutation performance', 'FontSize', cfg.fontSize, 'FontWeight','bold');
ylabel('Count', 'FontSize', cfg.fontSize, 'FontWeight','bold');
title(sprintf('Permutation test for the best decoding model\n%s', localBestCellText(bestInfo)), ...
      'FontSize', cfg.bigFontSize, 'FontWeight','bold', 'Interpreter','none');
set(ax, 'FontName', cfg.fontName, 'FontSize', cfg.fontSize, 'LineWidth', 1.2);
ylim(yl);

localSaveFigure(fig, cfg.outDir, '03_permutation_histogram_best_model', cfg);
close(fig);

end

%% ========================================================================
function [CM, classNames, sourceFile] = localFindConfusionMatrix(rootDir, cfg, bestInfo)

CM = [];
classNames = {};
sourceFile = '';

files = localRecursiveFiles(rootDir, {'.csv','.xlsx','.mat'});
bestScore = -Inf;

for i = 1:numel(files)
    f = files{i};
    [~, nm, ext] = fileparts(f);
    lowName = lower(nm);
    if ~(contains(lowName,'pred') || contains(lowName,'conf') || contains(lowName,'fold'))
        continue;
    end

    try
        tables = localReadTablesFromFile(f, ext);
    catch
        continue;
    end

    for j = 1:numel(tables)
        T = tables{j};
        if ~istable(T) || height(T) < 2; continue; end

        vars = string(T.Properties.VariableNames);
        trueCol = localFindColumn(vars, {'ytrue','true','actual','target','label'});
        predCol = localFindColumn(vars, {'ypred','predicted','prediction','pred'});

        if isempty(trueCol) || isempty(predCol) || strcmp(trueCol, predCol)
            [cm2, names2] = localCMFromCounts(T);
            if ~isempty(cm2)
                score = 5 + height(T);
                if score > bestScore
                    bestScore = score;
                    CM = cm2; classNames = names2; sourceFile = f;
                end
            end
            continue;
        end

        TT = T;
        ccols = localDetectPerformanceColumns(T, cfg);
        mask = true(height(TT),1);
        mask = localSameValueMask(mask, TT, ccols.run, bestInfo.run);
        mask = localSameValueMask(mask, TT, ccols.phase, bestInfo.phase);
        mask = localSameValueMask(mask, TT, ccols.contrast, bestInfo.contrast);
        mask = localSameValueMask(mask, TT, ccols.classifier, bestInfo.classifier);
        mask = localSameValueMask(mask, TT, ccols.K, bestInfo.K);

        if sum(mask) < 2
            mask = true(height(TT),1);
        end

        yTrue = localToString(TT.(trueCol));
        yPred = localToString(TT.(predCol));
        yTrue = yTrue(mask);
        yPred = yPred(mask);

        valid = ~(strcmp(yTrue,'') | strcmp(yPred,''));
        yTrue = yTrue(valid); yPred = yPred(valid);
        if numel(yTrue) < 2; continue; end

        labels = unique([yTrue; yPred], 'stable');
        cm = zeros(numel(labels));
        for a = 1:numel(yTrue)
            r = find(strcmp(labels, yTrue{a}),1);
            c = find(strcmp(labels, yPred{a}),1);
            cm(r,c) = cm(r,c) + 1;
        end

        score = 10 + numel(yTrue);
        if score > bestScore
            bestScore = score;
            CM = cm;
            classNames = labels;
            sourceFile = f;
        end
    end
end

end

%% ========================================================================
function localPlotConfusionMatrix(CM, classNames, bestInfo, cfg)

rowSums = sum(CM, 2);
CMpct = CM ./ max(rowSums, 1) * 100;

fig = figure('Color','w','Position',[160 160 850 760]);
ax = axes(fig);
imagesc(CMpct);
colormap(ax, localVividColormap(256));
caxis([0 100]);
cb = colorbar;
ylabel(cb, 'Row-normalized (%)', 'FontWeight','bold');

n = size(CM,1);
set(ax, 'XTick', 1:n, 'YTick', 1:n, ...
        'XTickLabel', classNames, 'YTickLabel', classNames, ...
        'TickLabelInterpreter','none', 'FontName', cfg.fontName, ...
        'FontSize', cfg.fontSize, 'LineWidth',1.2, 'Box','on');
xtickangle(35);
xlabel('Predicted class', 'FontSize', cfg.fontSize, 'FontWeight','bold');
ylabel('True class', 'FontSize', cfg.fontSize, 'FontWeight','bold');
title(sprintf('Confusion matrix for the best decoding model\n%s', localBestCellText(bestInfo)), ...
      'FontSize', cfg.bigFontSize, 'FontWeight','bold', 'Interpreter','none');

for r = 1:n
    for c = 1:n
        val = CMpct(r,c);
        if val > 50
            tc = [0 0 0];
        else
            tc = [1 1 1];
        end
        text(c, r, sprintf('%d\n%.1f%%', CM(r,c), val), ...
             'HorizontalAlignment','center', 'FontWeight','bold', ...
             'FontSize', cfg.fontSize, 'Color', tc, 'FontName', cfg.fontName);
    end
end

axis square;
localSaveFigure(fig, cfg.outDir, '04_confusion_matrix_best_model', cfg);
close(fig);

end

%% ========================================================================
function [Feat, sourceFile] = localFindSelectedFeatures(rootDir, cfg, bestInfo)

Feat = table();
sourceFile = '';

files = localRecursiveFiles(rootDir, {'.csv','.xlsx','.mat'});
bestScore = -Inf;
bestFeat = table();
bestFile = '';

for i = 1:numel(files)
    f = files{i};
    [~, nm, ext] = fileparts(f);
    lowName = lower(nm);
    if ~(contains(lowName,'feature') || contains(lowName,'selected') || contains(lowName,'stability') || contains(lowName,'fisher'))
        continue;
    end

    try
        tables = localReadTablesFromFile(f, ext);
    catch
        continue;
    end

    for j = 1:numel(tables)
        T = tables{j};
        if ~istable(T) || height(T) < 1; continue; end
        vars = string(T.Properties.VariableNames);

        featCol = localFindColumn(vars, {'featurename','feature','selectedfeature','unit','name'});
        countCol = localFindColumn(vars, {'selectioncount','selectedcount','count','frequency','freq','nfolds','foldcount','stability'});

        if isempty(featCol)
            continue;
        end

        ccols = localDetectPerformanceColumns(T, cfg);
        mask = true(height(T),1);
        mask = localSameValueMask(mask, T, ccols.run, bestInfo.run);
        mask = localSameValueMask(mask, T, ccols.phase, bestInfo.phase);
        mask = localSameValueMask(mask, T, ccols.contrast, bestInfo.contrast);
        mask = localSameValueMask(mask, T, ccols.classifier, bestInfo.classifier);
        mask = localSameValueMask(mask, T, ccols.K, bestInfo.K);

        if sum(mask) < 1
            mask = true(height(T),1);
        end

        F = table();
        F.Feature = localToString(T.(featCol));
        if ~isempty(countCol)
            F.Count = localToDouble(T.(countCol));
        else
            F.Count = ones(height(T),1);
        end
        F = F(mask,:);
        F = F(~strcmp(F.Feature,''),:);

        if isempty(F); continue; end

        F.Channel = localExtractChannelNumber(F.Feature);
        F.Family = localExtractFeatureFamily(F.Feature);

        score = 5 + height(F);
        if ~isempty(countCol); score = score + 3; end
        if score > bestScore
            bestScore = score;
            bestFeat = F;
            bestFile = f;
        end
    end
end

if ~isempty(bestFeat)
    Feat = localAggregateFeatures(bestFeat);
    sourceFile = bestFile;
    return;
end

if cfg.allowManualSelectedFeatureFallback
    Feat = table();
    Feat.Feature = {'ch26 RBP_delta'; 'ch10 RBP_alpha'; 'ch9 RBP_delta'; 'ch10 RBP_delta'};
    Feat.Count = [22; 22; 22; 22];
    Feat.Channel = [26; 10; 9; 10];
    Feat.Family = {'RBP delta'; 'RBP alpha'; 'RBP delta'; 'RBP delta'};
    sourceFile = 'Manual fallback from final primary report';
end

end

%% ========================================================================
function localPlotSelectedFeaturesAndChannels(Feat, cfg, bestInfo)

Feat = sortrows(Feat, 'Count', 'descend');
topN = min(12, height(Feat));
FeatTop = Feat(1:topN,:);

fig = figure('Color','w','Position',[80 80 1500 760]);

% ---- Left: feature stability bar chart
ax1 = subplot(1,2,1);
barh(ax1, FeatTop.Count, 'FaceColor',[0.2 0.55 0.95], 'EdgeColor','none');
set(ax1, 'YTick', 1:topN, 'YTickLabel', FeatTop.Feature, ...
         'YDir','reverse', 'TickLabelInterpreter','none', ...
         'FontName', cfg.fontName, 'FontSize', 10, 'LineWidth', 1.2);
xlabel('Selection count / stability', 'FontSize', cfg.fontSize, 'FontWeight','bold');
title('Most stable selected features', 'FontSize', cfg.bigFontSize, 'FontWeight','bold');
grid(ax1,'on'); box(ax1,'on');

% ---- Right: schematic channel map
ax2 = subplot(1,2,2); hold(ax2,'on'); axis(ax2,'equal'); axis(ax2,'off');

theta = linspace(0, 2*pi, 250);
plot(ax2, cos(theta), sin(theta), 'k-', 'LineWidth', 2.2);
plot(ax2, [0 0.12 -0.12 0], [1 1.14 1.14 1], 'k-', 'LineWidth', 2.0); % nose
plot(ax2, [-1 -1.08 -1], [0.18 0 -0.18], 'k-', 'LineWidth', 1.8);
plot(ax2, [1 1.08 1], [0.18 0 -0.18], 'k-', 'LineWidth', 1.8);

[x, y] = localSchematic64Layout();
scatter(ax2, x, y, 42, [0.70 0.70 0.70], 'filled', 'MarkerEdgeColor',[0.45 0.45 0.45]);

validCh = Feat.Channel(~isnan(Feat.Channel));
uCh = unique(validCh);
counts = zeros(size(uCh));
for i = 1:numel(uCh)
    counts(i) = sum(Feat.Count(Feat.Channel == uCh(i)), 'omitnan');
end
if ~isempty(uCh)
    maxCount = max(counts);
    for i = 1:numel(uCh)
        ch = uCh(i);
        if ch >= 1 && ch <= numel(x)
            sz = 120 + 330 * counts(i) / max(maxCount,1);
            scatter(ax2, x(ch), y(ch), sz, [0.95 0.20 0.20], 'filled', ...
                    'MarkerEdgeColor','k', 'LineWidth',1.3);
            text(ax2, x(ch), y(ch)+0.075, sprintf('ch%d', ch), ...
                 'HorizontalAlignment','center', 'FontWeight','bold', ...
                 'FontSize', 11, 'FontName', cfg.fontName, 'Color','k');
        end
    end
end

title(ax2, sprintf('Important channels\n%s', localBestCellText(bestInfo)), ...
      'FontSize', cfg.bigFontSize, 'FontWeight','bold', 'Interpreter','none');
text(ax2, 0, -1.20, 'Schematic sensor layout; highlighted size reflects feature-selection stability.', ...
     'HorizontalAlignment','center', 'FontSize', 10, 'FontName', cfg.fontName);

localSaveFigure(fig, cfg.outDir, '05_selected_features_and_important_channels', cfg);
close(fig);

% Save selected feature table
try
    writetable(Feat, fullfile(cfg.outDir, '05_selected_features_table.csv'));
catch
end

end

%% ============================ HELPERS ============================

function files = localRecursiveFiles(rootDir, exts)
files = {};
if ~exist(rootDir, 'dir'); return; end
d = dir(rootDir);
for i = 1:numel(d)
    nm = d(i).name;
    if d(i).isdir
        if strcmp(nm,'.') || strcmp(nm,'..'); continue; end
        if startsWith(nm, '.') || contains(lower(nm), 'private'); continue; end
        files = [files; localRecursiveFiles(fullfile(rootDir,nm), exts)]; %#ok<AGROW>
    else
        [~,~,ext] = fileparts(nm);
        if any(strcmpi(ext, exts))
            files{end+1,1} = fullfile(rootDir,nm); %#ok<AGROW>
        end
    end
end
end

function tables = localReadTablesFromFile(f, ext)
tables = {};
switch lower(ext)
    case '.csv'
        try
            opts = detectImportOptions(f);
            T = readtable(f, opts);
        catch
            T = readtable(f);
        end
        tables{end+1} = T;

    case '.xlsx'
        try
            T = readtable(f);
            tables{end+1} = T;
        catch
        end

    case '.mat'
        S = load(f);
        tables = localExtractTablesFromStruct(S, 0);
end
end

function tables = localExtractTablesFromStruct(S, depth)
tables = {};
if depth > 3; return; end
names = fieldnames(S);
for i = 1:numel(names)
    v = S.(names{i});
    if istable(v)
        tables{end+1} = v; %#ok<AGROW>
    elseif isstruct(v)
        if numel(v) > 1
            try
                T = struct2table(v);
                tables{end+1} = T; %#ok<AGROW>
            catch
            end
        else
            subTables = localExtractTablesFromStruct(v, depth+1);
            tables = [tables; subTables]; %#ok<AGROW>
        end
    end
end
end

function col = localFindColumn(vars, patterns)
col = '';
clean = lower(regexprep(string(vars), '[^a-z0-9]', ''));
for p = 1:numel(patterns)
    pat = lower(regexprep(patterns{p}, '[^a-z0-9]', ''));
    idx = find(strcmp(clean, pat) | contains(clean, pat), 1, 'first');
    if ~isempty(idx)
        col = char(vars(idx));
        return;
    end
end
end

function tf = localContainsAny(strs, patterns)
tf = false(size(strs));
for i = 1:numel(patterns)
    tf = tf | contains(strs, patterns{i});
end
end

function x = localToDouble(v)
if isnumeric(v)
    x = double(v);
elseif islogical(v)
    x = double(v);
elseif iscell(v)
    x = nan(numel(v),1);
    for i = 1:numel(v)
        if isnumeric(v{i}) && isscalar(v{i})
            x(i) = double(v{i});
        else
            x(i) = str2double(string(v{i}));
        end
    end
elseif isstring(v) || ischar(v) || iscategorical(v)
    x = str2double(string(v));
else
    try
        x = str2double(string(v));
    catch
        x = nan(numel(v),1);
    end
end
x = x(:);
end

function s = localToString(v)
if isstring(v)
    s = cellstr(v(:));
elseif iscellstr(v)
    s = v(:);
elseif iscategorical(v)
    s = cellstr(string(v(:)));
elseif isnumeric(v)
    s = cellstr(string(v(:)));
elseif iscell(v)
    s = cell(numel(v),1);
    for i = 1:numel(v)
        try
            s{i} = char(string(v{i}));
        catch
            s{i} = '';
        end
    end
else
    s = cellstr(string(v(:)));
end
end

function val = localGetCellValue(T, col, idx, defaultVal)
if isempty(col)
    val = defaultVal;
    return;
end
v = T.(col);
if isnumeric(v)
    val = v(idx);
else
    ss = localToString(v);
    val = ss{idx};
end
end

function labels = localMakeCellLabels(T, cols)
n = height(T);
r = repmat({''}, n, 1);
p = repmat({''}, n, 1);
c = repmat({''}, n, 1);

if ~isempty(cols.run); r = localToString(T.(cols.run)); end
if ~isempty(cols.phase); p = localToString(T.(cols.phase)); end
if ~isempty(cols.contrast); c = localToString(T.(cols.contrast)); end

labels = cell(n,1);
for i = 1:n
    labels{i} = sprintf('Run %s | %s | %s', r{i}, p{i}, c{i});
    labels{i} = strrep(labels{i}, '_', ' ');
end
end

function [labels, id] = localGetGroupLabels(T, col, defaultLabel)
if isempty(col)
    labels = repmat({defaultLabel}, height(T), 1);
else
    labels = localToString(T.(col));
    labels(strcmp(labels,'')) = {defaultLabel};
end
[~,~,id] = unique(labels, 'stable');
end

function mask = localSameValueMask(mask, T, col, target)
if isempty(col) || isempty(target)
    return;
end
v = T.(col);
if isnumeric(v) && isnumeric(target)
    mask = mask & (double(v(:)) == double(target));
else
    s = lower(localToString(v));
    targetS = lower(char(string(target)));
    if isempty(targetS); return; end
    mask = mask & strcmp(s, targetS);
end
end

function txt = localBestCellText(bestInfo)
txt = sprintf('Run %s | %s | %s | %s | K=%s', ...
    char(string(bestInfo.run)), char(string(bestInfo.phase)), ...
    char(string(bestInfo.contrast)), char(string(bestInfo.classifier)), ...
    char(string(bestInfo.K)));
txt = strrep(txt, '_', ' ');
end

function [sortIdx, sortedLabels] = localSortCellLabels(labels)
n = numel(labels);
runNum = nan(n,1);
phaseRank = nan(n,1);
for i = 1:n
    tok = regexp(labels{i}, 'Run\s*([0-9]+)', 'tokens', 'once');
    if ~isempty(tok); runNum(i) = str2double(tok{1}); end
    low = lower(labels{i});
    if contains(low,'stim')
        phaseRank(i) = 1;
    elseif contains(low,'maint')
        phaseRank(i) = 2;
    elseif contains(low,'retr')
        phaseRank(i) = 3;
    else
        phaseRank(i) = 9;
    end
end
[~, sortIdx] = sortrows([runNum phaseRank (1:n)']);
sortedLabels = labels(sortIdx);
end

function cm = localVividColormap(n)
% Blue -> cyan -> yellow -> red custom map, vivid but clean.
if nargin < 1; n = 256; end
base = [ ...
    0.10 0.10 0.45
    0.00 0.45 0.85
    0.00 0.75 0.75
    0.95 0.90 0.20
    0.95 0.25 0.15];
x = linspace(0,1,size(base,1));
xi = linspace(0,1,n);
cm = interp1(x, base, xi, 'linear');
cm = max(0, min(1, cm));
end

function localSaveFigure(fig, outDir, baseName, cfg)
if ~exist(outDir, 'dir'); mkdir(outDir); end
set(fig, 'PaperPositionMode','auto');

pngFile = fullfile(outDir, [baseName '.png']);
figFile = fullfile(outDir, [baseName '.fig']);
pdfFile = fullfile(outDir, [baseName '.pdf']);

try
    print(fig, pngFile, '-dpng', sprintf('-r%d', cfg.pngResolution));
catch ME
    warning('Could not save PNG: %s', ME.message);
end

try
    savefig(fig, figFile);
catch
end

try
    print(fig, pdfFile, '-dpdf', '-painters');
catch
end

fprintf('Saved: %s\n', pngFile);
end

function localWriteHeatmapCSV(M, rowLabels, colLabels, outDir)
try
    T = array2table(M, 'VariableNames', matlab.lang.makeValidName(colLabels));
    T.Cell = rowLabels(:);
    T = movevars(T, 'Cell', 'Before', 1);
    writetable(T, fullfile(outDir, '01_decoding_all_cells_heatmap_values.csv'));
catch
end
end

function localWriteMissingNote(outDir, fileName, msg)
fid = fopen(fullfile(outDir, fileName), 'w');
if fid > 0
    fprintf(fid, '%s\n', msg);
    fclose(fid);
end
end

function [cm, names] = localCMFromCounts(T)
cm = [];
names = {};
vars = string(T.Properties.VariableNames);
tp = localFindColumn(vars, {'tp','truepositive','truepositives'});
tn = localFindColumn(vars, {'tn','truenegative','truenegatives'});
fp = localFindColumn(vars, {'fp','falsepositive','falsepositives'});
fn = localFindColumn(vars, {'fn','falsenegative','falsenegatives'});
if isempty(tp) || isempty(tn) || isempty(fp) || isempty(fn); return; end

tpv = localToDouble(T.(tp)); tnv = localToDouble(T.(tn));
fpv = localToDouble(T.(fp)); fnv = localToDouble(T.(fn));
tpv = tpv(1); tnv = tnv(1); fpv = fpv(1); fnv = fnv(1);

cm = [tnv fpv; fnv tpv];
names = {'Class 1','Class 2'};
end

function ch = localExtractChannelNumber(featureNames)
ch = nan(numel(featureNames),1);
for i = 1:numel(featureNames)
    s = char(featureNames{i});
    tok = regexp(s, 'ch\s*0*([0-9]+)', 'tokens', 'once', 'ignorecase');
    if isempty(tok)
        tok = regexp(s, 'chan\s*0*([0-9]+)', 'tokens', 'once', 'ignorecase');
    end
    if ~isempty(tok)
        ch(i) = str2double(tok{1});
    end
end
end

function fam = localExtractFeatureFamily(featureNames)
fam = cell(numel(featureNames),1);
for i = 1:numel(featureNames)
    s = lower(char(featureNames{i}));
    if contains(s, 'rbp') || contains(s, 'relative')
        if contains(s,'alpha'); fam{i} = 'RBP alpha';
        elseif contains(s,'delta'); fam{i} = 'RBP delta';
        elseif contains(s,'theta'); fam{i} = 'RBP theta';
        elseif contains(s,'beta'); fam{i} = 'RBP beta';
        elseif contains(s,'gamma'); fam{i} = 'RBP gamma';
        else; fam{i} = 'Relative bandpower';
        end
    elseif contains(s, 'bp') || contains(s, 'bandpower')
        fam{i} = 'Bandpower';
    elseif contains(s, 'rms')
        fam{i} = 'RMS';
    elseif contains(s, 'lzc')
        fam{i} = 'LZC';
    elseif contains(s, 'tkeo')
        fam{i} = 'TKEO';
    elseif contains(s, 'linelen')
        fam{i} = 'Line length';
    else
        fam{i} = 'Other';
    end
end
end

function F = localAggregateFeatures(Fin)
% Aggregate repeated feature rows
feat = Fin.Feature;
[uFeat, ~, id] = unique(feat, 'stable');
cnt = zeros(numel(uFeat),1);
ch = nan(numel(uFeat),1);
fam = cell(numel(uFeat),1);

for i = 1:numel(uFeat)
    mask = id == i;
    cnt(i) = sum(Fin.Count(mask), 'omitnan');
    chVals = Fin.Channel(mask);
    chVals = chVals(~isnan(chVals));
    if ~isempty(chVals); ch(i) = chVals(1); end
    famVals = Fin.Family(mask);
    if ~isempty(famVals); fam{i} = famVals{1}; else; fam{i} = 'Other'; end
end

F = table(uFeat, cnt, ch, fam, 'VariableNames', {'Feature','Count','Channel','Family'});
end

function [x, y] = localSchematic64Layout()
% Schematic 64-channel layout. Not a source-localized anatomical map.
x = nan(64,1); y = nan(64,1);
rings = {1:4, 5:12, 13:24, 25:40, 41:56, 57:64};
radii = [0.18 0.36 0.54 0.72 0.88 0.98];
phaseShift = [pi/4 0 pi/12 0 pi/16 pi/8];
for rr = 1:numel(rings)
    idx = rings{rr};
    th = linspace(0, 2*pi, numel(idx)+1);
    th(end) = [];
    th = th + phaseShift(rr);
    x(idx) = radii(rr) * sin(th);
    y(idx) = radii(rr) * cos(th);
end
end

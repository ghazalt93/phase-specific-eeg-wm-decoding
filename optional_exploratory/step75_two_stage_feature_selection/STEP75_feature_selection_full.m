function [FSmodel, Ranking] = STEP75_feature_selection_full(X, y, featureNames, cfg)


method = lower(string(cfg.featureMethod));
featureNames = string(featureNames(:));
y = removecats(categorical(y));
cls = categories(y);

if numel(cls) ~= 2
    error('Feature selection requires binary labels.');
end

switch method

    case "ranksum"
        [Ranking, idxSel] = local_rank_ranksum(X, y, featureNames, cfg.Kfeat);

        FSmodel = struct();
        FSmodel.method = 'ranksum';
        FSmodel.selectedIdx = idxSel;
        FSmodel.selectedNames = cellstr(featureNames(idxSel));

    case "fisher"
        [Ranking, idxSel] = local_rank_fisher(X, y, featureNames, cfg.Kfeat);

        FSmodel = struct();
        FSmodel.method = 'fisher';
        FSmodel.selectedIdx = idxSel;
        FSmodel.selectedNames = cellstr(featureNames(idxSel));

    case "mrmr"
        yNum = double(y == cls{2});

        try
            [idx, scores] = fscmrmr(X, yNum);
        catch ME
            error('fscmrmr failed: %s', ME.message);
        end

        idx = idx(:);
        idx = idx(isfinite(idx) & idx >= 1 & idx <= size(X,2));
        K = min(cfg.Kfeat, numel(idx));
        idxSel = idx(1:K);

        scoreAll = NaN(numel(idx),1);
        if exist('scores','var') && ~isempty(scores)
            scores = scores(:);
            scoreAll(1:min(numel(scores), numel(scoreAll))) = scores(1:min(numel(scores), numel(scoreAll)));
        end

        Ranking = table((1:numel(idx))', featureNames(idx), scoreAll, ...
            'VariableNames', {'Rank','Feature','MRMRScore'});

        FSmodel = struct();
        FSmodel.method = 'mrmr';
        FSmodel.selectedIdx = idxSel;
        FSmodel.selectedNames = cellstr(featureNames(idxSel));

    case "relieff"
        yNum = double(y == cls{2});
        kNeighbors = min(10, max(1, sum(yNum==1)-1));
        kNeighbors = min(kNeighbors, max(1, sum(yNum==0)-1));

        try
            [idx, weights] = relieff(X, yNum, kNeighbors);
        catch ME
            error('relieff failed: %s', ME.message);
        end

        idx = idx(:);
        idx = idx(isfinite(idx) & idx >= 1 & idx <= size(X,2));
        K = min(cfg.Kfeat, numel(idx));
        idxSel = idx(1:K);

        weights = weights(:);
        wAll = NaN(numel(idx),1);
        wAll(1:min(numel(weights), numel(wAll))) = weights(1:min(numel(weights), numel(wAll)));

        Ranking = table((1:numel(idx))', featureNames(idx), wAll, ...
            'VariableNames', {'Rank','Feature','ReliefFWeight'});

        FSmodel = struct();
        FSmodel.method = 'relieff';
        FSmodel.selectedIdx = idxSel;
        FSmodel.selectedNames = cellstr(featureNames(idxSel));

    case "pca"
        [FSmodel, Ranking] = local_fit_pca_svd(X, featureNames, cfg);

    otherwise
        error('Unknown cfg.featureMethod: %s', cfg.featureMethod);
end

end

%% ===================== LOCAL METHODS =====================

function [Ranking, idxSel] = local_rank_ranksum(X, y, featureNames, Kfeat)

cls = categories(y);
gA = y == cls{1};
gB = y == cls{2};

nFeat = size(X,2);
pval = NaN(nFeat,1);
aucVal = NaN(nFeat,1);
score = NaN(nFeat,1);

for j = 1:nFeat

    x = X(:,j);
    ok = isfinite(x);

    xA = x(gA & ok);
    xB = x(gB & ok);

    if numel(xA) < 5 || numel(xB) < 5
        continue;
    end

    if std([xA; xB],0,'omitnan') <= eps
        continue;
    end

    try
        pval(j) = ranksum(xA, xB);

        xBoth = [xA; xB];
        r = tiedrank(xBoth);
        nA = numel(xA);
        nB = numel(xB);
        rA = r(1:nA);
        U_A = sum(rA) - nA*(nA+1)/2;
        aucA = U_A/(nA*nB);

        aucVal(j) = aucA;
        score(j) = abs(aucA - 0.5);
    catch
    end
end

idxAll = STEP75_topk(score, numel(score));

Ranking = table((1:numel(idxAll))', featureNames(idxAll), score(idxAll), aucVal(idxAll), pval(idxAll), ...
    'VariableNames', {'Rank','Feature','RankAUC_score','UnivariateAUC','Ranksum_p'});

K = min(Kfeat, numel(idxAll));
idxSel = idxAll(1:K);

end

function [Ranking, idxSel] = local_rank_fisher(X, y, featureNames, Kfeat)

cls = categories(y);
gA = y == cls{1};
gB = y == cls{2};

nFeat = size(X,2);
score = NaN(nFeat,1);

for j = 1:nFeat
    xA = X(gA,j);
    xB = X(gB,j);

    denom = var(xA,0,'omitnan') + var(xB,0,'omitnan');

    if ~isfinite(denom) || denom <= eps
        continue;
    end

    score(j) = (mean(xA,'omitnan') - mean(xB,'omitnan')).^2 ./ denom;
end

idxAll = STEP75_topk(score, numel(score));

Ranking = table((1:numel(idxAll))', featureNames(idxAll), score(idxAll), ...
    'VariableNames', {'Rank','Feature','FisherScore'});

K = min(Kfeat, numel(idxAll));
idxSel = idxAll(1:K);

end

function [FSmodel, Ranking] = local_fit_pca_svd(X, featureNames, cfg)

if isfield(cfg,'nPC') && ~isempty(cfg.nPC)
    K = cfg.nPC;
else
    K = cfg.Kfeat;
end

sd = std(X,0,1);
keepVar = isfinite(sd) & sd > 0;

X0 = X(:, keepVar);

if isempty(X0) || size(X0,2) < 2
    error('Too few valid features for PCA.');
end

try
    [~, S, V] = svd(X0, 'econ');
catch ME
    error('SVD failed in PCA: %s', ME.message);
end

eigvals = diag(S).^2 ./ max(size(X0,1)-1, 1);
explained = 100 * eigvals ./ sum(eigvals);

if isfield(cfg,'pcaVarianceToKeep') && ~isempty(cfg.pcaVarianceToKeep)
    nByVar = find(cumsum(explained) >= cfg.pcaVarianceToKeep, 1, 'first');
    if isempty(nByVar)
        nByVar = numel(explained);
    end
else
    nByVar = numel(explained);
end

nComp = min([K, nByVar, size(V,2), size(X0,1)-1, size(X0,2)]);

coeff = V(:, 1:nComp);

pcNames = strings(nComp,1);
for k = 1:nComp
    pcNames(k) = sprintf('PC%03d', k);
end

Ranking = table((1:nComp)', pcNames, explained(1:nComp), cumsum(explained(1:nComp)), ...
    'VariableNames', {'Component','PCName','ExplainedVariancePercent','CumulativeExplainedVariancePercent'});

FSmodel = struct();
FSmodel.method = 'pca_svd';
FSmodel.nComp = nComp;
FSmodel.keepVar = keepVar;
FSmodel.coeff = coeff;
FSmodel.explained = explained(1:nComp);
FSmodel.originalFeatureNamesAfterVarFilter = cellstr(featureNames(keepVar));

end

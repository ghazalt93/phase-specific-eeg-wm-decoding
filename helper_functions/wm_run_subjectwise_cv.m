function [PredT, FoldT, TopFeatureT] = wm_run_subjectwise_cv(D, meta, featVars, featCh, cfg, classKeep, classifierName)
    if nargin < 7 || isempty(classifierName)
        classifierName = 'svm_ecoc';
    end

    [Xraw, y0, subj0] = wm_get_cell_data(D, meta, featVars, cfg, classKeep);

    if isempty(Xraw)
        error('No rows selected for run=%d phase=%s.', cfg.run, cfg.phase);
    end

    Xnorm = wm_apply_subject_norm(Xraw, subj0, cfg.normMethod);

    uSub = unique(subj0, 'stable');
    nFold = min(cfg.nFolds, numel(uSub));
    if nFold < 2
        error('Too few subjects for CV.');
    end

    rng(cfg.randomSeed + cfg.run*1000);
    cv = cvpartition(numel(uSub), 'KFold', nFold);

    predRows = {};
    foldRows = {};
    featRows = {};
    pi = 0; fi = 0; ti = 0;

    for ff = 1:nFold
        testSub = uSub(test(cv,ff));
        te = ismember(subj0, testSub);
        tr = ~te;

        Xtr0 = Xnorm(tr,:);
        Xte0 = Xnorm(te,:);
        ytr0 = y0(tr);
        yte0 = y0(te);

        try
            [Xtr, Xte, keepCol] = wm_prep_train_test(Xtr0, Xte0);
            varsLocal = featVars(keepCol);
            chLocal = featCh(keepCol);

            F = wm_anovaF_multiclass(Xtr, ytr0);
            F(~isfinite(F)) = -Inf;
            [~, ord] = sort(F, 'descend');
            K = min(cfg.Kfeat, numel(ord));
            sel = ord(1:K);

            for jj = 1:numel(sel)
                ti = ti + 1;
                featRows(ti,:) = {ff, jj, string(varsLocal{sel(jj)}), chLocal(sel(jj)), F(sel(jj)), ...
                    string(wm_feature_family(varsLocal{sel(jj)}))}; %#ok<AGROW>
            end

            switch lower(classifierName)
                case 'svm_ecoc'
                    mdl = wm_train_svm_ecoc(Xtr(:,sel), ytr0);
                    yh = double(predict(mdl, Xte(:,sel)));

                case 'lda'
                    mdl = fitcdiscr(Xtr(:,sel), ytr0, 'DiscrimType','linear');
                    yh = double(predict(mdl, Xte(:,sel)));

                case 'logistic_ecoc'
                    try
                        t = templateLinear('Learner','logistic', 'Regularization','ridge', ...
                            'Lambda',1e-4, 'Solver','lbfgs');
                        mdl = fitcecoc(Xtr(:,sel), ytr0, 'Learners', t, 'Coding','onevsone', ...
                            'ClassNames', unique(ytr0(:))');
                        yh = double(predict(mdl, Xte(:,sel)));
                    catch ME
                        warning('Logistic ECOC failed in fold %d: %s', ff, ME.message);
                        yh = nan(size(yte0));
                    end

                otherwise
                    error('Unknown classifierName: %s', classifierName);
            end

            classes = unique([y0(:); yh(:)]);
            classes = classes(isfinite(classes));
            [acc, bal, rec, cm] = wm_acc_bal(yte0, yh, classes);

            fi = fi + 1;
            foldRows(fi,:) = {string(classifierName), ff, numel(ytr0), numel(yte0), numel(unique(subj0(tr))), ...
                numel(unique(subj0(te))), K, acc, bal, string(strjoin(testSub,';')), "ok", ""}; %#ok<AGROW>

            rowIdx = find(te);
            for ii = 1:numel(yte0)
                pi = pi + 1;
                pi2 = rowIdx(ii);
                piPred = NaN;
                if ii <= numel(yh), piPred = yh(ii); end
                predRows(pi,:) = {string(classifierName), ff, string(subj0{pi2}), yte0(ii), piPred, ...
                    cfg.run, string(cfg.phase)}; %#ok<AGROW>
            end

        catch ME
            fi = fi + 1;
            foldRows(fi,:) = {string(classifierName), ff, sum(tr), sum(te), numel(unique(subj0(tr))), ...
                numel(unique(subj0(te))), NaN, NaN, NaN, string(strjoin(testSub,';')), "fail", string(ME.message)}; %#ok<AGROW>
        end
    end

    PredT = cell2table(predRows, 'VariableNames', ...
        {'Classifier','Fold','Subject','TrueLabel','PredLabel','Run','Phase'});

    FoldT = cell2table(foldRows, 'VariableNames', ...
        {'Classifier','Fold','NTrain','NTest','NTrainSubjects','NTestSubjects','Kfeatures','Acc','BalAcc','TestSubjects','Status','Message'});

    TopFeatureT = cell2table(featRows, 'VariableNames', ...
        {'Fold','Rank','Feature','Channel','Score','Family'});
end

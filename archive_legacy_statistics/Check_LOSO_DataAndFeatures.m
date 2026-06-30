function Check_LOSO_DataAndFeatures()

dbclear if caught error
dbclear if error

results_root = '<RESULT_MINE_ROOT>/New Folder/no_ica/GroupResults_LOSO';
phases = {'stim','maint','retr'};

fprintf('\n=== CHECK: loading per-run feature mats ===\n');

for p = 1:numel(phases)
    phase = phases{p};
    fprintf('\n----- PHASE: %s -----\n', phase);

    % find all *_feat.mat for this phase
    files = dir(fullfile(results_root,'s*',sprintf('*_%s_feat.mat',phase)));
    fprintf('Found %d feat files\n', numel(files));
    if isempty(files), continue; end

    allTbl = table();

    for i=1:numel(files)
        f = fullfile(files(i).folder, files(i).name);
        S = load(f,'featAll','trialInfo','cleanInfo');
        if ~isfield(S,'featAll') || isempty(S.featAll) || height(S.featAll)==0
            fprintf('[EMPTY] %s\n', files(i).name);
            continue;
        end

        T = S.featAll;

        % basic meta check
        if ~ismember('Subject', T.Properties.VariableNames)
            warning('No Subject col in %s', files(i).name);
            continue;
        end
        if ~ismember('Condition', T.Properties.VariableNames)
            warning('No Condition col in %s', files(i).name);
            continue;
        end

        % append
        allTbl = [allTbl; T]; %#ok<AGROW>
    end

    if isempty(allTbl) || height(allTbl)==0
        fprintf('No rows collected for %s\n', phase);
        continue;
    end

    allTbl.Subject = string(allTbl.Subject);
    allTbl.Condition = double(allTbl.Condition);

    fprintf('Total rows: %d | Total subjects: %d | Total cols: %d\n', ...
        height(allTbl), numel(unique(allTbl.Subject)), width(allTbl));

    % --- label sanity ---
    bad = sum(~ismember(allTbl.Condition,[1 2 3]) | isnan(allTbl.Condition));
    fprintf('Bad/NaN labels: %d\n', bad);

    % --- class distribution overall ---
    c1 = sum(allTbl.Condition==1);
    c2 = sum(allTbl.Condition==2);
    c3 = sum(allTbl.Condition==3);
    fprintf('Class counts overall: [1]=%d [2]=%d [3]=%d\n', c1,c2,c3);

    % --- per-subject class coverage ---
    subs = unique(allTbl.Subject);
    miss3 = 0;
    for s=1:numel(subs)
        y = allTbl.Condition(allTbl.Subject==subs(s));
        u = unique(y);
        if numel(u) < 3
            miss3 = miss3 + 1;
        end
    end
    fprintf('Subjects missing some classes (<3): %d / %d\n', miss3, numel(subs));

    % --- numeric feature sanity ---
    metaVars = ["Subject","Condition"];
    vnames = string(allTbl.Properties.VariableNames);
    featMask = false(1,numel(vnames));
    for j=1:numel(vnames)
        if any(vnames(j)==metaVars), continue; end
        x = allTbl.(vnames(j));
        if isnumeric(x) || islogical(x)
            featMask(j) = true;
        end
    end
    X = allTbl{:,featMask};
    fprintf('Numeric feature cols: %d\n', size(X,2));

    % NaN ratio
    nanRatio = mean(isnan(X),1);
    fprintf('Median NaN ratio across features: %.3f | Max NaN ratio: %.3f\n', ...
        median(nanRatio), max(nanRatio));

    % zero variance
    vv = var(fillmissing(X,'constant',0),0,1);
    zv = sum(vv==0);
    fprintf('Zero-variance features: %d\n', zv);

    % --- leakage quick scan by name ---
    leakagePatterns = ["trial","run","start","end","sample","time","sec","row","trigger","latency","event"];
    leakCols = strings(0,1);
    for j=1:numel(vnames)
        for k=1:numel(leakagePatterns)
            if contains(lower(vnames(j)), leakagePatterns(k))
                leakCols(end+1,1) = vnames(j); %#ok<AGROW>
                break;
            end
        end
    end
    fprintf('Potential leakage-like columns (name-based): %d\n', numel(unique(leakCols)));

end

fprintf('\n=== DONE CHECK ===\n');
end

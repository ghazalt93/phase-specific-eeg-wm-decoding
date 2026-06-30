%% STEP107B_MAXPERM_primaryFamily_CHANNEL_binary_run2retr_OVC_noLeakage.m

clear; clc;
cfgBase = get_project_config();

cfg = struct();
cfg.ROOT = fullfile(cfgBase.outputRoot, 'new_article_outputs');

cfg.datasetPaths = { ...
    fullfile(cfg.ROOT,'_wm_ml','dataset.mat')
    fullfile(cfg.ROOT,'Subjects','_wm_ml','dataset.mat')
    fullfile(cfg.ROOT,'_wm_dataset','dataset.mat')
    fullfile(cfg.ROOT,'dataset.mat')
};

cfg.nPerm = 1000;
cfg.rngSeed = 2026;

cfg.Klist = 1:22;
cfg.classifiers = {'lda','svm','randomforest','decisiontree'};
cfg.nTrees = 100;

% Primary family
cfg.taskType = 'binary';
cfg.contrastName = 'orientation_vs_conjunction';
cfg.run = 2;
cfg.phase = 'retr';

cfg.labelMap.color       = {'1','color'};
cfg.labelMap.orientation = {'2','orientation'};
cfg.labelMap.conjunction = {'3','conjunction'};

cfg.channelsToUse = 1:64;
cfg.excludeConnectivityFeatures = true;
cfg.unitFeatureAggregate = 'median';
cfg.subjectConditionAggregate = 'median';
cfg.trainStandardize = true;

cfg.outDir = fullfile(cfg.ROOT,'_wm_STEP107B_MAXPERM_primaryFamily_CHANNEL_binary_run2retr_OVC_noLeakage');
if ~exist(cfg.outDir,'dir'), mkdir(cfg.outDir); end

diary(fullfile(cfg.outDir,'STEP107B_log.txt'));

fprintf('\n=== STEP107B max-permutation for primary channel model family ===\n');
fprintf('Started: %s\n', datestr(now));
fprintf('nPerm = %d\n', cfg.nPerm);
fprintf('Task: channel | %s | %s | run%d | %s\n', cfg.taskType, cfg.contrastName, cfg.run, cfg.phase);

%% ===================== LOAD DATA =====================
[T,col,fNames] = load_dataset(cfg);
subjAll  = cleanstr(T.(col.subj));
condAll  = cleanstr(T.(col.cond));
runAll   = double(T.(col.run));
phaseAll = cleanstr(T.(col.phase));

fprintf('Dataset rows=%d | numeric candidate features=%d | subjects=%d\n', height(T),numel(fNames),numel(unique(subjAll)));

%% ===================== FOLDS =====================
[GlobalFoldTable, globalFoldIDs] = create_global_subject_folds(subjAll);
writetable(GlobalFoldTable, fullfile(cfg.outDir,'STEP107B_internal_global_subject_LOO_folds.csv'));
fprintf('Global Subject-LOO folds=%d\n',numel(globalFoldIDs));

%% ===================== CHANNEL UNITS =====================
Finfo = parse_channel_feature_info(fNames,cfg);
ChanUnits = build_channel_units(Finfo);
fprintf('Detected channel units=%d\n',height(ChanUnits));
if isempty(ChanUnits), error('No channel-level features detected.'); end

%% ===================== TASK MATRIX =====================
[idxTrial,yTrial,subjTrial] = get_task_trials(subjAll,condAll,runAll,phaseAll,cfg);
fprintf('Task trial rows=%d\n',sum(idxTrial));
if sum(idxTrial) < 20, error('Too few task trials.'); end

Xraw = double(table2array(T(idxTrial,fNames)));
Xraw(~isfinite(Xraw)) = NaN;
XunitTrial = build_unit_matrix(Xraw,ChanUnits,cfg.unitFeatureAggregate);
[Xtask,yTask,subjTask] = build_subject_condition_samples(XunitTrial,yTrial,subjTrial,cfg);
Folds = make_global_folds_for_this_task(subjTask,GlobalFoldTable,globalFoldIDs);

fprintf('Task samples=%d | subjects=%d | classes=%s | folds=%d | units=%d\n', ...
    size(Xtask,1),numel(unique(subjTask)),strjoin(categories(removecats(yTask))',','),numel(Folds),size(Xtask,2));
if isempty(Folds), error('No valid folds.'); end

%% ===================== OBSERVED FAMILY SEARCH =====================
Observed = run_primary_family_search(Xtask,yTask,Folds,cfg);
Observed = sortrows(Observed,{'BalAcc','AUC','Accuracy'},{'descend','descend','descend'});
writetable(Observed, fullfile(cfg.outDir,'STEP107B_observed_primaryFamily_results.csv'));

fprintf('\nTop observed primary-family results:\n');
disp(Observed(1:min(20,height(Observed)),:));

%% ===================== MAX-PERMUTATION =====================
maxRows = cell(cfg.nPerm,1);

for pi = 1:cfg.nPerm
    rng(cfg.rngSeed + pi);
    yPerm = permute_labels_within_subject(yTask,subjTask,cfg.taskType);
    PermResults = run_primary_family_search(Xtask,yPerm,Folds,cfg);

    maxBal = max(PermResults.BalAcc,[],'omitnan');
    maxAcc = max(PermResults.Accuracy,[],'omitnan');
    maxAUC = max(PermResults.AUC,[],'omitnan');

    bestIdx = find(PermResults.BalAcc == maxBal,1,'first');
    if isempty(bestIdx)
        bestClassifier = "NA"; bestK = NaN;
    else
        bestClassifier = string(PermResults.Classifier(bestIdx));
        bestK = PermResults.K(bestIdx);
    end

    maxRows{pi} = table(pi,maxBal,maxAcc,maxAUC,bestClassifier,bestK, ...
        'VariableNames',{'Permutation','MaxBalAcc','MaxAccuracy','MaxAUC','BestNullClassifier','BestNullK'});

    if mod(pi,10)==0 || pi==cfg.nPerm
        MaxNullPartial = vertcat(maxRows{1:pi});
        writetable(MaxNullPartial, fullfile(cfg.outDir,'STEP107B_max_null_partial.csv'));
        fprintf('perm %d/%d | maxBalAcc=%.3f | best=%s K=%d\n', pi,cfg.nPerm,maxBal,bestClassifier,bestK);
    end
end

MaxNull = vertcat(maxRows{:});
writetable(MaxNull, fullfile(cfg.outDir,'STEP107B_max_null_by_permutation.csv'));

%% ===================== FAMILY-WISE P-VALUES =====================
Corrected = Observed;
pBal = NaN(height(Observed),1);
pAcc = NaN(height(Observed),1);
pAUC = NaN(height(Observed),1);

for i = 1:height(Observed)
    pBal(i) = (1 + sum(MaxNull.MaxBalAcc >= Observed.BalAcc(i))) / (cfg.nPerm + 1);
    pAcc(i) = (1 + sum(MaxNull.MaxAccuracy >= Observed.Accuracy(i))) / (cfg.nPerm + 1);
    if isfinite(Observed.AUC(i))
        pAUC(i) = (1 + sum(MaxNull.MaxAUC >= Observed.AUC(i))) / (cfg.nPerm + 1);
    end
end

Corrected.p_FWER_primaryFamily_BalAcc = pBal;
Corrected.p_FWER_primaryFamily_Accuracy = pAcc;
Corrected.p_FWER_primaryFamily_AUC = pAUC;
Corrected = sortrows(Corrected,{'p_FWER_primaryFamily_BalAcc','BalAcc','AUC'},{'ascend','descend','descend'});

writetable(Corrected, fullfile(cfg.outDir,'STEP107B_FWER_corrected_primaryFamily_results.csv'));
save(fullfile(cfg.outDir,'STEP107B_primaryFamily_maxperm_results.mat'), 'Observed','MaxNull','Corrected','cfg','-v7.3');

fprintf('\n=== STEP107B finished ===\n');
fprintf('Output folder:\n%s\n',cfg.outDir);
fprintf('\nTop FWER-corrected primary-family results:\n');
disp(Corrected(1:min(20,height(Corrected)),:));
fprintf('Finished: %s\n',datestr(now));
diary off;

%% ========================================================================
%% FUNCTIONS
%% ========================================================================
function Results = run_primary_family_search(X,y,Folds,cfg)
rows = {};
for ic = 1:numel(cfg.classifiers)
    clf = cfg.classifiers{ic};
    for ik = 1:numel(cfg.Klist)
        K = cfg.Klist(ik);
        [acc,bal,auc] = run_one_task_model(X,y,Folds,clf,K,cfg);
        rows{end+1,1} = table("channel","binary",string(cfg.contrastName),cfg.run,string(cfg.phase),string(clf),K,acc,bal,auc, ...
            'VariableNames',{'UnitLevel','TaskType','Contrast','Run','Phase','Classifier','K','Accuracy','BalAcc','AUC'}); %#ok<AGROW>
    end
end
Results = vertcat(rows{:});
end

function [Accuracy,BalAcc,AUC] = run_one_task_model(X,y,Folds,clf,K,cfg)
yt = categorical(); yp = categorical(); scoreAll = [];
for fi = 1:numel(Folds)
    tr = Folds(fi).train; te = Folds(fi).test;
    XtrainAll = X(tr,:); ytrain = removecats(y(tr));
    XtestAll = X(te,:); ytest = removecats(y(te));

    Scores = fisher_score_train_only(XtrainAll,ytrain);
    Scores = Scores(isfinite(Scores.FisherScore),:);
    if isempty(Scores), continue; end
    Scores = sortrows(Scores,'FisherScore','descend');
    Sel = Scores(1:min(K,height(Scores)),:);
    idxU = Sel.UnitIndex;

    Xtr = XtrainAll(:,idxU); Xte = XtestAll(:,idxU);
    [Xtr,Xte] = clean_train_test(Xtr,Xte);
    if cfg.trainStandardize, [Xtr,Xte] = train_z(Xtr,Xte); end

    try
        [pred,score] = classify_predict(Xtr,ytrain,Xte,clf,cfg);
    catch
        continue;
    end
    yt = [yt; ytest]; %#ok<AGROW>
    yp = [yp; pred]; %#ok<AGROW>
    scoreAll = [scoreAll; score(:)]; %#ok<AGROW>
end
if isempty(yt), Accuracy = NaN; BalAcc = NaN; AUC = NaN; return; end
M = compute_metrics(yt,yp,scoreAll,false);
Accuracy = M.Accuracy; BalAcc = M.BalAcc; AUC = M.AUC;
end

function [T,col,fNames] = load_dataset(cfg)
datasetPath = '';
for i = 1:numel(cfg.datasetPaths)
    if exist(cfg.datasetPaths{i},'file'), datasetPath = cfg.datasetPaths{i}; break; end
end
if isempty(datasetPath), error('dataset.mat not found.'); end
fprintf('Loading dataset:\n%s\n',datasetPath);
S = load(datasetPath);
if isfield(S,'DS') && istable(S.DS), T = S.DS;
elseif isfield(S,'T') && istable(S.T), T = S.T;
else, error('No table DS or T found.'); end
names = T.Properties.VariableNames;
col.subj  = findvar(names,{'Subject','subject','SubjectID','subjectID','subj','subjID','Subj','SubjID'});
col.run   = findvar(names,{'runNum','RunNum','run','Run'});
col.phase = findvar(names,{'phase','Phase'});
col.cond  = findvar(names,{'Condition','condition','yCondition','Label','label'});
if isempty(col.subj)||isempty(col.run)||isempty(col.phase)||isempty(col.cond), error('Missing required metadata columns.'); end
fNames = detect_features(T);
end

function v = findvar(names,cands)
v = '';
for i=1:numel(cands)
    idx = strcmp(names,cands{i});
    if any(idx), v = names{find(idx,1)}; return; end
end
end

function f = detect_features(T)
names = T.Properties.VariableNames; isNum = false(1,numel(names));
for j=1:numel(names), x = T.(names{j}); isNum(j) = isnumeric(x) || islogical(x); end
numNames = names(isNum);
meta = {'Subject','subject','SubjectID','subjectID','subj','subjID','Subj','SubjID','Run','run','runNum','RunNum','session','Session','Phase','phase','Condition','condition','yCondition','Label','label','Correct','correct','yCorrect','TrialNum','trialNum','TrialIndex','trialIndex','Trial','trial','PatternID','patternID','StartRow','EndRow','Second10Row','RetrRow','Fold','fold','CVFold'};
isMeta = false(size(numNames));
for j=1:numel(numNames), isMeta(j)=any(strcmpi(numNames{j},meta)); end
f = numNames(~isMeta);
end

function s = cleanstr(x)
if isnumeric(x)||islogical(x), s = string(x);
elseif iscell(x), s = string(x);
elseif iscategorical(x), s = string(x);
elseif isstring(x), s = x;
elseif ischar(x), s = string(cellstr(x));
else, s = string(x); end
s = lower(strtrim(s)); s = regexprep(s,'\s+',''); s = regexprep(s,'[^\w]','');
end

function [FoldTable,foldIDs] = create_global_subject_folds(subjAll)
subjects = unique(cleanstr(subjAll)); subjects = subjects(~ismissing(subjects)); nSubj = numel(subjects);
rows = cell(nSubj*nSubj,1); rr = 0;
for i=1:nSubj
    testSubject = subjects(i);
    for j=1:nSubj
        sid = subjects(j); rr = rr + 1;
        if sid == testSubject, role = "test"; else, role = "train"; end
        rows{rr} = table(i,sid,role,testSubject,'VariableNames',{'Fold','SubjectID','Role','LeftOutSubject'});
    end
end
FoldTable = vertcat(rows{:}); foldIDs = unique(FoldTable.Fold);
end

function Finfo = parse_channel_feature_info(fNames,cfg)
n = numel(fNames); FeatureIndex = (1:n)'; FeatureName = string(fNames(:));
Channel = NaN(n,1); ROI = strings(n,1); SpecificFeature = strings(n,1); FeatureFamily = strings(n,1);
IsConnectivity = false(n,1); IsChannelFeature = false(n,1);
for i=1:n
    fname = lower(char(FeatureName(i)));
    IsConnectivity(i) = is_connectivity_name(fname);
    ch = parse_single_channel(fname); Channel(i) = ch; ROI(i) = channel_to_roi(ch);
    SpecificFeature(i) = infer_specific_feature(fname); FeatureFamily(i) = infer_feature_family(fname);
    if isfinite(ch) && ismember(ch,cfg.channelsToUse) && ~(cfg.excludeConnectivityFeatures && IsConnectivity(i)), IsChannelFeature(i) = true; end
end
Finfo = table(FeatureIndex,FeatureName,Channel,ROI,SpecificFeature,FeatureFamily,IsConnectivity,IsChannelFeature);
Finfo = Finfo(Finfo.IsChannelFeature,:);
end

function tf = is_connectivity_name(fname)
s = string(fname);
tf = contains(s,"rie") || contains(s,"pli") || contains(s,"plv") || contains(s,"coh") || contains(s,"conn") || ~isempty(regexp(fname,'c0?\d{1,2}[_-]c0?\d{1,2}','once')) || ~isempty(regexp(fname,'ch0?\d{1,2}[_-]ch0?\d{1,2}','once'));
end

function ch = parse_single_channel(fname)
ch = NaN;
tok = regexp(fname,'(?:^|[_-])ch0?(\d{1,2})(?:[_-]|$)','tokens');
if isempty(tok), tok = regexp(fname,'(?:^|[_-])c0?(\d{1,2})(?:[_-]|$)','tokens'); end
if isempty(tok), tok = regexp(fname,'(?:ch|c)0?(\d{1,2})','tokens'); end
if ~isempty(tok), ch = str2double(tok{1}{1}); end
if ~(isfinite(ch)&&ch>=1&&ch<=128), ch = NaN; end
end

function sf = infer_specific_feature(fname)
s = string(fname);
if contains(s,"rbp") || contains(s,"relative")
    if contains(s,"delta"), sf = "relative_delta"; return; end
    if contains(s,"theta"), sf = "relative_theta"; return; end
    if contains(s,"alpha"), sf = "relative_alpha"; return; end
    if contains(s,"beta"),  sf = "relative_beta"; return; end
    if contains(s,"gamma"), sf = "relative_gamma"; return; end
    sf = "relative_bandpower"; return;
end
if contains(s,"delta"), sf = "delta"; return; end
if contains(s,"theta"), sf = "theta"; return; end
if contains(s,"alpha"), sf = "alpha"; return; end
if contains(s,"beta"),  sf = "beta"; return; end
if contains(s,"gamma"), sf = "gamma"; return; end
if contains(s,"rms"), sf = "rms"; return; end
if contains(s,"peakabs"), sf = "peak_abs"; return; end
if contains(s,"p2p"), sf = "peak_to_peak"; return; end
if contains(s,"mean"), sf = "mean"; return; end
if contains(s,"std"), sf = "std"; return; end
if contains(s,"var"), sf = "variance"; return; end
if contains(s,"tkeo"), sf = "tkeo"; return; end
if contains(s,"auc"), sf = "auc_abs"; return; end
if contains(s,"linelen") || contains(s,"line"), sf = "line_length"; return; end
if contains(s,"lzc") || contains(s,"lz"), sf = "lzc"; return; end
if contains(s,"entropy"), sf = "entropy"; return; end
if contains(s,"higuchi"), sf = "higuchi"; return; end
if contains(s,"fractal"), sf = "fractal"; return; end
sf = "unknown";
end

function ff = infer_feature_family(fname)
sf = infer_specific_feature(fname);
if any(sf == ["delta","theta","alpha","beta","gamma"]), ff = "bandpower";
elseif startsWith(sf,"relative"), ff = "relative_bandpower";
elseif any(sf == ["rms","peak_abs","peak_to_peak","mean","std","variance","tkeo","auc_abs","line_length"]), ff = "temporal_statistical";
elseif any(sf == ["lzc","entropy","higuchi","fractal"]), ff = "complexity";
else, ff = "unknown"; end
end

function roi = channel_to_roi(ch)
if ~isfinite(ch), roi = "unknown"; return; end
if ismember(ch,1:10), roi = "frontal";
elseif ismember(ch,11:20), roi = "frontocentral";
elseif ismember(ch,21:34), roi = "central";
elseif ismember(ch,35:48), roi = "centroparietal";
elseif ismember(ch,49:58), roi = "parietooccipital";
elseif ismember(ch,59:64), roi = "occipital";
else, roi = "unknown"; end
end

function Units = build_channel_units(Finfo)
[G,ch,roi,sf,ff] = findgroups(Finfo.Channel,Finfo.ROI,Finfo.SpecificFeature,Finfo.FeatureFamily);
FeatureIdx = splitapply(@(x){x},Finfo.FeatureIndex,G);
UnitID = "ch" + string(ch) + "_" + string(sf);
Units = table(UnitID,ch,string(roi),string(sf),string(ff),FeatureIdx,'VariableNames',{'UnitID','Channel','ROI','SpecificFeature','FeatureFamily','FeatureIdx'});
end

function [idx,y,subj] = get_task_trials(subjAll,condAll,runAll,phaseAll,cfg)
base = runAll==cfg.run & phaseAll==cleanstr({cfg.phase});
switch string(cfg.contrastName)
    case "color_vs_orientation", A=cfg.labelMap.color; B=cfg.labelMap.orientation; Aname="color"; Bname="orientation";
    case "color_vs_conjunction", A=cfg.labelMap.color; B=cfg.labelMap.conjunction; Aname="color"; Bname="conjunction";
    case "orientation_vs_conjunction", A=cfg.labelMap.orientation; B=cfg.labelMap.conjunction; Aname="orientation"; Bname="conjunction";
    otherwise, error('Unsupported binary contrast: %s',cfg.contrastName);
end
idxA = base & ismember(condAll,cleanstr(A)); idxB = base & ismember(condAll,cleanstr(B)); idx = idxA | idxB;
y = strings(sum(idx),1); y(idxA(idx)) = Aname; y(idxB(idx)) = Bname; y = categorical(y);
subj = subjAll(idx);
end

function XU = build_unit_matrix(Xraw,Units,agg)
XU = NaN(size(Xraw,1),height(Units));
for u=1:height(Units)
    Xi = Xraw(:,Units.FeatureIdx{u});
    if strcmpi(agg,'mean'), XU(:,u) = mean(Xi,2,'omitnan'); else, XU(:,u) = median(Xi,2,'omitnan'); end
end
end

function [Xs,ys,subjS] = build_subject_condition_samples(Xtrial,yTrial,subjTrial,cfg)
subjTrial = cleanstr(subjTrial); subs = unique(subjTrial); subs = subs(~ismissing(subs));
cats = categories(removecats(yTrial)); rows = {}; ys = categorical(); subjS = strings(0,1);
for s=1:numel(subs)
    ixS = subjTrial==subs(s);
    for c=1:numel(cats)
        cls = cats{c}; ix = ixS & yTrial==cls;
        if ~any(ix), continue; end
        Xi = Xtrial(ix,:); Xi(~isfinite(Xi)) = NaN;
        if strcmpi(cfg.subjectConditionAggregate,'mean'), row = mean(Xi,1,'omitnan'); else, row = median(Xi,1,'omitnan'); end
        rows{end+1,1}=row; ys=[ys; categorical(string(cls))]; subjS(end+1,1)=subs(s); %#ok<AGROW>
    end
end
Xs = vertcat(rows{:}); ys = removecats(ys);
end

function Folds = make_global_folds_for_this_task(subjTask,GlobalFoldTable,globalFoldIDs)
subjTask = cleanstr(subjTask); Folds = struct([]); k = 0;
for i=1:numel(globalFoldIDs)
    rows = GlobalFoldTable(GlobalFoldTable.Fold==globalFoldIDs(i),:);
    testSubj = string(rows.SubjectID(rows.Role=="test"));
    if isempty(testSubj), continue; end
    testSubj = cleanstr(testSubj(1)); te = subjTask==testSubj; tr = subjTask~=testSubj;
    if ~any(te) || ~any(tr), continue; end
    k = k + 1; Folds(k).subject = testSubj; Folds(k).foldID = globalFoldIDs(i); Folds(k).train = tr; Folds(k).test = te;
end
end

function yPerm = permute_labels_within_subject(y,subj,taskType)
yPerm = y; subj = cleanstr(subj); subs = unique(subj); subs = subs(~ismissing(subs));
for s=1:numel(subs)
    ix = find(subj==subs(s));
    if strcmpi(taskType,'binary')
        if numel(ix)==2 && rand < 0.5, yPerm(ix)=flipud(yPerm(ix)); else, yPerm(ix)=yPerm(ix(randperm(numel(ix)))); end
    else
        yPerm(ix)=yPerm(ix(randperm(numel(ix))));
    end
end
yPerm = removecats(categorical(yPerm));
end

function Scores = fisher_score_train_only(X,y)
X = double(X); X(~isfinite(X)) = NaN; y = removecats(categorical(y)); classes = categories(y);
nU = size(X,2); FS = NaN(nU,1);
for u=1:nU
    x = X(:,u); ok = isfinite(x);
    if sum(ok)<3, continue; end
    x = x(ok); yy = y(ok);
    if numel(categories(removecats(yy)))<2, continue; end
    muAll = mean(x,'omitnan'); between = 0; within = 0;
    for c=1:numel(classes)
        xc = x(yy==classes{c}); xc = xc(isfinite(xc));
        if isempty(xc), continue; end
        nc = numel(xc); muc = mean(xc,'omitnan'); vc = var(xc,0,'omitnan'); if ~isfinite(vc), vc = 0; end
        between = between + nc*(muc-muAll)^2; within = within + nc*vc;
    end
    if within>eps && isfinite(within), FS(u)=between/within; end
end
Scores = table((1:nU)',FS,'VariableNames',{'UnitIndex','FisherScore'});
end

function [Xtr,Xte] = clean_train_test(Xtr,Xte)
Xtr=double(Xtr); Xte=double(Xte); Xtr(~isfinite(Xtr))=NaN; Xte(~isfinite(Xte))=NaN;
med = median(Xtr,1,'omitnan'); med(~isfinite(med))=0;
for j=1:size(Xtr,2), Xtr(isnan(Xtr(:,j)),j)=med(j); Xte(isnan(Xte(:,j)),j)=med(j); end
sd = std(Xtr,0,1,'omitnan'); keep = isfinite(sd) & sd>eps;
Xtr = Xtr(:,keep); Xte = Xte(:,keep);
end

function [Xtr,Xte] = train_z(Xtr,Xte)
mu = mean(Xtr,1,'omitnan'); sd = std(Xtr,0,1,'omitnan'); mu(~isfinite(mu))=0; sd(~isfinite(sd)|sd<=eps)=1;
Xtr = (Xtr-mu)./sd; Xte = (Xte-mu)./sd;
end

function [pred,scoreOut] = classify_predict(Xtr,ytr,Xte,clf,cfg)
ytr = removecats(categorical(ytr)); cls = categories(ytr);
switch lower(string(clf))
    case "lda"
        M = fitcdiscr(Xtr,ytr,'DiscrimType','linear','ClassNames',cls); [pred,score] = predict(M,Xte);
    case "svm"
        M = fitcsvm(Xtr,ytr,'KernelFunction','linear','Standardize',false,'ClassNames',cls); [pred,score] = predict(M,Xte);
    case "randomforest"
        M = TreeBagger(cfg.nTrees,Xtr,ytr,'Method','classification'); [pc,score] = predict(M,Xte); pred = categorical(string(pc),cls);
    case "decisiontree"
        M = fitctree(Xtr,ytr,'ClassNames',cls); [pred,score] = predict(M,Xte);
    otherwise
        error('Unknown classifier: %s', clf);
end
if exist('score','var') && ~isempty(score) && size(score,2)>=2, scoreOut=score(:,2); else, scoreOut=double(pred==cls{min(2,numel(cls))}); end
scoreOut=double(scoreOut(:));
end

function M = compute_metrics(ytrue,ypred,score,isThree)
ytrue = removecats(categorical(ytrue)); ypred = categorical(ypred,categories(ytrue));
M = struct(); M.Accuracy = mean(ytrue==ypred); cats = categories(ytrue); recalls = NaN(numel(cats),1);
for c=1:numel(cats), recalls(c)=sum(ytrue==cats{c} & ypred==cats{c})/max(sum(ytrue==cats{c}),1); end
M.BalAcc = mean(recalls,'omitnan');
if isThree, M.AUC = NaN; else, try, [~,~,~,M.AUC] = perfcurve(ytrue==cats{2},score,true); catch, M.AUC = NaN; end; end
end

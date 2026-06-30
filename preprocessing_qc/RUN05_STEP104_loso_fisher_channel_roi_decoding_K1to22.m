%% RUN05_STEP104_loso_fisher_channel_roi_decoding_K1to22.m

clear; clc;
cfgBase = get_project_config();

%% ---------------- CONFIG ----------------
cfg = struct();
cfg.ROOT = fullfile(cfgBase.outputRoot, 'new_article_outputs');

cfg.datasetPaths = { ...
    fullfile(cfg.ROOT,'_wm_ml','dataset.mat')
    fullfile(cfg.ROOT,'Subjects','_wm_ml','dataset.mat')
    fullfile(cfg.ROOT,'_wm_dataset','dataset.mat')
    fullfile(cfg.ROOT,'dataset.mat')
};

cfg.Klist = 1:22;
cfg.classifiers = {'lda','svm','randomforest','decisiontree'};
cfg.nTrees = 100;

cfg.runAllTasks = true;
cfg.runList = [1 2 3];
cfg.phaseList = {'stim','maint','retr'};
cfg.binaryContrasts = {'color_vs_orientation','color_vs_conjunction','orientation_vs_conjunction'};
cfg.doBinary = true;
cfg.doThreeClass = true;

cfg.labelMap.color       = {'1','color'};
cfg.labelMap.orientation = {'2','orientation'};
cfg.labelMap.conjunction = {'3','conjunction'};

cfg.channelsToUse = 1:64;
cfg.excludeConnectivityFeatures = true;

cfg.unitFeatureAggregate = 'median';
cfg.subjectConditionAggregate = 'median';
cfg.trainStandardize = true;

cfg.outDir = fullfile(cfg.ROOT,'_wm_STEP104_FisherScore_CHANNEL_ROI_1to22_4classifiers_noLeakage');
if ~exist(cfg.outDir,'dir'), mkdir(cfg.outDir); end
diary(fullfile(cfg.outDir,'STEP104_log.txt'));

fprintf('\n=== STEP104 Fisher Score CHANNEL/ROI features 1..22 + 4 classifiers ===\n');
fprintf('Started: %s\n', datestr(now));

%% ---------------- LOAD DATA ----------------
[T,col,fNames] = load_dataset(cfg);

subjAll  = cleanstr(T.(col.subj));
condAll  = cleanstr(T.(col.cond));
runAll   = double(T.(col.run));
phaseAll = cleanstr(T.(col.phase));

fprintf('Dataset rows=%d | numeric candidate features=%d | subjects=%d\n', ...
    height(T), numel(fNames), numel(unique(subjAll)));

%% ---------------- GLOBAL SUBJECT FOLDS ----------------
[GlobalFoldTable, globalFoldIDs] = create_global_subject_folds(subjAll);
writetable(GlobalFoldTable, fullfile(cfg.outDir,'STEP104_internal_global_subject_LOO_folds.csv'));
fprintf('Global Subject-LOO folds=%d\n', numel(globalFoldIDs));

%% ---------------- CHANNEL / ROI UNITS ----------------
Finfo = parse_channel_feature_info(fNames,cfg);
ChanUnits = build_channel_units(Finfo);
ROIUnits  = build_roi_units(Finfo);

fprintf('Raw channel-level columns=%d\n',height(Finfo));
fprintf('Channel x SpecificFeature units=%d\n',height(ChanUnits));
fprintf('ROI x SpecificFeature units=%d\n',height(ROIUnits));

if isempty(ChanUnits)
    error('No channel-level features detected. Check feature name patterns.');
end

%% ---------------- TASK LIST ----------------
TaskList = build_task_list(cfg);
fprintf('Tasks=%d\n',height(TaskList));

%% ---------------- MAIN ----------------
allResults = {};
allPredictions = {};
allSelected = {};
allStatus = {};

unitSets = struct();
unitSets(1).Level = "channel";
unitSets(1).Units = ChanUnits;
unitSets(2).Level = "roi";
unitSets(2).Units = ROIUnits;

for ui = 1:numel(unitSets)
    unitLevel = unitSets(ui).Level;
    UnitsThis = unitSets(ui).Units;
    if isempty(UnitsThis), continue; end

    fprintf('\n######## UNIT LEVEL: %s | units=%d ########\n',unitLevel,height(UnitsThis));

    for ti = 1:height(TaskList)
        task = TaskList(ti,:);
        fprintf('\nTask %d/%d | %s | %s | run=%d | phase=%s | unit=%s\n', ...
            ti,height(TaskList),task.TaskType,task.Contrast,task.Run,task.Phase,unitLevel);

        try
            cfgTask = cfg;
            cfgTask.taskType = char(string(task.TaskType));
            cfgTask.contrastName = char(string(task.Contrast));
            cfgTask.run = task.Run;
            cfgTask.phase = char(string(task.Phase));

            [idxTrial,yTrial,subjTrial] = get_task_trials(subjAll,condAll,runAll,phaseAll,cfgTask);
            if sum(idxTrial) < 20
                allStatus{end+1,1} = make_status(unitLevel,task,NaN,"too_few_trials"); %#ok<AGROW>
                continue;
            end

            Xraw = double(table2array(T(idxTrial,fNames)));
            Xraw(~isfinite(Xraw)) = NaN;

            XunitTrial = build_unit_matrix(Xraw,UnitsThis,cfg.unitFeatureAggregate);
            [Xtask,yTask,subjTask] = build_subject_condition_samples(XunitTrial,yTrial,subjTrial,cfgTask);

            Folds = make_global_folds_for_this_task(subjTask,GlobalFoldTable,globalFoldIDs);
            if isempty(Folds)
                allStatus{end+1,1} = make_status(unitLevel,task,NaN,"no_valid_folds"); %#ok<AGROW>
                continue;
            end

            taskResults = {};
            taskPreds = {};
            taskSel = {};

            for ic = 1:numel(cfg.classifiers)
                clf = cfg.classifiers{ic};
                for ik = 1:numel(cfg.Klist)
                    K = cfg.Klist(ik);

                    [R,P,S] = run_one_model_K(Xtask,yTask,subjTask,Folds,UnitsThis,unitLevel,task,clf,K,cfgTask);

                    taskResults{end+1,1} = R; %#ok<AGROW>
                    taskPreds{end+1,1} = P; %#ok<AGROW>
                    taskSel{end+1,1} = S; %#ok<AGROW>

                    allResults{end+1,1} = R; %#ok<AGROW>
                    allPredictions{end+1,1} = P; %#ok<AGROW>
                    allSelected{end+1,1} = S; %#ok<AGROW>
                end
                fprintf('  done classifier=%s\n',clf);
            end

            TaskResults = sortrows(vertcat(taskResults{:}), {'BalAcc','AUC','Accuracy'}, {'descend','descend','descend'});
            TaskPreds = vertcat(taskPreds{:});
            TaskSel = vertcat(taskSel{:});

            tag = sprintf('%s_%s_%s_run%d_%s',unitLevel,task.TaskType,task.Contrast,task.Run,task.Phase);
            tag = regexprep(tag,'[^\w]','_');
            outSub = fullfile(cfg.outDir,tag);
            if ~exist(outSub,'dir'), mkdir(outSub); end

            writetable(TaskResults, fullfile(outSub,sprintf('STEP104_%s_results.csv',tag)));
            writetable(TaskPreds, fullfile(outSub,sprintf('STEP104_%s_predictions.csv',tag)));
            writetable(TaskSel, fullfile(outSub,sprintf('STEP104_%s_selected_features.csv',tag)));
            plot_curves(TaskResults,outSub,tag);

            allStatus{end+1,1} = make_status(unitLevel,task,height(TaskResults),"ok"); %#ok<AGROW>
            writetable(vertcat(allResults{:}), fullfile(cfg.outDir,'STEP104_partial_results.csv'));

        catch ME
            warning('Task failed: %s',ME.message);
            allStatus{end+1,1} = make_status(unitLevel,task,NaN,"failed_" + string(ME.message)); %#ok<AGROW>
        end
    end
end

%% ---------------- SAVE ----------------
if isempty(allResults), AllResults=table(); else, AllResults=sortrows(vertcat(allResults{:}), {'BalAcc','AUC','Accuracy'}, {'descend','descend','descend'}); end
if isempty(allPredictions), AllPredictions=table(); else, AllPredictions=vertcat(allPredictions{:}); end
if isempty(allSelected), AllSelected=table(); else, AllSelected=vertcat(allSelected{:}); end
Status = vertcat(allStatus{:});

writetable(AllResults, fullfile(cfg.outDir,'STEP104_FisherScore_CHANNEL_ROI_featureCountCurve_results.csv'));
writetable(AllPredictions, fullfile(cfg.outDir,'STEP104_FisherScore_CHANNEL_ROI_featureCountCurve_predictions.csv'));
writetable(AllSelected, fullfile(cfg.outDir,'STEP104_FisherScore_CHANNEL_ROI_selected_features_by_fold.csv'));
writetable(Status, fullfile(cfg.outDir,'STEP104_status.csv'));

save(fullfile(cfg.outDir,'STEP104_FisherScore_CHANNEL_ROI_results.mat'), ...
    'AllResults','AllPredictions','AllSelected','Status','TaskList','cfg','ChanUnits','ROIUnits','Finfo','-v7.3');

fprintf('\n=== STEP104 finished ===\n');
fprintf('Output folder:\n%s\n',cfg.outDir);
disp(AllResults(1:min(30,height(AllResults)),:));
diary off;

%% ================= FUNCTIONS =================

function TaskList = build_task_list(cfg)
TaskList = table();
if cfg.doBinary
    for r=cfg.runList
        for p=1:numel(cfg.phaseList)
            for c=1:numel(cfg.binaryContrasts)
                TaskList = [TaskList; table("binary",string(cfg.binaryContrasts{c}),r,string(cfg.phaseList{p}), ...
                    'VariableNames',{'TaskType','Contrast','Run','Phase'})]; %#ok<AGROW>
            end
        end
    end
end
if cfg.doThreeClass
    for r=cfg.runList
        for p=1:numel(cfg.phaseList)
            TaskList = [TaskList; table("threeclass","color_vs_orientation_vs_conjunction",r,string(cfg.phaseList{p}), ...
                'VariableNames',{'TaskType','Contrast','Run','Phase'})]; %#ok<AGROW>
        end
    end
end
end

function S = make_status(unitLevel,task,nRows,status)
S = table(string(unitLevel),string(task.TaskType),string(task.Contrast),task.Run,string(task.Phase),nRows,string(status), ...
    'VariableNames',{'UnitLevel','TaskType','Contrast','Run','Phase','N_rows','Status'});
end

function [T,col,fNames] = load_dataset(cfg)
datasetPath = '';
for i=1:numel(cfg.datasetPaths)
    if exist(cfg.datasetPaths{i},'file'), datasetPath=cfg.datasetPaths{i}; break; end
end
if isempty(datasetPath), error('dataset.mat not found.'); end
S=load(datasetPath);
if isfield(S,'DS') && istable(S.DS), T=S.DS;
elseif isfield(S,'T') && istable(S.T), T=S.T;
else, error('No DS/T table found.'); end

names=T.Properties.VariableNames;
col.subj=findvar(names,{'Subject','subject','SubjectID','subjectID','subj','subjID','Subj','SubjID'});
col.run=findvar(names,{'runNum','RunNum','run','Run'});
col.phase=findvar(names,{'phase','Phase'});
col.cond=findvar(names,{'Condition','condition','yCondition','Label','label'});
fNames=detect_features(T);
end

function v=findvar(names,cands)
v='';
for i=1:numel(cands)
    idx=strcmp(names,cands{i});
    if any(idx), v=names{find(idx,1)}; return; end
end
end

function f=detect_features(T)
names=T.Properties.VariableNames;
isNum=false(1,numel(names));
for j=1:numel(names)
    x=T.(names{j});
    isNum(j)=isnumeric(x)||islogical(x);
end
numNames=names(isNum);
meta={'Subject','subject','SubjectID','subjectID','subj','subjID','Subj','SubjID','Run','run','runNum','RunNum','session','Session','Phase','phase','Condition','condition','yCondition','Label','label','Correct','correct','yCorrect','TrialNum','trialNum','TrialIndex','trialIndex','Trial','trial','PatternID','patternID','StartRow','EndRow','Second10Row','RetrRow','Fold','fold','CVFold'};
isMeta=false(size(numNames));
for j=1:numel(numNames), isMeta(j)=any(strcmpi(numNames{j},meta)); end
f=numNames(~isMeta);
end

function s=cleanstr(x)
if isnumeric(x)||islogical(x), s=string(x);
elseif iscell(x), s=string(x);
elseif iscategorical(x), s=string(x);
elseif isstring(x), s=x;
elseif ischar(x), s=string(cellstr(x));
else, s=string(x); end
s=lower(strtrim(s)); s=regexprep(s,'\s+',''); s=regexprep(s,'[^\w]','');
end

function [FoldTable, foldIDs]=create_global_subject_folds(subjAll)
subjects=unique(cleanstr(subjAll)); subjects=subjects(~ismissing(subjects));
rows=cell(numel(subjects)^2,1); rr=0;
for i=1:numel(subjects)
    testSubject=subjects(i);
    for j=1:numel(subjects)
        sid=subjects(j); rr=rr+1;
        if sid==testSubject, role="test"; else, role="train"; end
        rows{rr}=table(i,sid,role,testSubject,'VariableNames',{'Fold','SubjectID','Role','LeftOutSubject'});
    end
end
FoldTable=vertcat(rows{:}); foldIDs=unique(FoldTable.Fold);
end

function Finfo=parse_channel_feature_info(fNames,cfg)
n=numel(fNames);
FeatureIndex=(1:n)';
FeatureName=string(fNames(:));
Channel=NaN(n,1); ROI=strings(n,1); SpecificFeature=strings(n,1); FeatureFamily=strings(n,1);
IsConnectivity=false(n,1); IsChannelFeature=false(n,1);

for i=1:n
    fname=lower(char(FeatureName(i)));
    IsConnectivity(i)=is_connectivity_name(fname);
    ch=parse_single_channel(fname);
    Channel(i)=ch;
    ROI(i)=channel_to_roi(ch);
    SpecificFeature(i)=infer_specific_feature(fname);
    FeatureFamily(i)=infer_feature_family(fname);

    if isfinite(ch) && ismember(ch,cfg.channelsToUse) && ~(cfg.excludeConnectivityFeatures && IsConnectivity(i))
        IsChannelFeature(i)=true;
    end
end

Finfo=table(FeatureIndex,FeatureName,Channel,ROI,SpecificFeature,FeatureFamily,IsConnectivity,IsChannelFeature);
Finfo=Finfo(Finfo.IsChannelFeature,:);
end

function tf=is_connectivity_name(fname)
s=string(fname);
tf=contains(s,"rie")||contains(s,"pli")||contains(s,"plv")||contains(s,"coh")||contains(s,"conn")|| ...
   ~isempty(regexp(fname,'c0?\d{1,2}[_-]c0?\d{1,2}','once')) || ...
   ~isempty(regexp(fname,'ch0?\d{1,2}[_-]ch0?\d{1,2}','once'));
end

function ch=parse_single_channel(fname)
ch=NaN;
tok=regexp(fname,'(?:^|[_-])ch0?(\d{1,2})(?:[_-]|$)','tokens');
if isempty(tok), tok=regexp(fname,'(?:^|[_-])c0?(\d{1,2})(?:[_-]|$)','tokens'); end
if isempty(tok), tok=regexp(fname,'(?:ch|c)0?(\d{1,2})','tokens'); end
if ~isempty(tok), ch=str2double(tok{1}{1}); end
if ~(isfinite(ch)&&ch>=1&&ch<=128), ch=NaN; end
end

function sf=infer_specific_feature(fname)
s=string(fname);
if contains(s,"rbp")||contains(s,"relative")
    if contains(s,"delta"), sf="relative_delta"; return; end
    if contains(s,"theta"), sf="relative_theta"; return; end
    if contains(s,"alpha"), sf="relative_alpha"; return; end
    if contains(s,"beta"), sf="relative_beta"; return; end
    if contains(s,"gamma"), sf="relative_gamma"; return; end
    sf="relative_bandpower"; return;
end
if contains(s,"delta"), sf="delta"; return; end
if contains(s,"theta"), sf="theta"; return; end
if contains(s,"alpha"), sf="alpha"; return; end
if contains(s,"beta"), sf="beta"; return; end
if contains(s,"gamma"), sf="gamma"; return; end
if contains(s,"rms"), sf="rms"; return; end
if contains(s,"peakabs"), sf="peak_abs"; return; end
if contains(s,"p2p"), sf="peak_to_peak"; return; end
if contains(s,"mean"), sf="mean"; return; end
if contains(s,"std"), sf="std"; return; end
if contains(s,"var"), sf="variance"; return; end
if contains(s,"tkeo"), sf="tkeo"; return; end
if contains(s,"auc"), sf="auc_abs"; return; end
if contains(s,"linelen")||contains(s,"line"), sf="line_length"; return; end
if contains(s,"lzc")||contains(s,"lz"), sf="lzc"; return; end
if contains(s,"entropy"), sf="entropy"; return; end
if contains(s,"higuchi"), sf="higuchi"; return; end
if contains(s,"fractal"), sf="fractal"; return; end
sf="unknown";
end

function ff=infer_feature_family(fname)
sf=infer_specific_feature(fname);
if any(sf==["delta","theta","alpha","beta","gamma"]), ff="bandpower";
elseif startsWith(sf,"relative"), ff="relative_bandpower";
elseif any(sf==["rms","peak_abs","peak_to_peak","mean","std","variance","tkeo","auc_abs","line_length"]), ff="temporal_statistical";
elseif any(sf==["lzc","entropy","higuchi","fractal"]), ff="complexity";
else, ff="unknown"; end
end

function roi=channel_to_roi(ch)
if ~isfinite(ch), roi="unknown"; return; end
if ismember(ch,1:10), roi="frontal";
elseif ismember(ch,11:20), roi="frontocentral";
elseif ismember(ch,21:34), roi="central";
elseif ismember(ch,35:48), roi="centroparietal";
elseif ismember(ch,49:58), roi="parietooccipital";
elseif ismember(ch,59:64), roi="occipital";
else, roi="unknown"; end
end

function Units=build_channel_units(Finfo)
[G,ch,roi,sf,ff]=findgroups(Finfo.Channel,Finfo.ROI,Finfo.SpecificFeature,Finfo.FeatureFamily);
FeatureIdx=splitapply(@(x){x},Finfo.FeatureIndex,G);
UnitID="ch"+string(ch)+"_"+string(sf);
Units=table(UnitID,ch,string(roi),string(sf),string(ff),FeatureIdx,'VariableNames',{'UnitID','Channel','ROI','SpecificFeature','FeatureFamily','FeatureIdx'});
end

function Units=build_roi_units(Finfo)
[G,roi,sf,ff]=findgroups(Finfo.ROI,Finfo.SpecificFeature,Finfo.FeatureFamily);
FeatureIdx=splitapply(@(x){x},Finfo.FeatureIndex,G);
UnitID=string(roi)+"_"+string(sf);
Units=table(UnitID,string(roi),string(sf),string(ff),FeatureIdx,'VariableNames',{'UnitID','ROI','SpecificFeature','FeatureFamily','FeatureIdx'});
end

function [idx,y,subj]=get_task_trials(subjAll,condAll,runAll,phaseAll,cfg)
base=runAll==cfg.run & phaseAll==cleanstr({cfg.phase});
if strcmpi(cfg.taskType,'threeclass')
    idxC=base & ismember(condAll,cleanstr(cfg.labelMap.color));
    idxO=base & ismember(condAll,cleanstr(cfg.labelMap.orientation));
    idxJ=base & ismember(condAll,cleanstr(cfg.labelMap.conjunction));
    idx=idxC|idxO|idxJ;
    y=strings(sum(idx),1); y(idxC(idx))="color"; y(idxO(idx))="orientation"; y(idxJ(idx))="conjunction"; y=categorical(y);
else
    switch string(cfg.contrastName)
        case "color_vs_orientation", A=cfg.labelMap.color; B=cfg.labelMap.orientation; Aname="color"; Bname="orientation";
        case "color_vs_conjunction", A=cfg.labelMap.color; B=cfg.labelMap.conjunction; Aname="color"; Bname="conjunction";
        case "orientation_vs_conjunction", A=cfg.labelMap.orientation; B=cfg.labelMap.conjunction; Aname="orientation"; Bname="conjunction";
        otherwise, error('Unsupported contrast');
    end
    idxA=base & ismember(condAll,cleanstr(A));
    idxB=base & ismember(condAll,cleanstr(B));
    idx=idxA|idxB;
    y=strings(sum(idx),1); y(idxA(idx))=Aname; y(idxB(idx))=Bname; y=categorical(y);
end
subj=subjAll(idx);
end

function XU=build_unit_matrix(Xraw,Units,agg)
XU=NaN(size(Xraw,1),height(Units));
for u=1:height(Units)
    Xi=Xraw(:,Units.FeatureIdx{u});
    if strcmpi(agg,'mean'), XU(:,u)=mean(Xi,2,'omitnan');
    else, XU(:,u)=median(Xi,2,'omitnan'); end
end
end

function [Xs,ys,subjS]=build_subject_condition_samples(Xtrial,yTrial,subjTrial,cfg)
subjTrial=cleanstr(subjTrial); subs=unique(subjTrial); subs=subs(~ismissing(subs));
cats=categories(removecats(yTrial));
rows={}; ys=categorical(); subjS=strings(0,1);
for s=1:numel(subs)
    ixS=subjTrial==subs(s);
    for c=1:numel(cats)
        cls=cats{c}; ix=ixS & yTrial==cls;
        if ~any(ix), continue; end
        Xi=Xtrial(ix,:); Xi(~isfinite(Xi))=NaN;
        if strcmpi(cfg.subjectConditionAggregate,'mean'), row=mean(Xi,1,'omitnan');
        else, row=median(Xi,1,'omitnan'); end
        rows{end+1,1}=row; ys=[ys; categorical(string(cls))]; subjS(end+1,1)=subs(s); %#ok<AGROW>
    end
end
Xs=vertcat(rows{:}); ys=removecats(ys);
end

function Folds=make_global_folds_for_this_task(subjTask,GlobalFoldTable,globalFoldIDs)
subjTask=cleanstr(subjTask); Folds=struct([]); k=0;
for i=1:numel(globalFoldIDs)
    rows=GlobalFoldTable(GlobalFoldTable.Fold==globalFoldIDs(i),:);
    testSubj=string(rows.SubjectID(rows.Role=="test"));
    if isempty(testSubj), continue; end
    testSubj=cleanstr(testSubj(1));
    te=subjTask==testSubj; tr=subjTask~=testSubj;
    if ~any(te)||~any(tr), continue; end
    k=k+1; Folds(k).subject=testSubj; Folds(k).foldID=globalFoldIDs(i); Folds(k).train=tr; Folds(k).test=te;
end
end

function [R,P,S]=run_one_model_K(X,y,subj,Folds,Units,unitLevel,task,clf,K,cfg)
yt=categorical(); yp=categorical(); scoreAll=[]; pRows={}; sRows={};
for fi=1:numel(Folds)
    tr=Folds(fi).train; te=Folds(fi).test;
    XtrainAll=X(tr,:); ytrain=removecats(y(tr));
    XtestAll=X(te,:); ytest=removecats(y(te));

    Scores=fisher_score_train_only(XtrainAll,ytrain,Units);
    Scores=Scores(isfinite(Scores.FisherScore),:);
    if isempty(Scores), continue; end
    Scores=sortrows(Scores,'FisherScore','descend');
    Sel=Scores(1:min(K,height(Scores)),:);
    idxU=Sel.UnitIndex;

    Xtr=XtrainAll(:,idxU); Xte=XtestAll(:,idxU);
    [Xtr,Xte]=clean_train_test(Xtr,Xte);
    if cfg.trainStandardize, [Xtr,Xte]=train_z(Xtr,Xte); end

    [pred,score]=classify_predict(Xtr,ytrain,Xte,clf,cfg);

    yt=[yt;ytest]; yp=[yp;pred]; scoreAll=[scoreAll;score(:)]; %#ok<AGROW>

    pRows{end+1,1}=table(repmat(string(unitLevel),sum(te),1),repmat(string(task.TaskType),sum(te),1),repmat(string(task.Contrast),sum(te),1),repmat(task.Run,sum(te),1),repmat(string(task.Phase),sum(te),1),repmat(string(clf),sum(te),1),repmat(K,sum(te),1),repmat(fi,sum(te),1),repmat(string(Folds(fi).subject),sum(te),1),string(ytest),string(pred),double(score(:)), ...
        'VariableNames',{'UnitLevel','TaskType','Contrast','Run','Phase','Classifier','N_features','Fold','LeftOutSubject','YTrue','YPred','Score'}); %#ok<AGROW>

    Sel.UnitLevel=repmat(string(unitLevel),height(Sel),1);
    Sel.TaskType=repmat(string(task.TaskType),height(Sel),1);
    Sel.Contrast=repmat(string(task.Contrast),height(Sel),1);
    Sel.Run=repmat(task.Run,height(Sel),1);
    Sel.Phase=repmat(string(task.Phase),height(Sel),1);
    Sel.Classifier=repmat(string(clf),height(Sel),1);
    Sel.N_features=repmat(K,height(Sel),1);
    Sel.Fold=repmat(fi,height(Sel),1);
    Sel.LeftOutSubject=repmat(string(Folds(fi).subject),height(Sel),1);
    Sel=movevars(Sel,{'UnitLevel','TaskType','Contrast','Run','Phase','Classifier','N_features','Fold','LeftOutSubject'},'Before',1);
    sRows{end+1,1}=Sel; %#ok<AGROW>
end
if isempty(yt), error('No predictions generated.'); end
isThree=string(task.TaskType)=="threeclass";
M=compute_metrics(yt,yp,scoreAll,isThree);
R=table(string(unitLevel),string(task.TaskType),string(task.Contrast),task.Run,string(task.Phase),string(clf),K,numel(Folds),numel(yt),M.Accuracy,M.BalAcc,M.AUC, ...
    'VariableNames',{'UnitLevel','TaskType','Contrast','Run','Phase','Classifier','N_features','N_folds','N_predictions','Accuracy','BalAcc','AUC'});
P=vertcat(pRows{:}); S=vertcat(sRows{:});
end

function Scores=fisher_score_train_only(X,y,Units)
X=double(X); X(~isfinite(X))=NaN; y=removecats(categorical(y)); classes=categories(y);
nU=size(X,2); FS=NaN(nU,1);
for u=1:nU
    x=X(:,u); ok=isfinite(x); x=x(ok); yy=y(ok);
    if numel(x)<3 || numel(categories(removecats(yy)))<2, continue; end
    muAll=mean(x,'omitnan'); between=0; within=0;
    for c=1:numel(classes)
        xc=x(yy==classes{c}); xc=xc(isfinite(xc));
        if isempty(xc), continue; end
        nc=numel(xc); muc=mean(xc,'omitnan'); vc=var(xc,0,'omitnan'); if ~isfinite(vc), vc=0; end
        between=between+nc*(muc-muAll)^2; within=within+nc*vc;
    end
    if within>eps && isfinite(within), FS(u)=between/within; end
end
Scores=base_unit_table(Units);
Scores.UnitIndex=(1:nU)';
Scores=movevars(Scores,'UnitIndex','Before',1);
Scores.FisherScore=FS;
end

function B=base_unit_table(Units)
B=table(); B.UnitID=Units.UnitID;
vars=Units.Properties.VariableNames;
if any(strcmp(vars,'Channel')), B.Channel=Units.Channel; end
if any(strcmp(vars,'ROI')), B.ROI=Units.ROI; end
if any(strcmp(vars,'SpecificFeature')), B.SpecificFeature=Units.SpecificFeature; end
if any(strcmp(vars,'FeatureFamily')), B.FeatureFamily=Units.FeatureFamily; end
end

function [Xtr,Xte]=clean_train_test(Xtr,Xte)
Xtr=double(Xtr); Xte=double(Xte); Xtr(~isfinite(Xtr))=NaN; Xte(~isfinite(Xte))=NaN;
med=median(Xtr,1,'omitnan'); med(~isfinite(med))=0;
for j=1:size(Xtr,2), Xtr(isnan(Xtr(:,j)),j)=med(j); Xte(isnan(Xte(:,j)),j)=med(j); end
sd=std(Xtr,0,1,'omitnan'); keep=isfinite(sd)&sd>eps;
Xtr=Xtr(:,keep); Xte=Xte(:,keep);
end

function [Xtr,Xte]=train_z(Xtr,Xte)
mu=mean(Xtr,1,'omitnan'); sd=std(Xtr,0,1,'omitnan'); mu(~isfinite(mu))=0; sd(~isfinite(sd)|sd<=eps)=1;
Xtr=(Xtr-mu)./sd; Xte=(Xte-mu)./sd;
end

function [pred,scoreOut]=classify_predict(Xtr,ytr,Xte,clf,cfg)
ytr=removecats(categorical(ytr)); cls=categories(ytr); isMulti=numel(cls)>2;
switch lower(string(clf))
    case "lda"
        M=fitcdiscr(Xtr,ytr,'DiscrimType','linear','ClassNames',cls); [pred,score]=predict(M,Xte);
    case "svm"
        if isMulti, t=templateSVM('KernelFunction','linear','Standardize',false); M=fitcecoc(Xtr,ytr,'Learners',t,'ClassNames',cls);
        else, M=fitcsvm(Xtr,ytr,'KernelFunction','linear','Standardize',false,'ClassNames',cls); end
        [pred,score]=predict(M,Xte);
    case "randomforest"
        M=TreeBagger(cfg.nTrees,Xtr,ytr,'Method','classification'); [pc,score]=predict(M,Xte); pred=categorical(string(pc),cls);
    case "decisiontree"
        M=fitctree(Xtr,ytr,'ClassNames',cls); [pred,score]=predict(M,Xte);
    otherwise
        error('Unknown classifier');
end
if exist('score','var') && ~isempty(score) && size(score,2)>=2
    if isMulti, scoreOut=max(score,[],2); else, scoreOut=score(:,2); end
else
    scoreOut=double(pred==cls{min(2,numel(cls))});
end
scoreOut=double(scoreOut(:));
end

function M=compute_metrics(ytrue,ypred,score,isThree)
ytrue=removecats(categorical(ytrue)); ypred=categorical(ypred,categories(ytrue));
M=struct(); M.Accuracy=mean(ytrue==ypred);
cats=categories(ytrue); recalls=NaN(numel(cats),1);
for c=1:numel(cats)
    recalls(c)=sum(ytrue==cats{c} & ypred==cats{c})/max(sum(ytrue==cats{c}),1);
end
M.BalAcc=mean(recalls,'omitnan');
if isThree, M.AUC=NaN;
else
    try, [~,~,~,M.AUC]=perfcurve(ytrue==cats{2},score,true);
    catch, M.AUC=NaN; end
end
end

function plot_curves(R,outSub,tag)
clfs=unique(string(R.Classifier));
for i=1:numel(clfs)
    Rc=R(string(R.Classifier)==clfs(i),:); Rc=sortrows(Rc,'N_features');
    figure('Color','w','Position',[100 100 900 600]);
    plot(Rc.N_features,Rc.Accuracy,'-o','LineWidth',2); hold on;
    plot(Rc.N_features,Rc.BalAcc,'-s','LineWidth',2);
    if any(isfinite(Rc.AUC)), plot(Rc.N_features,Rc.AUC,'-^','LineWidth',2); end
    grid on; xlabel('Number of selected features'); ylabel('Performance');
    title(sprintf('%s | %s | Fisher Score',strrep(tag,'_',' '),clfs(i)),'Interpreter','none');
    legend({'Accuracy','Balanced Accuracy','AUC'},'Location','best'); ylim([0 1]);
    figName=regexprep(sprintf('STEP104_%s_%s_curve',tag,clfs(i)),'[^\w]','_');
    savefig(fullfile(outSub,[figName '.fig']));
    saveas(gcf,fullfile(outSub,[figName '.png']));
    close(gcf);
    writetable(Rc(:,{'N_features','Accuracy','BalAcc','AUC'}),fullfile(outSub,[figName '_values.csv']));
end
end

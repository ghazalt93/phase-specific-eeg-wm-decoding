%% STEP100_channel_ROI_feature_effects_subjectMedian_stats.m

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

cfg.runList = [1 2 3];
cfg.phaseList = {'stim','maint','retr'};
cfg.binaryContrasts = {'color_vs_orientation','color_vs_conjunction','orientation_vs_conjunction'};

cfg.labelMap.color       = {'1','color'};
cfg.labelMap.orientation = {'2','orientation'};
cfg.labelMap.conjunction = {'3','conjunction'};

cfg.alphaFDR = 0.05;
cfg.topN = 50;
cfg.subjectConditionAggregate = 'median';

cfg.channelsToUse = 1:64;

cfg.excludeConnectivityFeatures = true;

cfg.outDir = fullfile(cfg.ROOT,'_wm_STEP100_channel_ROI_feature_effects_subjectMedian_stats');
if ~exist(cfg.outDir,'dir'), mkdir(cfg.outDir); end

diary(fullfile(cfg.outDir,'STEP100_log.txt'));

fprintf('\n=== STEP100 channel/ROI feature-level subject-median statistics ===\n');
fprintf('Started: %s\n', datestr(now));
fprintf('ROOT: %s\n', cfg.ROOT);

%% ===================== LOAD DATA =====================

[T,col,fNames] = local_load_dataset(cfg);

subjAll  = local_cleanstr(T.(col.subj));
condAll  = local_cleanstr(T.(col.cond));
runAll   = double(T.(col.run));
phaseAll = local_cleanstr(T.(col.phase));

fprintf('Dataset rows=%d | numeric features=%d | subjects=%d\n', ...
    height(T), numel(fNames), numel(unique(subjAll)));

%% ===================== PARSE CHANNEL FEATURES =====================

Finfo = local_parse_channel_features(fNames,cfg);
ChanUnits = local_build_channel_units(Finfo);
ROIUnits  = local_build_roi_units(Finfo);

fprintf('Channel-level raw features detected: %d\n', height(Finfo));
fprintf('Channel x specific-feature units: %d\n', height(ChanUnits));
fprintf('ROI x specific-feature units: %d\n', height(ROIUnits));

if isempty(ChanUnits)
    error('No channel-level units detected. Check feature-name patterns.');
end

%% ===================== MAIN LOOP =====================

AllBinChan = {};
AllBinROI  = {};
AllTriChan = {};
AllTriROI  = {};
AllStatus  = {};

for rr = cfg.runList
    for pp = 1:numel(cfg.phaseList)
        phaseName = cfg.phaseList{pp};

        fprintf('\n============================================================\n');
        fprintf('Run %d | Phase %s\n', rr, phaseName);
        fprintf('============================================================\n');

        %% ---------------- BINARY ----------------
        for cc = 1:numel(cfg.binaryContrasts)
            contrastName = cfg.binaryContrasts{cc};
            fprintf('\n--- Binary: %s | run=%d | phase=%s ---\n', contrastName, rr, phaseName);

            try
                task = struct();
                task.contrastName = contrastName;
                task.run = rr;
                task.phase = phaseName;

                [idx,yTrial,subjTrial] = local_binary_rows(subjAll,condAll,runAll,phaseAll,task,cfg);

                if sum(idx) < 10
                    AllStatus{end+1,1} = local_status('binary',contrastName,rr,phaseName,sum(idx),NaN,'too_few_trials'); %#ok<AGROW>
                    continue;
                end

                Xraw = double(table2array(T(idx,fNames)));
                Xraw(~isfinite(Xraw)) = NaN;

                XchanTrial = local_unit_matrix(Xraw,ChanUnits);
                [XchanSC,ySC,subjSC] = local_subject_condition_samples(XchanTrial,yTrial,subjTrial,cfg);
                Rchan = local_binary_paired_stats(XchanSC,ySC,subjSC,ChanUnits,cfg);
                Rchan.TaskType = repmat("binary",height(Rchan),1);
                Rchan.Contrast = repmat(string(contrastName),height(Rchan),1);
                Rchan.Run = repmat(rr,height(Rchan),1);
                Rchan.Phase = repmat(string(phaseName),height(Rchan),1);
                Rchan = movevars(Rchan,{'TaskType','Contrast','Run','Phase'},'Before',1);
                AllBinChan{end+1,1} = Rchan; %#ok<AGROW>

                XroiTrial = local_unit_matrix(Xraw,ROIUnits);
                [XroiSC,yRSC,subjRSC] = local_subject_condition_samples(XroiTrial,yTrial,subjTrial,cfg);
                Rroi = local_binary_paired_stats(XroiSC,yRSC,subjRSC,ROIUnits,cfg);
                Rroi.TaskType = repmat("binary",height(Rroi),1);
                Rroi.Contrast = repmat(string(contrastName),height(Rroi),1);
                Rroi.Run = repmat(rr,height(Rroi),1);
                Rroi.Phase = repmat(string(phaseName),height(Rroi),1);
                Rroi = movevars(Rroi,{'TaskType','Contrast','Run','Phase'},'Before',1);
                AllBinROI{end+1,1} = Rroi; %#ok<AGROW>

                AllStatus{end+1,1} = local_status('binary',contrastName,rr,phaseName,sum(idx),numel(unique(subjSC)),'ok'); %#ok<AGROW>

                fprintf('Done binary: channel units=%d | ROI units=%d | subjects=%d\n', ...
                    height(Rchan),height(Rroi),numel(unique(subjSC)));

            catch ME
                warning('Binary failed: %s', ME.message);
                AllStatus{end+1,1} = local_status('binary',contrastName,rr,phaseName,NaN,NaN,['failed_' ME.message]); %#ok<AGROW>
            end
        end

        %% ---------------- THREE CLASS ----------------
        fprintf('\n--- Three-class | run=%d | phase=%s ---\n', rr, phaseName);

        try
            task = struct(); task.run = rr; task.phase = phaseName;
            [idx,yTrial,subjTrial] = local_three_rows(subjAll,condAll,runAll,phaseAll,task,cfg);

            if sum(idx) < 10
                AllStatus{end+1,1} = local_status('threeclass','color_vs_orientation_vs_conjunction',rr,phaseName,sum(idx),NaN,'too_few_trials'); %#ok<AGROW>
                continue;
            end

            Xraw = double(table2array(T(idx,fNames)));
            Xraw(~isfinite(Xraw)) = NaN;

            XchanTrial = local_unit_matrix(Xraw,ChanUnits);
            [XchanSC,ySC,subjSC] = local_subject_condition_samples(XchanTrial,yTrial,subjTrial,cfg);
            Rchan3 = local_three_friedman_stats(XchanSC,ySC,subjSC,ChanUnits,cfg);
            Rchan3.TaskType = repmat("threeclass",height(Rchan3),1);
            Rchan3.Contrast = repmat("color_vs_orientation_vs_conjunction",height(Rchan3),1);
            Rchan3.Run = repmat(rr,height(Rchan3),1);
            Rchan3.Phase = repmat(string(phaseName),height(Rchan3),1);
            Rchan3 = movevars(Rchan3,{'TaskType','Contrast','Run','Phase'},'Before',1);
            AllTriChan{end+1,1} = Rchan3; %#ok<AGROW>

            XroiTrial = local_unit_matrix(Xraw,ROIUnits);
            [XroiSC,yRSC,subjRSC] = local_subject_condition_samples(XroiTrial,yTrial,subjTrial,cfg);
            Rroi3 = local_three_friedman_stats(XroiSC,yRSC,subjRSC,ROIUnits,cfg);
            Rroi3.TaskType = repmat("threeclass",height(Rroi3),1);
            Rroi3.Contrast = repmat("color_vs_orientation_vs_conjunction",height(Rroi3),1);
            Rroi3.Run = repmat(rr,height(Rroi3),1);
            Rroi3.Phase = repmat(string(phaseName),height(Rroi3),1);
            Rroi3 = movevars(Rroi3,{'TaskType','Contrast','Run','Phase'},'Before',1);
            AllTriROI{end+1,1} = Rroi3; %#ok<AGROW>

            AllStatus{end+1,1} = local_status('threeclass','color_vs_orientation_vs_conjunction',rr,phaseName,sum(idx),numel(unique(subjSC)),'ok'); %#ok<AGROW>

            fprintf('Done three-class: channel units=%d | ROI units=%d | subjects=%d\n', ...
                height(Rchan3),height(Rroi3),numel(unique(subjSC)));

        catch ME
            warning('Three-class failed: %s', ME.message);
            AllStatus{end+1,1} = local_status('threeclass','color_vs_orientation_vs_conjunction',rr,phaseName,NaN,NaN,['failed_' ME.message]); %#ok<AGROW>
        end
    end
end

%% ===================== SAVE =====================

BinaryChannel = local_vertcat_or_empty(AllBinChan);
BinaryROI     = local_vertcat_or_empty(AllBinROI);
ThreeChannel  = local_vertcat_or_empty(AllTriChan);
ThreeROI      = local_vertcat_or_empty(AllTriROI);
Status        = local_vertcat_or_empty(AllStatus);

BinaryChannel = local_sort_results(BinaryChannel,'binary');
BinaryROI     = local_sort_results(BinaryROI,'binary');
ThreeChannel  = local_sort_results(ThreeChannel,'three');
ThreeROI      = local_sort_results(ThreeROI,'three');

BinaryChannelSig = BinaryChannel(isfinite(BinaryChannel.q_FDR) & BinaryChannel.q_FDR < cfg.alphaFDR,:);
BinaryROISig     = BinaryROI(isfinite(BinaryROI.q_FDR) & BinaryROI.q_FDR < cfg.alphaFDR,:);
ThreeChannelSig  = ThreeChannel(isfinite(ThreeChannel.q_FDR) & ThreeChannel.q_FDR < cfg.alphaFDR,:);
ThreeROISig      = ThreeROI(isfinite(ThreeROI.q_FDR) & ThreeROI.q_FDR < cfg.alphaFDR,:);

BinaryChannelTop = local_top_by_task(BinaryChannel,cfg.topN);
BinaryROITop     = local_top_by_task(BinaryROI,cfg.topN);
ThreeChannelTop  = local_top_by_task(ThreeChannel,cfg.topN);
ThreeROITop      = local_top_by_task(ThreeROI,cfg.topN);

writetable(BinaryChannel, fullfile(cfg.outDir,'STEP100_binary_channel_specificFeature_all.csv'));
writetable(BinaryROI, fullfile(cfg.outDir,'STEP100_binary_ROI_specificFeature_all.csv'));
writetable(ThreeChannel, fullfile(cfg.outDir,'STEP100_threeclass_channel_specificFeature_all.csv'));
writetable(ThreeROI, fullfile(cfg.outDir,'STEP100_threeclass_ROI_specificFeature_all.csv'));

writetable(BinaryChannelSig, fullfile(cfg.outDir,'STEP100_binary_channel_specificFeature_significant.csv'));
writetable(BinaryROISig, fullfile(cfg.outDir,'STEP100_binary_ROI_specificFeature_significant.csv'));
writetable(ThreeChannelSig, fullfile(cfg.outDir,'STEP100_threeclass_channel_specificFeature_significant.csv'));
writetable(ThreeROISig, fullfile(cfg.outDir,'STEP100_threeclass_ROI_specificFeature_significant.csv'));

writetable(BinaryChannelTop, fullfile(cfg.outDir,'STEP100_binary_channel_specificFeature_topEffect.csv'));
writetable(BinaryROITop, fullfile(cfg.outDir,'STEP100_binary_ROI_specificFeature_topEffect.csv'));
writetable(ThreeChannelTop, fullfile(cfg.outDir,'STEP100_threeclass_channel_specificFeature_topEffect.csv'));
writetable(ThreeROITop, fullfile(cfg.outDir,'STEP100_threeclass_ROI_specificFeature_topEffect.csv'));
writetable(Status, fullfile(cfg.outDir,'STEP100_task_status.csv'));

Summary = table( ...
    ["binary_channel";"binary_ROI";"threeclass_channel";"threeclass_ROI"], ...
    [height(BinaryChannel);height(BinaryROI);height(ThreeChannel);height(ThreeROI)], ...
    [height(BinaryChannelSig);height(BinaryROISig);height(ThreeChannelSig);height(ThreeROISig)], ...
    'VariableNames',{'Analysis','N_all','N_significant_qFDR'});

writetable(Summary, fullfile(cfg.outDir,'STEP100_summary.csv'));

save(fullfile(cfg.outDir,'STEP100_results.mat'), ...
    'BinaryChannel','BinaryROI','ThreeChannel','ThreeROI', ...
    'BinaryChannelSig','BinaryROISig','ThreeChannelSig','ThreeROISig', ...
    'BinaryChannelTop','BinaryROITop','ThreeChannelTop','ThreeROITop', ...
    'Status','Summary','cfg','ChanUnits','ROIUnits','Finfo','-v7.3');

fprintf('\n=== STEP100 finished ===\n');
disp(Summary);
fprintf('Output folder:\n%s\n', cfg.outDir);
fprintf('Finished: %s\n', datestr(now));
diary off;

%% ========================================================================
%% LOCAL FUNCTIONS
%% ========================================================================

function [T,col,fNames] = local_load_dataset(cfg)
path = '';
for i=1:numel(cfg.datasetPaths)
    if exist(cfg.datasetPaths{i},'file'), path = cfg.datasetPaths{i}; break; end
end
if isempty(path), error('dataset.mat not found.'); end
fprintf('Loading dataset:\n%s\n', path);
S = load(path);
if isfield(S,'DS') && istable(S.DS)
    T = S.DS;
elseif isfield(S,'T') && istable(S.T)
    T = S.T;
else
    error('No DS or T table found.');
end
names = T.Properties.VariableNames;
col.subj  = local_findvar(names,{'Subject','subject','SubjectID','subjectID','subj','subjID','Subj','SubjID'});
col.run   = local_findvar(names,{'runNum','RunNum','run','Run'});
col.phase = local_findvar(names,{'phase','Phase'});
col.cond  = local_findvar(names,{'Condition','condition','yCondition','Label','label'});
if isempty(col.subj)||isempty(col.run)||isempty(col.phase)||isempty(col.cond)
    error('Missing metadata columns.');
end
fNames = local_detect_features(T);
end

function v = local_findvar(names,cands)
v = '';
for i=1:numel(cands)
    idx = strcmp(names,cands{i});
    if any(idx), v = names{find(idx,1)}; return; end
end
end

function f = local_detect_features(T)
names = T.Properties.VariableNames;
isNum = false(1,numel(names));
for j=1:numel(names)
    x = T.(names{j});
    isNum(j) = isnumeric(x) || islogical(x);
end
numNames = names(isNum);
meta = {'Subject','subject','SubjectID','subjectID','subj','subjID','Subj','SubjID','Run','run','runNum','RunNum','session','Session','Phase','phase','Condition','condition','yCondition','Label','label','Correct','correct','yCorrect','TrialNum','trialNum','TrialIndex','trialIndex','Trial','trial','PatternID','patternID','StartRow','EndRow','Second10Row','RetrRow','Fold','fold','CVFold'};
isMeta = false(size(numNames));
for j=1:numel(numNames), isMeta(j) = any(strcmpi(numNames{j},meta)); end
f = numNames(~isMeta);
end

function s = local_cleanstr(x)
if isnumeric(x)||islogical(x), s=string(x);
elseif iscell(x), s=string(x);
elseif iscategorical(x), s=string(x);
elseif isstring(x), s=x;
elseif ischar(x), s=string(cellstr(x));
else, s=string(x); end
s = lower(strtrim(s));
s = regexprep(s,'\s+','');
s = regexprep(s,'[^\w]','');
end

function Finfo = local_parse_channel_features(fNames,cfg)
n = numel(fNames);
FeatureIndex = (1:n)';
FeatureName = string(fNames(:));
Channel = NaN(n,1);
ROI = strings(n,1);
SpecificFeature = strings(n,1);
FeatureFamily = strings(n,1);
IsConnectivity = false(n,1);
IsChannelFeature = false(n,1);

for i=1:n
    fname = lower(char(FeatureName(i)));
    IsConnectivity(i) = local_is_connectivity(fname);
    ch = local_parse_channel(fname);
    Channel(i) = ch;
    ROI(i) = local_channel_to_roi(ch);
    SpecificFeature(i) = local_specific_feature(fname);
    FeatureFamily(i) = local_feature_family(fname);
    if isfinite(ch) && ismember(ch,cfg.channelsToUse)
        if ~(cfg.excludeConnectivityFeatures && IsConnectivity(i))
            IsChannelFeature(i) = true;
        end
    end
end

Finfo = table(FeatureIndex,FeatureName,Channel,ROI,SpecificFeature,FeatureFamily,IsConnectivity,IsChannelFeature);
Finfo = Finfo(Finfo.IsChannelFeature,:);
end

function tf = local_is_connectivity(fname)
s = string(fname);
tf = contains(s,"rie") || contains(s,"pli") || contains(s,"plv") || contains(s,"coh") || contains(s,"conn") || ...
     ~isempty(regexp(fname,'c0?\d{1,2}[_-]c0?\d{1,2}','once')) || ...
     ~isempty(regexp(fname,'ch0?\d{1,2}[_-]ch0?\d{1,2}','once'));
end

function ch = local_parse_channel(fname)
ch = NaN;
patterns = { ...
    '(?:^|[_-])ch0?(\d{1,2})(?:[_-]|$)', ...
    '(?:^|[_-])c0?(\d{1,2})(?:[_-]|$)', ...
    '(?:ch|c)0?(\d{1,2})(?:[_-]|$)' ...
};
for k=1:numel(patterns)
    tok = regexp(fname,patterns{k},'tokens');
    if ~isempty(tok)
        val = str2double(tok{1}{1});
        if isfinite(val) && val>=1 && val<=128
            ch = val; return;
        end
    end
end
end

function sf = local_specific_feature(fname)
s = string(fname);
% relative first
if contains(s,"rbp") || contains(s,"relative")
    if contains(s,"delta"), sf="relative_delta"; return; end
    if contains(s,"theta"), sf="relative_theta"; return; end
    if contains(s,"alpha"), sf="relative_alpha"; return; end
    if contains(s,"beta"),  sf="relative_beta";  return; end
    if contains(s,"gamma"), sf="relative_gamma"; return; end
    sf="relative_bandpower"; return;
end
if contains(s,"delta"), sf="delta"; return; end
if contains(s,"theta"), sf="theta"; return; end
if contains(s,"alpha"), sf="alpha"; return; end
if contains(s,"beta"),  sf="beta";  return; end
if contains(s,"gamma"), sf="gamma"; return; end
if contains(s,"rms"), sf="rms"; return; end
if contains(s,"peakabs"), sf="peak_abs"; return; end
if contains(s,"p2p"), sf="peak_to_peak"; return; end
if contains(s,"mean"), sf="mean"; return; end
if contains(s,"std"), sf="std"; return; end
if contains(s,"var"), sf="variance"; return; end
if contains(s,"deriv"), sf="derivative_rms"; return; end
if contains(s,"tkeo"), sf="tkeo"; return; end
if contains(s,"auc"), sf="auc_abs"; return; end
if contains(s,"linelen") || contains(s,"line"), sf="line_length"; return; end
if contains(s,"lzc") || contains(s,"lz"), sf="lzc"; return; end
if contains(s,"entropy"), sf="entropy"; return; end
if contains(s,"higuchi"), sf="higuchi"; return; end
if contains(s,"fractal"), sf="fractal"; return; end
sf = "unknown";
end

function ff = local_feature_family(fname)
sf = local_specific_feature(fname);
if any(sf == ["delta","theta","alpha","beta","gamma"])
    ff = "bandpower";
elseif startsWith(sf,"relative")
    ff = "relative_bandpower";
elseif any(sf == ["rms","peak_abs","peak_to_peak","mean","std","variance","derivative_rms","tkeo","auc_abs","line_length"])
    ff = "temporal_statistical";
elseif any(sf == ["lzc","entropy","higuchi","fractal"])
    ff = "complexity";
else
    ff = "unknown";
end
end

function roi = local_channel_to_roi(ch)
if ~isfinite(ch), roi="unknown"; return; end
if ismember(ch,1:10), roi="frontal";
elseif ismember(ch,11:20), roi="frontocentral";
elseif ismember(ch,21:34), roi="central";
elseif ismember(ch,35:48), roi="centroparietal";
elseif ismember(ch,49:58), roi="parietooccipital";
elseif ismember(ch,59:64), roi="occipital";
else, roi="unknown";
end
end

function Units = local_build_channel_units(Finfo)
[G,ch,roi,sf,ff] = findgroups(Finfo.Channel,Finfo.ROI,Finfo.SpecificFeature,Finfo.FeatureFamily);
FeatureIdx = splitapply(@(x){x},Finfo.FeatureIndex,G);
UnitID = "ch" + string(ch) + "_" + string(sf);
Units = table(UnitID,ch,string(roi),string(sf),string(ff),FeatureIdx, ...
    'VariableNames',{'UnitID','Channel','ROI','SpecificFeature','FeatureFamily','FeatureIdx'});
end

function Units = local_build_roi_units(Finfo)
[G,roi,sf,ff] = findgroups(Finfo.ROI,Finfo.SpecificFeature,Finfo.FeatureFamily);
FeatureIdx = splitapply(@(x){x},Finfo.FeatureIndex,G);
UnitID = string(roi) + "_" + string(sf);
Units = table(UnitID,string(roi),string(sf),string(ff),FeatureIdx, ...
    'VariableNames',{'UnitID','ROI','SpecificFeature','FeatureFamily','FeatureIdx'});
end

function [idx,y,subj] = local_binary_rows(subjAll,condAll,runAll,phaseAll,task,cfg)
switch string(task.contrastName)
    case "color_vs_orientation"
        A=cfg.labelMap.color; B=cfg.labelMap.orientation; Aname="color"; Bname="orientation";
    case "color_vs_conjunction"
        A=cfg.labelMap.color; B=cfg.labelMap.conjunction; Aname="color"; Bname="conjunction";
    case "orientation_vs_conjunction"
        A=cfg.labelMap.orientation; B=cfg.labelMap.conjunction; Aname="orientation"; Bname="conjunction";
    otherwise
        error('Unknown binary contrast.');
end
base = runAll==task.run & phaseAll==local_cleanstr({task.phase});
idxA = base & ismember(condAll,local_cleanstr(A));
idxB = base & ismember(condAll,local_cleanstr(B));
idx = idxA | idxB;
y = strings(sum(idx),1);
y(idxA(idx)) = Aname;
y(idxB(idx)) = Bname;
y = categorical(y);
subj = subjAll(idx);
end

function [idx,y,subj] = local_three_rows(subjAll,condAll,runAll,phaseAll,task,cfg)
base = runAll==task.run & phaseAll==local_cleanstr({task.phase});
idxC = base & ismember(condAll,local_cleanstr(cfg.labelMap.color));
idxO = base & ismember(condAll,local_cleanstr(cfg.labelMap.orientation));
idxJ = base & ismember(condAll,local_cleanstr(cfg.labelMap.conjunction));
idx = idxC | idxO | idxJ;
y = strings(sum(idx),1);
y(idxC(idx)) = "color";
y(idxO(idx)) = "orientation";
y(idxJ(idx)) = "conjunction";
y = categorical(y);
subj = subjAll(idx);
end

function XU = local_unit_matrix(Xraw,Units)
XU = NaN(size(Xraw,1),height(Units));
for u=1:height(Units)
    Xi = Xraw(:,Units.FeatureIdx{u});
    XU(:,u) = median(Xi,2,'omitnan');
end
end

function [Xs,ys,subjS] = local_subject_condition_samples(Xtrial,yTrial,subjTrial,cfg)
subjTrial = local_cleanstr(subjTrial);
subs = unique(subjTrial);
subs = subs(~ismissing(subs));
cats = categories(removecats(yTrial));
rows = {};
ys = categorical();
subjS = strings(0,1);
for s=1:numel(subs)
    ixS = subjTrial==subs(s);
    for c=1:numel(cats)
        cls = cats{c};
        ix = ixS & yTrial==cls;
        if ~any(ix), continue; end
        Xi = Xtrial(ix,:); Xi(~isfinite(Xi)) = NaN;
        if strcmpi(cfg.subjectConditionAggregate,'mean')
            row = mean(Xi,1,'omitnan');
        else
            row = median(Xi,1,'omitnan');
        end
        rows{end+1,1} = row; %#ok<AGROW>
        ys = [ys; categorical(string(cls))]; %#ok<AGROW>
        subjS(end+1,1) = subs(s); %#ok<AGROW>
    end
end
Xs = vertcat(rows{:});
ys = removecats(ys);
end

function R = local_binary_paired_stats(X,y,subj,Units,cfg)
cats = categories(removecats(y));
if numel(cats)~=2, error('Binary requires two classes.'); end
subs = unique(local_cleanstr(subj));
subs = subs(~ismissing(subs));
nU = size(X,2);
p = NaN(nU,1); eff = NaN(nU,1); medDiff = NaN(nU,1); nPairs = NaN(nU,1);
for u=1:nU
    d = [];
    for s=1:numel(subs)
        ixS = local_cleanstr(subj)==subs(s);
        xa = X(ixS & y==cats{1},u);
        xb = X(ixS & y==cats{2},u);
        if isempty(xa)||isempty(xb), continue; end
        xa=xa(1); xb=xb(1);
        if isfinite(xa)&&isfinite(xb), d(end+1,1)=xa-xb; end %#ok<AGROW>
    end
    nPairs(u)=numel(d);
    if numel(d)<5 || std(d,0,'omitnan')<=eps, continue; end
    try
        p(u)=signrank(d,0);
        medDiff(u)=median(d,'omitnan');
        eff(u)=medDiff(u)/max(mad(d,1),eps);
    catch
    end
end
q = local_bh_fdr(p);
R = local_base_unit_table(Units);
R.N_pairs = nPairs;
R.ClassA = repmat(string(cats{1}),height(R),1);
R.ClassB = repmat(string(cats{2}),height(R),1);
R.p_value = p;
R.q_FDR = q;
R.MedianDiff_AminusB = medDiff;
R.Effect = eff;
R.AbsEffect = abs(eff);
R.Significant_q05 = q < cfg.alphaFDR;
end

function R = local_three_friedman_stats(X,y,subj,Units,cfg)
wanted = ["color","orientation","conjunction"];
subs = unique(local_cleanstr(subj));
subs = subs(~ismissing(subs));
nU = size(X,2);
p = NaN(nU,1); eff = NaN(nU,1); nComplete = NaN(nU,1);
for u=1:nU
    M = [];
    for s=1:numel(subs)
        ixS = local_cleanstr(subj)==subs(s);
        row = NaN(1,3);
        for c=1:3
            val = X(ixS & string(y)==wanted(c),u);
            if ~isempty(val), row(c)=val(1); end
        end
        if all(isfinite(row)), M(end+1,:)=row; end %#ok<AGROW>
    end
    nComplete(u)=size(M,1);
    if size(M,1)<5 || std(M(:),0,'omitnan')<=eps, continue; end
    try
        [p(u),tbl]=friedman(M,1,'off');
        chi2=tbl{2,5};
        N=size(M,1); k=3;
        eff(u)=chi2/(N*(k-1));
    catch
    end
end
q = local_bh_fdr(p);
R = local_base_unit_table(Units);
R.N_subjects_complete = nComplete;
R.p_value = p;
R.q_FDR = q;
R.Effect = eff;
R.AbsEffect = abs(eff);
R.Significant_q05 = q < cfg.alphaFDR;
end

function R = local_base_unit_table(Units)
R = table();
R.UnitID = Units.UnitID;
if any(strcmp(Units.Properties.VariableNames,'Channel')), R.Channel = Units.Channel; end
if any(strcmp(Units.Properties.VariableNames,'ROI')), R.ROI = Units.ROI; end
if any(strcmp(Units.Properties.VariableNames,'SpecificFeature')), R.SpecificFeature = Units.SpecificFeature; end
if any(strcmp(Units.Properties.VariableNames,'FeatureFamily')), R.FeatureFamily = Units.FeatureFamily; end
end

function q = local_bh_fdr(p)
p = double(p(:)); q = NaN(size(p));
ok = isfinite(p) & p>=0 & p<=1;
p0 = p(ok); m = numel(p0);
if m==0, return; end
[ps,ord] = sort(p0,'ascend');
qs = ps .* m ./ (1:m)';
for i=m-1:-1:1, qs(i)=min(qs(i),qs(i+1)); end
qs(qs>1)=1;
tmp = NaN(m,1); tmp(ord)=qs; q(ok)=tmp;
end

function T = local_vertcat_or_empty(C)
if isempty(C), T = table(); else, T = vertcat(C{:}); end
end

function T = local_sort_results(T,mode)
if isempty(T), return; end
if any(strcmp(T.Properties.VariableNames,'q_FDR')) && any(strcmp(T.Properties.VariableNames,'AbsEffect'))
    T = sortrows(T,{'q_FDR','AbsEffect'},{'ascend','descend'});
elseif any(strcmp(T.Properties.VariableNames,'q_FDR'))
    T = sortrows(T,'q_FDR','ascend');
end
end

function Top = local_top_by_task(T,topN)
if isempty(T), Top=T; return; end
[G,~] = findgroups(T(:,{'TaskType','Contrast','Run','Phase'}));
rows = [];
for g=1:max(G)
    idx = find(G==g);
    Tg = T(idx,:);
    Tg2 = local_sort_results(Tg,'');
    % recover by matching UnitID positions in sorted table
    take = min(topN,height(Tg2));
    for k=1:take
        hit = idx(find(T.UnitID(idx)==Tg2.UnitID(k),1,'first'));
        rows(end+1,1)=hit; %#ok<AGROW>
    end
end
Top = T(rows,:);
end

function S = local_status(taskType,contrast,runNum,phaseName,nTrials,nSubjects,status)
S = table(string(taskType),string(contrast),runNum,string(phaseName),nTrials,nSubjects,string(status), ...
    'VariableNames',{'TaskType','Contrast','Run','Phase','N_trials','N_subjects','Status'});
end

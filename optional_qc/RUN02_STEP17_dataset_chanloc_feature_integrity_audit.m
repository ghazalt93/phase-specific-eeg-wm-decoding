function RUN02_STEP17_dataset_chanloc_feature_integrity_audit()
cfg = get_project_config();

clc;

%% ================= CONFIG =================
cfg.subjectRoot = fullfile(cfg.dataRoot, 'Subjects');
cfg.outBaseDir  = cfg.outputRoot;
cfg.outDir      = fullfile(cfg.outBaseDir, 'STEP17_dataset_chanloc_feature_integrity_audit');

if ~exist(cfg.outDir,'dir'), mkdir(cfg.outDir); end

cfg.dsCandidates = { ...
    fullfile(cfg.subjectRoot, '_wm_ml', 'dataset.mat'), ...
    fullfile(cfg.dataRoot, 'Features', 'Subjects', '_wm_ml', 'dataset.mat')};

cfg.chanlocsFile = '';

cfg.expectedLabels = {};

cfg.nEEGChan = 64;
cfg.eegChan  = 1:64;

cfg.workspaceRoots = { ...
    cfg.subjectRoot, ...
    fullfile(cfg.dataRoot, 'best', 'Subjects'), ...
    fullfile(cfg.dataRoot, 'Features', 'Subjects')};

cfg.excludeNameContains = {'eye','open','close','rest','baseline'};
cfg.validRuns = 1:3;

cfg.highNaNFeatureThr = 0.05;
cfg.constantStdThr    = 1e-12;
cfg.extremeAbsZThr    = 8;
cfg.highExtremeFracThr = 0.01;

cfg.familyRules.bandpower = {'bandpower','band_power','logbp','log_bp','bp_','_bp','power','psd','fft','delta','theta','alpha','beta','gamma'};
cfg.familyRules.morletTF  = {'morlet','wavelet','timefreq','time_freq','tf_','_tf','ersp','tfr'};
cfg.familyRules.temporalStat = {'temporal','mean','median','std','var','rms','skew','kurt','slope','peak','p2p','ptp','hjorth'};
cfg.familyRules.complexity = {'entropy','sampen','sampleentropy','perm_entropy','permentropy','higuchi','fractal','hurst','dfa','lz','lempel'};

fprintf('\n[STEP17] Dataset + chanloc + feature integrity audit\n');
fprintf('Output: %s\n', cfg.outDir);

%% ================= LOAD DATASET =================
cfg.DSpath = '';
for i=1:numel(cfg.dsCandidates)
    if exist(cfg.dsCandidates{i},'file')
        cfg.DSpath = cfg.dsCandidates{i};
        break;
    end
end
if isempty(cfg.DSpath)
    error('dataset.mat not found. Checked:\n%s', strjoin(cfg.dsCandidates, newline));
end

S = load(cfg.DSpath);
DS = get_table_from_mat(S);
meta = detect_meta(DS);
[yAll, runAll, phaseAll, subjAll] = get_core_vectors(DS, meta);
trialCol = detect_col(DS, {'trialnum','trial','trialindex','trial_idx','trialid'});
correctCol = detect_col(DS, {'ycorrect','correct','iscorrect','accuracy'});

[featVars, featCh] = detect_channel_features(DS, meta);
Family = detect_feature_families(featVars, cfg.familyRules);

fprintf('Dataset: %s\n', cfg.DSpath);
fprintf('Rows=%d | Vars=%d | Channel-coded numeric features=%d\n', height(DS), width(DS), numel(featVars));
fprintf('Meta: subject=%s | run=%s | phase=%s | label=%s | trial=%s | correct=%s\n', ...
    meta.subject, meta.run, meta.phase, meta.label, trialCol, correctCol);

%% ================= TABLES =================
OverviewT = make_overview(cfg, DS, meta, yAll, runAll, phaseAll, subjAll, featVars);

TrialCondT = make_trial_counts_condition(DS, yAll, runAll, phaseAll, subjAll, trialCol);
PhaseCountT = make_phase_counts(DS, yAll, runAll, phaseAll, subjAll, trialCol);

FeatureInvT = table(featVars(:), featCh(:), Family(:), ...
    'VariableNames', {'FeatureName','ChannelIndex','FeatureFamily'});
FeatureInvT = sortrows(FeatureInvT, {'ChannelIndex','FeatureFamily','FeatureName'});

FamilyCountT = summarize_family_counts(FeatureInvT);
ChannelFeatureCountT = summarize_channel_feature_counts(FeatureInvT, cfg);

FeatureQualityT = compute_feature_quality(DS, featVars, featCh, Family, cfg);

[ChanMapT, ChanLocSummaryT, chanMsg] = channel_location_audit(cfg);
FeatureFileInvT = feature_file_inventory(cfg);

%% ================= SAVE =================
writetable(OverviewT, fullfile(cfg.outDir, 'STEP17_dataset_overview.csv'));
writetable(TrialCondT, fullfile(cfg.outDir, 'STEP17_trial_counts_subject_run_phase_condition.csv'));
writetable(PhaseCountT, fullfile(cfg.outDir, 'STEP17_subject_run_phase_counts.csv'));
writetable(FeatureInvT, fullfile(cfg.outDir, 'STEP17_feature_inventory.csv'));
writetable(FamilyCountT, fullfile(cfg.outDir, 'STEP17_feature_family_counts.csv'));
writetable(ChannelFeatureCountT, fullfile(cfg.outDir, 'STEP17_channel_feature_counts.csv'));
writetable(FeatureQualityT, fullfile(cfg.outDir, 'STEP17_feature_quality_summary.csv'));

if ~isempty(ChanMapT)
    writetable(ChanMapT, fullfile(cfg.outDir, 'STEP17_channel_index_label_map.csv'));
end
if ~isempty(ChanLocSummaryT)
    writetable(ChanLocSummaryT, fullfile(cfg.outDir, 'STEP17_channel_location_quality.csv'));
end
if ~isempty(FeatureFileInvT)
    writetable(FeatureFileInvT, fullfile(cfg.outDir, 'STEP17_feature_file_inventory.csv'));
end

save(fullfile(cfg.outDir, 'STEP17_integrity_audit_results.mat'), ...
    'cfg','meta','OverviewT','TrialCondT','PhaseCountT','FeatureInvT','FamilyCountT', ...
    'ChannelFeatureCountT','FeatureQualityT','ChanMapT','ChanLocSummaryT','FeatureFileInvT','-v7.3');

write_report(fullfile(cfg.outDir, 'STEP17_integrity_report.txt'), cfg, DS, meta, ...
    OverviewT, TrialCondT, PhaseCountT, FeatureInvT, FamilyCountT, ...
    ChannelFeatureCountT, FeatureQualityT, ChanMapT, ChanLocSummaryT, FeatureFileInvT, chanMsg);

fprintf('\n[DONE STEP17]\n');
fprintf('Report:\n%s\n', fullfile(cfg.outDir, 'STEP17_integrity_report.txt'));
fprintf('\nMain files to send me:\n');
fprintf('  STEP17_integrity_report.txt\n');
fprintf('  STEP17_channel_index_label_map.csv\n');
fprintf('  STEP17_feature_quality_summary.csv\n');
fprintf('  STEP17_feature_family_counts.csv\n');
fprintf('  STEP17_channel_feature_counts.csv\n');

end

%% ================= BASIC HELPERS =================
function DS = get_table_from_mat(S)
if isfield(S,'DS'), DS=S.DS; return; end
if isfield(S,'dataset'), DS=S.dataset; return; end
if isfield(S,'T'), DS=S.T; return; end
fn = fieldnames(S);
for i=1:numel(fn)
    if istable(S.(fn{i}))
        DS = S.(fn{i});
        return;
    end
end
error('No table found in MAT file.');
end

function meta = detect_meta(T)
v = T.Properties.VariableNames;
l = lower(v);
meta.subject = pick(v,l,{'subject','subj','subjid','participant'});
meta.run     = pick(v,l,{'run','runnum','session'});
meta.phase   = pick(v,l,{'phase','phasename'});
meta.label   = pick(v,l,{'ycondition','condition','y','label'});
if isempty(meta.subject), error('Subject column not found.'); end
if isempty(meta.run), error('Run column not found.'); end
if isempty(meta.phase), error('Phase column not found.'); end
if isempty(meta.label), error('Condition/y label column not found.'); end
end

function out = pick(v,l,cands)
out = '';
for i=1:numel(cands)
    k = find(strcmp(l, lower(cands{i})), 1);
    if ~isempty(k), out = v{k}; return; end
end
for i=1:numel(cands)
    k = find(contains(l, lower(cands{i})), 1);
    if ~isempty(k), out = v{k}; return; end
end
end

function out = detect_col(T, cands)
v = T.Properties.VariableNames;
l = lower(v);
out = pick(v,l,cands);
end

function [y, runv, phasev, subjv] = get_core_vectors(DS, meta)
y = double(DS.(meta.label));
runv = double(DS.(meta.run));
phasev = lower(cellstr(string(DS.(meta.phase))));
subjv = cellstr(string(DS.(meta.subject)));
end

function [featVars, featCh] = detect_channel_features(DS, meta)
vars = DS.Properties.VariableNames;
n = height(DS);

exclude = false(1,numel(vars));
must = {meta.subject, meta.run, meta.phase, meta.label, ...
    'ycondition','condition','ycorrect','correct','trial','trialnum','patternid','rt','response','resp'};

for i=1:numel(vars)
    for j=1:numel(must)
        if strcmpi(vars{i}, must{j})
            exclude(i)=true;
        end
    end
end

featVars = {};
featCh = [];
for i=1:numel(vars)
    if exclude(i), continue; end
    ch = get_ch(vars{i});
    if isnan(ch), continue; end

    x = DS.(vars{i});
    if isnumeric(x) && isvector(x) && numel(x)==n
        featVars{end+1,1} = vars{i}; %#ok<AGROW>
        featCh(end+1,1) = ch; %#ok<AGROW>
    end
end
end

function ch = get_ch(name)
ch = NaN;
tok = regexp(name, 'ch[_-]?(\d{1,3})', 'tokens', 'once', 'ignorecase');
if ~isempty(tok)
    ch = str2double(tok{1});
    if ch<1 || ch>256, ch=NaN; end
end
end

%% ================= OVERVIEW / COUNTS =================
function T = make_overview(cfg, DS, meta, yAll, runAll, phaseAll, subjAll, featVars)
rows = {};
i=0;
i=i+1; rows(i,:)={'Dataset path', cfg.DSpath};
i=i+1; rows(i,:)={'Rows', height(DS)};
i=i+1; rows(i,:)={'Variables', width(DS)};
i=i+1; rows(i,:)={'Subjects', numel(unique(subjAll))};
i=i+1; rows(i,:)={'Subject list', strjoin(unique(subjAll,'stable'), ', ')};
i=i+1; rows(i,:)={'Runs', mat2str(unique(runAll(isfinite(runAll)))')};
i=i+1; rows(i,:)={'Phases', strjoin(unique(phaseAll,'stable'), ', ')};
i=i+1; rows(i,:)={'Condition codes', mat2str(unique(yAll(isfinite(yAll)))')};
i=i+1; rows(i,:)={'N condition 1', sum(yAll==1)};
i=i+1; rows(i,:)={'N condition 2', sum(yAll==2)};
i=i+1; rows(i,:)={'N condition 3', sum(yAll==3)};
i=i+1; rows(i,:)={'Channel-coded numeric features', numel(featVars)};
i=i+1; rows(i,:)={'Subject column', meta.subject};
i=i+1; rows(i,:)={'Run column', meta.run};
i=i+1; rows(i,:)={'Phase column', meta.phase};
i=i+1; rows(i,:)={'Label column', meta.label};
T = cell2table(rows, 'VariableNames', {'Item','Value'});
end

function Tcond = make_trial_counts_condition(DS, yAll, runAll, phaseAll, subjAll, trialCol)
[G,s,r,p,c] = findgroups(subjAll, runAll, phaseAll, yAll);
Nrows = splitapply(@numel, yAll, G);
NuniqueTrials = count_unique_trials(DS, trialCol, G, Nrows);
Tcond = table(s,r,p,c,Nrows,NuniqueTrials, ...
    'VariableNames', {'Subject','Run','Phase','Condition','Nrows','NuniqueTrials'});
Tcond = sortrows(Tcond, {'Subject','Run','Phase','Condition'});
end

function Tphase = make_phase_counts(DS, yAll, runAll, phaseAll, subjAll, trialCol)
[G,s,r,p] = findgroups(subjAll, runAll, phaseAll);
Nrows = splitapply(@numel, yAll, G);
NuniqueTrials = count_unique_trials(DS, trialCol, G, Nrows);
Ncond1 = splitapply(@(x) sum(x==1), yAll, G);
Ncond2 = splitapply(@(x) sum(x==2), yAll, G);
Ncond3 = splitapply(@(x) sum(x==3), yAll, G);
Tphase = table(s,r,p,Nrows,NuniqueTrials,Ncond1,Ncond2,Ncond3, ...
    'VariableNames', {'Subject','Run','Phase','Nrows','NuniqueTrials','Ncond1','Ncond2','Ncond3'});
Tphase = sortrows(Tphase, {'Subject','Run','Phase'});
end

function NuniqueTrials = count_unique_trials(DS, trialCol, G, fallback)
if isempty(trialCol)
    NuniqueTrials = fallback;
    return;
end
tv = DS.(trialCol);
if isnumeric(tv)
    NuniqueTrials = splitapply(@(x) numel(unique(x(isfinite(x)))), tv, G);
else
    tvs = cellstr(string(tv));
    NuniqueTrials = splitapply(@(x) numel(unique(x)), tvs, G);
end
end

%% ================= FEATURE FAMILY / QC =================
function Family = detect_feature_families(featVars, rules)
Family = repmat({'unknown'}, numel(featVars), 1);
fields = fieldnames(rules);
for i=1:numel(featVars)
    nm = lower(featVars{i});
    nm = regexprep(nm, 'ch[_-]?\d{1,3}', '');
    for f=1:numel(fields)
        fam = fields{f};
        pats = rules.(fam);
        hit = false;
        for p=1:numel(pats)
            if contains(nm, lower(pats{p}))
                hit = true; break;
            end
        end
        if hit
            Family{i}=fam; break;
        end
    end
end
end

function T = summarize_family_counts(FeatureInvT)
[G,fam] = findgroups(FeatureInvT.FeatureFamily);
Nfeatures = splitapply(@numel, FeatureInvT.FeatureName, G);
Nchannels = splitapply(@(x) numel(unique(x)), FeatureInvT.ChannelIndex, G);
T = table(fam,Nfeatures,Nchannels, 'VariableNames', {'FeatureFamily','Nfeatures','Nchannels'});
T = sortrows(T,'Nfeatures','descend');
end

function T = summarize_channel_feature_counts(FeatureInvT, cfg)
ch = (1:cfg.nEEGChan)';
Nfeatures = zeros(cfg.nEEGChan,1);
Nbandpower = zeros(cfg.nEEGChan,1);
NmorletTF = zeros(cfg.nEEGChan,1);
NtemporalStat = zeros(cfg.nEEGChan,1);
Ncomplexity = zeros(cfg.nEEGChan,1);
Nunknown = zeros(cfg.nEEGChan,1);

for i=1:cfg.nEEGChan
    idx = FeatureInvT.ChannelIndex==i;
    Nfeatures(i)=sum(idx);
    Nbandpower(i)=sum(idx & strcmp(FeatureInvT.FeatureFamily,'bandpower'));
    NmorletTF(i)=sum(idx & strcmp(FeatureInvT.FeatureFamily,'morletTF'));
    NtemporalStat(i)=sum(idx & strcmp(FeatureInvT.FeatureFamily,'temporalStat'));
    Ncomplexity(i)=sum(idx & strcmp(FeatureInvT.FeatureFamily,'complexity'));
    Nunknown(i)=sum(idx & strcmp(FeatureInvT.FeatureFamily,'unknown'));
end
MissingFromDataset = Nfeatures==0;
T = table(ch,Nfeatures,Nbandpower,NmorletTF,NtemporalStat,Ncomplexity,Nunknown,MissingFromDataset, ...
    'VariableNames', {'ChannelIndex','Nfeatures','Nbandpower','NmorletTF','NtemporalStat','Ncomplexity','Nunknown','MissingFromDataset'});
end

function T = compute_feature_quality(DS, featVars, featCh, Family, cfg)
n = numel(featVars);
NaNFrac = NaN(n,1); InfFrac = NaN(n,1); FiniteFrac = NaN(n,1);
MeanVal = NaN(n,1); MedianVal = NaN(n,1); StdVal = NaN(n,1); MADVal = NaN(n,1);
MinVal = NaN(n,1); MaxVal = NaN(n,1); ExtremeFrac = NaN(n,1);

for i=1:n
    x = double(DS.(featVars{i}));
    NaNFrac(i)=mean(isnan(x));
    InfFrac(i)=mean(isinf(x));
    FiniteFrac(i)=mean(isfinite(x));
    xf = x(isfinite(x));
    if isempty(xf), continue; end
    MeanVal(i)=mean(xf);
    MedianVal(i)=median(xf);
    StdVal(i)=std(xf);
    MADVal(i)=median(abs(xf-median(xf)));
    MinVal(i)=min(xf);
    MaxVal(i)=max(xf);
    rz = robust_z(xf);
    ExtremeFrac(i)=mean(abs(rz)>cfg.extremeAbsZThr);
end

Flag_HighNaN = (NaNFrac+InfFrac)>cfg.highNaNFeatureThr;
Flag_Constant = ~isfinite(StdVal) | StdVal<cfg.constantStdThr;
Flag_ExtremeFrac = ExtremeFrac>cfg.highExtremeFracThr;

T = table(featVars(:),featCh(:),Family(:),NaNFrac,InfFrac,FiniteFrac,MeanVal,MedianVal,StdVal,MADVal,MinVal,MaxVal,ExtremeFrac, ...
    Flag_HighNaN,Flag_Constant,Flag_ExtremeFrac, ...
    'VariableNames', {'FeatureName','ChannelIndex','FeatureFamily','NaNFrac','InfFrac','FiniteFrac','Mean','Median','Std','MAD','Min','Max','ExtremeFrac','Flag_HighNaN','Flag_Constant','Flag_ExtremeFrac'});
T = sortrows(T, {'Flag_HighNaN','Flag_Constant','Flag_ExtremeFrac','ExtremeFrac'}, {'descend','descend','descend','descend'});
end

function z = robust_z(x)
m = median(x);
mad0 = median(abs(x-m));
if ~isfinite(mad0) || mad0<eps
    s = std(x);
    if ~isfinite(s) || s<eps, s=1; end
    z = (x-m)/s;
else
    z = (x-m)/(1.4826*mad0);
end
end

%% ================= CHANNEL LOCATION AUDIT =================
function [ChanMapT, SummaryT, msg] = channel_location_audit(cfg)
ChanMapT = table(); SummaryT = table(); msg = '';

edfPath = find_example_edf(cfg);
if isempty(edfPath)
    msg = 'No task EDF found; channel-location audit skipped.';
    return;
end

if exist('pop_biosig','file') ~= 2
    try, eeglab nogui; catch, end
end
if exist('pop_biosig','file') ~= 2
    msg = 'pop_biosig not found; channel-location audit skipped.';
    return;
end

try
    EEG = pop_biosig(edfPath, 'channels', cfg.eegChan, 'importevent', 'off');
catch
    try
        EEG = pop_biosig(edfPath, 'channels', cfg.eegChan);
    catch ME
        msg = ['Could not load example EDF: ' ME.message];
        return;
    end
end

if ~isempty(cfg.chanlocsFile) && exist(cfg.chanlocsFile,'file') && exist('pop_chanedit','file')==2
    try
        EEG = pop_chanedit(EEG, 'load', {cfg.chanlocsFile, 'filetype', 'autodetect'});
        msg = ['Chanlocs loaded from cfg.chanlocsFile: ' cfg.chanlocsFile];
    catch ME
        msg = ['Could not apply cfg.chanlocsFile: ' ME.message];
    end
else
    msg = 'Using chanlocs from EDF/EEGLAB import. If labels are generic, set cfg.chanlocsFile.';
end

nCh = min(cfg.nEEGChan, size(EEG.data,1));
ChannelIndex = (1:nCh)';
ChannelLabel = strings(nCh,1);
X=NaN(nCh,1); Y=NaN(nCh,1); Z=NaN(nCh,1); Theta=NaN(nCh,1); Radius=NaN(nCh,1);

for i=1:nCh
    lab = sprintf('Ch%d', i);
    if isfield(EEG,'chanlocs') && numel(EEG.chanlocs)>=i
        if isfield(EEG.chanlocs,'labels') && ~isempty(EEG.chanlocs(i).labels)
            lab = EEG.chanlocs(i).labels;
        end
        if isfield(EEG.chanlocs,'X') && ~isempty(EEG.chanlocs(i).X), X(i)=EEG.chanlocs(i).X; end
        if isfield(EEG.chanlocs,'Y') && ~isempty(EEG.chanlocs(i).Y), Y(i)=EEG.chanlocs(i).Y; end
        if isfield(EEG.chanlocs,'Z') && ~isempty(EEG.chanlocs(i).Z), Z(i)=EEG.chanlocs(i).Z; end
        if isfield(EEG.chanlocs,'theta') && ~isempty(EEG.chanlocs(i).theta), Theta(i)=EEG.chanlocs(i).theta; end
        if isfield(EEG.chanlocs,'radius') && ~isempty(EEG.chanlocs(i).radius), Radius(i)=EEG.chanlocs(i).radius; end
    end
    ChannelLabel(i)=string(lab);
end

HasLabel = strlength(ChannelLabel)>0;
IsGenericLabel = startsWith(lower(ChannelLabel),"ch") | startsWith(lower(ChannelLabel),"chan");
HasXYZ = isfinite(X)&isfinite(Y)&isfinite(Z);
HasThetaRadius = isfinite(Theta)&isfinite(Radius);
HasAnyCoordinate = HasXYZ | HasThetaRadius;

DuplicateLabel = false(nCh,1);
for i=1:nCh
    DuplicateLabel(i)=sum(strcmpi(ChannelLabel,ChannelLabel(i)))>1;
end

ExpectedLabel = strings(nCh,1);
ExpectedProvided = false(nCh,1);
MatchesExpected = false(nCh,1);
if ~isempty(cfg.expectedLabels)
    for i=1:min(nCh,numel(cfg.expectedLabels))
        ExpectedLabel(i)=string(cfg.expectedLabels{i});
        ExpectedProvided(i)=true;
        MatchesExpected(i)=strcmpi(ChannelLabel(i),ExpectedLabel(i));
    end
end

ChanMapT = table(ChannelIndex,ChannelLabel,ExpectedLabel,ExpectedProvided,MatchesExpected, ...
    X,Y,Z,Theta,Radius,HasLabel,IsGenericLabel,DuplicateLabel,HasXYZ,HasThetaRadius,HasAnyCoordinate, ...
    'VariableNames', {'ChannelIndex','ChannelLabel','ExpectedLabel','ExpectedProvided','MatchesExpected','X','Y','Z','Theta','Radius','HasLabel','IsGenericLabel','DuplicateLabel','HasXYZ','HasThetaRadius','HasAnyCoordinate'});

rows = {};
r=0;
r=r+1; rows(r,:)={'Example EDF', edfPath};
r=r+1; rows(r,:)={'Chanloc status', msg};
r=r+1; rows(r,:)={'N channels checked', nCh};
r=r+1; rows(r,:)={'N generic labels', sum(IsGenericLabel)};
r=r+1; rows(r,:)={'N duplicate labels', sum(DuplicateLabel)};
r=r+1; rows(r,:)={'N channels with any coordinate', sum(HasAnyCoordinate)};
if any(ExpectedProvided)
    r=r+1; rows(r,:)={'Expected labels provided', sum(ExpectedProvided)};
    r=r+1; rows(r,:)={'Expected label matches', sum(MatchesExpected & ExpectedProvided)};
    r=r+1; rows(r,:)={'Expected label mismatches', sum(~MatchesExpected & ExpectedProvided)};
else
    r=r+1; rows(r,:)={'Expected label check', 'Skipped; cfg.expectedLabels is empty'};
end
SummaryT = cell2table(rows, 'VariableNames', {'Item','Value'});
end

function edfPath = find_example_edf(cfg)
edfPath = '';
WsT = index_task_workspaces_for_edf(cfg);
if ~isempty(WsT)
    edfPath = WsT.EDFPath{1};
end
end

function WsT = index_task_workspaces_for_edf(cfg)
wsPaths = {};
for r=1:numel(cfg.workspaceRoots)
    root=cfg.workspaceRoots{r};
    if ~exist(root,'dir'), continue; end
    d=dir(fullfile(root,'**','workspace*.mat'));
    for i=1:numel(d)
        wsPaths{end+1,1}=fullfile(d(i).folder,d(i).name); %#ok<AGROW>
    end
end
wsPaths = unique(wsPaths,'stable');
rows={}; ri=0;
for i=1:numel(wsPaths)
    wsPath=wsPaths{i};
    [subjDir,wsName]=fileparts(wsPath);
    if should_exclude_name(wsName,cfg.excludeNameContains), continue; end
    subj=infer_subject(wsPath);
    runNum=infer_run(wsName);
    if ~(isfinite(runNum)&&ismember(runNum,cfg.validRuns)), continue; end
    try
        info=whos('-file',wsPath); names={info.name};
        if ~any(strcmp(names,'T_phases')), continue; end
        W=load_minimal_workspace(wsPath);
        if ~isfield(W,'T_phases') || ~istable(W.T_phases), continue; end
        edf=find_edf_for_workspace(wsPath,subjDir,W,cfg);
        if isempty(edf) || ~isfile(edf), continue; end
        ri=ri+1;
        rows(ri,:)={char(subj),runNum,wsName,wsPath,edf}; %#ok<AGROW>
    catch
        continue;
    end
end
if isempty(rows)
    WsT=table();
else
    WsT=cell2table(rows,'VariableNames',{'Subject','Run','Workspace','WorkspacePath','EDFPath'});
    WsT=sortrows(WsT,{'Subject','Run'});
end
end

function W = load_minimal_workspace(wsPath)
varsWanted={'T_phases','filename1','filename2','filename','edfPath','EDFPath','fs'};
info=whos('-file',wsPath); names={info.name};
toLoad={};
for i=1:numel(varsWanted)
    if any(strcmp(names,varsWanted{i})), toLoad{end+1}=varsWanted{i}; end %#ok<AGROW>
end
W=load(wsPath,toLoad{:});
end

function edfPath = find_edf_for_workspace(wsPath, subjDir, W, cfg)
edfPath='';
candidateVars={'filename1','filename2','filename','edfPath','EDFPath'};
for cv=1:numel(candidateVars)
    v=candidateVars{cv};
    if isfield(W,v)
        candList=normalize_candidate_to_cell(W.(v));
        for ci=1:numel(candList)
            p=resolve_edf_candidate(strtrim(candList{ci}),subjDir);
            if ~isempty(p) && isfile(p) && endsWith(lower(p),'.edf') && ~should_exclude_name(p,cfg.excludeNameContains)
                edfPath=p; return;
            end
        end
    end
end

[~,wsBase]=fileparts(wsPath);
wsBase=regexprep(wsBase,'^workspace[_ -]?','');
wsClean=clean_name(wsBase);
wsNoRun=regexprep(wsClean,'\d+$','');
wsRun=infer_run(wsBase);

edfs=[dir(fullfile(subjDir,'*.edf')); dir(fullfile(subjDir,'*.EDF'))];
if isempty(edfs), edfs=[dir(fullfile(subjDir,'**','*.edf')); dir(fullfile(subjDir,'**','*.EDF'))]; end
if isempty(edfs), return; end

bestScore=-inf; bestFile='';
for i=1:numel(edfs)
    eName=edfs(i).name;
    ePath=fullfile(edfs(i).folder,edfs(i).name);
    if ~endsWith(lower(ePath),'.edf'), continue; end
    if should_exclude_name(eName,cfg.excludeNameContains), continue; end
    if startsWith(lower(eName),'workspace'), continue; end
    eBase=regexprep(eName,'\.edf$','','ignorecase');
    eClean=clean_name(eBase);
    eNoRun=regexprep(eClean,'\d+$','');
    eRun=infer_run(eBase);
    score=0;
    if strcmp(eClean,wsClean), score=score+300; end
    if strcmp(eNoRun,wsNoRun)&&isfinite(wsRun)&&isfinite(eRun)&&wsRun==eRun, score=score+180; end
    if contains(eClean,wsClean)||contains(wsClean,eClean), score=score+100; end
    if isfinite(wsRun)&&isfinite(eRun)&&wsRun==eRun, score=score+30; end
    if score>bestScore, bestScore=score; bestFile=ePath; end
end
if bestScore>0, edfPath=bestFile; end
end

function c = normalize_candidate_to_cell(x)
c={};
if isempty(x), return; end
if isstring(x), x=cellstr(x); end
if iscell(x)
    for i=1:numel(x)
        if ischar(x{i})||isstring(x{i}), c{end+1}=char(x{i}); end %#ok<AGROW>
    end
elseif ischar(x)
    c{1}=x;
end
end

function p = resolve_edf_candidate(cand, subjDir)
p='';
if isempty(cand), return; end
if endsWith(lower(cand),'.mat'), return; end
poss={cand, fullfile(subjDir,cand)};
if ~endsWith(lower(cand),'.edf')
    poss{end+1}=[cand '.edf'];
    poss{end+1}=fullfile(subjDir,[cand '.edf']);
end
for i=1:numel(poss)
    if isfile(poss{i}) && endsWith(lower(poss{i}),'.edf')
        p=poss{i}; return;
    end
end
end

%% ================= FEATURE FILE INVENTORY =================
function T = feature_file_inventory(cfg)
roots={fullfile(cfg.subjectRoot,'_wm_features'), ...
       fullfile(cfg.dataRoot, 'Features', 'Subjects', '_wm_features'), ...
       fullfile(cfg.dataRoot, 'Features', 'Subjects', '_logbp_blwin_features')};
rows={}; ri=0;
for r=1:numel(roots)
    root=roots{r};
    if ~exist(root,'dir'), continue; end
    d=dir(fullfile(root,'**','features.mat'));
    for i=1:numel(d)
        p=fullfile(d(i).folder,d(i).name);
        subj=char(infer_subject(p)); runNum=infer_run(p); phase=infer_phase(p);
        loadOK=false; nRows=NaN; nVars=NaN; vars=''; msg='';
        try
            info=whos('-file',p); vars=strjoin({info.name},', ');
            tmp=load(p); loadOK=true;
            fn=fieldnames(tmp);
            for k=1:numel(fn)
                obj=tmp.(fn{k});
                if istable(obj), nRows=height(obj); nVars=width(obj); break; end
                if isnumeric(obj)&&ismatrix(obj), nRows=size(obj,1); nVars=size(obj,2); end
            end
        catch ME
            msg=ME.message;
        end
        ri=ri+1;
        rows(ri,:)={root,p,subj,runNum,phase,loadOK,nRows,nVars,vars,msg}; %#ok<AGROW>
    end
end
if isempty(rows)
    T=table();
else
    T=cell2table(rows,'VariableNames',{'FeatureRoot','FeaturePath','Subject','Run','Phase','LoadOK','NrowsOrDim1','NvarsOrDim2','VariablesInMat','Message'});
    T=sortrows(T,{'Subject','Run','Phase'});
end
end

function ph = infer_phase(p)
lp=lower(p);
if contains(lp,'/stim/')||contains(lp,'\stim\'), ph='stim';
elseif contains(lp,'/maint/')||contains(lp,'\maint\'), ph='maint';
elseif contains(lp,'/retr/')||contains(lp,'\retr\'), ph='retr';
else, ph=''; end
end

%% ================= REPORT =================
function write_report(reportPath, cfg, DS, meta, OverviewT, TrialCondT, PhaseCountT, FeatureInvT, FamilyCountT, ChannelFeatureCountT, FeatureQualityT, ChanMapT, ChanLocSummaryT, FeatureFileInvT, chanMsg)
fid=fopen(reportPath,'w');

fprintf(fid,'STEP17 DATASET + CHANLOC + FEATURE INTEGRITY AUDIT\n');
fprintf(fid,'=================================================\n\n');
fprintf(fid,'Dataset: %s\n', cfg.DSpath);
fprintf(fid,'Rows=%d | Vars=%d\n', height(DS), width(DS));
fprintf(fid,'Meta: subject=%s | run=%s | phase=%s | label=%s\n\n', meta.subject, meta.run, meta.phase, meta.label);

fprintf(fid,'Dataset overview:\n');
for i=1:height(OverviewT)
    fprintf(fid,'  %s: %s\n', OverviewT.Item{i}, string(OverviewT.Value(i)));
end

fprintf(fid,'\nChannel-location audit:\n');
fprintf(fid,'  %s\n', chanMsg);
if isempty(ChanMapT)
    fprintf(fid,'  Channel map not created.\n');
else
    fprintf(fid,'  N channels mapped: %d\n', height(ChanMapT));
    fprintf(fid,'  Generic labels: %d\n', sum(ChanMapT.IsGenericLabel));
    fprintf(fid,'  Duplicate labels: %d\n', sum(ChanMapT.DuplicateLabel));
    fprintf(fid,'  With any coordinate: %d\n', sum(ChanMapT.HasAnyCoordinate));
    if any(ChanMapT.ExpectedProvided)
        fprintf(fid,'  Expected label mismatches: %d\n', sum(~ChanMapT.MatchesExpected & ChanMapT.ExpectedProvided));
    else
        fprintf(fid,'  Expected label check skipped because cfg.expectedLabels is empty.\n');
    end
end

fprintf(fid,'\nImportant note:\n');
fprintf(fid,'  ch10 in feature names means EEG.data row/index 10 during feature extraction.\n');
fprintf(fid,'  It maps to Fz/Cz/etc only if EEG.chanlocs(10).labels is correct and chanloc ordering matches EDF order.\n\n');

fprintf(fid,'Feature inventory:\n');
fprintf(fid,'  Total numeric channel-coded features: %d\n', height(FeatureInvT));
fprintf(fid,'  Channels represented: %d\n', numel(unique(FeatureInvT.ChannelIndex)));
missing = ChannelFeatureCountT.ChannelIndex(ChannelFeatureCountT.MissingFromDataset);
if isempty(missing)
    fprintf(fid,'  No missing channels among 1:%d in feature names.\n', cfg.nEEGChan);
else
    fprintf(fid,'  Missing channels among 1:%d: %s\n', cfg.nEEGChan, mat2str(missing'));
end

fprintf(fid,'\nFeature family counts:\n');
for i=1:height(FamilyCountT)
    fprintf(fid,'  %s: %d features across %d channels\n', FamilyCountT.FeatureFamily{i}, FamilyCountT.Nfeatures(i), FamilyCountT.Nchannels(i));
end

fprintf(fid,'\nFeature quality warnings:\n');
fprintf(fid,'  High NaN/Inf features: %d\n', sum(FeatureQualityT.Flag_HighNaN));
fprintf(fid,'  Constant features: %d\n', sum(FeatureQualityT.Flag_Constant));
fprintf(fid,'  High extreme-value fraction features: %d\n', sum(FeatureQualityT.Flag_ExtremeFrac));

warn = FeatureQualityT(FeatureQualityT.Flag_HighNaN | FeatureQualityT.Flag_Constant | FeatureQualityT.Flag_ExtremeFrac,:);
if ~isempty(warn)
    fprintf(fid,'\nTop warning features:\n');
    for i=1:min(30,height(warn))
        fprintf(fid,'  %s | ch%d | %s | NaN+Inf=%.4f | Std=%.4g | ExtremeFrac=%.4f\n', ...
            warn.FeatureName{i}, warn.ChannelIndex(i), warn.FeatureFamily{i}, ...
            warn.NaNFrac(i)+warn.InfFrac(i), warn.Std(i), warn.ExtremeFrac(i));
    end
end

fprintf(fid,'\nTrial-count summary:\n');
fprintf(fid,'  Subject x Run x Phase x Condition rows: %d\n', height(TrialCondT));
fprintf(fid,'  Subject x Run x Phase rows: %d\n', height(PhaseCountT));
fprintf(fid,'  Min rows per Subject x Run x Phase: %d\n', min(PhaseCountT.Nrows));
fprintf(fid,'  Median rows per Subject x Run x Phase: %.1f\n', median(PhaseCountT.Nrows));
fprintf(fid,'  Max rows per Subject x Run x Phase: %d\n', max(PhaseCountT.Nrows));

fprintf(fid,'\nFeature file inventory:\n');
fprintf(fid,'  features.mat files found: %d\n', height(FeatureFileInvT));
if ~isempty(FeatureFileInvT)
    fprintf(fid,'  Load failures: %d\n', sum(~FeatureFileInvT.LoadOK));
end

fprintf(fid,'\nRecommended actions:\n');
fprintf(fid,'  1) If labels are generic or coordinates are missing, set cfg.chanlocsFile and rerun.\n');
fprintf(fid,'  2) If expected channel order is known, fill cfg.expectedLabels and rerun.\n');
fprintf(fid,'  3) If many features are NaN/constant/extreme, inspect feature extraction and run STEP15.\n');
fprintf(fid,'  4) If channels are missing from feature names, inspect feature extraction completeness.\n');

fclose(fid);
end

%% ================= PATH HELPERS =================
function tf = should_exclude_name(name,pats)
tf=false; ln=lower(char(name));
for i=1:numel(pats)
    if contains(ln,lower(pats{i})), tf=true; return; end
end
end

function subj = infer_subject(pathstr)
parts=split(string(pathstr),filesep);
subj="UNKNOWN";
for i=1:numel(parts)
    x=char(parts(i));
    if ~isempty(regexp(x,'^s\d+$','once')), subj=string(x); return; end
end
tok=regexp(pathstr,'(s\d+)','tokens','once','ignorecase');
if ~isempty(tok), subj=string(tok{1}); end
end

function runNum = infer_run(name)
runNum=NaN; name=char(name);
tok=regexp(name,'(?:run|session|sess|block)[_-]?(\d+)','tokens','once','ignorecase');
if ~isempty(tok), runNum=str2double(tok{1}); return; end
tok=regexp(name,'(\d+)$','tokens','once');
if ~isempty(tok), runNum=str2double(tok{1}); return; end
tok=regexp(name,'(\d+)','tokens');
if ~isempty(tok), runNum=str2double(tok{end}{1}); end
end

function out = clean_name(x)
out=lower(char(x));
out=regexprep(out,'\.edf$','');
out=regexprep(out,'_filtered$','');
out=regexprep(out,'-filtered$','');
out=regexprep(out,'[^a-z0-9]+','_');
out=regexprep(out,'^_+|_+$','');
end

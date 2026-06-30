function RUN03_STEP18_s21_duplicate_missingness_audit()

cfg = get_project_config();

clc;

%% ================= CONFIG =================
cfg.targetSubject = 's21';

cfg.subjectRoot = fullfile(cfg.dataRoot, 'Subjects');

cfg.dsCandidates = { ...
    fullfile(cfg.subjectRoot, '_wm_ml', 'dataset.mat'), ...
    fullfile(cfg.dataRoot, 'Features', 'Subjects', '_wm_ml', 'dataset.mat')};

cfg.outBaseDir = cfg.outputRoot;
cfg.outDir = fullfile(cfg.outBaseDir, 'STEP18_s21_duplicate_missingness_audit');

if ~exist(cfg.outBaseDir,'dir'), mkdir(cfg.outBaseDir); end
if ~exist(cfg.outDir,'dir'), mkdir(cfg.outDir); end

% Feature detection
cfg.excludeBehavior = true;
cfg.expectedEEGChannels = 1:64;

% QC thresholds
cfg.highFeatureNaNRate_subject = 0.30;  
cfg.highFeatureNaNRate_group   = 0.30;  
cfg.highFeatureNaNRate_feature = 0.05;   
cfg.constantStdThr = 1e-12;
cfg.extremeAbsRobustZ = 8;
cfg.highExtremeFracThr = 0.01;

% Duplicate similarity thresholds
cfg.identicalTolerance = 1e-10;
cfg.nearlyIdenticalCorrelation = 0.9999;

fprintf('\n[STEP18] s21 DUPLICATE / MISSINGNESS AUDIT\n');
fprintf('Target subject: %s\n', cfg.targetSubject);
fprintf('Output: %s\n', cfg.outDir);

%% ================= LOAD DATASET =================
cfg.DSpath = '';
for i=1:numel(cfg.dsCandidates)
    if exist(cfg.dsCandidates{i}, 'file')
        cfg.DSpath = cfg.dsCandidates{i};
        break;
    end
end

if isempty(cfg.DSpath)
    error('Dataset not found. Checked:\n%s', strjoin(cfg.dsCandidates, newline));
end

fprintf('Dataset: %s\n', cfg.DSpath);

S = load(cfg.DSpath);
DS = get_table_from_mat(S);
meta = detect_meta(DS);

[yAll, runAll, phaseAll, subjAll] = get_core_vectors(DS, meta);

trialCol = detect_col(DS, {'trialnum','trial','trialindex','trial_idx','trialid'});
if isempty(trialCol)
    warning('TrialNum column not found. Duplicate trial audit will use row index fallback.');
end

correctCol = detect_col(DS, {'ycorrect','correct','iscorrect','accuracy'});

[featVars, featCh] = detect_features(DS, meta, cfg);

fprintf('\n[DATASET]\n');
fprintf('Rows=%d | Vars=%d | EEG features=%d\n', height(DS), width(DS), numel(featVars));
fprintf('Meta: subject=%s | run=%s | phase=%s | condition=%s | trial=%s | correct=%s\n', ...
    meta.subject, meta.run, meta.phase, meta.label, trialCol, correctCol);

%% ================= BASIC OVERVIEW =================
isTarget = strcmp(subjAll, cfg.targetSubject);
if ~any(isTarget)
    error('Subject %s was not found in dataset.', cfg.targetSubject);
end

OverviewT = make_overview(DS, cfg, meta, isTarget, yAll, runAll, phaseAll, subjAll, featVars);

%% ================= COUNTS =================
CountT = make_run_phase_condition_counts(DS, cfg, isTarget, yAll, runAll, phaseAll, trialCol);

%% ================= DUPLICATE TRIAL GROUPS =================
[DupGroupT, DupSimilarityT] = audit_duplicates(DS, cfg, isTarget, yAll, runAll, phaseAll, trialCol, featVars);

%% ================= MISSINGNESS =================
MissingT = compute_missingness_by_subject_run_phase(DS, yAll, runAll, phaseAll, subjAll, featVars);
SubjectMissingT = summarize_subject_missingness(MissingT);

%% ================= FEATURE QUALITY FOR s21 =================
FeatureQualityT = compute_target_feature_quality(DS, cfg, isTarget, featVars, featCh);

%% ================= DECISION FLAGS =================
Decision = make_decision_struct(cfg, OverviewT, CountT, DupGroupT, DupSimilarityT, MissingT, SubjectMissingT, FeatureQualityT);

%% ================= SAVE =================
f_overview = fullfile(cfg.outDir, 'STEP18_s21_overview.csv');
f_count = fullfile(cfg.outDir, 'STEP18_s21_run_phase_condition_counts.csv');
f_dup = fullfile(cfg.outDir, 'STEP18_s21_duplicate_trial_groups.csv');
f_sim = fullfile(cfg.outDir, 'STEP18_s21_duplicate_similarity.csv');
f_miss = fullfile(cfg.outDir, 'STEP18_missingness_by_subject_run_phase.csv');
f_submiss = fullfile(cfg.outDir, 'STEP18_subject_missingness_summary.csv');
f_featq = fullfile(cfg.outDir, 'STEP18_s21_feature_quality.csv');
f_mat = fullfile(cfg.outDir, 'STEP18_s21_audit_results.mat');
f_report = fullfile(cfg.outDir, 'STEP18_s21_decision_report.txt');

writetable(OverviewT, f_overview);
writetable(CountT, f_count);
writetable(DupGroupT, f_dup);
writetable(DupSimilarityT, f_sim);
writetable(MissingT, f_miss);
writetable(SubjectMissingT, f_submiss);
writetable(FeatureQualityT, f_featq);

save(f_mat, 'cfg', 'meta', 'OverviewT', 'CountT', 'DupGroupT', 'DupSimilarityT', ...
    'MissingT', 'SubjectMissingT', 'FeatureQualityT', 'Decision', '-v7.3');

write_report(f_report, cfg, OverviewT, CountT, DupGroupT, DupSimilarityT, ...
    MissingT, SubjectMissingT, FeatureQualityT, Decision);

fprintf('\n[DONE STEP18]\n');
fprintf('Overview:     %s\n', f_overview);
fprintf('Counts:       %s\n', f_count);
fprintf('Duplicates:   %s\n', f_dup);
fprintf('Similarity:   %s\n', f_sim);
fprintf('Missingness:  %s\n', f_miss);
fprintf('Subject miss: %s\n', f_submiss);
fprintf('Feature QC:   %s\n', f_featq);
fprintf('Report:       %s\n\n', f_report);

disp(OverviewT);
fprintf('\n[DECISION SUMMARY]\n');
disp(struct2table(Decision, 'AsArray', true));

end

%% ========================================================================
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
if isempty(meta.label), error('Condition label column not found.'); end
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

function [featVars, featCh] = detect_features(DS, meta, cfg)
vars = DS.Properties.VariableNames;
low  = lower(vars);
n = height(DS);

exclude = false(1,numel(vars));
must = {meta.subject, meta.run, meta.phase, meta.label, ...
    'ycondition','condition','ycorrect','correct','trial','trialnum','patternid'};

for i=1:numel(vars)
    for j=1:numel(must)
        if strcmpi(vars{i}, must{j})
            exclude(i) = true;
        end
    end

    if cfg.excludeBehavior
        pats = {'rt','reaction','correct','accuracy','error','resp','response','missing','ismissing'};
        for p=1:numel(pats)
            if contains(low{i}, pats{p})
                exclude(i) = true;
            end
        end
    end
end

featVars0 = {};
featCh0 = [];

for i=1:numel(vars)
    if exclude(i), continue; end
    x = DS.(vars{i});
    if ~(isnumeric(x) && isvector(x) && numel(x)==n), continue; end

    ch = get_ch(vars{i});
    if isnan(ch), continue; end

    featVars0{end+1,1} = vars{i}; %#ok<AGROW>
    featCh0(end+1,1) = ch; %#ok<AGROW>
end

featVars = featVars0;
featCh = featCh0;
end

function ch = get_ch(name)
ch = NaN;
tok = regexp(name, 'ch[_-]?(\d{1,3})', 'tokens', 'once', 'ignorecase');
if ~isempty(tok)
    ch = str2double(tok{1});
    if ch<1 || ch>128
        ch = NaN;
    end
end
end

%% ================= OVERVIEW =================
function OverviewT = make_overview(DS, cfg, meta, isTarget, yAll, runAll, phaseAll, subjAll, featVars)
rows = {};
i = 0;

i=i+1; rows(i,:) = {'Dataset path', cfg.DSpath};
i=i+1; rows(i,:) = {'Target subject', cfg.targetSubject};
i=i+1; rows(i,:) = {'Total dataset rows', height(DS)};
i=i+1; rows(i,:) = {'Target subject rows', sum(isTarget)};
i=i+1; rows(i,:) = {'Other subject rows', sum(~isTarget)};
i=i+1; rows(i,:) = {'Total subjects', numel(unique(subjAll))};
i=i+1; rows(i,:) = {'Runs in target', mat2str(unique(runAll(isTarget))')};
i=i+1; rows(i,:) = {'Phases in target', strjoin(unique(phaseAll(isTarget),'stable'), ', ')};
i=i+1; rows(i,:) = {'Condition codes in target', mat2str(unique(yAll(isTarget))')};
i=i+1; rows(i,:) = {'Condition 1 rows in target', sum(isTarget & yAll==1)};
i=i+1; rows(i,:) = {'Condition 2 rows in target', sum(isTarget & yAll==2)};
i=i+1; rows(i,:) = {'Condition 3 rows in target', sum(isTarget & yAll==3)};
i=i+1; rows(i,:) = {'EEG channel-coded features', numel(featVars)};
i=i+1; rows(i,:) = {'Subject column', meta.subject};
i=i+1; rows(i,:) = {'Run column', meta.run};
i=i+1; rows(i,:) = {'Phase column', meta.phase};
i=i+1; rows(i,:) = {'Condition column', meta.label};

OverviewT = cell2table(rows, 'VariableNames', {'Item','Value'});
end

%% ================= COUNTS =================
function CountT = make_run_phase_condition_counts(DS, cfg, isTarget, yAll, runAll, phaseAll, trialCol)
idx = isTarget;
subj = repmat({cfg.targetSubject}, sum(idx), 1);
runv = runAll(idx);
phasev = phaseAll(idx);
condv = yAll(idx);

if isempty(trialCol)
    trialv = (1:sum(idx))';
else
    trialvRaw = DS.(trialCol);
    trialv = trialvRaw(idx);
end

[G, r, p, c] = findgroups(runv, phasev, condv);
Nrows = splitapply(@numel, condv, G);

if isnumeric(trialv)
    NuniqueTrialNum = splitapply(@(x) numel(unique(x(isfinite(x)))), trialv, G);
    MinTrialNum = splitapply(@(x) min_or_nan(x), trialv, G);
    MaxTrialNum = splitapply(@(x) max_or_nan(x), trialv, G);
else
    trialStr = cellstr(string(trialv));
    NuniqueTrialNum = splitapply(@(x) numel(unique(x)), trialStr, G);
    MinTrialNum = NaN(size(Nrows));
    MaxTrialNum = NaN(size(Nrows));
end

RowsPerUniqueTrial = Nrows ./ max(NuniqueTrialNum,1);
PotentialTrialNumReuse = Nrows > NuniqueTrialNum;

CountT = table(repmat({cfg.targetSubject}, numel(Nrows),1), r, p, c, Nrows, ...
    NuniqueTrialNum, RowsPerUniqueTrial, MinTrialNum, MaxTrialNum, PotentialTrialNumReuse, ...
    'VariableNames', {'Subject','Run','Phase','Condition','Nrows','NuniqueTrialNum', ...
    'RowsPerUniqueTrial','MinTrialNum','MaxTrialNum','PotentialTrialNumReuse'});

CountT = sortrows(CountT, {'Run','Phase','Condition'});
end

function y = min_or_nan(x)
x = x(isfinite(x));
if isempty(x), y = NaN; else, y = min(x); end
end

function y = max_or_nan(x)
x = x(isfinite(x));
if isempty(x), y = NaN; else, y = max(x); end
end

%% ================= DUPLICATES =================
function [DupGroupT, DupSimilarityT] = audit_duplicates(DS, cfg, isTarget, yAll, runAll, phaseAll, trialCol, featVars)
idxTarget = find(isTarget);
nTarget = numel(idxTarget);

if isempty(trialCol)
    trialv = (1:nTarget)';
else
    tvRaw = DS.(trialCol);
    trialv = tvRaw(isTarget);
end

runv = runAll(isTarget);
phasev = phaseAll(isTarget);
condv = yAll(isTarget);

if isnumeric(trialv)
    trialKey = cellstr(string(trialv));
else
    trialKey = cellstr(string(trialv));
end

[G, r, p, c, tr] = findgroups(runv, phasev, condv, trialKey);
Nrows = splitapply(@numel, condv, G);
FirstDatasetRow = splitapply(@(x) x(1), idxTarget, G);

DupGroupT = table(repmat({cfg.targetSubject}, numel(Nrows),1), r, p, c, tr, Nrows, FirstDatasetRow, ...
    Nrows>1, ...
    'VariableNames', {'Subject','Run','Phase','Condition','TrialKey','Nrows','FirstDatasetRow','IsDuplicateGroup'});
DupGroupT = sortrows(DupGroupT, {'Run','Phase','Condition','TrialKey'});

dupGroups = find(Nrows > 1);
simRows = {};
si = 0;

if isempty(dupGroups)
    DupSimilarityT = table();
    return;
end

% Use EEG feature matrix for exact/nearly duplicate check
Xall = double(table2array(DS(isTarget, featVars)));

for dg = dupGroups(:)'
    localRows = find(G == dg);
    datasetRows = idxTarget(localRows);

    X = Xall(localRows, :);

    % Pairwise comparisons among duplicate rows
    pairCount = 0;
    maxAbsDiffList = [];
    meanAbsDiffList = [];
    corrList = [];
    identicalList = [];
    nearlyIdenticalList = [];

    for a=1:size(X,1)-1
        for b=a+1:size(X,1)
            xa = X(a,:);
            xb = X(b,:);
            finite = isfinite(xa) & isfinite(xb);
            if sum(finite) < 10
                maxAbsDiff = NaN;
                meanAbsDiff = NaN;
                rr = NaN;
            else
                d = abs(xa(finite)-xb(finite));
                maxAbsDiff = max(d);
                meanAbsDiff = mean(d);
                rr = corr_safe(xa(finite)', xb(finite)');
            end

            isIdent = isfinite(maxAbsDiff) && maxAbsDiff <= cfg.identicalTolerance;
            isNear = isfinite(rr) && rr >= cfg.nearlyIdenticalCorrelation;

            pairCount = pairCount + 1;
            maxAbsDiffList(end+1) = maxAbsDiff; %#ok<AGROW>
            meanAbsDiffList(end+1) = meanAbsDiff; %#ok<AGROW>
            corrList(end+1) = rr; %#ok<AGROW>
            identicalList(end+1) = isIdent; %#ok<AGROW>
            nearlyIdenticalList(end+1) = isNear; %#ok<AGROW>
        end
    end

    si = si + 1;
    simRows(si,:) = {cfg.targetSubject, r(dg), p{dg}, c(dg), tr{dg}, Nrows(dg), ...
        strjoin(cellstr(string(datasetRows)), ','), pairCount, ...
        nanmax_local(maxAbsDiffList), nanmean_local(meanAbsDiffList), nanmean_local(corrList), ...
        sum(identicalList), sum(nearlyIdenticalList), ...
        mean(identicalList), mean(nearlyIdenticalList)};
end

DupSimilarityT = cell2table(simRows, 'VariableNames', { ...
    'Subject','Run','Phase','Condition','TrialKey','Nrows','DatasetRows','Npairs', ...
    'MaxAbsDiffAcrossPairs','MeanAbsDiffAcrossPairs','MeanCorrAcrossPairs', ...
    'NidenticalPairs','NnearlyIdenticalPairs','FracIdenticalPairs','FracNearlyIdenticalPairs'});

DupSimilarityT = sortrows(DupSimilarityT, {'Run','Phase','Condition','TrialKey'});
end

function r = corr_safe(x,y)
if numel(x) < 3 || std(x) < eps || std(y) < eps
    r = NaN;
    return;
end
C = corrcoef(x,y);
r = C(1,2);
end

function y = nanmax_local(x)
x = x(isfinite(x));
if isempty(x), y = NaN; else, y = max(x); end
end

function y = nanmean_local(x)
x = x(isfinite(x));
if isempty(x), y = NaN; else, y = mean(x); end
end

%% ================= MISSINGNESS =================
function MissingT = compute_missingness_by_subject_run_phase(DS, yAll, runAll, phaseAll, subjAll, featVars)
X = double(table2array(DS(:, featVars)));
isBad = ~isfinite(X);

[G, s, r, p] = findgroups(subjAll, runAll, phaseAll);
Nrows = splitapply(@numel, yAll, G);

FeatureNaNRate_allCells = splitapply(@(idx) mean(isBad(idx,:),'all'), (1:height(DS))', G);
FractionRows_anyFeatureNaN = splitapply(@(idx) mean(any(isBad(idx,:),2)), (1:height(DS))', G);
FractionRows_allFeatureFinite = splitapply(@(idx) mean(all(~isBad(idx,:),2)), (1:height(DS))', G);
FractionRows_allFeatureNaN = splitapply(@(idx) mean(all(isBad(idx,:),2)), (1:height(DS))', G);

MissingT = table(s, r, p, Nrows, FeatureNaNRate_allCells, ...
    FractionRows_anyFeatureNaN, FractionRows_allFeatureFinite, FractionRows_allFeatureNaN, ...
    'VariableNames', {'Subject','Run','Phase','Nrows','FeatureNaNRate_allCells', ...
    'FractionRows_anyFeatureNaN','FractionRows_allFeatureFinite','FractionRows_allFeatureNaN'});

MissingT = sortrows(MissingT, {'Subject','Run','Phase'});
end

function SubjectMissingT = summarize_subject_missingness(MissingT)
[G, s] = findgroups(MissingT.Subject);
Ngroups = splitapply(@numel, MissingT.Nrows, G);
TotalRows = splitapply(@sum, MissingT.Nrows, G);
MeanFeatureNaNRate = splitapply(@(x) mean(x,'omitnan'), MissingT.FeatureNaNRate_allCells, G);
MaxFeatureNaNRate = splitapply(@(x) max(x), MissingT.FeatureNaNRate_allCells, G);
MeanAnyFeatureNaNRows = splitapply(@(x) mean(x,'omitnan'), MissingT.FractionRows_anyFeatureNaN, G);
MinAllFeatureFiniteRows = splitapply(@(x) min(x), MissingT.FractionRows_allFeatureFinite, G);

SubjectMissingT = table(s, Ngroups, TotalRows, MeanFeatureNaNRate, MaxFeatureNaNRate, ...
    MeanAnyFeatureNaNRows, MinAllFeatureFiniteRows, ...
    'VariableNames', {'Subject','Ngroups','TotalRows','MeanFeatureNaNRate','MaxFeatureNaNRate', ...
    'MeanAnyFeatureNaNRows','MinAllFeatureFiniteRows'});

SubjectMissingT = sortrows(SubjectMissingT, 'MeanFeatureNaNRate', 'descend');
end

%% ================= FEATURE QUALITY =================
function FeatureQualityT = compute_target_feature_quality(DS, cfg, isTarget, featVars, featCh)
X = double(table2array(DS(isTarget, featVars)));
nFeat = numel(featVars);

NaNFrac = mean(isnan(X),1)';
InfFrac = mean(isinf(X),1)';
FiniteFrac = mean(isfinite(X),1)';
MeanVal = NaN(nFeat,1);
MedianVal = NaN(nFeat,1);
StdVal = NaN(nFeat,1);
MADVal = NaN(nFeat,1);
MinVal = NaN(nFeat,1);
MaxVal = NaN(nFeat,1);
ExtremeFrac = NaN(nFeat,1);

for j=1:nFeat
    x = X(:,j);
    xf = x(isfinite(x));
    if isempty(xf), continue; end

    MeanVal(j) = mean(xf);
    MedianVal(j) = median(xf);
    StdVal(j) = std(xf);
    MADVal(j) = median(abs(xf - median(xf)));
    MinVal(j) = min(xf);
    MaxVal(j) = max(xf);

    rz = robust_z_vec(xf);
    ExtremeFrac(j) = mean(abs(rz) > cfg.extremeAbsRobustZ);
end

Flag_HighNaN = (NaNFrac + InfFrac) > cfg.highFeatureNaNRate_feature;
Flag_Constant = ~isfinite(StdVal) | StdVal < cfg.constantStdThr;
Flag_ExtremeFrac = ExtremeFrac > cfg.highExtremeFracThr;

FeatureQualityT = table(featVars(:), featCh(:), NaNFrac, InfFrac, FiniteFrac, ...
    MeanVal, MedianVal, StdVal, MADVal, MinVal, MaxVal, ExtremeFrac, ...
    Flag_HighNaN, Flag_Constant, Flag_ExtremeFrac, ...
    'VariableNames', {'FeatureName','ChannelIndex','NaNFrac','InfFrac','FiniteFrac', ...
    'Mean','Median','Std','MAD','Min','Max','ExtremeFrac', ...
    'Flag_HighNaN','Flag_Constant','Flag_ExtremeFrac'});

FeatureQualityT = sortrows(FeatureQualityT, {'Flag_HighNaN','Flag_Constant','Flag_ExtremeFrac','ExtremeFrac'}, ...
    {'descend','descend','descend','descend'});
end

function z = robust_z_vec(x)
m = median(x);
mad0 = median(abs(x-m));
if ~isfinite(mad0) || mad0 < eps
    s = std(x);
    if ~isfinite(s) || s < eps, s = 1; end
    z = (x-m)/s;
else
    z = (x-m)/(1.4826*mad0);
end
end

%% ================= DECISION =================
function Decision = make_decision_struct(cfg, OverviewT, CountT, DupGroupT, DupSimilarityT, MissingT, SubjectMissingT, FeatureQualityT)

targetMiss = SubjectMissingT(strcmp(SubjectMissingT.Subject, cfg.targetSubject), :);

nTargetRows = get_overview_number(OverviewT, 'Target subject rows');
nDupGroups = sum(DupGroupT.IsDuplicateGroup);
nDupRows = sum(DupGroupT.Nrows(DupGroupT.IsDuplicateGroup));
meanRowsPerUnique = mean(CountT.RowsPerUniqueTrial, 'omitnan');
maxRowsPerUnique = max(CountT.RowsPerUniqueTrial);

if isempty(DupSimilarityT)
    fracIdenticalPairs = 0;
    fracNearPairs = 0;
    meanCorr = NaN;
else
    fracIdenticalPairs = mean(DupSimilarityT.FracIdenticalPairs, 'omitnan');
    fracNearPairs = mean(DupSimilarityT.FracNearlyIdenticalPairs, 'omitnan');
    meanCorr = mean(DupSimilarityT.MeanCorrAcrossPairs, 'omitnan');
end

meanMiss = targetMiss.MeanFeatureNaNRate;
maxMiss = targetMiss.MaxFeatureNaNRate;
rankMiss = find(strcmp(SubjectMissingT.Subject, cfg.targetSubject), 1);

nHighNaNFeatures = sum(FeatureQualityT.Flag_HighNaN);
nConstantFeatures = sum(FeatureQualityT.Flag_Constant);
nExtremeFeatures = sum(FeatureQualityT.Flag_ExtremeFrac);

% Decision logic
likelyTrueDuplicate = (nDupGroups > 0) && (fracIdenticalPairs > 0.80 || fracNearPairs > 0.90);
likelyTrialNumReuse = (nDupGroups > 0) && (meanRowsPerUnique >= 1.5) && ~likelyTrueDuplicate;
missingnessConcern = meanMiss > cfg.highFeatureNaNRate_subject || maxMiss > cfg.highFeatureNaNRate_group;
featureQualityConcern = nHighNaNFeatures > 0 || nConstantFeatures > 0 || nExtremeFeatures > 0;

if likelyTrueDuplicate
    recommendation = 'Strong QC concern: duplicated rows appear identical/nearly identical. Prefer correcting duplicates or excluding duplicated entries before excluding entire subject.';
elseif missingnessConcern && likelyTrialNumReuse
    recommendation = 'QC candidate: TrialNum reuse plus high missingness. Audit source sessions/workspaces; consider excluding s21 only if source error or missingness criterion is confirmed.';
elseif missingnessConcern
    recommendation = 'QC candidate due to high feature missingness. Run sensitivity excluding s21 and define an independent missingness threshold.';
elseif likelyTrialNumReuse
    recommendation = 'Likely repeated session / TrialNum reset rather than duplicate. Do not exclude automatically; add session/block ID or verify source workspaces.';
else
    recommendation = 'No strong independent reason to exclude s21 from this audit alone. Keep in main analysis and optionally report sensitivity.';
end

Decision = struct();
Decision.TargetSubject = cfg.targetSubject;
Decision.TargetRows = nTargetRows;
Decision.DuplicateGroups = nDupGroups;
Decision.RowsInDuplicateGroups = nDupRows;
Decision.MeanRowsPerUniqueTrial = meanRowsPerUnique;
Decision.MaxRowsPerUniqueTrial = maxRowsPerUnique;
Decision.FractionIdenticalDuplicatePairs = fracIdenticalPairs;
Decision.FractionNearlyIdenticalDuplicatePairs = fracNearPairs;
Decision.MeanDuplicatePairCorrelation = meanCorr;
Decision.MeanFeatureNaNRate = meanMiss;
Decision.MaxFeatureNaNRate = maxMiss;
Decision.MissingnessRankAmongSubjects = rankMiss;
Decision.NHighNaNFeatures = nHighNaNFeatures;
Decision.NConstantFeatures = nConstantFeatures;
Decision.NExtremeFeatures = nExtremeFeatures;
Decision.LikelyTrueDuplicate = likelyTrueDuplicate;
Decision.LikelyTrialNumReuse = likelyTrialNumReuse;
Decision.MissingnessConcern = missingnessConcern;
Decision.FeatureQualityConcern = featureQualityConcern;
Decision.Recommendation = recommendation;
end

function val = get_overview_number(OverviewT, item)
idx = strcmp(OverviewT.Item, item);
if ~any(idx), val = NaN; return; end
val = str2double(string(OverviewT.Value(idx)));
end

%% ================= REPORT =================
function write_report(reportPath, cfg, OverviewT, CountT, DupGroupT, DupSimilarityT, ...
    MissingT, SubjectMissingT, FeatureQualityT, Decision)

fid = fopen(reportPath, 'w');

fprintf(fid, 'STEP18 s21 DUPLICATE / MISSINGNESS AUDIT REPORT\n');
fprintf(fid, '==============================================\n\n');

fprintf(fid, 'Target subject: %s\n', cfg.targetSubject);
fprintf(fid, 'Dataset: %s\n\n', cfg.DSpath);

fprintf(fid, 'Overview:\n');
for i=1:height(OverviewT)
    fprintf(fid, '  %s: %s\n', OverviewT.Item{i}, string(OverviewT.Value(i)));
end

fprintf(fid, '\nRun×Phase×Condition counts:\n');
for i=1:height(CountT)
    fprintf(fid, '  run%d | %s | cond%d | rows=%d | uniqueTrialNum=%d | rows/unique=%.2f | reuse=%d\n', ...
        CountT.Run(i), CountT.Phase{i}, CountT.Condition(i), CountT.Nrows(i), ...
        CountT.NuniqueTrialNum(i), CountT.RowsPerUniqueTrial(i), CountT.PotentialTrialNumReuse(i));
end

fprintf(fid, '\nDuplicate groups:\n');
fprintf(fid, '  Total duplicate groups: %d\n', sum(DupGroupT.IsDuplicateGroup));
fprintf(fid, '  Rows in duplicate groups: %d\n', sum(DupGroupT.Nrows(DupGroupT.IsDuplicateGroup)));
if ~isempty(DupSimilarityT)
    fprintf(fid, '  Mean fraction identical pairs: %.4f\n', mean(DupSimilarityT.FracIdenticalPairs,'omitnan'));
    fprintf(fid, '  Mean fraction nearly identical pairs: %.4f\n', mean(DupSimilarityT.FracNearlyIdenticalPairs,'omitnan'));
    fprintf(fid, '  Mean duplicate-pair correlation: %.6f\n', mean(DupSimilarityT.MeanCorrAcrossPairs,'omitnan'));
else
    fprintf(fid, '  No duplicate groups detected.\n');
end

targetMiss = SubjectMissingT(strcmp(SubjectMissingT.Subject, cfg.targetSubject), :);
fprintf(fid, '\nMissingness:\n');
if ~isempty(targetMiss)
    fprintf(fid, '  %s MeanFeatureNaNRate: %.4f\n', cfg.targetSubject, targetMiss.MeanFeatureNaNRate);
    fprintf(fid, '  %s MaxFeatureNaNRate: %.4f\n', cfg.targetSubject, targetMiss.MaxFeatureNaNRate);
    fprintf(fid, '  Missingness rank among subjects: %d of %d (1 = highest missingness)\n', ...
        Decision.MissingnessRankAmongSubjects, height(SubjectMissingT));
end

fprintf(fid, '\nTop subjects by missingness:\n');
n = min(10, height(SubjectMissingT));
for i=1:n
    fprintf(fid, '  %s | MeanNaN=%.4f | MaxNaN=%.4f | TotalRows=%d\n', ...
        SubjectMissingT.Subject{i}, SubjectMissingT.MeanFeatureNaNRate(i), ...
        SubjectMissingT.MaxFeatureNaNRate(i), SubjectMissingT.TotalRows(i));
end

fprintf(fid, '\nFeature-level warnings for %s:\n', cfg.targetSubject);
fprintf(fid, '  High NaN/Inf features: %d\n', sum(FeatureQualityT.Flag_HighNaN));
fprintf(fid, '  Constant features: %d\n', sum(FeatureQualityT.Flag_Constant));
fprintf(fid, '  Extreme-fraction features: %d\n', sum(FeatureQualityT.Flag_ExtremeFrac));

warn = FeatureQualityT(FeatureQualityT.Flag_HighNaN | FeatureQualityT.Flag_Constant | FeatureQualityT.Flag_ExtremeFrac, :);
if ~isempty(warn)
    fprintf(fid, '\nTop warning features:\n');
    n = min(30, height(warn));
    for i=1:n
        fprintf(fid, '  %s | ch%d | NaN=%.4f | Std=%.4g | ExtremeFrac=%.4f | HighNaN=%d | Const=%d | Extreme=%d\n', ...
            warn.FeatureName{i}, warn.ChannelIndex(i), warn.NaNFrac(i)+warn.InfFrac(i), ...
            warn.Std(i), warn.ExtremeFrac(i), warn.Flag_HighNaN(i), ...
            warn.Flag_Constant(i), warn.Flag_ExtremeFrac(i));
    end
end

fprintf(fid, '\nDecision flags:\n');
fprintf(fid, '  LikelyTrueDuplicate: %d\n', Decision.LikelyTrueDuplicate);
fprintf(fid, '  LikelyTrialNumReuse: %d\n', Decision.LikelyTrialNumReuse);
fprintf(fid, '  MissingnessConcern: %d\n', Decision.MissingnessConcern);
fprintf(fid, '  FeatureQualityConcern: %d\n', Decision.FeatureQualityConcern);
fprintf(fid, '\nRecommendation:\n  %s\n', Decision.Recommendation);

fprintf(fid, '\nHow to use this in the manuscript:\n');
fprintf(fid, '  - If true duplicates are confirmed, correct duplicated rows or exclude duplicated entries using an independent QC rule.\n');
fprintf(fid, '  - If TrialNum reuse reflects extra sessions, do not delete s21 only for having 200 rows; document extra sessions or add session ID.\n');
fprintf(fid, '  - If high missingness is confirmed, define a subject-level missingness threshold and report an exclusion/sensitivity analysis.\n');
fprintf(fid, '  - Never exclude s21 solely because model accuracy improves; exclusion must be justified by QC independent of labels/performance.\n');

fclose(fid);
end

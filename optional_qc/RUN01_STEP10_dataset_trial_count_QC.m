function RUN01_STEP10_dataset_trial_count_QC()

cfg = get_project_config();

clc;

%% ================= CONFIG =================
cfg.ROOT = fullfile(cfg.dataRoot, 'Subjects');

dsCandidates = { ...
    fullfile(cfg.ROOT, '_wm_ml', 'dataset.mat'), ...
    fullfile(cfg.dataRoot, 'Features', 'Subjects', '_wm_ml', 'dataset.mat')};

cfg.DSpath = '';
for di = 1:numel(dsCandidates)
    if exist(dsCandidates{di}, 'file')
        cfg.DSpath = dsCandidates{di};
        break;
    end
end

if isempty(cfg.DSpath)
    error('Dataset not found. Checked:\n  %s\n  %s', dsCandidates{1}, dsCandidates{2});
end

cfg.outBaseDir = cfg.outputRoot;
cfg.outDir = fullfile(cfg.outBaseDir, 'STEP10_dataset_trial_count_QC');

if ~exist(cfg.outBaseDir,'dir'), mkdir(cfg.outBaseDir); end
if ~exist(cfg.outDir,'dir'), mkdir(cfg.outDir); end

% Raw-onset flagged runs from final audit
cfg.flagged.Subject = {'s7','s12','s14','s15','s17','s19','s20'}';
cfg.flagged.Run     = [  1,   1,    1,    1,    1,    1,    1]';

fprintf('\n[STEP10] DATASET / TRIAL-COUNT QC\n');
fprintf('Dataset: %s\n', cfg.DSpath);
fprintf('Output:  %s\n', cfg.outDir);

%% ================= LOAD DATASET =================
S = load(cfg.DSpath);
DS = get_table_from_mat(S);
meta = detect_meta(DS);

[yAll, runAll, phaseAll, subjAll] = get_core_vectors(DS, meta);

trialCol = detect_trial_col(DS);
correctCol = detect_correct_col(DS);

[featVars, ~] = detect_features(DS, meta);
fprintf('\n[DATASET]\n');
fprintf('Rows=%d | Vars=%d | EEG channel-coded features=%d\n', height(DS), width(DS), numel(featVars));
fprintf('Meta: subject=%s | run=%s | phase=%s | label=%s | trial=%s | correct=%s\n', ...
    meta.subject, meta.run, meta.phase, meta.label, trialCol, correctCol);

%% ================= ADD QC FLAGS =================
flagMask = false(height(DS),1);
for i=1:numel(cfg.flagged.Subject)
    flagMask = flagMask | (strcmp(subjAll, cfg.flagged.Subject{i}) & runAll == cfg.flagged.Run(i));
end

%% ================= TABLE 1: Subject×Run×Phase×Condition =================
[G1, gSubj, gRun, gPhase, gCond] = findgroups(subjAll, runAll, phaseAll, yAll);
Nrows = splitapply(@numel, yAll, G1);

if ~isempty(trialCol)
    trialVals = DS.(trialCol);
    if isnumeric(trialVals)
        NuniqueTrials = splitapply(@(x) numel(unique(x(isfinite(x)))), trialVals, G1);
    else
        trialValsStr = cellstr(string(trialVals));
        NuniqueTrials = splitapply(@(x) numel(unique(x)), trialValsStr, G1);
    end
else
    NuniqueTrials = Nrows;
end

FlaggedRawOnsetRun = splitapply(@(x) any(x), flagMask, G1);

Tcond = table(gSubj, gRun, gPhase, gCond, Nrows, NuniqueTrials, FlaggedRawOnsetRun, ...
    'VariableNames', {'Subject','Run','Phase','Condition','Nrows','NuniqueTrials','FlaggedRawOnsetRun'});
Tcond = sortrows(Tcond, {'Subject','Run','Phase','Condition'});

%% ================= TABLE 2: Subject×Run×Phase totals =================
[G2, s2, r2, p2] = findgroups(subjAll, runAll, phaseAll);
Nrows2 = splitapply(@numel, yAll, G2);
Ncond1 = splitapply(@(x) sum(x==1), yAll, G2);
Ncond2 = splitapply(@(x) sum(x==2), yAll, G2);
Ncond3 = splitapply(@(x) sum(x==3), yAll, G2);

if ~isempty(trialCol)
    trialVals = DS.(trialCol);
    if isnumeric(trialVals)
        NuniqueTrials2 = splitapply(@(x) numel(unique(x(isfinite(x)))), trialVals, G2);
    else
        trialValsStr = cellstr(string(trialVals));
        NuniqueTrials2 = splitapply(@(x) numel(unique(x)), trialValsStr, G2);
    end
else
    NuniqueTrials2 = Nrows2;
end

Flagged2 = splitapply(@(x) any(x), flagMask, G2);

if ~isempty(correctCol) && isnumeric(DS.(correctCol))
    c = double(DS.(correctCol));
    NcorrectKnown = splitapply(@(x) sum(isfinite(x)), c, G2);
    Ncorrect1 = splitapply(@(x) sum(x==1), c, G2);
    Ncorrect0 = splitapply(@(x) sum(x==0), c, G2);
    CorrectRate = Ncorrect1 ./ max(NcorrectKnown,1);
else
    NcorrectKnown = NaN(size(Nrows2));
    Ncorrect1 = NaN(size(Nrows2));
    Ncorrect0 = NaN(size(Nrows2));
    CorrectRate = NaN(size(Nrows2));
end

Tphase = table(s2, r2, p2, Nrows2, NuniqueTrials2, Ncond1, Ncond2, Ncond3, ...
    NcorrectKnown, Ncorrect1, Ncorrect0, CorrectRate, Flagged2, ...
    'VariableNames', {'Subject','Run','Phase','Nrows','NuniqueTrials','Ncond1','Ncond2','Ncond3', ...
    'NcorrectKnown','Ncorrect1','Ncorrect0','CorrectRate','FlaggedRawOnsetRun'});
Tphase = sortrows(Tphase, {'Subject','Run','Phase'});

%% ================= TABLE 3: Wide Subject×Run =================
subjects = unique(subjAll, 'stable');
runs = unique(runAll(isfinite(runAll)));
phases = {'stim','maint','retr'};

wideRows = {};
wi = 0;
for si = 1:numel(subjects)
    for ri = 1:numel(runs)
        ss = subjects{si};
        rr = runs(ri);
        idxSR = strcmp(subjAll, ss) & runAll==rr;
        if ~any(idxSR), continue; end

        wi = wi + 1;
        vals = {};
        vals{end+1} = ss;
        vals{end+1} = rr;
        vals{end+1} = any(flagMask(idxSR));

        for pp=1:numel(phases)
            idx = idxSR & strcmp(phaseAll, phases{pp});
            vals{end+1} = sum(idx);
            vals{end+1} = sum(idx & yAll==1);
            vals{end+1} = sum(idx & yAll==2);
            vals{end+1} = sum(idx & yAll==3);
        end
        wideRows(wi,:) = vals; %#ok<AGROW>
    end
end

wideNames = {'Subject','Run','FlaggedRawOnsetRun'};
for pp=1:numel(phases)
    ph = phases{pp};
    wideNames{end+1} = [ph '_N']; %#ok<AGROW>
    wideNames{end+1} = [ph '_cond1']; %#ok<AGROW>
    wideNames{end+1} = [ph '_cond2']; %#ok<AGROW>
    wideNames{end+1} = [ph '_cond3']; %#ok<AGROW>
end
Twide = cell2table(wideRows, 'VariableNames', wideNames);
Twide = sortrows(Twide, {'Subject','Run'});

%% ================= TABLE 4: Feature missingness =================
missRows = {};
mi = 0;
if ~isempty(featVars)
    Xall = table2array(DS(:, featVars));
    for gi=1:max(G2)
        idx = G2==gi;
        Xi = Xall(idx,:);
        mi = mi + 1;
        missRows(mi,:) = {s2{gi}, r2(gi), p2{gi}, sum(idx), ...
            mean(~isfinite(Xi(:))), ...
            mean(all(~isfinite(Xi),2)), ...
            mean(all(isfinite(Xi),2)), ...
            Flagged2(gi)};
    end
    Tmiss = cell2table(missRows, 'VariableNames', ...
        {'Subject','Run','Phase','Nrows','FeatureNaNRate_allCells','FractionRows_allFeatureNaN','FractionRows_allFeatureFinite','FlaggedRawOnsetRun'});
else
    Tmiss = table();
end

%% ================= TABLE 5: Subject overall =================
[Gs, ss0] = findgroups(subjAll);
Ntotal = splitapply(@numel, yAll, Gs);
Nrun1 = splitapply(@(r) sum(r==1), runAll, Gs);
Nrun2 = splitapply(@(r) sum(r==2), runAll, Gs);
Nrun3 = splitapply(@(r) sum(r==3), runAll, Gs);
Nflag = splitapply(@(x) sum(x), flagMask, Gs);
Ncond1tot = splitapply(@(x) sum(x==1), yAll, Gs);
Ncond2tot = splitapply(@(x) sum(x==2), yAll, Gs);
Ncond3tot = splitapply(@(x) sum(x==3), yAll, Gs);

Tsubj = table(ss0, Ntotal, Nrun1, Nrun2, Nrun3, Ncond1tot, Ncond2tot, Ncond3tot, Nflag, ...
    'VariableNames', {'Subject','NtotalRows','Nrun1Rows','Nrun2Rows','Nrun3Rows','Ncond1','Ncond2','Ncond3','NflaggedRows'});
Tsubj = sortrows(Tsubj, 'Subject');

%% ================= SAVE =================
f1 = fullfile(cfg.outDir, 'STEP10_trials_by_subject_run_phase_condition.csv');
f2 = fullfile(cfg.outDir, 'STEP10_trials_by_subject_run_phase.csv');
f3 = fullfile(cfg.outDir, 'STEP10_trials_wide_subject_run.csv');
f4 = fullfile(cfg.outDir, 'STEP10_feature_missingness_by_subject_run_phase.csv');
f5 = fullfile(cfg.outDir, 'STEP10_subject_overall_summary.csv');
matPath = fullfile(cfg.outDir, 'STEP10_dataset_trial_count_QC_results.mat');

writetable(Tcond, f1);
writetable(Tphase, f2);
writetable(Twide, f3);
if ~isempty(Tmiss), writetable(Tmiss, f4); end
writetable(Tsubj, f5);
save(matPath, 'cfg', 'Tcond', 'Tphase', 'Twide', 'Tmiss', 'Tsubj', 'meta');

%% ================= REPORT =================
reportPath = fullfile(cfg.outDir, 'STEP10_QC_report.txt');
fid = fopen(reportPath, 'w');
fprintf(fid, 'STEP10 DATASET / TRIAL-COUNT QC REPORT\n');
fprintf(fid, '=====================================\n\n');
fprintf(fid, 'Dataset: %s\n', cfg.DSpath);
fprintf(fid, 'Rows: %d\n', height(DS));
fprintf(fid, 'Variables: %d\n', width(DS));
fprintf(fid, 'Subjects: %d\n', numel(unique(subjAll)));
fprintf(fid, 'Runs: %s\n', mat2str(unique(runAll(isfinite(runAll)))'));
fprintf(fid, 'Phases: %s\n', strjoin(unique(phaseAll,'stable'), ', '));
fprintf(fid, 'EEG channel-coded features: %d\n\n', numel(featVars));
fprintf(fid, 'Raw-onset flagged runs:\n');
for i=1:numel(cfg.flagged.Subject)
    fprintf(fid, '  %s run%d\n', cfg.flagged.Subject{i}, cfg.flagged.Run(i));
end
fprintf(fid, 'Rows belonging to flagged runs: %d / %d\n\n', sum(flagMask), height(DS));
fprintf(fid, 'Condition counts overall:\n');
fprintf(fid, '  Condition 1: %d\n', sum(yAll==1));
fprintf(fid, '  Condition 2: %d\n', sum(yAll==2));
fprintf(fid, '  Condition 3: %d\n\n', sum(yAll==3));
fprintf(fid, 'Minimum rows per Subject×Run×Phase group: %d\n', min(Tphase.Nrows));
fprintf(fid, 'Maximum rows per Subject×Run×Phase group: %d\n', max(Tphase.Nrows));
fprintf(fid, 'Median rows per Subject×Run×Phase group: %.1f\n\n', median(Tphase.Nrows));
fprintf(fid, 'Files saved:\n');
fprintf(fid, '  %s\n', f1);
fprintf(fid, '  %s\n', f2);
fprintf(fid, '  %s\n', f3);
fprintf(fid, '  %s\n', f4);
fprintf(fid, '  %s\n', f5);
fprintf(fid, '  %s\n', matPath);
fclose(fid);

fprintf('\n[DONE STEP10]\n');
fprintf('Output: %s\n', cfg.outDir);
fprintf('Report: %s\n\n', reportPath);
disp(Tphase(1:min(30,height(Tphase)),:));

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

function [y, runv, phasev, subjv] = get_core_vectors(DS, meta)
y = double(DS.(meta.label));
runv = double(DS.(meta.run));
phasev = lower(cellstr(string(DS.(meta.phase))));
subjv = cellstr(string(DS.(meta.subject)));
end

function trialCol = detect_trial_col(DS)
v = DS.Properties.VariableNames;
l = lower(v);
trialCol = pick(v,l,{'trialnum','trial','trialindex','trial_idx','trialid'});
end

function correctCol = detect_correct_col(DS)
v = DS.Properties.VariableNames;
l = lower(v);
correctCol = pick(v,l,{'ycorrect','correct','iscorrect','accuracy'});
end

function [featVars, featCh] = detect_features(DS, meta)
vars = DS.Properties.VariableNames;
low  = lower(vars);
n = height(DS);
exclude = false(1,numel(vars));
must = {meta.subject, meta.run, meta.phase, meta.label, ...
    'ycondition','condition','ycorrect','correct','trial','trialnum','patternid','rt','response','resp'};
for i=1:numel(vars)
    for j=1:numel(must)
        if strcmpi(vars{i}, must{j})
            exclude(i) = true;
        end
    end
end
isNum = false(1,numel(vars));
for i=1:numel(vars)
    x = DS.(vars{i});
    isNum(i) = isnumeric(x) && isvector(x) && numel(x)==n;
end
idx = find(isNum & ~exclude);
featVars0 = vars(idx);
featCh0 = NaN(1,numel(featVars0));
for i=1:numel(featVars0)
    featCh0(i) = get_ch(featVars0{i});
end
keep = ~isnan(featCh0);
featVars = featVars0(keep);
featCh = featCh0(keep);
end

function ch = get_ch(name)
ch = NaN;
tok = regexp(name, 'ch[_-]?(\d{1,2})', 'tokens', 'once', 'ignorecase');
if ~isempty(tok)
    ch = str2double(tok{1});
    if ch<1 || ch>64
        ch = NaN;
    end
end
end

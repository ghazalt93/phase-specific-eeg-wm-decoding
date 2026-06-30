function wm_subjectwise_channel_qc()
cfg = get_project_config();


clc; close all;

%% ================= CONFIG =================
cfg.projectDir = fullfile(cfg.dataRoot, 'best');
cfg.outBaseDir = cfg.outputRoot;
cfg.outDir = fullfile(cfg.outBaseDir, 'SUBJECTWISE_CHANNEL_QC');

if ~exist(cfg.outBaseDir, 'dir'), mkdir(cfg.outBaseDir); end
if ~exist(cfg.outDir, 'dir'), mkdir(cfg.outDir); end

candidateDirs = { ...
    fullfile(cfg.projectDir, 'Subjects', 'DERIVED_RESULTS'), ...
    fullfile(cfg.projectDir, 'DERIVED_RESULTS')};

cfg.derivedDir = '';
for i = 1:numel(candidateDirs)
    if exist(candidateDirs{i}, 'dir')
        d = dir(fullfile(candidateDirs{i}, '**', '*__DERIVED.mat'));
        fprintf('[PATH-CHECK] %s | derived files=%d\n', candidateDirs{i}, numel(d));
        if ~isempty(d)
            cfg.derivedDir = candidateDirs{i};
            break;
        end
    end
end

if isempty(cfg.derivedDir)
    error('No folder containing *__DERIVED.mat was found.');
end

cfg.phases = {'Stim','Maint','Retr'};

cfg.maxTrialsPerPhase = 40;
cfg.maxFiles = Inf;

cfg.flatStdRatioThr = 0.08;     
cfg.robustZthr = 4.5;           
cfg.badScoreThr = 2;           

cfg.subjectBadFractionThr = 0.20;

cfg.globalBadSubjectFractionThr = 0.20;

fprintf('\n[SUBJECT-WISE CHANNEL QC]\n');
fprintf('derivedDir: %s\n', cfg.derivedDir);
fprintf('outDir:     %s\n', cfg.outDir);

%% ================= FIND FILES =================
files = dir(fullfile(cfg.derivedDir, '**', '*__DERIVED.mat'));
if isfinite(cfg.maxFiles)
    files = files(1:min(cfg.maxFiles, numel(files)));
end
fprintf('Files to scan: %d\n', numel(files));

%% ================= SESSION/RUN/PHASE QC =================
rows = {};
ri = 0;

for f = 1:numel(files)
    filePath = fullfile(files(f).folder, files(f).name);
    [~, runTag] = fileparts(files(f).name);
    subj = infer_subject(files(f).folder, files(f).name);

    try
        S = load(filePath, 'SEG');
        if ~isfield(S, 'SEG')
            fprintf('[SKIP] no SEG: %s\n', filePath);
            continue;
        end
        SEG = S.SEG;
    catch ME
        fprintf('[SKIP-LOAD] %s | %s\n', filePath, ME.message);
        continue;
    end

    fprintf('\n[%03d/%03d] %s | %s\n', f, numel(files), subj, runTag);

    for p = 1:numel(cfg.phases)
        ph = cfg.phases{p};
        A = get_phase_array(SEG, ph);

        if isempty(A) || ndims(A) ~= 3
            fprintf('  [SKIP] %s missing/not 3D\n', ph);
            continue;
        end

        nCh = size(A,1);
        nTime = size(A,2);
        nTrOrig = size(A,3);

        if isfinite(cfg.maxTrialsPerPhase) && nTrOrig > cfg.maxTrialsPerPhase
            useTr = round(linspace(1, nTrOrig, cfg.maxTrialsPerPhase));
            A = A(:,:,useTr);
        end

        M = compute_basic_qc(A);
        F = flag_bad(M, cfg);

        fprintf('  %s | nCh=%d nTime=%d nTrial=%d | bad=%s\n', ...
            ph, nCh, nTime, nTrOrig, mat2str(find(F.IsBad)'));

        for ch = 1:nCh
            ri = ri + 1;
            rows(ri,:) = { ...
                subj, runTag, ph, ch, nTime, nTrOrig, ...
                M.Std(ch), M.P2P(ch), M.DerivSTD(ch), M.CorrRef(ch), ...
                F.Z_STD(ch), F.Z_P2P(ch), F.Z_Deriv(ch), F.Z_CorrRef(ch), ...
                F.FlagFlat(ch), F.FlagHighSTD(ch), F.FlagHighP2P(ch), F.FlagHighDeriv(ch), F.FlagLowCorr(ch), ...
                F.BadScore(ch), F.IsBad(ch), F.IsSuspicious(ch), filePath};
        end
    end
end

if isempty(rows)
    error('No QC rows created.');
end

Tsession = cell2table(rows, 'VariableNames', { ...
    'Subject','RunTag','Phase','Channel','Ntime','NtrialsOriginal', ...
    'STD','P2P','DerivSTD','CorrRef', ...
    'Z_STD','Z_P2P','Z_DerivSTD','Z_CorrRef', ...
    'FlagFlat','FlagHighSTD','FlagHighP2P','FlagHighDerivSTD','FlagLowCorr', ...
    'BadScore','IsBad','IsSuspicious','SourceFile'});

sessionCsv = fullfile(cfg.outDir, 'subjectwise_channel_qc_session.csv');
writetable(Tsession, sessionCsv);

%% ================= SUBJECT-LEVEL SUMMARY =================
[Gsubch, subjG, chG] = findgroups(Tsession.Subject, Tsession.Channel);

TotalCount = splitapply(@numel, Tsession.IsBad, Gsubch);
BadCount = splitapply(@sum, double(Tsession.IsBad), Gsubch);
SuspiciousCount = splitapply(@sum, double(Tsession.IsSuspicious), Gsubch);

BadFraction = BadCount ./ max(TotalCount,1);
SuspiciousFraction = SuspiciousCount ./ max(TotalCount,1);

Tsubject = table(subjG, chG, TotalCount, BadCount, BadFraction, SuspiciousCount, SuspiciousFraction, ...
    'VariableNames', {'Subject','Channel','TotalCount','BadCount','BadFraction','SuspiciousCount','SuspiciousFraction'});

Tsubject.SubjectBad = Tsubject.BadFraction >= cfg.subjectBadFractionThr;
Tsubject.SubjectGood = ~Tsubject.SubjectBad;

Tsubject = sortrows(Tsubject, {'Subject','BadFraction'}, {'ascend','descend'});
subjectCsv = fullfile(cfg.outDir, 'subjectwise_channel_qc_subject_summary.csv');
writetable(Tsubject, subjectCsv);

%% ================= SUBJECT GOOD/BAD LISTS =================
subjects = unique(Tsubject.Subject, 'stable');

subjRows = {};
si = 0;

for i = 1:numel(subjects)
    s = subjects{i};
    idx = strcmp(Tsubject.Subject, s);

    badCh = Tsubject.Channel(idx & Tsubject.SubjectBad)';
    goodCh = setdiff(1:64, badCh);

    si = si + 1;
    subjRows(si,:) = {s, mat2str(badCh), mat2str(goodCh), numel(badCh), numel(goodCh)};
end

Tlists = cell2table(subjRows, 'VariableNames', ...
    {'Subject','BadChannels','GoodChannels','Nbad','Ngood'});

listsCsv = fullfile(cfg.outDir, 'subjectwise_good_bad_channel_lists.csv');
writetable(Tlists, listsCsv);

%% ================= GLOBAL SUMMARY ACROSS SUBJECTS =================
[Gch, ch2] = findgroups(Tsubject.Channel);
Nsubjects = splitapply(@numel, Tsubject.SubjectBad, Gch);
NsubjectBad = splitapply(@sum, double(Tsubject.SubjectBad), Gch);
SubjectBadFraction = NsubjectBad ./ max(Nsubjects,1);

Tglobal = table(ch2, Nsubjects, NsubjectBad, SubjectBadFraction, ...
    'VariableNames', {'Channel','Nsubjects','NsubjectBad','SubjectBadFraction'});

Tglobal.GlobalBad = Tglobal.SubjectBadFraction >= cfg.globalBadSubjectFractionThr;
Tglobal.GlobalGood = ~Tglobal.GlobalBad;

Tglobal = sortrows(Tglobal, 'SubjectBadFraction', 'descend');
globalCsv = fullfile(cfg.outDir, 'subjectwise_channel_qc_global_summary.csv');
writetable(Tglobal, globalCsv);

globalBadChannels = sort(Tglobal.Channel(Tglobal.GlobalBad))';
globalGoodChannels = setdiff(1:64, globalBadChannels);

matPath = fullfile(cfg.outDir, 'subjectwise_channel_qc_results.mat');
save(matPath, 'cfg', 'Tsession', 'Tsubject', 'Tlists', 'Tglobal', ...
    'globalBadChannels', 'globalGoodChannels');

txtPath = fullfile(cfg.outDir, 'subjectwise_channel_qc_recommendation.txt');
fid = fopen(txtPath, 'w');
fprintf(fid, 'Subject-level bad threshold: BadFraction >= %.3f\n', cfg.subjectBadFractionThr);
fprintf(fid, 'Global bad threshold: SubjectBadFraction >= %.3f\n\n', cfg.globalBadSubjectFractionThr);
fprintf(fid, 'Global bad channels = %s\n', mat2str(globalBadChannels));
fprintf(fid, 'Global good channels = %s\n\n', mat2str(globalGoodChannels));
fprintf(fid, 'Important recommendation:\n');
fprintf(fid, 'For final across-subject modeling, do NOT select different feature dimensions per test subject.\n');
fprintf(fid, 'Best option: interpolate each subject/run bad channels on raw EEG, then re-extract same 64-channel features.\n');
fprintf(fid, 'Fast option: remove globally bad channels from all subjects.\n');
fprintf(fid, 'Exploratory option: subject-specific masking is possible but should be described as exploratory.\n');
fclose(fid);

%% ================= PLOTS =================
make_plots(Tglobal, Tsubject, cfg);

fprintf('\n[DONE SUBJECT-WISE QC]\n');
fprintf('Session CSV: %s\n', sessionCsv);
fprintf('Subject summary CSV: %s\n', subjectCsv);
fprintf('Subject lists CSV: %s\n', listsCsv);
fprintf('Global summary CSV: %s\n', globalCsv);
fprintf('MAT: %s\n', matPath);
fprintf('TXT: %s\n', txtPath);

fprintf('\nGlobal bad channels: %s\n', mat2str(globalBadChannels));
fprintf('Global good channels: %s\n', mat2str(globalGoodChannels));
disp(Tglobal);

end

%% ========================================================================
function subj = infer_subject(folder, fname)
[~, parent] = fileparts(folder);
if startsWith(lower(parent), 's') && numel(parent) <= 5
    subj = parent;
else
    tok = regexp(fname, '(s\d+)', 'tokens', 'once', 'ignorecase');
    if ~isempty(tok)
        subj = tok{1};
    else
        subj = parent;
    end
end
end

function A = get_phase_array(SEG, ph)
A = [];
names = fieldnames(SEG);
for i = 1:numel(names)
    if strcmpi(names{i}, ph)
        A = SEG.(names{i});
        return;
    end
end
end

function M = compute_basic_qc(A)
A = double(A);
nCh = size(A,1);

Std = NaN(nCh,1);
P2P = NaN(nCh,1);
DerivSTD = NaN(nCh,1);
CorrRef = NaN(nCh,1);

Xall = reshape(A, nCh, []);
ref = median(Xall, 1, 'omitnan');

for ch = 1:nCh
    x = Xall(ch,:);
    xgood = x(isfinite(x));

    if isempty(xgood), continue; end

    Std(ch) = std(xgood);

    try
        P2P(ch) = prctile(xgood,99.5) - prctile(xgood,0.5);
    catch
        P2P(ch) = max(xgood) - min(xgood);
    end

    xd = squeeze(diff(A(ch,:,:),1,2));
    xd = xd(:);
    xd = xd(isfinite(xd));
    if ~isempty(xd)
        DerivSTD(ch) = std(xd);
    end

    good = isfinite(x) & isfinite(ref);
    if sum(good) > 20 && std(x(good)) > eps && std(ref(good)) > eps
        C = corrcoef(double(x(good)), double(ref(good)));
        CorrRef(ch) = C(1,2);
    end
end

M.Std = Std;
M.P2P = P2P;
M.DerivSTD = DerivSTD;
M.CorrRef = CorrRef;
end

function F = flag_bad(M, cfg)
Z_STD = robust_z(M.Std);
Z_P2P = robust_z(M.P2P);
Z_Deriv = robust_z(M.DerivSTD);
Z_CorrRef = robust_z(M.CorrRef);

medStd = median(M.Std(isfinite(M.Std)));
if isempty(medStd) || ~isfinite(medStd), medStd = 0; end

FlagFlat = M.Std < max(1e-12, cfg.flatStdRatioThr * medStd);
FlagHighSTD = abs(Z_STD) > cfg.robustZthr;
FlagHighP2P = Z_P2P > cfg.robustZthr;
FlagHighDeriv = Z_Deriv > cfg.robustZthr;
FlagLowCorr = Z_CorrRef < -cfg.robustZthr;

BadScore = double(FlagFlat) + double(FlagHighSTD) + double(FlagHighP2P) + ...
    double(FlagHighDeriv) + double(FlagLowCorr);

IsBad = BadScore >= cfg.badScoreThr;
IsSuspicious = BadScore == 1;

F.Z_STD = Z_STD;
F.Z_P2P = Z_P2P;
F.Z_Deriv = Z_Deriv;
F.Z_CorrRef = Z_CorrRef;

F.FlagFlat = FlagFlat;
F.FlagHighSTD = FlagHighSTD;
F.FlagHighP2P = FlagHighP2P;
F.FlagHighDeriv = FlagHighDeriv;
F.FlagLowCorr = FlagLowCorr;

F.BadScore = BadScore;
F.IsBad = IsBad;
F.IsSuspicious = IsSuspicious;
end

function z = robust_z(x)
x = double(x(:));
good = isfinite(x);
z = zeros(size(x));
if ~any(good), return; end
med = median(x(good));
mad0 = median(abs(x(good)-med));
if ~isfinite(mad0) || mad0 < eps
    mad0 = std(x(good));
end
if ~isfinite(mad0) || mad0 < eps
    mad0 = 1;
end
z = 0.6745 * (x-med) ./ mad0;
z(~isfinite(z)) = 0;
end

function make_plots(Tglobal, Tsubject, cfg)
T = sortrows(Tglobal, 'Channel', 'ascend');

fig = figure('Color','w','Position',[100 100 1200 500]);
bar(T.Channel, T.SubjectBadFraction);
hold on;
yline(cfg.globalBadSubjectFractionThr, '--r', 'Global bad threshold');
xlabel('Channel');
ylabel('Fraction of subjects where channel is bad');
title('Subject-wise bad-channel frequency');
grid on;
saveas(fig, fullfile(cfg.outDir, 'PLOT_global_subject_bad_fraction.png'));
close(fig);

subjects = unique(Tsubject.Subject, 'stable');
channels = 1:64;
H = NaN(numel(subjects), numel(channels));

for i=1:numel(subjects)
    for ch=channels
        idx = strcmp(Tsubject.Subject, subjects{i}) & Tsubject.Channel == ch;
        if any(idx)
            H(i,ch) = Tsubject.BadFraction(find(idx,1));
        end
    end
end

fig = figure('Color','w','Position',[100 100 1400 700]);
imagesc(channels, 1:numel(subjects), H);
colorbar;
xlabel('Channel');
ylabel('Subject');
set(gca, 'YTick', 1:numel(subjects), 'YTickLabel', subjects);
title('Bad fraction per subject and channel');
saveas(fig, fullfile(cfg.outDir, 'PLOT_subject_by_channel_bad_fraction_heatmap.png'));
close(fig);
end

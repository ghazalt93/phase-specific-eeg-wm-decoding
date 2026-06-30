function strict_channel_qc_psd()

cfg = get_project_config();

clc; close all;

%% ================= CONFIG =================
cfg.projectDir = fullfile(cfg.dataRoot, 'best');
cfg.outBaseDir = cfg.outputRoot;
cfg.outDir = fullfile(cfg.outBaseDir, 'STRICT_CHANNEL_QC_PSD');

if ~exist(cfg.outBaseDir,'dir'), mkdir(cfg.outBaseDir); end
if ~exist(cfg.outDir,'dir'), mkdir(cfg.outDir); end

candidateDirs = { ...
    fullfile(cfg.projectDir, 'DERIVED_RESULTS'), ...
    fullfile(cfg.projectDir, 'Subjects', 'DERIVED_RESULTS')};

cfg.derivedDir = '';
for i = 1:numel(candidateDirs)
    if exist(candidateDirs{i}, 'dir')
        cfg.derivedDir = candidateDirs{i};
        break;
    end
end
if isempty(cfg.derivedDir)
    error('DERIVED_RESULTS folder not found. Check cfg.projectDir.');
end

cfg.fs = 250;
cfg.phases = {'Stim','Maint','Retr'};

cfg.maxTrialsPerPhase = 40;
cfg.maxFiles = Inf;

cfg.psd.totalBand = [1 45];
cfg.psd.deltaBand = [1 4];
cfg.psd.thetaBand = [4 8];
cfg.psd.alphaBand = [8 13];
cfg.psd.betaBand  = [13 30];
cfg.psd.hfBand    = [30 45];  
cfg.psd.slopeBand = [2 40];

cfg.zHighAmpThr     = 4.0;
cfg.zHighP2PThr     = 4.0;
cfg.zHighDerivThr   = 4.0;
cfg.zHighHFRatioThr = 3.0;
cfg.zLowHFRatioThr  = -4.0;
cfg.zTotalPowThr    = 4.0;
cfg.zSlopeThr       = 4.0;

cfg.flatStdRatioThr = 0.10;  

cfg.minCorrAbs = 0.05;      
cfg.zLowCorrThr = -3.5;

cfg.veryStrictExcellentFraction = 0.90;
cfg.strictExcellentFraction     = 0.80;
cfg.moderateExcellentFraction   = 0.70;

cfg.maxBadFractionForStrict = 0.05;
cfg.maxBadFractionForModerate = 0.10;

fprintf('\n[STRICT CHANNEL QC + PSD]\n');
fprintf('derivedDir: %s\n', cfg.derivedDir);
fprintf('outDir:     %s\n', cfg.outDir);

files = dir(fullfile(cfg.derivedDir, '**', '*__DERIVED.mat'));
if isempty(files)
    error('No *__DERIVED.mat found under %s', cfg.derivedDir);
end
if isfinite(cfg.maxFiles)
    files = files(1:min(cfg.maxFiles, numel(files)));
end
fprintf('Files to scan: %d\n', numel(files));

%% ================= SCAN =================
rows = {};
ri = 0;

for f = 1:numel(files)
    filePath = fullfile(files(f).folder, files(f).name);
    [~, runTag] = fileparts(files(f).name);
    subj = infer_subject(files(f).folder, files(f).name);

    try
        S = load(filePath, 'SEG');
        if ~isfield(S,'SEG')
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
        nTrials = size(A,3);

        if isfinite(cfg.maxTrialsPerPhase) && nTrials > cfg.maxTrialsPerPhase
            useTr = round(linspace(1, nTrials, cfg.maxTrialsPerPhase));
            A = A(:,:,useTr);
        end

        M = compute_metrics_psd(A, cfg);
        F = flag_channels_strict(M, cfg);

        badCh = find(F.IsBad);
        excellentCh = find(F.IsExcellent);

        fprintf('  %s | ch=%d time=%d trials=%d | bad=%s | excellent n=%d\n', ...
            ph, nCh, nTime, nTrials, mat2str(badCh'), numel(excellentCh));

        for ch = 1:nCh
            ri = ri + 1;
            rows(ri,:) = { ...
                subj, runTag, ph, ch, nTime, nTrials, ...
                M.Std(ch), M.MAD(ch), M.P2P(ch), M.DerivSTD(ch), M.AbsMean(ch), ...
                M.TotalPower(ch), M.LogTotalPower(ch), M.DeltaRatio(ch), M.ThetaRatio(ch), ...
                M.AlphaRatio(ch), M.BetaRatio(ch), M.HFRatio(ch), M.PSDSlope(ch), M.CorrRef(ch), ...
                F.Z_STD(ch), F.Z_P2P(ch), F.Z_Deriv(ch), F.Z_HFRatio(ch), F.Z_TotalPower(ch), F.Z_Slope(ch), F.Z_CorrRef(ch), ...
                F.FlagFlat(ch), F.FlagHighAmp(ch), F.FlagHighP2P(ch), F.FlagHighDeriv(ch), ...
                F.FlagHighHF(ch), F.FlagAbnormalTotalPower(ch), F.FlagAbnormalSlope(ch), F.FlagLowCorr(ch), ...
                F.BadScore(ch), F.SuspiciousScore(ch), F.IsBad(ch), F.IsExcellent(ch), ...
                filePath};
        end
    end
end

if isempty(rows)
    error('No QC rows created.');
end

Tsession = cell2table(rows, 'VariableNames', { ...
    'Subject','RunTag','Phase','Channel','Ntime','NtrialsOriginal', ...
    'STD','MAD','P2P','DerivSTD','AbsMean', ...
    'TotalPower','LogTotalPower','DeltaRatio','ThetaRatio','AlphaRatio','BetaRatio','HFRatio','PSDSlope','CorrRef', ...
    'Z_STD','Z_P2P','Z_DerivSTD','Z_HFRatio','Z_TotalPower','Z_PSDSlope','Z_CorrRef', ...
    'FlagFlat','FlagHighAmp','FlagHighP2P','FlagHighDerivSTD', ...
    'FlagHighHF','FlagAbnormalTotalPower','FlagAbnormalSlope','FlagLowCorr', ...
    'BadScore','SuspiciousScore','IsBad','IsExcellent','SourceFile'});

sessionCsv = fullfile(cfg.outDir, 'STRICT_PSD_QC_session.csv');
writetable(Tsession, sessionCsv);

%% ================= GLOBAL SUMMARY =================
[G, ch] = findgroups(Tsession.Channel);

TotalCount = splitapply(@numel, Tsession.IsBad, G);
BadCount = splitapply(@sum, double(Tsession.IsBad), G);
ExcellentCount = splitapply(@sum, double(Tsession.IsExcellent), G);

BadFraction = BadCount ./ max(TotalCount,1);
ExcellentFraction = ExcellentCount ./ max(TotalCount,1);

STD_median = splitapply(@(x) mean_omitnan(x), Tsession.STD, G);
HFRatio_median = splitapply(@(x) mean_omitnan(x), Tsession.HFRatio, G);
CorrRef_median = splitapply(@(x) mean_omitnan(x), Tsession.CorrRef, G);
TotalPower_median = splitapply(@(x) mean_omitnan(x), Tsession.TotalPower, G);
Slope_median = splitapply(@(x) mean_omitnan(x), Tsession.PSDSlope, G);

Tglobal = table(ch, TotalCount, BadCount, BadFraction, ExcellentCount, ExcellentFraction, ...
    STD_median, TotalPower_median, HFRatio_median, Slope_median, CorrRef_median, ...
    'VariableNames', {'Channel','TotalCount','BadCount','BadFraction','ExcellentCount','ExcellentFraction', ...
    'STD_mean','TotalPower_mean','HFRatio_mean','PSDSlope_mean','CorrRef_mean'});

Tglobal.SelectVeryStrict = Tglobal.ExcellentFraction >= cfg.veryStrictExcellentFraction & ...
    Tglobal.BadFraction <= cfg.maxBadFractionForStrict;

Tglobal.SelectStrict = Tglobal.ExcellentFraction >= cfg.strictExcellentFraction & ...
    Tglobal.BadFraction <= cfg.maxBadFractionForStrict;

Tglobal.SelectModerate = Tglobal.ExcellentFraction >= cfg.moderateExcellentFraction & ...
    Tglobal.BadFraction <= cfg.maxBadFractionForModerate;

Tglobal = sortrows(Tglobal, {'ExcellentFraction','BadFraction'}, {'descend','ascend'});

excellentChannels_veryStrict = sort(Tglobal.Channel(Tglobal.SelectVeryStrict))';
excellentChannels_strict     = sort(Tglobal.Channel(Tglobal.SelectStrict))';
excellentChannels_moderate   = sort(Tglobal.Channel(Tglobal.SelectModerate))';

nTopList = [16 24 32 48];
topChannels = struct();
for i = 1:numel(nTopList)
    nTop = min(nTopList(i), height(Tglobal));
    fn = sprintf('top%d_byQuality', nTopList(i));
    topChannels.(fn) = sort(Tglobal.Channel(1:nTop))';
end

globalCsv = fullfile(cfg.outDir, 'STRICT_PSD_QC_global.csv');
writetable(Tglobal, globalCsv);

matPath = fullfile(cfg.outDir, 'STRICT_PSD_QC_channels.mat');
save(matPath, 'cfg', 'Tsession', 'Tglobal', ...
    'excellentChannels_veryStrict', 'excellentChannels_strict', 'excellentChannels_moderate', 'topChannels');

txtPath = fullfile(cfg.outDir, 'STRICT_PSD_QC_selected_channels.txt');
fid = fopen(txtPath, 'w');
fprintf(fid, 'Very strict threshold: ExcellentFraction >= %.2f and BadFraction <= %.2f\n', ...
    cfg.veryStrictExcellentFraction, cfg.maxBadFractionForStrict);
fprintf(fid, 'excellentChannels_veryStrict = %s\n\n', mat2str(excellentChannels_veryStrict));

fprintf(fid, 'Strict threshold: ExcellentFraction >= %.2f and BadFraction <= %.2f\n', ...
    cfg.strictExcellentFraction, cfg.maxBadFractionForStrict);
fprintf(fid, 'excellentChannels_strict = %s\n\n', mat2str(excellentChannels_strict));

fprintf(fid, 'Moderate threshold: ExcellentFraction >= %.2f and BadFraction <= %.2f\n', ...
    cfg.moderateExcellentFraction, cfg.maxBadFractionForModerate);
fprintf(fid, 'excellentChannels_moderate = %s\n\n', mat2str(excellentChannels_moderate));

fns = fieldnames(topChannels);
for i=1:numel(fns)
    fprintf(fid, '%s = %s\n', fns{i}, mat2str(topChannels.(fns{i})));
end
fclose(fid);

%% ================= PLOTS =================
make_plots(Tglobal, cfg);

fprintf('\n[DONE STRICT PSD QC]\n');
fprintf('Session CSV: %s\n', sessionCsv);
fprintf('Global CSV:  %s\n', globalCsv);
fprintf('MAT:         %s\n', matPath);
fprintf('TXT:         %s\n', txtPath);

fprintf('\nSelected channels:\n');
fprintf('  veryStrict n=%d: %s\n', numel(excellentChannels_veryStrict), mat2str(excellentChannels_veryStrict));
fprintf('  strict     n=%d: %s\n', numel(excellentChannels_strict), mat2str(excellentChannels_strict));
fprintf('  moderate   n=%d: %s\n', numel(excellentChannels_moderate), mat2str(excellentChannels_moderate));
disp(Tglobal(:, {'Channel','BadFraction','ExcellentFraction','SelectVeryStrict','SelectStrict','SelectModerate'}));

end

%% ========================================================================
function subj = infer_subject(folder, fname)
[~, parent] = fileparts(folder);
if startsWith(lower(parent), 's') && numel(parent) <= 5
    subj = parent;
else
    tok = regexp(fname, '(s\d+)', 'tokens', 'once', 'ignorecase');
    if ~isempty(tok), subj = tok{1}; else, subj = parent; end
end
end

function A = get_phase_array(SEG, ph)
A = [];
names = fieldnames(SEG);
for i=1:numel(names)
    if strcmpi(names{i}, ph)
        A = SEG.(names{i});
        return;
    end
end
end

function M = compute_metrics_psd(A, cfg)
A = double(A);
nCh = size(A,1);
nT = size(A,2);
nTr = size(A,3);

Std = NaN(nCh,1);
MAD = NaN(nCh,1);
P2P = NaN(nCh,1);
DerivSTD = NaN(nCh,1);
AbsMean = NaN(nCh,1);

TotalPower = NaN(nCh,1);
LogTotalPower = NaN(nCh,1);
DeltaRatio = NaN(nCh,1);
ThetaRatio = NaN(nCh,1);
AlphaRatio = NaN(nCh,1);
BetaRatio = NaN(nCh,1);
HFRatio = NaN(nCh,1);
PSDSlope = NaN(nCh,1);
CorrRef = NaN(nCh,1);

XcatAll = reshape(A, nCh, []);
refAll = median(XcatAll, 1, 'omitnan');

for ch=1:nCh
    X = squeeze(A(ch,:,:)); % time x trials
    xcat = X(:);
    xcat = xcat(isfinite(xcat));

    if isempty(xcat), continue; end

    Std(ch) = std(xcat);
    MAD(ch) = median(abs(xcat - median(xcat)));
    AbsMean(ch) = abs(mean(xcat));

    try
        P2P(ch) = prctile(xcat,99.5) - prctile(xcat,0.5);
    catch
        P2P(ch) = max(xcat) - min(xcat);
    end

    d = diff(X,1,1);
    d = d(:);
    d = d(isfinite(d));
    if ~isempty(d)
        DerivSTD(ch) = std(d);
    end

    xfull = XcatAll(ch,:);
    good = isfinite(xfull) & isfinite(refAll);
    if sum(good) > 20 && std(xfull(good)) > eps && std(refAll(good)) > eps
        C = corrcoef(double(xfull(good)), double(refAll(good)));
        CorrRef(ch) = C(1,2);
    end

    % PSD per trial and average
    Pacc = [];
    Freq = [];
    for tr=1:nTr
        x = X(:,tr);
        x = x(isfinite(x));
        if numel(x) < 16, continue; end
        x = x - mean(x);

        wlen = min(numel(x), 256);
        if wlen < 16, continue; end
        nover = floor(wlen/2);
        nfft = max(256, 2^nextpow2(wlen));

        try
            [pxx, f] = pwelch(x, hamming(wlen), nover, nfft, cfg.fs);
        catch
            [pxx, f] = periodogram(x, [], nfft, cfg.fs);
        end

        if isempty(Pacc)
            Pacc = zeros(size(pxx));
            Freq = f;
        end
        Pacc = Pacc + pxx;
    end

    if isempty(Pacc), continue; end
    P = Pacc ./ nTr;
    f = Freq;

    total = band_int(f, P, cfg.psd.totalBand);
    if ~isfinite(total) || total <= 0, continue; end

    TotalPower(ch) = total;
    LogTotalPower(ch) = log10(total + eps);
    DeltaRatio(ch) = band_int(f, P, cfg.psd.deltaBand) / total;
    ThetaRatio(ch) = band_int(f, P, cfg.psd.thetaBand) / total;
    AlphaRatio(ch) = band_int(f, P, cfg.psd.alphaBand) / total;
    BetaRatio(ch)  = band_int(f, P, cfg.psd.betaBand)  / total;
    HFRatio(ch)    = band_int(f, P, cfg.psd.hfBand)    / total;

    idxSlope = f >= cfg.psd.slopeBand(1) & f <= cfg.psd.slopeBand(2) & P > 0;
    if sum(idxSlope) >= 5
        xx = log10(f(idxSlope));
        yy = log10(P(idxSlope));
        pp = polyfit(xx(:), yy(:), 1);
        PSDSlope(ch) = pp(1);
    end
end

M.Std = Std; M.MAD = MAD; M.P2P = P2P; M.DerivSTD = DerivSTD; M.AbsMean = AbsMean;
M.TotalPower = TotalPower; M.LogTotalPower = LogTotalPower;
M.DeltaRatio = DeltaRatio; M.ThetaRatio = ThetaRatio; M.AlphaRatio = AlphaRatio;
M.BetaRatio = BetaRatio; M.HFRatio = HFRatio; M.PSDSlope = PSDSlope; M.CorrRef = CorrRef;
end

function val = band_int(f, P, band)
idx = f >= band(1) & f <= band(2) & isfinite(P);
if sum(idx) < 2
    val = NaN;
else
    val = trapz(f(idx), P(idx));
end
end

function F = flag_channels_strict(M, cfg)
Z_STD = robust_z(M.Std);
Z_P2P = robust_z(M.P2P);
Z_Deriv = robust_z(M.DerivSTD);
Z_HFRatio = robust_z(M.HFRatio);
Z_TotalPower = robust_z(M.LogTotalPower);
Z_Slope = robust_z(M.PSDSlope);
Z_CorrRef = robust_z(M.CorrRef);

medStd = median(M.Std(isfinite(M.Std)));
if isempty(medStd) || ~isfinite(medStd), medStd = 0; end

FlagFlat = M.Std < max(1e-12, cfg.flatStdRatioThr * medStd);
FlagHighAmp = abs(Z_STD) > cfg.zHighAmpThr;
FlagHighP2P = Z_P2P > cfg.zHighP2PThr;
FlagHighDeriv = Z_Deriv > cfg.zHighDerivThr;
FlagHighHF = Z_HFRatio > cfg.zHighHFRatioThr;
FlagAbnormalTotalPower = abs(Z_TotalPower) > cfg.zTotalPowThr;
FlagAbnormalSlope = abs(Z_Slope) > cfg.zSlopeThr;
FlagLowCorr = (M.CorrRef < cfg.minCorrAbs) | (Z_CorrRef < cfg.zLowCorrThr);

BadScore = double(FlagFlat) + double(FlagHighAmp) + double(FlagHighP2P) + ...
    double(FlagHighDeriv) + double(FlagHighHF) + double(FlagAbnormalTotalPower) + ...
    double(FlagAbnormalSlope) + double(FlagLowCorr);

SuspiciousScore = BadScore;

IsBad = BadScore >= 2 | (FlagHighHF & FlagHighDeriv);

IsExcellent = BadScore == 0;

F.Z_STD = Z_STD; F.Z_P2P = Z_P2P; F.Z_Deriv = Z_Deriv; F.Z_HFRatio = Z_HFRatio;
F.Z_TotalPower = Z_TotalPower; F.Z_Slope = Z_Slope; F.Z_CorrRef = Z_CorrRef;
F.FlagFlat = FlagFlat; F.FlagHighAmp = FlagHighAmp; F.FlagHighP2P = FlagHighP2P;
F.FlagHighDeriv = FlagHighDeriv; F.FlagHighHF = FlagHighHF;
F.FlagAbnormalTotalPower = FlagAbnormalTotalPower; F.FlagAbnormalSlope = FlagAbnormalSlope;
F.FlagLowCorr = FlagLowCorr;
F.BadScore = BadScore; F.SuspiciousScore = SuspiciousScore;
F.IsBad = IsBad; F.IsExcellent = IsExcellent;
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
z = 0.6745 * (x - med) ./ mad0;
z(~isfinite(z)) = 0;
end

function m = mean_omitnan(x)
x = x(isfinite(x));
if isempty(x), m = NaN; else, m = mean(x); end
end

function make_plots(Tglobal, cfg)
% sort back by channel
T = sortrows(Tglobal, 'Channel', 'ascend');

fig = figure('Color','w','Position',[100 100 1200 500]);
bar(T.Channel, T.ExcellentFraction);
hold on;
yline(cfg.strictExcellentFraction, '--r', 'strict');
yline(cfg.veryStrictExcellentFraction, '--k', 'very strict');
xlabel('Channel'); ylabel('Excellent fraction');
title('Strict PSD-QC: Excellent fraction per channel');
grid on;
saveas(fig, fullfile(cfg.outDir, 'PLOT_excellent_fraction_per_channel.png'));
close(fig);

fig = figure('Color','w','Position',[100 100 1200 500]);
bar(T.Channel, T.BadFraction);
hold on;
yline(cfg.maxBadFractionForStrict, '--r', 'max bad strict');
xlabel('Channel'); ylabel('Bad fraction');
title('Strict PSD-QC: Bad fraction per channel');
grid on;
saveas(fig, fullfile(cfg.outDir, 'PLOT_bad_fraction_per_channel.png'));
close(fig);

fig = figure('Color','w','Position',[100 100 1200 500]);
bar(T.Channel, T.HFRatio_mean);
xlabel('Channel'); ylabel('Mean HF ratio 30-45 / 1-45');
title('High-frequency ratio per channel');
grid on;
saveas(fig, fullfile(cfg.outDir, 'PLOT_HF_ratio_per_channel.png'));
close(fig);

fig = figure('Color','w','Position',[100 100 1200 500]);
bar(T.Channel, T.CorrRef_mean);
xlabel('Channel'); ylabel('Mean corr with median reference');
title('Correlation with robust median reference');
grid on;
saveas(fig, fullfile(cfg.outDir, 'PLOT_corr_reference_per_channel.png'));
close(fig);
end

%% STEP28B_FEATURE_FAMILY_RECLASSIFY_AND_QC

cfg = get_project_config();

clear; clc;

ROOT = cfg.outputRoot;
DSpath = fullfile(ROOT, 'Subjects', '_wm_ml', 'dataset.mat');
OUTDIR = fullfile(ROOT, '_wm_reports', 'STEP28B_feature_family_reclassified');
if exist(OUTDIR,'dir') ~= 7, mkdir(OUTDIR); end

robustZ_thr = 8;
high_naninf_thr = 0.05;
high_extreme_thr = 0.01;

fprintf('\n=== STEP28B FEATURE FAMILY RECLASSIFY AND QC ===\n');
fprintf('Dataset: %s\nOutput: %s\n', DSpath, OUTDIR);

S = load(DSpath);
D = get_table(S);
vars = D.Properties.VariableNames;
low = lower(vars);

subjCol = pick_col(vars, low, {'subject','subj','subjid','participant'});
runCol = pick_col(vars, low, {'run','runnum','session'});
phaseCol = pick_col(vars, low, {'phase','phasename'});
labelCol = pick_col(vars, low, {'ycondition','condition','cond','y','label','class'});
trialCol = pick_col(vars, low, {'trialnum','trial','trialid'});

isNum = false(1,numel(vars));
for i = 1:numel(vars)
    x = D.(vars{i});
    isNum(i) = isnumeric(x) && isvector(x) && numel(x)==height(D);
end

isMeta = false(1,numel(vars));
metaNames = {subjCol, runCol, phaseCol, labelCol, trialCol, ...
    'correct','ycorrect','rt','reactiontime','response','patternid','session'};
for i = 1:numel(vars)
    for j = 1:numel(metaNames)
        if ~isempty(metaNames{j}) && strcmpi(vars{i}, metaNames{j})
            isMeta(i) = true;
        end
    end
    if contains(low{i}, 'rt') || contains(low{i}, 'reaction') || contains(low{i}, 'response') || ...
       contains(low{i}, 'correct') || contains(low{i}, 'accuracy')
        isMeta(i) = true;
    end
end

featIdxAll = find(isNum & ~isMeta);
featVarsAll = string(vars(featIdxAll));
Channel = nan(numel(featVarsAll),1);
Family = strings(numel(featVarsAll),1);

for i = 1:numel(featVarsAll)
    Channel(i) = get_ch(featVarsAll(i));
    Family(i) = family2(featVarsAll(i));
end

hasCh = ~isnan(Channel);
featVars = featVarsAll(hasCh);
featIdx = featIdxAll(hasCh);
Channel = Channel(hasCh);
Family = Family(hasCh);

fprintf('Rows=%d | Vars=%d | channel-coded numeric features=%d\n', height(D), width(D), numel(featVars));

X = double(table2array(D(:, cellstr(featVars))));

% Inventory
[GF, Fam] = findgroups(Family);
NFeatures = splitapply(@numel, Family, GF);
NChannels = splitapply(@(x) numel(unique(x)), Channel, GF);
Examples = strings(numel(Fam),1);
for i = 1:numel(Fam)
    idx = find(Family == Fam(i));
    Examples(i) = strjoin(cellstr(featVars(idx(1:min(6,numel(idx))))), ', ');
end
InvFamily = table(Fam, NFeatures, NChannels, Examples, ...
    'VariableNames', {'Family','NFeatures','NChannels','Examples'});
InvFamily = sortrows(InvFamily, 'NFeatures','descend');
writetable(InvFamily, fullfile(OUTDIR, 'STEP28B_feature_inventory_by_family.csv'));

[GFC, Fam2, Ch2] = findgroups(Family, Channel);
NFC = splitapply(@numel, Family, GFC);
InvFC = table(Fam2, Ch2, NFC, 'VariableNames', {'Family','Channel','NFeatures'});
writetable(InvFC, fullfile(OUTDIR, 'STEP28B_feature_inventory_family_channel.csv'));

% Quality audit
naninfRate = mean(~isfinite(X),1)';
zeroVar = false(numel(featVars),1);
for j = 1:numel(featVars)
    x = X(:,j); x = x(isfinite(x));
    zeroVar(j) = isempty(x) || std(x) < 1e-12;
end

med = nanmedian_cols(X);
mad0 = nanmad_cols(X, med);
scale = 1.4826 * mad0;
scale(~isfinite(scale) | scale < 1e-12) = NaN;
extremeRate = nan(numel(featVars),1);
for j = 1:numel(featVars)
    if isfinite(scale(j))
        z = abs((X(:,j)-med(j))./scale(j));
        extremeRate(j) = mean(isfinite(z) & z > robustZ_thr);
    end
end

QualityT = table(featVars(:), Family(:), Channel(:), naninfRate, extremeRate, zeroVar, ...
    naninfRate>=high_naninf_thr, extremeRate>=high_extreme_thr, ...
    'VariableNames', {'Feature','Family','Channel','NaNInfRate','ExtremeRate','ZeroVariance','HighNaNInfFlag','HighExtremeFlag'});
writetable(QualityT, fullfile(OUTDIR, 'STEP28B_feature_quality_audit.csv'));

QualitySummary = table( ...
    ["N_channel_coded_features";"N_high_NaNInf_features";"N_high_extreme_features";"N_zero_variance_features";"Median_NaNInf_rate";"Median_extreme_rate"], ...
    [numel(featVars); sum(QualityT.HighNaNInfFlag); sum(QualityT.HighExtremeFlag); sum(QualityT.ZeroVariance); median(naninfRate,'omitnan'); median(extremeRate,'omitnan')], ...
    'VariableNames', {'Metric','Value'});
writetable(QualitySummary, fullfile(OUTDIR, 'STEP28B_feature_quality_summary.csv'));

writetable(sortrows(QualityT(QualityT.HighNaNInfFlag,:), 'NaNInfRate','descend'), fullfile(OUTDIR, 'STEP28B_high_naninf_features.csv'));
writetable(sortrows(QualityT(QualityT.HighExtremeFlag,:), 'ExtremeRate','descend'), fullfile(OUTDIR, 'STEP28B_high_extreme_features.csv'));
writetable(QualityT(QualityT.ZeroVariance,:), fullfile(OUTDIR, 'STEP28B_constant_features.csv'));

% Figures
try
    fig = figure('Color','w','Position',[100 100 850 520]);
    bar(InvFamily.NFeatures);
    set(gca,'XTick',1:height(InvFamily),'XTickLabel',InvFamily.Family,'XTickLabelRotation',25);
    ylabel('Number of features');
    title('Feature inventory by family, reclassified');
    grid on;
    saveas(fig, fullfile(OUTDIR, 'FIG_STEP28B_feature_family_counts.png'));
    savefig(fig, fullfile(OUTDIR, 'FIG_STEP28B_feature_family_counts.fig'));
    close(fig);
catch ME
    warning('Family count plot failed: %s', ME.message);
end

try
    fig = figure('Color','w','Position',[100 100 700 480]);
    cats = categorical({'High NaN/Inf','High extreme','Zero variance'});
    vals = [sum(QualityT.HighNaNInfFlag), sum(QualityT.HighExtremeFlag), sum(QualityT.ZeroVariance)];
    bar(cats, vals);
    ylabel('Number of features');
    title('Feature-quality audit summary, reclassified');
    grid on;
    saveas(fig, fullfile(OUTDIR, 'FIG_STEP28B_feature_quality_summary.png'));
    savefig(fig, fullfile(OUTDIR, 'FIG_STEP28B_feature_quality_summary.fig'));
    close(fig);
catch ME
    warning('Quality plot failed: %s', ME.message);
end

txt = fullfile(OUTDIR, 'STEP28B_INTERPRETATION_SUMMARY.txt');
fid = fopen(txt,'w');
fprintf(fid, 'STEP28B FEATURE FAMILY RECLASSIFICATION SUMMARY\n');
fprintf(fid, '==============================================\n\n');
fprintf(fid, 'Dataset rows: %d\n', height(D));
fprintf(fid, 'Dataset variables: %d\n', width(D));
fprintf(fid, 'Channel-coded numeric features: %d\n\n', numel(featVars));
for i = 1:height(InvFamily)
    fprintf(fid, '%s: %d features across %d channels. Examples: %s\n', InvFamily.Family(i), InvFamily.NFeatures(i), InvFamily.NChannels(i), InvFamily.Examples(i));
end
fprintf(fid, '\nQuality summary:\n');
for i = 1:height(QualitySummary)
    fprintf(fid, '%s = %.6g\n', QualitySummary.Metric(i), QualitySummary.Value(i));
end
fclose(fid);

fprintf('\nDONE STEP28B.\n');
disp(InvFamily)
disp(QualitySummary)

%% LOCAL FUNCTIONS
function D = get_table(S)
    if isfield(S,'DS'), D = S.DS; return; end
    if isfield(S,'dataset'), D = S.dataset; return; end
    if isfield(S,'D'), D = S.D; return; end
    if isfield(S,'T'), D = S.T; return; end
    fn = fieldnames(S);
    for i = 1:numel(fn)
        if istable(S.(fn{i})), D = S.(fn{i}); return; end
    end
    error('No table found in MAT.');
end

function out = pick_col(vars, low, candidates)
    out = '';
    for i = 1:numel(candidates)
        k = find(strcmp(low, lower(candidates{i})), 1);
        if ~isempty(k), out = vars{k}; return; end
    end
    for i = 1:numel(candidates)
        k = find(contains(low, lower(candidates{i})), 1);
        if ~isempty(k), out = vars{k}; return; end
    end
end

function ch = get_ch(name)
    ch = NaN; s = char(name);
    tok = regexp(s, 'ch[_-]?(\d{1,2})', 'tokens', 'once', 'ignorecase');
    if isempty(tok), tok = regexp(s, 'chan(?:nel)?[_-]?(\d{1,2})', 'tokens', 'once', 'ignorecase'); end
    if ~isempty(tok)
        ch = str2double(tok{1});
        if ch<1 || ch>64, ch = NaN; end
    end
end

function fam = family2(name)
    s = lower(string(name));
    % waveform descriptors first, because line length / AUC are waveform not generic temporal
    if contains(s,"aucabs") || contains(s,"linelen") || contains(s,"line_len") || contains(s,"line_length") || contains(s,"tmaxabs")
        fam = "waveformDescriptor";
        elseif contains(s,"erpptp") || contains(s,"erp_ptp") || contains(s,"erp") || contains(s,"ptp_0_100") || contains(s,"ptp_100_200") || contains(s,"ptp_200_300") || contains(s,"ptp_300_500")
    fam = "ERPwindow";
    elseif startsWith(s,"bp_") || startsWith(s,"rbp_") || contains(s,"bandpower") || contains(s,"pow_") || startsWith(s,"pow") || ...
           contains(s,"psd") || contains(s,"delta") || contains(s,"theta") || contains(s,"alpha") || contains(s,"beta") || contains(s,"gamma")
        fam = "bandpower";
    elseif contains(s,"lzc") || contains(s,"dfa") || contains(s,"hurst") || contains(s,"hurst_rs") || ...
       contains(s,"entropy") || contains(s,"sampen") || contains(s,"apen") || ...
       contains(s,"fractal") || contains(s,"higuchi") || contains(s,"katz") || contains(s,"complex")
    fam = "complexity";
    elseif contains(s,"rms") || contains(s,"peak") || contains(s,"mean") || contains(s,"std") || contains(s,"p2p") || ...
           contains(s,"deriv") || contains(s,"tkeo") || contains(s,"var") || contains(s,"skew") || contains(s,"kurt") || ...
           contains(s,"median") || contains(s,"iqr") || contains(s,"mad")
        fam = "temporalStat";
    else
        fam = "other";
    end
end

function m = nanmedian_cols(X)
    m = NaN(1,size(X,2));
    for j = 1:size(X,2)
        x = X(:,j); x = x(isfinite(x));
        if ~isempty(x), m(j) = median(x); end
    end
end

function m = nanmad_cols(X, med)
    m = NaN(1,size(X,2));
    for j = 1:size(X,2)
        x = X(:,j); x = x(isfinite(x));
        if ~isempty(x) && isfinite(med(j)), m(j) = median(abs(x-med(j))); end
    end
end

%% STEP85_collect_all_significant_STEP84_STEP84B.m

cfg = get_project_config();

clear; clc;

%% ===================== CONFIG =====================

cfg = struct();

cfg.ROOT = cfg.outputRoot;

cfg.step84Dir  = fullfile(cfg.ROOT, '_wm_STEP84_channel_featuregroup_subjectlevel_stats');
cfg.step84BDir = fullfile(cfg.ROOT, '_wm_STEP84B_channel_featuregroup_KW_ranksum_subjectAggregated');

cfg.outDir = fullfile(cfg.ROOT, '_wm_STEP85_collected_significant_effects');

cfg.alphaFDR = 0.05;

% top rows for compact output
cfg.topN = 200;

if ~exist(cfg.outDir, 'dir')
    mkdir(cfg.outDir);
end

diary(fullfile(cfg.outDir, 'STEP85_collect_all_significant_log.txt'));

fprintf('\n=== STEP85 collect all significant STEP84/STEP84B effects ===\n');
fprintf('Started: %s\n', datestr(now));
fprintf('FDR alpha = %.3f\n', cfg.alphaFDR);

%% ===================== INPUT FILES =====================

inputs = {};

% STEP84: paired/repeated-measures version
inputs{end+1,1} = local_input('STEP84', 'SignedRank_binary', 'ChannelFeatureGroup', ...
    fullfile(cfg.step84Dir, 'STEP84_binary_channel_featuregroup_stats.csv'));

inputs{end+1,1} = local_input('STEP84', 'Friedman_3class', 'ChannelFeatureGroup', ...
    fullfile(cfg.step84Dir, 'STEP84_threeclass_channel_featuregroup_stats.csv'));

inputs{end+1,1} = local_input('STEP84', 'SignedRank_binary', 'EdgeFeatureGroup', ...
    fullfile(cfg.step84Dir, 'STEP84_binary_edge_featuregroup_stats.csv'));

inputs{end+1,1} = local_input('STEP84', 'Friedman_3class', 'EdgeFeatureGroup', ...
    fullfile(cfg.step84Dir, 'STEP84_threeclass_edge_featuregroup_stats.csv'));

% STEP84B: requested KW/ranksum version
inputs{end+1,1} = local_input('STEP84B', 'Ranksum_binary', 'ChannelFeatureGroup', ...
    fullfile(cfg.step84BDir, 'STEP84B_binary_channel_featuregroup_RANKSUM_stats.csv'));

inputs{end+1,1} = local_input('STEP84B', 'KruskalWallis_3class', 'ChannelFeatureGroup', ...
    fullfile(cfg.step84BDir, 'STEP84B_threeclass_channel_featuregroup_KW_stats.csv'));

inputs{end+1,1} = local_input('STEP84B', 'Ranksum_binary', 'EdgeFeatureGroup', ...
    fullfile(cfg.step84BDir, 'STEP84B_binary_edge_featuregroup_RANKSUM_stats.csv'));

inputs{end+1,1} = local_input('STEP84B', 'KruskalWallis_3class', 'EdgeFeatureGroup', ...
    fullfile(cfg.step84BDir, 'STEP84B_threeclass_edge_featuregroup_KW_stats.csv'));

%% ===================== READ AND STANDARDIZE =====================

allSig = {};
readSummary = {};

for i = 1:numel(inputs)

    I = inputs{i};

    fprintf('\n[%d/%d] %s | %s | %s\n', i, numel(inputs), I.SourceStep, I.TestFamily, I.UnitType);
    fprintf('File: %s\n', I.File);

    if ~exist(I.File, 'file')
        warning('File not found. Skipping.');
        readSummary{end+1,1} = table(string(I.SourceStep), string(I.TestFamily), string(I.UnitType), string(I.File), ...
            0, 0, string("missing"), ...
            'VariableNames', {'SourceStep','TestFamily','UnitType','File','N_rows','N_sig','Status'}); %#ok<AGROW>
        continue;
    end

    T = readtable(I.File);

    if isempty(T)
        readSummary{end+1,1} = table(string(I.SourceStep), string(I.TestFamily), string(I.UnitType), string(I.File), ...
            0, 0, string("empty"), ...
            'VariableNames', {'SourceStep','TestFamily','UnitType','File','N_rows','N_sig','Status'}); %#ok<AGROW>
        continue;
    end

    if ~any(strcmp(T.Properties.VariableNames, 'q_FDR'))
        warning('No q_FDR column. Skipping.');
        readSummary{end+1,1} = table(string(I.SourceStep), string(I.TestFamily), string(I.UnitType), string(I.File), ...
            height(T), 0, string("no_q_FDR"), ...
            'VariableNames', {'SourceStep','TestFamily','UnitType','File','N_rows','N_sig','Status'}); %#ok<AGROW>
        continue;
    end

    sig = isfinite(T.q_FDR) & T.q_FDR < cfg.alphaFDR;
    Tsig = T(sig, :);

    fprintf('Rows=%d | Significant=%d\n', height(T), height(Tsig));

    readSummary{end+1,1} = table(string(I.SourceStep), string(I.TestFamily), string(I.UnitType), string(I.File), ...
        height(T), height(Tsig), string("ok"), ...
        'VariableNames', {'SourceStep','TestFamily','UnitType','File','N_rows','N_sig','Status'}); %#ok<AGROW>

    if isempty(Tsig)
        continue;
    end

    S = local_standardize(Tsig, I.SourceStep, I.TestFamily, I.UnitType, I.File);

    allSig{end+1,1} = S; %#ok<AGROW>
end

if isempty(allSig)
    Combined = table();
else
    Combined = vertcat(allSig{:});
end

ReadSummary = vertcat(readSummary{:});

%% ===================== SORT / SUMMARIZE =====================

if ~isempty(Combined)

    Combined.AbsEffect = abs(Combined.EffectValue);

    Combined = sortrows(Combined, ...
        {'SourceStep','q_FDR','AbsEffect'}, ...
        {'ascend','ascend','descend'});

    TopSignificant = sortrows(Combined, {'q_FDR','AbsEffect'}, {'ascend','descend'});
    TopSignificant = TopSignificant(1:min(cfg.topN,height(TopSignificant)), :);

    SummaryByCondition = local_summary_by_condition(Combined);
    SummaryByUnitFamily = local_summary_by_unit_family(Combined);

else
    TopSignificant = table();
    SummaryByCondition = table();
    SummaryByUnitFamily = table();
end

%% ===================== SAVE =====================

outCombined = fullfile(cfg.outDir, 'STEP85_all_significant_units_combined.csv');
outTop = fullfile(cfg.outDir, 'STEP85_top_significant_units.csv');
outCond = fullfile(cfg.outDir, 'STEP85_significant_summary_by_condition.csv');
outUnit = fullfile(cfg.outDir, 'STEP85_significant_summary_by_unit_family.csv');
outRead = fullfile(cfg.outDir, 'STEP85_input_read_summary.csv');
outMAT = fullfile(cfg.outDir, 'STEP85_all_significant_results.mat');

writetable(Combined, outCombined);
writetable(TopSignificant, outTop);
writetable(SummaryByCondition, outCond);
writetable(SummaryByUnitFamily, outUnit);
writetable(ReadSummary, outRead);

save(outMAT, 'Combined', 'TopSignificant', 'SummaryByCondition', 'SummaryByUnitFamily', 'ReadSummary', 'cfg', '-v7.3');

fprintf('\n=== Saved STEP85 outputs ===\n');
fprintf('Combined significant units:\n  %s\n', outCombined);
fprintf('Top significant units:\n  %s\n', outTop);
fprintf('Summary by condition:\n  %s\n', outCond);
fprintf('Summary by unit/family:\n  %s\n', outUnit);
fprintf('Input read summary:\n  %s\n', outRead);
fprintf('MAT:\n  %s\n', outMAT);

fprintf('\nFinished: %s\n', datestr(now));
diary off;

%% ========================================================================
%% FUNCTIONS
%% ========================================================================

function I = local_input(sourceStep, testFamily, unitType, filePath)
    I = struct();
    I.SourceStep = sourceStep;
    I.TestFamily = testFamily;
    I.UnitType = unitType;
    I.File = filePath;
end

function S = local_standardize(T, sourceStep, testFamily, unitType, sourceFile)

    n = height(T);

    SourceStep = repmat(string(sourceStep), n, 1);
    TestFamily = repmat(string(testFamily), n, 1);
    UnitType = repmat(string(unitType), n, 1);
    SourceFile = repmat(string(sourceFile), n, 1);

    % Required/common fields
    Contrast = local_get_string_col(T, 'Contrast', n, "unknown");
    Run = local_get_numeric_col(T, 'Run', n, NaN);
    Phase = local_get_string_col(T, 'Phase', n, "unknown");

    UnitID = local_get_string_col(T, 'UnitID', n, "");
    Level = local_get_string_col(T, 'Level', n, "");
    FeatureGroup = local_get_string_col(T, 'FeatureGroup', n, "unknown");

    Channel = local_get_numeric_col(T, 'Channel', n, NaN);
    ChannelLabel = local_get_string_col(T, 'ChannelLabel', n, "");

    Edge = local_get_string_col(T, 'Edge', n, "");
    EdgeLabel = local_get_string_col(T, 'EdgeLabel', n, "");

    N_features_in_unit = local_get_numeric_col(T, 'N_features_in_unit', n, NaN);
    q_FDR = T.q_FDR;

    % p-value column depends on test
    p_value = NaN(n,1);

    if any(strcmp(T.Properties.VariableNames, 'p_signrank'))
        p_value = T.p_signrank;
    elseif any(strcmp(T.Properties.VariableNames, 'p_friedman'))
        p_value = T.p_friedman;
    elseif any(strcmp(T.Properties.VariableNames, 'p_ranksum'))
        p_value = T.p_ranksum;
    elseif any(strcmp(T.Properties.VariableNames, 'p_KW'))
        p_value = T.p_KW;
    end

    % effect column depends on test
    EffectValue = NaN(n,1);
    EffectName = strings(n,1);

    if any(strcmp(T.Properties.VariableNames, 'RankBiserial_paired'))
        EffectValue = T.RankBiserial_paired;
        EffectName(:) = "RankBiserial_paired";
    elseif any(strcmp(T.Properties.VariableNames, 'KendallW'))
        EffectValue = T.KendallW;
        EffectName(:) = "KendallW";
    elseif any(strcmp(T.Properties.VariableNames, 'RankBiserial_A_vs_B'))
        EffectValue = T.RankBiserial_A_vs_B;
        EffectName(:) = "RankBiserial_A_vs_B";
    elseif any(strcmp(T.Properties.VariableNames, 'Epsilon2_KW'))
        EffectValue = T.Epsilon2_KW;
        EffectName(:) = "Epsilon2_KW";
    end

    Direction = strings(n,1);

    for i = 1:n
        if contains(string(testFamily), "binary", 'IgnoreCase', true) || contains(string(testFamily), "Ranksum", 'IgnoreCase', true) || contains(string(testFamily), "SignedRank", 'IgnoreCase', true)
            if isfinite(EffectValue(i))
                if EffectValue(i) > 0
                    Direction(i) = "A_greater_than_B";
                elseif EffectValue(i) < 0
                    Direction(i) = "B_greater_than_A";
                else
                    Direction(i) = "mixed_or_zero";
                end
            else
                Direction(i) = "unknown";
            end
        else
            Direction(i) = "three_class_no_binary_direction";
        end
    end

    % Optional medians
    Median_A = local_get_numeric_col(T, 'Median_A_subjectAgg', n, NaN);
    Median_B = local_get_numeric_col(T, 'Median_B_subjectAgg', n, NaN);
    Median_Diff_AminusB = local_get_numeric_col(T, 'Median_Diff_AminusB', n, NaN);

    Median_color = local_get_numeric_col(T, 'Median_color_subjectAgg', n, NaN);
    Median_orientation = local_get_numeric_col(T, 'Median_orientation_subjectAgg', n, NaN);
    Median_conjunction = local_get_numeric_col(T, 'Median_conjunction_subjectAgg', n, NaN);

    % Optional N
    N_subjects = NaN(n,1);

    if any(strcmp(T.Properties.VariableNames, 'N_subjects'))
        N_subjects = T.N_subjects;
    elseif any(strcmp(T.Properties.VariableNames, 'N_subjects_per_class'))
        N_subjects = T.N_subjects_per_class;
    elseif any(strcmp(T.Properties.VariableNames, 'N_subjects_feature'))
        N_subjects = T.N_subjects_feature;
    end

    S = table(SourceStep, TestFamily, UnitType, Contrast, Run, Phase, ...
        UnitID, Level, FeatureGroup, ...
        Channel, ChannelLabel, Edge, EdgeLabel, ...
        N_subjects, N_features_in_unit, ...
        p_value, q_FDR, EffectName, EffectValue, Direction, ...
        Median_A, Median_B, Median_Diff_AminusB, ...
        Median_color, Median_orientation, Median_conjunction, ...
        SourceFile);
end

function x = local_get_string_col(T, name, n, defaultVal)

    if any(strcmp(T.Properties.VariableNames, name))
        x = string(T.(name));
    else
        x = repmat(string(defaultVal), n, 1);
    end
end

function x = local_get_numeric_col(T, name, n, defaultVal)

    if any(strcmp(T.Properties.VariableNames, name))
        x = double(T.(name));
    else
        x = repmat(defaultVal, n, 1);
    end
end

function Summary = local_summary_by_condition(T)

    [G, keys] = findgroups(T(:, {'SourceStep','TestFamily','UnitType','Contrast','Run','Phase'}));

    N_sig = splitapply(@numel, T.q_FDR, G);
    Min_q = splitapply(@(x) min(x, [], 'omitnan'), T.q_FDR, G);
    Min_p = splitapply(@(x) min(x, [], 'omitnan'), T.p_value, G);
    MaxAbsEffect = splitapply(@(x) max(abs(x), [], 'omitnan'), T.EffectValue, G);
    MeanAbsEffect = splitapply(@(x) mean(abs(x), 'omitnan'), T.EffectValue, G);

    Summary = [keys, table(N_sig, Min_p, Min_q, MeanAbsEffect, MaxAbsEffect)];

    Summary = sortrows(Summary, {'N_sig','Min_q'}, {'descend','ascend'});
end

function Summary = local_summary_by_unit_family(T)

    [G, keys] = findgroups(T(:, {'SourceStep','TestFamily','UnitType','Contrast','Run','Phase','FeatureGroup'}));

    N_sig = splitapply(@numel, T.q_FDR, G);
    Min_q = splitapply(@(x) min(x, [], 'omitnan'), T.q_FDR, G);
    Min_p = splitapply(@(x) min(x, [], 'omitnan'), T.p_value, G);
    MaxAbsEffect = splitapply(@(x) max(abs(x), [], 'omitnan'), T.EffectValue, G);
    MeanAbsEffect = splitapply(@(x) mean(abs(x), 'omitnan'), T.EffectValue, G);

    Summary = [keys, table(N_sig, Min_p, Min_q, MeanAbsEffect, MaxAbsEffect)];

    Summary = sortrows(Summary, {'N_sig','Min_q'}, {'descend','ascend'});
end

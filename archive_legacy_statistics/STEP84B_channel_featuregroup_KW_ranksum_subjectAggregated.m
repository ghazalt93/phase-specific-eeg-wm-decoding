%% STEP84B_channel_featuregroup_KW_ranksum_subjectAggregated.m

cfg = get_project_config();

clear; clc;

%% ===================== CONFIG =====================

cfg = struct();

cfg.ROOT = cfg.outputRoot;
cfg.datasetPath = '';

cfg.candidatePaths = { ...
    fullfile(cfg.ROOT,'_wm_ml','dataset.mat')
    fullfile(cfg.ROOT,'Subjects','_wm_ml','dataset.mat')
    fullfile(cfg.ROOT,'_wm_dataset','dataset.mat')
    fullfile(cfg.ROOT,'dataset.mat')
};

cfg.montageFile = fullfile(cfg.ROOT,'ipm2.ced.tsv');
cfg.outDir = fullfile(cfg.ROOT,'_wm_STEP84B_channel_featuregroup_KW_ranksum_subjectAggregated');

cfg.runList = [1 2 3];
cfg.phaseList = {'stim','maint','retr'};

cfg.trialAggregate = 'median';    % aggregate trials within Subject x Condition
cfg.featureAggregate = 'median';  % aggregate features within Channel/Edge x FeatureGroup

cfg.alphaFDR = 0.05;
cfg.minSubjectsPerClass = 8;

cfg.doChannelFeatureGroups = true;
cfg.doConnectivityEdges = true;

cfg.labelMap.color       = {'1','color'};
cfg.labelMap.orientation = {'2','orientation'};
cfg.labelMap.conjunction = {'3','conjunction'};

contrasts = struct([]);
contrasts(1).name = 'orientation_vs_conjunction';
contrasts(1).A = cfg.labelMap.orientation;
contrasts(1).B = cfg.labelMap.conjunction;

contrasts(2).name = 'color_vs_orientation';
contrasts(2).A = cfg.labelMap.color;
contrasts(2).B = cfg.labelMap.orientation;

contrasts(3).name = 'color_vs_conjunction';
contrasts(3).A = cfg.labelMap.color;
contrasts(3).B = cfg.labelMap.conjunction;

%% ===================== PREPARE =====================

if ~exist(cfg.outDir,'dir')
    mkdir(cfg.outDir);
end

diary(fullfile(cfg.outDir,'STEP84B_KW_ranksum_subjectAggregated_log.txt'));

fprintf('\n=== STEP84B KW/ranksum on subject-aggregated Channel/Edge x FeatureGroup ===\n');
fprintf('Started: %s\n', datestr(now));
fprintf('ROOT: %s\n', cfg.ROOT);
fprintf('trialAggregate=%s | featureAggregate=%s\n', cfg.trialAggregate, cfg.featureAggregate);

%% ===================== LOAD DATA =====================

[T,col,featureNames] = local_load_data(cfg);

subjVals  = cleanstr(T.(col.subj));
condVals  = cleanstr(T.(col.cond));
runVals   = double(T.(col.run));
phaseVals = cleanstr(T.(col.phase));

labelMap = load_montage(cfg.montageFile);

fprintf('Loaded data: rows=%d | features=%d | subjects=%d\n', ...
    height(T), numel(featureNames), numel(unique(subjVals)));

%% ===================== PARSE FEATURES =====================

Finfo = parse_feature_info(featureNames, labelMap);

fprintf('\nFeature groups detected:\n');
disp(groupsummary(Finfo,'FeatureGroup'));

fprintf('\nFeature level detected:\n');
disp(groupsummary(Finfo,'Level'));

%% ===================== MAIN LOOP =====================

allRS_Channel = {};
allKW_Channel = {};
allRS_Edge = {};
allKW_Edge = {};
allJobs = {};

for ir = 1:numel(cfg.runList)

    runNum = cfg.runList(ir);

    for ip = 1:numel(cfg.phaseList)

        phaseName = cfg.phaseList{ip};
        idxBase = runVals == runNum & phaseVals == cleanstr({phaseName});

        fprintf('\n\n====================================================\n');
        fprintf('Run %d | Phase %s | rows=%d\n', runNum, phaseName, sum(idxBase));
        fprintf('====================================================\n');

        if sum(idxBase) < 20
            fprintf('Too few rows. Skipping.\n');
            continue;
        end

        %% ---------- Channel x FeatureGroup ----------
        if cfg.doChannelFeatureGroups

            UnitsCh = build_channel_units(Finfo);
            fprintf('Channel x FeatureGroup units: %d\n', height(UnitsCh));

            if height(UnitsCh) > 0

                % 3-class Kruskal-Wallis
                [X3, subj3] = build_subject_condition_unit_matrix( ...
                    T, idxBase, subjVals, condVals, featureNames, UnitsCh, ...
                    {cfg.labelMap.color,cfg.labelMap.orientation,cfg.labelMap.conjunction}, cfg);

                if numel(subj3) >= cfg.minSubjectsPerClass
                    R = run_kruskalwallis_units(X3, UnitsCh);
                    R.Test = repmat("KruskalWallis_3class_ChannelFeatureGroup",height(R),1);
                    R.Contrast = repmat("color_vs_orientation_vs_conjunction",height(R),1);
                    R.Run = repmat(runNum,height(R),1);
                    R.Phase = repmat(string(phaseName),height(R),1);
                    R.N_subjects_per_class = repmat(numel(subj3),height(R),1);
                    R = movevars(R,{'Test','Contrast','Run','Phase','N_subjects_per_class'},'Before',1);

                    allKW_Channel{end+1,1} = R; %#ok<AGROW>

                    nSig = sum(R.q_FDR < cfg.alphaFDR,'omitnan');
                    minq = min(R.q_FDR,[],'omitnan');
                    maxe = max(R.Epsilon2_KW,[],'omitnan');

                    allJobs{end+1,1}=table(string("KruskalWallis_3class_ChannelFeatureGroup"),string("color_vs_orientation_vs_conjunction"),runNum,string(phaseName),numel(subj3),height(R),nSig,minq,maxe, ...
                        'VariableNames',{'Test','Contrast','Run','Phase','N_subjects_per_class','N_tests','N_sig_FDR','Min_q','MaxEffect'}); %#ok<AGROW>

                    save_unit_outputs(cfg,'ChannelFeatureGroup','KruskalWallis_3class','color_vs_orientation_vs_conjunction',runNum,phaseName,R);
                    fprintf('ChannelGroup KW: subjects/class=%d | Nsig=%d | min q=%.3g\n', numel(subj3), nSig, minq);
                end

                % Binary rank-sum
                for ic = 1:numel(contrasts)

                    C = contrasts(ic);

                    [X2, subj2] = build_subject_condition_unit_matrix( ...
                        T, idxBase, subjVals, condVals, featureNames, UnitsCh, ...
                        {C.A,C.B}, cfg);

                    if numel(subj2) < cfg.minSubjectsPerClass
                        fprintf('ChannelGroup ranksum %s skipped: subjects/class=%d\n', C.name, numel(subj2));
                        continue;
                    end

                    R = run_ranksum_units(X2, UnitsCh);
                    R.Test = repmat("Ranksum_binary_ChannelFeatureGroup",height(R),1);
                    R.Contrast = repmat(string(C.name),height(R),1);
                    R.Run = repmat(runNum,height(R),1);
                    R.Phase = repmat(string(phaseName),height(R),1);
                    R.N_subjects_per_class = repmat(numel(subj2),height(R),1);
                    R = movevars(R,{'Test','Contrast','Run','Phase','N_subjects_per_class'},'Before',1);

                    allRS_Channel{end+1,1} = R; %#ok<AGROW>

                    nSig = sum(R.q_FDR < cfg.alphaFDR,'omitnan');
                    minq = min(R.q_FDR,[],'omitnan');
                    maxe = max(abs(R.RankBiserial_A_vs_B),[],'omitnan');

                    allJobs{end+1,1}=table(string("Ranksum_binary_ChannelFeatureGroup"),string(C.name),runNum,string(phaseName),numel(subj2),height(R),nSig,minq,maxe, ...
                        'VariableNames',{'Test','Contrast','Run','Phase','N_subjects_per_class','N_tests','N_sig_FDR','Min_q','MaxEffect'}); %#ok<AGROW>

                    save_unit_outputs(cfg,'ChannelFeatureGroup','Ranksum_binary',C.name,runNum,phaseName,R);
                    fprintf('ChannelGroup ranksum %s: subjects/class=%d | Nsig=%d | min q=%.3g\n', C.name, numel(subj2), nSig, minq);
                end
            end
        end

        %% ---------- Edge x ConnectivityGroup ----------
        if cfg.doConnectivityEdges

            UnitsEdge = build_edge_units(Finfo);
            fprintf('Connectivity Edge units: %d\n', height(UnitsEdge));

            if height(UnitsEdge) > 0

                % 3-class Kruskal-Wallis
                [X3e, subj3e] = build_subject_condition_unit_matrix( ...
                    T, idxBase, subjVals, condVals, featureNames, UnitsEdge, ...
                    {cfg.labelMap.color,cfg.labelMap.orientation,cfg.labelMap.conjunction}, cfg);

                if numel(subj3e) >= cfg.minSubjectsPerClass
                    R = run_kruskalwallis_units(X3e, UnitsEdge);
                    R.Test = repmat("KruskalWallis_3class_EdgeFeatureGroup",height(R),1);
                    R.Contrast = repmat("color_vs_orientation_vs_conjunction",height(R),1);
                    R.Run = repmat(runNum,height(R),1);
                    R.Phase = repmat(string(phaseName),height(R),1);
                    R.N_subjects_per_class = repmat(numel(subj3e),height(R),1);
                    R = movevars(R,{'Test','Contrast','Run','Phase','N_subjects_per_class'},'Before',1);

                    allKW_Edge{end+1,1} = R; %#ok<AGROW>

                    nSig = sum(R.q_FDR < cfg.alphaFDR,'omitnan');
                    minq = min(R.q_FDR,[],'omitnan');
                    maxe = max(R.Epsilon2_KW,[],'omitnan');

                    allJobs{end+1,1}=table(string("KruskalWallis_3class_EdgeFeatureGroup"),string("color_vs_orientation_vs_conjunction"),runNum,string(phaseName),numel(subj3e),height(R),nSig,minq,maxe, ...
                        'VariableNames',{'Test','Contrast','Run','Phase','N_subjects_per_class','N_tests','N_sig_FDR','Min_q','MaxEffect'}); %#ok<AGROW>

                    save_unit_outputs(cfg,'EdgeFeatureGroup','KruskalWallis_3class','color_vs_orientation_vs_conjunction',runNum,phaseName,R);
                    fprintf('EdgeGroup KW: subjects/class=%d | Nsig=%d | min q=%.3g\n', numel(subj3e), nSig, minq);
                end

                % Binary rank-sum
                for ic = 1:numel(contrasts)

                    C = contrasts(ic);

                    [X2e, subj2e] = build_subject_condition_unit_matrix( ...
                        T, idxBase, subjVals, condVals, featureNames, UnitsEdge, ...
                        {C.A,C.B}, cfg);

                    if numel(subj2e) < cfg.minSubjectsPerClass
                        fprintf('EdgeGroup ranksum %s skipped: subjects/class=%d\n', C.name, numel(subj2e));
                        continue;
                    end

                    R = run_ranksum_units(X2e, UnitsEdge);
                    R.Test = repmat("Ranksum_binary_EdgeFeatureGroup",height(R),1);
                    R.Contrast = repmat(string(C.name),height(R),1);
                    R.Run = repmat(runNum,height(R),1);
                    R.Phase = repmat(string(phaseName),height(R),1);
                    R.N_subjects_per_class = repmat(numel(subj2e),height(R),1);
                    R = movevars(R,{'Test','Contrast','Run','Phase','N_subjects_per_class'},'Before',1);

                    allRS_Edge{end+1,1} = R; %#ok<AGROW>

                    nSig = sum(R.q_FDR < cfg.alphaFDR,'omitnan');
                    minq = min(R.q_FDR,[],'omitnan');
                    maxe = max(abs(R.RankBiserial_A_vs_B),[],'omitnan');

                    allJobs{end+1,1}=table(string("Ranksum_binary_EdgeFeatureGroup"),string(C.name),runNum,string(phaseName),numel(subj2e),height(R),nSig,minq,maxe, ...
                        'VariableNames',{'Test','Contrast','Run','Phase','N_subjects_per_class','N_tests','N_sig_FDR','Min_q','MaxEffect'}); %#ok<AGROW>

                    save_unit_outputs(cfg,'EdgeFeatureGroup','Ranksum_binary',C.name,runNum,phaseName,R);
                    fprintf('EdgeGroup ranksum %s: subjects/class=%d | Nsig=%d | min q=%.3g\n', C.name, numel(subj2e), nSig, minq);
                end
            end
        end
    end
end

%% ===================== SAVE GLOBAL =====================

RS_Channel = cat_tables(allRS_Channel);
KW_Channel = cat_tables(allKW_Channel);
RS_Edge = cat_tables(allRS_Edge);
KW_Edge = cat_tables(allKW_Edge);

if isempty(allJobs)
    JobSummary = table();
else
    JobSummary = sortrows(vertcat(allJobs{:}), {'N_sig_FDR','Min_q'}, {'descend','ascend'});
end

writetable(RS_Channel, fullfile(cfg.outDir,'STEP84B_binary_channel_featuregroup_RANKSUM_stats.csv'));
writetable(KW_Channel, fullfile(cfg.outDir,'STEP84B_threeclass_channel_featuregroup_KW_stats.csv'));
writetable(RS_Edge, fullfile(cfg.outDir,'STEP84B_binary_edge_featuregroup_RANKSUM_stats.csv'));
writetable(KW_Edge, fullfile(cfg.outDir,'STEP84B_threeclass_edge_featuregroup_KW_stats.csv'));
writetable(JobSummary, fullfile(cfg.outDir,'STEP84B_job_summary.csv'));

save(fullfile(cfg.outDir,'STEP84B_KW_ranksum_subjectAggregated.mat'), ...
    'RS_Channel','KW_Channel','RS_Edge','KW_Edge','JobSummary','cfg','-v7.3');

fprintf('\n=== STEP84B finished ===\n');
fprintf('Main summary:\n%s\n', fullfile(cfg.outDir,'STEP84B_job_summary.csv'));
fprintf('Finished: %s\n', datestr(now));

diary off;

%% ========================================================================
%% FUNCTIONS
%% ========================================================================

function [T,col,featureNames] = local_load_data(cfg)

    datasetPath = cfg.datasetPath;

    if isempty(datasetPath)
        for i=1:numel(cfg.candidatePaths)
            if exist(cfg.candidatePaths{i},'file')
                datasetPath = cfg.candidatePaths{i};
                break;
            end
        end
    end

    if isempty(datasetPath)
        error('dataset.mat not found.');
    end

    fprintf('Loading dataset:\n%s\n', datasetPath);

    S = load(datasetPath);

    if isfield(S,'DS') && istable(S.DS)
        T = S.DS;
    elseif isfield(S,'T') && istable(S.T)
        T = S.T;
    else
        error('Could not find table DS or T in dataset.mat');
    end

    names = T.Properties.VariableNames;

    col.subj  = findvar(names, {'Subject','subject','SubjectID','subjectID','subj','subjID','Subj','SubjID'});
    col.run   = findvar(names, {'runNum','RunNum','run','Run'});
    col.phase = findvar(names, {'phase','Phase'});
    col.cond  = findvar(names, {'Condition','condition','yCondition','Label','label'});

    if isempty(col.subj) || isempty(col.run) || isempty(col.phase) || isempty(col.cond)
        error('Missing metadata columns.');
    end

    featureNames = detect_feature_columns(T);
end

function v = findvar(names,candidates)
    v = '';
    for i=1:numel(candidates)
        idx = strcmp(names,candidates{i});
        if any(idx)
            v = names{find(idx,1)};
            return;
        end
    end
end

function featureNames = detect_feature_columns(T)

    names = T.Properties.VariableNames;
    isNum = false(1,numel(names));

    for j=1:numel(names)
        x = T.(names{j});
        isNum(j) = isnumeric(x) || islogical(x);
    end

    numNames = names(isNum);

    meta = {'Subject','subject','SubjectID','subjectID','subj','subjID','Subj','SubjID', ...
        'Run','run','runNum','RunNum','session','Session','Phase','phase', ...
        'Condition','condition','yCondition','Label','label','Correct','correct','yCorrect', ...
        'TrialNum','trialNum','TrialIndex','trialIndex','Trial','trial','PatternID','patternID', ...
        'StartRow','EndRow','Second10Row','RetrRow','Fold','fold','CVFold'};

    isMeta = false(size(numNames));

    for j=1:numel(numNames)
        isMeta(j) = any(strcmpi(numNames{j}, meta));
    end

    featureNames = numNames(~isMeta);
end

function s = cleanstr(x)

    if isnumeric(x) || islogical(x)
        s = string(x);
    elseif iscell(x)
        s = string(x);
    elseif iscategorical(x)
        s = string(x);
    elseif isstring(x)
        s = x;
    elseif ischar(x)
        s = string(cellstr(x));
    else
        s = string(x);
    end

    s = lower(strtrim(s));
    s = regexprep(s, '\s+', '');
    s = regexprep(s, '[^\w]', '');
end

function labelMap = load_montage(file)

    labelMap = containers.Map('KeyType','double','ValueType','char');

    if ~exist(file,'file')
        warning('Montage file not found. Labels will use Cxx.');
        return;
    end

    try
        M = readtable(file,'FileType','text','Delimiter','\t');
        vars = M.Properties.VariableNames;

        if any(strcmp(vars,'Number')), ncol='Number'; else, ncol=vars{1}; end
        if any(strcmp(vars,'labels')), lcol='labels'; else, lcol=vars{2}; end

        nums = double(M.(ncol));
        labs = string(M.(lcol));

        for i=1:numel(nums)
            if isfinite(nums(i))
                labelMap(nums(i)) = char(labs(i));
            end
        end
    catch ME
        warning('Could not read montage: %s', ME.message);
    end
end

function out = lab(labelMap,ch)
    if isKey(labelMap,double(ch))
        out = string(labelMap(double(ch)));
    else
        out = "C" + compose("%02d",ch);
    end
end

function Finfo = parse_feature_info(featureNames,labelMap)

    n = numel(featureNames);

    Feature = string(featureNames(:));
    FeatureGroup = strings(n,1);
    Level = strings(n,1);
    Channel = NaN(n,1);
    ChannelLabel = strings(n,1);
    ChannelA = NaN(n,1);
    ChannelB = NaN(n,1);
    Edge = strings(n,1);
    EdgeLabel = strings(n,1);

    for i=1:n
        f = lower(char(Feature(i)));

        fg = infer_group(f);
        FeatureGroup(i) = fg;

        chs = parse_channels(f);

        if fg == "connectivity" && numel(chs) >= 2
            Level(i) = "edge";
            a = chs(1); b = chs(2);
            ChannelA(i) = a;
            ChannelB(i) = b;
            aa = min(a,b); bb = max(a,b);
            Edge(i) = sprintf('C%02d-C%02d',aa,bb);
            EdgeLabel(i) = lab(labelMap,a) + "-" + lab(labelMap,b);
        elseif numel(chs) >= 1
            Level(i) = "channel";
            Channel(i) = chs(1);
            ChannelLabel(i) = lab(labelMap,chs(1));
        else
            Level(i) = "unknown";
        end
    end

    Finfo = table(Feature,FeatureGroup,Level,Channel,ChannelLabel,ChannelA,ChannelB,Edge,EdgeLabel);
end

function fg = infer_group(f)

    s = string(f);

    if contains(s,"rie") || contains(s,"pli") || contains(s,"plv") || contains(s,"coh") || contains(s,"conn") || ~isempty(regexp(f,'c0?\d{1,2}[_-]c0?\d{1,2}','once'))
        fg = "connectivity";
    elseif contains(s,"rbp") || contains(s,"relative")
        fg = "relative_bandpower";
    elseif contains(s,"bp_") || contains(s,"bandpower") || contains(s,"alpha") || contains(s,"beta") || contains(s,"theta") || contains(s,"delta") || contains(s,"gamma")
        fg = "bandpower";
    elseif contains(s,"morlet") || contains(s,"tf_") || contains(s,"wavelet") || contains(s,"tfr")
        fg = "time_frequency";
    elseif contains(s,"lz") || contains(s,"lzc") || contains(s,"entropy") || contains(s,"samp") || contains(s,"perm") || contains(s,"higuchi") || contains(s,"fractal")
        fg = "complexity";
    elseif contains(s,"skew") || contains(s,"kurt") || contains(s,"hjorth") || contains(s,"line") || contains(s,"auc") || contains(s,"rms") || contains(s,"mean") || contains(s,"std") || contains(s,"var") || contains(s,"max") || contains(s,"min")
        fg = "temporal_statistical";
    else
        fg = "unknown";
    end
end

function chs = parse_channels(f)

    chs = [];

    tok = regexp(f,'c0?(\d{1,2})[_-]c0?(\d{1,2})','tokens');

    if ~isempty(tok)
        chs = [str2double(tok{1}{1}), str2double(tok{1}{2})];
    else
        tok = regexp(f,'ch0?(\d{1,2})[_-]ch0?(\d{1,2})','tokens');
        if ~isempty(tok)
            chs = [str2double(tok{1}{1}), str2double(tok{1}{2})];
        end
    end

    if isempty(chs)
        toks = regexp(f,'(?:^|[_\-\s])c(?:h)?0?(\d{1,2})(?=[_\-\s]|$)','tokens');
        vals = [];
        for k=1:numel(toks)
            v = str2double(toks{k}{1});
            if isfinite(v) && v>=1 && v<=128
                vals(end+1) = v; %#ok<AGROW>
            end
        end
        if ~isempty(vals)
            chs = vals;
        end
    end

    chs = unique(chs(isfinite(chs) & chs>=1 & chs<=128),'stable');
end

function Units = build_channel_units(Finfo)

    idx = Finfo.Level == "channel" & isfinite(Finfo.Channel) & Finfo.FeatureGroup ~= "connectivity";

    T = Finfo(idx,:);

    if isempty(T)
        Units = table();
        return;
    end

    [G,ch,fg] = findgroups(T.Channel,T.FeatureGroup);

    FeatureIdx = splitapply(@(x){x}, find(idx), G);
    UnitID = "CH" + compose("%02d",ch) + "_" + string(fg);

    Channel = ch;
    FeatureGroup = string(fg);
    ChannelLabel = strings(numel(ch),1);

    for i=1:numel(ch)
        ii = find(T.Channel == ch(i),1);
        ChannelLabel(i) = T.ChannelLabel(ii);
    end

    Level = repmat("channel_featuregroup",numel(ch),1);
    Edge = strings(numel(ch),1);
    EdgeLabel = strings(numel(ch),1);

    Units = table(UnitID,Level,Channel,ChannelLabel,Edge,EdgeLabel,FeatureGroup,FeatureIdx);
end

function Units = build_edge_units(Finfo)

    idx = Finfo.Level == "edge" & strlength(Finfo.Edge)>0;

    T = Finfo(idx,:);

    if isempty(T)
        Units = table();
        return;
    end

    [G,edge,fg] = findgroups(T.Edge,T.FeatureGroup);

    FeatureIdx = splitapply(@(x){x}, find(idx), G);
    UnitID = string(edge) + "_" + string(fg);

    Level = repmat("edge_featuregroup",numel(edge),1);
    Channel = NaN(numel(edge),1);
    ChannelLabel = strings(numel(edge),1);
    Edge = string(edge);
    EdgeLabel = strings(numel(edge),1);
    FeatureGroup = string(fg);

    for i=1:numel(edge)
        ii = find(T.Edge == edge(i),1);
        EdgeLabel(i) = T.EdgeLabel(ii);
    end

    Units = table(UnitID,Level,Channel,ChannelLabel,Edge,EdgeLabel,FeatureGroup,FeatureIdx);
end

function [Xunit,subjKeep] = build_subject_condition_unit_matrix(T,idxBase,subjVals,condVals,featureNames,Units,classSets,cfg)

    subs = unique(subjVals(idxBase));
    subs = subs(~ismissing(subs));

    nC = numel(classSets);
    nU = height(Units);

    Xlist = {};
    subjKeep = strings(0,1);

    for s=1:numel(subs)
        idxS = idxBase & subjVals == subs(s);
        Xi = NaN(nC,nU);
        okSubj = true;

        for c=1:nC
            idxC = idxS & ismember(condVals, cleanstr(classSets{c}));

            if sum(idxC) < 1
                okSubj = false;
                break;
            end

            XrawAll = double(table2array(T(idxC,featureNames)));
            XrawAll(~isfinite(XrawAll)) = NaN;

            for u=1:nU
                featIdx = Units.FeatureIdx{u};

                if isempty(featIdx)
                    continue;
                end

                Xu = XrawAll(:,featIdx);

                switch lower(string(cfg.featureAggregate))
                    case "mean"
                        perTrial = mean(Xu,2,'omitnan');
                    otherwise
                        perTrial = median(Xu,2,'omitnan');
                end

                switch lower(string(cfg.trialAggregate))
                    case "mean"
                        Xi(c,u) = mean(perTrial,'omitnan');
                    otherwise
                        Xi(c,u) = median(perTrial,'omitnan');
                end
            end
        end

        if okSubj
            Xlist{end+1,1} = Xi; %#ok<AGROW>
            subjKeep(end+1,1) = subs(s); %#ok<AGROW>
        end
    end

    Xunit = NaN(numel(Xlist),nC,nU);

    for s=1:numel(Xlist)
        Xunit(s,:,:) = Xlist{s};
    end
end

function R = run_ranksum_units(Xunit,Units)

    nU = size(Xunit,3);

    p = NaN(nU,1);
    z = NaN(nU,1);
    U_A = NaN(nU,1);
    aucA = NaN(nU,1);
    rb = NaN(nU,1);
    medA = NaN(nU,1);
    medB = NaN(nU,1);
    nA = NaN(nU,1);
    nB = NaN(nU,1);

    for u=1:nU
        A = squeeze(Xunit(:,1,u));
        B = squeeze(Xunit(:,2,u));

        A = A(isfinite(A));
        B = B(isfinite(B));

        nA(u) = numel(A);
        nB(u) = numel(B);

        if numel(A) < 5 || numel(B) < 5
            continue;
        end

        if std([A;B],0,'omitnan') <= eps
            continue;
        end

        try
            [p(u),~,stats] = ranksum(A,B);
            if isfield(stats,'zval'), z(u) = stats.zval; end

            x = [A;B];
            ranks = tiedrank(x);
            na = numel(A);
            nb = numel(B);

            U = sum(ranks(1:na)) - na*(na+1)/2;
            U_A(u) = U;
            aucA(u) = U/(na*nb);
            rb(u) = 2*aucA(u) - 1;

            medA(u) = median(A,'omitnan');
            medB(u) = median(B,'omitnan');
        catch
        end
    end

    q = bh_fdr(p);

    R = Units(:,{'UnitID','Level','Channel','ChannelLabel','Edge','EdgeLabel','FeatureGroup'});
    R.N_features_in_unit = cellfun(@numel,Units.FeatureIdx);
    R.p_ranksum = p;
    R.q_FDR = q;
    R.z_ranksum = z;
    R.U_A = U_A;
    R.AUC_A = aucA;
    R.RankBiserial_A_vs_B = rb;
    R.Median_A_subjectAgg = medA;
    R.Median_B_subjectAgg = medB;
    R.N_A_subjects = nA;
    R.N_B_subjects = nB;
end

function R = run_kruskalwallis_units(Xunit,Units)

    nU = size(Xunit,3);

    p = NaN(nU,1);
    chi2 = NaN(nU,1);
    df = NaN(nU,1);
    eps2 = NaN(nU,1);

    med1 = NaN(nU,1);
    med2 = NaN(nU,1);
    med3 = NaN(nU,1);

    n1 = NaN(nU,1);
    n2 = NaN(nU,1);
    n3 = NaN(nU,1);

    for u=1:nU
        A = squeeze(Xunit(:,1,u));
        B = squeeze(Xunit(:,2,u));
        C = squeeze(Xunit(:,3,u));

        A = A(isfinite(A));
        B = B(isfinite(B));
        C = C(isfinite(C));

        n1(u)=numel(A);
        n2(u)=numel(B);
        n3(u)=numel(C);

        if min([numel(A),numel(B),numel(C)]) < 5
            continue;
        end

        x = [A;B;C];

        if std(x,0,'omitnan') <= eps
            continue;
        end

        g = [repmat("color",numel(A),1); repmat("orientation",numel(B),1); repmat("conjunction",numel(C),1)];
        g = categorical(g);

        try
            [p(u),tbl] = kruskalwallis(x,g,'off');

            chi2(u) = tbl{2,5};
            df(u) = tbl{2,3};

            N = numel(x);
            k = 3;
            eps2(u) = max(0,(chi2(u)-k+1)/max(N-k,1));

            med1(u) = median(A,'omitnan');
            med2(u) = median(B,'omitnan');
            med3(u) = median(C,'omitnan');
        catch
        end
    end

    q = bh_fdr(p);

    R = Units(:,{'UnitID','Level','Channel','ChannelLabel','Edge','EdgeLabel','FeatureGroup'});
    R.N_features_in_unit = cellfun(@numel,Units.FeatureIdx);
    R.p_KW = p;
    R.q_FDR = q;
    R.KW_chi2 = chi2;
    R.KW_df = df;
    R.Epsilon2_KW = eps2;
    R.Median_color_subjectAgg = med1;
    R.Median_orientation_subjectAgg = med2;
    R.Median_conjunction_subjectAgg = med3;
    R.N_color_subjects = n1;
    R.N_orientation_subjects = n2;
    R.N_conjunction_subjects = n3;
end

function q = bh_fdr(p)

    p = double(p(:));
    q = NaN(size(p));

    ok = isfinite(p) & p>=0 & p<=1;
    p0 = p(ok);
    m = numel(p0);

    if m==0
        return;
    end

    [ps,ord] = sort(p0,'ascend');
    qs = ps .* m ./ (1:m)';

    for i=m-1:-1:1
        qs(i) = min(qs(i),qs(i+1));
    end

    qs(qs>1) = 1;

    tmp = NaN(m,1);
    tmp(ord) = qs;
    q(ok) = tmp;
end

function T = cat_tables(cells)

    if isempty(cells)
        T = table();
    else
        T = vertcat(cells{:});
    end
end

function save_unit_outputs(cfg,unitType,testName,contrastName,runNum,phaseName,R)

    sig = isfinite(R.q_FDR) & R.q_FDR < cfg.alphaFDR;
    Rsig = R(sig,:);

    if isempty(Rsig)
        return;
    end

    if any(strcmp(Rsig.Properties.VariableNames,'RankBiserial_A_vs_B'))
        Rsig.AbsEffect = abs(Rsig.RankBiserial_A_vs_B);
        Rsig = sortrows(Rsig, {'q_FDR','AbsEffect'}, {'ascend','descend'});
    elseif any(strcmp(Rsig.Properties.VariableNames,'Epsilon2_KW'))
        Rsig = sortrows(Rsig, {'q_FDR','Epsilon2_KW'}, {'ascend','descend'});
    else
        Rsig = sortrows(Rsig,'q_FDR','ascend');
    end

    tag = regexprep(sprintf('%s_%s_%s_run%d_%s',unitType,testName,contrastName,runNum,phaseName),'[^\w]','_');

    outSub = fullfile(cfg.outDir,tag);

    if ~exist(outSub,'dir')
        mkdir(outSub);
    end

    writetable(Rsig, fullfile(outSub,sprintf('STEP84B_%s_significant_units.csv',tag)));

    if contains(unitType,'Channel')
        writetable(local_channel_unit_summary(Rsig), fullfile(outSub,sprintf('STEP84B_%s_channel_summary.csv',tag)));
    else
        writetable(local_edge_unit_summary(Rsig), fullfile(outSub,sprintf('STEP84B_%s_edge_summary.csv',tag)));
    end
end

function S = local_channel_unit_summary(T)

    if isempty(T)
        S = table();
        return;
    end

    [G,ch,fg] = findgroups(T.Channel,T.FeatureGroup);

    N_sig_units = splitapply(@numel,T.q_FDR,G);
    Min_q = splitapply(@(x)min(x,[],'omitnan'),T.q_FDR,G);

    if any(strcmp(T.Properties.VariableNames,'RankBiserial_A_vs_B'))
        eff = abs(T.RankBiserial_A_vs_B);
    else
        eff = T.Epsilon2_KW;
    end

    MaxEffect = splitapply(@(x)max(x,[],'omitnan'),eff,G);
    MeanEffect = splitapply(@(x)mean(x,'omitnan'),eff,G);

    labels = strings(size(ch));
    for i=1:numel(ch)
        ii = find(T.Channel==ch(i) & T.FeatureGroup==string(fg(i)),1);
        labels(i) = T.ChannelLabel(ii);
    end

    S = table(ch,labels,string(fg),N_sig_units,Min_q,MeanEffect,MaxEffect, ...
        'VariableNames',{'Channel','ChannelLabel','FeatureGroup','N_sig_units','Min_q','MeanEffect','MaxEffect'});

    S = sortrows(S,{'N_sig_units','Min_q'},{'descend','ascend'});
end

function S = local_edge_unit_summary(T)

    if isempty(T)
        S = table();
        return;
    end

    [G,edge,fg] = findgroups(T.Edge,T.FeatureGroup);

    N_sig_units = splitapply(@numel,T.q_FDR,G);
    Min_q = splitapply(@(x)min(x,[],'omitnan'),T.q_FDR,G);

    if any(strcmp(T.Properties.VariableNames,'RankBiserial_A_vs_B'))
        eff = abs(T.RankBiserial_A_vs_B);
        rawEff = T.RankBiserial_A_vs_B;
    else
        eff = T.Epsilon2_KW;
        rawEff = T.Epsilon2_KW;
    end

    MaxEffect = splitapply(@(x)max(x,[],'omitnan'),eff,G);
    MeanEffect = splitapply(@(x)mean(x,'omitnan'),eff,G);
    MeanRawEffect = splitapply(@(x)mean(x,'omitnan'),rawEff,G);

    labels = strings(size(edge));
    direction = strings(size(edge));

    for i=1:numel(edge)
        ii = find(T.Edge==string(edge(i)) & T.FeatureGroup==string(fg(i)),1);
        labels(i) = T.EdgeLabel(ii);

        if any(strcmp(T.Properties.VariableNames,'RankBiserial_A_vs_B'))
            if MeanRawEffect(i)>0
                direction(i)="A_greater_than_B";
            elseif MeanRawEffect(i)<0
                direction(i)="B_greater_than_A";
            else
                direction(i)="mixed_or_zero";
            end
        else
            direction(i)="three_class_no_binary_direction";
        end
    end

    S = table(string(edge),labels,string(fg),N_sig_units,Min_q,MeanEffect,MaxEffect,direction, ...
        'VariableNames',{'Edge','EdgeLabel','FeatureGroup','N_sig_units','Min_q','MeanEffect','MaxEffect','Direction'});

    S = sortrows(S,{'N_sig_units','Min_q'},{'descend','ascend'});
end

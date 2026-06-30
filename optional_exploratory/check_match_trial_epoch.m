%% ===================== WS/EDF QA REPORT (3 tests) =====================
clc; clear; close all;

% --------- PATHS (EDIT THESE) ---------
root_in = '<RESULT_MINE_ROOT>/New Folder/no_ica';
subjects_root = fullfile(root_in,'Subjects');

results_root = fullfile(root_in,'GroupResults_LOSO_new'); 

out_dir = fullfile(results_root,'_debug');
if ~exist(out_dir,'dir'), mkdir(out_dir); end

% ---------- Filters ----------
wsRejectTokens  = ["eye","backup","preproc","filtered","(copy)"," copy"];
edfRejectTokens = ["eye"]; 

% ---------- Run ----------
subj_dirs = dir(fullfile(subjects_root,'s*'));
subj_dirs = subj_dirs([subj_dirs.isdir]);

% numeric sort s1,s2,...
subj_nums = nan(1,numel(subj_dirs));
for k=1:numel(subj_dirs)
    tok = regexp(subj_dirs(k).name,'s(\d+)','tokens','once');
    if ~isempty(tok), subj_nums(k)=str2double(tok{1}); end
end
[~,ord]=sort(subj_nums);
subj_dirs=subj_dirs(ord);

rows = {};  % will become table at end

fprintf('=== WS/EDF QA REPORT ===\nSubjects root: %s\nNsubjects=%d\n', subjects_root, numel(subj_dirs));

for isub = 1:numel(subj_dirs)
    subj = string(subj_dirs(isub).name);
    subj_dir = fullfile(subjects_root, subj);
    
    wsList = dir(fullfile(subj_dir,'workspace_*.mat'));
    if isempty(wsList)
        fprintf('[%s] no workspace files\n', subj);
        continue;
    end
    
    wsNames = lower(string({wsList.name}));
    wsKeep  = ~contains(wsNames, wsRejectTokens);
    wsList  = wsList(wsKeep);
    if isempty(wsList)
        fprintf('[%s] no usable workspace files after filter\n', subj);
        continue;
    end
    
    edfList = [dir(fullfile(subj_dir,'*.edf')); dir(fullfile(subj_dir,'*.EDF'))];
    if isempty(edfList)
        fprintf('[%s] no EDF files\n', subj);
        continue;
    end
    edfNames = lower(string({edfList.name}));
    edfList  = edfList(~contains(edfNames, edfRejectTokens));
    
    fprintf('\n[%s] workspaces=%d | edf=%d\n', subj, numel(wsList), numel(edfList));
    
    % sort workspaces by run number if possible
    runNums = nan(1,numel(wsList));
    for i=1:numel(wsList), runNums(i)=extract_run_number(wsList(i).name); end
    if any(isfinite(runNums))
        [~,o2]=sort(runNums);
        wsList=wsList(o2);
    end
    
    for iws = 1:numel(wsList)
        ws_name = string(wsList(iws).name);
        ws_file = fullfile(wsList(iws).folder, wsList(iws).name);
        
        % ---------- Load needed vars ----------
        [TrialsN, TriggersN, PhasesN, c123, hasTph, badOrderN, nanRateSM, maxPhSample, maxTrigSample] = inspect_workspace(ws_file);
        
        % ---------- Match EDF ----------
        [edf_file, okMatch, runTag, whyMatch] = matchEdfToWorkspace_BASIC(ws_name, edfList, subj_dir);
        
        % ---------- EDF alignment check ----------
        edfNSamp = NaN;
        mismatch = false;
        if okMatch
            [edfNSamp, mismatch] = check_ws_edf_alignment_numbers(edf_file, max([maxPhSample,maxTrigSample]));
        end
        
        % ---------- Print line ----------
        fprintf('  %s | runTag=%s | okEDF=%d | why=%s | Trials=%d Trig=%d Ph=%d | C=[%d %d %d] | badOrder=%d | nanRateSM=%.3f | mismatch=%d\n', ...
            ws_name, runTag, okMatch, whyMatch, TrialsN, TriggersN, PhasesN, c123(1),c123(2),c123(3), badOrderN, nanRateSM, mismatch);
        
        % ---------- Store row ----------
        rows(end+1,:) = { ...
            char(subj), ...
            char(ws_name), ...
            char(runTag), ...
            char(string(edf_file)), ...
            okMatch, char(string(whyMatch)), ...
            TrialsN, TriggersN, PhasesN, ...
            c123(1), c123(2), c123(3), ...
            hasTph, badOrderN, nanRateSM, ...
            maxPhSample, maxTrigSample, edfNSamp, mismatch ...
            }; %#ok<SAGROW>
    end
end

% ---------- Build report table ----------
varNames = { ...
    'Subject','Workspace','RunTag','EDF_File','EDF_Matched','MatchWhy', ...
    'TrialsN','TriggersN','PhasesN', ...
    'Cond1','Cond2','Cond3', ...
    'HasTphases','BadOrderN','NanRateStimMaint', ...
    'MaxPhaseSample','MaxTriggerSample','EDF_TotalSamples','WS_EDF_Mismatch'};

R = cell2table(rows, 'VariableNames', varNames);

% ---------- Summary ----------
fprintf('\n=== SUMMARY ===\n');
fprintf('Total runs checked: %d\n', height(R));
fprintf('EDF matched:        %d (%.1f%%)\n', sum(R.EDF_Matched), 100*mean(R.EDF_Matched));
fprintf('Mismatch flagged:   %d\n', sum(R.WS_EDF_Mismatch==true));
fprintf('Runs w/o T_phases:  %d\n', sum(R.HasTphases==false));
fprintf('Runs badOrder>0:    %d\n', sum(R.BadOrderN>0));
fprintf('Runs nanRate>0:     %d\n', sum(R.NanRateStimMaint>0));

% show worst offenders
bad = R(R.WS_EDF_Mismatch==true | R.HasTphases==false | R.BadOrderN>0 | R.NanRateStimMaint>0, :);
if ~isempty(bad)
    fprintf('\n--- PROBLEM RUNS (first 30) ---\n');
    disp(bad(1:min(30,height(bad)), {'Subject','Workspace','RunTag','EDF_Matched','MatchWhy','WS_EDF_Mismatch','BadOrderN','NanRateStimMaint'}));
end

% ---------- Save ----------
stamp = datestr(now,'yyyymmdd_HHMMSS');
csvFile = fullfile(out_dir, ['ws_edf_report_' stamp '.csv']);
matFile = fullfile(out_dir, ['ws_edf_report_' stamp '.mat']);
try
    writetable(R, csvFile);
    save(matFile,'R');
    fprintf('\n[SAVED]\nCSV: %s\nMAT: %s\n', csvFile, matFile);
catch ME
    warning('Save failed: %s', ME.message);
end

%% 

root_in = '<RESULT_MINE_ROOT>/New Folder/no_ica';
subjects_root = fullfile(root_in,'Subjects');
out_csv = fullfile(root_in,'QA_Tphases_report.csv');

subj_dirs = dir(fullfile(subjects_root,'s*'));
subj_dirs = subj_dirs([subj_dirs.isdir]);

R = table();

for si=1:numel(subj_dirs)
    subj = subj_dirs(si).name;
    sdir = fullfile(subjects_root, subj);
    ws = dir(fullfile(sdir,'workspace_*.mat'));

    for wi=1:numel(ws)
        wfile = fullfile(sdir, ws(wi).name);

        vars = {whos('-file', wfile).name};
        hasPh = ismember('T_phases', vars);
        hasTr = ismember('T_trials_p5', vars);
        hasTg = ismember('T_triggers', vars);

        nPh = NaN; badOrder = NaN; nanRate = NaN;
        c1=0; c2=0; c3=0;

        if hasPh
            S = load(wfile,'T_phases');
            T = S.T_phases;
            if istable(T) && ~isempty(T) ...
               && all(ismember({'StimSample','MaintSample','RetrSample','Condition'}, T.Properties.VariableNames))

                nPh = height(T);
                stim = T.StimSample; maint = T.MaintSample; retr = T.RetrSample;

                badOrder = sum(~(stim < maint & maint < retr));
                nanRate  = mean(isnan(stim) | isnan(maint));  % Stim/Maint

                y = double(T.Condition);
                c1 = sum(y==1); c2 = sum(y==2); c3 = sum(y==3);
            end
        end

        R = [R; table(string(subj), string(ws(wi).name), hasPh, hasTr, hasTg, ...
            nPh, badOrder, nanRate, c1, c2, c3, ...
            'VariableNames', {'Subject','Workspace','hasTphases','hasTrials','hasTriggers',...
            'nPhases','badOrderN','nanRateStimMaint','nC1','nC2','nC3'})]; %#ok<AGROW>
    end
end

writetable(R, out_csv);
disp(R);
fprintf('\nSaved: %s\n', out_csv);

%% 

results_root = '<RESULT_MINE_ROOT>/New Folder/no_ica/GroupResults_LOSO';
phase = 'maint';

ff = dir(fullfile(results_root,'s*',sprintf('*_%s_feat.mat',phase)));
assert(~isempty(ff),'No feat mats found');

rep = table();

for i=1:numel(ff)
    S = load(fullfile(ff(i).folder, ff(i).name),'featAll');
    T = S.featAll;
    if ~istable(T) || height(T)<10, continue; end

    v = T.Properties.VariableNames;
    numIdx = find(varfun(@isnumeric,T,'OutputFormat','uniform'));
    numIdx = numIdx(~ismember(string(v(numIdx)), ["Condition"]));
    if isempty(numIdx), continue; end
    f = v{numIdx(1)};  

    y = double(T.Condition);
    m1 = mean(T.(f)(y==1),'omitnan');
    m2 = mean(T.(f)(y==2),'omitnan');
    m3 = mean(T.(f)(y==3),'omitnan');

    rep = [rep; table(string(ff(i).name), string(f), m1,m2,m3, ...
        'VariableNames', {'File','Feature','meanC1','meanC2','meanC3'})]; %#ok<AGROW>
end

disp(rep(1:min(20,height(rep)),:));

%% ===================== LOCAL FUNCTIONS =====================

function [TrialsN, TriggersN, PhasesN, c123, hasTph, badOrderN, nanRateSM, maxPhSample, maxTrigSample] = inspect_workspace(ws_file)
    TrialsN = NaN; TriggersN = NaN; PhasesN = NaN;
    c123 = [0 0 0];
    hasTph = false;
    badOrderN = NaN;
    nanRateSM = NaN;
    maxPhSample = NaN;
    maxTrigSample = NaN;
    
    v = {whos('-file', ws_file).name};
    
    % trials
    if any(strcmp(v,'T_trials_p5'))
        S = load(ws_file,'T_trials_p5');
        if istable(S.T_trials_p5), TrialsN = height(S.T_trials_p5); end
    end
    
    % triggers
    if any(strcmp(v,'T_triggers'))
        S = load(ws_file,'T_triggers');
        if istable(S.T_triggers)
            TriggersN = height(S.T_triggers);
            maxTrigSample = max_any_sample_col(S.T_triggers, ["sample"]);
        end
    end
    
    % phases
    if any(strcmp(v,'T_phases'))
        S = load(ws_file,'T_phases');
        T = S.T_phases;
        if istable(T) && ~isempty(T)
            hasTph = true;
            PhasesN = height(T);
            
            % class counts from Condition
            if ismember('Condition', T.Properties.VariableNames)
                y = double(T.Condition);
                c123(1) = sum(y==1);
                c123(2) = sum(y==2);
                c123(3) = sum(y==3);
            end
            
            % sanity: order + nanrate
            req = {'StimSample','MaintSample','RetrSample'};
            if all(ismember(req, T.Properties.VariableNames))
                a = double(T.StimSample);
                b = double(T.MaintSample);
                c = double(T.RetrSample);
                
                nanRateSM = mean(isnan(a) | isnan(b));
                
                good = isfinite(a) & isfinite(b) & isfinite(c);
                if any(good)
                    badOrder = ~(a(good) < b(good) & b(good) < c(good));
                    badOrderN = sum(badOrder);
                else
                    badOrderN = NaN;
                end
                
                maxPhSample = max([max(a,[],'omitnan'), max(b,[],'omitnan'), max(c,[],'omitnan')]);
            else
                % if column names differ, still try any "sample" columns
                maxPhSample = max_any_sample_col(T, ["stimsample","maintsample","retrsample","sample"]);
                badOrderN = NaN;
                nanRateSM = NaN;
            end
        end
    end
end

function mx = max_any_sample_col(T, preferredNames)
    if nargin<2, preferredNames = ["sample"]; end
    v = lower(string(T.Properties.VariableNames));
    mx = -inf;
    
    % preferred pass
    for p = 1:numel(preferredNames)
        hit = find(contains(v, lower(preferredNames(p))));
        for k = hit(:)'
            x = T{:,k};
            if isnumeric(x) || islogical(x)
                x = double(x);
                mx = max(mx, max(x,[],'omitnan'));
            end
        end
        if isfinite(mx), return; end
    end
    
    % fallback: any "sample"
    hit = find(contains(v,"sample"));
    for k = hit(:)'
        x = T{:,k};
        if isnumeric(x) || islogical(x)
            x = double(x);
            mx = max(mx, max(x,[],'omitnan'));
        end
    end
    
    if ~isfinite(mx), mx = NaN; end
end

function [edfNSamp, mismatch] = check_ws_edf_alignment_numbers(edf_file, maxWsSample)
    edfNSamp = NaN;
    mismatch = false;
    if ~isfinite(maxWsSample)
        return;
    end
    
    try
        info = edfinfo(edf_file);
        % total samples (use first signal as reference)
        edfNSamp = info.NumSamples(1) * info.NumDataRecords;
        
        if isfinite(edfNSamp) && maxWsSample > edfNSamp
            mismatch = true;
        end
    catch
        edfNSamp = NaN;
        mismatch = false; % unknown
    end
end

function r = extract_run_number(wsName)
    s = lower(string(wsName));
    s = regexprep(s,'\s+',' ');
    s = regexprep(s,'\s*\.mat\s*$','');
    tok = regexp(s, '(\d+)\D*$', 'tokens', 'once');
    if isempty(tok), r = NaN;
    else, r = str2double(tok{1});
    end
end

function [edf_file, ok, runTag, why] = matchEdfToWorkspace_BASIC(wsName, edfList, subj_dir)
    ok = false; why = "noMatch";
    edf_file = "";
    runTag = "run?";
    
    [wsCore, wsBase, wsRun, okWs] = parse_core_run(wsName, true);
    if ~okWs
        why = "badWSname";
        return;
    end
    runTag = regexprep(wsBase + "_" + wsRun, '[^a-zA-Z0-9_]+', '_');
    
    nE = numel(edfList);
    edCore = strings(nE,1); edBase = strings(nE,1); edRun = strings(nE,1); okEd = false(nE,1);
    for i=1:nE
        [c,b,r,ok1] = parse_core_run(edfList(i).name, false);
        edCore(i)=c; edBase(i)=b; edRun(i)=r; okEd(i)=ok1;
    end
    
    % PASS1: exact core
    hit = find(okEd & (edCore==wsCore));
    if numel(hit)==1
        edf_file = fullfile(subj_dir, edfList(hit).name);
        ok = true; why = "ok_core";
        return;
    elseif numel(hit)>1
        % pick shortest name
        lens = arrayfun(@(k) strlength(string(edfList(k).name)), hit);
        [~,ii]=min(lens);
        pick = hit(ii);
        edf_file = fullfile(subj_dir, edfList(pick).name);
        ok = true; why = "multi_core_pickShortest";
        return;
    end
    
    % PASS2: same run + base exact
    hit = find(okEd & (edRun==wsRun) & (edBase==wsBase));
    if numel(hit)==1
        edf_file = fullfile(subj_dir, edfList(hit).name);
        ok = true; why = "ok_baseRun";
        return;
    elseif numel(hit)>1
        why = "ambiguous_baseRun";
        return;
    end
    
    % PASS3: fallback within same run (score contains)
    cand = find(okEd & (edRun==wsRun));
    if isempty(cand)
        why = "noRunMatch";
        return;
    end
    
    sc = zeros(numel(cand),1);
    for k=1:numel(cand)
        sc(k) = score_name_match(wsBase, edBase(cand(k)), edfList(cand(k)).name);
    end
    [best,ib]=max(sc);
    if best < 150
        why = "lowConfidence";
        return;
    end
    edf_file = fullfile(subj_dir, edfList(cand(ib)).name);
    ok = true; why = "ok_scored";
end

function [coreNorm, baseNorm, runStr, ok] = parse_core_run(fname, isMat)
    ok = false;
    coreNorm = ""; baseNorm = ""; runStr = "";
    
    s = lower(string(fname));
    s = strtrim(s);
    s = regexprep(s, '\s+', ' ');
    s = regexprep(s, '\s*\.(mat|edf)\s*$', '');
    
    if isMat
        s = regexprep(s, '^workspace[_\-\s]*', '');
    end
    
    % clean junk
    s = regexprep(s, '\(copy\)', '');
    s = regexprep(s, 'copy', '');
    s = regexprep(s, 'filtered', '');
    s = regexprep(s, 'backup', '');
    s = regexprep(s, 'preproc', '');
    s = strtrim(s);
    
    tok = regexp(s, '(\d+)\D*$', 'tokens', 'once');
    if isempty(tok), return; end
    runStr = string(tok{1});
    
    baseRaw = regexprep(s, '(\d+)\D*$', '');
    baseRaw = strtrim(baseRaw);
    
    baseNorm = regexprep(baseRaw, '[^a-z0-9]+', '');
    if baseNorm=="" || runStr=="", return; end
    
    coreNorm = baseNorm + runStr; % e.g., masume1
    ok = true;
end

function sc = score_name_match(wsBase, edBase, edfName)
    sc = 0;
    if edBase==wsBase, sc=sc+500; end
    if contains(edBase, wsBase), sc=sc+250+2*strlength(wsBase); end
    if contains(wsBase, edBase), sc=sc+200+2*strlength(edBase); end
    
    n = lower(string(edfName));
    if contains(n,"eye") || contains(n,"filtered") || contains(n,"ica") || contains(n,"clean") || ...
       contains(n,"preproc") || contains(n,"backup") || contains(n,"copy")
        sc = sc - 200;
    end
    sc = sc - 0.2*strlength(n);
end

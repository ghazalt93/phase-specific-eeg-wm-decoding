%% CHECK_SUBJECTS_CHAN_SELECTED

cfg = get_project_config();

clear; clc;

%% ================= CONFIG =================
cfg.root = cfg.outputRoot;

cfg.targetSubjects = {'s1','s3','s5','s6','s7','s8','s9', ...
                      's11','s12','s13','s14','s15','s16','s17','s18','s19','s20','s21','s22','s23','s24','s25'};

cfg.outDir = fullfile(cfg.root, 'SUBJECT_AUDIT_CHAN_SELECTED');
if ~exist(cfg.outDir, 'dir'); mkdir(cfg.outDir); end

reportFile = fullfile(cfg.outDir, 'SUBJECT_AUDIT_CHAN_SELECTED_Report.txt');
if exist(reportFile, 'file'); delete(reportFile); end
diary(reportFile);

fprintf('=== SUBJECT AUDIT: CHAN SELECTED PATH ===\n');
fprintf('Generated: %s\n', datestr(now));
fprintf('Root path: %s\n\n', cfg.root);

target = local_sort_subjects(string(cfg.targetSubjects(:)));
fprintf('Target list currently set to %d subjects:\n%s\n\n', numel(target), strjoin(target, ', '));

%% ================= SCAN ALL CANDIDATE FILES =================
Rows = {};
AllSubjectsRows = {};

% 1) scan folders
folderSubjects = local_subjects_from_folders(cfg.root);
[missing, extra, status] = local_compare_subjects(target, folderSubjects);
Rows(end+1,:) = {'Folders named s* under root', cfg.root, numel(folderSubjects), ...
    char(strjoin(folderSubjects, ', ')), char(strjoin(missing, ', ')), char(strjoin(extra, ', ')), status}; %#ok<AGROW>
local_print_block('Folders named s* under root', cfg.root, folderSubjects, missing, extra, status);

% 2) scan clean epoch mat files
matSubjects = local_subjects_from_mat_files(cfg.root);
[missing, extra, status] = local_compare_subjects(target, matSubjects);
Rows(end+1,:) = {'MAT files with subject/run names', cfg.root, numel(matSubjects), ...
    char(strjoin(matSubjects, ', ')), char(strjoin(missing, ', ')), char(strjoin(extra, ', ')), status}; %#ok<AGROW>
local_print_block('MAT files with subject/run names', cfg.root, matSubjects, missing, extra, status);

% 3) scan important CSV tables by likely names
patterns = { ...
    '**/STEP208_SubjectMedianConnectivityTable.csv', ...
    '**/*SubjectMedian*Connectivity*.csv', ...
    '**/*LOSO*Predictions*.csv', ...
    '**/*SelectedFeatures*.csv', ...
    '**/*FeatureStability*.csv', ...
    '**/*PermutationSummary*.csv', ...
    '**/*BestModels*.csv', ...
    '**/*Subject*.csv', ...
    '**/*subject*.csv'};

seenFiles = strings(0,1);
for ip = 1:numel(patterns)
    files = dir(fullfile(cfg.root, patterns{ip}));
    for i = 1:numel(files)
        fpath = string(fullfile(files(i).folder, files(i).name));
        if any(seenFiles == fpath)
            continue;
        end
        seenFiles(end+1,1) = fpath; %#ok<AGROW>

        subj = local_subjects_from_csv(char(fpath));
        [missing, extra, status] = local_compare_subjects(target, subj);

        Rows(end+1,:) = {char(files(i).name), char(fpath), numel(subj), ...
            char(strjoin(subj, ', ')), char(strjoin(missing, ', ')), char(strjoin(extra, ', ')), status}; %#ok<AGROW>

        if ~isempty(subj)
            for s = 1:numel(subj)
                AllSubjectsRows(end+1,:) = {char(subj(s)), char(files(i).name), char(fpath)}; %#ok<AGROW>
            end
        end

        fprintf('\n--- CSV: %s ---\n', files(i).name);
        fprintf('Path: %s\n', fpath);
        fprintf('N subjects: %d\n', numel(subj));
        if isempty(subj)
            fprintf('Subjects: none detected\n');
        else
            fprintf('Subjects: %s\n', strjoin(subj, ', '));
        end
        if ~isempty(extra)
            fprintf('Extra not in target: %s\n', strjoin(extra, ', '));
        end
        if ~isempty(missing)
            fprintf('Missing from target: %s\n', strjoin(missing, ', '));
        end
        fprintf('Status: %s\n', status);
    end
end

%% ================= SAVE SUMMARY =================
Audit = cell2table(Rows, 'VariableNames', ...
    {'Source','Path','NSubjects','SubjectList','MissingFromTarget','ExtraNotInTarget','Status'});
writetable(Audit, fullfile(cfg.outDir, 'SUBJECT_AUDIT_CHAN_SELECTED_Table.csv'));

if ~isempty(AllSubjectsRows)
    Occ = cell2table(AllSubjectsRows, 'VariableNames', {'Subject','FileName','Path'});
    writetable(Occ, fullfile(cfg.outDir, 'SUBJECT_AUDIT_CHAN_SELECTED_AllOccurrences.csv'));
else
    Occ = table();
end

%% ================= PRINT FINAL CHECK =================
fprintf('\n\n=== FINAL SUMMARY ===\n');
disp(Audit(:, {'Source','NSubjects','MissingFromTarget','ExtraNotInTarget','Status'}));

fprintf('\nFiles with 23 subjects:\n');
m23 = Audit.NSubjects == 23;
if any(m23)
    disp(Audit(m23, {'Source','Path','SubjectList','ExtraNotInTarget'}));
else
    fprintf('No source with exactly 23 detected in scanned files.\n');
end

fprintf('\nFiles with extra subjects outside target list:\n');
mExtra = string(Audit.ExtraNotInTarget) ~= "";
if any(mExtra)
    disp(Audit(mExtra, {'Source','Path','NSubjects','ExtraNotInTarget'}));
else
    fprintf('No extra subjects outside target list detected.\n');
end

fprintf('\nSaved report:\n%s\n', reportFile);
fprintf('Saved table:\n%s\n', fullfile(cfg.outDir, 'SUBJECT_AUDIT_CHAN_SELECTED_Table.csv'));
fprintf('Saved occurrences:\n%s\n', fullfile(cfg.outDir, 'SUBJECT_AUDIT_CHAN_SELECTED_AllOccurrences.csv'));
fprintf('\n=== AUDIT COMPLETE ===\n');

diary off;

%% ================= LOCAL FUNCTIONS =================
function subj = local_subjects_from_folders(rootPath)
    subj = strings(0,1);
    if exist(rootPath, 'dir') ~= 7
        return;
    end
    d = dir(fullfile(rootPath, '**', 's*'));
    for i = 1:numel(d)
        if d(i).isdir
            name = string(d(i).name);
            if ~isempty(regexp(char(name), '^s\d+$', 'once'))
                subj(end+1,1) = name; %#ok<AGROW>
            end
        end
    end
    subj = unique(local_sort_subjects(local_normalize_subject(subj)), 'stable');
end

function subj = local_subjects_from_mat_files(rootPath)
    subj = strings(0,1);
    if exist(rootPath, 'dir') ~= 7
        return;
    end
    files = dir(fullfile(rootPath, '**', '*.mat'));
    for i = 1:numel(files)
        nm = string(files(i).name);
        tok = regexp(char(nm), '(s\d+)', 'tokens', 'once');
        if ~isempty(tok)
            subj(end+1,1) = string(tok{1}); %#ok<AGROW>
        end
    end
    subj = unique(local_sort_subjects(local_normalize_subject(subj)), 'stable');
end

function subj = local_subjects_from_csv(file)
    subj = strings(0,1);
    if exist(file, 'file') ~= 2
        return;
    end

    try
        T = readtable(file, 'PreserveVariableNames', true);
    catch
        try
            T = readtable(file);
        catch
            return;
        end
    end

    cand = {'Subject','TestSubject','subj','participant','Participant','SubjectID','SubID','ID'};
    col = local_find_col(T.Properties.VariableNames, cand);

    if ~isempty(col)
        subj = local_normalize_subject(local_to_string(T.(col)));
        subj = subj(subj ~= "");
        subj = unique(local_sort_subjects(subj), 'stable');
        return;
    end

    % If no subject column, try scanning first few text-like columns
    maxCols = min(width(T), 12);
    tmp = strings(0,1);
    for j = 1:maxCols
        try
            v = local_to_string(T.(T.Properties.VariableNames{j}));
            for k = 1:numel(v)
                tok = regexp(char(v(k)), '(s\d+)', 'tokens');
                for t = 1:numel(tok)
                    tmp(end+1,1) = string(tok{t}{1}); %#ok<AGROW>
                end
            end
        catch
        end
    end
    subj = unique(local_sort_subjects(local_normalize_subject(tmp)), 'stable');
end

function [missing, extra, status] = local_compare_subjects(target, subj)
    target = local_sort_subjects(string(target(:)));
    subj = local_sort_subjects(string(subj(:)));

    if isempty(subj)
        missing = target;
        extra = strings(0,1);
        status = 'NO_SUBJECTS_DETECTED';
        return;
    end

    missing = setdiff(target, subj, 'stable');
    extra = setdiff(subj, target, 'stable');

    if isempty(missing) && isempty(extra) && numel(subj) == numel(target)
        status = 'OK';
    elseif ~isempty(extra) && isempty(missing)
        status = 'HAS_EXTRA_SUBJECTS';
    elseif isempty(extra) && ~isempty(missing)
        status = 'MISSING_TARGET_SUBJECTS';
    else
        status = 'MISMATCH';
    end
end

function local_print_block(label, pathValue, subj, missing, extra, status)
    fprintf('\n--- %s ---\n', label);
    fprintf('Path: %s\n', pathValue);
    fprintf('N subjects: %d\n', numel(subj));
    if isempty(subj)
        fprintf('Subjects: none detected\n');
    else
        fprintf('Subjects: %s\n', strjoin(subj, ', '));
    end
    if isempty(missing)
        fprintf('Missing from target: none\n');
    else
        fprintf('Missing from target: %s\n', strjoin(missing, ', '));
    end
    if isempty(extra)
        fprintf('Extra not in target: none\n');
    else
        fprintf('Extra not in target: %s\n', strjoin(extra, ', '));
    end
    fprintf('Status: %s\n', status);
end

function col = local_find_col(vnames, keys)
    col = '';
    vnLower = lower(regexprep(vnames, '[^a-zA-Z0-9]', ''));
    for k = 1:numel(keys)
        key = lower(regexprep(keys{k}, '[^a-zA-Z0-9]', ''));
        idx = find(strcmp(vnLower, key), 1);
        if ~isempty(idx)
            col = vnames{idx};
            return;
        end
    end
    for k = 1:numel(keys)
        key = lower(regexprep(keys{k}, '[^a-zA-Z0-9]', ''));
        idx = find(contains(vnLower, key), 1);
        if ~isempty(idx)
            col = vnames{idx};
            return;
        end
    end
end

function s = local_to_string(x)
    if iscell(x)
        s = string(x);
    elseif iscategorical(x)
        s = string(x);
    elseif isnumeric(x) || islogical(x)
        s = string(x);
    elseif isstring(x)
        s = x;
    elseif ischar(x)
        s = string(cellstr(x));
    else
        try
            s = string(x);
        catch
            s = repmat("", size(x,1), 1);
        end
    end
    s = s(:);
    s(ismissing(s)) = "";
    s = strtrim(s);
end

function subj = local_normalize_subject(s)
    subj = lower(strtrim(string(s)));
    subj = regexprep(subj, '\s+', '');
    subj = regexprep(subj, '^subj', 's');
    subj = regexprep(subj, '^subject', 's');
    subj = regexprep(subj, '^participant', 's');

    for i = 1:numel(subj)
        nums = regexp(char(subj(i)), '\d+', 'match');
        if ~isempty(nums)
            subj(i) = "s" + string(str2double(nums{1}));
        end
    end
end

function s2 = local_sort_subjects(s)
    s = string(s(:));
    s = s(s ~= "");
    if isempty(s)
        s2 = strings(0,1);
        return;
    end

    nums = nan(numel(s),1);
    for i = 1:numel(s)
        tok = regexp(char(s(i)), '\d+', 'match');
        if ~isempty(tok)
            nums(i) = str2double(tok{1});
        else
            nums(i) = inf;
        end
    end
    T = table(nums, s);
    T = sortrows(T, {'nums','s'});
    s2 = T.s;
end

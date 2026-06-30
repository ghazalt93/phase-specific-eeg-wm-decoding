function summaryTbl = check_edf_epochs_trials(rootDir)

    if nargin < 1 || isempty(rootDir)
        rootDir = uigetdir(pwd, 'Select root folder of Subjects (s1, s2, ...)');
        if rootDir == 0
            error('No folder selected.');
        end
    end

    baseRow = struct( ...
        'Subject',        '', ...
        'EDFfile',        '', ...
        'FullPath',       '', ...
        'Type',           '', ...  % TASK / REST
        'Status',         '', ...  
        'nbchan',         NaN, ...
        'nsamples',       NaN, ...
        'fs',             NaN, ...
        'nTriggers',      NaN, ...
        'nTrialsTotal',   NaN, ...
        'nCond1_40',      NaN, ...
        'nCond2_60',      NaN, ...
        'nCond3_16',      NaN, ...
        'nEpochsCond1',   NaN, ...
        'nEpochsCond2',   NaN, ...
        'nEpochsCond3',   NaN, ...
        'Message',        ''  ...
    );

    results = baseRow([]);  % 0x1 struct

    preTime  = 0.5;   
    postTime = 1.5;   
    
    subjDirs = dir(fullfile(rootDir, 's*'));
    subjDirs = subjDirs([subjDirs.isdir]);

    fprintf('Found %d subject folders in %s\n', numel(subjDirs), rootDir);

    for isub = 1:numel(subjDirs)
        subjName = subjDirs(isub).name;
        subjPath = fullfile(rootDir, subjName);

        edfFiles = dir(fullfile(subjPath, '*.edf'));
        fprintf('--------------------------------------------------\n');
        fprintf('Subject %d/%d : %s  (%d EDF files)\n', ...
            isub, numel(subjDirs), subjName, numel(edfFiles));

        for ifile = 1:numel(edfFiles)
            edfName = edfFiles(ifile).name;
            edfPath = fullfile(subjPath, edfName);

            fprintf('  [%d/%d] %s\n', ifile, numel(edfFiles), edfName);

            R = baseRow;
            R.Subject  = subjName;
            R.EDFfile  = edfName;
            R.FullPath = edfPath;

            lname = lower(edfName);
            if contains(lname, 'eye') && (contains(lname, 'open') || contains(lname, 'close')) ...
                    || contains(lname, 'rest')
                R.Type   = 'REST';
                R.Status = 'SKIP_REST';
                R.Message = 'REST file (eye open/close or rest) - skipped.';
                results(end+1) = R; %#ok<AGROW>
                fprintf('    -> REST file, skipped.\n');
                continue;
            else
                R.Type = 'TASK';
            end

            try
                try
                   
                    EEG = pop_biosig(edfPath);
                    data = double(EEG.data);  % [nbchan x nsamples]
                    fs   = EEG.srate;
                catch ME1
                    hdr = sopen(edfPath, 'r', 'OVERFLOWDETECTION:OFF');
                    [S, hdr] = sread(hdr);
                    sclose(hdr);
                    data = double(S'); 
                    if isfield(hdr, 'SampleRate')
                        fs = hdr.SampleRate(1);
                    elseif isfield(hdr, 'EVENT') && isfield(hdr.EVENT, 'SampleRate')
                        fs = hdr.EVENT.SampleRate;
                    else
                        fs = NaN;
                    end
                end

                [nbchan, nsamples] = size(data);
                R.nbchan   = nbchan;
                R.nsamples = nsamples;
                R.fs       = fs;

                if nbchan < 65
                    R.Status  = 'NO_TRIGGER_CHANNEL';
                    R.Message = sprintf('nbchan=%d < 65, cannot use last channel as trigger.', nbchan);
                    results(end+1) = R; %#ok<AGROW>
                    fprintf('    -> %s\n', R.Message);
                    continue;
                end

                trig = data(nbchan, :);

                [nTriggers, trialInfo, epochInfo, msg] = local_analyse_triggers( ...
                    trig, fs, preTime, postTime);

                R.nTriggers      = nTriggers;
                R.nTrialsTotal   = trialInfo.nTrialsTotal;
                R.nCond1_40      = trialInfo.nCond1;
                R.nCond2_60      = trialInfo.nCond2;
                R.nCond3_16      = trialInfo.nCond3;
                R.nEpochsCond1   = epochInfo.nEpochsCond1;
                R.nEpochsCond2   = epochInfo.nEpochsCond2;
                R.nEpochsCond3   = epochInfo.nEpochsCond3;
                R.Status         = trialInfo.status;
                R.Message        = msg;

                fprintf('    -> %s\n', msg);

            catch ME
                R.Status  = 'ERROR';
                R.Message = sprintf('%s: %s', ME.identifier, ME.message);
                fprintf('    !! ERROR: %s\n', R.Message);
            end

            results(end+1) = R; %#ok<AGROW>
        end
    end

    if isempty(results)
        warning('No EDF files processed.');
        summaryTbl = table();
    else
        summaryTbl = struct2table(results);
    end

    outCsv = fullfile(rootDir, 'EDF_check_summary.csv');
    try
        writetable(summaryTbl, outCsv);
        fprintf('\nSummary table saved to:\n  %s\n', outCsv);
    catch
        fprintf('\nCould not save CSV automatically. You can save summaryTbl manually.\n');
    end
end

function [nTriggers, trialInfo, epochInfo, msg] = local_analyse_triggers(trig, fs, preTime, postTime)


    trig = trig(:);                     
    trigBin = (trig ~= 0);             
    onsets = find(diff([0; trigBin]) == 1);  
    codes  = round(trig(onsets));     

    nTriggers = numel(onsets);

    trialInfo = struct('nTrialsTotal', 0, ...
                       'nCond1',      0, ...
                       'nCond2',      0, ...
                       'nCond3',      0, ...
                       'status',      'NO_TRIALS');
    epochInfo = struct('nEpochsCond1', 0, ...
                       'nEpochsCond2', 0, ...
                       'nEpochsCond3', 0);

    if nTriggers == 0
        msg = 'No triggers detected on last channel.';
        return;
    end
    if isempty(fs) || isnan(fs)
        msg = sprintf('Triggers=%d, but fs is unknown. Trials counted, epochs skipped.', nTriggers);
    end

    is10        = (codes == 10);
    isMid       = ismember(codes, [20 30 40]);
    isRetrieval = ismember(codes, [40 60 16]);

    respSamples = [];
    respConds   = [];

    for k = 1:(numel(codes)-3)
        if is10(k) && isMid(k+1) && is10(k+2) && isRetrieval(k+3)
            respSamples(end+1,1) = onsets(k+3);   %#ok<AGROW>
            respConds(end+1,1)   = codes(k+3);   %#ok<AGROW>
        end
    end

    nTrials = numel(respSamples);
    if nTrials == 0
        msg = sprintf('Triggers=%d, but no valid 10-20/30/40-10-40/60/16 patterns found.', nTriggers);
        return;
    end

    trialInfo.nTrialsTotal = nTrials;
    trialInfo.nCond1       = sum(respConds == 40);
    trialInfo.nCond2       = sum(respConds == 60);
    trialInfo.nCond3       = sum(respConds == 16);
    trialInfo.status       = 'OK';

    if isempty(fs) || isnan(fs)
        epochInfo.nEpochsCond1 = NaN;
        epochInfo.nEpochsCond2 = NaN;
        epochInfo.nEpochsCond3 = NaN;
        msg = sprintf('Triggers=%d, Trials=%d (C1=%d,C2=%d,C3=%d). fs unknown → epochs not checked.', ...
            nTriggers, nTrials, trialInfo.nCond1, trialInfo.nCond2, trialInfo.nCond3);
        return;
    end

    preS  = round(preTime  * fs);
    postS = round(postTime * fs);
    N     = numel(trig);

    validMask = (respSamples - preS >= 1) & (respSamples + postS <= N);
    validConds = respConds(validMask);

    epochInfo.nEpochsCond1 = sum(validConds == 40);
    epochInfo.nEpochsCond2 = sum(validConds == 60);
    epochInfo.nEpochsCond3 = sum(validConds == 16);

    msg = sprintf(['Triggers=%d, Trials=%d (C1=%d,C2=%d,C3=%d), ' ...
                   'Epochs=(%d,%d,%d) in [%.1f %.1f] s window.'], ...
        nTriggers, nTrials, ...
        trialInfo.nCond1, trialInfo.nCond2, trialInfo.nCond3, ...
        epochInfo.nEpochsCond1, epochInfo.nEpochsCond2, epochInfo.nEpochsCond3, ...
        -preTime, postTime);
end

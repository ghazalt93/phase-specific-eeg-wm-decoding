function wmFeat = compute_wm_features_for_phase(ep, t, fs, trialInfo, useChanIdx, cfg, outDirWM, baseName, verbose)


if nargin < 5, useChanIdx = []; end
if nargin < 6 || isempty(cfg), cfg = struct(); end
if nargin < 7 || isempty(outDirWM), outDirWM = pwd; end
if nargin < 8 || isempty(baseName), baseName = 'wm'; end
if nargin < 9 || isempty(verbose), verbose = true; end

% ---- normalize types ----
if isstring(outDirWM), outDirWM = char(outDirWM); end
if isstring(baseName), baseName = char(baseName); end
outDirWM = strtrim(outDirWM);
if exist(outDirWM,'dir') ~= 7
    mkdir(outDirWM);
end

% ---- validate ep/t ----
if ~isnumeric(ep) || ndims(ep) ~= 3
    error('compute_wm_features_for_phase:BadEp','ep must be numeric [chan x time x trial].');
end
if ~isnumeric(t) || ~isvector(t)
    error('compute_wm_features_for_phase:BadT','t must be numeric vector.');
end
t = double(t(:)'); % row

% auto-fix dimension order if needed
[n1,n2,n3] = size(ep);
if numel(t) == n1 && numel(t) ~= n2
    % likely [time x chan x trial] -> [chan x time x trial]
    ep = permute(ep, [2 1 3]);
    [n1,n2,n3] = size(ep);
end
if numel(t) ~= n2
    error('compute_wm_features_for_phase:TimeMismatch','Length(t)=%d but size(ep,2)=%d', numel(t), n2);
end

nTrials = n3;
if istable(trialInfo)
    if height(trialInfo) ~= nTrials
        % align by truncation
        m = min(height(trialInfo), nTrials);
        trialInfo = trialInfo(1:m,:);
        ep = ep(:,:,1:m);
        nTrials = m;
        if verbose
            fprintf('[WM][WARN] trialInfo/trials mismatch -> truncated to %d trials\n', m);
        end
    end
else
    trialInfo = table();
end

if isempty(useChanIdx)
    useChanIdx = 1:size(ep,1);
end

% ensure double
ep = double(ep);
fs = double(fs);

% ---- drop channels that are entirely non-finite (prevents all-NaN features) ----
finiteByChan = false(size(ep,1),1);
for ch = 1:size(ep,1)
    finiteByChan(ch) = any(isfinite(ep(ch,:,:)), 'all');
end
useChanIdx = useChanIdx(:)';
useChanIdx = useChanIdx(useChanIdx>=1 & useChanIdx<=size(ep,1));
useChanIdx = useChanIdx(finiteByChan(useChanIdx));

if isempty(useChanIdx)
    % Save a small debug file to outDirWM to show the issue
    try
        dbg.epFiniteFrac = squeeze(mean(isfinite(ep), [2 3]));
        dbg.fs = fs;
        dbg.t = t;
        save(fullfile(outDirWM, [baseName '_DEBUG_allNaN.mat']), 'dbg', '-v7');
    catch
    end
    error('compute_wm_features_for_phase:AllNaN', ...
        'All selected channels are non-finite (NaN/Inf) in ep. This usually means EDF channel labels did not match expected labels, so data was padded with NaN.');
end


wmFeat = table();

% function list (override-able)
funcs = { ...
    'WM_Compute_FrequencyFeatures', ...
    'WM_BandpowerFeatures_FFT_FAST', ...
    'WM_Compute_TemporalFeatures', ...
    'WM_Compute_StatisticalFeatures', ...
    'WM_Compute_ComplexityFeatures', ...
    'WM_Compute_TFeatures_Morlet', ...  % may be mismatched inside -> handled
    'WM_Compute_BehavioralFeatures' ...
    };

if isfield(cfg,'wm_funcs') && ~isempty(cfg.wm_funcs)
    funcs = cfg.wm_funcs;
end

% ---- PSD plot/save (summary) ----
try
    doPSD = true;
    if isfield(cfg,'psd') && isfield(cfg.psd,'enable'), doPSD = logical(cfg.psd.enable); end
    if doPSD
        save_psd_summary(ep, t, fs, useChanIdx, outDirWM, [baseName '_PSD'], verbose);
    end
catch ME
    if verbose
        fprintf('[PSD][FAIL] %s\n', ME.message);
    end
end

% ---- run WM functions ----
for i = 1:numel(funcs)
    fn = funcs{i};
    if ~exist(fn,'file')
        if verbose, fprintf('[WM][MISS] %s not on path\n', fn); end
        continue;
    end

    try
        out = [];
        switch fn
            case 'WM_BandpowerFeatures_FFT_FAST'
                out = feval(fn, ep, fs, t, useChanIdx, cfg);

            case 'WM_Compute_BehavioralFeatures'
                out = feval(fn, trialInfo, cfg);

            case 'WM_Compute_TFeatures_Morlet'
                % file you uploaded has a name mismatch: function is WM_Compute_TFFeatures_Morlet
                cfgM = cfg;
                cfgM.fs = fs;
                try
                    out = feval('WM_Compute_TFeatures_Morlet', ep, t, useChanIdx, cfgM);
                catch
                    if exist('WM_Compute_TFFeatures_Morlet','file') == 2
                        out = feval('WM_Compute_TFFeatures_Morlet', ep, t, useChanIdx, cfgM);
                    else
                        rethrow(lasterror); %#ok<LERR>
                    end
                end

            otherwise
                % standard signature: (ep, t, useChanIdx, cfg)
                out = feval(fn, ep, t, useChanIdx, cfg);
        end

        % Save raw output for debugging
        try
            outPath = fullfile(outDirWM, [baseName '_' fn '.mat']);
            save(outPath, 'out', '-v7');
        catch
        end

        T = to_table(out, nTrials);
        if isempty(T) || ~istable(T) || height(T) ~= nTrials || width(T)==0
            if verbose
                fprintf('[WM][SKIP] %s returned empty or wrong size\n', fn);
            end
            continue;
        end

        % Prefix columns to avoid collisions and to be easily countable
        T = prefix_vars(T, ['WM_' fn '_']);

        wmFeat = merge_tables_horiz(wmFeat, T);

        if verbose
            fprintf('[WM][OK] %s -> +%d cols\n', fn, width(T));
        end

    catch ME
        if verbose
            fprintf('[WM][FAIL] %s: %s\n', fn, ME.message);
        end
    end
end

end

% ================= helpers =================

function T = to_table(out, nTrials)
T = table();
if isempty(out), return; end

if istable(out)
    T = out;
    return;
end

if isnumeric(out)
    if size(out,1) == nTrials
        T = array2table(out);
        T.Properties.VariableNames = make_names('X', size(out,2));
        return;
    elseif size(out,2) == nTrials
        T = array2table(out');
        T.Properties.VariableNames = make_names('X', size(out,1));
        return;
    end
end

if isstruct(out)
    % scalar struct with numeric vectors/matrices
    f = fieldnames(out);
    cols = {};
    data = [];
    for i=1:numel(f)
        v = out.(f{i});
        if isnumeric(v)
            if isvector(v) && numel(v)==nTrials
                cols{end+1} = f{i}; %#ok<AGROW>
                data(:,end+1) = v(:); %#ok<AGROW>
            elseif ismatrix(v) && size(v,1)==nTrials
                for k=1:size(v,2)
                    cols{end+1} = sprintf('%s_%d', f{i}, k); %#ok<AGROW>
                end
                data = [data, v]; %#ok<AGROW>
            end
        end
    end
    if ~isempty(cols)
        T = array2table(data);
        T.Properties.VariableNames = make_valid(cols);
        return;
    end
end

try
    T0 = struct2table(out);
    if height(T0)==nTrials
        T = T0;
    end
catch
end
end

function T = prefix_vars(T, prefix)
vars = T.Properties.VariableNames;
for i=1:numel(vars)
    vars{i} = [prefix vars{i}];
end
vars = make_valid(vars);
vars = make_unique(vars);
T.Properties.VariableNames = vars;
end

function vars = make_names(base, n)
vars = cell(1,n);
for i=1:n
    vars{i} = sprintf('%s_%03d', base, i);
end
end

function out = make_valid(in)
out = in;
for i=1:numel(out)
    s = in{i};
    s = regexprep(s,'[^a-zA-Z0-9_]','_');
    if isempty(s) || ~isletter(s(1))
        s = ['x_' s];
    end
    out{i} = s;
end
end

function out = make_unique(in)
out = in;
seen = containers.Map('KeyType','char','ValueType','double');
for i=1:numel(out)
    k = out{i};
    if ~isKey(seen,k)
        seen(k) = 1;
    else
        seen(k) = seen(k) + 1;
        out{i} = sprintf('%s_%d', k, seen(k));
    end
end
end

function C = merge_tables_horiz(A, B)
if isempty(A)
    C = B;
    return;
end
if height(A) ~= height(B)
    m = min(height(A), height(B));
    A = A(1:m,:);
    B = B(1:m,:);
end
bvars = B.Properties.VariableNames;
avars = A.Properties.VariableNames;
collide = ismember(bvars, avars);
if any(collide)
    bvars(collide) = make_unique(bvars(collide));
    B.Properties.VariableNames = bvars;
end
C = [A B];
end

function save_psd_summary(ep, t, fs, useChanIdx, outDir, base, verbose)
if isempty(useChanIdx)
    useChanIdx = 1:size(ep,1);
end
useChanIdx = useChanIdx(:)';
useChanIdx = useChanIdx(useChanIdx>=1 & useChanIdx<=size(ep,1));
if isempty(useChanIdx), return; end

nChan = numel(useChanIdx);
nTrial = size(ep,3);
maxSeg = 100;

Pacc = [];
F = [];

segCount = 0;
for tr=1:nTrial
    for ci=1:nChan
        x = squeeze(ep(useChanIdx(ci), :, tr));
        x = x(:);
        x = x(~isnan(x));
        if numel(x) < 16, continue; end
        [Pxx,Ftmp] = pwelch_fallback(x, fs);
        if isempty(Pxx), continue; end
        if isempty(Pacc)
            Pacc = zeros(numel(Pxx),1);
            F = Ftmp;
        end
        if numel(Pxx)==numel(Pacc)
            Pacc = Pacc + Pxx(:);
            segCount = segCount + 1;
        end
        if segCount >= maxSeg, break; end
    end
    if segCount >= maxSeg, break; end
end
if isempty(Pacc) || segCount==0, return; end
Pmean = Pacc / segCount;

try
    save(fullfile(outDir, [base '.mat']), 'F', 'Pmean', 'fs', 't');
end

try
    h = figure('Visible','off');
    plot(F, 10*log10(Pmean));
    grid on;
    xlabel('Hz'); ylabel('Power (dB/Hz)');
    title(strrep(base,'_','\_'));
    xlim([0 60]);
    pngPath = fullfile(outDir, [base '.png']);
    figPath = fullfile(outDir, [base '.fig']);
    saveas(h, pngPath);
    saveas(h, figPath);
    close(h);
    if verbose, fprintf('[PSD][SAVE] %s\n', pngPath); end
catch ME
    if verbose, fprintf('[PSD][FAIL] plot/save: %s\n', ME.message); end
end
end

function [Pxx,F] = pwelch_fallback(x, fs)
Pxx = []; F = [];
try
    if exist('pwelch','file')==2
        n = numel(x);
        wlen = max(64, 2^nextpow2(min(n, 512)));
        nover = round(0.5*wlen);
        nfft = max(wlen, 256);
        [Pxx,F] = pwelch(x, wlen, nover, nfft, fs);
    else
        nfft = max(256, 2^nextpow2(numel(x)));
        [Pxx,F] = periodogram(x, [], nfft, fs);
    end
catch
    Pxx = []; F = [];
end
end

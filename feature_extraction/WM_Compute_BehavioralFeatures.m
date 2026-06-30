function behFeat = WM_Compute_BehavioralFeatures(trialInfo, cfg)

if nargin < 2 || isempty(cfg), cfg = struct(); end

% ---- pick phase if trialInfo is struct ----
trialInfo = pickTrialInfoIfStruct(trialInfo, cfg);

if isempty(trialInfo) || ~istable(trialInfo)
    behFeat = table();
    return;
end

n = height(trialInfo);
behFeat = table();
vnames = lower(string(trialInfo.Properties.VariableNames));

% defaults
if ~isfield(cfg,'rt_min') || isempty(cfg.rt_min), cfg.rt_min = -Inf; end
if ~isfield(cfg,'rt_max') || isempty(cfg.rt_max), cfg.rt_max = +Inf; end

% --- helper: find likely columns ---
findCol = @(patterns) find(any(contains(vnames, patterns), 2), 1, 'first');

ixRT   = findCol(["rt","reaction","responsetime","response_time","latency"]);
ixAcc  = findCol(["acc","correct","accuracy","hit","iscorrect"]);
ixResp = findCol(["resp","response","answer","key","button","choice"]);
ixMiss = findCol(["miss","omission","noresponse","nr"]);

% ---- RT ----
if ~isempty(ixRT)
    rt = toNumericCol(trialInfo.(trialInfo.Properties.VariableNames{ixRT}), n);
else
    rt = nan(n,1);
end
rt(rt < cfg.rt_min | rt > cfg.rt_max) = NaN;
behFeat.RT = rt;

% ---- Correct / Accuracy ----
if ~isempty(ixAcc)
    accRaw = trialInfo.(trialInfo.Properties.VariableNames{ixAcc});
    behFeat.Correct = toBinaryCol(accRaw, n);
else
    behFeat.Correct = nan(n,1);
end

% ---- Missing response flag (optional) ----
if ~isempty(ixMiss)
    missRaw = trialInfo.(trialInfo.Properties.VariableNames{ixMiss});
    behFeat.IsMissing = toBinaryCol(missRaw, n);
else
    behFeat.IsMissing = double(~isfinite(behFeat.RT)); % fallback
end

% ---- Response mapping ----
% RespCode: Left=1, Right=2, Up=3, Down=4, Other=0, Missing=NaN
if ~isempty(ixResp)
    respRaw = trialInfo.(trialInfo.Properties.VariableNames{ixResp});
    [respCode, flags] = mapResponse(respRaw, n);
else
    respCode = nan(n,1);
    flags = struct('Left',zeros(n,1),'Right',zeros(n,1),'Up',zeros(n,1), ...
                   'Down',zeros(n,1),'Other',zeros(n,1),'Missing',ones(n,1));
end

behFeat.RespCode    = respCode;
behFeat.RespLeft    = flags.Left;
behFeat.RespRight   = flags.Right;
behFeat.RespUp      = flags.Up;
behFeat.RespDown    = flags.Down;
behFeat.RespOther   = flags.Other;
behFeat.RespMissing = flags.Missing;

behFeat.TrialIndex = (1:n)';

end

%% ================= helpers =================
function TOut = pickTrialInfoIfStruct(TIn, cfg)
TOut = TIn;
if ~isstruct(TIn), return; end

ph = "";
if isfield(cfg,'phase') && ~isempty(cfg.phase), ph = lower(string(cfg.phase)); end

if ph ~= "" && isfield(TIn, ph)
    TOut = TIn.(ph);
    return;
end

f = fieldnames(TIn);
for i=1:numel(f)
    v = TIn.(f{i});
    if istable(v) && ~isempty(v)
        TOut = v;
        return;
    end
end

TOut = [];
end

function x = toNumericCol(col, n)
try
    if isnumeric(col) || islogical(col)
        x = double(col(:));
        x = padOrCrop(x,n);
        return;
    end

    if iscategorical(col), col = string(col); end

    if ischar(col)
        s = string(col);
        x = parseNumberFromString(s);
        x = repmat(x, n, 1);
        return;
    end

    if isstring(col)
        s = col(:);
        x = str2double(s);
        if all(isnan(x))
            x = parseNumberFromString(s);
        end
        x = padOrCrop(x,n);
        return;
    end

    if iscell(col)
        s = string(col(:));
        x = str2double(s);
        if all(isnan(x))
            x = parseNumberFromString(s);
        end
        x = padOrCrop(x,n);
        return;
    end
catch
end
x = nan(n,1);
end

function b = toBinaryCol(col, n)
x = toNumericCol(col, n);
if any(isfinite(x))
    b = x;
    b(b>1) = 1;
    b(b<0) = 0;
    return;
end

try
    if iscategorical(col), col = string(col); end
    if iscell(col), col = string(col); end
    if ischar(col), col = repmat(string(col), n, 1); end
    if isstring(col)
        s = lower(strtrim(col(:)));
        b = nan(n,1);
        b(ismember(s, ["1","true","t","yes","y","correct","hit","ok","success"])) = 1;
        b(ismember(s, ["0","false","f","no","n","incorrect","miss","wrong","error","fail"])) = 0;
        b = padOrCrop(b,n);
        return;
    end
catch
end
b = nan(n,1);
end

function [code, flags] = mapResponse(col, n)
code = toNumericCol(col, n);

flags = struct();
flags.Left    = zeros(n,1);
flags.Right   = zeros(n,1);
flags.Up      = zeros(n,1);
flags.Down    = zeros(n,1);
flags.Other   = zeros(n,1);
flags.Missing = zeros(n,1);

if any(isfinite(code))
    flags.Missing(~isfinite(code)) = 1;

    flags.Left(code==1)  = 1;
    flags.Right(code==2) = 1;
    flags.Up(code==3)    = 1;
    flags.Down(code==4)  = 1;

    other = isfinite(code) & ~(code==1 | code==2 | code==3 | code==4);
    flags.Other(other) = 1;
    code(other) = 0;

    return;
end

try
    if iscategorical(col), col = string(col); end
    if iscell(col), col = string(col); end
    if ischar(col), col = repmat(string(col), n, 1); end

    s = lower(strtrim(string(col(:))));
    code = nan(n,1);

    isMiss = (s=="" | s=="nan" | s=="none" | s=="null" | s=="missing");
    flags.Missing(isMiss) = 1;

    isL = contains(s, ["left","arrowleft","←"]) | ismember(s, ["l","lt","a","1"]);
    isR = contains(s, ["right","arrowright","→"]) | ismember(s, ["r","rt","d","2"]);
    isU = contains(s, ["up","arrowup","↑"]) | ismember(s, ["u","w","3"]);
    isD = contains(s, ["down","arrowdown","↓"]) | ismember(s, ["dn","s","4"]);

    flags.Left(isL & ~isMiss)  = 1;
    flags.Right(isR & ~isMiss) = 1;
    flags.Up(isU & ~isMiss)    = 1;
    flags.Down(isD & ~isMiss)  = 1;

    code(flags.Left==1)  = 1;
    code(flags.Right==1) = 2;
    code(flags.Up==1)    = 3;
    code(flags.Down==1)  = 4;

    other = ~isMiss & ~isL & ~isR & ~isU & ~isD;
    flags.Other(other) = 1;
    code(other) = 0;

catch
    code = nan(n,1);
    flags.Missing = ones(n,1);
end
end

function x = parseNumberFromString(s)
s = string(s(:));
x = nan(numel(s),1);
for i=1:numel(s)
    tok = regexp(s(i), '[-+]?\d*\.?\d+([eE][-+]?\d+)?', 'match', 'once');
    if ~isempty(tok)
        x(i) = str2double(tok);
    end
end
end

function y = padOrCrop(x,n)
y = nan(n,1);
m = min(n, numel(x));
y(1:m) = x(1:m);
end

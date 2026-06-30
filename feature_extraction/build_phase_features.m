function [featAll, trialInfo, meta] = build_phase_features(varargin)

nIn = numel(varargin);
if nIn < 11
    error('build_phase_features:NotEnoughInputs', 'Expected at least 11 inputs, got %d', nIn);
end

T             = varargin{1};
XttEEG        = varargin{2};
fs            = varargin{3};
chanLabels    = varargin{4};
phase         = varargin{5};
BANDS         = varargin{6};
baselineSec   = varargin{7};
epochSec      = varargin{8};
useHamming    = varargin{9};
REGIONS       = varargin{10};
addRegionFeats= varargin{11};

% optional cfg
cfg = struct();
if nIn >= 12 && isstruct(varargin{12})
    cfg = varargin{12};
else
    for k = 12:nIn
        if isstruct(varargin{k})
            cfg = varargin{k};
        end
    end
end


cfg = fill_defaults(cfg);
% --- validate T_phases ---
vn = T.Properties.VariableNames;
reqCols = {'TrialNum','Condition','StimSample','MaintSample','RetrSample'};
if ~all(ismember(reqCols, vn))
    missing = reqCols(~ismember(reqCols, vn));
    error('T_phases missing: %s', strjoin(cellstr(missing), ', '));
end

if isstring(phase), phase = char(phase); end
phase = lower(phase);
switch lower(phase)
    case 'stim',  startS = col_to_double(T.StimSample);
    case 'maint', startS = col_to_double(T.MaintSample);
    case 'retr',  startS = col_to_double(T.RetrSample);
    otherwise, error('Bad phase=%s', phase);
end

N = numel(startS);
labels = string(chanLabels(:));
Nc = numel(labels);
bandNames = fieldnames(BANDS);
Nb = numel(bandNames);

% --- EDF table -> matrix time x chan ---
X = edf_table_to_matrix(XttEEG);
X = double(X);
Tlen = size(X,1);
  % --- normalize labels to match number of data channels ---
  NcX = size(X,2);

  % If input was table/timetable and labels mismatch, try to use variable names of kept EEG columns
  if (istable(XttEEG) || istimetable(XttEEG)) && (numel(labels) ~= NcX)
      try
          vnames = string(XttEEG.Properties.VariableNames);
          nv = numel(vnames);
          keep = false(1,nv);
          for j=1:nv
              nm = lower(char(vnames(j)));
              if any(strcmp(nm, {'time','timestep','timestamp','datetime','date','times'}))
                  keep(j)=false; continue;
              end
              col = XttEEG{:,j};
              if isnumeric(col) || islogical(col)
                  keep(j)=true;
              elseif iscell(col)
                  k0 = find(~cellfun('isempty',col), 1, 'first');
                  if ~isempty(k0)
                      x0 = col{k0};
                      keep(j) = (isnumeric(x0) || islogical(x0));
                  end
              end
          end
          vnames = vnames(keep);
          if numel(vnames) == NcX
              labels = vnames(:);
          end
      catch
          
      end
  end

  if numel(labels) < NcX
      labels(end+1:NcX) = string(arrayfun(@(k) sprintf('ch%03d', k), (numel(labels)+1):NcX, 'UniformOutput', false));
  elseif numel(labels) > NcX
      labels = labels(1:NcX);
  end
  Nc = NcX;


trialInfo = T(:, {'TrialNum','Condition','StimSample','MaintSample','RetrSample'});
trialInfo.epochType = repmat(string(phase), height(trialInfo), 1);

meta = struct();
meta.fs = fs;
meta.bandDefs = BANDS;
meta.chanLabels = labels;
meta.baselineSec = baselineSec;
meta.useHamming  = useHamming;
meta.cfg = cfg;


  regionNames = string(REGIONS.names);
  regionOfChan = strings(Nc,1);
  for cc=1:Nc
      regOut = [];
      try
          regOut = REGIONS.fn(char(labels(cc)));
      catch
          regOut = [];
      end
      if isempty(regOut)
          regionOfChan(cc) = "Unknown";
      else
          regOut = string(regOut);
          if numel(regOut) > 1
              regionOfChan(cc) = join(regOut, "_");
          else
              regionOfChan(cc) = regOut;
          end
      end
  end
  meta.regionOfChan = regionOfChan;
  meta.regionNames  = regionNames;


% --- allocate blocks ---
epsPow = 1e-12;

freqMat = nan(N, Nb*Nc);
freqNames = strings(1, Nb*Nc);

statFns = {'mean','std','skew','kurt','median','iqr','mad','min','max','rms','lineLength'};
Ns = numel(statFns);
statMat = nan(N, Ns*Nc); statNames = strings(1, Ns*Nc);

tempFns = {'zcr','ssc','hj_activity','hj_mobility','hj_complexity'};
Nt = numel(tempFns);
tempMat = nan(N, Nt*Nc); tempNames = strings(1, Nt*Nc);

compFns = {'permEntropy3','lzComplexity'};
Nk = numel(compFns);
compMat = nan(N, Nk*Nc); compNames = strings(1, Nk*Nc);

tfBins = cfg.TF_BINS;
tfMat = []; tfNames = strings(1,0);
if cfg.DO_TF
    tfMat = nan(N, tfBins*Nb*Nc);
    tfNames = strings(1, tfBins*Nb*Nc);
end

% PSD accumulators (mean over trials)
psdFreq = [];
psdPxxSum = [];
psdCount = 0;

% --- loop trials ---
for i=1:N
    s0 = round(startS(i));
    if isnan(s0) || s0<1 || s0>Tlen
        continue;
    end

    if isscalar(epochSec)
        epSec = double(epochSec);
    else
        epSec = double(epochSec(i));
    end
    m = max(8, round(epSec * fs));

    b0 = max(1, s0 - round(baselineSec*fs));
    b1 = max(b0, s0-1);
    e1 = min(Tlen, s0 + m - 1);

    segRaw  = X(s0:e1, :);
    baseRaw = X(b0:b1, :);

    % baseline subtraction
    seg = segRaw;
    if ~isempty(baseRaw)
        seg = segRaw - nanmean_local(baseRaw, 1);
    end

    if useHamming
        seg = seg .* hamming(size(seg,1));
    end

    % -------- Frequency log bandpower --------
    col = 1;
    for bb=1:Nb
        bn = bandNames{bb};
        fr = BANDS.(bn);
        for cc=1:Nc
            val = bandpower_fft(seg(:,cc), fs, fr(1), fr(2));
            freqMat(i,col) = log10(val + epsPow);
            if i==1
                freqNames(col) = string(matlab.lang.makeValidName(['logbpblwin_' bn '_' char(labels(cc))]));
            end
            col = col + 1;
        end
    end

    % -------- Statistical --------
    if cfg.DO_STAT
        S = stat_block(seg);
        c = 1;
        for ss=1:Ns
            for cc=1:Nc
                statMat(i,c) = S(ss,cc);
                if i==1
                    statNames(c) = string(matlab.lang.makeValidName(['stat_' statFns{ss} '_' char(labels(cc))]));
                end
                c = c+1;
            end
        end
    end

    % -------- Temporal --------
    if cfg.DO_TEMP
        Tb = temp_block(seg);
        c = 1;
        for tt=1:Nt
            for cc=1:Nc
                tempMat(i,c) = Tb(tt,cc);
                if i==1
                    tempNames(c) = string(matlab.lang.makeValidName(['temp_' tempFns{tt} '_' char(labels(cc))]));
                end
                c = c+1;
            end
        end
    end

    % -------- Complexity --------
    if cfg.DO_COMP
        c = 1;
        for kk=1:Nk
            fn = compFns(kk);
            for cc=1:Nc
                x = seg(:,cc);
                if strcmp(fn,'permEntropy3')
                    compMat(i,c) = perm_entropy(x, 3, 1);
                else
                    compMat(i,c) = lz_complexity_binary(x);
                end
                if i==1
                    compNames(c) = string(matlab.lang.makeValidName(['comp_' fn '_' char(labels(cc))]));
                end
                c = c+1;
            end
        end
    end

    % -------- TF (optional) --------
    if cfg.DO_TF
        % baseline for normalization
        if ~isempty(baseRaw) && size(baseRaw,1) >= 8
            baseSig = baseRaw;
        else
            nBL = min(size(segRaw,1), max(8, round(baselineSec*fs)));
            baseSig = segRaw(1:nBL,:);
        end

        nT = size(segRaw,1);
        edges = round(linspace(1, nT+1, tfBins+1));

        colTF = 1;
        for bb=1:Nb
            bn = bandNames{bb};
            fr = BANDS.(bn);
            f0 = mean(fr);
            for cc=1:Nc
                p_db = morlet_power_db(segRaw(:,cc), baseSig(:,cc), fs, f0, cfg.TF_NCYCLES);
                for bi=1:tfBins
                    idx = edges(bi):(edges(bi+1)-1);
                    tfMat(i,colTF) = nanmean_local(p_db(idx));
                    if i==1
                        tfNames(colTF) = string(matlab.lang.makeValidName(['tf_morlet_' bn '_bin' num2str(bi) '_' char(labels(cc))]));
                    end
                    colTF = colTF + 1;
                end
            end
        end
    end

    % -------- PSD accumulate (optional) --------
    if cfg.SAVE_PSD
        [f, Pxx] = psd_multichan(seg, fs, cfg);
        if isempty(psdFreq)
            psdFreq = f;
            psdPxxSum = zeros(numel(f), Nc);
        end
        % align (just in case)
        L = min(numel(psdFreq), numel(f));
        psdPxxSum(1:L,:) = psdPxxSum(1:L,:) + Pxx(1:L,:);
        psdCount = psdCount + 1;
    end
end

% -------- assemble table --------
featAll = array2table(freqMat, 'VariableNames', cellstr(freqNames));
if cfg.DO_STAT
    featAll = [featAll, array2table(statMat, 'VariableNames', cellstr(statNames))];
end
if cfg.DO_TEMP
    featAll = [featAll, array2table(tempMat, 'VariableNames', cellstr(tempNames))];
end
if cfg.DO_COMP
    featAll = [featAll, array2table(compMat, 'VariableNames', cellstr(compNames))];
end
if cfg.DO_TF
    featAll = [featAll, array2table(tfMat, 'VariableNames', cellstr(tfNames))];
end

% -------- region averages --------
% -------- region averages --------
if addRegionFeats
    regNames = regionNames;
    % freq per band -> region mean
    featAll = [featAll, add_region_block(freqMat, 'logbpblwin', bandNames, regionOfChan, regNames, Nc)];
    if cfg.DO_STAT
        featAll = [featAll, add_region_block(statMat, 'stat', statFns, regionOfChan, regNames, Nc)];
    end
    if cfg.DO_TEMP
        featAll = [featAll, add_region_block(tempMat, 'temp', tempFns, regionOfChan, regNames, Nc)];
    end
    if cfg.DO_COMP
        featAll = [featAll, add_region_block(compMat, 'comp', compFns, regionOfChan, regNames, Nc)];
    end
    if cfg.DO_TF
        sub = cell(1, Nb*tfBins);
        k=1;
        for bb=1:Nb
            for bi=1:tfBins
                sub{k} = [bandNames{bb} '_bin' num2str(bi)];
                k=k+1;
            end
        end
        featAll = [featAll, add_region_block(tfMat, 'tf_morlet', sub, regionOfChan, regNames, Nc)];
    end
end

% -------- save PSD outputs in meta + png --------
if logical(cfg.SAVE_PSD) && psdCount > 0
    PxxMean = psdPxxSum ./ psdCount;
    meta.psd = struct();
    meta.psd.f = psdFreq;
    meta.psd.Pxx_mean = PxxMean;              % [freq x chan]
    meta.psd.method = cfg.PSD_METHOD;
    meta.psd.countTrials = psdCount;

    try
        save_psd_plot(psdFreq, PxxMean, labels, regionOfChan, regionNames, cfg);
    catch
        
    end
end

end

% ==================== helpers ====================
function cfg = fill_defaults(cfg)
def = struct('DO_STAT',true,'DO_TEMP',true,'DO_COMP',true,'DO_TF',false,...
             'TF_NCYCLES',6,'TF_BINS',3,'SAVE_PSD',false,'PSD_FMAX',80,'PSD_METHOD','pwelch',...
             'PSD_OUTDIR','','PSD_PREFIX','PSD');
fn = fieldnames(def);
for i=1:numel(fn)
    if ~isfield(cfg, fn{i}) || isempty(cfg.(fn{i}))
        cfg.(fn{i}) = def.(fn{i});
    end
end
end

function v = col_to_double(col)
if isnumeric(col) || islogical(col)
    v = double(col(:)); return;
end
if iscell(col)
    v = nan(numel(col),1);
    for i=1:numel(col)
        x = col{i};
        if isnumeric(x) || islogical(x)
            v(i) = double(x);
        else
            v(i) = str2double(char(x));
        end
    end
    return;
end
if isstring(col), col = cellstr(col); end
if ischar(col), col = cellstr(col); end
v = double(str2double(col(:)));
end

function X = edf_table_to_matrix(Xtt)

if isnumeric(Xtt) || islogical(Xtt)
    X = double(Xtt);
    return;
end

if istimetable(Xtt)
    try
        Xtt = timetable2table(Xtt, 'ConvertRowTimes', false);
    catch
        Xtt = timetable2table(Xtt);
    end
end

if ~istable(Xtt)
    X = double(Xtt);
    return;
end

vnames = string(Xtt.Properties.VariableNames);
nv = width(Xtt);
keep = false(1,nv);

for j=1:nv
    nm = lower(char(vnames(j)));
    if any(strcmp(nm, {'time','timestep','timestamp','datetime','date','times'}))
        keep(j) = false; continue;
    end

    col = Xtt{:,j};
    if isnumeric(col) || islogical(col)
        keep(j) = true;
    elseif iscell(col)
        k0 = find(~cellfun('isempty',col), 1, 'first');
        if ~isempty(k0)
            x0 = col{k0};
            keep(j) = (isnumeric(x0) || islogical(x0));
        end
    end
end

idxKeep = find(keep);
if isempty(idxKeep)
    error('build_phase_features:BadXtt', 'No numeric EEG columns found in XttEEG table/timetable.');
end

nChan = numel(idxKeep);
Xcols = cell(1, nChan);

for ii = 1:nChan
    j = idxKeep(ii);
    v = Xtt{:,j};

    if isnumeric(v) || islogical(v)
        Xcols{ii} = double(v(:));
    elseif iscell(v)
        parts = cell(numel(v),1);
        for k=1:numel(v)
            parts{k} = cell_elem_to_col(v{k});
        end
        try
            Xcols{ii} = double(vertcat(parts{:}));
        catch
            Xcols{ii} = double(cell2mat(parts));
            Xcols{ii} = Xcols{ii}(:);
        end
    else
        Xcols{ii} = nan(height(Xtt),1);
    end
end

L = cellfun(@numel, Xcols);
Lmax = max(L);
X = nan(Lmax, nChan);
for ii=1:nChan
    x = Xcols{ii};
    if isempty(x), continue; end
    X(1:numel(x), ii) = x(:);
end

if istable(Xtt) || istimetable(Xtt)
    try
        nv = width(Xtt);
        keep = false(1,nv);
        for j=1:nv
            try
                keep(j) = isnumeric(Xtt{:,j});
            catch
                keep(j) = false;
            end
        end
        if any(keep)
            X = Xtt{:, keep};
        else
            X = Xtt{:,:};
        end
    catch
        X = Xtt{:,:};
    end
    X = double(X);
    return;
end

error('build_phase_features:BadXtt', 'XttEEG must be numeric, table, or timetable.');

if ~istable(Xtt)
    X = double(Xtt); return;
end

nChan = width(Xtt);
Xcols = cell(1, nChan);

for c = 1:nChan
    v = Xtt{:,c};

    if isnumeric(v)
        Xcols{c} = double(v(:));
    elseif iscell(v)
        parts = cell(numel(v),1);
        for k=1:numel(v)
            parts{k} = cell_elem_to_col(v{k});
        end
        try
            Xcols{c} = double(vertcat(parts{:}));
        catch
            Xcols{c} = double(cell2mat(parts));
            Xcols{c} = Xcols{c}(:);
        end
    else
        try
            Xcols{c} = double(v(:));
        catch
            Xcols{c} = nan(height(Xtt),1);
        end
    end
end

L = cellfun(@numel, Xcols);
Lmax = max(L);
X = nan(Lmax, nChan);
for c=1:nChan
    x = Xcols{c};
    if isempty(x), continue; end
    X(1:numel(x), c) = x(:);
end
end

function y = cell_elem_to_col(x)
while iscell(x) && numel(x)==1
    x = x{1};
end
if isnumeric(x) || islogical(x)
    y = double(x(:));
elseif iscell(x)
    tmp = cell(numel(x),1);
    for i=1:numel(x)
        tmp{i} = cell_elem_to_col(x{i});
    end
    y = double(vertcat(tmp{:}));
else
    if isstring(x), x = cellstr(x); end
if ischar(x), x = cellstr(x); end
y = double(str2double(x(:)));
    y = y(:);
end
end

function S = stat_block(seg)
x = double(seg); x(isnan(x)) = 0;
mu = mean(x,1);
sd = std(x,0,1);
sk = skewness_simple(x);
ku = kurtosis_simple(x);
med = median(x,1);
iq = iqr_simple(x);
md = mad_simple(x);
mn = min(x,[],1);
mx = max(x,[],1);
rmsv = sqrt(mean(x.^2,1));
ll = sum(abs(diff(x,1,1)), 1);
S = [mu; sd; sk; ku; med; iq; md; mn; mx; rmsv; ll];
end

function Tb = temp_block(seg)
x = double(seg); x(isnan(x)) = 0;
zcr = zcr_simple(x);
ssc = ssc_simple(x);
[a,m,c] = hjorth_simple(x);
Tb = [zcr; ssc; a; m; c];
end

% ---------- simple stats (toolbox-free) ----------
function sk = skewness_simple(x)
mu = mean(x,1);
xc = x - mu;
m2 = mean(xc.^2,1) + eps;
m3 = mean(xc.^3,1);
sk = m3 ./ (m2.^(3/2));
end

function ku = kurtosis_simple(x)
mu = mean(x,1);
xc = x - mu;
m2 = mean(xc.^2,1) + eps;
m4 = mean(xc.^4,1);
ku = m4 ./ (m2.^2) - 3;
end

function iq = iqr_simple(x)
xs = sort(x,1);
n = size(xs,1);
q1 = xs(max(1, round(0.25*(n+1))), :);
q3 = xs(max(1, round(0.75*(n+1))), :);
iq = q3 - q1;
end

function md = mad_simple(x)
med = median(x,1);
md = median(abs(x - med), 1);
end

% ---------- temporal ----------
function z = zcr_simple(x)
s = sign(x); s(s==0) = 1;
z = sum(abs(diff(s,1,1))>0, 1) ./ max(1, size(x,1)-1);
end

function ssc = ssc_simple(x)
dx = diff(x,1,1);
ssc = sum((dx(1:end-1,:).*dx(2:end,:))<0, 1) ./ max(1, size(dx,1)-1);
end

function [a,m,c] = hjorth_simple(x)
dx = diff(x,1,1);
ddx = diff(dx,1,1);
a = var(x,0,1);
m = sqrt(var(dx,0,1) ./ (var(x,0,1) + eps));
m2 = sqrt(var(ddx,0,1) ./ (var(dx,0,1) + eps));
c = m2 ./ (m + eps);
end

% ---------- complexity ----------
function pe = perm_entropy(x, m, tau)
x = double(x(:)); x(isnan(x)) = 0;
n = numel(x) - (m-1)*tau;
if n<=0, pe = NaN; return; end
patterns = zeros(n, m);
for i=1:m
    patterns(:,i) = x( (1:n) + (i-1)*tau );
end
[~, ord] = sort(patterns, 2);
base = (m.^(0:m-1));
keys = sum((ord-1).*base, 2);
[~,~,ic] = unique(keys);
cnt = accumarray(ic, 1);
p = cnt / sum(cnt);
pe = -sum(p .* log(p+eps));
pe = pe / log(numel(cnt)+eps);
end

function c = lz_complexity_binary(x)
x = double(x(:)); x(isnan(x)) = 0;
thr = median(x);
b = x > thr;
s = char(b.' + '0'); %#ok<CHARTEN>
n = length(s);
i = 1; k = 1; l = 1; c0 = 1;
while true
    if i+k > n
        c0 = c0 + 1;
        break;
    end
    sub1 = s(l:l+k-1);
    sub2 = s(i:i+k-1);
    if strcmp(sub1, sub2)
        k = k + 1;
        if l + k - 1 > i
            l = l + 1;
            k = 1;
        end
    else
        i = i + 1;
        if i == l
            c0 = c0 + 1;
            l = l + k;
            if l > n
                break;
            end
            i = 1; k = 1;
        end
    end
end
c = c0 * log2(n) / n;
end

% ---------- TF ----------
function p_db = morlet_power_db(sig, baseSig, fs, f0, nCycles)
sig = double(sig(:)); sig(isnan(sig)) = 0;
baseSig = double(baseSig(:)); baseSig(isnan(baseSig)) = 0;

t = (-2:1/fs:2);
s = nCycles/(2*pi*f0);
wave = exp(2*1i*pi*f0.*t) .* exp(-(t.^2)./(2*s^2));

nTime = numel(sig);
nWave = numel(wave);
nConv = nTime + nWave - 1;

sigX = fft(sig, nConv);
waveX = fft(wave, nConv);
convRes = ifft(sigX .* waveX);
convRes = convRes(floor(nWave/2)+1 : floor(nWave/2)+nTime);
pow = abs(convRes).^2;

nB = numel(baseSig);
nConvB = nB + nWave - 1;
baseX = fft(baseSig, nConvB);
waveXB = fft(wave, nConvB);
convB = ifft(baseX .* waveXB);
convB = convB(floor(nWave/2)+1 : floor(nWave/2)+nB);
powB = nanmean_local(abs(convB).^2);

p_db = 10*log10( pow ./ (powB + eps) );
end

% ---------- PSD helpers ----------
function [f, Pxx] = psd_multichan(seg, fs, cfg)
% returns f [F x 1], Pxx [F x Nc]
Nc = size(seg,2);
fmax = cfg.PSD_FMAX;

if strcmpi(cfg.PSD_METHOD,'pwelch') && exist('pwelch','file')==2
    % choose window
    n = size(seg,1);
    win = hamming(min(n, max(64, round(n/2))));
    nover = floor(numel(win)/2);
    nfft = max(256, 2^nextpow2(numel(win)));
    Pxx = [];
    for c=1:Nc
        [P,f] = pwelch(seg(:,c), win, nover, nfft, fs);
        if isempty(Pxx), Pxx = zeros(numel(P),Nc); end
        Pxx(:,c) = P;
    end
else
    % FFT periodogram
    n = size(seg,1);
    nfft = max(256, 2^nextpow2(n));
    f = fs*(0:floor(nfft/2))/nfft;
    Pxx = zeros(numel(f), Nc);
    for c=1:Nc
        x = seg(:,c); x = x - nanmean_local(x); x(isnan(x)) = 0;
        X = fft(x, nfft);
        P2 = abs(X/nfft).^2;
        P1 = P2(1:floor(nfft/2)+1);
        Pxx(:,c) = P1;
    end
end

% trim fmax
idx = f <= fmax;
f = f(idx);
Pxx = Pxx(idx,:);
end

function save_psd_plot(f, Pxx, labels, regionOfChan, regionNames, cfg)
% Plot region-mean PSD (dB/Hz) and save png
outDir = cfg.PSD_OUTDIR; if isstring(outDir), outDir = char(outDir); end
if isempty(outDir), return; end
prefix = cfg.PSD_PREFIX; if isstring(prefix), prefix = char(prefix); end
if isempty(prefix), prefix = 'PSD'; end

% region mean
Nr = numel(regionNames);
Preg = nan(numel(f), Nr);
for r=1:Nr
    idx = regionOfChan==regionNames(r);
    if any(idx)
        Preg(:,r) = nanmean_local(Pxx(:,idx), 2);
    end
end

h = figure('Visible','off');
plot(f, 10*log10(Preg + eps));
grid on;
xlabel('Hz'); ylabel('PSD (dB/Hz)');
title([prefix ' | Region-mean PSD']);
legend(cellstr(regionNames), 'Location','northeastoutside');

outPng = fullfile(outDir, [prefix '_PSD.png']);
try
    exportgraphics(h, outPng, 'Resolution', 200);
catch
    saveas(h, outPng);
end
close(h);
end

% ---------- bandpower via FFT ----------
function p = bandpower_fft(x, fs, f1, f2)
x = double(x(:));
if numel(x) < 4, p = 0; return; end
x = x - nanmean_local(x);
x(isnan(x)) = 0;
n = numel(x);
nfft = 2^nextpow2(n);
X = fft(x, nfft);
P2 = abs(X/nfft).^2;
P1 = P2(1:floor(nfft/2)+1);
f = fs*(0:floor(nfft/2))/nfft;
idx = (f>=f1) & (f<=f2);
if ~any(idx), p = 0; else, p = trapz(f(idx), P1(idx)); end
end

function regTbl = add_region_block(blockMat, basePrefix, subNames, regionOfChan, regNames, Nc)
  % Add region-averaged versions of per-channel blocks.
  % blockMat: [nEpoch x (Nc*nSub)] with per-channel features concatenated by subfeature.
  basePrefix   = char(basePrefix);
  regionOfChan = string(regionOfChan);
  regNames     = string(regNames);
  nEpoch = size(blockMat,1);
  nSub = floor(size(blockMat,2) / Nc);
  regTbl = table();
  for kk2 = 1:nSub
      cols = (kk2-1)*Nc + (1:Nc);
      if cols(end) > size(blockMat,2)
          break; % safety
      end
      subMat = blockMat(:, cols);
      if isstring(subNames)
          subName = char(subNames(kk2));
      else
          subName = subNames{kk2};
      end
      for rr = 1:numel(regNames)
          rname = regNames(rr);
          idx = (regionOfChan == rname);
          if any(idx)
              regFeat = nanmean(subMat(:, idx), 2);
          else
              regFeat = nan(nEpoch, 1);
          end
          vname = matlab.lang.makeValidName([basePrefix '_' subName '_reg_' char(rname)]);
          regTbl.(vname) = regFeat;
      end
  end
end

function m = nanmean_local(x, dim)
% toolbox-free nanmean replacement
if nargin < 2 || isempty(dim)
    % first non-singleton
    dim = find(size(x)~=1, 1, 'first');
    if isempty(dim), dim = 1; end
end
x = double(x);
mask = ~isnan(x);
x(~mask) = 0;
cnt = sum(mask, dim);
cnt(cnt==0) = NaN;
m = sum(x, dim) ./ cnt;
end

function v = nanvar_local(x, dim)
% toolbox-free nanvar (normalization by N-1 like var)
if nargin < 2 || isempty(dim)
    dim = find(size(x)~=1, 1, 'first');
    if isempty(dim), dim = 1; end
end
m = nanmean_local(x, dim);
% expand m to x size
rep = ones(1, ndims(x)); rep(dim) = size(x, dim);
mexp = repmat(m, rep);
xc = x - mexp;
xc(isnan(xc)) = 0;
mask = ~isnan(x);
cnt = sum(mask, dim);
cntm1 = cnt - 1;
cntm1(cntm1<=0) = NaN;
v = sum(xc.^2, dim) ./ cntm1;
end

function s = nanstd_local(x, dim)
s = sqrt(nanvar_local(x, dim));
end
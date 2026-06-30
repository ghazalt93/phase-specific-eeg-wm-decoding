function freqFeat = WM_BandpowerFeatures_FFT_FAST(ep, fs, t, useChanIdx, cfg)


% bands
if isfield(cfg,'freqBands') && ~isempty(cfg.freqBands)
    bands = cfg.freqBands;
else
    bands = struct('theta',[4 7],'alpha',[8 12],'beta',[13 30]);
end
bnames = fieldnames(bands);
nB = numel(bnames);

nTrial = size(ep,3);
useChanIdx = useChanIdx(:)';
useChanIdx = useChanIdx(useChanIdx>=1 & useChanIdx<=size(ep,1));
nCh = numel(useChanIdx);

% use post-stim only
idxUse = (t >= 0);
if ~any(idxUse), idxUse = true(size(t)); end

N = nnz(idxUse);
if N < 16
    % too short -> all NaN
    out = nan(nTrial, nB*nCh);
    varNames = makeNames(bnames, useChanIdx);
    freqFeat = array2table(out, 'VariableNames', varNames);
    return
end

% freq axis (one-sided)
nHalf = floor(N/2);
f = fs*(0:nHalf)/N;

% frequency resolution (constant for this FFT grid)
df = fs / N;

% band indices precompute
idxF = cell(nB,1);
for ib=1:nB
    fr = bands.(bnames{ib});
    idxF{ib} = (f >= fr(1)) & (f <= fr(2));
end

out = nan(nTrial, nB*nCh);

for tr = 1:nTrial
    x = double(ep(useChanIdx, idxUse, tr));  % [nCh x N]
    x(~isfinite(x)) = 0;
    x = x - mean(x,2);  % detrend mean per channel

    X = fft(x, [], 2);
    P2 = (abs(X)/N).^2;           % [nCh x N]
    P1 = P2(:, 1:nHalf+1);        % one-sided
    if size(P1,2) > 2
        P1(:,2:end-1) = 2*P1(:,2:end-1);
    end

    for ib=1:nB
        msk = idxF{ib};
        if ~any(msk)
            bp = nan(nCh,1);
        else
            % FIX: discrete band "area" over FFT grid
            % (works even if sum(msk)==1)
            bp = sum(P1(:,msk), 2) * df;   % [nCh x 1]
        end
        col0 = (ib-1)*nCh;
        out(tr, col0+(1:nCh)) = bp(:);
    end
end

varNames = makeNames(bnames, useChanIdx);
freqFeat = array2table(out, 'VariableNames', varNames);

end

function varNames = makeNames(bnames, useChanIdx)
nB = numel(bnames);
nCh = numel(useChanIdx);
varNames = cell(1, nB*nCh);
k = 0;
for ib=1:nB
    for ic=1:nCh
        k = k+1;
        varNames{k} = sprintf('Pow_%s_ch%02d', bnames{ib}, useChanIdx(ic));
    end
end
varNames = matlab.lang.makeUniqueStrings(varNames);
end

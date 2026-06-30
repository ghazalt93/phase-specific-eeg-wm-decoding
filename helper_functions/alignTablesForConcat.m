function [A,B] = alignTablesForConcat(A,B)
% alignTablesForConcat
% Ensures A and B have identical variables (names + compatible types)
% so you can safely do: A = [A; B];

varsA = string(A.Properties.VariableNames);
varsB = string(B.Properties.VariableNames);

% ---- add missing columns into B (based on A types) ----
missInB = setdiff(varsA, varsB);
for v = missInB
    refCol = A.(v);
    B.(v)  = makeFillerLike(refCol, height(B));
end

% ---- add missing columns into A (based on B types) ----
missInA = setdiff(varsB, varsA);
for v = missInA
    refCol = B.(v);
    A.(v)  = makeFillerLike(refCol, height(A));
end

% ---- match column order ----
B = B(:, A.Properties.VariableNames);

end

% =====================================================================
function filler = makeFillerLike(refCol, nRows)
% Create a filler column matching the type of refCol

% string
if isstring(refCol)
    filler = repmat(missing, nRows, 1);

% cellstr / cell
elseif iscell(refCol)
    filler = repmat({''}, nRows, 1);

% categorical
elseif iscategorical(refCol)
    filler = categorical(repmat({''}, nRows, 1));

% logical
elseif islogical(refCol)
    filler = false(nRows, 1);

% numeric
elseif isnumeric(refCol)
    filler = nan(nRows, 1, 'like', refCol);

% datetime
elseif isdatetime(refCol)
    filler = repmat(NaT, nRows, 1);

% duration
elseif isduration(refCol)
    filler = repmat(duration(NaN,0,0), nRows, 1);

else
    % fallback
    filler = nan(nRows,1);
end

end

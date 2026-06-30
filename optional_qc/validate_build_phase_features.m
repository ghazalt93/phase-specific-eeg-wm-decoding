function validate_build_phase_features()
%VALIDATE_BUILD_PHASE_FEATURES Quick sanity checks for build_phase_features.m on MATLAB path.
p = which('build_phase_features');
fprintf('build_phase_features path: %s\n', p);

% show first non-empty non-comment line
fid = fopen(p,'r');
firstLine = "";
while true
    t = fgetl(fid);
    if ~ischar(t), break; end
    tt = strtrim(t);
    if isempty(tt), continue; end
    if startsWith(tt,'%'), continue; end
    firstLine = t; break;
end
fclose(fid);

fprintf('First non-empty line: %s\n', firstLine);
try
    n = nargin('build_phase_features');
    fprintf('nargin(build_phase_features) = %d\n', n);
catch ME
    fprintf('Could not evaluate nargin: %s\n', ME.message);
end

% show any nested function definitions before the first "end" (common cause of MATLAB:m_function_def)
txt = fileread(p);
lines = splitlines(txt);
idxEnd = find(strcmp(strtrim(lines),'end'), 1, 'first');
if isempty(idxEnd), idxEnd = numel(lines); end
bad = [];
for i=2:idxEnd-1
    if startsWith(strtrim(lines{i}), 'function ')
        bad(end+1) = i; %#ok<AGROW>
    end
end
if ~isempty(bad)
    fprintf('[WARN] Found function definitions before first end at lines: %s\n', mat2str(bad));
    fprintf('       This will trigger "Function definition is misplaced".\n');
else
    fprintf('[OK] No nested function definitions before the first end.\n');
end
end

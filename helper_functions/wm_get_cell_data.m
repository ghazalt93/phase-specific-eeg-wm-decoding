function [X, y, subj, rows] = wm_get_cell_data(D, meta, featVars, cfg, classKeep)
    yAll = double(D.(meta.label));
    runRaw = D.(meta.run);

    if isnumeric(runRaw)
        runAll = double(runRaw);
    else
        runStr = string(runRaw);
        runAll = NaN(size(runStr));
        for i = 1:numel(runStr)
            tok = regexp(char(runStr(i)), '\d+', 'match');
            if ~isempty(tok)
                runAll(i) = str2double(tok{end});
            end
        end
    end

    phaseAll = lower(string(D.(meta.phase)));
    subjAll  = string(D.(meta.subject));

    rows = runAll == cfg.run & phaseAll == lower(string(cfg.phase)) & isfinite(yAll);

    if nargin >= 6 && ~isempty(classKeep)
        rows = rows & ismember(yAll, classKeep);
    else
        rows = rows & ismember(yAll, [1 2 3]);
    end

    if isfield(cfg,'excludeSubjects') && ~isempty(cfg.excludeSubjects)
        for i = 1:numel(cfg.excludeSubjects)
            rows = rows & ~strcmpi(subjAll, string(cfg.excludeSubjects{i}));
        end
    end

    X = double(table2array(D(rows, featVars)));
    y = yAll(rows);
    subj = cellstr(subjAll(rows));
end

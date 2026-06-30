function [featVars, featCh] = wm_detect_features(D, meta, cfg)
    vars = D.Properties.VariableNames;
    low  = lower(vars);
    n = height(D);

    exclude = false(1,numel(vars));
    must = {meta.subject, meta.run, meta.phase, meta.label, ...
        'ycondition','condition','cond','ycorrect','correct','trial','trialnum','patternid'};

    for i = 1:numel(vars)
        for j = 1:numel(must)
            if ~isempty(must{j}) && strcmpi(vars{i}, must{j})
                exclude(i) = true;
            end
        end

        if cfg.excludeBehavior
            pats = {'rt','reaction','correct','accuracy','error','resp','response','missing','ismissing','behavior','behav'};
            for p = 1:numel(pats)
                if contains(low{i}, pats{p})
                    exclude(i) = true;
                end
            end
        end
    end

    isNum = false(1,numel(vars));
    for i = 1:numel(vars)
        x = D.(vars{i});
        isNum(i) = isnumeric(x) && isvector(x) && numel(x)==n;
    end

    idx = find(isNum & ~exclude);
    featVars = vars(idx);
    featCh = NaN(1,numel(featVars));

    for i = 1:numel(featVars)
        featCh(i) = wm_get_ch(featVars{i});
    end

    keep = ~isnan(featCh);
    featVars = featVars(keep);
    featCh = featCh(keep);
end

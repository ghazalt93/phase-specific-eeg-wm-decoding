function out = wm_pick_col(v,l,cands)
    out = '';
    for i = 1:numel(cands)
        k = find(strcmp(l, lower(cands{i})), 1);
        if ~isempty(k), out = v{k}; return; end
    end
    for i = 1:numel(cands)
        k = find(contains(l, lower(cands{i})), 1);
        if ~isempty(k), out = v{k}; return; end
    end
end

function D = wm_get_table_from_mat(S)
    if isfield(S,'DS'), D = S.DS; return; end
    if isfield(S,'dataset'), D = S.dataset; return; end
    if isfield(S,'D'), D = S.D; return; end
    if isfield(S,'T'), D = S.T; return; end

    fn = fieldnames(S);
    for i = 1:numel(fn)
        if istable(S.(fn{i}))
            D = S.(fn{i});
            return;
        end
    end
    error('No dataset table found in MAT file.');
end

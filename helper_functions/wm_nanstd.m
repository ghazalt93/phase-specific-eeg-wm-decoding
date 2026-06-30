function s = wm_nanstd(X)
    s = NaN(1,size(X,2));
    for j = 1:size(X,2)
        x = X(:,j); x = x(isfinite(x));
        if ~isempty(x), s(j) = std(x); end
    end
end

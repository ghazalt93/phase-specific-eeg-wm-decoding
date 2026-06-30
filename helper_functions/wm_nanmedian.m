function m = wm_nanmedian(X)
    m = NaN(1,size(X,2));
    for j = 1:size(X,2)
        x = X(:,j); x = x(isfinite(x));
        if ~isempty(x), m(j) = median(x); end
    end
end

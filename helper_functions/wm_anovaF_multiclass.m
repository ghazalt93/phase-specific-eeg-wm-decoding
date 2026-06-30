function F = wm_anovaF_multiclass(X,y)
    classes = unique(y(isfinite(y)));
    n = size(X,1);
    p = size(X,2);
    grand = mean(X,1);

    SSb = zeros(1,p);
    SSw = zeros(1,p);
    k = numel(classes);

    for ci = 1:k
        idx = y == classes(ci);
        ni = sum(idx);
        if ni < 1, continue; end
        mu = mean(X(idx,:),1);
        SSb = SSb + ni*(mu-grand).^2;
        D = bsxfun(@minus, X(idx,:), mu);
        SSw = SSw + sum(D.^2,1);
    end

    F = (SSb/max(k-1,1)) ./ max(SSw/max(n-k,1), eps);
end

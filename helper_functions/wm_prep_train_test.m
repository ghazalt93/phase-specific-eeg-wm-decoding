function [Xtr, Xte, keepCol] = wm_prep_train_test(Xtr0, Xte0)
    med = wm_nanmedian(Xtr0);

    for j = 1:size(Xtr0,2)
        if ~isfinite(med(j)), med(j) = 0; end
        bad = ~isfinite(Xtr0(:,j)); if any(bad), Xtr0(bad,j) = med(j); end
        bad = ~isfinite(Xte0(:,j)); if any(bad), Xte0(bad,j) = med(j); end
    end

    sd0 = std(Xtr0,0,1);
    keep = isfinite(sd0) & sd0 > 1e-12;
    keepCol = find(keep);

    Xtr0 = Xtr0(:,keep);
    Xte0 = Xte0(:,keep);

    mu = mean(Xtr0,1);
    sd = std(Xtr0,0,1);
    sd(~isfinite(sd) | sd < 1e-12) = 1;

    Xtr = bsxfun(@rdivide, bsxfun(@minus, Xtr0, mu), sd);
    Xte = bsxfun(@rdivide, bsxfun(@minus, Xte0, mu), sd);

    Xtr(~isfinite(Xtr)) = 0;
    Xte(~isfinite(Xte)) = 0;
end

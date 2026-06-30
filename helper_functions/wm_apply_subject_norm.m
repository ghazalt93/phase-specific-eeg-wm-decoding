function Xout = wm_apply_subject_norm(X, subj, method)
    Xout = X;
    subjects = unique(subj, 'stable');

    for i = 1:numel(subjects)
        idx = strcmp(subj, subjects{i});
        Xi = X(idx,:);

        switch lower(method)
            case 'subjectz'
                mu = wm_nanmean(Xi);
                sd = wm_nanstd(Xi);
                sd(~isfinite(sd) | sd < 1e-12) = 1;
                Xout(idx,:) = bsxfun(@rdivide, bsxfun(@minus, Xi, mu), sd);
            case 'none'
                Xout(idx,:) = Xi;
            otherwise
                error('Unknown normalization method: %s', method);
        end
    end

    Xout(~isfinite(Xout)) = 0;
end

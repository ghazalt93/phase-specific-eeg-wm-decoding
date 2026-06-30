function CI = wm_bootstrap_ci(x, nBoot, alpha)
    x = x(isfinite(x));
    if isempty(x)
        CI = [NaN NaN NaN];
        return;
    end
    rng(20260507);
    boot = NaN(nBoot,1);
    n = numel(x);
    for b = 1:nBoot
        idx = randi(n, [n 1]);
        boot(b) = mean(x(idx));
    end
    CI = [mean(x), prctile(boot, 100*alpha/2), prctile(boot, 100*(1-alpha/2))];
end

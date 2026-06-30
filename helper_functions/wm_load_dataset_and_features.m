function [D, meta, featVars, featCh] = wm_load_dataset_and_features(cfg)
    if exist(cfg.outDir,'dir') ~= 7
        mkdir(cfg.outDir);
    end

    if exist(cfg.DSpath,'file') ~= 2
        error('Dataset not found: %s', cfg.DSpath);
    end

    S = load(cfg.DSpath);
    D = wm_get_table_from_mat(S);
    meta = wm_detect_meta(D);
    [featVars, featCh] = wm_detect_features(D, meta, cfg);

    fmask = ismember(featCh, cfg.channels);
    featVars = featVars(fmask);
    featCh   = featCh(fmask);

    if isempty(featVars)
        error('No feature variables selected. Check cfg.channels and feature names.');
    end
end

function cfg = wm_default_cfg()
cfg = get_project_config();
    cfg = struct();

    % ===== EDIT THESE PATHS IF NEEDED =====
    cfg.ROOT   = fullfile(cfg.dataRoot, 'Subjects');
    cfg.DSpath = fullfile(cfg.ROOT, '_wm_ml', 'dataset.mat');

    cfg.outDir = fullfile(cfg.outputRoot, '_wm_reports', 'posthoc_steps_3_6');

    % Main cell to summarize in the paper.
    cfg.run   = 3;
    cfg.phase = 'retr';

    cfg.nFolds = 10;
    cfg.Kfeat  = 100;
    cfg.channels = 5:60;
    cfg.channelSetName = 'noEdge_5_60';
    cfg.normMethod = 'subjectZ';
    cfg.excludeBehavior = true;
    cfg.randomSeed = 20260507;

    % If primary permutation included s21, keep this empty.
    % If final analysis must exclude s21, set {'s21'}.
    cfg.excludeSubjects = {};

    cfg.binaryPairs = {
        1, 2, 'color_vs_orientation'
        1, 3, 'color_vs_conjunction'
        2, 3, 'orientation_vs_conjunction'
    };

    cfg.nBoot = 5000;
    cfg.alpha = 0.05;
end

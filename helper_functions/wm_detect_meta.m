function meta = wm_detect_meta(T)
    v = T.Properties.VariableNames;
    l = lower(v);

    meta.subject = wm_pick_col(v,l,{'subject','subj','subjid','participant','participantid'});
    meta.run     = wm_pick_col(v,l,{'run','runnum','runid','session','sessionid','sess'});
    meta.phase   = wm_pick_col(v,l,{'phase','phasename','epochphase'});
    meta.label   = wm_pick_col(v,l,{'ycondition','condition','cond','y','label','class'});
    meta.trial   = wm_pick_col(v,l,{'trialnum','trial','trialid'});

    if isempty(meta.subject), error('Subject column not found.'); end
    if isempty(meta.run), error('Run column not found.'); end
    if isempty(meta.phase), error('Phase column not found.'); end
    if isempty(meta.label), error('Condition label column not found.'); end
end

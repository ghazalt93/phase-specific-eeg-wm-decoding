function fam = wm_feature_family(name)
    s = lower(string(name));

    if startsWith(s,"bp_") || startsWith(s,"rbp_") || contains(s,"bandpower")
        fam = "bandpower";
    elseif contains(s,"lzc") || contains(s,"entropy") || contains(s,"complex")
        fam = "complexity";
    elseif contains(s,"aucabs") || contains(s,"linelen") || contains(s,"tmaxabs")
        fam = "waveformDescriptor";
    elseif contains(s,"rms") || contains(s,"peak") || contains(s,"mean") || contains(s,"std") || ...
            contains(s,"p2p") || contains(s,"deriv") || contains(s,"tkeo") || contains(s,"var") || ...
            contains(s,"skew") || contains(s,"kurt")
        fam = "temporalStat";
    else
        fam = "other";
    end
end

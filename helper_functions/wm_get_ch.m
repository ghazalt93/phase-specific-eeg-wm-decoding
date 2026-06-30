function ch = wm_get_ch(name)
    ch = NaN;
    tok = regexp(name, 'ch[_-]?(\d{1,2})', 'tokens', 'once', 'ignorecase');
    if isempty(tok)
        tok = regexp(name, 'chan(?:nel)?[_-]?(\d{1,2})', 'tokens', 'once', 'ignorecase');
    end
    if ~isempty(tok)
        ch = str2double(tok{1});
        if ch < 1 || ch > 64
            ch = NaN;
        end
    end
end

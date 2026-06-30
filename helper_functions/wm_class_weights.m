function W = wm_class_weights(y)
    classes = unique(y(:))';
    W = zeros(size(y));
    for i = 1:numel(classes)
        idx = y == classes(i);
        W(idx) = 1 / max(sum(idx),1);
    end
    W = W / mean(W);
end

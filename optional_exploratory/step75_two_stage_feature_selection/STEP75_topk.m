function idxSel = STEP75_topk(score, Kfeat)
score2 = double(score(:));
score2(~isfinite(score2)) = -Inf;

[~, ord] = sort(score2, 'descend');
ord = ord(isfinite(score2(ord)));

K = min(Kfeat, numel(ord));
idxSel = ord(1:K);
end

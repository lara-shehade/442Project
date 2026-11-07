function xq = quan(x, thr, lvl)
% Quantizes x with thresholds thr (M-1 values, increasing) and levels lvl (M values).
% x <= thr(1) -> lvl(1);  thr(j-1) < x <= thr(j) -> lvl(j);  x > thr(end) -> lvl(end).
idx = ones(size(x));
for j = 1:numel(thr)
    idx = idx + (x > thr(j));
end
xq = reshape(lvl(idx), size(x));
end

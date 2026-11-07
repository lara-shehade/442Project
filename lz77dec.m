function s = lz77dec(tok)
s = [];
for i = 1:size(tok, 1)
    for k = 1:tok(i,2), s(end+1) = s(end - tok(i,1) + 1); end
    s(end+1) = tok(i,3);
end
end

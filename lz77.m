function tok = lz77(s, W, L)
% LZ77: each token = [offset, length, next symbol]. The longest match of the
% upcoming symbols is searched in the last W symbols (match length at most L).
N = numel(s);  pos = 1;  tok = zeros(0, 3);
while pos <= N
    Lmax = min(L, N - pos);                  % keep one symbol for the "next symbol" field
    best = 0;  off = 0;
    j = max(1, pos - W) : pos - 1;           % candidate starts in the window
    if Lmax > 0 && ~isempty(j)
        idx = j(:) + (0:Lmax-1);             % candidate substrings (rows)
        eq  = s(idx) == s(pos + (0:Lmax-1)); % compare with the upcoming symbols
        len = sum(cumprod(eq, 2), 2);        % length of the match of each candidate
        [best, b] = max(len);
        if best > 0, off = pos - j(b); end
    end
    tok(end+1, :) = [off, best, s(pos + best)];
    pos = pos + best + 1;
end
end

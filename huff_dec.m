function s = huff_dec(bits, codes)
% Decode a bit string with a prefix-free code: read bits until they form a codeword.
lut = containers.Map(codes, num2cell(1:numel(codes)));
s = zeros(1, numel(bits));  n = 0;  cur = '';
for b = bits
    cur = [cur b];
    if isKey(lut, cur)
        n = n + 1;  s(n) = lut(cur);  cur = '';
    end
end
s = s(1:n);
end

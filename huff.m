function codes = huff(p)
% Huffman code for probabilities p: codes{i} is the codeword (char '0'/'1') of symbol i.
% Merge the two least likely nodes, prepend 0 / 1 to their symbols, repeat.
n = numel(p);
codes = repmat({''}, 1, n);
if n == 1, codes = {'0'}; return; end
groups = num2cell(1:n);
w = p(:).';
while numel(w) > 1
    [w, o] = sort(w);  groups = groups(o);
    for i = groups{1}, codes{i} = ['0' codes{i}]; end
    for i = groups{2}, codes{i} = ['1' codes{i}]; end
    groups = [{[groups{1} groups{2}]}, groups(3:end)];
    w = [w(1) + w(2), w(3:end)];
end
end

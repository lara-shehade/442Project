clear; clc; close all;
%% 3 - Source coding of the quantized sequences
load('quantized.mat');                 % from part2.m: xs, fs, thr, lvl, names, mse_s, mse_r
N  = numel(xs);
nq = numel(names);
K  = [2 3 4];                          % block lengths
W  = 255;  Lz_list = [15 31 63 127];  % LZ77 window and maximum match lengths tried

for q = 1:nq
    M = numel(lvl{q});
    s = quan(xs, thr{q}, 1:M);         % symbol stream: index of the level, alphabet A = {1..M}

    % (a) no model: fixed length
    Lfix(q) = ceil(log2(M));

    % (b) empirical distribution and entropy
    cnt = accumarray(s(:), 1, [M 1]).';  cnts{q} = cnt;   % count of each symbol
    p = cnt / N;
    H(q) = -sum(p(p > 0) .* log2(p(p > 0)));
    same(q) = mean(s(2:end) == s(1:end-1));                % (h) next symbol equal to the current one

    % (c) Huffman code, encode, decode, check
    codes = huff(p);
    bits  = [codes{s}];
    ok(q) = isequal(huff_dec(bits, codes), s);
    Lh(q) = numel(bits) / N;           % = sum p(a) l(a)

    % (e)-(g) blocks of k symbols: fixed length k*Lfix, Huffman on p(a^k)
    for j = 1:numel(K)
        k = K(j);
        blk = reshape(s, k, []).';                 % one block per row
        [ub, ~, b] = unique(blk, 'rows');          % different blocks and block index
        b  = b(:).';
        pb = accumarray(b(:), 1).' / numel(b);
        Hk(q,j) = -sum(pb .* log2(pb)) / k;        % H_k / k
        cb = huff(pb);
        bb = [cb{b}];
        okb(q,j) = isequal(huff_dec(bb, cb), b);
        Lk(q,j) = numel(bb) / N;                   % bits per symbol
        nblk(q,j) = size(ub, 1);                   % (h) number of different blocks
        once(q,j) = sum(accumarray(b(:), 1) == 1); % (h) blocks seen only once
        if q == 1 && k == 2, ex_ub = ub; ex_pb = pb; ex_cb = cb; end   % (f) example
    end

    % (k) LZ77, no model: try several maximum match lengths L, keep the best
    for i = 1:numel(Lz_list)
        tok = lz77(s, W, Lz_list(i));
        okz(q,i) = isequal(lz77dec(tok), s);
        Lz_all(q,i) = size(tok,1) * (ceil(log2(W+1)) + ceil(log2(Lz_list(i)+1)) + Lfix(q)) / N;
    end
    [Llz(q), ib] = min(Lz_all(q,:));   Lbest(q) = Lz_list(ib);
end

%% Tables
fprintf('\n(a)-(d) Symbol-wise coding (bits/symbol)\n');
fprintf('%-15s %4s %7s %8s %8s %8s %9s\n', 'Quantizer', 'M', 'Fixed', 'H(A)', 'Huffman', 'H/L', 'lossless');
for q = 1:nq
    fprintf('%-15s %4d %7d %8.4f %8.4f %8.4f %9d\n', names{q}, numel(lvl{q}), Lfix(q), H(q), Lh(q), H(q)/Lh(q), ok(q));
end

fprintf('\n(e)-(g) Block coding (bits/symbol): H_k/k and Huffman L_k/k\n');
fprintf('%-15s %6s', 'Quantizer', 'Fixed');
for k = K, fprintf('  H_%d/%d  Huff_%d', k, k, k); end
fprintf('\n');
for q = 1:nq
    fprintf('%-15s %6d', names{q}, Lfix(q));
    for j = 1:numel(K), fprintf(' %7.4f %7.4f', Hk(q,j), Lk(q,j)); end
    fprintf('\n');
end
fprintf('All block codes lossless: %d\n', all(okb(:)));

fprintf('\n(b) Example, %s: counts %s, p(a) = %s\n', names{2}, mat2str(cnts{2}), mat2str(cnts{2}/N, 3));
fprintf('\n(f) Example, %s, k = 2:\n', names{1});
for i = 1:size(ex_ub, 1)
    fprintf('  block (%d,%d): count %4d, p = %.3f, Huffman codeword %s\n', ex_ub(i,1), ex_ub(i,2), ...
            round(ex_pb(i)*N/2), ex_pb(i), ex_cb{i});
end
fprintf('\n(h) Next symbol equal to the current one: %s\n', sprintf('%.0f%% ', 100*same));
fprintf('(h) k = 4: different blocks %s, seen once %s (out of %d blocks)\n', ...
        mat2str(nblk(:,end).'), mat2str(once(:,end).'), N/4);

fprintf('\n(i), (j) Full chain: N = %d symbols, fs = %.0f Hz (bits/symbol unless stated)\n', N, fs);
fprintf('%-15s %10s %6s %8s %7s %7s %7s %10s %10s\n', 'Quantizer', 'MSE', 'Fixed', 'Huffman', ...
        'k=2', 'k=3', 'k=4', 'Huff bits', 'Huff kb/s');
for q = 1:nq
    fprintf('%-15s %10.3e %6d %8.4f %7.4f %7.4f %7.4f %10d %10.2f\n', names{q}, mse_r(q), Lfix(q), ...
            Lh(q), Lk(q,1), Lk(q,2), Lk(q,3), round(Lh(q)*N), Lh(q)*fs/1e3);
end

fprintf('\n(k) LZ77 bits/symbol, W = %d, for L = %s\n', W, mat2str(Lz_list));
for q = 1:nq
    fprintf('%-15s %s   best L = %3d: %.4f\n', names{q}, sprintf('%8.4f', Lz_all(q,:)), Lbest(q), Llz(q));
end
fprintf('LZ77 lossless (all L): %d\n', all(okz(:)));

fprintf('\n(i) Latency (algorithmic delay = samples the encoder waits for / fs)\n');
fprintf('  Fixed length:            1 sample      = %.2f ms\n', 1e3/fs);
fprintf('  Huffman (p from all N):  %d samples = %.2f s\n', N, N/fs);
fprintf('  LZ77 (look-ahead L+1):   %s ms\n', sprintf('%.1f ', 1e3*(Lbest+1)/fs));

%% Figure: distortion versus rate, with and without Huffman coding
figure('Name', 'Distortion vs rate with Huffman');
semilogy(Lh(2:4), mse_r(2:4), 'o-', 'LineWidth', 1.2); hold on;
semilogy(Lh(5:7), mse_r(5:7), 's-', 'LineWidth', 1.2);
semilogy(Lfix(2:4), mse_r(2:4), 'o--', 'Color', [0.6 0.75 0.9]);
semilogy(Lfix(5:7), mse_r(5:7), 's--', 'Color', [0.95 0.7 0.6]);
semilogy(Lk(2:4,end), mse_r(2:4), 'o:', 'Color', [0 0.447 0.741], 'LineWidth', 1.2);
semilogy(Lk(5:7,end), mse_r(5:7), 's:', 'Color', [0.85 0.325 0.098], 'LineWidth', 1.2);
semilogy(Lh(1), mse_r(1), 'k^', 'MarkerFaceColor', 'k'); grid on;
xlim([0.5 5.5]); ylim([1e-4 1e-1]);
xlabel('Bits per symbol'); ylabel('MSE of the reconstructed message');
legend('uniform + Huffman', 'Lloyd-Max + Huffman', 'uniform, fixed length', ...
       'Lloyd-Max, fixed length', 'uniform + Huffman (blocks of 4)', ...
       'Lloyd-Max + Huffman (blocks of 4)', 'two-level');
title('Distortion versus rate, with and without Huffman coding');

%% Save all figures
if ~exist('figs_part3', 'dir'), mkdir figs_part3; end
fh = findobj('Type', 'figure');
for k = 1:numel(fh)
    name = regexprep(fh(k).Name, '[^\w]+', '_');
    exportgraphics(fh(k), fullfile('figs_part3', sprintf('Fig%02d_%s.png', fh(k).Number, name)), 'Resolution', 200);
end
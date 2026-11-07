clear; clc; close all;
%% 2 - Quantization of the sampled message

% message and its samples (as in section 1.2)
fs_c = 44100;
t = 0:1/fs_c:1.24;
m = piano_signal(t);
Mf = abs(fft(m)).^2;  Mf = Mf(1:floor(end/2));
f  = (0:numel(Mf)-1) * fs_c/numel(m);
B  = f(find(cumsum(Mf)/sum(Mf) >= 0.9999, 1));
fs = 4*B;
[ts, xs] = sample(t, m, fs);
mid = (t > 0.05) & (t < t(end)-0.05);          % skip the two ends of the reconstruction
zs = (ts >= 0.60) & (ts <= 0.62);              % 20 ms zoom
xmin = min(xs);  xmax = max(xs);

%% Quantizers: thr{q}, lvl{q}
names = {'Two-level', 'Uniform M=4', 'Uniform M=16', 'Uniform M=32', ...
         'Lloyd-Max M=4', 'Lloyd-Max M=16', 'Lloyd-Max M=32'};

% (a) two-level (optimal, M = 2): levels = conditional means, threshold = their midpoint
t2 = 0;
for it = 1:100
    l2 = [mean(xs(xs <= t2)), mean(xs(xs > t2))];
    t2 = mean(l2);
end
thr = {t2};
lvl = {l2};

% (c) uniform: M equal intervals over [min, max], levels at the midpoints
for M = [4 16 32]
    D = (xmax - xmin)/M;
    thr{end+1} = xmin + D*(1:M-1);
    lvl{end+1} = xmin + D*((1:M) - 0.5);
end

% (f) Lloyd-Max: start from the uniform levels, then repeat
%     thresholds = midpoints of the levels, levels = mean of the samples in each cell
for M = [4 16 32]
    l = xmin + (xmax - xmin)/M * ((1:M) - 0.5);
    for it = 1:2000
        idx = quan(xs, (l(1:end-1) + l(2:end))/2, 1:M);   % cell of every sample
        new = l;
        for k = 1:M
            if any(idx == k), new(k) = mean(xs(idx == k)); end
        end
        done = max(abs(new - l)) < 1e-10;
        l = new;
        if done, break; end
    end
    thr{end+1} = (l(1:end-1) + l(2:end))/2;
    lvl{end+1} = l;
end

%% Quantize, reconstruct, measure
fprintf('\n');
for q = 1:7
    xq{q} = quan(xs, thr{q}, lvl{q});
    mr{q} = reconstruct(t, xq{q}, fs);                       % rest of the chain ideal
    mse_s(q) = mean((xs - xq{q}).^2);                        % error of the samples
    mse_r(q) = mean((m(mid) - mr{q}(mid)).^2);               % error of the reconstructed message
    bits(q)  = log2(numel(lvl{q}));
    fprintf('%s\n  thresholds: %s\n  levels:     %s\n', names{q}, ...   % (d)
            sprintf('%.4f ', thr{q}), sprintf('%.4f ', lvl{q}));
end

fprintf('\n%-15s %5s %14s %16s\n', 'Quantizer', 'bits', 'MSE (samples)', 'MSE (reconstr.)');   % (e), (h)
for q = 1:7
    fprintf('%-15s %5.0f %14.3e %16.3e\n', names{q}, bits(q), mse_s(q), mse_r(q));
end

%% Perceived quality: ViSQOL, audio mode (needs the Audio Toolbox); MOS from 1 (bad) to 5 (excellent)
% m(t) is a formula, so the reference and the reconstructions are evaluated on a 48 kHz grid.
if exist('visqol', 'file')
    fs48 = 48000;
    t48  = 0:1/fs48:1.24;
    m48  = piano_signal(t48);
    mos0 = visqol(reconstruct(t48, xs, fs).', m48.', fs48);   % sampling only, no quantization
    fprintf('\nViSQOL (MOS, 1-5)\n%-15s %5.2f\n', 'No quantizer', mos0);
    for q = 1:7
        mos(q) = visqol(reconstruct(t48, xq{q}, fs).', m48.', fs48);
        fprintf('%-15s %5.2f\n', names{q}, mos(q));
    end
else
    disp('ViSQOL skipped: the Audio Toolbox is not installed.');
end

%% Figures
figure('Name', 'Two-level');                                 % (b)
plot(ts(zs), xs(zs), 'k.-'); hold on; stairs(ts(zs), xq{1}(zs), 'r'); grid on;
xlabel('Time (s)'); ylabel('Amplitude'); legend('original', 'two-level');
title(sprintf('Two-level quantization, 20 ms zoom (MSE = %.3e)', mse_s(1)));

figure('Name', 'Uniform');                                   % (c)
for k = 1:3
    subplot(3,1,k);
    plot(ts(zs), xs(zs), 'k.-'); hold on; stairs(ts(zs), xq{1+k}(zs), 'r'); grid on;
    ylabel('Amplitude'); legend('original', 'quantized');
    title(sprintf('%s (MSE = %.3e)', names{1+k}, mse_s(1+k)));
end
xlabel('Time (s)');

figure('Name', 'Lloyd-Max');                                 % (f)
[cnt, ctr] = hist(xs, 100);
dens = cnt / (sum(cnt) * (ctr(2) - ctr(1)));
for k = 1:3
    subplot(3,1,k);
    bar(ctr, dens, 1, 'FaceColor', [0.8 0.8 0.8]); hold on;
    for th = thr{4+k}, plot([th th], [0 max(dens)], 'r--'); end
    plot(lvl{4+k}, zeros(size(lvl{4+k})), 'go', 'MarkerFaceColor', 'g');
    ylabel('Density');
    title(sprintf('%s: thresholds (red) and levels (green) on the histogram of the samples', names{4+k}));
end
xlabel('Amplitude');

figure('Name', 'MSE vs rate');                               % (e), (h), (i)
semilogy(bits(2:4), mse_r(2:4), 'o-', 'LineWidth', 1.2); hold on;
semilogy(bits(5:7), mse_r(5:7), 's-', 'LineWidth', 1.2);
semilogy(bits(1), mse_r(1), 'k^', 'MarkerFaceColor', 'k'); grid on;
xlim([0.5 5.5]); ylim([1e-4 1e-1]); xticks(1:5);
xticklabels({'1 (M=2)', '2 (M=4)', '3', '4 (M=16)', '5 (M=32)'});
xlabel('Rate R = log_2 M (bits/sample)'); ylabel('MSE of the reconstructed message');
legend('uniform', 'Lloyd-Max', 'two-level'); title('Distortion versus rate');

save('quantized.mat', 'xs', 'fs', 'thr', 'lvl', 'names', 'mse_s', 'mse_r');   % used by part3.m

%% Listen: original, sampled + reconstructed, then every quantizer
mr0 = reconstruct(t, xs, fs);
clips  = [{m, mr0}, mr];
labels = [{'original', 'sampled + reconstructed'}, names];
for k = 1:numel(clips)
    fprintf('Playing: %s\n', labels{k});
    sound(clips{k}, fs_c);
    pause(numel(t)/fs_c + 0.5);
end

%% Save all figures
if ~exist('figs_part2', 'dir'), mkdir figs_part2; end
fh = findobj('Type', 'figure');
for k = 1:numel(fh)
    name = regexprep(fh(k).Name, '[^\w]+', '_');
    exportgraphics(fh(k), fullfile('figs_part2', sprintf('Fig%02d_%s.png', fh(k).Number, name)), 'Resolution', 200);
end
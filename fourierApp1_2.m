clear; clc; close all;
%% 1.2 - Message signal, Fourier series and DTFT

%% Message signal m(t)  (44.1 kHz grid stands in for continuous time)
fs_c = 44100;
t = 0:1/fs_c:1.24;
m = piano_signal(t);

%% Frequency content: bandwidth B = frequency below which 99.99% of the energy lies
Mf = abs(fft(m)).^2;  Mf = Mf(1:floor(end/2));
f  = (0:numel(Mf)-1) * fs_c/numel(m);
B  = f(find(cumsum(Mf)/sum(Mf) >= 0.9999, 1));
fN = 2*B;                         % Nyquist rate
fs = 2*fN;                        % sampling rate used from now on
fprintf('B = %.0f Hz, fN = %.0f Hz, fs = %.0f Hz\n', B, fN, fs);

% Second convention (x-dB bandwidth): last frequency where |M(f)|^2 is within
% X dB of its peak. Only for comparison; B above is the bandwidth used.
dBspec = 10*log10(Mf/max(Mf));
f40 = f(find(dBspec >= -40, 1, 'last'));
f50 = f(find(dBspec >= -50, 1, 'last'));
f_max = f50;                      % highest harmonic of the message (used for n in (c), (e))
fprintf('x-dB bandwidth: %.0f Hz (-40 dB line), %.0f Hz (-50 dB line); 2 x %.0f = %.0f Hz < fs\n', ...
        f40, f50, f50, 2*f50);

figure('Name', 'Message signal');
subplot(3,1,1); plot(t, m); grid on;
xlabel('Time (s)'); ylabel('m(t)'); title('Message signal m(t)');
subplot(3,1,2); z = (t >= 0.60) & (t <= 0.62); plot(t(z), m(z)); grid on;
xlabel('Time (s)'); ylabel('m(t)'); title('Zoom: 20 ms of the third note (E5)');
subplot(3,1,3); plot(f, 10*log10(Mf/max(Mf))); hold on; grid on;
xline(B, 'r--', sprintf('B = %.0f Hz', B), 'LabelOrientation', 'horizontal');
xlim([0 6000]); ylim([-100 0]);
xlabel('Frequency (Hz)'); ylabel('|M(f)|^2 (dB)'); title('Energy spectral density of m(t), relative to its peak');

%% Check: sample and reconstruct the whole message
[ts, xs] = sample(t, m, fs);
mr = reconstruct(t, xs, fs);
mid = (t > 0.05) & (t < t(end)-0.05);          % skip the two ends
fprintf('Samples: %d | reconstruction MSE: %.2e\n', numel(xs), mean((m(mid)-mr(mid)).^2));

figure('Name', 'Sampling check');
z  = (t >= 0.600) & (t <= 0.605);
zs = (ts >= 0.600) & (ts <= 0.605);
plot(t(z), m(z), 'k', 'LineWidth', 1.5); hold on; grid on;
plot(t(z), mr(z), 'r--', 'LineWidth', 1.2);
stem(ts(zs), xs(zs), 'filled');
legend('m(t)', 'reconstructed', 'samples');
xlabel('Time (s)'); ylabel('Amplitude');
title(sprintf('Sampling at f_s = %.0f Hz and reconstruction (5 ms zoom)', fs));
% sound(m, fs_c); pause(1.5); sound(mr, fs_c);   % listen: original, then rebuilt

%% Chunk x(t) of the message
L  = 0.05;                                     % chunk length (s)
t1 = 0.60;                                     % chunk start (inside the third note)
c  = (t >= t1) & (t < t1 + L);
tc = t(c);  x = m(c);

%% (a), (b) Fourier series for several n, with T = chunk length
T = L;
n_list = [10 20 33 66 100 126 168 200];   % chosen to reach the harmonics (20k Hz each)
n_show = [20 33 66 168];                 % approximations shown in the plot
E_n = zeros(size(n_list));
xh_show = zeros(numel(n_show), numel(x));
for i = 1:numel(n_list)
    xhat = ffs(x, tc, n_list(i), T);
    E_n(i) = trapz(tc, (x - xhat).^2);
    j = find(n_show == n_list(i));
    if ~isempty(j), xh_show(j,:) = xhat; end
end

% (a) largest Fourier coefficients: |c_k| peaks at the harmonics, k = f*T
[~, ck] = ffs(x, tc, 200, T);                  % c_k for k = -200..200
for k = [31 33 66 100]
    fprintf('|c_%d| = %.4f  (k/T = %.0f Hz)\n', k, abs(ck(200+1+k)), k/T);
end

figure('Name', 'Fourier series vs n');
for p = 1:2
    subplot(2,1,p);
    plot(tc, x, 'k', 'LineWidth', 1.5); hold on; grid on;
    plot(tc, xh_show, 'LineWidth', 1);
    xlabel('Time (s)'); ylabel('Amplitude');
    legend('x(t)', 'n = 20', 'n = 33', 'n = 66', 'n = 168');
end
subplot(2,1,1); title(sprintf('Fourier series of the %.0f ms chunk (T = %.0f ms)', 1e3*L, 1e3*T));
subplot(2,1,2); xlim([t1 t1+0.01]); title('Zoom: first 10 ms of the chunk');

%% (d) Error energy versus n (fixed T)
figure('Name', 'E vs n');
semilogy(n_list, E_n, 'o-', 'LineWidth', 1.2); grid on;
xlabel('n'); ylabel('E(n,T)'); title(sprintf('Error energy vs n (T = %.0f ms)', 1e3*T));
disp(table(n_list(:), E_n(:), 'VariableNames', {'n','E'}));

% (d) where the remaining error is (n = 168): near the two ends of the chunk
xhat = ffs(x, tc, 168, T);   e2 = (x - xhat).^2;
sA = (tc - tc(1)) < 0.003;   sB = (tc(end) - tc) < 0.003;   sM = ~(sA | sB);
Etot = trapz(tc, e2);
fprintf('n = 168: %.1f%% of E within 3 ms of the chunk ends (start %.1f%%, end %.1f%%), E in the middle = %.2e\n', ...
        100*(trapz(tc(sA), e2(sA)) + trapz(tc(sB), e2(sB)))/Etot, ...
        100*trapz(tc(sA), e2(sA))/Etot, 100*trapz(tc(sB), e2(sB))/Etot, trapz(tc(sM), e2(sM)));
fprintf('Chunk ends: x(start) = %.3f, x(end) = %.3f\n', x(1), x(end));
% (d) the step at the fundamentals: D#5 tail (k = 31), then E5 (k = 33)
for nn = [30 31 32 33]
    xhat = ffs(x, tc, nn, T);
    fprintf('n = %d: E / chunk energy = %.1f%%\n', nn, 100*trapz(tc, (x - xhat).^2)/trapz(tc, x.^2));
end

%% (c), (e) Error energy versus T (fixed, large n)
T_list = L * [0.6 0.8 0.9 1 1.2 1.5 2 2.5];      % 30 ... 125 ms
n_big  = ceil(f_max * max(T_list));            % sufficiently large for every T: n/T reaches f_max
T_show = L * [0.8 1 2];                        % approximations shown in the plot
E_T = zeros(size(T_list));
xh_show = zeros(numel(T_show), numel(x));
for j = 1:numel(T_list)
    xhat = ffs(x, tc, n_big, T_list(j));
    E_T(j) = trapz(tc, (x - xhat).^2);
    q = find(abs(T_show - T_list(j)) < 1e-12);
    if ~isempty(q), xh_show(q,:) = xhat; end
end

figure('Name', 'Fourier series vs T');
for p = 1:2
    subplot(2,1,p);
    plot(tc, x, 'k', 'LineWidth', 1.5); hold on; grid on;
    plot(tc, xh_show, 'LineWidth', 1);
    xlabel('Time (s)'); ylabel('Amplitude');
    legend('x(t)', 'T = 40 ms', 'T = 50 ms', 'T = 100 ms');
end
subplot(2,1,1); title(sprintf('Fourier series for different T (n = %d)', n_big));
subplot(2,1,2); xlim([t1+0.035 t1+0.05]); title('Zoom: last 15 ms of the chunk');

figure('Name', 'E vs T');
semilogy(1e3*T_list, E_T, 'o-', 'LineWidth', 1.2); grid on;
xlabel('T (ms)'); ylabel('E(n,T)'); title(sprintf('Error energy vs T (n = %d)', n_big));
disp(table(1e3*T_list(:), E_T(:), 'VariableNames', {'T_ms','E'}));
fprintf('Energy of the chunk: %.3e\n', trapz(tc, x.^2));

% (e) E(n,T) over both n and T: the trade-off between n and T
Tg = (30:5:125)/1e3;   ng = 0:10:500;
G = zeros(numel(ng), numel(Tg));
for j = 1:numel(Tg)
    for i = 1:numel(ng)
        xh = ffs(x, tc, ng(i), Tg(j));
        G(i,j) = trapz(tc, (x - xh).^2);
    end
end
figure('Name', 'E map');
imagesc(1e3*Tg, ng, log10(G)); axis xy; clim([-7 -2]); hold on;
ylabel(colorbar, 'log_{10} E(n,T)');
hl = plot(1e3*Tg, f_max*Tg, 'w--', 'LineWidth', 1.4);
xline(1e3*L, 'w:', 'T = chunk length', 'LineWidth', 1.2, 'LabelVerticalAlignment', 'top');
xlabel('T (ms)'); ylabel('n'); title('Error energy E(n,T)');
legend(hl, 'n = f_{max} T', 'TextColor', 'w', 'Color', 'none');

%% DTFT on frequency segments (third sampling method)
% The sampled message xs (N samples at fs) is described by K = N samples of its
% DTFT X(f_k), f_k = k*fs/K, computed directly from the definition:
%   X(f) = sum_n x[n] e^{-j 2 pi f n/fs}
% and recovered by   x[n] = (1/K) sum_k X(f_k) e^{j 2 pi f_k n/fs}.
N  = numel(xs);
K  = N;
fk = (0:K-1) * fs/K;
n  = 0:N-1;
Xd = zeros(1, K);
for k = 1:K
    Xd(k) = sum(xs .* exp(-1j*2*pi*fk(k)*n/fs));
end

figure('Name', 'DTFT');
half = fk <= fs/2;                                     % one-sided view: 0 .. fs/2
plot(fk(half), 20*log10(abs(Xd(half))/max(abs(Xd)))); grid on; hold on;
edges = [1000 1600 2300];                              % band edges (Hz)
for e = edges, xline(e, 'r--'); end
ylim([-100 0]); xlabel('Frequency (Hz)'); ylabel('|X(f)| (dB)');
title(sprintf('DTFT of the sampled message, K = %d frequency samples', K));

% Segment the frequency axis into bands and rebuild each band
fa = min(fk, fs - fk);                                 % f and fs-f belong to the same band
band_lims = [0 edges; edges fs/2]';                    % one row per band: [low high]
names = {'0-1000 Hz (fundamentals)', '1000-1600 Hz (2nd harmonics)', ...
         '1600-2300 Hz (3rd harmonics)', '2300 Hz - fs/2 (4th, 5th, rest)'};
nb = size(band_lims, 1);
xb = zeros(nb, N);
fprintf('\n');
for b = 1:nb
    in = fa >= band_lims(b,1) & fa < band_lims(b,2);
    if b == nb, in = fa >= band_lims(b,1); end         % last band includes fs/2
    acc = zeros(1, N);
    for k = find(in)
        acc = acc + Xd(k) * exp(1j*2*pi*fk(k)*n/fs);
    end
    xb(b,:) = real(acc) / K;                           % x[n] is real
    fprintf('Band %d, %-32s energy share: %8.4f %%\n', b, names{b}, ...
            100*sum(xb(b,:).^2)/sum(xs.^2));
end
fprintf('Sum of all bands vs x[n]: MSE = %.3e\n', mean((xs - sum(xb,1)).^2));

figure('Name', 'Bands');
zs = (ts >= 0.60) & (ts <= 0.62);
for b = 1:nb
    subplot(nb,1,b); plot(ts(zs), xb(b,zs)); grid on;
    ylabel('Amplitude'); title(['Band ' num2str(b) ': ' names{b}]);
end
xlabel('Time (s)');

figure('Name', 'Sum of bands');
plot(ts(zs), xs(zs), 'k', 'LineWidth', 1.5); hold on; grid on;
plot(ts(zs), sum(xb(:,zs),1), 'r--', 'LineWidth', 1.2);
legend('x[n]', 'sum of all bands');
xlabel('Time (s)'); ylabel('Amplitude');
title('Verification: x[n] vs sum of the DTFT bands (20 ms zoom)');

figure('Name', 'ESD linear'); plot(f, Mf/max(Mf)); xlim([0 6000]); grid on;
xlabel('Frequency (Hz)'); ylabel('|M(f)|^2 (linear, normalised)');
title('Energy spectral density of m(t), linear scale');

%% Save all figures and zip them
if ~exist('figs', 'dir'), mkdir figs; end
fh = findobj('Type', 'figure');
for k = 1:numel(fh)
    name = regexprep(fh(k).Name, '[^\w]+', '_');
    exportgraphics(fh(k), fullfile('figs', sprintf('Fig%02d_%s.png', fh(k).Number, name)), 'Resolution', 200);
end
zip('figs.zip', 'figs');
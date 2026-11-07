clear; clc; close all;
%% 1.1 - Generic Sampling and Reconstruction (test with a cosine)

f0       = 250;                 % Signal frequency (Hz)
T_total  = 0.04;                % Duration (s), 10 cycles
fs_dense = 100*f0;              % Dense grid standing in for continuous time
t  = 0:1/fs_dense:T_total;
xt = cos(2*pi*f0*t);

fN      = 2*f0;                 % Nyquist rate
fs_low  = 0.5*fN;               % Below Nyquist
fs_high = 2*fN;                 % Above Nyquist

% Sampling
[t_low,  x_low]  = sample(t, xt, fs_low);
[t_high, x_high] = sample(t, xt, fs_high);

% Reconstruction
xr_low  = reconstruct(t, x_low,  fs_low);
xr_high = reconstruct(t, x_high, fs_high);

%% Plotting
originalColor = [0.15 0.15 0.15];
reconColor    = [0.05 0.35 0.75];
sampleColor   = [0.80 0.15 0.12];

figure('Name','Sampling and Reconstruction', 'Position',[100 100 1150 750]);

% ----- Below Nyquist -----
subplot(2,1,1);
h1 = plot(t, xt, '-', 'Color', originalColor, 'LineWidth', 1.4); hold on;
h2 = plot(t, xr_low, '--', 'Color', reconColor, 'LineWidth', 2.2);
h3 = stem(t_low, x_low, 'filled', 'Color', sampleColor, ...
    'MarkerFaceColor', sampleColor, 'MarkerEdgeColor', sampleColor, ...
    'LineWidth', 1.3, 'MarkerSize', 8);
grid on;
title({sprintf('Below Nyquist: f_s = %.0f Hz = 0.5 f_N', fs_low), ...
       sprintf('Signal Frequency: f_0 = %.0f Hz | Nyquist Rate: f_N = %.0f Hz', f0, fN)});
xlabel('Time (s)'); ylabel('Amplitude');
legend([h1 h2 h3], 'Original x(t)', 'Reconstruction', 'Samples', 'Location', 'northeast');
xlim([0 T_total]); ylim([-1.3 1.3]);

% ----- Above Nyquist -----
subplot(2,1,2);
h1 = plot(t, xt, '-', 'Color', originalColor, 'LineWidth', 1.4); hold on;
h2 = plot(t, xr_high, '--', 'Color', reconColor, 'LineWidth', 2.2);
h3 = stem(t_high, x_high, 'filled', 'Color', sampleColor, ...
    'MarkerFaceColor', sampleColor, 'MarkerEdgeColor', sampleColor, ...
    'LineWidth', 1.3, 'MarkerSize', 7);
grid on;
title({sprintf('Above Nyquist: f_s = %.0f Hz = 2 f_N', fs_high), ...
       sprintf('Signal Frequency: f_0 = %.0f Hz | Nyquist Rate: f_N = %.0f Hz', f0, fN)});
xlabel('Time (s)'); ylabel('Amplitude');
legend([h1 h2 h3], 'Original x(t)', 'Reconstruction', 'Samples', 'Location', 'northeast');
xlim([0 T_total]); ylim([-1.3 1.3]);

%% Reconstruction error (MSE), middle half only:
% near the ends the sinc pulses of samples outside the window are missing
mid = (t > 0.25*T_total) & (t < 0.75*T_total);
mse_low  = mean((xt(mid) - xr_low(mid)).^2);
mse_high = mean((xt(mid) - xr_high(mid)).^2);

fprintf('\n=== Sampling and Reconstruction Results ===\n');
fprintf('Signal frequency f0: %.0f Hz\n', f0);
fprintf('Nyquist rate fN: %.0f Hz\n\n', fN);
fprintf('Below Nyquist (fs = %.0f Hz): MSE = %.6e\n', fs_low,  mse_low);
fprintf('Above Nyquist (fs = %.0f Hz): MSE = %.6e\n', fs_high, mse_high);
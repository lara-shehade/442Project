function m = piano_signal(t)
% Message signal m(t): first five notes of Fuer Elise (Beethoven), modelled on a real
% recording (S. Steppeler, Wikimedia Commons). Each note is a sum of decaying
% sinusoids whose frequencies, amplitudes and decay rates were measured from
% the recording. Defined by a formula, so it can be evaluated at any t (s).
%
% note i:  sum_k A_k * exp(-a_k*tau) * sin(2*pi*f_k*tau),  tau = t - t0(i)
% times a 2 ms attack, and a damper (extra decay r) when the key is released.

names = {'E5', 'D#5', 'E5', 'D#5', 'E5'};
t0    = [0.051 0.313 0.548 0.778 1.016];   % note onsets (s)
toff  = [0.272 0.245 0.240 0.248 0.244];   % key release after onset (s)
r     = 10;                      % damper decay rate after release (1/s)

% per note: rows = [frequency (Hz); amplitude; decay rate (1/s)]
P = cell(1, 5);
P{1} = [661.0 1324.0 1992.5 3355.0;   % E5
        0.1006 0.01205 0.009277 0.0004075;
        2.71 4.12 10.6 4.19];
P{2} = [623.9 1247.8 1879.0 2515.9 3163.1;   % D#5
        0.08762 0.01066 0.00174 0.001479 0.0005717;
        0.3 3.99 1.28 4.14 9.62];
P{3} = [661.0 1324.0 1992.2 3355.6;   % E5
        0.1499 0.02072 0.02239 0.001766;
        3.31 4.24 12 6.88];
P{4} = [623.6 1247.6 1878.9 2515.6 3162.7;   % D#5
        0.1437 0.02141 0.004923 0.004863 0.001983;
        0.441 3.82 1.91 4.37 7.27];
P{5} = [660.8 1324.1 1992.8 3354.9;   % E5
        0.125 0.01413 0.00799 0.0006339;
        2.95 3.63 5.02 3.67];

m = zeros(size(t));
for i = 1:numel(t0)
    tau = t - t0(i);
    on  = tau >= 0;
    ti  = tau(on);
    g   = (1 - exp(-ti/0.002)) .* exp(-r*max(ti - toff(i), 0));   % attack and damper
    f = P{i}(1,:);  A = P{i}(2,:);  a = P{i}(3,:);
    s = zeros(size(ti));
    for k = 1:numel(f)
        s = s + A(k) * exp(-a(k)*ti) .* sin(2*pi*f(k)*ti);
    end
    m(on) = m(on) + g .* s;
end
m = m / 0.2513;                    % scale so that max|m| = 1
end
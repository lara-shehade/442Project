function xrcon = reconstruct(t, x_sample, fs)
% Rebuilds the signal on t from its samples by ideal sinc interpolation.
if iscolumn(t),        t = t';               end
if iscolumn(x_sample), x_sample = x_sample'; end

Ts = 1/fs;
t_sample = min(t) + (0:length(x_sample)-1)*Ts;   % same instants as sample.m

xrcon = zeros(size(t));
for n = 1:length(x_sample)
    arg = (t - t_sample(n))/Ts;
    s = sin(pi*arg)./(pi*arg);
    s(arg == 0) = 1;                              % sinc(0) = 1
    xrcon = xrcon + x_sample(n)*s;
end
end
